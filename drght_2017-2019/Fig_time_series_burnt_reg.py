import os
import cartopy
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from scipy.interpolate import griddata, interp1d
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, OCEAN
from common_utils import *

def regrid_to_fire_map_resolution(fire_path, var_in, lat_in, lon_in, loc_lat=None, loc_lon=None, burn=0):

    # =========== Read in fire data ============
    fire_file  = Dataset(fire_path, mode='r')
    Burn_Date  = fire_file.variables['Burn_Date'][0:8,:,:]  # 2019-07 - 2020-02
    lat_out    = fire_file.variables['lat'][:]
    lon_out    = fire_file.variables['lon'][:]

    var_regrid = regrid_data(lat_in, lon_in, lat_out, lon_out, var_in, method='nearest')

    # burnt region from 2019-07 to 2020-02
    burn_area  = np.where( Burn_Date[0,:,:] + Burn_Date[1,:,:] + Burn_Date[2,:,:] + Burn_Date[3,:,:] +
                           Burn_Date[4,:,:] + Burn_Date[5,:,:] + Burn_Date[6,:,:] + Burn_Date[7,:,:] > 0, 1, Burn_Date[0,:,:])
    if burn == 1:
        # burnt region
        var_regrid = np.where(burn_area==1, var_regrid, np.nan )
    elif burn == 0:
        # all region
        var_regrid = var_regrid
    elif burn == -1:
        # unburnt region
        var_regrid = np.where(burn_area==0, var_regrid, np.nan )

    if loc_lat !=None:
        lons_2D, lats_2D = np.meshgrid(lon_out, lat_out)
        var_regrid = np.where(np.all(( lats_2D>loc_lat[0],
                                       lats_2D<loc_lat[1],
                                       lons_2D>loc_lon[0],
                                       lons_2D<loc_lon[1]), axis=0),
                                       var_regrid, np.nan)
        lat_out    = lats_2D
        lon_out    = lons_2D

    return var_regrid, lat_out, lon_out

def output_time_series_burn_region(file_out, fire_path, wrf_path, file_name, land_ctl_path, land_sen_path,
                          time_s=None, time_e=None, loc_lats=None, loc_lons=None,
                          lat_name=None, lon_name=None, message=None, burn=0):

    # Read all pixels in the dataset
    time_ctl, Tmax_ctl = read_var_multi_file([land_ctl_path +"Tair_f_inst/"+ file_name], "Tair_f_inst",
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)
    time_sen, Tmax_sen = read_var_multi_file([land_sen_path +"Tair_f_inst/"+ file_name], "Tair_f_inst",
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)

    time_ctl, LAI_ctl = read_var_multi_file([land_ctl_path +"LAI_inst/"+ file_name], 'LAI_inst',
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)
    time_sen, LAI_sen = read_var_multi_file([land_sen_path +"LAI_inst/"+ file_name], 'LAI_inst',
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)

    time_ctl, ALB_ctl = read_var_multi_file([land_ctl_path +"Albedo_inst/"+ file_name], 'Albedo_inst',
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)
    time_sen, ALB_sen = read_var_multi_file([land_sen_path +"Albedo_inst/"+ file_name], 'Albedo_inst',
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)

    time_ctl, SM_ctl = read_var_multi_file([land_ctl_path +"SoilMoist_inst/"+ file_name], 'SoilMoist_inst',
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)
    time_sen, SM_sen = read_var_multi_file([land_sen_path +"SoilMoist_inst/"+ file_name], 'SoilMoist_inst',
                                             loc_lat=None, loc_lon=None, lat_name=lat_name, lon_name=lon_name)

    SMtop_ctl    = (SM_ctl[:,0,:,:]*0.022 + SM_ctl[:,1,:,:]*0.058 + SM_ctl[:,2,:,:]*0.154 + SM_ctl[:,3,:,:]*0.266)/0.5
    SMtop_sen    = (SM_sen[:,0,:,:]*0.022 + SM_sen[:,1,:,:]*0.058 + SM_sen[:,2,:,:]*0.154 + SM_sen[:,3,:,:]*0.266)/0.5

    # Calculate daily values
    Tmax_daily_ctl = time_clip_to_day_max(time_ctl, Tmax_ctl, time_s, time_e, seconds=None)
    Tmax_daily_sen = time_clip_to_day_max(time_sen, Tmax_sen, time_s, time_e, seconds=None)

    LAI_daily_ctl  = time_clip_to_day(time_ctl, LAI_ctl, time_s, time_e, seconds=None)
    LAI_daily_sen  = time_clip_to_day(time_sen, LAI_sen, time_s, time_e, seconds=None)

    ALB_daily_ctl  = time_clip_to_day(time_ctl, ALB_ctl, time_s, time_e, seconds=None)
    ALB_daily_sen  = time_clip_to_day(time_sen, ALB_sen, time_s, time_e, seconds=None)

    SMtop_daily_ctl = time_clip_to_day(time_ctl, SMtop_ctl, time_s, time_e, seconds=None)
    SMtop_daily_sen = time_clip_to_day(time_sen, SMtop_sen, time_s, time_e, seconds=None)


    # Set lat and lon input
    wrf     = Dataset(wrf_path,  mode='r')
    lon_in  = wrf.variables['XLONG'][0,:,:]
    lat_in  = wrf.variables['XLAT'][0,:,:]

    ntime   = np.shape(Tmax_daily_ctl)[0]
    print("ntime =",ntime)

    for i in np.arange(ntime):
        print("i=",i)
        # regrid to burn map resolution ~ 400m
        if i == 0:
            Tmax_regrid_ctl_tmp, lats, lons  = regrid_to_fire_map_resolution(fire_path, Tmax_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            Tmax_regrid_sen_tmp, lats, lons  = regrid_to_fire_map_resolution(fire_path, Tmax_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            LAI_regrid_ctl_tmp, lats, lons   = regrid_to_fire_map_resolution(fire_path, LAI_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            LAI_regrid_sen_tmp, lats, lons   = regrid_to_fire_map_resolution(fire_path, LAI_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            ALB_regrid_ctl_tmp, lats, lons   = regrid_to_fire_map_resolution(fire_path, ALB_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            ALB_regrid_sen_tmp, lats, lons   = regrid_to_fire_map_resolution(fire_path, ALB_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            SMtop_regrid_ctl_tmp, lats, lons = regrid_to_fire_map_resolution(fire_path, SMtop_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            SMtop_regrid_sen_tmp, lats, lons = regrid_to_fire_map_resolution(fire_path, SMtop_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            # Set up array
            nlat = np.shape(Tmax_regrid_ctl_tmp)[0]
            nlon = np.shape(Tmax_regrid_ctl_tmp)[1]

            Tmax_regrid_ctl = np.zeros((ntime, nlat, nlon))
            Tmax_regrid_sen = np.zeros((ntime, nlat, nlon))

            LAI_regrid_ctl = np.zeros((ntime, nlat, nlon))
            LAI_regrid_sen = np.zeros((ntime, nlat, nlon))

            ALB_regrid_ctl = np.zeros((ntime, nlat, nlon))
            ALB_regrid_sen = np.zeros((ntime, nlat, nlon))

            SMtop_regrid_ctl = np.zeros((ntime, nlat, nlon))
            SMtop_regrid_sen = np.zeros((ntime, nlat, nlon))

            # Assign the first time step value
            Tmax_regrid_ctl[i,:,:]  = Tmax_regrid_ctl_tmp
            Tmax_regrid_sen[i,:,:]  = Tmax_regrid_sen_tmp

            LAI_regrid_ctl[i,:,:]   = LAI_regrid_ctl_tmp
            LAI_regrid_sen[i,:,:]   = LAI_regrid_sen_tmp

            ALB_regrid_ctl[i,:,:]   = ALB_regrid_ctl_tmp
            ALB_regrid_sen[i,:,:]   = ALB_regrid_sen_tmp

            SMtop_regrid_ctl[i,:,:] = SMtop_regrid_ctl_tmp
            SMtop_regrid_sen[i,:,:] = SMtop_regrid_sen_tmp

        else:
            Tmax_regrid_ctl[i,:,:], lats, lons = regrid_to_fire_map_resolution(fire_path, Tmax_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            Tmax_regrid_sen[i,:,:], lats, lons = regrid_to_fire_map_resolution(fire_path, Tmax_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            LAI_regrid_ctl[i,:,:], lats, lons  = regrid_to_fire_map_resolution(fire_path, LAI_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            LAI_regrid_sen[i,:,:], lats, lons  = regrid_to_fire_map_resolution(fire_path, LAI_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            ALB_regrid_ctl[i,:,:], lats, lons  = regrid_to_fire_map_resolution(fire_path, ALB_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            ALB_regrid_sen[i,:,:], lats, lons  = regrid_to_fire_map_resolution(fire_path, ALB_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

            SMtop_regrid_ctl[i,:,:], lats, lons = regrid_to_fire_map_resolution(fire_path, SMtop_daily_ctl[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)
            SMtop_regrid_sen[i,:,:], lats, lons = regrid_to_fire_map_resolution(fire_path, SMtop_daily_sen[i,:,:], lat_in, lon_in, loc_lat=None, loc_lon=None, burn=burn)

    # ===== Make masks for three regions =====
    # make fire lats and lons into 2 D
    lons_2D, lats_2D = np.meshgrid(lons, lats)
    mask_val         = np.zeros((3,np.shape(lons_2D)[0],np.shape(lons_2D)[1]),dtype=bool)

    for i in np.arange(3):
        mask_val[i,:,:]  = np.all(( lats_2D>loc_lats[i][0],lats_2D<loc_lats[i][1],
                                    lons_2D>loc_lons[i][0],lons_2D<loc_lons[i][1]), axis=0)

    # Extend the 3D mask (nreg, nlat, nlon) to 4D (nreg, ntime, nlat, nlon)
    mask_val_4D      = np.expand_dims(mask_val,axis=1).repeat(ntime,axis=1)

    print('np.shape(mask_val_4D)',np.shape(mask_val_4D))

    # Set up the output variables
    nreg          = 3
    Tmax_mean_ctl = np.zeros((nreg,ntime))
    Tmax_std_ctl  = np.zeros((nreg,ntime))
    Tmax_mean_sen = np.zeros((nreg,ntime))
    Tmax_std_sen  = np.zeros((nreg,ntime))

    LAI_mean_ctl = np.zeros((nreg,ntime))
    LAI_std_ctl  = np.zeros((nreg,ntime))
    LAI_mean_sen = np.zeros((nreg,ntime))
    LAI_std_sen  = np.zeros((nreg,ntime))

    ALB_mean_ctl = np.zeros((nreg,ntime))
    ALB_std_ctl  = np.zeros((nreg,ntime))
    ALB_mean_sen = np.zeros((nreg,ntime))
    ALB_std_sen  = np.zeros((nreg,ntime))

    SMtop_mean_ctl = np.zeros((nreg,ntime))
    SMtop_std_ctl  = np.zeros((nreg,ntime))
    SMtop_mean_sen = np.zeros((nreg,ntime))
    SMtop_std_sen  = np.zeros((nreg,ntime))

    # Mask out three regions
    for i in np.arange(3):
        print("process reg",i)

        Tmax_masked_ctl  = np.where( mask_val_4D[i,:,:,:], Tmax_regrid_ctl, np.nan)
        Tmax_masked_sen  = np.where( mask_val_4D[i,:,:,:], Tmax_regrid_sen, np.nan)

        LAI_masked_ctl   = np.where( mask_val_4D[i,:,:,:], LAI_regrid_ctl, np.nan)
        LAI_masked_sen   = np.where( mask_val_4D[i,:,:,:], LAI_regrid_sen, np.nan)

        ALB_masked_ctl   = np.where( mask_val_4D[i,:,:,:], ALB_regrid_ctl, np.nan)
        ALB_masked_sen   = np.where( mask_val_4D[i,:,:,:], ALB_regrid_sen, np.nan)

        SMtop_masked_ctl = np.where( mask_val_4D[i,:,:,:], SMtop_regrid_ctl, np.nan)
        SMtop_masked_sen = np.where( mask_val_4D[i,:,:,:], SMtop_regrid_sen, np.nan)

        Tmax_mean_ctl[i,:] = np.nanmean(Tmax_masked_ctl,axis=(1,2))
        Tmax_std_ctl[i,:]  = np.nanstd(Tmax_masked_ctl, axis=(1,2))
        Tmax_mean_sen[i,:] = np.nanmean(Tmax_masked_sen,axis=(1,2))
        Tmax_std_sen[i,:]  = np.nanstd(Tmax_masked_sen, axis=(1,2))

        LAI_mean_ctl[i,:]  = np.nanmean(LAI_masked_ctl,axis=(1,2))
        LAI_std_ctl[i,:]   = np.nanstd(LAI_masked_ctl, axis=(1,2))
        LAI_mean_sen[i,:]  = np.nanmean(LAI_masked_sen,axis=(1,2))
        LAI_std_sen[i,:]   = np.nanstd(LAI_masked_sen, axis=(1,2))

        ALB_mean_ctl[i,:]  = np.nanmean(ALB_masked_ctl,axis=(1,2))
        ALB_std_ctl[i,:]   = np.nanstd(ALB_masked_ctl, axis=(1,2))
        ALB_mean_sen[i,:]  = np.nanmean(ALB_masked_sen,axis=(1,2))
        ALB_std_sen[i,:]   = np.nanstd(ALB_masked_sen, axis=(1,2))

        SMtop_mean_ctl[i,:] = np.nanmean(SMtop_masked_ctl,axis=(1,2))
        SMtop_std_ctl[i,:]  = np.nanstd(SMtop_masked_ctl, axis=(1,2))
        SMtop_mean_sen[i,:] = np.nanmean(SMtop_masked_sen,axis=(1,2))
        SMtop_std_sen[i,:]  = np.nanstd(SMtop_masked_sen, axis=(1,2))

        if 1:
            fig1, ax1    = plt.subplots(nrows=1, ncols=1, figsize=[5,4],sharex=True, sharey=True, squeeze=True,
                                        subplot_kw={'projection': ccrs.PlateCarree()})

            states= NaturalEarthFeature(category="cultural", scale="50m",
                                                facecolor="none",
                                                name="admin_1_states_provinces_shp")

            # ======================= Set colormap =======================
            cmap    = plt.cm.BrBG
            cmap.set_bad(color='lightgrey')
            ax1.coastlines(resolution="50m",linewidth=1)
            ax1.set_extent([135,155,-39,-23])
            ax1.add_feature(states, linewidth=.5, edgecolor="black")

            ctl_avg = np.nanmean(ctl_regrid, axis=0)

            plot1  = ax1.contourf( lons, lats, np.nanmean(Tmax_masked_ctl,axis=0), transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
            cbar1  = plt.colorbar(plot1, ax=ax1, ticklocation="right", pad=0.08, orientation="horizontal",
                                aspect=40, shrink=0.6)
            cbar1.ax.tick_params(labelsize=8, labelrotation=45)

            plt.savefig('./plots/spatial_map_check_burn_region_'+str(i)+'.png',dpi=300)


    # ================== make output file ==================

    # create file and write global attributes
    f = nc.Dataset(file_out, 'w', format='NETCDF4')

    ### Create nc file ###
    f.history           = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date     = "%s" % (datetime.now())
    f.description       = '201910-202002 ΔTmax, LAI, Albedo and SM top 50cm in three burnt regions, made by MU Mengyuan'
    f.Conventions       = "CF-1.0"

    # set dimensions
    f.createDimension('region', 3)
    f.createDimension('time',  ntime)

    # Set cooridiates
    region               = f.createVariable('region', 'S', ('region'))
    region.standard_name = "Burnt regions"
    region.long_name     = "Name of the burnt regions"
    region[:]            = ['North','Central','South']

    time                 = f.createVariable('time', 'f4', ('time'))
    time.standard_name   = "time"
    time.units           = "days since 2019-10-01 00:00:00"
    time[:]              = np.arange((time_e-time_s).days)
    print("np.arange((time_e-time_s).days)",np.arange((time_e-time_s).days))

    tmax_mean_ctl               = f.createVariable('Tmax_mean_ctl', 'f4', ('region','time'))
    tmax_mean_ctl.standard_name = "Tmax in ctl"
    tmax_mean_ctl.units         = "deg C"
    tmax_mean_ctl[:]            = Tmax_mean_ctl

    tmax_std_ctl               = f.createVariable('Tmax_std_ctl', 'f4', ('region','time'))
    tmax_std_ctl.standard_name = "standard deviation of Tmax in burnt region in ctl"
    tmax_std_ctl.units         = "deg C"
    tmax_std_ctl[:]            = Tmax_std_ctl

    lai_mean_ctl               = f.createVariable('LAI_mean_ctl', 'f4', ('region','time'))
    lai_mean_ctl.standard_name = "LAI in ctl"
    lai_mean_ctl.units         = "m2 m-2"
    lai_mean_ctl[:]            = LAI_mean_ctl

    lai_std_ctl               = f.createVariable('LAI_std_ctl', 'f4', ('region','time'))
    lai_std_ctl.standard_name = "standard deviation of LAI in burnt region in ctl"
    lai_std_ctl.units         = "m2 m-2"
    lai_std_ctl[:]            = LAI_std_ctl

    alb_mean_ctl               = f.createVariable('ALB_mean_ctl', 'f4', ('region','time'))
    alb_mean_ctl.standard_name = "Albedo in ctl"
    alb_mean_ctl.units         = "-"
    alb_mean_ctl[:]            = ALB_mean_ctl

    alb_std_ctl               = f.createVariable('ALB_std_ctl', 'f4', ('region','time'))
    alb_std_ctl.standard_name = "standard deviation of albedo in burnt region in ctl"
    alb_std_ctl.units         = "-"
    alb_std_ctl[:]            = ALB_std_ctl

    smtop_mean_ctl               = f.createVariable('SMtop_mean_ctl', 'f4', ('region','time'))
    smtop_mean_ctl.standard_name = "Soil moisture at top 50cm in ctl"
    smtop_mean_ctl.units         = "m3 m-3"
    smtop_mean_ctl[:]            = SMtop_mean_ctl

    smtop_std_ctl               = f.createVariable('SMtop_std_ctl', 'f4', ('region','time'))
    smtop_std_ctl.standard_name = "standard deviation of top 50 cm soil moisutre in burnt region in ctl"
    smtop_std_ctl.units         = "m3 m-3"
    smtop_std_ctl[:]            = SMtop_std_ctl


    tmax_mean_sen               = f.createVariable('Tmax_mean_sen', 'f4', ('region','time'))
    tmax_mean_sen.standard_name = "Tmax in sen"
    tmax_mean_sen.units         = "deg C"
    tmax_mean_sen[:]            = Tmax_mean_sen

    tmax_std_sen               = f.createVariable('Tmax_std_sen', 'f4', ('region','time'))
    tmax_std_sen.standard_name = "standard deviation of Tmax in burnt region in sen"
    tmax_std_sen.units         = "deg C"
    tmax_std_sen[:]            = Tmax_std_sen

    lai_mean_sen               = f.createVariable('LAI_mean_sen', 'f4', ('region','time'))
    lai_mean_sen.standard_name = "LAI in sen"
    lai_mean_sen.units         = "m2 m-2"
    lai_mean_sen[:]            = LAI_mean_sen

    lai_std_sen               = f.createVariable('LAI_std_sen', 'f4', ('region','time'))
    lai_std_sen.standard_name = "standard deviation of LAI in burnt region in sen"
    lai_std_sen.units         = "m2 m-2"
    lai_std_sen[:]            = LAI_std_sen

    alb_mean_sen               = f.createVariable('ALB_mean_sen', 'f4', ('region','time'))
    alb_mean_sen.standard_name = "Albedo in sen"
    alb_mean_sen.units         = "-"
    alb_mean_sen[:]            = ALB_mean_sen

    alb_std_sen               = f.createVariable('ALB_std_sen', 'f4', ('region','time'))
    alb_std_sen.standard_name = "standard deviation of albedo in burnt region in sen"
    alb_std_sen.units         = "d"
    alb_std_sen[:]            = ALB_std_sen

    smtop_mean_sen               = f.createVariable('SMtop_mean_sen', 'f4', ('region','time'))
    smtop_mean_sen.standard_name = "Soil moisture at top 50cm in sen"
    smtop_mean_sen.units         = "m3 m-3"
    smtop_mean_sen[:]            = SMtop_mean_sen

    smtop_std_sen               = f.createVariable('SMtop_std_sen', 'f4', ('region','time'))
    smtop_std_sen.standard_name = "standard deviation of top 50 cm soil moisutre in burnt region in sen"
    smtop_std_sen.units         = "m3 m-3"
    smtop_std_sen[:]            = SMtop_std_sen

    f.close()


def plot_time_series_burn_region(file_out, time_s=None, time_e=None, message=None):
    # Check whether mask right

    f              = Dataset(file_out,  mode='r')
    time           = f.variables['time']
    Tmax_mean_ctl  = f.variables['Tmax_mean_ctl']
    LAI_mean_ctl   = f.variables['LAI_mean_ctl']
    LAI_std_ctl    = f.variables['LAI_std_ctl']
    ALB_mean_ctl   = f.variables['ALB_mean_ctl']
    ALB_std_ctl    = f.variables['ALB_std_ctl']
    SMtop_mean_ctl = f.variables['SMtop_mean_ctl']
    SMtop_std_ctl  = f.variables['SMtop_std_ctl']

    Tmax_mean_sen  = f.variables['Tmax_mean_sen']
    LAI_mean_sen   = f.variables['LAI_mean_sen']
    LAI_std_sen    = f.variables['LAI_std_sen']
    ALB_mean_sen   = f.variables['ALB_mean_sen']
    ALB_std_sen    = f.variables['ALB_std_sen']
    SMtop_mean_sen = f.variables['SMtop_mean_sen']
    SMtop_std_sen  = f.variables['SMtop_std_sen']

    df_reg1                  = pd.DataFrame({'Tmax_diff': Tmax_mean_sen[0,:] - Tmax_mean_ctl[0,:]})
    df_reg1['LAI_ctl_mean']  = LAI_mean_ctl[0,:]
    df_reg1['LAI_sen_mean']  = LAI_mean_sen[0,:]
    df_reg1['ALB_ctl_mean']  = ALB_mean_ctl[0,:]
    df_reg1['ALB_sen_mean']  = ALB_mean_sen[0,:]
    df_reg1['SMtop_ctl_mean']  = SMtop_mean_ctl[0,:]
    df_reg1['SMtop_sen_mean']  = SMtop_mean_sen[0,:]
    df_reg1['LAI_ctl_low']   = LAI_mean_ctl[0,:] - LAI_std_ctl[0,:]
    df_reg1['LAI_ctl_high']  = LAI_mean_ctl[0,:] + LAI_std_ctl[0,:]
    df_reg1['LAI_sen_low']   = LAI_mean_sen[0,:] - LAI_std_sen[0,:]
    df_reg1['LAI_sen_high']  = LAI_mean_sen[0,:] + LAI_std_sen[0,:]
    df_reg1['ALB_ctl_low']   = ALB_mean_ctl[0,:] - ALB_std_ctl[0,:]
    df_reg1['ALB_ctl_high']  = ALB_mean_ctl[0,:] + ALB_std_ctl[0,:]
    df_reg1['ALB_sen_low']   = ALB_mean_sen[0,:] - ALB_std_sen[0,:]
    df_reg1['ALB_sen_high']  = ALB_mean_sen[0,:] + ALB_std_sen[0,:]
    df_reg1['SMtop_ctl_low'] = SMtop_mean_ctl[0,:] - SMtop_std_ctl[0,:]
    df_reg1['SMtop_ctl_high']= SMtop_mean_ctl[0,:] + SMtop_std_ctl[0,:]
    df_reg1['SMtop_sen_low'] = SMtop_mean_sen[0,:] - SMtop_std_sen[0,:]
    df_reg1['SMtop_sen_high']= SMtop_mean_sen[0,:] + SMtop_std_sen[0,:]

    df_reg2                  = pd.DataFrame({'Tmax_diff': Tmax_mean_sen[1,:] - Tmax_mean_ctl[1,:]})
    df_reg2['LAI_ctl_mean']  = LAI_mean_ctl[1,:]
    df_reg2['LAI_sen_mean']  = LAI_mean_sen[1,:]
    df_reg2['ALB_ctl_mean']  = ALB_mean_ctl[1,:]
    df_reg2['ALB_sen_mean']  = ALB_mean_sen[1,:]
    df_reg2['SMtop_ctl_mean']  = SMtop_mean_ctl[1,:]
    df_reg2['SMtop_sen_mean']  = SMtop_mean_sen[1,:]
    df_reg2['LAI_ctl_low']   = LAI_mean_ctl[1,:] - LAI_std_ctl[1,:]
    df_reg2['LAI_ctl_high']  = LAI_mean_ctl[1,:] + LAI_std_ctl[1,:]
    df_reg2['LAI_sen_low']   = LAI_mean_sen[1,:] - LAI_std_sen[1,:]
    df_reg2['LAI_sen_high']  = LAI_mean_sen[1,:] + LAI_std_sen[1,:]
    df_reg2['ALB_ctl_low']   = ALB_mean_ctl[1,:] - ALB_std_ctl[1,:]
    df_reg2['ALB_ctl_high']  = ALB_mean_ctl[1,:] + ALB_std_ctl[1,:]
    df_reg2['ALB_sen_low']   = ALB_mean_sen[1,:] - ALB_std_sen[1,:]
    df_reg2['ALB_sen_high']  = ALB_mean_sen[1,:] + ALB_std_sen[1,:]
    df_reg2['SMtop_ctl_low'] = SMtop_mean_ctl[1,:] - SMtop_std_ctl[1,:]
    df_reg2['SMtop_ctl_high']= SMtop_mean_ctl[1,:] + SMtop_std_ctl[1,:]
    df_reg2['SMtop_sen_low'] = SMtop_mean_sen[1,:] - SMtop_std_sen[1,:]
    df_reg2['SMtop_sen_high']= SMtop_mean_sen[1,:] + SMtop_std_sen[1,:]

    df_reg3                  = pd.DataFrame({'Tmax_diff': Tmax_mean_sen[2,:] - Tmax_mean_ctl[2,:]})
    df_reg3['LAI_ctl_mean']  = LAI_mean_ctl[2,:]
    df_reg3['LAI_sen_mean']  = LAI_mean_sen[2,:]
    df_reg3['ALB_ctl_mean']  = ALB_mean_ctl[2,:]
    df_reg3['ALB_sen_mean']  = ALB_mean_sen[2,:]
    df_reg3['SMtop_ctl_mean']  = SMtop_mean_ctl[2,:]
    df_reg3['SMtop_sen_mean']  = SMtop_mean_sen[2,:]
    df_reg3['LAI_ctl_low']   = LAI_mean_ctl[2,:] - LAI_std_ctl[2,:]
    df_reg3['LAI_ctl_high']  = LAI_mean_ctl[2,:] + LAI_std_ctl[2,:]
    df_reg3['LAI_sen_low']   = LAI_mean_sen[2,:] - LAI_std_sen[2,:]
    df_reg3['LAI_sen_high']  = LAI_mean_sen[2,:] + LAI_std_sen[2,:]
    df_reg3['ALB_ctl_low']   = ALB_mean_ctl[2,:] - ALB_std_ctl[2,:]
    df_reg3['ALB_ctl_high']  = ALB_mean_ctl[2,:] + ALB_std_ctl[2,:]
    df_reg3['ALB_sen_low']   = ALB_mean_sen[2,:] - ALB_std_sen[2,:]
    df_reg3['ALB_sen_high']  = ALB_mean_sen[2,:] + ALB_std_sen[2,:]
    df_reg3['SMtop_ctl_low'] = SMtop_mean_ctl[2,:] - SMtop_std_ctl[2,:]
    df_reg3['SMtop_ctl_high']= SMtop_mean_ctl[2,:] + SMtop_std_ctl[2,:]
    df_reg3['SMtop_sen_low'] = SMtop_mean_sen[2,:] - SMtop_std_sen[2,:]
    df_reg3['SMtop_sen_high']= SMtop_mean_sen[2,:] + SMtop_std_sen[2,:]

    cleaner_dates = ["Oct 2019", "Nov 2019", "Dec 2019", "Jan 2020", "Feb 2020"]
    xtickslocs    = [         0,         31,         61,         92,       123 ]

    # ===================== Plotting =====================
    fig, axs = plt.subplots(nrows=4, ncols=3, figsize=[10,10], sharex=True,
                sharey=True, squeeze=True)
    plt.subplots_adjust(wspace=0., hspace=0.)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 12
    plt.rcParams['font.size']       = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    almost_black                    = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')

    time_steps = np.arange(len(time))

    Tmax_bot_val = -0.5
    Tmax_up_val  = 1.
    LAI_bot_val = 0
    LAI_up_val  = 6
    Albedo_bot_val = 0.05
    Albedo_up_val  = 0.15
    SM_bot_val = 0.
    SM_up_val  = 0.4

    # Tmax
    axs[0,0].axhline(y=0, color='grey', linestyle='--')
    axs[0,0].plot(df_reg1['Tmax_diff'].rolling(window=5).mean(), c = almost_black, lw=1., alpha=1)

    axs[0,1].axhline(y=0, color='grey', linestyle='--')
    axs[0,1].plot(df_reg2['Tmax_diff'].rolling(window=5).mean(), c = almost_black, lw=1., alpha=1)

    axs[0,2].axhline(y=0, color='grey', linestyle='--')
    axs[0,2].plot(df_reg3['Tmax_diff'].rolling(window=5).mean(), c = almost_black, lw=1., alpha=1)

    # LAI
    axs[1,0].fill_between(time_steps, df_reg1['LAI_ctl_low'].rolling(window=5).mean(), df_reg1['LAI_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[1,0].fill_between(time_steps, df_reg1['LAI_sen_low'].rolling(window=5).mean(), df_reg1['LAI_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[1,0].plot(df_reg1['LAI_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[1,0].plot(df_reg1['LAI_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs[1,1].fill_between(time_steps, df_reg2['LAI_ctl_low'].rolling(window=5).mean(), df_reg2['LAI_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[1,1].fill_between(time_steps, df_reg2['LAI_sen_low'].rolling(window=5).mean(), df_reg2['LAI_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[1,1].plot(df_reg2['LAI_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[1,1].plot(df_reg2['LAI_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs[1,2].fill_between(time_steps, df_reg3['LAI_ctl_low'].rolling(window=5).mean(), df_reg3['LAI_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[1,2].fill_between(time_steps, df_reg3['LAI_sen_low'].rolling(window=5).mean(), df_reg3['LAI_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[1,2].plot(df_reg3['LAI_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[1,2].plot(df_reg3['LAI_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    # ALB
    axs[2,0].fill_between(time_steps, df_reg1['ALB_ctl_low'].rolling(window=5).mean(), df_reg1['ALB_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[2,0].fill_between(time_steps, df_reg1['ALB_sen_low'].rolling(window=5).mean(), df_reg1['ALB_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[2,0].plot(df_reg1['ALB_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[2,0].plot(df_reg1['ALB_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs[2,1].fill_between(time_steps, df_reg2['ALB_ctl_low'].rolling(window=5).mean(), df_reg2['ALB_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[2,1].fill_between(time_steps, df_reg2['ALB_sen_low'].rolling(window=5).mean(), df_reg2['ALB_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[2,1].plot(df_reg2['ALB_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[2,1].plot(df_reg2['ALB_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs[2,2].fill_between(time_steps, df_reg3['ALB_ctl_low'].rolling(window=5).mean(), df_reg3['ALB_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[2,2].fill_between(time_steps, df_reg3['ALB_sen_low'].rolling(window=5).mean(), df_reg3['ALB_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[2,2].plot(df_reg3['ALB_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[2,2].plot(df_reg3['ALB_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    # SMtop
    axs[3,0].fill_between(time_steps, df_reg1['SMtop_ctl_low'].rolling(window=5).mean(), df_reg1['SMtop_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[3,0].fill_between(time_steps, df_reg1['SMtop_sen_low'].rolling(window=5).mean(), df_reg1['SMtop_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[3,0].plot(df_reg1['SMtop_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[3,0].plot(df_reg1['SMtop_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs[3,1].fill_between(time_steps, df_reg2['SMtop_ctl_low'].rolling(window=5).mean(), df_reg2['SMtop_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[3,1].fill_between(time_steps, df_reg2['SMtop_sen_low'].rolling(window=5).mean(), df_reg2['SMtop_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[3,1].plot(df_reg2['SMtop_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[3,1].plot(df_reg2['SMtop_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs[3,2].fill_between(time_steps, df_reg3['SMtop_ctl_low'].rolling(window=5).mean(), df_reg3['SMtop_ctl_high'].rolling(window=5).mean(),
                    color="green", edgecolor="none", alpha=0.3)
    axs[3,2].fill_between(time_steps, df_reg3['SMtop_sen_low'].rolling(window=5).mean(), df_reg3['SMtop_sen_high'].rolling(window=5).mean(),
                    color="orange", edgecolor="none", alpha=0.3)
    axs[3,2].plot(df_reg3['SMtop_ctl_mean'].rolling(window=5).mean(), label="ctl", c = "green", lw=0.5, alpha=1)
    axs[3,2].plot(df_reg3['SMtop_sen_mean'].rolling(window=5).mean(), label="exp", c = "orange", lw=0.5, alpha=1)

    axs.set_xlim(0,time_steps[-1])
    axs[0].set_ylim(Tmax_bot_val,Tmax_up_val)
    axs[1].set_ylim(LAI_bot_val,LAI_up_val)
    axs[2].set_ylim(Albedo_bot_val,Albedo_up_val)
    axs[3].set_ylim(SM_bot_val,SM_up_val)

    # Set top titles
    ax[0,0].set_title("North")
    ax[0,1].set_title("Central")
    ax[0,2].set_title("South")

    # Set left labels
    ax[0,0].set_ylabel("ΔT$\mathregular{_{max}}$ ($\mathregular{^{o}}$C)", fontsize=12)
    ax[1,0].set_ylabel("ΔLAI (m$\mathregular{^{2}}$ m$\mathregular{^{-2}}$)", fontsize=12)
    ax[2,0].set_ylabel("Δ$α$ (-)", fontsize=12)
    ax[3,0].set_ylabel("ΔSM (m$\mathregular{^{3}}$ m$\mathregular{^{-3}}$)", fontsize=12)

    axs.set(xticks=xtickslocs, xticklabels=cleaner_dates)
    # ax.legend()
    fig.tight_layout()

    plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # plot burnt region time series
    case_ctl       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
    case_sen       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"

    file_name      = 'LIS.CABLE.201701-202002.nc'
    fire_path      = '/g/data/w97/mm3972/data/MODIS/MODIS_fire/MCD64A1.061_500m_aid0001.nc'
    wrf_path       = "/scratch/w97/mm3972/model/NUWRF/Tinderbox_drght_LAI_ALB/output/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2/WRF_output/wrfout_d01_2017-02-01_06:00:00"

    land_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/LIS_output/"
    land_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/LIS_output/"

    file_out       = "/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/nc_files/times_series_201910-202002.nc"

    time_s         = datetime(2019,10,1,0,0,0,0)
    time_e         = datetime(2019,10,3,0,0,0,0)
    # time_e         = datetime(2020,3,1,0,0,0,0)

    #                   North ,        Central,       South
    loc_lats       = [[-32,-28.5],   [-34.5,-32.5], [-38,-34.5]]
    loc_lons       = [[151.5,153.5], [149.5,151.5], [146.5,151]]

    message        = "time_series_201910-202002_burnt"

    output_time_series_burn_region(file_out, fire_path, wrf_path, file_name, land_ctl_path, land_sen_path,
                        time_s=time_s, time_e=time_e, loc_lats=loc_lats, loc_lons=loc_lons,
                        lat_name="lat", lon_name="lon", message=message, burn=1)
