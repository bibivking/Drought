#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

'''
Functions:
1.

Climdex indices: https://www.climdex.org/learn/indices/
'''

from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import scipy.stats as stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.cm import get_cmap
from sklearn.metrics import mean_squared_error
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from wrf import (getvar, interplevel, get_cartopy, cartopy_xlim,
                 cartopy_ylim, to_np, latlon_coords)
from common_utils import *

def regrid_EHI_4_WRF_domain(EHI_path,wrf_path,EHI_out):

    # regridding EHF index
    time, Var  = read_var(EHI_path, 'event', loc_lat=None, loc_lon=None, lat_name='lat', lon_name='lon')
    time, lats = read_var(EHI_path, 'lat', loc_lat=None, loc_lon=None, lat_name='lat', lon_name='lon')
    time, lons = read_var(EHI_path, 'lon', loc_lat=None, loc_lon=None, lat_name='lat', lon_name='lon')

    wrf        = Dataset(wrf_path,  mode='r')
    lats_out   = wrf.variables['XLAT'][0,:,:]
    lons_out   = wrf.variables['XLONG'][0,:,:]

    nlat       = np.shape(lats_out)[0]
    nlon       = np.shape(lats_out)[1]
    ntime      = np.shape(Var)[0]
    var_regrid = np.zeros([ntime,nlat,nlon])
    for i in np.arange(ntime):
        var_regrid[i,:,:]= regrid_data(lats, lons, lats_out, lons_out, Var,method="nearest")

    # Make EHI output
    f_org = nc.Dataset(EHI_path, 'r', format='NETCDF4')

    # create file and write global attributes
    f     = nc.Dataset(EHI_out, 'w', format='NETCDF4')
    f.history           = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date     = "%s" % (datetime.now())
    f.description       = 'wrf output content '+var_name+', created by MU Mengyuan'

    # Copy global attributes from old file
    f_org = nc.Dataset(EHI_path, 'r', format='NETCDF4')
    for attr_name in f_org.ncattrs():
        attr_value = getattr(f_org, attr_name)
        setattr(f, attr_name, attr_value)
    # f_org.close()
    f.Conventions       = "CF-1.0"

    # set dimensions
    f.createDimension('time', None)
    f.createDimension('lat', nlat)
    f.createDimension('lon', nlon)

    time                = f.createVariable('time', 'f4', ('time'))
    time.standard_name  = "time"
    time.units          = "days since 1970-01-01"
    time[:]             = f_org.var


    lat                = f.createVariable('lat', 'f4', ('lat','lon'))
    lat.standard_name  = "latitude"
    lat.long_name      = "Latitude"
    lat.units          = "degrees_north"
    lat[:]             = lats_out

    lon                = f.createVariable('lon', 'f4', ('lat','lon'))
    lon.standard_name  = "longitude"
    lon.long_name      = "Longitude"
    lon.units          = "degrees_east"
    lon[:]             = lons_out


    event               = f.createVariable('event', 'f4', ('time','lat','lon'))
    event.FillValue     = 9.96921e+36
    event.missing_value = -999.99
    event.long_name     = "Event indicator"
    event.description   = "Indicates whether a summer heatwave is happening on that day"
    event[:]            = var_regrid

    f.close()
    f_org.close()

def read_spatial_data(land_ctl_path, land_sen_path, var_name, time_s=None,
                      time_e=None, lat_names="lat", lon_names="lon",loc_lat=None,
                      loc_lon=None, wrf_path=None):

    '''
    Read ctl and sen data
    '''

    print("var_name= "+var_name)

    if var_name in ["Tmax","Tmin",]:
        land_ctl_files= [land_ctl_path+'Tair_f_inst/LIS.CABLE.201701-201912.nc']
        land_sen_files= [land_sen_path+'Tair_f_inst/LIS.CABLE.201701-201912.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
    elif var_name in ["VegTmax","VegTmin"]:
        land_ctl_files= [land_ctl_path+'VegT_tavg/LIS.CABLE.201701-201912.nc']
        land_sen_files= [land_sen_path+'VegT_tavg/LIS.CABLE.201701-201912.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
    elif var_name in ["SurfTmax","SurfTmin"]:
        land_ctl_files= [land_ctl_path+'AvgSurfT_tavg/LIS.CABLE.201701-201912.nc']
        land_sen_files= [land_sen_path+'AvgSurfT_tavg/LIS.CABLE.201701-201912.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)

    # time-step into daily
    if var_name in ["SurfTmax","Tmax","VegTmax"]:
        # average of daily max
        ctl_in       = time_clip_to_day_max(time,Ctl_tmp,time_s,time_e)
        sen_in       = time_clip_to_day_max(time,Sen_tmp,time_s,time_e)
    elif var_name in ["SurfTmin","Tmin","VegTmin"]:
        # average of daily min
        ctl_in       = time_clip_to_day_min(time,Ctl_tmp,time_s,time_e)
        sen_in       = time_clip_to_day_min(time,Sen_tmp,time_s,time_e)

    return ctl_in, sen_in

def calc_heat_length():

def calc_heatwave_maganitude(EHI_out, land_ctl_path, land_sen_path, var_name, time_s=None, time_e=None, 
                             lat_names="lat", lon_names="lon",loc_lat=None, loc_lon=None):
    
    time_ehi, hw_event = read_var_multi_file(EHI_out, 'event', loc_lat, loc_lon, lat_names, lon_names)

    if var_name in ["Tmax","Tmin","T_DR"]:
        land_ctl_files= [land_ctl_path+'Tair_f_inst/LIS.CABLE.201701-201912.nc']
        land_sen_files= [land_sen_path+'Tair_f_inst/LIS.CABLE.201701-201912.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
    elif var_name in ["VegTmax","VegTmin","Veg_DR"]:
        land_ctl_files= [land_ctl_path+'VegT_tavg/LIS.CABLE.201701-201912.nc']
        land_sen_files= [land_sen_path+'VegT_tavg/LIS.CABLE.201701-201912.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
    elif var_name in ["SurfTmax","SurfTmin","Surf_DR"]:
        land_ctl_files= [land_ctl_path+'AvgSurfT_tavg/LIS.CABLE.201701-201912.nc']
        land_sen_files= [land_sen_path+'AvgSurfT_tavg/LIS.CABLE.201701-201912.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)

    # time-step into daily
    if var_name in ["SurfTmax","Tmax","VegTmax"]:
        # average of daily max
        ctl_in       = time_clip_to_day_max(time,Ctl_tmp,time_s,time_e)
        sen_in       = time_clip_to_day_max(time,Sen_tmp,time_s,time_e)
    elif var_name in ["SurfTmin","Tmin","VegTmin"]:
        # average of daily min
        ctl_in       = time_clip_to_day_min(time,Ctl_tmp,time_s,time_e)
        sen_in       = time_clip_to_day_min(time,Sen_tmp,time_s,time_e)
    elif var_name in ["Surf_DR","T_DR","Veg_DR"]:
        ctl_in       = time_clip_to_day_max(time,Ctl_tmp,time_s,time_e) - time_clip_to_day_min(time,Ctl_tmp,time_s,time_e)
        sen_in       = time_clip_to_day_max(time,Sen_tmp,time_s,time_e) - time_clip_to_day_min(time,Sen_tmp,time_s,time_e)

    Var_diff = ctl_in - sen_in
    
    # calculate HW periods values
    time_cood    = time_mask(time_ehi, time_s, time_e, seconds=None)
    hw_event_new = hw_event[time_cood,:,:]
    Var_diff     = np.where(hw_event_new==1, Var_diff, np.nan)
    var_diff     = np.nanmean(Var_diff,axis=0)

    return var_diff

def plot_spatial_map_hw_magnitude(land_ctl_path, land_sen_path, var_names, time_s=None,time_e=None, lat_names="lat", lon_names="lon",
                    loc_lat=None, loc_lon=None, wrf_path=None, message=None):
    
    hw_tmax_diff = calc_heatwave_maganitude(EHI_out, land_ctl_path, land_sen_path, var_names[0], time_s, time_e, 
                                            lat_names, lon_names,loc_lat, loc_lon)
    hw_tmin_diff = calc_heatwave_maganitude(EHI_out, land_ctl_path, land_sen_path, var_names[1], time_s, time_e, 
                                            lat_names, lon_names,loc_lat, loc_lon)
    hw_dr_diff   = calc_heatwave_maganitude(EHI_out, land_ctl_path, land_sen_path, var_names[2], time_s, time_e, 
                                            lat_names, lon_names,loc_lat, loc_lon)

    time, lats   = read_var(land_ctl_path, 'lat', loc_lat, loc_lon, lat_name='lat', lon_name='lon')
    time, lons   = read_var(land_ctl_path, 'lon', loc_lat, loc_lon, lat_name='lat', lon_name='lon')

    # ================== Start Plotting =================
    fig = plt.figure(figsize=(6,5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

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
    # choose colormap

    # =============== CHANGE HERE ===============
    cmap  = plt.cm.seismic_r
    clevs = np.arange(-1,1.1,0.1)

    # start plotting
    if loc_lat == None:
        ax.set_extent([135,155,-40,-25])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    ax.coastlines(resolution="50m",linewidth=1)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top   = False
    gl.ylabels_right = False
    gl.xlines        = False
    gl.xlines        = False

    if loc_lat == None:
        gl.xlocator  = mticker.FixedLocator([135,140,145,150,155])
        gl.ylocator  = mticker.FixedLocator([-40,-35,-30,-25])
    else:
        gl.xlocator  = mticker.FixedLocator(loc_lon)
        gl.ylocator  = mticker.FixedLocator(loc_lat)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}

    plt.contourf(lons, lats, hw_tmax_diff, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #

    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    cb.ax.tick_params(labelsize=10)
    plt.title(message, size=16)

    plt.savefig('./plots/spatial_map_hw_'+message+'.png',dpi=300)

if __name__ == "__main__":

    # ======================= Option =======================
    region = "SE Aus" #"SE Aus" #"CORDEX" #"SE Aus"

    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-25]
        loc_lon    = [135,155]
    elif region == "CORDEX":
        loc_lat    = [-52.36,3.87]
        loc_lon    = [89.25,180]

    # #################################
    # Plot WRF-CABLE vs AWAP temperal metrics
    # #################################

    EHI_path   = '/g/data/w97/mm3972/scripts/ehfheatwaves/nc_file'

    if 1:
        '''
        Test WRF-CABLE output
        '''

        case_name      = "ALB-CTL_new" #"bl_pbl2_mp4_sf_sfclay2" #
        case_ctl       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
        case_sen       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"
        land_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/LIS_output/"
        land_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/LIS_output/"
