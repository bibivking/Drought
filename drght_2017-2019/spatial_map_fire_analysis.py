import sys
import cartopy
import numpy as np
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

def plot_LAI_MODIS(LAI_MODIS_path):

    # =========== Read in data ===========
    LAI_file   = Dataset(LAI_MODIS_path, mode='r')
    LAI_input  = LAI_file.variables['LAI'][:]
    lat        = LAI_file.variables['latitude'][:]
    lon        = LAI_file.variables['longitude'][:]
    time       = nc.num2date( LAI_file.variables['time'][:], LAI_file.variables['time'].units,
                              only_use_cftime_datetimes=False, only_use_python_datetimes=True )

    # =========== Plotting ============
    # for i in np.arange( 1, len(time) ):

    #     print('time[i]', time[i])

    fig, ax    = plt.subplots(nrows=1, ncols=1, figsize=[5,4],sharex=True, sharey=True, squeeze=True,
                                subplot_kw={'projection': ccrs.PlateCarree()})

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
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

    states= NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    # ======================= Set colormap =======================
    cmap    = plt.cm.BrBG

    ax.coastlines(resolution="50m",linewidth=1)
    ax.set_extent([135,155,-39,-23])
    ax.add_feature(states, linewidth=.5, edgecolor="black")

    # Add gridlines
    gl              = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator(np.arange(125,160,1))
    gl.ylocator     = mticker.FixedLocator(np.arange(-40,-20,1))
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':12, 'color':almost_black}

    gl.xlabels_bottom = True
    gl.ylabels_left   = True

    clevs_diff     = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]

    # left: LAI_obs_mean
    plot           = ax.contourf( lon, lat, LAI_input[9] - LAI_input[1], levels=clevs_diff,
                                    transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
    # ax.text(0.02, 0.15, time[i], transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
    cbar = plt.colorbar(plot, ax=ax, ticklocation="right", pad=0.08, orientation="horizontal",
                        aspect=40, shrink=0.6)

    cbar.ax.tick_params(labelsize=8, labelrotation=45)

    plt.savefig('./plots/spatial_map_LAI_'+ str(time[9]) +'.png',dpi=300)

    return

def plot_fire_map(fire_path):

    day_in_month = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

    # =========== Read in data ===========
    fire_file  = Dataset(fire_path, mode='r')
    Burn_Date  = fire_file.variables['Burn_Date'][:]
    First_Day  = fire_file.variables['First_Day'][:]
    Last_Day   = fire_file.variables['Last_Day'][:]
    lat        = fire_file.variables['lat'][:]
    lon        = fire_file.variables['lon'][:]
    time       = nc.num2date( fire_file.variables['time'][:], fire_file.variables['time'].units,
                              only_use_cftime_datetimes=False, only_use_python_datetimes=True )

    Burn_Date  = np.where(Burn_Date<=0, np.nan, Burn_Date)
    First_Day  = np.where(First_Day<=0, np.nan, First_Day)
    Last_Day   = np.where(Last_Day<=0, np.nan, Last_Day)

    # ============= Plotting ==============
    for i in np.arange( len(time) ):

        print('time[i].month', time[i].month)

        fig, ax    = plt.subplots(nrows=1, ncols=1, figsize=[5,4],sharex=True, sharey=True, squeeze=True,
                                    subplot_kw={'projection': ccrs.PlateCarree()})

        plt.rcParams['text.usetex']     = False
        plt.rcParams['font.family']     = "sans-serif"
        plt.rcParams['font.serif']      = "Helvetica"
        plt.rcParams['axes.linewidth']  = 1.5
        plt.rcParams['axes.labelsize']  = 14
        plt.rcParams['font.size']       = 14
        plt.rcParams['legend.fontsize'] = 10
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14

        almost_black = '#262626'
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

        states= NaturalEarthFeature(category="cultural", scale="50m",
                                            facecolor="none",
                                            name="admin_1_states_provinces_shp")

        # ======================= Set colormap =======================
        cmap    = plt.cm.BrBG
        cmap.set_bad(color='lightgrey')

        ax.coastlines(resolution="50m",linewidth=1)
        ax.set_extent([135,155,-39,-23])
        ax.add_feature(states, linewidth=.5, edgecolor="black")

        # Add gridlines
        gl              = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
        gl.xlabels_top  = False
        gl.ylabels_right= False
        gl.xlines       = False
        gl.ylines       = False
        gl.xlocator     = mticker.FixedLocator(np.arange(125,160,1))
        gl.ylocator     = mticker.FixedLocator(np.arange(-40,-20,1))
        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
        gl.ylabel_style = {'size':12, 'color':almost_black}

        gl.xlabels_bottom = True
        gl.ylabels_left   = True

        month_int = int(time[i].month)

        # clevs = np.arange(day_in_month[month_int-1]+1, day_in_month[month_int]+1)
        # clevs_diff = [-5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]

        plot  = ax.contourf( lon, lat, Burn_Date[i,:,:]-day_in_month[month_int-1],transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #levels=clevs, # levels=clevs_diff,
        # ax.text(0.02, 0.15, time[i], transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax.add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
        cbar  = plt.colorbar(plot, ax=ax, ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.6)

        cbar.ax.tick_params(labelsize=8, labelrotation=45)

        plt.savefig('./plots/spatial_map_fire_Burn_Date_'+ str(time[i]) +'.png',dpi=300)

    return

def plot_LAI_fire_map(fire_path,LAI_MODIS_path):

    # =========== Read in fire data ============
    fire_file  = Dataset(fire_path, mode='r')
    Burn_Date  = fire_file.variables['Burn_Date'][4:8,:,:]  # 2019-11 - 2020-02
    lat_out    = fire_file.variables['lat'][:]
    lon_out    = fire_file.variables['lon'][:]

    time_fire  = nc.num2date( fire_file.variables['time'][:], fire_file.variables['time'].units,
                              only_use_cftime_datetimes=False, only_use_python_datetimes=True )
    nlat       = len(lat_out)
    nlon       = len(lon_out)

    # =========== Read in MODIS data ===========
    LAI_file   = Dataset(LAI_MODIS_path, mode='r')
    LAI_in     = LAI_file.variables['LAI'][:]
    lat_in     = LAI_file.variables['latitude'][:]
    lon_in     = LAI_file.variables['longitude'][:]
    time       = nc.num2date( LAI_file.variables['time'][:], LAI_file.variables['time'].units,
                              only_use_cftime_datetimes=False, only_use_python_datetimes=True )
    ntime      = len(time)
    LAI_regrid = np.zeros((ntime,nlat,nlon))

    for i in np.arange(ntime):
        LAI_regrid[i,:,:] = regrid_data(lat_in, lon_in, lat_out, lon_out, LAI_in[i,:,:], method='nearest',threshold=0)

    # ============= Plotting ==============
    fig, ax    = plt.subplots(nrows=1, ncols=1, figsize=[5,4],sharex=True, sharey=True, squeeze=True,
                                subplot_kw={'projection': ccrs.PlateCarree()})

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
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

    states= NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    # ======================= Set colormap =======================
    cmap    = plt.cm.BrBG
    cmap.set_bad(color='lightgrey')

    ax.coastlines(resolution="50m",linewidth=1)
    ax.set_extent([135,155,-39,-23])
    ax.add_feature(states, linewidth=.5, edgecolor="black")

    # Add gridlines
    gl              = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.xlabels_top  = False
    gl.ylabels_right= False
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator(np.arange(125,160,1))
    gl.ylocator     = mticker.FixedLocator(np.arange(-40,-20,1))
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':12, 'color':almost_black}

    gl.xlabels_bottom = True
    gl.ylabels_left   = True

    clevs_diff     = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5]

    # Nov changes
    # LAI_plot       = np.where(Burn_Date[0,:,:] + Burn_Date[1,:,:]> 0, LAI_regrid[9,:,:] - LAI_regrid[0,:,:], np.nan)
    LAI_plot       = np.where( Burn_Date[1,:,:]> 0, LAI_regrid[9,:,:] - LAI_regrid[4,:,:], np.nan)
    plot           = ax.contourf( lon_out, lat_out, LAI_plot, levels=clevs_diff, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #

    # ax.text(0.02, 0.15, time[i], transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    ax.add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
    cbar = plt.colorbar(plot, ax=ax, ticklocation="right", pad=0.08, orientation="horizontal",
                        aspect=40, shrink=0.6)

    cbar.ax.tick_params(labelsize=8, labelrotation=45)

    plt.savefig('./plots/spatial_map_fire_Burnt_LAI_Dec.png',dpi=300)

def read_LIS_diff(var_name,land_ctl_path,land_sen_path, lat_names, lon_names):

    print("plotting "+var_name)

    if var_name in ["Tmax","Tmin","TDR"]:
        land_ctl_files= [land_ctl_path+'Tair_f_inst/LIS.CABLE.201912-202002.nc']
        land_sen_files= [land_sen_path+'Tair_f_inst/LIS.CABLE.201912-202002.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
        Ctl_tmp       = Ctl_tmp -273.15
        Sen_tmp       = Sen_tmp -273.15
    elif var_name in ["VegTmax","VegTmin","VegTDR"]:
        land_ctl_files= [land_ctl_path+'VegT_tavg/LIS.CABLE.201912-202002.nc']
        land_sen_files= [land_sen_path+'VegT_tavg/LIS.CABLE.201912-202002.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        Ctl_tmp       = Ctl_tmp -273.15
        Sen_tmp       = Sen_tmp -273.15
    elif var_name in ["SurfTmax","SurfTmin","SurfTDR"]:
        land_ctl_files= [land_ctl_path+'AvgSurfT_tavg/LIS.CABLE.201912-202002.nc']
        land_sen_files= [land_sen_path+'AvgSurfT_tavg/LIS.CABLE.201912-202002.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)
        Ctl_tmp       = Ctl_tmp -273.15
        Sen_tmp       = Sen_tmp -273.15
    elif var_name in ["Rnet",]:
        land_ctl_files= [land_ctl_path+'Lwnet_tavg/LIS.CABLE.201912-202002.nc']
        land_sen_files= [land_sen_path+'Lwnet_tavg/LIS.CABLE.201912-202002.nc']
        time, Ctl_Lwnet_tmp = read_var_multi_file(land_ctl_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_Lwnet_tmp = read_var_multi_file(land_sen_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        land_ctl_files= [land_ctl_path+'Swnet_tavg/LIS.CABLE.201912-202002.nc']
        land_sen_files= [land_sen_path+'Swnet_tavg/LIS.CABLE.201912-202002.nc']
        time, Ctl_Swnet_tmp = read_var_multi_file(land_ctl_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_Swnet_tmp = read_var_multi_file(land_sen_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        Ctl_tmp = Ctl_Lwnet_tmp+Ctl_Swnet_tmp
        Sen_tmp = Sen_Lwnet_tmp+Sen_Swnet_tmp
    elif var_name in ["SM_top50cm",]:
        land_ctl_files = [land_ctl_path+'SoilMoist_inst/LIS.CABLE.201912-202002.nc']
        land_sen_files = [land_sen_path+'SoilMoist_inst/LIS.CABLE.201912-202002.nc']
        time, Ctl_temp = read_var_multi_file(land_ctl_files, 'SoilMoist_inst', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_temp = read_var_multi_file(land_sen_files, 'SoilMoist_inst', loc_lat, loc_lon, lat_names, lon_names)
        # [.022, .058, .154, .409, 1.085, 2.872]
        Ctl_tmp    = Ctl_temp[:,0,:,:]*0.022 + Ctl_temp[:,1,:,:]*0.058 + Ctl_temp[:,2,:,:]*0.154 + Ctl_temp[:,3,:,:]*0.266
        Sen_tmp    = Sen_temp[:,0,:,:]*0.022 + Sen_temp[:,1,:,:]*0.058 + Sen_temp[:,2,:,:]*0.154 + Sen_temp[:,3,:,:]*0.266
    elif var_name in ['VPD','VPDmax','VPDmin']:
        tair_ctl_files= [land_ctl_path+'Tair_f_inst/LIS.CABLE.201912-202002.nc']
        tair_sen_files= [land_sen_path+'Tair_f_inst/LIS.CABLE.201912-202002.nc']
        qair_ctl_files= [land_ctl_path+'Qair_f_inst/LIS.CABLE.201912-202002.nc']
        qair_sen_files= [land_sen_path+'Qair_f_inst/LIS.CABLE.201912-202002.nc']
        pres_ctl_files= ['/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2/'
                            +'WRF_output/slp/wrfout_201701-202002.nc']
        pres_sen_files= ['/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB/'
                            +'WRF_output/slp/wrfout_201701-202002.nc']

        time, Tair_ctl    = read_var_multi_file(tair_ctl_files, "Tair_f_inst", loc_lat, loc_lon, lat_names, lon_names)
        time, Tair_sen    = read_var_multi_file(tair_sen_files, "Tair_f_inst", loc_lat, loc_lon, lat_names, lon_names)
        time, Qair_ctl    = read_var_multi_file(qair_ctl_files, "Qair_f_inst", loc_lat, loc_lon, lat_names, lon_names)
        time, Qair_sen    = read_var_multi_file(qair_sen_files, "Qair_f_inst", loc_lat, loc_lon, lat_names, lon_names)
        time_wrf, Pres_ctl_tmp= read_var_multi_file(pres_ctl_files, "slp", loc_lat, loc_lon, lat_names, lon_names)
        time_wrf, Pres_sen_tmp= read_var_multi_file(pres_sen_files, "slp", loc_lat, loc_lon, lat_names, lon_names)

        time_in = []
        time_out= []
        for t in time_wrf:
            time_in.append(t.total_seconds())
        for t in time:
            time_out.append(t.total_seconds())

        f_ctl             = interp1d(np.array(time_in), Pres_ctl_tmp[:], kind='linear',fill_value='extrapolate', axis=0)
        f_sen             = interp1d(np.array(time_in), Pres_sen_tmp[:],kind='linear', fill_value='extrapolate', axis=0)
        Pres_ctl          = f_ctl(np.array(time_out))
        Pres_sen          = f_sen(np.array(time_out))
        Ctl_tmp           = qair_to_vpd(Qair_ctl, Tair_ctl, Pres_ctl)
        Sen_tmp           = qair_to_vpd(Qair_sen, Tair_sen, Pres_sen)
    else:
        land_ctl_files= [land_ctl_path+var_name+'/LIS.CABLE.201912-202002.nc']
        land_sen_files= [land_sen_path+var_name+'/LIS.CABLE.201912-202002.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, var_name, loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, var_name, loc_lat, loc_lon, lat_names, lon_names)

    if 'max' in var_name:
        # average of daily max
        ctl_in       = spital_var_max(time,Ctl_tmp,time_s,time_e)
        sen_in       = spital_var_max(time,Sen_tmp,time_s,time_e)
    elif 'min' in var_name:
        # average of daily min
        ctl_in       = spital_var_min(time,Ctl_tmp,time_s,time_e)
        sen_in       = spital_var_min(time,Sen_tmp,time_s,time_e)
    elif 'TDR' in var_name:
        # average of daily min
        ctl_in_max   = spital_var_max(time,Ctl_tmp,time_s,time_e)
        sen_in_max   = spital_var_max(time,Sen_tmp,time_s,time_e)
        ctl_in_min   = spital_var_min(time,Ctl_tmp,time_s,time_e)
        sen_in_min   = spital_var_min(time,Sen_tmp,time_s,time_e)
        ctl_in       = ctl_in_max - ctl_in_min
        sen_in       = sen_in_max - sen_in_min
    else:
        ctl_in       = spital_var(time,Ctl_tmp,time_s,time_e)
        sen_in       = spital_var(time,Sen_tmp,time_s,time_e)

    if var_name in ['WaterTableD_tavg','WatTable']:
        ctl_in     = ctl_in/1000.
        sen_in     = sen_in/1000.
    if var_name in ['ESoil_tavg','Evap_tavg',"ECanop_tavg",'TVeg_tavg',"Rainf_tavg","Snowf_tavg","Qs_tavg","Qsb_tavg"]:
        t_s        = time_s - datetime(2000,1,1,0,0,0,0)
        t_e        = time_e - datetime(2000,1,1,0,0,0,0)
        ctl_in     = ctl_in*3600*24*(t_e.days - t_s.days)
        sen_in     = sen_in*3600*24*(t_e.days - t_s.days)
    if var_name in ['Qair_f_inst']:
        ctl_in     = ctl_in*1000
        sen_in     = sen_in*1000
    if var_name in ['GPP_tavg','NPP_tavg']:
        t_s        = time_s - datetime(2000,1,1,0,0,0,0)
        t_e        = time_e - datetime(2000,1,1,0,0,0,0)
        s2d        = 3600*24.          # s-1 to d-1
        GPP_scale  = -0.000001*12*s2d   # umol s-1 to g d-1
        ctl_in     = ctl_in*GPP_scale*(t_e.days - t_s.days)
        sen_in     = sen_in*GPP_scale*(t_e.days - t_s.days)

    var_diff     = sen_in - ctl_in


    # Select clevs
    cmap  = plt.cm.BrBG
    if var_name in ['WaterTableD_tavg','WatTable']:
        clevs = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]
    elif var_name in ['GWwb_tavg','GWMoist']:
        clevs = [-0.05,-0.04,-0.03,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.03,0.04,0.05]
    elif  var_name in ["Qair_f_inst"]:
        clevs = [-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]
        #clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2]
    elif var_name in ['SoilMoist_inst','SoilMoist',"SM_top50cm"]:
        clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]
    elif var_name in ['ESoil_tavg','Evap_tavg',"ECanop_tavg",'TVeg_tavg',"Rainf_tavg","Snowf_tavg","Qs_tavg","Qsb_tavg"]:
        # clevs = [-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30]
        # clevs = [-5.,-4.5,-4.,-3.5,-3.,-2.5,-2,-1.5,-1,-0.5,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
        clevs = [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,-5,5,10,20.,30,40,50,60,70,80,90,100]
        # clevs = [-140,-120,-100,-80,-60,-40,-20,20,40,60,80,100,120,140]
    elif var_name in ["GPP_tavg","NPP_tavg",]:
        clevs = [-200,-190,-180,-170,-160,-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,
                    -5,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]
        # clevs = [-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
        clevs_percentage =  [-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70]
        cmap  = plt.cm.BrBG
    elif var_name in ["CanopInt_inst","SnowCover_inst"]:
        clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.]
        # clevs = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]
    elif var_name in ["Qle_tavg","Qh_tavg","Qg_tavg","Rnet",]:
        clevs = [-140,-120,-100,-80,-60,-40,-20,-10,-5,5,10,20,40,60,80,100,120,140]
    elif var_name in ["Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst","Rnet"]:
        clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
    elif var_name in ["Wind_f_inst",]:
        clevs = [-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4]
    elif var_name in ["Psurf_f_inst"]:
        clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
    elif var_name in ["Tair_f_inst","Tmax","Tmin","VegT_tavg","VegTmax","VegTmin",
                        "AvgSurfT_tavg","SurfTmax","SurfTmin","SoilTemp_inst",'TDR','VegTDR','SurfTDR']:
        # clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.]
        clevs = [-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2]
        # clevs = [-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1.,1.5,2,2.5,3.]
        cmap  = plt.cm.coolwarm
    elif var_name in ["Wind_f_inst",]:
        clevs = [-2.,-1.5,-1,-0.5,-0.1,0.1,0.5,1.,1.5,2.]
    elif var_name in ["FWsoil_tavg","SmLiqFrac_inst","SmFrozFrac_inst"]:
        clevs = [-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5]
    elif var_name in ["LAI_inst"]:
        clevs = [-2,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2]
        clevs_percentage =  [-70,-60,-50,-40,-30,-20,-10,-5,5,10,20,30,40,50,60,70]
        cmap  = plt.cm.BrBG
    elif var_name in ["VPD","VPDmax","VPDmin",]:
        clevs = [-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
        cmap  = plt.cm.BrBG
    elif var_name in ["Albedo_inst"]:
        clevs = [-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08]
        clevs_percentage =   [-70,-60,-50,-40,-30,-20,-10,-5,5,10,20,30,40,50,60,70]
        cmap  = plt.cm.BrBG_r
    else:
        clevs = [-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5]

    return sen_in, var_diff, clevs, cmap

def regrid_to_fire_map_resolution(fire_path, var_in, lat_in, lon_in, burn=0):

    # =========== Read in fire data ============
    fire_file  = Dataset(fire_path, mode='r')
    Burn_Date  = fire_file.variables['Burn_Date'][0:8,:,:]  # 2019-07 - 2020-02
    lat_out    = fire_file.variables['lat'][:]
    lon_out    = fire_file.variables['lon'][:]

    var_regrid = regrid_data(lat_in, lon_in, lat_out, lon_out, var_in, method='nearest')

    # burnt region from 2019-07 to 2020-01
    burn_area  = np.where( Burn_Date[0,:,:] + Burn_Date[1,:,:] + Burn_Date[2,:,:] + Burn_Date[3,:,:] +
                           Burn_Date[4,:,:] + Burn_Date[5,:,:] + Burn_Date[6,:,:] >0, 1, np.nan)
    if burn == 1:
        # burnt region
        var_regrid = np.where(burn_area==1, var_regrid, np.nan )
    elif burn == 0:
        # all region
        var_regrid = var_regrid
    elif burn == -1:
        # unburnt region
        var_regrid = np.where(np.isnan(burn_area), var_regrid, np.nan )

    return var_regrid, lat_out, lon_out

def plot_LIS_diff(fire_path, land_ctl_path, land_sen_path, var_names, time_s=None, time_e=None,
                  lat_names="lat", lon_names="lon",loc_lat=None, loc_lon=None, wrf_path=None, message=None, burn=0):

    '''
    plot LIS variables in burnt / unburnt / all regions
    '''

    # Read in WRF lat and lon
    wrf            = Dataset(wrf_path,  mode='r')
    lon_in         = wrf.variables['XLONG'][0,:,:]
    lat_in         = wrf.variables['XLAT'][0,:,:]

    for var_name in var_names:

        # read in var
        sen_in, var_diff, clevs, cmap = read_LIS_diff(var_name, land_ctl_path, land_sen_path, lat_names, lon_names)

        # regrid to burn map resolution ~ 400m
        var_regrid, lats, lons = regrid_to_fire_map_resolution(fire_path, var_diff, lat_in, lon_in, burn=burn)
        sen_regrid, lats, lons = regrid_to_fire_map_resolution(fire_path, sen_in, lat_in, lon_in, burn=burn)

        # ================== Start Plotting =================
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=[10,6],sharex=True,
                    sharey=True, squeeze=True, subplot_kw={'projection': ccrs.PlateCarree()})

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

        states= NaturalEarthFeature(category="cultural", scale="50m",
                                            facecolor="none",
                                            name="admin_1_states_provinces_shp")

        # =============== CHANGE HERE ===============

        for i in np.arange(2):

            axs[i].coastlines(resolution="50m",linewidth=1)
            axs[i].set_extent([135,155,-39,-23])
            axs[i].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = axs[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            # gl.xlocator     = mticker.FixedLocator(np.arange(125,160,1))
            # gl.ylocator     = mticker.FixedLocator(np.arange(-40,-20,1))
            gl.xlocator     = mticker.FixedLocator([130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25,-20])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            gl.xlabels_bottom = True
            gl.ylabels_left   = True

        plot1 = axs[0].contourf(lons, lats, var_regrid, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
        cbar = plt.colorbar(plot1, ax=axs[0], ticklocation="right", pad=0.08, orientation="horizontal",
                aspect=40, shrink=1) # cax=cax,
        cbar.ax.tick_params(labelsize=10,labelrotation=45)

        plot2 = axs[1].contourf(lons, lats, sen_regrid, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') # clevs_percentage,
        cbar  = plt.colorbar(plot2, ax=axs[1], ticklocation="right", pad=0.08, orientation="horizontal",
                aspect=40, shrink=1) # cax=cax,
        cbar.ax.tick_params(labelsize=10,labelrotation=45)

        plt.savefig('./plots/spatial_map_' +message + '_' + var_name + '.png',dpi=300)

if __name__ == "__main__":
    region = "Aus" #"SE Aus" #"CORDEX" #"SE Aus"

    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-25]
        loc_lon    = [135,155]
    elif region == "CORDEX":
        loc_lat    = [-52.36,3.87]
        loc_lon    = [89.25,180]

    if 0:
        LAI_MODIS_path = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/MCD15A3H_c61_bigWRFroi_LAI_for_WRF_5000m_201911_202002.nc"
        plot_LAI_MODIS(LAI_MODIS_path)

    if 0:
        fire_path = '/g/data/w97/mm3972/data/MODIS/MODIS_fire/MCD64A1.061_500m_aid0001.nc'
        plot_fire_map(fire_path)

    if 0:
        fire_path = '/g/data/w97/mm3972/data/MODIS/MODIS_fire/MCD64A1.061_500m_aid0001.nc'
        LAI_MODIS_path = "/g/data/w97/mm3972/data/MODIS/MODIS_LAI/MCD15A3H_c61_bigWRFroi_LAI_for_WRF_5000m_201911_202002.nc"
        plot_LAI_fire_map(fire_path,LAI_MODIS_path)

    if 1:
        '''
        Difference plot yearly

        '''
        case_ctl       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
        case_sen       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"
        land_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/LIS_output/"
        land_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/LIS_output/"

        fire_path = '/g/data/w97/mm3972/data/MODIS/MODIS_fire/MCD64A1.061_500m_aid0001.nc'
        var_names  = [  "Tmax","Tmin",
                        #"FWsoil_tavg",
                        #"SM_top50cm",
                        # "VPDmax","VPDmin"
                        # "LAI_inst","Albedo_inst",
                        # "VPD",
                        # "GPP_tavg",#"NPP_tavg",
                        # "SurfTmax","SurfTmin",
                        # "Tmax","Tmin",
                        # "VegT_tavg","Tair_f_inst","AvgSurfT_tavg",
                        # # "VegTmax","VegTmin",
                        ]

        burn_message = "_unburnt"
        burn         = -1

        message    = "201920_summer"+burn_message
        time_s     = datetime(2019,12,1,0,0,0,0)
        time_e     = datetime(2020,3,1,0,0,0,0)

        plot_LIS_diff(fire_path, land_ctl_path, land_sen_path, var_names, time_s=time_s, time_e=time_e, lat_names="lat",
                      lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,
                      message=message, burn=burn)

        burn_message = "_burnt"
        burn         = 1

        message    = "201920_summer"+burn_message
        time_s     = datetime(2019,12,1,0,0,0,0)
        time_e     = datetime(2020,3,1,0,0,0,0)

        plot_LIS_diff(fire_path, land_ctl_path, land_sen_path, var_names, time_s=time_s, time_e=time_e, lat_names="lat",
                      lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,
                      message=message, burn=burn)

        # message    = "HW_20191220-20191223"+burn_message
        # time_s     = datetime(2019,12,20,0,0,0,0)
        # time_e     = datetime(2019,12,24,0,0,0,0)
        #
        # plot_LIS_diff(fire_path, land_ctl_path, land_sen_path, var_names, time_s=time_s, time_e=time_e, lat_names="lat",
        #               lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,
        #               message=message, burn=burn)
        #
        # message    = "HW_20191230-20200101"+burn_message
        # time_s     = datetime(2019,12,30,0,0,0,0)
        # time_e     = datetime(2020,1,2,0,0,0,0)
        #
        # plot_LIS_diff(fire_path, land_ctl_path, land_sen_path, var_names, time_s=time_s, time_e=time_e, lat_names="lat",
        #               lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,
        #               message=message, burn=burn)
        #
        # message    = "HW_20200131-20200203"+burn_message
        # time_s     = datetime(2020,1,31,0,0,0,0)
        # time_e     = datetime(2020,2,4,0,0,0,0)
        #
        # plot_LIS_diff(fire_path, land_ctl_path, land_sen_path, var_names, time_s=time_s, time_e=time_e, lat_names="lat",
        #               lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,
        #               message=message, burn=burn)
