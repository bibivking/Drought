#!/usr/bin/python

'''
Functions:
1. Analyze weather situation during heatwaves
2. Process ERAI, AWAP, offline CABLE and LIS-CABLE data
'''

from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from wrf import (getvar, interplevel, get_cartopy, cartopy_xlim,
                 cartopy_ylim, to_np, latlon_coords)
from common_utils import *
# from spatial_wrf_hgt_var import plot_spatial_map_hgt

def plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=None, loc_lon=None, lat_names=None, lon_names=None, message=None, diff=False):

    print("======== In plot_spital_map =========")
    print(var_names)
    if diff == False:
        # Open the NetCDF4 file (add a directory path if necessary) for reading:
        time1, Var1  = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        print(time1)
        print("The number of pixels are Nan: ", np.count_nonzero(np.isnan(Var1)))
        time1, lats1 = read_var(file_paths[0], lat_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, lons1 = read_var(file_paths[0], lon_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])

        if var_names[0] in ['tas','Tair','Tair_f_inst']:
            var1         = spital_var(time1,Var1,time_s,time_e)-273.15
            print(var1)
        elif var_names[0] in ['tp']:
            scale        = get_scale(var_names[0])
            var1         = spital_ERAI_tp(time1,Var1,time_s,time_e)*scale
        elif var_names[0] in ['Rainf','Rainf_tavg']:
            var1         = spital_var(time1,Var1,time_s,time_e)*24*60*60.
            print(var1)
        elif var_names[0] in ['Wind']:
            # !!!!!!!!! Note that !!!!!!!!!!!
            # Wind speeds is at 2 meter height in AWAP while 10 meters in WRF-CABLE
            # So here converting AWAP 2m wind speed to 10m wind speed by multipling 2
            var1         = spital_var(time1,Var1,time_s,time_e)*2.
            print(var1)
        else:
            scale        = get_scale(var_names[0])
            var1         = spital_var(time1,Var1,time_s,time_e)*scale

        if len(file_paths) > 1:
            time2, Var2  = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
            time2, lats2 = read_var(file_paths[1], lat_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
            time2, lons2 = read_var(file_paths[1], lon_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
            scale        = get_scale(var_names[1])
            var2         = spital_var(time2,Var2,time_s,time_e)*scale

        if len(file_paths) > 2:
            time3, Var3  = read_var(file_paths[2], var_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
            time3, lats3 = read_var(file_paths[2], lat_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
            time3, lons3 = read_var(file_paths[2], lon_names[2], loc_lat, loc_lon, lat_names[2], lon_names[2])
            scale        = get_scale(var_names[2])
            var3         = spital_var(time3,Var3,time_s,time_e)*scale

        if len(file_paths) > 3:
            time4, Var4  = read_var(file_paths[3], var_names[3], loc_lat, loc_lon, lat_names[3], lon_names[3])
            time4, lats4 = read_var(file_paths[3], lat_names[3], loc_lat, loc_lon, lat_names[3], lon_names[3])
            time4, lons4 = read_var(file_paths[3], lon_names[3], loc_lat, loc_lon, lat_names[3], lon_names[3])
            scale        = get_scale(var_names[3])
            var4         = spital_var(time4,Var4,time_s,time_e)*scale

    elif diff == True:
        time1, Var1  = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, lats1 = read_var(file_paths[0], lat_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, lons1 = read_var(file_paths[0], lon_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])

        time2, Var2  = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        time2, lats2 = read_var(file_paths[1], lat_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        time2, lons2 = read_var(file_paths[1], lon_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])

        print("np.shape(Var1) ", np.shape(Var1))
        print("np.shape(Var2) ", np.shape(Var2))

        if var_names[0] in ['tas','Tair','Tair_f_inst']:
            var1         = spital_var(time1,Var1,time_s,time_e)-273.15 # AWAP
            var2         = spital_var(time2,Var2,time_s,time_e)-273.15 # WRF
            # regrid_data(lat_in, lon_in, lat_out, lon_out, input_data)
            var1_regrid  = regrid_data(lats1, lons1, lats2, lons2, var1)
            var1         = var2 - var1_regrid
        elif var_names[0] in ['tp']:
            scale        = get_scale(var_names[0])
            var1         = spital_ERAI_tp(time1,Var1,time_s,time_e)*scale
            var2         = spital_ERAI_tp(time2,Var2,time_s,time_e)*scale
            # regrid_data(lat_in, lon_in, lat_out, lon_out, input_data)
            var1_regrid  = regrid_data(lats1, lons1, lats2, lons2, var1)
            var1         = var2 - var1_regrid
        elif var_names[0] in ['Rainf','Rainf_tavg']:
            var1         = spital_var(time1,Var1,time_s,time_e)*24*60*60.
            var2         = spital_var(time2,Var2,time_s,time_e)*24*60*60.
            # regrid_data(lat_in, lon_in, lat_out, lon_out, input_data)
            var1_regrid  = regrid_data(lats1, lons1, lats2, lons2, var1)
            var1         = var2 - var1_regrid
        else:
            scale        = get_scale(var_names[0])
            var1         = spital_var(time1,Var1,time_s,time_e)*scale
            var2         = spital_var(time2,Var2,time_s,time_e)*scale
            # regrid_data(lat_in, lon_in, lat_out, lon_out, input_data)
            var1_regrid  = regrid_data(lats1, lons1, lats2, lons2, var1)
            var1         = var2 - var1_regrid

    # ================== Start Plotting =================
    fig = plt.figure(figsize=(6,5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # start plotting
    if loc_lat == None:
        # ax.set_extent([140,154,-40,-28])
        ax.set_extent([135,155,-40,-25])
    else:
        ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

    ax.coastlines(resolution="50m",linewidth=1)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top   = False
    gl.ylabels_right = False
    gl.xlines        = True

    if loc_lat == None:
        # gl.xlocator = mticker.FixedLocator([140,145,150])
        # gl.ylocator = mticker.FixedLocator([-40,-35,-30])
        gl.xlocator     = mticker.FixedLocator([135,140,145,150,155])
        gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25])
    else:
        gl.xlocator = mticker.FixedLocator(loc_lon)
        gl.ylocator = mticker.FixedLocator(loc_lat)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}


    # plot Var1
    if var_names[0] in ['tas','Tair','Tair_f_inst']:
        clevs = np.arange( 15.,40.,1.) #np.linspace( 15.,45., num=31)
        cmap  = plt.cm.RdYlBu_r
    elif var_names[0] in ['Rainf','Rainf_tavg']:
        clevs = np.arange( 0.,22.,2.) #np.linspace( 15.,45., num=31)
        cmap  = plt.cm.Blues
    elif var_names[0] in ['Wind','Wind_f_inst']:
        clevs = np.arange( 0,10.,0.5) #np.linspace( 15.,45., num=31)
        cmap  = plt.cm.Blues
    elif var_names[0] in ['LWdown','LWdown_f_inst','SWdown','SWdown_f_inst']:
        clevs = np.arange( 80.,500.,20.) #np.linspace( 15.,45., num=31)
        cmap  = plt.cm.RdYlBu_r
    elif var_names[0] in ['Qair','Qair_f_inst']:
        # kg kg-1
        clevs = np.arange( 0.,0.02, 0.001) #np.linspace( 15.,45., num=31)
        cmap  = plt.cm.RdYlBu_r
    else:
        # clevs = np.linspace( 0.,120., num=13)
        clevs = np.linspace( 0.,5., num=11)
        cmap  = plt.cm.GnBu # BrBG
    # print(var1)
    # np.savetxt("./"+message,var1)
    plt.contourf(lons1, lats1, var1, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #,#bwr)#coolwarm)#cm.BrBG) # clevs,

    plt.title(var_names[0], size=16)
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    # cb.set_label(units,size=14,rotation=270,labelpad=15)
    cb.ax.tick_params(labelsize=10)

    # plot Var2
    if len(file_paths) > 1 and var_names[1] == 'ps':
        clevs = np.arange( -100.,100., 20.)
        cs = plt.contour(lons2, lats2, var2-1010., clevs, transform=ccrs.PlateCarree(),linewidths=0.8,colors="darkgray") #, ,cmap=plt.cm.hot_r)#bwr)#coolwarm)#cm.BrBG) # clevs,
        cl = plt.clabel(cs, inline=True, fmt="%4d",fontsize=6) #manual=True)

    # plot Var3, Var4
    if len(file_paths) > 3 and var_names[2] == 'uas' and var_names[3] == 'vas':
        qv = plt.quiver(lons1[::3,::3], lats1[::3,::3], var3[::3,::3], var4[::3,::3], scale=300, color='k')

    if message == None:
        message = var_names[0]
    else:
        message = message + "_" + var_names[0]
    if diff:
        message = message + "_diff"
    plt.savefig('./plots/weather/spatial_map_weather_analysis_'+message+'.png',dpi=300)

def plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=None, loc_lon=None, lat_names=None, lon_names=None, message=None):

    print("======== In plot_spital_map_multi =========")

    # ================== Plot setting ==================
    case_sum = np.shape(file_paths)[0]
    print(case_sum//2)
    fig, ax = plt.subplots( nrows=(case_sum // 2), ncols=2, figsize=[12,10],sharex=True, sharey=True, squeeze=True,
                            subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=0, hspace=0.2) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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

    for i, file_path in enumerate(file_paths):

        row = i // 2
        col = i % 2
        print("row & col ", row, col)
        # ================== Reading data =================
        time, Var  = read_var(file_path, var_names[i], loc_lat, loc_lon, lat_names[i], lon_names[i])
        time, lats = read_var(file_path, lat_names[i], loc_lat, loc_lon, lat_names[i], lon_names[i])
        time, lons = read_var(file_path, lon_names[i], loc_lat, loc_lon, lat_names[i], lon_names[i])

        if var_names[i] in ['2t','tas','Tair','Tair_f_inst']:
            if i == 0:
                clevs    = np.linspace( 15.,45., num=31)
                cmap     = plt.cm.RdYlBu_r
            else:
                clevs    = [-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5] # np.linspace( -5., 5., num=11)
                cmap     = plt.cm.RdYlBu_r
            var      = spital_var(time,Var,time_s,time_e)-273.15
        elif var_names[i] in ['Rainf','Rainf_tavg','tp']:
            if i == 0:
                clevs    = np.linspace( 0., 300., num=16)
                cmap     = plt.cm.BrBG
            else:
                clevs    = [-180,-160,-140,-120,-100,-80,-60,-40,-20,-10,10,20,40,60,80,100,120,140,160,180]
                cmap     = plt.cm.BrBG #RdYlBu_r
            var      = spital_var(time,Var,time_s,time_e)*24*60*60.*30
        elif var_names[i] in ['LWdown','LWdown_f_inst','SWdown','SWdown_f_inst']:
            if i == 0:
                clevs = np.arange( 80.,520.,20.) #np.linspace( 15.,45., num=31)
            else:
                clevs = np.arange( -90.,100.,10.)
            cmap  = plt.cm.BrBG_r
            scale = get_scale(var_names[i])
            var   = spital_var(time,Var,time_s,time_e)*scale
        elif var_names[i] in ['Wind','Wind_f_inst']:
            if i == 0:
                clevs = np.arange( 0,10.5,0.5) #np.linspace( 15.,45., num=31)
                var   = spital_var(time,Var,time_s,time_e)*2
            else:
                clevs = np.arange( -5,5.5,0.5)
                var   = spital_var(time,Var,time_s,time_e)
            cmap  = plt.cm.BrBG
        elif var_names[i] in ['Qair','Qair_f_inst']:
            # kg kg-1
            if i == 0:
                clevs = np.arange( 0.,0.02, 0.001) #np.linspace( 15.,45., num=31)
            else:
                clevs = np.arange( -0.006,0.007, 0.001)
            cmap  = plt.cm.BrBG
            scale = get_scale(var_names[i])
            var   = spital_var(time,Var,time_s,time_e)*scale
        else:
            if i == 0:
                clevs = np.linspace( 0.,5., num=11)
            else:
                clevs = np.linspace( -5.,5., num=11)
            cmap  = plt.cm.GnBu # BrBG
            scale = get_scale(var_names[i])
            var   = spital_var(time,Var,time_s,time_e)*scale

        if i == 0:
            # save AWAP dataset
            lat_AWAP   = lats
            lon_AWAP   = lons
            var_AWAP   = var
        elif i == 1:
            wrf = Dataset(wrf_path,  mode='r')
            lons_out = wrf.variables['XLONG'][0,:,:]
            lats_out = wrf.variables['XLAT'][0,:,:]
            AWAP_regrid  = regrid_data(lat_AWAP, lon_AWAP, lats_out, lons_out, var_AWAP)
            var          = var - AWAP_regrid
        else:
            var          = var - AWAP_regrid

        # =============== setting plots ===============
        if loc_lat == None:
            # ax.set_extent([140,154,-40,-28])
            ax[row,col].set_extent([135,155,-40,-25])
        else:
            ax[row,col].set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

        ax[row,col].coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl               = ax[row,col].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')

        if loc_lat == None:
            gl.xlocator = mticker.FixedLocator([135,140,145,150,155])
            gl.ylocator = mticker.FixedLocator([-40,-35,-30,-25])
        else:
            gl.xlocator = mticker.FixedLocator([135,140,144.2,150,155])
            gl.ylocator = mticker.FixedLocator([-40,-35,-31.8,-25])

        gl.xformatter   = LONGITUDE_FORMATTER
        gl.yformatter   = LATITUDE_FORMATTER
        gl.xlabel_style = {'size':10, 'color':'black'}
        gl.ylabel_style = {'size':10, 'color':'black'}
        gl.xlabels_bottom= True
        gl.xlabels_top   = False
        gl.ylabels_left  = True
        gl.ylabels_right = False
        gl.xlines        = True
        gl.ylines        = True
        plot1 = ax[row, col].contourf(lons, lats, var, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #,#bwr)#coolwarm)#cm.BrBG) # clevs,
        cb    = plt.colorbar(plot1, ax=ax[row, col], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)
        ax[row, col].set_title(names[i], size=12)
    # cb.set_label(units,size=14,rotation=270,labelpad=15)

    if message == None:
        message = var_names[0]
    else:
        message = message + "_" + var_names[0]

    plt.savefig('./plots/weather/spatial_map_weather_analysis_'+message+'_multi.png',dpi=300)


if __name__ == "__main__":


    # ======================= Option =======================
    region = "SE Aus" #"SE Aus" #"CORDEX" #"SE Aus"

    # ====================== Pre-load =======================
    ERAI_path    = '/g/data/ub4/erai/netcdf/3hr/atmos/oper_fc_sfc/v01'
    ERAI_T_file  = ERAI_path + '/tas/tas_3hrs_ERAI_historical_fc-sfc_20170101_20170131.nc' # air temperature
    ERAI_P_file  = ERAI_path + '/ps/ps_3hrs_ERAI_historical_fc-sfc_20170101_20170131.nc'   # surface pressure
    ERAI_U_file  = ERAI_path + '/uas/uas_3hrs_ERAI_historical_fc-sfc_20170101_20170131.nc' # 10 m wind speed
    ERAI_V_file  = ERAI_path + '/vas/vas_3hrs_ERAI_historical_fc-sfc_20170101_20170131.nc' # 10 m wind speed
    ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20170101_20170131.nc' # Total rainfall

    ERA5_path    = "/g/data/rt52/era5/single-levels/reanalysis"
    ERA5_T_file  = ERA5_path + '/2t/2017/2t_era5_oper_sfc_20170101-20170131.nc' # air temperature
    ERA5_P_file  = ERA5_path + '/sp/2017/sp_era5_oper_sfc_20170101-20170131.nc'   # surface pressure
    ERA5_U_file  = ERA5_path + '/10u/2017/10u_era5_oper_sfc_20170101-20170131.nc' # 10 m wind speed
    ERA5_V_file  = ERA5_path + '/10v/2017/10v_era5_oper_sfc_20170101-20170131.nc' # 10 m wind speed
    ERA5_R_file  = ERA5_path + '/tp/2017/tp_era5_oper_sfc_20170101-20170131.nc' # Total rainfall

    AWAP_path    = '/g/data/w97/W35_GDATA_MOVED/Shared_data/AWAP_3h_v1'
    AWAP_T_file  = AWAP_path + '/Tair/AWAP.Tair.3hr.2017.nc'     # air temperature
    AWAP_R_file  = AWAP_path + '/Rainf/AWAP.Rainf.3hr.2017.nc'   # Daily rainfall
    AWAP_LW_file = AWAP_path + '/LWdown/AWAP.LWdown.3hr.2017.nc'   # Downward Longwave Radiation
    AWAP_SW_file = AWAP_path + '/SWdown/AWAP.SWdown.3hr.2017.nc'  # Downward Shortwave Radiation
    AWAP_W_file  = AWAP_path + '/Wind/AWAP.Wind.3hr.2017.nc'     # Near surface wind speed
    AWAP_Q_file  = AWAP_path + '/Qair/AWAP.Qair.3hr.2017.nc'    # Near surface specific humidity

    cpl_land_path  = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2013_15Oct/gw_jan2019/LIS_output'
    cpl_land_file  = cpl_land_path + '/LIS.CABLE.201901-201901.d01.nc'  # land output of wrf-cable run

    cpl_atmo_file     = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/hw2009_15Oct/ensemble_avg'
    cpl_atmo_file_gw  = cpl_atmo_file + '/wrfout_20090122-20090213_gw'  # atmo output of wrf-cable run
    cpl_atmo_file_fd  = cpl_atmo_file + '/wrfout_20090122-20090213_fd'  # atmo output of wrf-cable run

    # datetime(year, month, day, hour, minute, second, microsecond)
    # 2017 heatwave 9-12 Jan
    hw_s    = [ datetime(2017,1,9,0,0,0,0),
                ]
    hw_e    = [ datetime(2017,1,10,0,0,0,0),
                ]
                # datetime(2014,1,18,0,0,0,0),
                # datetime(2017,2,12,0,0,0,0),

    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-25]
        loc_lon    = [135,155]
    elif region == "CORDEX":
        loc_lat    = [-52.36,3.87]
        loc_lon    = [89.25,180]


    # # =================== Opreation ====================

    # # #################################
    # # Plot AWAP
    # # #################################


    # # Plot AWAP Tair
    # file_paths  = [AWAP_T_file]

    # var_names   = ['Tair']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # # for i in np.arange(0,30):
    # time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    # message = "AWAP_2017-01"#+str(i+2)
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)

    # # Plot AWAP rainfall
    # file_paths  = [AWAP_R_file]

    # var_names   = ['Rainf']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # # for i in np.arange(0,30):
    # time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    # message = "AWAP_2017-01"#+str(i+2)
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)

    # '''
    # Plot AWAP Near surface specific humidity
    # '''

    # file_paths  = [AWAP_Q_file]

    # var_names   = ['Qair']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # # for i in np.arange(0,30):
    # time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    # message = "AWAP_2017-01"#+str(i+2)
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)


    # '''
    # Plot AWAP Wind
    # '''

    # file_paths  = [AWAP_W_file]

    # var_names   = ['Wind']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # # for i in np.arange(0,30):
    # time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    # message = "AWAP_2017-01"#+str(i+2)
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)

    # '''
    # Plot AWAP LWdown
    # '''

    # file_paths  = [AWAP_LW_file]

    # var_names   = ['LWdown']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # # for i in np.arange(0,30):
    # time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    # message = "AWAP_2017-01"#+str(i+2)
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)

    # # Plot AWAP SWdown
    # file_paths  = [AWAP_SW_file]

    # var_names   = ['SWdown']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]

    # # for i in np.arange(0,30):
    # time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    # message = "AWAP_2017-01"#+str(i+2)
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)

    #
    # # #################################
    # # Plot WRF-CABLE land var
    # # #################################
    # case_names     = [
    #                   "drght_2017_2019_bl_pbl5_mp6_sf_sfclay1",
    #                   "drght_2017_2019_bl_pbl7_mp8_sf_sfclay1" ]
    #                 #   "drght_2017_2019_bl_pbl5_mp8_sf_sfclay1",
    #                 # "drght_2017_2019",
    #                 #   "drght_2017_2019_bl_pbl1_mp6_sf_sfclay1",
    #     # "drght_2017_2019",
    # for case_name in case_names:
    #     # case_name      = "drght_2017_2019"
    #     cpl_land_path  = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/'+case_name+'/LIS_output'
    #     cpl_land_file  = cpl_land_path + '/LIS.CABLE.201701-201701.d01.nc'  # land output of wrf-cable run
    #
    #     file_paths  = [cpl_land_file]
    #
    #     # var_names   = ['Rainf_tavg'] # ,'Qair_f_inst','Rainf_tavg','Tair_f_inst'
    #     lat_names   = ["lat"]
    #     lon_names   = ["lon"]
    #
    #     # for i in np.arange(0,30):
    #     time_s = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    #     time_e = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    #     print(time_s)
    #     message     = case_name+"_2017-01"
    #     var_names   = ['Rainf_tavg']
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)
    #     var_names   = ['Tair_f_inst']
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)
    #
    #     var_names   = ['Qair_f_inst']
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)
    #
    #     var_names   = ['SWdown_f_inst']
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)
    #
    #     var_names   = ['LWdown_f_inst']
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)
    #
    #     var_names   = ['Wind_f_inst']
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)


    # # #################################
    # # Plot WRF-CABLE - AWAP
    # # #################################
    # message    = "WRF-AWAP_2017-01"
    # case_names = ["drght_2017_2019",
    #               "drght_2017_2019_bl_pbl1_mp6_sf_sfclay1",
    #               "drght_2017_2019_bl_pbl5_mp6_sf_sfclay1",
    #               "drght_2017_2019",
    #               "drght_2017_2019_bl_pbl7_mp8_sf_sfclay1"]
    #               # "drght_2017_2019_bl_pbl5_mp8_sf_sfclay1",
    # names = [ "AWAP",
    #           "bl_pbl2_mp4_sf_sfclay2",
    #           "bl_pbl1_mp6_sf_sfclay1",
    #           "bl_pbl5_mp6_sf_sfclay1",
    #           "drght_2017_2019",
    #           "bl_pbl7_mp8_sf_sfclay1"]
    # wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[0]+"/WRF_output/wrfout_d01_2017-01-01_11:00:00"
    # time_s     = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    # time_e     = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))
    #
    # # ============= Tair =============
    # file_paths  = [ AWAP_T_file]
    # var_names   = ['Tair']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('Tair_f_inst')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    #     print(file_paths)
    #
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)
    #
    # # ============= Rainf =============
    # file_paths  = [ AWAP_R_file]
    # var_names   = ['Rainf']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('Rainf_tavg')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    #
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)
    #
    # # ============= Wind =============
    # file_paths  = [ AWAP_W_file]
    # var_names   = ['Wind']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('Wind_f_inst')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    #
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)
    #
    #
    # # ============= LWdown =============
    # file_paths  = [ AWAP_LW_file]
    # var_names   = ['LWdown']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('LWdown_f_inst')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    #
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)
    #
    # # ============= SWdown =============
    # file_paths  = [ AWAP_SW_file]
    # var_names   = ['SWdown']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('SWdown_f_inst')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    #
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)
    #
    # # ============= Qair =============
    # file_paths  = [ AWAP_Q_file]
    # var_names   = ['Qair']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('Qair_f_inst')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    #
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)

    # #################################
    # Plot WRF-CABLE - ERA 5
    # #################################
    message    = "WRF-ERA5_2017-01"
    case_names = ["drght_2017_2019",
                  "drght_2017_2019_bl_pbl1_mp6_sf_sfclay1",
                  "drght_2017_2019_bl_pbl5_mp6_sf_sfclay1",
                  "drght_2017_2019",
                  "drght_2017_2019_bl_pbl7_mp8_sf_sfclay1"]
                  # "drght_2017_2019_bl_pbl5_mp8_sf_sfclay1",
    names = [ "ERA5",
              "bl_pbl2_mp4_sf_sfclay2",
              "bl_pbl1_mp6_sf_sfclay1",
              "bl_pbl5_mp6_sf_sfclay1",
              "drght_2017_2019",
              "bl_pbl7_mp8_sf_sfclay1"]
    wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_names[0]+"/WRF_output/wrfout_d01_2017-01-01_11:00:00"
    time_s     = datetime(2017,1,2,0,0,0,0) #+ timedelta(days=int(i))
    time_e     = datetime(2017,1,31,23,59,0,0) #+ timedelta(days=int(i))

    # ============= Tair =============
    file_paths  = [ ERA5_T_file]
    var_names   = ['2t']
    lat_names   = ['latitude']
    lon_names   = ['longitude']
    for case_name in case_names:
        file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
        var_names.append('Tair_f_inst')
        lat_names.append('lat')
        lon_names.append('lon')
    plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
                          loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)

    # # ============= Rain =============
    # !!!!!!!!!! NCI didn't download rainfall data
    # file_paths  = [ ERA5_R_file]
    # var_names   = ['tp']
    # lat_names   = ['lat']
    # lon_names   = ['lon']
    # for case_name in case_names:
    #     file_paths.append("/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/LIS.CABLE.201701-201701.d01.nc")
    #     var_names.append('Rainf_tavg')
    #     lat_names.append('lat')
    #     lon_names.append('lon')
    # plot_spital_map_multi(wrf_path, names, file_paths, var_names, time_s, time_e, loc_lat=loc_lat,
    #                       loc_lon=loc_lon, lat_names=lat_names,lon_names=lon_names,message=message)

    # # #################################
    # # Plot WRF-CABLE atmo var
    # # #################################
    #
    # file_paths  = [cpl_atmo_file_fd, cpl_atmo_file_gw]
    #
    # var_name    = "temp"
    # var_unit    = "degC"
    # height      = 1000
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]
    #
    # time_s = datetime(2009,1,22,0,0,0,0)
    # time_e = datetime(2009,2,14,0,0,0,0)
    #
    # for i in np.arange(0,23):
    #
    #     time_s = datetime(2009,1,22,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2009,1,22,23,59,0,0) + timedelta(days=int(i))
    #     message = "Couple_GW-FD"+str(time_s)
    #     plot_spatial_map_hgt(file_paths, var_name, var_unit, height, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lat, message=message)


    # ##################################
    # # Plot ERAI air temp + 10 m wind speed + pressure
    # ##################################

    # var_names   = ['tas','ps','uas','vas']
    # lat_names   = ["lat","lat","lat","lat"]
    # lon_names   = ["lon","lon","lon","lon"]

    # ERAI_period  = "20090101_20090131"
    # ERAI_T_file  = ERAI_path + '/tas/tas_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # air temperature
    # ERAI_P_file  = ERAI_path + '/ps/ps_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc'   # surface pressure
    # ERAI_U_file  = ERAI_path + '/uas/uas_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # 10 m wind speed
    # ERAI_V_file  = ERAI_path + '/vas/vas_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # 10 m wind speed
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # Total rainfall
    # file_paths   = [ERAI_T_file, ERAI_P_file, ERAI_U_file, ERAI_V_file]

    # for i in np.arange(0,31):
    #     time_s = datetime(2009,1,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2009,1,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "ERAI_2009-01-"+str(i+1)
    #     print(message)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names,message=message)

    # ERAI_period  = "20090201_20090228"
    # ERAI_T_file  = ERAI_path + '/tas/tas_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # air temperature
    # ERAI_P_file  = ERAI_path + '/ps/ps_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc'   # surface pressure
    # ERAI_U_file  = ERAI_path + '/uas/uas_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # 10 m wind speed
    # ERAI_V_file  = ERAI_path + '/vas/vas_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # 10 m wind speed
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_'+ERAI_period+'.nc' # Total rainfall
    # file_paths   = [ERAI_T_file, ERAI_P_file, ERAI_U_file, ERAI_V_file]

    # for i in np.arange(0,28):
    #     time_s = datetime(2009,2,1,0,0,0,0) + timedelta(days=int(i))
    #     time_e = datetime(2009,2,2,0,0,0,0) + timedelta(days=int(i))
    #     message = "ERAI_2009-02-"+str(i+1)
    #     print(message)
    #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                     lon_names=lon_names, message=message)

    #################################
    #    Plot ERAI rainfall         #
    #################################

    ### Month before HW ###
    # # 2009
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20090101_20090131.nc' # Total rainfall
    # time_s = datetime(2009,1,1,0,0,0,0)
    # time_e = datetime(2009,2,1,0,0,0,0)
    # message = "ERAI_2009-1_Rainfall"

    # # 2011
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20110101_20110131.nc' # Total rainfall
    # time_s = datetime(2011,1,1,0,0,0,0)
    # time_e = datetime(2011,2,1,0,0,0,0)
    # message = "ERAI_2011-1_Rainfall"

    # # 2013
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20121201_20121231.nc' # Total rainfall
    # time_s = datetime(2012,12,1,0,0,0,0)
    # time_e = datetime(2013,1,1,0,0,0,0)
    # message = "ERAI_2012-12_Rainfall"

    # # 2019
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20181201_20181231.nc' # Total rainfall
    # time_s = datetime(2018,12,1,0,0,0,0)
    # time_e = datetime(2019,1,1,0,0,0,0)
    # message = "ERAI_2018-12_Rainfall"

    # hw_e    = [ datetime(2009,2,8,0,0,0,0),
    #             datetime(2011,2,6,0,0,0,0),
    #             datetime(2013,1,9,0,0,0,0),
    #             datetime(2019,1,26,0,0,0,0)
    #             ]

    # ### HW's Month ###
    # # # 2009
    # # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20090201_20090228.nc' # Total rainfall
    # # time_s = datetime(2009,2,1,0,0,0,0)
    # # time_e = hw_e[0]
    # # message = "ERAI_2009hw_Rainfall"
    #
    # # # 2011
    # # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20110201_20110228.nc' # Total rainfall
    # # time_s = datetime(2011,2,1,0,0,0,0)
    # # time_e = hw_e[1]
    # # message = "ERAI_2011hw_Rainfall"
    #
    # # # 2013
    # # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20130101_20130131.nc' # Total rainfall
    # # time_s = datetime(2013,1,1,0,0,0,0)
    # # time_e = hw_e[2]
    # # message = "ERAI_2013hw_Rainfall"
    #
    # # 2019
    # ERAI_R_file  = ERAI_path + '/tp/tp_3hrs_ERAI_historical_fc-sfc_20190101_20190131.nc' # Total rainfall
    # time_s = datetime(2019,1,1,0,0,0,0)
    # time_e = datetime(2019,1,12,0,0,0,0)#hw_e[3]
    # message = "ERAI_20190101-12_Rainfall"
    #
    # file_paths  = [ERAI_R_file]
    #
    # var_names   = ['tp']
    # lat_names   = ["lat"]
    # lon_names   = ["lon"]
    #
    # # for i in np.arange(0,30):
    # #     time_s = datetime(2019,1,1,0,0,0,0) + timedelta(days=int(i))
    # #     time_e = datetime(2019,1,2,0,0,0,0) + timedelta(days=int(i))
    # #     message = "ERAI_2019-01-"+str(i+1)
    # #     plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    # #                     lon_names=lon_names,message=message)
    #
    # plot_spital_map(file_paths, var_names, time_s, time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_names=lat_names,
    #                 lon_names=lon_names,message=message)
