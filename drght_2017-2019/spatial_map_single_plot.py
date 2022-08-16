#!/usr/bin/env python
"""
Produce the Australia map with the label of EucFACE site - Figure 1
"""

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"


#!/usr/bin/python

import sys
import cartopy
import numpy as np
from netCDF4 import Dataset,num2date
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
import scipy.ndimage as ndimage
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def spatial_map_single_plot(file_path, var_name, time_s, time_e, lat_name, lon_name,
                            loc_lat=None, loc_lon=None, wrf_path=None,message=None):

    '''
    plot a single spatial map
    '''

    time, Var  = read_var(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
    print(time)
    var        = spital_var(time,Var,time_s,time_e)

    if 'LIS' in file_path:
        wrf        = Dataset(wrf_path,  mode='r')
        lons       = wrf.variables['XLONG'][0,:,:]
        lats       = wrf.variables['XLAT'][0,:,:]
        var        =  var/1000.
    else:
        time, lats = read_var(file_path, lat_name, loc_lat, loc_lon, lat_name, lon_name)
        time, lons = read_var(file_path, lon_name, loc_lat, loc_lon, lat_name, lon_name)


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
    # clevs = np.linspace( 0.,10., num=11)
    clevs = np.linspace( 0.,600., num=31)
    cmap  = plt.cm.GnBu_r # BrBG

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
    gl.xlines        = True

    if loc_lat == None:
        gl.xlocator     = mticker.FixedLocator([135,140,145,150,155])
        gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25])
    else:
        gl.xlocator = mticker.FixedLocator(loc_lon)
        gl.ylocator = mticker.FixedLocator(loc_lat)

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':10, 'color':'black'}
    gl.ylabel_style = {'size':10, 'color':'black'}

    plt.contourf(lons, lats, var,clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #

    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    cb.ax.tick_params(labelsize=10)
    plt.title(var_name, size=16)

    if message == None:
        message = var_name
    else:
        message = message + "_" + var_name

    plt.savefig('./plots/weather/spatial_map_'+message+'.png',dpi=300)

def spatial_map_single_plot_diff(file_paths, var_names, time_s=None, time_e=None, lat_names="lat",
                                 lon_names="lon",loc_lat=None, loc_lon=None, wrf_path=None,message=None):

    '''
    plot a single spatial map
    '''

    # WRF-CABLE
    if 'LIS_HIST_' in file_paths[0]:
        var_file     = Dataset(file_paths[0], mode='r')
        var1         = var_file.variables[var_names[0]][:,:]
        lats         = var_file.variables[lat_names[0]][:,:]
        lons         = var_file.variables[lon_names[0]][:,:]
    else:
        time1, Var1  = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, lats  = read_var(file_paths[0], lat_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        time1, lons  = read_var(file_paths[0], lon_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
        var1         = spital_var(time1,Var1,time_s,time_e)

    if 'LIS' in file_paths[0] and var_names[0] in ['WaterTableD_tavg','WatTable']:
        var1     = var1/1000.
    if var_names[0] in ['ESoil_tavg','Evap_tavg','TVeg_tavg']:
        var1     = var1*3600

    # read lat and lon outs
    wrf          = Dataset(wrf_path,  mode='r')
    lons_out     = wrf.variables['XLONG'][0,:,:]
    lats_out     = wrf.variables['XLAT'][0,:,:]

    # offline-CABLE

    if 'cable_out' in file_paths[1]:
        # offline sim
        time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        var2        = spital_var(time2,Var2,time_s,time_e)
    elif '-' in file_paths[1]:
        # lis-cable hist
        time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        var2        = Var2[0,:,:]
    elif 'LIS_HIST_' in file_paths[1]:
        var_file2   = Dataset(file_paths[1], mode='r')
        var2        = var_file2.variables[var_names[1]][:,:]
    elif 'LIS' in file_paths[0]:
        # lis restart
        time2, Var2 = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        var2        = Var2[-1]

    if 'LIS' in file_paths[1] and var_names[1] in ['WaterTableD_tavg','WatTable']:
        var2     = var2/1000.
    if var_names[1] in ['ESoil_tavg','Evap_tavg','TVeg_tavg']:
        var2     = var2*3600

    if 'cable_out' in file_paths[1] :
        time, lats_in= read_var(file_paths[1], lat_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
        time, lons_in= read_var(file_paths[1], lon_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])

        if var_names[1] in ['SoilMoist','SoilMoist_inst','ssat','sucs']:
            nlat   = len(lats_out[:,0])
            nlon   = len(lats_out[0,:])
            var2_regrid  = np.zeros((6,nlat,nlon))
            for j in np.arange(6):
                var2_regrid[j,:,:]  = regrid_data(lats_in, lons_in, lats_out, lons_out, var2[j,:,:])
        else:
            var2_regrid  = regrid_data(lats_in, lons_in, lats_out, lons_out, var2)

        if var_names[1] in ['ssat','sucs']:
            nlat   = len(lats_out[:,0])
            nlon   = len(lats_out[0,:])
            var    = np.zeros((6,nlat,nlon))
            for j in np.arange(6):
                var[j,:,:] = var1-var2_regrid[j,:,:]
        else:
            var    = var1-var2_regrid

    elif 'LIS' in file_paths[1]:
        var          = var1-var2

    if var_names[0] in ['WaterTableD_tavg','WatTable']:
        clevs = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]
    elif var_names[0] in ['GWwb_tavg','GWMoist']:
        clevs = [-0.05,-0.04,-0.03,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.03,0.04,0.05]
    elif var_names[0] in ['SoilMoist_inst','SoilMoist']:
        clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]
    elif var_names[0] in ['Sucs_inst','sucs_inst']:
        clevs = [-100,-80,-60,-40,-20,-10,10,20,40,60,80,100]
    elif var_names[0] in ['ESoil_tavg','Evap_tavg','TVeg_tavg']:
        clevs = [-1,-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8,1]
    else:
        # clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]
        clevs = [-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5]

    print("len(np.shape(var))",len(np.shape(var)))

    if len(np.shape(var)) >=3:

        for j in np.arange(6):

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
            cmap  = plt.cm.seismic

            #hot_r # BrBG

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
            gl.xlines        = True

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
            plt.contourf(lons, lats, var[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #

            cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
            cb.ax.tick_params(labelsize=10)
            plt.title(var_names[0], size=16)

            if j == 0:
                if message == None:
                    message = var_names[0]
                else:
                    message = message + "_" + var_names[0]

            plt.savefig('./plots/WTD_sudden_change/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
            cb = None
            gl = None
            ax = None
            fig= None

    elif len(np.shape(var)) ==2:
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
        cmap  = plt.cm.seismic

        #hot_r # BrBG

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
        gl.xlines        = True

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

        plt.contourf(lons, lats, var, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #
        print(var)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)
        plt.title(var_names[0], size=16)

        if message == None:
            message = var_names[0]
        else:
            message = message + "_" + var_names[0]

        plt.savefig('./plots/WTD_sudden_change/spatial_map_'+message+'.png',dpi=300)

def spatial_map_total_soil_water_diff(file_paths, lis_path, loc_lat=None, loc_lon=None, wrf_path=None,message=None):

    '''
    calculate the total soil water change to estimate water table depth changes
    '''

    # read WRF-CABLE
    var_file   = Dataset(file_paths[0], mode='r')
    wb1        = var_file.variables['SoilMoist_inst'][0,:,:,:]
    gwwb1      = var_file.variables['GWwb_tavg'][0,:,:]

    # read LIS_RST
    rst_file   = Dataset(file_paths[1], mode='r')
    wb2        = rst_file.variables['SoilMoist_inst'][-1,:,:,:]
    gwwb2      = rst_file.variables['GWwb_tavg'][-1,:,:]

    # calculate soil water storage
    lis_file   = Dataset(lis_path, mode='r')
    dtb        = lis_file.variables['DTB']

    sws1       = wb1[0,:,:]*0.005 + wb1[1,:,:]*0.075 + wb1[2,:,:]*0.154 + wb1[3,:,:]*0.409 + \
                 wb1[4,:,:]*1.085 + wb1[5,:,:]*2.872 + gwwb1*dtb

    sws2       = wb2[0,:,:]*0.005 + wb2[1,:,:]*0.075 + wb2[2,:,:]*0.154 + wb2[3,:,:]*0.409 + \
                 wb2[4,:,:]*1.085 + wb2[5,:,:]*2.872 + gwwb2*dtb

    sws_diff   = sws1-sws2

    # read lat and lon outs
    wrf        = Dataset(wrf_path,  mode='r')
    lons       = wrf.variables['XLONG'][0,:,:]
    lats       = wrf.variables['XLAT'][0,:,:]

    clevs      = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]

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
    gl.xlines        = True

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

    plt.contourf(lons, lats, sws_diff, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #

    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    cb.ax.tick_params(labelsize=10)
    plt.title(message, size=16)

    if message == None:
        message = 'Soil_watar_storage_diff'
    else:
        message = message + "_" + 'Soil_watar_storage_diff'

    plt.savefig('./plots/weather/spatial_map_'+message+'.png',dpi=300)

def spatial_map_land_info_diff(file_paths, loc_lat=None, loc_lon=None, wrf_path=None,message=None):

    '''
    compare land information
    '''

    var_names = ['SiltFrac_inst','silt'] #['ClayFrac_inst','clay']
    #['SandFrac_inst','sand'] #['sand_vec','SAND'] #['Sucs_inst','sucs'] # ['SoilSat_inst', 'ssat']

    if var_names[0] in ['Sucs_inst','sucs','sand_vec','SAND']:
        loop = 6
    else:
        loop = 1

    # read WRF-CABLE
    var_file   = Dataset(file_paths[0], mode='r')
    var1       = var_file.variables[var_names[0]][0,:,:] #

    # read off-CABLE
    off_file   = Dataset(file_paths[1], mode='r')
    var2       = off_file.variables[var_names[1]]#[:,:,:]
    print(file_paths[1])
    print(var2)
    print(off_file)
    lats_in    = off_file.variables['latitude'][:]
    lons_in    = off_file.variables['longitude'][:]

    # read lat and lon outs
    if wrf_path == None:
        lats       = var_file.variables['lat']
        lons       = var_file.variables['lon']
        lons_out, lats_out = np.meshgrid(lons, lats)
    else:
        wrf        = Dataset(wrf_path,  mode='r')
        lons_out   = wrf.variables['XLONG'][0,:,:]
        lats_out   = wrf.variables['XLAT'][0,:,:]

    nlat       = len(lons_out[:,0])
    nlon       = len(lons_out[0,:])


    var2_regrid = regrid_data(lats_in, lons_in, lats_out, lons_out, var2) #[j,:,:]
    if loop > 1:
        var_diff   = np.zeros((loop,nlat,nlon))
        for j in np.arange(loop):
            var_diff[j,:,:] = var1[j,:,:] - var2_regrid
    else:
        var_diff = var1 - var2_regrid

    if var_names[0] in ['SoilSat_inst', 'ssat']:
        clevs      = [-0.03,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.03]
    elif var_names[0] in ['SAND','sand_vec','SandFrac_inst','ClayFrac_inst','SiltFrac_inst']:
        clevs      = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]
    elif var_names[0] in ['Sucs_inst','sucs']:
        clevs      = [-100.,-80,-60,-40,-20,-10,10,20,40,60,80,100]

    # ================== Start Plotting =================
    for j in np.arange(loop):
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
        gl.xlines        = True

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

        if loop >1:
            plt.contourf(lons_out, lats_out, var_diff[j,:,:],  clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
        else:
            plt.contourf(lons_out, lats_out, var_diff,  clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)
        plt.title(message, size=16)
        if j == 0:
            if message == None:
                message = var_names[0]+'_diff'
            else:
                message = message +"_"+ var_names[0]+'_diff'
        if loop > 1:
            plt.savefig('./plots/weather/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
        else:
            plt.savefig('./plots/weather/spatial_map_'+message+'.png',dpi=300)

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


    #######################################################
    # Decks to run:
    #    plot a single map
    #######################################################
    # message    = "LIS-CABLE"
    # var_name   = 'WaterTableD_tavg'
    # lat_name   = "lat"
    # lon_name   = "lon"
    # file_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc"
    # wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"

    # message    = "off-CABLE"
    # var_name   = 'WatTable'
    # lat_name   = "latitude"
    # lon_name   = "longitude"
    # file_path  = "/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc"
    # wrf_path   = None

    # time_s     = datetime(2017,1,1,0,0,0,0)
    # time_e     = datetime(2017,1,1,23,59,0,0)
    # # # spatial_map_single_plot(file_path, var_name, time_s, time_e, lat_name, lon_name, loc_lat=loc_lat,loc_lon=loc_lon, wrf_path=wrf_path,message=message)

    # message    = "LIS-OFF"
    # file_paths = ["/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc",
    #               "/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc"]
    # wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"
    # var_names  = ['SoilMoist_inst','SoilMoist']
    # # var_names  = ['GWwb_tavg','GWMoist']
    # # var_names  = ['WaterTableD_tavg','WatTable']
    # lat_names  = ['lat','latitude']
    # lon_names  = ['lon','longitude']
    # spatial_map_single_plot_diff(file_paths, var_names, time_s, time_e, lat_names, lon_names,
    #                             loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)

    # time_s     = datetime(2017,1,1,0,0,0,0)
    # time_e     = datetime(2017,1,1,23,59,0,0)

    # message    = "LIS-LIS_rst"
    # file_paths = ["/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc",
    #               "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/LIS_output/LIS.CABLE.20170101110000.d01.nc"]
    # wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"

    # var_names  = ['SoilMoist_inst','SoilMoist_inst']
    # # var_names  = ['GWwb_tavg','GWwb_tavg']
    # # var_names  = ['WaterTableD_tavg','WaterTableD_tavg']
    # lat_names  = ['lat','lat']
    # lon_names  = ['lon','lon']
    # spatial_map_single_plot_diff(file_paths, var_names, time_s, time_e, lat_names, lon_names,
    #                             loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)


    if 0:
        message    = "LIS_depth_varying-LIS_rst_2016-12-31"
        file_paths = ['/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run/depth_varying/OUTPUT/SURFACEMODEL/LIS_HIST_201701011200.d01.nc',
                      '/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc']
        wrf_path   = "/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run_old/Jan_unfinish/wrfout_d01_2017-01-01_11:00:00"
        lat_names  = ['lat','latitude']
        lon_names  = ['lon','longitude']
        time_s     = datetime(2016,12,31,0,0,0,0)
        time_e     = datetime(2016,12,31,23,59,0,0)

        var_names  = ['WaterTableD_tavg','WatTable']
        spatial_map_single_plot_diff(file_paths, var_names, time_s=time_s, time_e = time_e, lat_names=lat_names, lon_names=lon_names,
                                    loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)
        var_names  = ['SoilMoist_inst','SoilMoist']
        spatial_map_single_plot_diff(file_paths, var_names, time_s=time_s, time_e = time_e, lat_names=lat_names, lon_names=lon_names,
                                    loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)

    if 1:
        message    = "LIS_depth_varying_ts1-LIS_rst_2016-12-31"
        file_paths = ['/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run/depth_varying/OUTPUT/SURFACEMODEL/LIS_HIST_201701011200.d01.nc',
                      '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/offline_rst_output/output_1719_drght/LIS.CABLE.20170101110000-day1.d01.nc',]
                    # '/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run/depth_varying/OUTPUT/SURFACEMODEL/LIS_HIST_201701021200.d01.nc',
                    #  '/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run/depth_varying/OUTPUT/SURFACEMODEL/LIS_HIST_201701011200.d01.nc']
                    #   '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght1719_bdy_data/LIS_output/LIS.CABLE.20170101110000.d01.nc']
        wrf_path   = "/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run_old/Jan_unfinish/wrfout_d01_2017-01-01_11:00:00"
        lat_names  = ['lat','lat']#'latitude']
        lon_names  = ['lon','lon']#'longitude']


        # var_names  = ['GWwb_tavg','GWwb_tavg']#'GWMoist']
        # var_names  = ['WaterTableD_tavg','WaterTableD_tavg']#'WatTable']
        # var_names  = ['SoilMoist_inst','SoilMoist_inst']#'SoilMoist']
        # var_names  = ['WaterTableD_tavg','WaterTableD_tavg']
        var_names = ['Evap_tavg','Evap_tavg']
        # var_names  = ['SoilTemp_inst','SoilTemp_inst']
        spatial_map_single_plot_diff(file_paths, var_names, lat_names=lat_names, lon_names=lon_names,
                                    loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)
        # var_names  = ['Tair_f_inst','Tair_f_inst']
        # var_names  = ['SoilMoist_inst','SoilMoist_inst']#'SoilMoist']
        var_names  = ['TVeg_tavg','TVeg_tavg']
        spatial_map_single_plot_diff(file_paths, var_names, lat_names=lat_names, lon_names=lon_names,
                                    loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)
        var_names  = ['ESoil_tavg','ESoil_tavg']
        spatial_map_single_plot_diff(file_paths, var_names, lat_names=lat_names, lon_names=lon_names,
                                    loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)
        var_names  = ['FWsoil_tavg','FWsoil_tavg']
        spatial_map_single_plot_diff(file_paths, var_names, lat_names=lat_names, lon_names=lon_names,
                                    loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)


    if 0:
        message    = "LIS-LIS_rst"
        file_paths = ['/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc',
                    '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/LIS_output/LIS.CABLE.20170101110000.d01.nc']
        lis_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/bdy_data/lis_input.d01.nc"
        wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"

        spatial_map_total_soil_water_diff(file_paths, lis_path, loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)


    if 0:
        message    = "LIS-LIS_rst"
        file_paths = ['/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc',
                      '/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc']
        wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"
        spatial_map_land_info_diff(file_paths, loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)

    if 0:
        message    = "LIS_new-off_gridinfo"
        # file_paths = ['/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc',
        file_paths = ['/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run/OUTPUT/LIS_HIST_201701011200.d01.nc',
                      '/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc']
        wrf_path   = "/scratch/w97/mm3972/model/NUWRF/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/coupled_run/wrfout_d01_2017-01-01_11:00:00"
        spatial_map_land_info_diff(file_paths, loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)

    if 0:
        message    = "gridinfo-LIS_input"
        file_paths = ['/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix_10km_6layer_uniform.nc',
                      '/g/data/w97/mm3972/scripts/process_data/make_LIS_landinfo/nc_file/Openlandmap_soilcomposition_CORDEX_180E.nc']
        spatial_map_land_info_diff(file_paths, loc_lat=loc_lat, loc_lon=loc_lon, message=message)
