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

def spatial_map_single_plot_diff(file_paths, var_names, time_s, time_e, lat_names, lon_names, 
                                loc_lat=None, loc_lon=None, wrf_path=None,message=None):

    ''' 
    plot a single spatial map
    '''
    
    # WRF-CABLE
    time1, Var1  = read_var(file_paths[0], var_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, lats  = read_var(file_paths[0], lat_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, lons  = read_var(file_paths[0], lon_names[0], loc_lat, loc_lon, lat_names[0], lon_names[0])
    var1         = spital_var(time1,Var1,time_s,time_e)    
    if 'LIS' in file_paths[0] and var_names[0] in ['WaterTableD_tavg','WatTable']:
        var1     = var1/1000.
        
    # read lat and lon outs
    wrf          = Dataset(wrf_path,  mode='r')
    lons_out     = wrf.variables['XLONG'][0,:,:]
    lats_out     = wrf.variables['XLAT'][0,:,:]

    # offline-CABLE
    time2, Var2  = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])

    if 'cable_out' in file_paths[1]:
        var2     = spital_var(time2,Var2,time_s,time_e)
    elif 'LIS' in file_paths[1] :
        var2     = Var2[-1]
        
    time, lats_in= read_var(file_paths[1], lat_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
    time, lons_in= read_var(file_paths[1], lon_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
    if 'LIS' in file_paths[1] and var_names[1] in ['WaterTableD_tavg','WatTable']:
        var2     = var2/1000.
        
    print("shape lats_in=", np.shape(lats_in))
    print("shape lons_in=", np.shape(lons_in))
    print("shape lats_out=",np.shape(lats_out))
    print("shape lons_out=",np.shape(lons_out))
    print("shape var2=",    np.shape(var2))
    
    if 'cable_out' in file_paths[1] :
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
    else:
        clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]       
    
    if var_names[1] in ['SoilMoist','SoilMoist_inst','ssat','sucs']:
    
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

            plt.savefig('./plots/weather/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
            cb = None
            gl = None
            ax = None
            fig= None
    else:
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

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)
        plt.title(var_names[0], size=16)

        if message == None:
            message = var_names[0]
        else:
            message = message + "_" + var_names[0]

        plt.savefig('./plots/weather/spatial_map_'+message+'.png',dpi=300)

def spatial_map_total_soil_water_diff(file_paths, var_names, time_s, time_e, lat_names, lon_names, 
                                       loc_lat=None, loc_lon=None, wrf_path=None,message=None):

    ''' 
    calculate the total soil water change to estimate water table depth changes
    '''
    
    # WRF-CABLE
    time1, Wb1   = read_var(file_paths[0], 'SoilMoist_inst', loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, GWwb1 = read_var(file_paths[0], 'GWMoist_inst', loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, lats  = read_var(file_paths[0], 'lat', loc_lat, loc_lon, lat_names[0], lon_names[0])
    time1, lons  = read_var(file_paths[0], 'lon', loc_lat, loc_lon, lat_names[0], lon_names[0])
    wb1         = spital_var(time1,wb1,time_s,time_e)    
    gWwb1         = spital_var(time1,GWwb1,time_s,time_e)    
    
    if 'LIS' in file_paths[0] and var_names[0] in ['WaterTableD_tavg','WatTable']:
        var1     = var1/1000.
        
    # read lat and lon outs
    wrf          = Dataset(wrf_path,  mode='r')
    lons_out     = wrf.variables['XLONG'][0,:,:]
    lats_out     = wrf.variables['XLAT'][0,:,:]

    # offline-CABLE
    time2, Var2  = read_var(file_paths[1], var_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])

    if 'cable_out' in file_paths[1]:
        var2     = spital_var(time2,Var2,time_s,time_e)
    elif 'LIS' in file_paths[1] :
        var2     = Var2[-1]
        
    time, lats_in= read_var(file_paths[1], lat_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
    time, lons_in= read_var(file_paths[1], lon_names[1], loc_lat, loc_lon, lat_names[1], lon_names[1])
    if 'LIS' in file_paths[1] and var_names[1] in ['WaterTableD_tavg','WatTable']:
        var2     = var2/1000.
        
    print("shape lats_in=", np.shape(lats_in))
    print("shape lons_in=", np.shape(lons_in))
    print("shape lats_out=",np.shape(lats_out))
    print("shape lons_out=",np.shape(lons_out))
    print("shape var2=",    np.shape(var2))
    
    if 'cable_out' in file_paths[1] :
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
    else:
        clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]       
    
    if var_names[1] in ['SoilMoist','SoilMoist_inst','ssat','sucs']:
    
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

            plt.savefig('./plots/weather/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
            cb = None
            gl = None
            ax = None
            fig= None
    else:
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

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)
        plt.title(var_names[0], size=16)

        if message == None:
            message = var_names[0]
        else:
            message = message + "_" + var_names[0]

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

    # time_s     = datetime(2017,1,1,0,0,0,0)
    # time_e     = datetime(2017,1,1,23,59,0,0)        

    # message    = "LIS-OFF"
    # file_paths = ["/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/LIS_output/LIS.CABLE.201701-201701.d01.nc",
    #               "/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc"]
    # wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"
    # var_names  = ['SoilSat_inst','ssat']
    # # var_names  = ['SoilMoist_inst','SoilMoist']           
    # # var_names  = ['GWwb_tavg','GWMoist']    
    # # var_names  = ['WaterTableD_tavg','WatTable']
    # lat_names  = ['lat','latitude']
    # lon_names  = ['lon','longitude']
    # spatial_map_single_plot_diff(file_paths, var_names, time_s, time_e, lat_names, lon_names, 
    #                             loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)

    # var_names  = ['Sucs_inst','sucs']
    # # var_names  = ['SoilMoist_inst','SoilMoist']           
    # # var_names  = ['GWwb_tavg','GWMoist']    
    # # var_names  = ['WaterTableD_tavg','WatTable']
    # lat_names  = ['lat','latitude']
    # lon_names  = ['lon','longitude']
    # spatial_map_single_plot_diff(file_paths, var_names, time_s, time_e, lat_names, lon_names, 
    #                             loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path,message=message)