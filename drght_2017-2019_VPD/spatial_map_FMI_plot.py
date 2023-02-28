#!/usr/bin/env python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

'''
History:
1. 8 Feb 2023: cp ./drght_2017-2019_VPD/spatial_map_single_plot.py ./drght_2017-2019_VPD/
'''

import sys
import cartopy
import numpy as np
import shapefile as shp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import scipy.ndimage as ndimage
from netCDF4 import Dataset,num2date
from datetime import datetime, timedelta
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)

from metpy.calc import relative_humidity_from_specific_humidity
from metpy.units import units

from common_utils import *

def plot_spatial_map(filename, var, varname, clevs=None, message="plot"):
    
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
    cmap  = plt.cm.RdBu
            #plt.cm.GnBu
            #plt.cm.seismic #hot_r # BrBG
    
    latlon = Dataset(filename,  mode='r')
    lons   = latlon.variables['longitude'][:,:]
    lats   = latlon.variables['latitude'][:,:]

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
    if np.all(clevs) == None:
        plt.contourf(lons, lats, var, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') 
    else:
        plt.contourf(lons, lats, var, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') 
        
    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    cb.ax.tick_params(labelsize=10)
    plt.title(varname, size=16)

    plt.savefig('./plots/spatial_map_'+message+'.png',dpi=300)

def calc_relative_humidity(filenames, time_s=None, time_e=None, lat_names="latitude", lon_names="longitude",
                           loc_lat=None,loc_lon=None, plot=False):
    
    time, PSurf_tmp = read_var_multi_file(filenames, "PSurf", loc_lat, loc_lon, lat_names, lon_names)
    time, Tair_tmp  = read_var_multi_file(filenames, "Tair", loc_lat, loc_lon, lat_names, lon_names)
    time, Qair_tmp  = read_var_multi_file(filenames, "Qair", loc_lat, loc_lon, lat_names, lon_names)
    
    PSurf           = PSurf_tmp
    Tair            = Tair_tmp - 273.15
    Qair            = Qair_tmp

    RH_tmp          = relative_humidity_from_specific_humidity(PSurf * units.hPa, Tair * units.degC, Qair).to('percent')
    print("RH_tmp",RH_tmp)
    
    rh              = spital_var(time,RH_tmp.magnitude,time_s,time_e)
    tair            = spital_var(time,Tair,time_s,time_e)
    
    print('rh',rh)
    print('tair',tair)
    
    if plot:
        plot_spatial_map(filenames[0], rh, "RH", np.arange(0,100,10), "RH")

    return rh, tair
   
def calc_fuel_moisture_index(land_ctl_files, land_sen_files, time_s=None, time_e=None, lat_names="latitude", lon_names="longitude",loc_lat=None,
                            loc_lon=None, shape_path=None, message=None, plot=False):
    
    rh_ctl, tair_ctl = calc_relative_humidity(land_ctl_files, time_s, time_e, lat_names, lon_names, loc_lat, loc_lon, plot)
    rh_sen, tair_sen = calc_relative_humidity(land_sen_files, time_s, time_e, lat_names, lon_names, loc_lat, loc_lon, plot)
    
    print("np.nanmean(rh_ctl)",np.nanmean(rh_ctl))
    print("np.nanmean(rh_sen)",np.nanmean(rh_sen))
    
    FMI_ctl = 10 - 0.25*(tair_ctl - rh_ctl)
    FMI_sen = 10 - 0.25*(tair_sen - rh_sen)
    
    print("np.nanmean(tair_ctl - rh_ctl)",np.nanmean(30 - rh_ctl))
    print("np.nanmean(tair_sen - rh_sen)",np.nanmean(30 - rh_sen))
    
    FMI_diff = FMI_sen - FMI_ctl
    
    if plot:
        plot_spatial_map(land_ctl_files[0], FMI_ctl, varname="FMI", clevs=None, message=message+"_ctl")
        
    if plot:
        plot_spatial_map(land_ctl_files[0], FMI_diff, varname="FMI", clevs=None, message=message+"_diff")

if __name__ == "__main__":

    # ======================= Option =======================
    # 2017-2019 drought polygon shapefile
    shape_path = "/g/data/w97/ad9701/drought_2017to2020/drought_focusArea/smooth_polygon_drought_focusArea.shp"

    region     = "SE Aus" #"SE Aus" #"CORDEX" #"SE Aus"

    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-23.6]
        loc_lon    = [134,155]
    elif region == "CORDEX":
        loc_lat    = [-52.36,3.87]
        loc_lon    = [89.25,180]

    #######################################################
    # Decks to run:
    #    plot a single map
    #######################################################
    if 1:

        case_name      = "FMI"
        case_ctl       = "default"
        case_sen       = "90th"

        land_sen_path  = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/"+case_sen+"/outputs/"
        land_ctl_path  = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/"+case_ctl+"/outputs/"

        land_sen_files = [ land_sen_path+"cable_out_2017.nc"]
        land_ctl_files = [ land_ctl_path+"cable_out_2017.nc"]

        period         = "2017"
        time_s         = datetime(2017,1,1,0,0,0,0)
        time_e         = datetime(2018,1,1,0,0,0,0)
        message        = case_name+"_"+period # +"_"+case_sen+"-"+case_ctl
        
        calc_fuel_moisture_index(land_ctl_files, land_sen_files, time_s=time_s, time_e=time_e, lat_names="latitude", lon_names="longitude",
                                 loc_lat=loc_lat, loc_lon=loc_lon, shape_path=None, message=message, plot=True)
        
        period         = "Jan2017"
        time_s         = datetime(2017,1,1,0,0,0,0)
        time_e         = datetime(2017,2,1,0,0,0,0)
        message        = case_name+"_"+period # +"_"+case_sen+"-"+case_ctl
        
        calc_fuel_moisture_index(land_ctl_files, land_sen_files, time_s=time_s, time_e=time_e, lat_names="latitude", lon_names="longitude",
                                 loc_lat=loc_lat, loc_lon=loc_lon, shape_path=None, message=message, plot=True)