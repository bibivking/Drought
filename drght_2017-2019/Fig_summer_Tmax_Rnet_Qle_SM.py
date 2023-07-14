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
import shapefile as shp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
import scipy.ndimage as ndimage
from scipy.interpolate import griddata, interp1d
from netCDF4 import Dataset,num2date
from datetime import datetime, timedelta
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                        cartopy_ylim, latlon_coords, ALL_TIMES)
from common_utils import *

def read_LIS_var(file_name, land_ctl_path, land_sen_path, var_name, loc_lat, loc_lon, lat_names, lon_names, time_s, time_e):

    if var_name in ["Tmax","Tmin","TDR"]:
        land_ctl_files= [land_ctl_path+'Tair_f_inst/'+file_name]
        land_sen_files= [land_sen_path+'Tair_f_inst/'+file_name]
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
        Ctl_tmp       = Ctl_tmp - 273.15
        Sen_tmp       = Sen_tmp - 273.15
    elif var_name in ["Rnet",]:
        land_ctl_files= [land_ctl_path+'Lwnet_tavg/'+file_name]
        land_sen_files= [land_sen_path+'Lwnet_tavg/'+file_name]
        time, Ctl_Lwnet_tmp = read_var_multi_file(land_ctl_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_Lwnet_tmp = read_var_multi_file(land_sen_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        land_ctl_files= [land_ctl_path+'Swnet_tavg/'+file_name]
        land_sen_files= [land_sen_path+'Swnet_tavg/'+file_name]
        time, Ctl_Swnet_tmp = read_var_multi_file(land_ctl_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_Swnet_tmp = read_var_multi_file(land_sen_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
        Ctl_tmp = Ctl_Lwnet_tmp+Ctl_Swnet_tmp
        Sen_tmp = Sen_Lwnet_tmp+Sen_Swnet_tmp
    elif var_name in ["SM_top50cm",]:
        land_ctl_files = [land_ctl_path+'SoilMoist_inst/'+file_name]
        land_sen_files = [land_sen_path+'SoilMoist_inst/'+file_name]
        time, Ctl_temp = read_var_multi_file(land_ctl_files, 'SoilMoist_inst', loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_temp = read_var_multi_file(land_sen_files, 'SoilMoist_inst', loc_lat, loc_lon, lat_names, lon_names)
        # [.022, .058, .154, .409, 1.085, 2.872]
        Ctl_tmp    = (Ctl_temp[:,0,:,:]*0.022 + Ctl_temp[:,1,:,:]*0.058 + Ctl_temp[:,2,:,:]*0.154 + Ctl_temp[:,3,:,:]*0.266)/0.5
        Sen_tmp    = (Sen_temp[:,0,:,:]*0.022 + Sen_temp[:,1,:,:]*0.058 + Sen_temp[:,2,:,:]*0.154 + Sen_temp[:,3,:,:]*0.266)/0.5
    else:
        land_ctl_files= [land_ctl_path+var_name+'/LIS.CABLE.201701-202002.nc']
        land_sen_files= [land_sen_path+var_name+'/LIS.CABLE.201701-202002.nc']
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, var_name, loc_lat, loc_lon, lat_names, lon_names)
        time, Sen_tmp = read_var_multi_file(land_sen_files, var_name, loc_lat, loc_lon, lat_names, lon_names)

    if 'max' in var_name:
        # average of daily max
        ctl_in       = spatial_var_max(time,Ctl_tmp,time_s,time_e)
        sen_in       = spatial_var_max(time,Sen_tmp,time_s,time_e)
    else:
        ctl_in       = spatial_var(time,Ctl_tmp,time_s,time_e)
        sen_in       = spatial_var(time,Sen_tmp,time_s,time_e)
    var_diff  = sen_in - ctl_in

    # =============== CHANGE HERE ===============
    cmap  = plt.cm.BrBG
    if var_name in ['SoilMoist_inst','SoilMoist',"SM_top50cm"]:
        clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]
    elif var_name in ["Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst","Qle_tavg","Qh_tavg","Qg_tavg"]:
        clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
    elif var_name in ["Rnet"]:
        clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
        cmap  = plt.cm.BrBG_r
    elif var_name in ["Tair_f_inst","Tmax","Tmin","VegT_tavg","VegTmax","VegTmin",
                        "AvgSurfT_tavg","SurfTmax","SurfTmin","SoilTemp_inst",'TDR','VegTDR','SurfTDR']:
        clevs = [-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2]
        cmap  = plt.cm.seismic
    elif var_name in ["FWsoil_tavg","SmLiqFrac_inst","SmFrozFrac_inst"]:
        clevs = [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35]
    else:
        clevs = [-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5]

    return var_diff, cmap, clevs

def spatial_map_summer_Tmax_Rnet_Qle_SM(file_name, land_ctl_path, land_sen_path, time_ss=None,
                              time_es=None, lat_names="lat", lon_names="lon",loc_lat=None,
                              loc_lon=None, wrf_path=None,  message=None):

    '''
    plot a single spatial map
    '''

    # Read lat and lon outs
    wrf            = Dataset(wrf_path,  mode='r')
    lons           = wrf.variables['XLONG'][0,:,:]
    lats           = wrf.variables['XLAT'][0,:,:]

    # ================== Start Plotting =================
    fig, axs = plt.subplots(nrows=4, ncols=3, figsize=[10,12],sharex=True,
                sharey=True, squeeze=True, subplot_kw={'projection': ccrs.PlateCarree()})

    plt.subplots_adjust(wspace=0, hspace=0.)

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

    states= NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")

    # WRF-CABLE
    for i in np.arange(3):

        Tmax_diff, cmap, clevs = read_LIS_var(file_name, land_ctl_path, land_sen_path, 'Tmax', loc_lat, loc_lon, lat_names, lon_names, time_ss[i], time_es[i])
        Rnet_diff, cmap, clevs = read_LIS_var(file_name, land_ctl_path, land_sen_path, 'Rnet', loc_lat, loc_lon, lat_names, lon_names, time_ss[i], time_es[i])
        Qle_diff, cmap, clevs  = read_LIS_var(file_name, land_ctl_path, land_sen_path, 'Qle', loc_lat, loc_lon, lat_names, lon_names, time_ss[i], time_es[i])
        SM_diff, cmap, clevs   = read_LIS_var(file_name, land_ctl_path, land_sen_path, 'SM_top50cm', loc_lat, loc_lon, lat_names, lon_names, time_ss[i], time_es[i])

        for i in np.arange(4):

            axs[i,j].coastlines(resolution="50m",linewidth=1)
            axs[i,j].set_extent([135,155,-39,-24])
            axs[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = axs[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
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

        if i == 3:
            gl.xlabels_bottom = True
        else:
            gl.xlabels_bottom = False

        if j == 0:
            gl.ylabels_left   = True
        else:
            gl.ylabels_left   = False

        plot1 = axs[0,j].contourf(lons, lats, Tmax_diff[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
        plot2 = axs[1,j].contourf(lons, lats, Rnet_diff[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
        plot3 = axs[2,j].contourf(lons, lats, Qle_diff[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
        plot4 = axs[3,j].contourf(lons, lats, SM_diff[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')

    # Adding colormap
    cbar = plt.colorbar(plot1, ax=axs[0], ticklocation="right", pad=0.08, orientation="horizontal",
        aspect=40, shrink=1) # cax=cax,
    cbar.set_label('$\mathregular{^{o}}$C', loc='center',size=12)# rotation=270,

    cbar = plt.colorbar(plot2, ax=axs[1], ticklocation="right", pad=0.08, orientation="horizontal",
        aspect=40, shrink=1) # cax=cax,
    cbar.set_label('W m$\mathregular{^{-2}}$', loc='center',size=12)# rotation=270,

    cbar = plt.colorbar(plot3, ax=axs[2], ticklocation="right", pad=0.08, orientation="horizontal",
        aspect=40, shrink=1) # cax=cax,
    cbar.set_label('W m$\mathregular{^{-2}}$', loc='center',size=12)# rotation=270,

    cbar = plt.colorbar(plot4, ax=axs[3], ticklocation="right", pad=0.08, orientation="horizontal",
        aspect=40, shrink=1) # cax=cax,
    cbar.set_label('m$\mathregular{^{3}}$ m$\mathregular{^{-3}}$', loc='center',size=12)# rotation=270,

    cbar.ax.tick_params(labelsize=12)#,labelrotation=45)

    # Add titles
    axs[0,0].set_title("2017-18", fontsize=12)
    axs[0,1].set_title("2018-19", fontsize=12)
    axs[0,2].set_title("2019-20", fontsize=12)

    axs[0,0].set_ylabel("ΔT$\mathregular{_{max}}$")
    axs[1,0].set_ylabel("ΔR$\mathregular{_{net}}$")
    axs[2,0].set_ylabel("ΔLH")
    axs[3,0].set_ylabel("ΔSM$\mathregular{_{top}}$")

    # Apply tight layout
    # plt.tight_layout()
    plt.savefig('./plots/spatial_map_summer_Tmax_Rnet_Qle_SM.png',dpi=300)


if __name__ == "__main__":

    # ======================= Option =======================
    region = "SE Aus" #"SE Aus" #"CORDEX" #"SE Aus"

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
        '''
        Test WRF-CABLE LIS output
        '''

        case_name      = "ALB-CTL_new" #"bl_pbl2_mp4_sf_sfclay2" #
        case_ctl       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
        case_sen       = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/WRF_output/wrfout_d01_2019-12-01_01:00:00"
        land_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/LIS_output/"
        land_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/LIS_output/"
        atmo_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/WRF_output/"
        atmo_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/WRF_output/"

        time_ss    = [datetime(2017,12,1,0,0,0,0),
                      datetime(2018,12,1,0,0,0,0),
                      datetime(2019,12,1,0,0,0,0)]

        time_es    = [datetime(2018,3,1,0,0,0,0),
                      datetime(2019,3,1,0,0,0,0),
                      datetime(2020,3,1,0,0,0,0)]

        message    = "Summer_Tmax_Rnet_Qle_SM"

        file_name  = "LIS.CABLE.201701-202002.nc"

        spatial_map_summer_Tmax_Rnet_Qle_SM(file_name, land_ctl_path, land_sen_path, time_ss=time_ss, time_es=time_es, lat_names="lat",
                            lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path, message=message)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              