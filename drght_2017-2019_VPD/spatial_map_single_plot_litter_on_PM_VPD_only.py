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
from common_utils import *

def qair_to_vpd(qair, tair, press):
    '''
    calculate vpd
    '''
    DEG_2_KELVIN = 273.15
    PA_TO_KPA    = 0.001
    PA_TO_HPA    = 0.01

    # convert back to Pa
    press        /= PA_TO_HPA
    tair         -= DEG_2_KELVIN

    # saturation vapor pressure
    es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

    # vapor pressure
    ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair)

    vpd = (es - ea) * PA_TO_KPA
    vpd = np.where(vpd < 0.05, 0.05, vpd)

    return vpd

def read_LIS_vars(var_type):

    '''
    List the variables in a LIS file
    '''

    if var_type == "var_CABLE_3D":
        var_names  =  [ "Rnet","SWnet","LWnet","Qle","Qh","Rainf","Evap","Qs","Qsb","VegT","RadT","ECanop","TVeg","Fwsoil","ESoil","Qair"] # "Albedo","Tair", "Wind",
        # "Rainf_f_inst","Albedo_inst",
        # "GPP_tavg","Qg_tavg","Snowf_tavg","SWE_inst","SnowDepth_inst","SoilWet_inst","CanopInt_inst","SnowCover_inst",
    elif var_type == "var_3D":
        var_names  =  [ "Albedo_inst","Tair_f_inst","Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg",
                        "Rainf_tavg","Evap_tavg","Qs_tavg","Qsb_tavg","VegT_tavg","AvgSurfT_tavg",
                        "ECanop_tavg","TVeg_tavg",
                        "FWsoil_tavg","ESoil_tavg","Wind_f_inst",
                        "Qair_f_inst","Psurf_f_inst","SWdown_f_inst","LWdown_f_inst"]
        # "Rainf_f_inst","Albedo_inst",
        # "GPP_tavg","Qg_tavg","Snowf_tavg","SWE_inst","SnowDepth_inst","SoilWet_inst","CanopInt_inst","SnowCover_inst",
    elif var_type == "var_landinfo_3D":
        var_names  =  [ "Landmask_inst","Landcover_inst","Soiltype_inst","SandFrac_inst","ClayFrac_inst","SiltFrac_inst",
                        "SoilFieldCap_inst","SoilSat_inst","SoilWiltPt_inst","Hyds_inst","Bch_inst","Sucs_inst",
                        "Elevation_inst","LAI_inst"]
    elif var_type == "var_4D":
        var_names  =  ["RelSMC_inst","SoilMoist_inst","SoilTemp_inst","SmLiqFrac_inst","SmFrozFrac_inst"]
    elif var_type == "var_3D_basic":
        var_names  = ["Tair_f_inst",'Evap_tavg',"ESoil_tavg","ECanop_tavg",'TVeg_tavg',"FWsoil_tavg","Qle_tavg","Qh_tavg",
                      "Qg_tavg","VegT_tavg","WaterTableD_tavg","Rainf_tavg","Qs_tavg","Qsb_tavg",]
    elif var_type == "var_4D_basic":
        var_names  =  ["SoilMoist_inst","SoilTemp_inst"]
    elif var_type == "var_energy":
        var_names  = ["Tair_f_inst"] #["Swnet_tavg","Lwnet_tavg","Qle_tavg","Qh_tavg","Qg_tavg","Qair_f_inst","Rnet","EF"] #
    elif var_type == "var_albedo":
        var_names  = ["Albedo_inst"]
    elif var_type == "var_wrf_hgt":
        var_names  = [
                    'cape_3d',# 3D CAPE and CIN
                    'p',    # Full Model Pressure
                    'avo',    # Absolute Vorticity
                    'eth',    # Equivalent Potential Temperature
                    'dbz',    # Reflectivity
                    'geopt',  # Geopotential for the Mass Grid
                    'omg',  # Omega
                    'pvo',  # Potential Vorticity
                    'rh',   # Relative Humidity
                    'td',   # Dew Point Temperature
                    'tc',   # Temperature in Celsius
                    'th',   # Potential Temperature
                    'temp', # Temperature (in specified units)
                    'tv',   # Virtual Temperature
                    'twb',  # Wet Bulb Temperature
                    'ua',   # U-component of Wind on Mass Points
                    'va',   # V-component of Wind on Mass Points
                    'wa',   # W-component of Wind on Mass Points
                    'z',    # Model Height for Mass Grid
                    ]
    elif var_type == "var_wrf_hgt":
        var_names  = [
                    'cape_3d',# 3D CAPE and CIN
                    'p',    # Full Model Pressure
                    'avo',    # Absolute Vorticity
                    'eth',    # Equivalent Potential Temperature
                    'dbz',    # Reflectivity
                    'geopt',  # Geopotential for the Mass Grid
                    'omg',  # Omega
                    'pvo',  # Potential Vorticity
                    'rh',   # Relative Humidity
                    'td',   # Dew Point Temperature
                    'tc',   # Temperature in Celsius
                    'th',   # Potential Temperature
                    'temp', # Temperature (in specified units)
                    'tv',   # Virtual Temperature
                    'twb',  # Wet Bulb Temperature
                    'ua',   # U-component of Wind on Mass Points
                    'va',   # V-component of Wind on Mass Points
                    'wa',   # W-component of Wind on Mass Points
                    'z',    # Model Height for Mass Grid
                    ]
    elif var_type == "var_wrf_surf":
        var_names = [
                    'cloudfrac', # Cloud Fraction
                    'td2',  # 2m Dew Point Temperature
                    'rh2',  # 2m Relative Humidity
                    'T2',   # 2m Temperature
                    'slp',  # Sea Level Pressure
                    'ter',  # Model Terrain Height
                    'updraft_helicity', # Updraft Helicity
                    'helicity',        # Storm Relative Helicity
                    'ctt',  # Cloud Top Temperature
                    'mdbz', # Maximum Reflectivity
                    'td2',  # 2m Dew Point Temperature
                    'rh2',  # 2m Relative Humidity
                    'T2',   # 2m Temperature
                    'slp',  # Sea Level Pressure
                    'pw',   # Precipitable Water
                    'cape_2d', # 2D CAPE (MCAPE/MCIN/LCL/LFC)
                    'cloudfrac', # Cloud Fraction
                ]
    elif var_type == "var_wrf_surf_basic":
        var_names = [
                    'cloudfrac', # Cloud Fraction
                    'td2',  # 2m Dew Point Temperature
                    'rh2',  # 2m Relative Humidity
                    'T2',   # 2m Temperature
                    'slp',  # Sea Level Pressure
                    'td2',  # 2m Dew Point Temperature
                    'rh2',  # 2m Relative Humidity
                    'T2',   # 2m Temperature
                    'slp',  # Sea Level Pressure
                    'pw',   # Precipitable Water
                    'cape_2d', # 2D CAPE (MCAPE/MCIN/LCL/LFC)
                ]
    elif var_type == "var_wrf_surf_other":
        var_names  = [  "SWDNB", # INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM
                        "LWDNB", # INSTANTANEOUS DOWNWELLING LONGWAVE FLUX AT BOTTOM
                        "SWUPB", # INSTANTANEOUS UPWELLING SHORTWAVE FLUX AT BOTTOM
                        "LWUPB", # INSTANTANEOUS UPWELLING LONGWAVE FLUX AT BOTTOM
                        ]
        # ['RAINC','RAINNC','PSFC','U10','V10','TSK','PBLH']
    return var_names

def spatial_map_single_plot_diff_multifile(land_ctl_files, land_sen_files, var_names, time_s=None,
                                     time_e=None, lat_names="latitude", lon_names="longitude",loc_lat=None,
                                     loc_lon=None, shape_path=None, message=None, case_sen=None):

    '''
    plot a single spatial map
    '''

    # CABLE
    for var_name in var_names:
        if var_name in ["VPD","VPD_rate"]:
            time, Tair_ctl = read_var_multi_file(land_ctl_files, "Tair", loc_lat, loc_lon, lat_names, lon_names)
            time, Qair_ctl = read_var_multi_file(land_ctl_files, "Qair", loc_lat, loc_lon, lat_names, lon_names)
            time, Pres_ctl = read_var_multi_file(land_ctl_files, "PSurf", loc_lat, loc_lon, lat_names, lon_names)
            Ctl_tmp        = qair_to_vpd(Qair_ctl, Tair_ctl, Pres_ctl)

            time, Tair_sen = read_var_multi_file(land_sen_files, "Tair", loc_lat, loc_lon, lat_names, lon_names)
            time, Qair_sen = read_var_multi_file(land_sen_files, "Qair", loc_lat, loc_lon, lat_names, lon_names)
            time, Pres_sen = read_var_multi_file(land_sen_files, "PSurf", loc_lat, loc_lon, lat_names, lon_names)
            Sen_tmp        = qair_to_vpd(Qair_sen, Tair_sen, Pres_sen)
        else:
            time, Ctl_tmp = read_var_multi_file(land_ctl_files, var_name, loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_tmp = read_var_multi_file(land_sen_files, var_name, loc_lat, loc_lon, lat_names, lon_names)

        Ctl_var       = spital_var(time,Ctl_tmp,time_s,time_e)
        Sen_var       = spital_var(time,Sen_tmp,time_s,time_e)

        latlon        = Dataset(land_ctl_files[0],  mode='r')
        lons          = latlon.variables['longitude'][:,:]
        lats          = latlon.variables['latitude'][:,:]

        if var_name in ['WaterTableD_tavg','WatTable']:
            Ctl_var   = Ctl_var/1000.
            Sen_var   = Sen_var/1000.
        if var_name in ['ESoil_tavg','Evap_tavg',"ECanop_tavg",'TVeg_tavg',"Rainf_tavg","Snowf_tavg","Qs_tavg","Qsb_tavg",
                        'ESoil','Evap',"ECanop",'TVeg',"Rainf","Snowf","Qs","Qsb","PotEvap"]:
            Ctl_var   = Ctl_var*3600*24*(time_e-time_s).days
            Sen_var   = Sen_var*3600*24*(time_e-time_s).days
            if (time_e-time_s).days > 365:
                Ctl_var = Ctl_var/3.
                Sen_var = Sen_var/3.
        if var_name in ['Qair_f_inst','Qair']:
            Ctl_var   = Ctl_var*1000
            Sen_var   = Sen_var*1000
        if var_name in ["VegT_tavg","AvgSurfT_tavg","CanopInt_inst","SnowCover_inst", "SoilTemp_inst",
                        "VegT","AvgSurfT","CanopInt","SnowCover", "SoilTemp","RadT"]:
            Ctl_var   = Ctl_var -273.15
            Sen_var   = Sen_var -273.15
        if var_name in ["GPP","NPP"]:
            s2d       = 3600*24. # s-1 to d-1
            GPP_scale = 0.000001*12*s2d # umol s-1 to g d-1
            Ctl_var   = Ctl_var*GPP_scale*(time_e-time_s).days
            Sen_var   = Sen_var*GPP_scale*(time_e-time_s).days
            if (time_e-time_s).days > 365:
                Ctl_var = Ctl_var/3.
                Sen_var = Sen_var/3.
        if var_name in ["Gs"]:
            s2d       = 3600*24. # s-1 to d-1
            Ctl_var   = Ctl_var*s2d
            Sen_var   = Sen_var*s2d      
        if var_name in ["VPD_rate"]:  
            var_diff  = (Sen_var - Ctl_var)/Ctl_var*100
        else:    
            var_diff  = Sen_var - Ctl_var

        # read lat and lon outs
        cmap  = plt.cm.RdBu_r
        #plt.cm.GnBu
        #plt.cm.seismic #hot_r # BrBG

        if var_name in ['WaterTableD_tavg','WatTable']:
            clevs = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]
        elif var_name in ['GWwb_tavg','GWMoist']:
            clevs = [-0.05,-0.04,-0.03,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.03,0.04,0.05]
        elif  var_name in ["Qair_f_inst","Qair"]:
            clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2]
            # clevs = [-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14,16,18,20]
            # clevs = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        elif var_name in ['SoilMoist_inst','SoilMoist']:
            clevs = [-0.1,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
        elif var_name in ['ESoil_tavg','Evap_tavg',"ECanop_tavg",'TVeg_tavg',"Rainf_tavg","Snowf_tavg","Qs_tavg","Qsb_tavg",
                          'ESoil','Evap',"ECanop",'TVeg',"Rainf","Snowf","Qs","Qsb"]:
            # clevs = [-5.,-4.5,-4.,-3.5,-3.,-2.5,-2,-1.5,-1,-0.5,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.]
            # clevs = [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,-5,5,10,20.,30,40,50,60,70,80,90,100]
            # clevs = [-140,-120,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20.,30,40,50,60,70,80,90,100,120,140]
            # clevs = [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,-5,5,10,20.,30,40,50,60,70,80,90,100]
            clevs = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20.,25,30,35,40,45,50]
            cmap  = plt.cm.RdBu_r #plt.cm.GnBu
            # clevs = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5, 5,10,15,20,25,30,35,40,45,50]
        elif var_name in ["VPD_rate"]:  
            clevs = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20.,25,30,35,40,45,50]
        elif var_name in ["PotEvap"]:
            clevs = [-200,-180,-160,-140,-120,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20.,30,40,50,60,70,80,90,100,120,140,160,180,200]
        elif var_name in ["Qle_tavg","Qh_tavg","Qg_tavg","Qle","Qh","Qg"]:
            # clevs = [-140,-120,-100,-80,-60,-40,-20,-10,-5,5,10,20,40,60,80,100,120,140]
            # clevs = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5, 5,10,15,20,25,30,35,40,45,50]
            clevs = [-20,-18,-16,-14,-12,-10,-8,-6,-4,-2, 2,4,6,8,10,12,14,16,18,20]
            # clevs = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 1,2,3,4,5,6,7,8,9,10]
        elif var_name in ["Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst","Rnet","SWnet","LWnet","SWdown","LWdown"]:
            # clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
            clevs = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 1,2,3,4,5,6,7,8,9,10]
        elif var_name in ["Tair_f_inst","Tair","RadT","VegT","VegT_tavg","AvgSurfT_tavg","CanopInt_inst","SnowCover_inst", "SoilTemp_inst",
                          "AvgSurfT","CanopInt","SnowCover", "SoilTemp"]:
            clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.]
            cmap  = plt.cm.RdBu_r
        elif var_name in ["Wind_f_inst","Wind"]:
            clevs = [-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4]
        elif var_name in ["Psurf_f_inst","Psurf"]:
            clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
        elif var_name in ["Tair_f_inst","Tair",]:
            # clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.]
            clevs = [-1,-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8,1.]
            # clevs = [-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1.,1.5,2,2.5,3.]
        elif var_name in ["FWsoil_tavg","SmLiqFrac_inst","SmFrozFrac_inst","Fwsoil"]:
            clevs = [-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45]
        elif  var_name in ["GPP","NPP"]:
            clevs =  [-140,-120,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20.,30,40,50,60,70,80,90,100,120,140]
        elif  var_name in ["Gs"]:
            clevs = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20.,25,30,35,40,45,50]
        else:
            clevs = [-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5]

        print("len(np.shape(var_diff))",len(np.shape(var_diff)))

        if len(np.shape(var_diff)) >=3:

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
                plt.contourf(lons, lats, var_diff[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #

                cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
                cb.ax.tick_params(labelsize=10)
                plt.title(var_name, size=16)

                if shape_path != None:
                    # Requires the pyshp package
                    sf = shp.Reader(shape_path)

                    for shape in sf.shapeRecords():
                        x = [i[0] for i in shape.shape.points[:]]
                        y = [i[1] for i in shape.shape.points[:]]
                        plt.plot(x,y,c="black")

                if j == 0:
                    if message == None:
                        message = var_name
                    else:
                        message = message + "_" + var_name

                plt.savefig('./plots/'+case_sen+'/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
                cb = None
                gl = None
                ax = None
                fig= None

        elif len(np.shape(var_diff)) ==2:
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

            plt.contourf(lons, lats, var_diff, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') #
            print(var_diff)
            cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
            cb.ax.tick_params(labelsize=10)
            plt.title(var_name, size=16)

            if shape_path != None:
                # Requires the pyshp package
                sf = shp.Reader(shape_path)

                for shape in sf.shapeRecords():
                    x = [i[0] for i in shape.shape.points[:]]
                    y = [i[1] for i in shape.shape.points[:]]
                    plt.plot(x,y,c="black")
            if case_sen == None:
                plt.savefig('./plots/spatial_map_'+message + "_" + var_name+'.png',dpi=300)
            else:
                plt.savefig('./plots/'+case_sen+'/spatial_map_'+message + "_" + var_name+'.png',dpi=300)

def spatial_map_single_plot_multifile(land_ctl_files, var_names, time_s=None,
                                     time_e=None, lat_names="latitude", lon_names="longitude",loc_lat=None,
                                     loc_lon=None, shape_path=None, message=None):

    '''
    plot a single spatial map
    '''

    # CABLE
    for var_name in var_names:
        time, Ctl_tmp = read_var_multi_file(land_ctl_files, var_name, loc_lat, loc_lon, lat_names, lon_names)
        Ctl_var       = spital_var(time,Ctl_tmp,time_s,time_e)
        latlon        = Dataset(land_ctl_files[0],  mode='r')
        lons          = latlon.variables['longitude'][:,:]
        lats          = latlon.variables['latitude'][:,:]

        if var_name in ['WaterTableD_tavg','WatTable']:
            Ctl_var   = Ctl_var/1000.
        if var_name in ['ESoil_tavg','Evap_tavg',"ECanop_tavg",'TVeg_tavg',"Rainf_tavg","Snowf_tavg","Qs_tavg","Qsb_tavg",
                        'ESoil','Evap',"ECanop",'TVeg',"Rainf","Snowf","Qs","Qsb","PotEvap"]:
            Ctl_var   = Ctl_var*3600*24*(time_e-time_s).days
        if var_name in ['Qair_f_inst','Qair']:
            Ctl_var   = Ctl_var*1000
        if var_name in ["VegT_tavg","AvgSurfT_tavg","CanopInt_inst","SnowCover_inst", "SoilTemp_inst",
                          "VegT","AvgSurfT","CanopInt","SnowCover", "SoilTemp","RadT"]:
            Ctl_var   = Ctl_var -273.15
        if var_name in ["GPP","NPP"]:
            s2d       = 3600*24. # s-1 to d-1
            GPP_scale = 0.000001*12*s2d # umol s-1 to g d-1
            Ctl_var   = Ctl_var*GPP_scale*(time_e-time_s).days
        if var_name in ["Gs"]:
            s2d       = 3600*24. # s-1 to d-1
            Ctl_var   = Ctl_var*s2d
        # read lat and lon outs

        cmap  = plt.cm.GnBu
        if var_name in ['WaterTableD_tavg','WatTable']:
            clevs = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]
        elif var_name in ['GWwb_tavg','GWMoist']:
            clevs = [-0.05,-0.04,-0.03,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.03,0.04,0.05]
        elif  var_name in ["Qair_f_inst","Qair"]:
            # clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2]
            clevs = [-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14,16,18,20]
            # clevs = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
        elif var_name in ['SoilMoist_inst','SoilMoist']:
            clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]
        elif var_name in ['ESoil_tavg','Evap_tavg',"ECanop_tavg",'TVeg_tavg',"Rainf_tavg","Snowf_tavg","Qs_tavg","Qsb_tavg",
                          'ESoil','Evap',"ECanop",'TVeg',"Rainf","Snowf","Qs","Qsb","PotEvap"]:
            clevs = np.arange(0,1000,20)
            cmap  = plt.cm.GnBu
        elif var_name in ["Qle_tavg","Qh_tavg","Qg_tavg","Qle","Qh","Qg"]:
            clevs = [-40,-20,-10,-5,5,10,20,40,60,80,100,120,140,160,180,200,220,240]
            # clevs = [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5, 5,10,15,20,25,30,35,40,45,50]
            # clevs = [-20,-18,-16,-14,-12,-10,-8,-6,-4,-2, 2,4,6,8,10,12,14,16,18,20]
            # clevs = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 1,2,3,4,5,6,7,8,9,10]
        elif var_name in ["Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst","Rnet","SWnet","LWnet","SWdown","LWdown"]:
            # clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
            clevs = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 1,2,3,4,5,6,7,8,9,10]
        elif var_name in ["VegT_tavg","AvgSurfT_tavg","CanopInt_inst","SnowCover_inst", "SoilTemp_inst",
                          "VegT","AvgSurfT","CanopInt","SnowCover", "SoilTemp",'RadT']:
            clevs = np.arange(10,40,2)
            cmap  = plt.cm.RdBu_r
        elif var_name in ["Wind_f_inst","Wind"]:
            clevs = [-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4]
        elif var_name in ["GPP","NPP"]:
            clevs = np.arange(0,2000,20)
        elif var_name in ["Psurf_f_inst","Psurf"]:
            clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
        elif var_name in ["Tair_f_inst","Tair",]:
            # clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.]
            clevs = [-1,-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8,1.]
            # clevs = [-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1.,1.5,2,2.5,3.]
        elif var_name in ["FWsoil_tavg","SmLiqFrac_inst","SmFrozFrac_inst","Fwsoil"]:
            clevs = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]
        elif var_name in ["Gs"]:
            clevs = np.arange(0,400,20)
        else:
            clevs = [-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4,0.5]

        print("len(np.shape(Ctl_var))",len(np.shape(Ctl_var)))

        if len(np.shape(Ctl_var)) >=3:

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
                plt.contourf(lons, lats, Ctl_var[j,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
                cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
                cb.ax.tick_params(labelsize=10)
                plt.title(var_name, size=16)

                if shape_path != None:
                    # Requires the pyshp package
                    sf = shp.Reader(shape_path)

                    for shape in sf.shapeRecords():
                        x = [i[0] for i in shape.shape.points[:]]
                        y = [i[1] for i in shape.shape.points[:]]
                        plt.plot(x,y,c="black")

                if j == 0:
                    if message == None:
                        message = var_name
                    else:
                        message = message + "_" + var_name

                plt.savefig('./plots/ctl/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
                cb = None
                gl = None
                ax = None
                fig= None

        elif len(np.shape(Ctl_var)) ==2:
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

            plt.contourf(lons, lats, Ctl_var, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
            print(Ctl_var)
            cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
            cb.ax.tick_params(labelsize=10)
            plt.title(var_name, size=16)

            if shape_path != None:
                # Requires the pyshp package
                sf = shp.Reader(shape_path)

                for shape in sf.shapeRecords():
                    x = [i[0] for i in shape.shape.points[:]]
                    y = [i[1] for i in shape.shape.points[:]]
                    plt.plot(x,y,c="black")

            plt.savefig('./plots/ctl/spatial_map_'+message + "_" + var_name+'.png',dpi=300)

def sen_vs_ctl(case_ctl,case_sen,loc_lat=None,loc_lon=None,sen_vs_ctl=None):

    land_ctl_path  = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/"+case_ctl+"/"+config_set[0]+"/outputs/"
    land_sen_path  = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/"+case_sen+"/"+config_set[1]+"/outputs/"

    land_ctl_files = [ land_ctl_path+"cable_out_2017.nc",
                       land_ctl_path+"cable_out_2018.nc",
                       land_ctl_path+"cable_out_2019.nc",
                       ]

    land_sen_files = [ land_sen_path+"cable_out_2017.nc",
                       land_sen_path+"cable_out_2018.nc",
                       land_sen_path+"cable_out_2019.nc",
                       ]

    var_names      =  [ "Gs","GPP","Qle","Qh","RadT","Rnet","VegT","Evap","TVeg",
                        "ESoil","Fwsoil","Tair","Qair","VPD","VPD_rate"] #"SoilMoist",

    period         = "2019_fire"
    time_s         = datetime(2019,9,1,0,0,0,0)
    time_e         = datetime(2020,1,1,0,0,0,0)
    message        = period+"-"+case_sen+"-"+case_ctl+"_"+config_set[1]
    spatial_map_single_plot_diff_multifile(land_ctl_files, land_sen_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                    lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message, case_sen =case_sen+"/"+config_set[1])

    # period         = "2017"
    # time_s         = datetime(2017,1,1,0,0,0,0)
    # time_e         = datetime(2018,1,1,0,0,0,0)
    # message        = period+"-"+case_sen+"-"+case_ctl+"_"+config_set[1]
    # spatial_map_single_plot_diff_multifile(land_ctl_files, land_sen_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
    #                                 lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message,case_sen=case_sen)

    # period         = "2018"
    # time_s         = datetime(2018,1,1,0,0,0,0)
    # time_e         = datetime(2019,1,1,0,0,0,0)
    # message        =  period+"-"+case_sen+"-"+case_ctl+"_"+config_set[1]
    # spatial_map_single_plot_diff_multifile(land_ctl_files, land_sen_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
    #                                 lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message,case_sen=case_sen)

    period         = "2019"
    time_s         = datetime(2019,1,1,0,0,0,0)
    time_e         = datetime(2020,1,1,0,0,0,0)
    message        =  period+"-"+case_sen+"-"+case_ctl+"_"+config_set[1]
    spatial_map_single_plot_diff_multifile(land_ctl_files, land_sen_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                    lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message, case_sen =case_sen+"/"+config_set[1])


    period         = "2017-2019"
    time_s         = datetime(2017,1,1,0,0,0,0)
    time_e         = datetime(2020,1,1,0,0,0,0)
    message        = period+"-"+case_sen+"-"+case_ctl+"_"+config_set[1]
    spatial_map_single_plot_diff_multifile(land_ctl_files, land_sen_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                    lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message, case_sen =case_sen+"/"+config_set[1])

if __name__ == "__main__":

    # ======================= Option =======================
    # 2017-2019 drought polygon shapefile
    shape_path = "/g/data/w97/ad9701/drought_2017to2020/drought_focusArea/smooth_polygon_drought_focusArea.shp"

    region     = "SE Aus" #"Aus" #"SE Aus" #"CORDEX" #"SE Aus"

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


    if 0:
        # ctl sim
        case_ctl       = "ctl"
        land_ctl_path  = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/"+case_ctl+"/outputs/"
        land_ctl_files = [ land_ctl_path+"cable_out_2017.nc",
                           land_ctl_path+"cable_out_2018.nc",
                           land_ctl_path+"cable_out_2019.nc"]

        var_names      = ["VPD"]
        #["Gs","GPP","NPP","Qle","Qh","RadT","VegT","Evap","TVeg","ESoil","SoilMoist","Fwsoil","Tair","Qair",]

        period         = "2019_fire"
        time_s         = datetime(2019,9,1,0,0,0,0)
        time_e         = datetime(2020,1,1,0,0,0,0)
        message        = case_ctl+"_"+period
        spatial_map_single_plot_multifile(land_ctl_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                        lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message)

        period         = "2017"
        time_s         = datetime(2017,1,1,0,0,0,0)
        time_e         = datetime(2018,1,1,0,0,0,0)
        message        = case_ctl+"_"+period
        spatial_map_single_plot_multifile(land_ctl_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                        lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message)

        period         = "2018"
        time_s         = datetime(2018,1,1,0,0,0,0)
        time_e         = datetime(2019,1,1,0,0,0,0)
        message        = case_ctl+"_"+period
        spatial_map_single_plot_multifile(land_ctl_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                        lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message)

        period         = "2019"
        time_s         = datetime(2019,1,1,0,0,0,0)
        time_e         = datetime(2020,1,1,0,0,0,0)
        message        = case_ctl+"_"+period
        spatial_map_single_plot_multifile(land_ctl_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                        lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message)

        period         = "2017-2019"
        time_s         = datetime(2017,1,1,0,0,0,0)
        time_e         = datetime(2020,1,1,0,0,0,0)
        message        = case_ctl+"_"+period
        spatial_map_single_plot_multifile(land_ctl_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                        lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message)

    if 1:
        config_set     = ["litter_on_PM","litter_on_PM_VPD_only"]
        #["litter_off_HDM","litter_off_HDM"]
        
        case_ctl       = "ctl"
        case_sen       = "T_Q_LWdown_detrend_2017_2019"

        sen_vs_ctl(case_ctl,case_sen,loc_lat,loc_lon,config_set) #

        # case_sen       = "Q_detrend_2000_2019"
        # sen_vs_ctl(case_ctl,case_sen)

        # case_sen       = "T_detrend_2000_2019"
        # sen_vs_ctl(case_ctl,case_sen)

        #case_sen       = "T_Q_detrend_2000_2019"
        #sen_vs_ctl(case_ctl,case_sen)

        # case_sen       = "Q_detrend_2017_2019"
        # sen_vs_ctl(case_ctl,case_sen)

        # case_sen       = "T_detrend_2017_2019"
        # sen_vs_ctl(case_ctl,case_sen)

    if 0:
        # compare old VPD sims
        case_name      = "VPD"
        case_ctl       = "default"
        case_sen       = "80th"

        land_ctl_path  = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/"+case_ctl+"/outputs/"

        land_ctl_files = [ land_ctl_path+"cable_out_2017.nc"]
        land_sen_files = [ land_sen_path+"cable_out_2017.nc"]

        message    = case_name+"_"+case_sen+"-"+case_ctl+"_"+period
        spatial_map_single_plot_multifile(land_ctl_files, land_sen_files, var_names, time_s=time_s, time_e=time_e, lat_names="latitude",
                                        lon_names="longitude",loc_lat=loc_lat, loc_lon=loc_lon, shape_path=shape_path, message=message)
