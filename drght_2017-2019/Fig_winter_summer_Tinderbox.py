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

def read_LIS_vars(var_type):

    '''
    List the variables in a LIS file
    '''

    if var_type == "var_3D":
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

def spatial_map_single_plot_LIS_diff(land_ctl_path, land_sen_path, var_names, time_ss=None,
                                     time_es=None, lat_names="lat", lon_names="lon",loc_lat=None,
                                     loc_lon=None, wrf_path=None, shape_path=None, message=None):

    '''
    plot a single spatial map
    '''

    # WRF-CABLE
    for var_name in var_names:
        print("plotting "+var_name)

        if var_name in ["Tmax","Tmin","TDR"]:
            land_ctl_files= [land_ctl_path+'Tair_f_inst/LIS.CABLE.201701-202002.nc']
            land_sen_files= [land_sen_path+'Tair_f_inst/LIS.CABLE.201701-202002.nc']
            time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_tmp = read_var_multi_file(land_sen_files, 'Tair_f_inst', loc_lat, loc_lon, lat_names, lon_names)
            Ctl_tmp       = Ctl_tmp -273.15
            Sen_tmp       = Sen_tmp -273.15
        elif var_name in ["VegTmax","VegTmin","VegTDR"]:
            land_ctl_files= [land_ctl_path+'VegT_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files= [land_sen_path+'VegT_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_tmp = read_var_multi_file(land_sen_files, 'VegT_tavg', loc_lat, loc_lon, lat_names, lon_names)
            Ctl_tmp       = Ctl_tmp -273.15
            Sen_tmp       = Sen_tmp -273.15
        elif var_name in ["SurfTmax","SurfTmin","SurfTDR"]:
            land_ctl_files= [land_ctl_path+'AvgSurfT_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files= [land_sen_path+'AvgSurfT_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_tmp = read_var_multi_file(land_ctl_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_tmp = read_var_multi_file(land_sen_files, 'AvgSurfT_tavg', loc_lat, loc_lon, lat_names, lon_names)
            Ctl_tmp       = Ctl_tmp -273.15
            Sen_tmp       = Sen_tmp -273.15
        elif var_name in ["Rnet",]:
            land_ctl_files= [land_ctl_path+'Lwnet_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files= [land_sen_path+'Lwnet_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_Lwnet_tmp = read_var_multi_file(land_ctl_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_Lwnet_tmp = read_var_multi_file(land_sen_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            land_ctl_files= [land_ctl_path+'Swnet_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files= [land_sen_path+'Swnet_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_Swnet_tmp = read_var_multi_file(land_ctl_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_Swnet_tmp = read_var_multi_file(land_sen_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            Ctl_tmp = Ctl_Lwnet_tmp+Ctl_Swnet_tmp
            Sen_tmp = Sen_Lwnet_tmp+Sen_Swnet_tmp
        elif var_name in ["EF",]:
            land_ctl_files      = [land_ctl_path+'Lwnet_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files      = [land_sen_path+'Lwnet_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_Lwnet_tmp = read_var_multi_file(land_ctl_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_Lwnet_tmp = read_var_multi_file(land_sen_files, 'Lwnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            land_ctl_files      = [land_ctl_path+'Swnet_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files      = [land_sen_path+'Swnet_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_Swnet_tmp = read_var_multi_file(land_ctl_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_Swnet_tmp = read_var_multi_file(land_sen_files, 'Swnet_tavg', loc_lat, loc_lon, lat_names, lon_names)
            land_ctl_files      = [land_ctl_path+'Qle_tavg/LIS.CABLE.201701-202002.nc']
            land_sen_files      = [land_sen_path+'Qle_tavg/LIS.CABLE.201701-202002.nc']
            time, Ctl_Qle_tmp   = read_var_multi_file(land_ctl_files, 'Qle_tavg', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_Qle_tmp   = read_var_multi_file(land_sen_files, 'Qle_tavg', loc_lat, loc_lon, lat_names, lon_names)
            ctl_Rnet            = Ctl_Lwnet_tmp + Ctl_Swnet_tmp
            sen_Rnet            = Sen_Lwnet_tmp + Sen_Swnet_tmp
            Ctl_tmp             = np.where(abs(ctl_Rnet)>0.01, Ctl_Qle_tmp/ctl_Rnet,np.nan)
            Sen_tmp             = np.where(abs(sen_Rnet)>0.01, Sen_Qle_tmp/sen_Rnet,np.nan)
        elif var_name in ["SM_top50cm",]:
            land_ctl_files = [land_ctl_path+'SoilMoist_inst/LIS.CABLE.201701-202002.nc']
            land_sen_files = [land_sen_path+'SoilMoist_inst/LIS.CABLE.201701-202002.nc']
            time, Ctl_temp = read_var_multi_file(land_ctl_files, 'SoilMoist_inst', loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_temp = read_var_multi_file(land_sen_files, 'SoilMoist_inst', loc_lat, loc_lon, lat_names, lon_names)
            # [.022, .058, .154, .409, 1.085, 2.872]
            Ctl_tmp    = Ctl_temp[:,0,:,:]*0.022 + Ctl_temp[:,1,:,:]*0.058 + Ctl_temp[:,2,:,:]*0.154 + Ctl_temp[:,3,:,:]*0.266
            Sen_tmp    = Sen_temp[:,0,:,:]*0.022 + Sen_temp[:,1,:,:]*0.058 + Sen_temp[:,2,:,:]*0.154 + Sen_temp[:,3,:,:]*0.266
        elif var_name in ['VPD','VPDmax','VPDmin']:
            tair_ctl_files= [land_ctl_path+'Tair_f_inst/LIS.CABLE.201701-202002.nc']
            tair_sen_files= [land_sen_path+'Tair_f_inst/LIS.CABLE.201701-202002.nc']
            qair_ctl_files= [land_ctl_path+'Qair_f_inst/LIS.CABLE.201701-202002.nc']
            qair_sen_files= [land_sen_path+'Qair_f_inst/LIS.CABLE.201701-202002.nc']
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

            # print("type(time_in)",type(time_in),"time=",time_in)
            # print("type(time_out)",type(time_out),"time_out=",time_out)

            f_ctl             = interp1d(np.array(time_in), Pres_ctl_tmp[:], kind='linear',fill_value='extrapolate', axis=0)
            f_sen             = interp1d(np.array(time_in), Pres_sen_tmp[:],kind='linear', fill_value='extrapolate', axis=0)
            Pres_ctl          = f_ctl(np.array(time_out))
            Pres_sen          = f_sen(np.array(time_out))
            Ctl_tmp           = qair_to_vpd(Qair_ctl, Tair_ctl, Pres_ctl)
            Sen_tmp           = qair_to_vpd(Qair_sen, Tair_sen, Pres_sen)
        else:
            land_ctl_files= [land_ctl_path+var_name+'/LIS.CABLE.201701-202002.nc']
            land_sen_files= [land_sen_path+var_name+'/LIS.CABLE.201701-202002.nc']
            time, Ctl_tmp = read_var_multi_file(land_ctl_files, var_name, loc_lat, loc_lon, lat_names, lon_names)
            time, Sen_tmp = read_var_multi_file(land_sen_files, var_name, loc_lat, loc_lon, lat_names, lon_names)

        wrf            = Dataset(wrf_path,  mode='r')
        lons           = wrf.variables['XLONG'][0,:,:]
        lats           = wrf.variables['XLAT'][0,:,:]

        for j in np.arange(6):

            if 'max' in var_name:
                # average of daily max
                ctl_in       = spital_var_max(time,Ctl_tmp,time_ss[j],time_es[j])
                sen_in       = spital_var_max(time,Sen_tmp,time_ss[j],time_es[j])
            elif 'min' in var_name:
                # average of daily min
                ctl_in       = spital_var_min(time,Ctl_tmp,time_ss[j],time_es[j])
                sen_in       = spital_var_min(time,Sen_tmp,time_ss[j],time_es[j])
            elif 'TDR' in var_name:
                # average of daily min
                ctl_in_max   = spital_var_max(time,Ctl_tmp,time_ss[j],time_es[j])
                sen_in_max   = spital_var_max(time,Sen_tmp,time_ss[j],time_es[j])
                ctl_in_min   = spital_var_min(time,Ctl_tmp,time_ss[j],time_es[j])
                sen_in_min   = spital_var_min(time,Sen_tmp,time_ss[j],time_es[j])
                ctl_in       = ctl_in_max - ctl_in_min
                sen_in       = sen_in_max - sen_in_min
            else:
                ctl_in       = spital_var(time,Ctl_tmp,time_ss[j],time_es[j])
                sen_in       = spital_var(time,Sen_tmp,time_ss[j],time_es[j])
            var_diff     = sen_in - ctl_in

            if var_name in ['WaterTableD_tavg','WatTable']:
                var_diff     = var_diff/1000.
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


        # read lat and lon outs

        # =============== CHANGE HERE ===============
        cmap  = plt.cm.BrBG
        clevs_percentage =  [-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,-4,-3,-2,-1,1,2,3,4,5,10,15,20,25,30,35,40,45,50]

        if var_name in ['WaterTableD_tavg','WatTable']:
            clevs = [-4,-3,-2,-1.5,-1,-0.5,0.5,1,1.5,2,3,4]
        elif var_name in ['EF']:
            clevs = [-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]
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
        # elif var_name in ["Qle_tavg","Qh_tavg","Qg_tavg"]:
        #     clevs = [-140,-120,-100,-80,-60,-40,-20,-10,-5,5,10,20,40,60,80,100,120,140]
        elif var_name in ["Swnet_tavg","Lwnet_tavg","SWdown_f_inst","LWdown_f_inst","Qle_tavg","Qh_tavg","Qg_tavg"]:
            clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
        elif var_name in ["Rnet"]:
            clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
            cmap  = plt.cm.BrBG_r
        elif var_name in ["Wind_f_inst",]:
            clevs = [-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4]
        elif var_name in ["Psurf_f_inst"]:
            clevs = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
        elif var_name in ["Tair_f_inst","Tmax","Tmin","VegT_tavg","VegTmax","VegTmin",
                          "AvgSurfT_tavg","SurfTmax","SurfTmin","SoilTemp_inst",'TDR','VegTDR','SurfTDR']:
            # clevs = [-2.,-1.8,-1.6,-1.4,-1.2,-1.,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.]
            clevs = [-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2]
            # clevs = [-3,-2.5,-2,-1.5,-1,-0.5,-0.1,0.1,0.5,1.,1.5,2,2.5,3.]
            cmap  = plt.cm.seismic
        elif var_name in ["Wind_f_inst",]:
            clevs = [-2.,-1.5,-1,-0.5,-0.1,0.1,0.5,1.,1.5,2.]
        elif var_name in ["FWsoil_tavg","SmLiqFrac_inst","SmFrozFrac_inst"]:
            clevs = [-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3,0.35]
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

                #hot_r # BrBG

                # start plotting
                if loc_lat == None:
                    ax.set_extent([135,155,-40,-25])
                else:
                    ax.set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

                ax.coastlines(resolution="50m",linewidth=1)

                # # Add gridlines
                gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
                gl.xlabels_top   = False
                gl.ylabels_right = False
                gl.xlines        = False
                gl.ylines        = False

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
                cb.ax.tick_params(labelsize=10, labelrotation=45)
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

                plt.savefig('./plots/spatial_map_'+message+'_layer='+str(j)+'.png',dpi=300)
                cb = None
                gl = None
                ax = None
                fig= None

        elif len(np.shape(var_diff)) ==2:
            # ================== Start Plotting =================
            # fig = plt.figure(figsize=(6,5))
            # ax = plt.axes(projection=ccrs.PlateCarree())

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
            # choose colormap
            # cmap  = plt.cm.seismic

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

            # print("any(not np.isnan(var_diff))",any(not np.isnan(var_diff)))
            plot1 = axs[0].contourf(lons, lats, var_diff, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
            cbar = plt.colorbar(plot1, ax=axs[0], ticklocation="right", pad=0.08, orientation="horizontal",
                    aspect=40, shrink=1) # cax=cax,
            cbar.ax.tick_params(labelsize=10,labelrotation=45)
            # plt.title(var_name, size=16)
            # if shape_path != None:
            #     # Requires the pyshp package
            #     sf = shp.Reader(shape_path)
            #
            #     for shape in sf.shapeRecords():
            #         x = [i[0] for i in shape.shape.points[:]]
            #         y = [i[1] for i in shape.shape.points[:]]
            #         plt.plot(x,y,c="black")

            rate = np.where( ctl_in != 0, var_diff/ctl_in, np.nan)
            plot2 = axs[1].contourf(lons, lats, rate*100., clevs_percentage, transform=ccrs.PlateCarree(), cmap=cmap, extend='both')
            cbar  = plt.colorbar(plot2, ax=axs[1], ticklocation="right", pad=0.08, orientation="horizontal",
                    aspect=40, shrink=1) # cax=cax,
            cbar.ax.tick_params(labelsize=10,labelrotation=45)
            # plt.title(var_name+"_diff_percentage", size=16)

            # if shape_path != None:
            #     # Requires the pyshp package
            #     sf = shp.Reader(shape_path)
            #
            #     for shape in sf.shapeRecords():
            #         x = [i[0] for i in shape.shape.points[:]]
            #         y = [i[1] for i in shape.shape.points[:]]
            #         plt.plot(x,y,c="black")

            plt.savefig('./plots/spatial_map_'+message + "_" + var_name+'.png',dpi=300)


if __name__ == "__main__":


    # ======================= Option =======================
    # 2017-2019 drought polygon shapefile
    shape_path = "/g/data/w97/ad9701/drought_2017to2020/drought_focusArea/smooth_polygon_drought_focusArea.shp"

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

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"
        land_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/LIS_output/"
        land_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/LIS_output/"
        atmo_sen_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_sen+"/WRF_output/"
        atmo_ctl_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_ctl+"/WRF_output/"

        if 1:
            '''
            Difference plot yearly
            '''
            var_names  = [  "TDR"
                            #"FWsoil_tavg","EF"
                            #"FWsoil_tavg",
                            #"SM_top50cm",
                            #"VPDmax","VPDmin"
                            # "LAI_inst","Albedo_inst",
                            # "VPD",
                            # "GPP_tavg",#"NPP_tavg",
                            # "SurfTmax","SurfTmin",
                            # "Tmax","Tmin",
                            # "VegT_tavg","Tair_f_inst","AvgSurfT_tavg",
                            # # "VegTmax","VegTmin",
                            ]

            time_ss    = [datetime(2017,6,1,0,0,0,0), datetime(2017,12,1,0,0,0,0),
                          datetime(2018,6,1,0,0,0,0), datetime(2018,12,1,0,0,0,0),
                          datetime(2019,6,1,0,0,0,0), datetime(2019,12,1,0,0,0,0)]
            time_ss    = [datetime(2017,9,1,0,0,0,0), datetime(2018,3,1,0,0,0,0),
                          datetime(2018,9,1,0,0,0,0), datetime(2019,3,1,0,0,0,0),
                          datetime(2019,9,1,0,0,0,0), datetime(2020,3,1,0,0,0,0)]

            message    = "Winter_Summer"
            spatial_map_winter_summer(land_ctl_path, land_sen_path, var_names, time_ss=time_ss, time_es=time_es, lat_names="lat",
                                lon_names="lon",loc_lat=loc_lat, loc_lon=loc_lon, wrf_path=wrf_path, shape_path=shape_path,
                                message=message)
