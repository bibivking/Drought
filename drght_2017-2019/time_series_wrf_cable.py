#!/usr/bin/python

'''
Functions:
1. Compare LIS-CABLE with GRACE, GLEAM, & DOLCE
2. GW vs FD
3. plot time-series and spitial (difference) map
'''

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from common_utils import *

def plot_time_series(file_paths, var_name, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None, message=None):

    print("======== In plot_time_series =========")
    time,Var = read_var_multi_file(file_paths, var_name, loc_lat, loc_lon, lat_name, lon_name)
    
    if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg"]:          
        var = var*3600.

    if len(np.shape(Var)) == 3:
        
        print(var_name, "dimension=",len(np.shape(Var)))
        
        var = np.nanmean(Var,axis=(1,2))
        fig, ax = plt.subplots()
        ax.plot(np.arange(len(var)), var, c = "blue", label=var_name, alpha=0.5)
        # ax.set_title(var_name)
        ax.legend()
        fig.tight_layout()
        if message == None:
            message = var_name
        else:
            message = message + "_" + var_name

        plt.savefig('./plots/time_series_'+message+'.png',dpi=300)
        
    elif len(np.shape(Var)) == 4:
        
        print(var_name, "dimension=",len(np.shape(Var)))
        
        var = np.nanmean(Var,axis=(2,3))
    
        fig, ax = plt.subplots()
        ax.plot(np.arange(len(var[:,0])), var, label="lyr="+str(np.arange(1,7,1)), alpha=0.5) #c = "blue",
        # ax.set_title(var_name)
        ax.legend()
        fig.tight_layout()
        if message == None:
            message = var_name
        else:
            message = message + "_" + var_name

        plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

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

    ####################################
    #         plot_time_series         #
    ####################################

    case_name       = "drght_2017_2019_bl_pbl2_mp4_sf_sfclay2"
                     #"drght_2017_2019_bl_pbl5_mp6_sf_sfclay1" 
                     #"drght_2017_2019_bl_pbl2_mp4_sf_sfclay2" #"spinup_30yrx3_rerun" #
    LIS_path        = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"
    # SP4_off_path = "/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/"+case_name+"/outputs_4th/"

    file_paths = [ LIS_path+"LIS.CABLE.201701-201701.d01.nc",
                   LIS_path+"LIS.CABLE.201702-201702.d01.nc",
                   LIS_path+"LIS.CABLE.201703-201703.d01.nc",
                   LIS_path+"LIS.CABLE.201704-201704.d01.nc",
                   LIS_path+"LIS.CABLE.201705-201705.d01.nc",
                   LIS_path+"LIS.CABLE.201706-201706.d01.nc",
                   LIS_path+"LIS.CABLE.201707-201707.d01.nc",
                   LIS_path+"LIS.CABLE.201708-201708.d01.nc",
                   LIS_path+"LIS.CABLE.201709-201709.d01.nc",
                   LIS_path+"LIS.CABLE.201710-201710.d01.nc",
                   LIS_path+"LIS.CABLE.201711-201711.d01.nc",
                   LIS_path+"LIS.CABLE.201712-201712.d01.nc",
                   LIS_path+"LIS.CABLE.201801-201801.d01.nc",
                   LIS_path+"LIS.CABLE.201802-201802.d01.nc",
                   LIS_path+"LIS.CABLE.201803-201803.d01.nc",
                   LIS_path+"LIS.CABLE.201804-201804.d01.nc",
                   LIS_path+"LIS.CABLE.201805-201805.d01.nc",
                   LIS_path+"LIS.CABLE.201806-201806.d01.nc",
                   LIS_path+"LIS.CABLE.201807-201807.d01.nc",
                   LIS_path+"LIS.CABLE.201808-201808.d01.nc",
                   LIS_path+"LIS.CABLE.201809-201809.d01.nc",
                   LIS_path+"LIS.CABLE.201810-201810.d01.nc",
                   LIS_path+"LIS.CABLE.201811-201811.d01.nc",
                   LIS_path+"LIS.CABLE.201812-201812.d01.nc"]

    print(file_paths)

    message   = case_name
    var_names = ["SoilMoist_inst"]
    # var_names = ["GWwb_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
    #               "WaterTableD_tavg","Qle_tavg","Qh_tavg","Qg_tavg","AvgSurfT_tavg","VegT_tavg","FWsoil_tavg"
    #               "SoilMoist_inst"]
    lat_name  = "lat"
    lon_name  = "lon"
    for var_name in var_names:
        plot_time_series(file_paths, var_name, loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message)
