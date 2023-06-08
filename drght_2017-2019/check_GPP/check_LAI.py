#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"


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
import pandas as pd
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature, OCEAN
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from common_utils import *

def plot_time_series(file_paths_ctl, file_paths_sen, file_paths_sen_2=None,var_name=None,
                          time_s=None, time_e=None, loc_lat=None, loc_lon=None,
                          lat_name=None, lon_name=None, message=None, multi=None,iveg=None):


    # read data
    time_ctl, Var_ctl = read_var_multi_file(file_paths_ctl, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time_sen, Var_sen = read_var_multi_file(file_paths_sen, var_name, loc_lat, loc_lon, lat_name, lon_name)
    print("(~np.isnan(Var_ctl[0,:,:])).sum()",(~np.isnan(Var_ctl[0,:,:])).sum())
    Var_ctl           = np.nanmean(Var_ctl,axis=(1,2))
    Var_sen           = np.nanmean(Var_sen,axis=(1,2))
    print(message," Var_ctl",Var_ctl)
    print("Var_sen",Var_sen)

    if file_paths_sen_2 != None:
        if message == "ALB_CTL_AU-Cum":
            file_sen_2    = Dataset(file_paths_sen_2[0],  mode='r')
            LAI_sen_2     = file_sen_2.variables['LAI'][4:6]
        elif message == "ALB_CTL_AU-Tum":
            file_sen_2    = Dataset(file_paths_sen_2[0],  mode='r')
            LAI_sen_2     = file_sen_2.variables['LAI'][15]
        print("LAI_sen_2",LAI_sen_2)

if __name__ == "__main__":

    # ======================= Option =======================
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


    # small box
    # loc_lat    = [-33,-29]
    # loc_lon    = [147,149]

    # # east coast
    # loc_lat    = [-33,-27]
    # loc_lon    = [152,154]
    # PFT        = False


    ####################################
    #         plot_time_series         #
    ####################################
    if 1:
        lat_name       = "lat"
        lon_name       = "lon"
        iveg           = None

        case_name_ctl  = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
        case_name_sen  = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_name_ctl+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"

        LIS_path_ctl   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_name_ctl+"/LIS_output/"
        LIS_path_sen   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/"+case_name_sen+"/LIS_output/"


        var_name       = "LAI_inst"
        file_paths_ctl = [ LIS_path_ctl+var_name+'/LIS.CABLE.201701-201912_yearmean.nc' ]
        file_paths_sen = [ LIS_path_sen+var_name+'/LIS.CABLE.201701-201912_yearmean.nc' ]


        # AU-Cum
        message = "ALB_CTL_AU-Cum"
        loc_lat = [-33.6133-0.00000001,-33.6133+0.00000001] #000000
        loc_lon = [150.7225-0.00000001,150.7225+0.00000001]
        time_s  = datetime(2017,1,1,0,0,0,0)
        time_e  = datetime(2019,1,1,0,0,0,0)
        file_paths_sen_2 = ['/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/check_GPP/AU-Cum_2013-2018_OzFlux_Met_yearmean.nc']
        plot_time_series(file_paths_ctl,file_paths_sen,file_paths_sen_2, var_name,
                            time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message,
                            multi=True, iveg=iveg)


        # AU-Tum
        message = "ALB_CTL_AU-Tum"
        loc_lat = [-35.6566-0.00000001,-35.6566+0.00000001]
        loc_lon = [148.1517-0.00000001,148.1517+0.00000001]
        time_s  = datetime(2017,1,1,0,0,0,0)
        time_e  = datetime(2018,1,1,0,0,0,0)
        file_paths_sen_2 = ['/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/check_GPP/AU-Tum_2002-2017_OzFlux_Met_yearmean.nc']
        plot_time_series(file_paths_ctl,file_paths_sen,file_paths_sen_2, var_name,
                            time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message,
                            multi=True, iveg=iveg)
