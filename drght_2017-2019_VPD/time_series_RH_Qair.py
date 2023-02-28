#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"


'''
Functions:
1. Compare LIS-CABLE with GRACE, GLEAM, & DOLCE
2. GW vs FD
3. plot time-series and spitial (difference) map
'''

'''
History:
'''

from netCDF4 import Dataset
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
from scipy.stats import linregress
import matplotlib.ticker as mticker
from common_utils import *


def plot_time_series(file_paths_ctl, file_paths_sen, file_paths_sen_2=None,var_names=None,
                          time_s=None, time_e=None, loc_lat=None, loc_lon=None,
                          lat_name=None, lon_name=None, message=None):

    print("======== In plot_time_series =========")
    time_ctl, Var_ctl = read_var_multi_file(file_paths_ctl, var_names[0], loc_lat, loc_lon, lat_name, lon_name)
    if file_paths_sen != None:
        time_sen, Var_sen = read_var_multi_file(file_paths_sen, var_names[1], loc_lat, loc_lon, lat_name, lon_name)
    if file_paths_sen_2 != None:
        time_sen_2, Var_sen_2 = read_var_multi_file(file_paths_sen_2, var_names[2], loc_lat, loc_lon, lat_name, lon_name)
    print("time_ctl ",time_ctl)
    print("time_sen ",time_sen)
    print("time_sen_2 ",time_sen_2)
    
    var_ctl   = np.nanmean(Var_ctl[:,:,:],axis=(1,2))
    var_sen   = np.nanmean(Var_sen[:,:,:],axis=(1,2))
    var_sen_2 = np.nanmean(Var_sen_2[:,:,:],axis=(1,2))

    time_ctl_day   = []
    time_sen_day   = []
    time_sen_2_day = []
    
    for t in time_ctl:
        time_ctl_day.append(t.days) 
    for t in time_sen:
        time_sen_day.append(t.days) 
    for t in time_sen_2:
        time_sen_2_day.append(t.days) 
    
    cleaner_dates = ["2017","2018", "2019", "2020" ]
    xtickslocs    = [0,     365,      730,   1095  ]

    fig, ax       = plt.subplots()
    df_ctl        = pd.DataFrame({'ctl': var_ctl})
    df_sen        = pd.DataFrame({'sen': var_sen})
    
    ctl_smooth    = df_ctl['ctl'].rolling(window=30).mean().values
    sen_smooth    = df_sen['sen'].rolling(window=30).mean().values
    
    print("ctl_smooth", ctl_smooth)
    print("sen_smooth", sen_smooth)
    
    mask_ctl      = np.array(~np.isnan(ctl_smooth))
    mask_sen      = np.array(~np.isnan(sen_smooth))
    mask_sen2     = np.array(~np.isnan(var_sen_2))
    
    print("mask_ctl", mask_ctl)
    print("mask_sen", mask_sen)
    print("mask_sen2", mask_sen2)
    
    print("ctl_smooth[mask_ctl]", ctl_smooth[mask_ctl])
    
    result_ctl    = linregress(time_ctl_day[mask_ctl==True],ctl_smooth[mask_ctl==True])
    result_sen    = linregress(time_sen_day[mask_sen==True],sen_smooth[mask_sen==True])
    result_sen2   = linregress(time_sen_2_day[mask_sen2==True],var_sen_2[mask_sen2==True])
    
    print("result_ctl", result_ctl)
    print("result_ctl.intercept", result_ctl.intercept)
    print("result_ctl.intercept", result_ctl.slope)
    
    ax.plot(time_ctl_day, df_ctl['ctl'].rolling(window=30).mean(), c = "red",  label=var_names[0], alpha=0.5)
    ax.plot(time_sen_day, df_sen['sen'].rolling(window=30).mean(), c = "blue", label=var_names[1], alpha=0.5)
    ax.plot(time_sen_2_day, var_sen_2, c = "green", label=var_names[2], alpha=0.5)
    
    pre_ctl  = result_ctl.intercept + result_ctl.slope*time_ctl_day
    pre_sen  = result_sen.intercept + result_sen.slope*time_sen_day
    pre_sen2 = result_sen2.intercept + result_sen2.slope*time_sen_2_day
    
    ax.plot(time_ctl_day,  pre_ctl, c = "red",  label=var_names[0], alpha=1)
    ax.plot(time_sen_day,  pre_sen, c = "blue", label=var_names[1], alpha=1)
    ax.plot(time_sen_2_day,pre_sen2,c = "green",label=var_names[2], alpha=1)

    ax.legend()
    fig.tight_layout()

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
    
    if 1:
        message               = "VPH_vs"
        lat_name              = "latitude"
        lon_name              = "longitude"

        time_s                = datetime(1970,1,1,0,0,0,0)
        time_e                = datetime(2020,1,1,0,0,0,0)

        file_path_AWAP_vph09  = ["/g/data/w97/Shared_data/Observations/AWAP_all_variables/daily/vph09/AWAP_daily_vph09_1970_2019.nc"]
        file_path_AWAP_vph15  = ["/g/data/w97/Shared_data/Observations/AWAP_all_variables/daily/vph15/AWAP_daily_vph15_1970_2019.nc"]
        file_path_HadISDH_vph = ["/g/data/w97/mm3972/data/HadISDH/HadISDH.lande.4.4.0.2021f_FLATgridHOM5by5_anoms9120.nc"]

        var_names             = ["vph09","vph15", "e_abs"]

        plot_time_series(file_path_AWAP_vph09, file_path_AWAP_vph15, file_path_HadISDH_vph, var_names,
                         time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                         lat_name=lat_name, lon_name=lon_name, message=message)

