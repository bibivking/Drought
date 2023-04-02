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
    
    if 'AWAP_3h_v1' in file_paths_ctl[0]:
        time_ctl, Var_ctl    = read_var_multi_file(file_paths_ctl, var_names[0], loc_lat, loc_lon, "lat", "lon")
        Var_ctl,time_ctl_day = time_clip_to_day_dates(time_ctl, Var_ctl, time_s, time_e, seconds=None)
        if "Qair" in file_paths_ctl[0]:
            Var_ctl              = Var_ctl*1000.
    else:
        time_ctl, Var_ctl    = read_var_multi_file(file_paths_ctl, var_names[0], loc_lat, loc_lon, lat_name, lon_name)
        time_ctl_day         = []
        for t in time_ctl:
            time_ctl_day.append(t.days) 
        time_ctl_day  = np.array(time_ctl_day)
    var_ctl           = np.nanmean(Var_ctl[:,:,:],axis=(1,2))

    if file_paths_sen != None:
        time_sen, Var_sen = read_var_multi_file(file_paths_sen, var_names[1], loc_lat, loc_lon, lat_name, lon_name)
        var_sen   = np.nanmean(Var_sen[:,:,:],axis=(1,2))
        time_sen_day   = []
        for t in time_sen:
            time_sen_day.append(t.days) 
        time_sen_day = np.array(time_sen_day)

    if file_paths_sen_2 != None:
        time_sen_2, Var_sen_2 = read_var_multi_file(file_paths_sen_2, var_names[2], loc_lat, loc_lon, lat_name, lon_name)
        var_sen_2 = np.nanmean(Var_sen_2[:,:,:],axis=(1,2))
        time_sen_2_day = []
        for t in time_sen_2:
            time_sen_2_day.append(t.days) 
        time_sen_2_day = np.array(time_sen_2_day)

    # ============== Plotting ==============
    fig, ax       = plt.subplots()
    
    df_ctl        = pd.DataFrame({'ctl': var_ctl})
    ctl_smooth    = df_ctl['ctl'].rolling(window=30).mean().values
    mask_ctl      = np.array(~np.isnan(ctl_smooth))
    result_ctl    = linregress(time_ctl_day[mask_ctl],ctl_smooth[mask_ctl])
    
    ax.plot(time_ctl_day, df_ctl['ctl'].rolling(window=30).mean(), c = "red",  label=var_names[0], alpha=0.5)
    pre_ctl  = result_ctl.intercept + result_ctl.slope*time_ctl_day
    ax.plot(time_ctl_day,  pre_ctl, c = "red",  label=var_names[0], alpha=1)
    
    if file_paths_sen != None:
        df_sen        = pd.DataFrame({'sen': var_sen})
        sen_smooth    = df_sen['sen'].rolling(window=30).mean().values
        mask_sen      = np.array(~np.isnan(sen_smooth))
        result_sen    = linregress(time_sen_day[mask_sen==True],sen_smooth[mask_sen])
        
        ax.plot(time_sen_day, df_sen['sen'].rolling(window=30).mean(), c = "blue", label=var_names[1], alpha=0.5)
        pre_sen  = result_sen.intercept + result_sen.slope*time_sen_day
        ax.plot(time_sen_day,  pre_sen, c = "blue", label=var_names[1], alpha=1)
        
    if file_paths_sen_2 != None:
        mask_sen2     = np.array(~np.isnan(var_sen_2))
        result_sen2   = linregress(time_sen_2_day[mask_sen2==True],var_sen_2[mask_sen2==True])
        
        #ax.plot(time_sen_2_day, var_sen_2, c = "green", label=var_names[2], alpha=0.5)
        pre_sen2 = result_sen2.intercept + result_sen2.slope*time_sen_2_day
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
        
        lat_name              = "latitude"
        lon_name              = "longitude"

        time_s                = datetime(1970,1,1,0,0,0,0)
        time_e                = datetime(2020,1,1,0,0,0,0)
        file_path_AWAP        = "/g/data/w97/Shared_data/AWAP_3h_v1/"
        file_path_AWAP_detrend= "/g/data/w97/mm3972/data/AWAP_detrend/"

        file_path_AWAP_vph09  = ["/g/data/w97/Shared_data/Observations/AWAP_all_variables/daily/vph09/AWAP_daily_vph09_1970_2019.nc"]
        file_path_AWAP_vph15  = ["/g/data/w97/Shared_data/Observations/AWAP_all_variables/daily/vph15/AWAP_daily_vph15_1970_2019.nc"]
        file_path_HadISDH_vph = ["/g/data/w97/mm3972/data/HadISDH/HadISDH.lande.4.4.0.2021f_FLATgridHOM5by5_anoms9120.nc"]
        file_path_HadISDH_q   = ["/g/data/w97/mm3972/data/HadISDH/HadISDH.landq.4.4.0.2021f_FLATgridHOM5by5_anoms9120.nc"]
        file_path_HadISDH_T   = ["/g/data/w97/mm3972/data/HadISDH/HadISDH.landT.4.4.0.2021f_FLATgridHOM5by5_anoms9120.nc"]
        file_path_AWAP_tmax   = ["/g/data/w97/Shared_data/Observations/AWAP_all_variables/daily/tmax/AWAP_daily_tmax_1970_2019.nc"]
        file_path_AWAP_tmin   = ["/g/data/w97/Shared_data/Observations/AWAP_all_variables/daily/tmin/AWAP_daily_tmin_1970_2019.nc"]

        file_path_detrend_tmax_daily   = ["/g/data/w97/mm3972/data/AWAP_detrend/tmax/AWAP_daily_tmax_1970_2019_daily.nc"]
        file_path_detrend_tmax_monthly = ["/g/data/w97/mm3972/data/AWAP_detrend/tmax/AWAP_daily_tmax_1970_2019_monthly.nc"]
        file_path_detrend_tmin         = ["/g/data/w97/mm3972/data/AWAP_detrend/tmin/AWAP_daily_tmin_1970_2019.nc"]
        file_path_detrend_vph09        = ["/g/data/w97/mm3972/data/AWAP_detrend/vph09_detrend/AWAP_daily_vph09_1970_2019.nc"]
        file_path_detrend_vph15        = ["/g/data/w97/mm3972/data/AWAP_detrend/vph15_detrend/AWAP_daily_vph15_1970_2019.nc"]
        
        if 0:
            file_path_AWAP_q           = []
            file_path_AWAP_T           = []
            file_path_AWAP_q_detrend   = []
            file_path_AWAP_T_detrend   = []
            for y in np.arange(1970,2020,1):
                file_path_AWAP_q.append(file_path_AWAP + "Qair/AWAP.Qair.3hr."+str(y)+".nc")
                file_path_AWAP_T.append(file_path_AWAP + "Tair/AWAP.Tair.3hr."+str(y)+".nc")
                file_path_AWAP_q.append(file_path_AWAP_detrend + "Qair/AWAP.Qair.3hr."+str(y)+".nc")
                file_path_AWAP_T.append(file_path_AWAP_detrend + "Tair/AWAP.Tair.3hr."+str(y)+".nc")
                
            message               = "Qair_vs_detrend"
            var_names             = ["Qair","Qair"]

            plot_time_series(file_path_AWAP_q, file_path_AWAP_q_detrend, None, var_names,
                             time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                             lat_name=lat_name, lon_name=lon_name, message=message)


            message               = "Tair_vs_detrend"
            var_names             = ["Tair","Tair"]

            plot_time_series(file_path_AWAP_T, file_path_AWAP_T_detrend, None, var_names,
                             time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                             lat_name=lat_name, lon_name=lon_name, message=message)
        
                
        if 0:
            message               = "tmin_vs_detrend_daily"
            var_names             = ["tmin","tmin"]

            #plot_time_series(file_path_AWAP_tmin, file_path_detrend_tmin, None, var_names,
            #                time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
            #                lat_name=lat_name, lon_name=lon_name, message=message)

            message               = "vph09_vs_detrend_daily"
            var_names             = ["vph09","vph09"]

            plot_time_series(file_path_AWAP_vph09, file_path_detrend_vph09, None, var_names,
                            time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message)

            message               = "vph15_vs_detrend_daily"
            var_names             = ["vph15","vph15"]

            plot_time_series(file_path_AWAP_vph15, file_path_detrend_vph15, None, var_names,
                            time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message)




        if 1:
            file_path_out_detrend = ["/g/data/w97/mm3972/model/cable/runs/VPD_drought/100th_check/outputs/cable_out_2017.nc"]
            file_path_out_default = ["/g/data/w97/mm3972/model/cable/runs/VPD_drought/detrended_Tair_VPD/outputs/cable_out_2017.nc"]
                
            message               = "detrend_vs_100th_check"
            var_names             = ["Qle","Qle"]

            plot_time_series(file_path_out_default, file_path_out_detrend, None, var_names,
                             time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                             lat_name=lat_name, lon_name=lon_name, message=message)
