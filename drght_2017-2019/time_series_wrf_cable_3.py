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
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
from common_utils import *

def plot_time_series(file_paths, var_name, time_s=None, time_e=None, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None, message=None):

    print("======== In plot_time_series =========")
    time, Var = read_var_multi_file(file_paths, var_name, loc_lat, loc_lon, lat_name, lon_name)

    if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg"]:
        Var_daily = time_clip_to_day_sum(time, Var, time_s, time_e, seconds=None)
        Var_daily = Var_daily*3600.
    elif var_name in ["WaterTableD_tavg"]:
        Var_daily = time_clip_to_day(time, Var, time_s, time_e, seconds=None)
        Var_daily = Var_daily/1000.
    else:
        Var_daily = time_clip_to_day(time, Var, time_s, time_e, seconds=None)

    colors        = [ "forestgreen", "yellowgreen","orange","red","black",]
                    # [ "black", "grey","lightcoral","red","orange","gold",
                    #  "yellow","yellowgreen","forestgreen","aquamarine","skyblue",
                    #  "blue","blueviolet","violet","purple","crimson","pink"]
    cleaner_dates = ["2017","2018", "2019", "2020" ]
    xtickslocs    = [0,     365,      730,   1095  ]

    if len(np.shape(Var_daily)) == 3:

        print(var_name, "dimension=",np.shape(Var_daily))

        var = np.nanmean(Var_daily[1:,:,:],axis=(1,2))
        fig, ax = plt.subplots(figsize=[9,9])
        ax.plot(np.arange(len(var)), var, c = "blue", label=var_name, alpha=0.5)
        # ax.set_title(var_name)
        # ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
        ax.legend()
        fig.tight_layout()
        if message == None:
            message = var_name
        else:
            message = message + "_" + var_name

        plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

    elif len(np.shape(Var_daily)) == 4:

        print(var_name, "dimension=",np.shape(Var_daily))

        labels = ["lyr1","lyr2","lyr3","lyr4","lyr5","lyr6"]
        var = np.nanmean(Var_daily[1:,:,:,:],axis=(2,3))

        fig, ax = plt.subplots(figsize=[9,9])
        ax.plot(np.arange(len(var[:,0])), var, label=labels, alpha=0.5) #c = "blue",

        # ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
        ax.legend()
        fig.tight_layout()
        if message == None:
            message = var_name
        else:
            message = message + "_" + var_name

        plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

def plot_time_series_diff(file_paths_ctl, file_paths_sen, file_paths_sen_2=None,var_name=None,
                          time_s=None, time_e=None, loc_lat=None, loc_lon=None,
                          lat_name=None, lon_name=None, message=None, multi=None,iveg=None):

    print("======== In plot_time_series =========")
    time_ctl, Var_ctl = read_var_multi_file(file_paths_ctl, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time_sen, Var_sen = read_var_multi_file(file_paths_sen, var_name, loc_lat, loc_lon, lat_name, lon_name)
    if file_paths_sen_2 != None:
        time_sen_2, Var_sen_2 = read_var_multi_file(file_paths_sen_2, var_name, loc_lat, loc_lon, lat_name, lon_name)

    if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg"]:
        Var_daily_ctl = time_clip_to_day_sum(time_ctl, Var_ctl, time_s, time_e, seconds=None)
        Var_daily_ctl = Var_daily_ctl*3600.
    elif var_name in ["WaterTableD_tavg"]:
        Var_daily_ctl = time_clip_to_day(time_ctl, Var_ctl, time_s, time_e, seconds=None)
        Var_daily_ctl = Var_daily_ctl/1000.
    else:
        Var_daily_ctl = time_clip_to_day(time_ctl, Var_ctl, time_s, time_e, seconds=None)

    if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg"]:
        Var_daily_sen = time_clip_to_day_sum(time_sen, Var_sen, time_s, time_e, seconds=None)
        Var_daily_sen = Var_daily_sen*3600.
    elif var_name in ["WaterTableD_tavg"]:
        Var_daily_sen = time_clip_to_day(time_sen, Var_sen, time_s, time_e, seconds=None)
        Var_daily_sen = Var_daily_sen/1000.
    else:
        Var_daily_sen = time_clip_to_day(time_sen, Var_sen, time_s, time_e, seconds=None)

    if file_paths_sen_2 != None:
        if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg"]:
            Var_daily_sen_2 = time_clip_to_day_sum(time_sen_2, Var_sen_2, time_s, time_e, seconds=None)
            Var_daily_sen_2 = Var_daily_sen_2*3600.
        elif var_name in ["WaterTableD_tavg"]:
            Var_daily_sen_2 = time_clip_to_day(time_sen_2, Var_sen_2, time_s, time_e, seconds=None)
            Var_daily_sen_2 = Var_daily_sen_2/1000.
        else:
            Var_daily_sen_2 = time_clip_to_day(time_sen_2, Var_sen_2, time_s, time_e, seconds=None)

    if multi == None:
        Var_daily_diff = Var_daily_sen - Var_daily_ctl

    if iveg != None:
        LC_file        = ["/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/LIS_output/LIS.CABLE.201701-202006_ALB_LAI.nc"]
        time_lc, LC    = read_var_multi_file(LC_file, "Landcover_inst", loc_lat, loc_lon, lat_name, lon_name)
        landcover      = time_clip_to_day(time_lc, LC, time_s, time_e, seconds=None)

    colors        = [ "forestgreen", "yellowgreen","orange","red","black",]
                    # [ "black", "grey","lightcoral","red","orange","gold",
                    #  "yellow","yellowgreen","forestgreen","aquamarine","skyblue",
                    #  "blue","blueviolet","violet","purple","crimson","pink"]
    cleaner_dates = ["2017","2018", "2019", "2020" ]
    xtickslocs    = [0,     365,      730,   1095  ]

    if len(np.shape(Var_daily_ctl)) == 3:

        print(var_name, "dimension=",np.shape(Var_daily_ctl))
        fig, ax = plt.subplots(figsize=[9,9])

        if multi == None:
            if iveg == None:
                var = np.nanmean(Var_daily_diff[:,:,:],axis=(1,2)) # 1
                ax.plot(np.arange(len(var)), var, c = "blue", label="Δ"+var_name, alpha=0.5)
            elif np.shape(iveg)[0]>1:
                for i,pft in enumerate(iveg):
                    var = np.nanmean(np.where(landcover == pft, Var_daily_diff, np.nan),axis=(1,2))
                    ax.plot(np.arange(len(var)), var, c = colors[i], label="Δ"+var_name+" PFT="+str(pft), alpha=0.5)
            else:
                var = np.nanmean(np.where(landcover == iveg, Var_daily_diff, np.nan),axis=(1,2))
                ax.plot(np.arange(len(var)), var, c = "blue", label="Δ"+var_name, alpha=0.5)

        if multi == True:
            if iveg == None:
                var_ctl = np.nanmean(Var_daily_ctl[:,:,:],axis=(1,2))
                var_sen = np.nanmean(Var_daily_sen[:,:,:],axis=(1,2))
                ax.plot(np.arange(len(var_ctl)), var_ctl, c = "red",  label=var_name+"_ctl", alpha=0.5)
                ax.plot(np.arange(len(var_sen)), var_sen, c = "blue", label=var_name+"_sen", alpha=0.5)
            elif np.shape(iveg)[0]>1:
                for i,pft in enumerate(iveg):
                    var_ctl = np.nanmean(np.where(landcover == pft, Var_daily_ctl, np.nan),axis=(1,2))
                    var_sen = np.nanmean(np.where(landcover == pft, Var_daily_sen, np.nan),axis=(1,2))
                    df_ctl  = pd.DataFrame({'ctl': var_ctl})
                    df_sen  = pd.DataFrame({'sen': var_sen})
                    ax.plot( df_ctl['ctl'].rolling(window=30).mean(), c = colors[i],  label=var_name+"_ctl PFT="+str(pft), alpha=0.8)
                    ax.plot( df_sen['sen'].rolling(window=30).mean(), c = colors[i], label=var_name+"_sen PFT="+str(pft), alpha=0.5)
            else:
                var_ctl = np.nanmean(np.where(landcover == iveg, Var_daily_ctl, np.nan),axis=(1,2))
                var_sen = np.nanmean(np.where(landcover == iveg, Var_daily_sen, np.nan),axis=(1,2))
                ax.plot(np.arange(len(var_ctl)), var_ctl, c = "red",  label=var_name+"_ctl", alpha=0.5)
                ax.plot(np.arange(len(var_sen)), var_sen, c = "blue", label=var_name+"_sen", alpha=0.5)

            if file_paths_sen_2 != None:
                if iveg == None:
                    var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:],axis=(1,2))
                    ax.plot(np.arange(len(var_sen_2)), var_sen_2, c = "green", label=var_name+"_sen_2", alpha=0.5)
                elif np.shape(iveg)[0]>1:
                    for i,pft in enumerate(iveg):
                        var_sen_2 = np.nanmean(np.where(landcover == pft, Var_daily_sen_2, np.nan),axis=(1,2))
                        ax.plot(np.arange(len(var_sen_2)), var_sen_2, ls="-.", c = colors[i], label=var_name+"_sen_2 PFT="+str(pft), alpha=0.5)
                else:
                    var_sen_2 = np.nanmean(np.where(landcover == iveg, Var_daily_sen_2, np.nan),axis=(1,2))
                    ax.plot(np.arange(len(var_sen_2)), var_sen_2, c = "green", label=var_name+"_sen_2", alpha=0.5)

        # ax.set_title(var_name)
        # ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
        ax.legend()
        fig.tight_layout()
        if message == None:
            message = var_name
        else:
            message = message + "_" + var_name
        if multi == True:
            message = message + "_multi"
        if iveg != None and np.shape(iveg)[0] >1:
            message = message + "_iveg="+str(iveg)
        elif iveg != None:
            message = message + "_iveg="+str(iveg[0])+"-"+str(iveg[-1])

        plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

    elif len(np.shape(Var_daily_ctl)) == 4:
        # note that 4-D var doesn't support PFT lines
        print(var_name, "dimension=",np.shape(Var_daily_ctl))

        fig, ax = plt.subplots(figsize=[9,9])
        labels = ["lyr1","lyr2","lyr3","lyr4","lyr5","lyr6"]

        if multi == None:
            if iveg == None:
                var = np.nanmean(Var_daily_diff[:,:,:,:],axis=(2,3))
                ax.plot(np.arange(len(var[:,0])), var, label=labels, alpha=0.5) #c = "blue",
            else:
                var = np.nanmean(Var_daily_diff[:,:,:,:],axis=(2,3))
                for i in np.arange(6):
                    var[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_diff[:,i,:,:], np.nan),axis=(1,2))
        if multi == True:
            if iveg == None:
                var_ctl = np.nanmean(Var_daily_ctl[:,:,:,:],axis=(2,3))
                var_sen = np.nanmean(Var_daily_sen[:,:,:,:],axis=(2,3))
                ax.plot(np.arange(len(var_ctl[:,0])), var_ctl, ls = "-",  label=labels+"_ctl", alpha=0.5) #c = "blue",
                ax.plot(np.arange(len(var_sen[:,0])), var_sen, ls = "-.", label=labels+"_sen", alpha=0.5) #c = "blue",
            else:
                var_ctl = np.nanmean(Var_daily_ctl[:,:,:,:],axis=(2,3))
                var_sen = np.nanmean(Var_daily_sen[:,:,:,:],axis=(2,3))
                for i in np.arange(6):
                    var_ctl[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_ctl[:,i,:,:], np.nan),axis=(1,2))
                    var_sen[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_sen[:,i,:,:], np.nan),axis=(1,2))
                ax.plot(np.arange(len(var_ctl[:,0])), var_ctl, ls = "-",  label=labels+"_ctl", alpha=0.5) #c = "blue",
                ax.plot(np.arange(len(var_sen[:,0])), var_sen, ls = "-.", label=labels+"_sen", alpha=0.5) #c = "blue",

            if file_paths_sen_2 != None:
                if iveg == None:
                    var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:,:],axis=(2,3))
                    ax.plot(np.arange(len(var_sen_2[:,0])), var_sen_2, ls = "--", label=var_name+"_sen_2", alpha=0.5)
                else:
                    var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:,:],axis=(2,3))
                    for i in np.arange(6):
                        var_sen_2[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_sen_2[:,i,:,:], np.nan),axis=(1,2))
                    ax.plot(np.arange(len(var_sen_2[:,0])), var_sen_2, ls = "-.", label=labels+"_sen", alpha=0.5) #c = "blue",

        # ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
        ax.legend()
        fig.tight_layout()

        if message == None:
            message = var_name
        else:
            message = message + "_" + var_name
        if multi == True:
            message = message + "_multi"
        if iveg != None:
            message = message + "_iveg="+str(iveg)
        plt.savefig('./plots/time_series_diff_'+message+'.png',dpi=300)

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
    if 0:
        case_names = [ "drght_2017_2019_bl_pbl2_mp4_sf_sfclay2",
                       "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"]

        time_s     = datetime(2017,2,1,0,0,0,0)
        time_e     = datetime(2017,7,31,23,59,0,0)

        for case_name in case_names:
            wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"
            LIS_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"

            file_paths = [
                          LIS_path+"LIS.CABLE.201701-201701.d01.nc",
                           LIS_path+"LIS.CABLE.201702-201702.d01.nc",
                           LIS_path+"LIS.CABLE.201703-201703.d01.nc",
                           LIS_path+"LIS.CABLE.201704-201704.d01.nc",
                           LIS_path+"LIS.CABLE.201705-201705.d01.nc",
                           LIS_path+"LIS.CABLE.201706-201706.d01.nc",
                           LIS_path+"LIS.CABLE.201707-201707.d01.nc",]

            message    = case_name+"_Jan_2017-Jul_2017"
            var_names  = [ "FWsoil_tavg", "SoilMoist_inst","GWwb_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                           "WaterTableD_tavg","Qle_tavg","Qh_tavg", "Qg_tavg","AvgSurfT_tavg","VegT_tavg",]
            lat_name   = "lat"
            lon_name   = "lon"
            for var_name in var_names:
                plot_time_series(file_paths, var_name, time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message)

    if 1:
        message        = "bl_pbl2_mp4_sf_sfclay2_CTL_ALB"
        lat_name       = "lat"
        lon_name       = "lon"
        iveg           = [2,5,6,9,14] #2
        
        case_name_ctl  = "drght_2017_2019_bl_pbl2_mp4_sf_sfclay2"
        case_name_sen  = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
        case_name_sen_2= "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR"

        time_s         = datetime(2017,1,1,0,0,0,0)
        time_e         = datetime(2020,7,1,0,0,0,0)

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_ctl+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"
        LIS_path_ctl   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_ctl+"/LIS_output/"
        LIS_path_sen   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_sen+"/LIS_output/"
        LIS_path_sen_2 = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_sen_2+"/LIS_output/"

        file_paths_ctl  = [
                            LIS_path_ctl+"LIS.CABLE.201701-202006_met.nc"
                            ]

        file_paths_sen = [
                            LIS_path_sen+"LIS.CABLE.201701-202006_met.nc"
                            ]

        file_paths_sen_2 = None
                        #    [LIS_path_sen_2+"LIS.CABLE.201701-201701.d01.nc",]

        var_names  = ['Tair_f_inst','Wind_f_inst']
                    #    [ "Evap_tavg","TVeg_tavg","ESoil_tavg"]
                    #    "WaterTableD_tavg","Qle_tavg","Qh_tavg", "Qg_tavg","AvgSurfT_tavg","VegT_tavg",'Albedo_inst',
                    #    'Tair_f_inst',"SoilMoist_inst", "FWsoil_tavg", "GWwb_tavg",,"ECanop_tavg","Qs_tavg","Qsb_tavg",]

        for var_name in var_names:
            plot_time_series_diff(file_paths_ctl,file_paths_sen,file_paths_sen_2, var_name,
                                  time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                                  lat_name=lat_name, lon_name=lon_name, message=message,
                                  multi=True, iveg=iveg)

    if 0:
        case_name       =  "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI"
                         # "drght_2017_2019_bl_pbl2_mp4_sf_sfclay2"
                         # "drght_2017_2019_bl_pbl5_mp6_sf_sfclay1"
        LIS_path        = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"

        time_s          = datetime(2017,1,1,0,0,0,0)
        time_e          = datetime(2017,4,30,23,59,0,0)
        # time_e          = datetime(2020,6,30,23,59,0,0)

        file_paths = [ LIS_path+"LIS.CABLE.201701-201701.d01.nc",
                       LIS_path+"LIS.CABLE.201702-201702.d01.nc",
                       LIS_path+"LIS.CABLE.201703-201703.d01.nc",
                       LIS_path+"LIS.CABLE.201704-201704.d01.nc",]

        message   = case_name
        var_names = [ "FWsoil_tavg", "SoilMoist_inst", "GWwb_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                      "WaterTableD_tavg","Qle_tavg","Qh_tavg", "Qg_tavg","AvgSurfT_tavg","VegT_tavg",]
        lat_name  = "lat"
        lon_name  = "lon"
        for var_name in var_names:
            plot_time_series(file_paths, var_name, time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message)

    if 0:
        # difference plots for drought-summer
        case_names = [ "drght_2017_2019_bl_pbl5_mp6_sf_sfclay1_1201",
                       "drght_2017_2019_bl_pbl5_mp6_sf_sfclay1_climatology"]

        time_s     = datetime(2018,12,1,0,0,0,0)
        time_e     = datetime(2019,2,28,23,59,0,0)

        for case_name in case_names:
            wrf_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/WRF_output/wrfout_d01_2018-12-01_01:00:00"
            LIS_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name+"/LIS_output/"

            file_paths = [ LIS_path+"LIS.CABLE.201812-201812.d01.nc",
                           LIS_path+"LIS.CABLE.201901-201901.d01.nc",
                           LIS_path+"LIS.CABLE.201902-201902.d01.nc",]

            message    = case_name+"_Dec_2018-Jan_2019"
            var_names  = [ "FWsoil_tavg", "SoilMoist_inst","GWwb_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                           "WaterTableD_tavg","Qle_tavg","Qh_tavg", "Qg_tavg","AvgSurfT_tavg","VegT_tavg",]
            lat_name   = "lat"
            lon_name   = "lon"
            for var_name in var_names:
                plot_time_series(file_paths, var_name, time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon, lat_name=lat_name, lon_name=lon_name, message=message)
