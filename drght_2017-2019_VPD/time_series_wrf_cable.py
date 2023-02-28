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
1. 8 Feb 2023: cp ./drght_2017-2019_VPD/time_series_wrf_cable.py ./drght_2017-2019_VPD/
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

    cleaner_dates = ["2017","2018", "2019", "2020" ]
    xtickslocs    = [0,     365,      730,   1095  ]

    if len(np.shape(Var_daily)) == 3:

        print(var_name, "dimension=",np.shape(Var_daily))

        var = np.nanmean(Var_daily[1:,:,:],axis=(1,2))
        fig, ax = plt.subplots()
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

        fig, ax = plt.subplots()
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
        time_lc, LC    = read_var_multi_file(file_paths_ctl, "Landcover_inst", loc_lat, loc_lon, lat_name, lon_name)
        landcover      = time_clip_to_day(time_lc, LC, time_s, time_e, seconds=None)

    cleaner_dates = ["2017","2018", "2019", "2020" ]
    xtickslocs    = [0,     365,      730,   1095  ]

    if len(np.shape(Var_daily_ctl)) == 3:

        print(var_name, "dimension=",np.shape(Var_daily_ctl))
        fig, ax = plt.subplots()

        if multi == None:
            if iveg == None:
                var = np.nanmean(Var_daily_diff[:,:,:],axis=(1,2)) # 1
            else:
                var = np.nanmean(np.where(landcover == iveg, Var_daily_diff, np.nan),axis=(1,2))
            ax.plot(np.arange(len(var)), var, c = "blue", label=var_name+"_diff", alpha=0.5)

        if multi == True:
            if iveg == None:
                var_ctl = np.nanmean(Var_daily_ctl[:,:,:],axis=(1,2))
                var_sen = np.nanmean(Var_daily_sen[:,:,:],axis=(1,2))
            else:
                var_ctl = np.nanmean(np.where(landcover == iveg, Var_daily_ctl, np.nan),axis=(1,2))
                var_sen = np.nanmean(np.where(landcover == iveg, Var_daily_sen, np.nan),axis=(1,2))

            ax.plot(np.arange(len(var_ctl)), var_ctl, c = "red",  label=var_name+"_ctl", alpha=0.5)
            ax.plot(np.arange(len(var_sen)), var_sen, c = "blue", label=var_name+"_sen", alpha=0.5)

            if file_paths_sen_2 != None:
                if iveg == None:
                    var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:],axis=(1,2))
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
        if iveg != None:
            message = message + "_iveg="+str(iveg)

        plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

    elif len(np.shape(Var_daily_ctl)) == 4:

        print(var_name, "dimension=",np.shape(Var_daily_ctl))

        fig, ax = plt.subplots()
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
        message        = "drght_2017_2019_bl_pbl2_mp4_sf_sfclay2_VS"
        lat_name       = "lat"
        lon_name       = "lon"
        iveg           = 2

        case_name_ctl  = "drght_2017_2019_bl_pbl2_mp4_sf_sfclay2"
                            #"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
        case_name_sen  = "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
                             #"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR"

        case_name_sen_2= "drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR"

        time_s         = datetime(2017,1,1,0,0,0,0)
        time_e         = datetime(2019,1,1,0,0,0,0)

        wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_ctl+"/WRF_output/wrfout_d01_2017-02-01_06:00:00"
        LIS_path_ctl   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_ctl+"/LIS_output/"
        LIS_path_sen   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_sen+"/LIS_output/"
        LIS_path_sen_2 = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/"+case_name_sen_2+"/LIS_output/"

        file_paths_ctl  = [
                            LIS_path_ctl+"LIS.CABLE.201701-201701.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201702-201702.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201703-201703.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201704-201704.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201705-201705.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201706-201706.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201707-201707.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201708-201708.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201709-201709.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201710-201710.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201711-201711.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201712-201712.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201801-201801.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201802-201802.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201803-201803.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201804-201804.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201805-201805.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201806-201806.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201807-201807.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201808-201808.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201809-201809.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201810-201810.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201811-201811.d01.nc",
                            LIS_path_ctl+"LIS.CABLE.201812-201812.d01.nc",
                            ]

        file_paths_sen = [
                            LIS_path_sen+"LIS.CABLE.201701-201701.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201702-201702.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201703-201703.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201704-201704.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201705-201705.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201706-201706.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201707-201707.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201708-201708.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201709-201709.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201710-201710.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201711-201711.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201712-201712.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201801-201801.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201802-201802.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201803-201803.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201804-201804.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201805-201805.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201806-201806.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201807-201807.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201808-201808.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201809-201809.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201810-201810.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201811-201811.d01.nc",
                            LIS_path_sen+"LIS.CABLE.201812-201812.d01.nc"
                            ]

        file_paths_sen_2 = [
                            LIS_path_sen_2+"LIS.CABLE.201701-201701.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201702-201702.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201703-201703.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201704-201704.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201705-201705.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201706-201706.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201707-201707.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201708-201708.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201709-201709.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201710-201710.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201711-201711.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201712-201712.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201801-201801.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201802-201802.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201803-201803.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201804-201804.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201805-201805.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201806-201806.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201807-201807.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201808-201808.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201809-201809.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201810-201810.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201811-201811.d01.nc",
                            LIS_path_sen_2+"LIS.CABLE.201812-201812.d01.nc"
                            ]

        var_names  = [ "GWwb_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                       "WaterTableD_tavg","Qle_tavg","Qh_tavg", "Qg_tavg","AvgSurfT_tavg","VegT_tavg",'Albedo_inst',
                       'Tair_f_inst',"SoilMoist_inst", "FWsoil_tavg", ]

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
                       # LIS_path+"LIS.CABLE.201705-201705.d01.nc",
                       # LIS_path+"LIS.CABLE.201706-201706.d01.nc",
                       # LIS_path+"LIS.CABLE.201707-201707.d01.nc",]
                       # LIS_path+"LIS.CABLE.201708-201708.d01.nc",
                       # LIS_path+"LIS.CABLE.201709-201709.d01.nc",
                       # LIS_path+"LIS.CABLE.201710-201710.d01.nc",
                       # LIS_path+"LIS.CABLE.201711-201711.d01.nc",
                       # LIS_path+"LIS.CABLE.201712-201712.d01.nc",
                       # LIS_path+"LIS.CABLE.201801-201801.d01.nc",
                       # LIS_path+"LIS.CABLE.201802-201802.d01.nc",
                       # LIS_path+"LIS.CABLE.201803-201803.d01.nc",
                       # LIS_path+"LIS.CABLE.201804-201804.d01.nc",
                       # LIS_path+"LIS.CABLE.201805-201805.d01.nc",
                       # LIS_path+"LIS.CABLE.201806-201806.d01.nc",
                       # LIS_path+"LIS.CABLE.201807-201807.d01.nc",
                       # LIS_path+"LIS.CABLE.201808-201808.d01.nc",
                       # LIS_path+"LIS.CABLE.201809-201809.d01.nc",
                       # LIS_path+"LIS.CABLE.201810-201810.d01.nc",
                       # LIS_path+"LIS.CABLE.201811-201811.d01.nc",
                       # LIS_path+"LIS.CABLE.201812-201812.d01.nc",
                       # LIS_path+"LIS.CABLE.201901-201901.d01.nc",
                       # LIS_path+"LIS.CABLE.201902-201902.d01.nc",
                       # LIS_path+"LIS.CABLE.201903-201903.d01.nc",
                       # LIS_path+"LIS.CABLE.201904-201904.d01.nc",
                       # LIS_path+"LIS.CABLE.201905-201905.d01.nc",
                       # LIS_path+"LIS.CABLE.201906-201906.d01.nc",
                       # LIS_path+"LIS.CABLE.201907-201907.d01.nc",
                       # LIS_path+"LIS.CABLE.201908-201908.d01.nc",
                       # LIS_path+"LIS.CABLE.201909-201909.d01.nc",
                       # LIS_path+"LIS.CABLE.201910-201910.d01.nc",
                       # LIS_path+"LIS.CABLE.201911-201911.d01.nc",
                       # LIS_path+"LIS.CABLE.201912-201912.d01.nc",
                       # LIS_path+"LIS.CABLE.202001-202001.d01.nc",
                       # LIS_path+"LIS.CABLE.202002-202002.d01.nc",
                       # LIS_path+"LIS.CABLE.202003-202003.d01.nc",
                       # LIS_path+"LIS.CABLE.202004-202004.d01.nc",]
                    #    LIS_path+"LIS.CABLE.202005-202005.d01.nc",
                    #    LIS_path+"LIS.CABLE.202006-202006.d01.nc",]

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
