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

    # GPP scaling
    s2d               = 3600               # s-1 to d-1
    GPP_scale         = -0.000001*12*s2d   # umol s-1 to g d-1

    # read data
    time_ctl, Var_ctl = read_var_multi_file(file_paths_ctl, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time_sen, Var_sen = read_var_multi_file(file_paths_sen, var_name, loc_lat, loc_lon, lat_name, lon_name)
    Var_ctl           = np.nanmean(Var_ctl,axis=(1,2))*GPP_scale
    Var_sen           = np.nanmean(Var_sen,axis=(1,2))*GPP_scale

    print(message," Var_ctl",Var_ctl)
    print("Var_sen",Var_sen)

    # Var_daily_ctl     = time_clip_to_day(time_ctl, Var_ctl, time_s, time_e, seconds=None)
    # Var_daily_sen     = time_clip_to_day(time_sen, Var_sen, time_s, time_e, seconds=None)
    # Var_daily_ctl     = Var_daily_ctl*GPP_scale
    # Var_daily_sen     = Var_daily_sen*GPP_scale

    if file_paths_sen_2 != None:
        if message == "ALB_CTL_AU-Cum":
            file_sen_2    = Dataset(file_paths_sen_2[0],  mode='r')
            GPP_sen_2     = file_sen_2.variables['GPP'][4:6]
            GPP_sen_2_LL  = file_sen_2.variables['GPP_LL'][4:6]
            GPP_sen_2_SOLO= file_sen_2.variables['GPP_SOLO'][4:6]
            GPP_sen_2     = GPP_sen_2*0.000001*12*1800
            GPP_sen_2_LL  = GPP_sen_2_LL*0.000001*12*1800
            GPP_sen_2_SOLO= GPP_sen_2_SOLO*0.000001*12*1800
        elif message == "ALB_CTL_AU-Tum":
            file_sen_2    = Dataset(file_paths_sen_2[0],  mode='r')
            GPP_sen_2     = file_sen_2.variables['GPP'][15]
            GPP_sen_2_LL  = file_sen_2.variables['GPP_LL'][15]
            GPP_sen_2_SOLO= file_sen_2.variables['GPP_SOLO'][15]
            GPP_sen_2     = GPP_sen_2*0.000001*12*3600
            GPP_sen_2_LL  = GPP_sen_2_LL*0.000001*12*3600
            GPP_sen_2_SOLO= GPP_sen_2_SOLO*0.000001*12*3600
        print("GPP_sen_2",GPP_sen_2,"GPP_sen_2_LL",GPP_sen_2_LL,"GPP_sen_2_SOLO",GPP_sen_2_SOLO)

        # # interpolate observation to WRF domain
        # # read latitude and longitude from observation file
        # time_sen_2, lats_sen_2    = read_var(file_paths_sen_2[0], 'latitude', loc_lat, loc_lon, 'latitude', 'longitude')
        # time_sen_2, lons_sen_2    = read_var(file_paths_sen_2[0], 'longitude', loc_lat, loc_lon, 'latitude', 'longitude')
        # print("np.shape(lats_sen_2)",np.shape(lats_sen_2))
        # print("np.shape(lons_sen_2)",np.shape(lons_sen_2))
        # print("np.shape(Var_daily_sen_2)",np.shape(Var_daily_sen_2))

        # # read latitude and longitude from lis-cable file
        # lis_cable         = Dataset(file_paths_ctl[0],  mode='r')
        # lats_out          = lis_cable.variables['lat'][:,:]
        # lons_out          = lis_cable.variables['lon'][:,:]
        # sen_2_regrid      = np.zeros((len(Var_daily_sen_2[:,0,0]),len(lats_out[:,0]),len(lats_out[0,:])))
       
        # for i in np.arange(len(Var_daily_sen_2[:,0,0])):
        #     sen_2_regrid[i,:,:]    = regrid_data(lats_sen_2, lons_sen_2, lats_out, lons_out, Var_daily_sen_2[i,:,:], threshold=0)
        # print('np.shape(sen_2_regrid)',np.shape(sen_2_regrid))


    # colors        = [ "forestgreen", "yellowgreen","orange","red","black",]
    #                 # [ "black", "grey","lightcoral","red","orange","gold",
    #                 #  "yellow","yellowgreen","forestgreen","aquamarine","skyblue",
    #                 #  "blue","blueviolet","violet","purple","crimson","pink"]
    # cleaner_dates = ["2017","2018", "2019", "2020" ]
    # xtickslocs    = [0,        365,    730,  1095  ]

    # print(var_name, "dimension=",np.shape(Var_daily_ctl))
    # fig, ax = plt.subplots(figsize=[9,9])

    # if multi == None:
    #     if iveg == None:
    #         var = np.nanmean(Var_daily_diff[:,:,:],axis=(1,2)) # 1
    #         ax.plot(np.arange(len(var)), var, c = "blue", label="Δ"+var_name, alpha=0.5)
    #     elif np.shape(iveg)[0]>1:
    #         for i,pft in enumerate(iveg):
    #             var = np.nanmean(np.where(landcover == pft, Var_daily_diff, np.nan),axis=(1,2))
    #             ax.plot(np.arange(len(var)), var, c = colors[i], label="Δ"+var_name+" PFT="+str(pft), alpha=0.5)
    #     else:
    #         var = np.nanmean(np.where(landcover == iveg, Var_daily_diff, np.nan),axis=(1,2))
    #         ax.plot(np.arange(len(var)), var, c = "blue", label="Δ"+var_name, alpha=0.5)

    # if multi == True:
    #     if iveg == None:
    #         var_ctl = np.nanmean(Var_daily_ctl[:,:,:],axis=(1,2))
    #         var_sen = np.nanmean(Var_daily_sen[:,:,:],axis=(1,2))
    #         df_ctl  = pd.DataFrame({'ctl': var_ctl})
    #         df_sen  = pd.DataFrame({'sen': var_sen})
    #         ax.plot(df_ctl['ctl'].rolling(window=30).mean(), c = "red",  label=var_name+"_ctl", alpha=0.5) # np.arange(len(var_ctl)),
    #         ax.plot(df_sen['sen'].rolling(window=30).mean(), c = "blue", label=var_name+"_sen", alpha=0.5) # np.arange(len(var_sen)),
    #         if file_paths_sen_2 != None:
    #             var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:],axis=(1,2))
    #             df_sen_2  = pd.DataFrame({'sen_2': var_sen_2})
    #             ax.plot(df_sen_2['sen_2'].rolling(window=30).mean(), c = "green", label=var_name+"_sen_2", alpha=0.5)
    #     elif np.shape(iveg)[0]>1:
    #         for i,pft in enumerate(iveg):
    #             var_ctl = np.nanmean(np.where(landcover == pft, Var_daily_ctl, np.nan),axis=(1,2))
    #             var_sen = np.nanmean(np.where(landcover == pft, Var_daily_sen, np.nan),axis=(1,2))
    #             df_ctl  = pd.DataFrame({'ctl': var_ctl})
    #             df_sen  = pd.DataFrame({'sen': var_sen})
    #             ax.plot( df_ctl['ctl'].rolling(window=30).mean(), c = colors[i], label=var_name+"_ctl PFT="+str(pft), alpha=0.8) # .rolling(window=30).mean()
    #             ax.plot( df_sen['sen'].rolling(window=30).mean(), c = colors[i], label=var_name+"_sen PFT="+str(pft), alpha=0.5) # .rolling(window=30).mean()
    #             print("iveg = ",pft)
    #             if file_paths_sen_2 != None:
    #                 sen_2_tmp = sen_2_regrid*1.0
    #                 for j in np.arange(len(sen_2_regrid[:,0,0])):
    #                     sen_2_tmp[j,:,:] = np.where(landcover[0,:,:] == pft, sen_2_regrid[j,:,:], np.nan)
    #                 var_sen_2 = np.nanmean(sen_2_tmp,axis=(1,2))
    #                 df_sen_2  = pd.DataFrame({'sen_2': var_sen_2})
    #                 time_sen_2 = [6224,6255,6283,6314,6344,6375,6405,6436,6467,6497,6528,6558,6589,6620,6648,
    #                               6679,6709,6740,6770,6801,6832,6862,6893,6923,6954,6985,7013,7044,7074,7105,
    #                               7135,7166,7197,7227,7258,7288,7319,7350]
    #                 ax.plot(time_sen_2, df_sen_2['sen_2'], c = colors[i], label=var_name+"_sen_2 PFT="+str(pft), alpha=0.3) # .rolling(window=30).mean()
    #             for j in np.arange(len(df_ctl['ctl'])):
    #                 print("day = ", j, ", values = ", df_ctl['ctl'][j], " & ", df_sen['sen'][j])
    #     else:
    #         var_ctl = np.nanmean(np.where(landcover == iveg, Var_daily_ctl, np.nan),axis=(1,2))
    #         var_sen = np.nanmean(np.where(landcover == iveg, Var_daily_sen, np.nan),axis=(1,2))
    #         ax.plot(np.arange(len(var_ctl)), var_ctl, c = "red",  label=var_name+"_ctl", alpha=0.5)
    #         ax.plot(np.arange(len(var_sen)), var_sen, c = "blue", label=var_name+"_sen", alpha=0.5)
    #         if file_paths_sen_2 != None:
    #             var_sen_2 = np.nanmean(np.where(landcover == iveg, Var_daily_sen_2, np.nan),axis=(1,2))
    #             ax.plot(np.arange(len(var_sen_2)), var_sen_2, c = "green", label=var_name+"_sen_2", alpha=0.5)  



    # ax.set_title(var_name)
    # ax.set(xticks=xtickslocs, xticklabels=cleaner_dates)
    # ax.legend()
    # fig.tight_layout()
    # if message == None:
    #     message = var_name
    # else:
    #     message = message + "_" + var_name
    # if multi == True:
    #     message = message + "_multi"
    # if iveg != None and np.shape(iveg)[0] >1:
    #     message = message + "_iveg="+str(iveg)
    # elif iveg != None:
    #     message = message + "_iveg="+str(iveg[0])+"-"+str(iveg[-1])

    # plt.savefig('./plots/time_series_'+message+'.png',dpi=300)

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


        var_name       = "GPP_tavg"
        file_paths_ctl = [ LIS_path_ctl+var_name+'/LIS.CABLE.201701-201912_yearly.nc' ]
        file_paths_sen = [ LIS_path_sen+var_name+'/LIS.CABLE.201701-201912_yearly.nc' ]


        # AU-Cum
        message = "ALB_CTL_AU-Cum"
        loc_lat = [-33.6133-0.02,-33.6133+0.02]
        loc_lon = [150.7225-0.02,150.7225+0.02]
        time_s  = datetime(2017,1,1,0,0,0,0)
        time_e  = datetime(2019,1,1,0,0,0,0)
        file_paths_sen_2 = ['/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/check_GPP/AU-Cum_2013-2018_OzFlux_Flux_yearly.nc']
        plot_time_series(file_paths_ctl,file_paths_sen,file_paths_sen_2, var_name,
                            time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message,
                            multi=True, iveg=iveg)


        # AU-Tum
        message = "ALB_CTL_AU-Tum"
        loc_lat = [-35.6566-0.02,-35.6566+0.02]
        loc_lon = [148.1517-0.02,148.1517+0.02]
        time_s  = datetime(2017,1,1,0,0,0,0)
        time_e  = datetime(2018,1,1,0,0,0,0)
        file_paths_sen_2 = ['/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/check_GPP/AU-Tum_2002-2017_OzFlux_Flux_yearly.nc']
        plot_time_series(file_paths_ctl,file_paths_sen,file_paths_sen_2, var_name,
                            time_s=time_s,time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message,
                            multi=True, iveg=iveg)
