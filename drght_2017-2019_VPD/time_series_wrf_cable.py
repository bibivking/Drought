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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import griddata
import matplotlib.ticker as mticker
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

def plot_time_series(file_paths, var_name, time_s=None, time_e=None, loc_lat=None, loc_lon=None, lat_name=None, lon_name=None, message=None):

    print("======== In plot_time_series =========")
    time, Var = read_var_multi_file(file_paths, var_name, loc_lat, loc_lon, lat_name, lon_name)

    if var_name in ["Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                    "Evap","TVeg","ESoil","ECanop","Qs","Qsb"]:
        Var_daily = time_clip_to_day_sum(time, Var, time_s, time_e, seconds=None)
        Var_daily = Var_daily*3600.
    elif var_name in ["WaterTableD_tavg","WatTable"]:
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
    if var_name == "EF":
        time_ctl, Var_ctl_Qle = read_var_multi_file(file_paths_ctl, "Qle_tavg", loc_lat, loc_lon, lat_name, lon_name)
        time_ctl, Var_ctl_Qh  = read_var_multi_file(file_paths_ctl, "Qh_tavg", loc_lat, loc_lon, lat_name, lon_name)

        time_sen, Var_sen_Qle = read_var_multi_file(file_paths_sen, "Qle_tavg", loc_lat, loc_lon, lat_name, lon_name)
        time_sen, Var_sen_Qh  = read_var_multi_file(file_paths_sen, "Qh_tavg", loc_lat, loc_lon, lat_name, lon_name)

        ctl_QleQh = Var_ctl_Qle+Var_ctl_Qh
        sen_QleQh = Var_sen_Qle+Var_sen_Qh
        Var_ctl = np.where(abs(ctl_QleQh)>0.01, Var_ctl_Qle/ctl_QleQh,np.nan)
        Var_sen = np.where(abs(sen_QleQh)>0.01, Var_sen_Qle/sen_QleQh,np.nan)
    elif var_name == "VPD":
        time_ctl, Tair_ctl = read_var_multi_file(file_paths_ctl, "Tair", loc_lat, loc_lon, lat_name, lon_name)
        time_ctl, Qair_ctl = read_var_multi_file(file_paths_ctl, "Qair", loc_lat, loc_lon, lat_name, lon_name)
        time_ctl, Pres_ctl = read_var_multi_file(file_paths_ctl, "PSurf", loc_lat, loc_lon, lat_name, lon_name)
        Var_ctl            = qair_to_vpd(Qair_ctl, Tair_ctl, Pres_ctl)

        time_sen, Tair_sen = read_var_multi_file(file_paths_sen, "Tair", loc_lat, loc_lon, lat_name, lon_name)
        time_sen, Qair_sen = read_var_multi_file(file_paths_sen, "Qair", loc_lat, loc_lon, lat_name, lon_name)
        time_sen, Pres_sen = read_var_multi_file(file_paths_sen, "PSurf", loc_lat, loc_lon, lat_name, lon_name)
        Var_sen            = qair_to_vpd(Qair_sen, Tair_sen, Pres_sen)
    elif var_name == "SM35":
        time_ctl, Var_ctl_SM = read_var_multi_file(file_paths_ctl, "SoilMoist", loc_lat, loc_lon, lat_name, lon_name)
        time_sen, Var_sen_SM = read_var_multi_file(file_paths_sen, "SoilMoist", loc_lat, loc_lon, lat_name, lon_name)
        Var_ctl              = (Var_ctl_SM[:,0,:,:]*0.022+Var_ctl_SM[:,1,:,:]*0.058+Var_ctl_SM[:,2,:,:]*0.154+Var_ctl_SM[:,3,:,:]*0.116)/0.35
        Var_sen              = (Var_sen_SM[:,0,:,:]*0.022+Var_sen_SM[:,1,:,:]*0.058+Var_sen_SM[:,2,:,:]*0.154+Var_sen_SM[:,3,:,:]*0.116)/0.35
    else:
        time_ctl, Var_ctl = read_var_multi_file(file_paths_ctl, var_name, loc_lat, loc_lon, lat_name, lon_name)
        time_sen, Var_sen = read_var_multi_file(file_paths_sen, var_name, loc_lat, loc_lon, lat_name, lon_name)

    if file_paths_sen_2 != None:

        if var_name == "EF":
            time_sen_2, Var_sen_2_Qle = read_var_multi_file(file_paths_sen_2, "Qle_tavg", loc_lat, loc_lon, lat_name, lon_name)
            time_sen_2, Var_sen_2_Qh  = read_var_multi_file(file_paths_sen_2, "Qh_tavg", loc_lat, loc_lon, lat_name, lon_name)
            sen_2_QleQh = Var_sen_2_Qle+Var_sen_2_Qh
            Var_sen_2 = np.where(abs(sen_2_QleQh)>0.01, Var_sen_2_Qle/sen_2_QleQh,np.nan)
        elif var_name == "VPD":
            time_sen_2, Tair_sen_2 = read_var_multi_file(file_paths_sen_2, "Tair", loc_lat, loc_lon, lat_name, lon_name)
            time_sen_2, Qair_sen_2 = read_var_multi_file(file_paths_sen_2, "Qair", loc_lat, loc_lon, lat_name, lon_name)
            time_sen_2, Pres_sen_2 = read_var_multi_file(file_paths_sen_2, "PSurf", loc_lat, loc_lon, lat_name, lon_name)
            Var_sen_2              = qair_to_vpd(Qair_sen_2, Tair_sen_2, Pres_sen_2)
        elif var_name == "SM35":
            time_sen_2, Var_sen_2_SM = read_var_multi_file(file_paths_sen_2, "SoilMoist", loc_lat, loc_lon, lat_name, lon_name)
            Var_sen_2                = (Var_sen_2_SM[:,0,:,:]*0.022+Var_sen_2_SM[:,1,:,:]*0.058+Var_sen_2_SM[:,2,:,:]*0.154+Var_sen_2_SM[:,3,:,:]*0.116)/0.35
        else:
            time_sen_2, Var_sen_2 = read_var_multi_file(file_paths_sen_2, var_name, loc_lat, loc_lon, lat_name, lon_name)

    if var_name in ["Rainf_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                    "Rainf","Evap","TVeg","ESoil","ECanop","Qs","Qsb"]:
        Var_daily_ctl = time_clip_to_day_sum(time_ctl, Var_ctl, time_s, time_e, seconds=None)
        Var_daily_ctl = Var_daily_ctl*3600.

        Var_daily_sen = time_clip_to_day_sum(time_sen, Var_sen, time_s, time_e, seconds=None)
        Var_daily_sen = Var_daily_sen*3600.

    elif var_name in ["WaterTableD_tavg","WatTable"]:
        Var_daily_ctl = time_clip_to_day(time_ctl, Var_ctl, time_s, time_e, seconds=None)
        Var_daily_ctl = Var_daily_ctl/1000.

        Var_daily_sen = time_clip_to_day(time_sen, Var_sen, time_s, time_e, seconds=None)
        Var_daily_sen = Var_daily_sen/1000.

    else:
        Var_daily_ctl = time_clip_to_day(time_ctl, Var_ctl, time_s, time_e, seconds=None)
        Var_daily_sen = time_clip_to_day(time_sen, Var_sen, time_s, time_e, seconds=None)

    if file_paths_sen_2 != None:
        if var_name in ["Rainf_tavg","Evap_tavg","TVeg_tavg","ESoil_tavg","ECanop_tavg","Qs_tavg","Qsb_tavg",
                        "Rainf","Evap","TVeg","ESoil","ECanop","Qs","Qsb"]:
            Var_daily_sen_2 = time_clip_to_day_sum(time_sen_2, Var_sen_2, time_s, time_e, seconds=None)
            Var_daily_sen_2 = Var_daily_sen_2*3600.
        elif var_name in ["WaterTableD_tavg","WatTable"]:
            Var_daily_sen_2 = time_clip_to_day(time_sen_2, Var_sen_2, time_s, time_e, seconds=None)
            Var_daily_sen_2 = Var_daily_sen_2/1000.
        else:
            Var_daily_sen_2 = time_clip_to_day(time_sen_2, Var_sen_2, time_s, time_e, seconds=None)

    if multi == None:
        Var_daily_diff = Var_daily_sen - Var_daily_ctl

    if iveg != None:
        LC_file        = nc.Dataset(file_paths_ctl[0], mode='r')
        LC             = LC_file.variables["iveg"][:,:]
        mask           = mask_by_lat_lon(file_paths_ctl[0], loc_lat, loc_lon, lat_name, lon_name)
        LC_temp        = np.where(mask,LC,np.nan)
        ntime          = np.shape(Var_daily_ctl)[0]
        nlat           = np.shape(Var_daily_ctl)[-2]
        nlon           = np.shape(Var_daily_ctl)[-1]
        print("ntime = ", ntime)
        landcover      = np.zeros((ntime,nlat,nlon))

        for i in np.arange(ntime):
            landcover[i,:,:]= LC_temp[:,:]

    # _____________________ Plotting _____________________
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
                    ax.plot( df_ctl['ctl'].rolling(window=90).mean(), c = colors[i],  label=var_name+"_ctl PFT="+str(pft), alpha=0.8)
                    ax.plot( df_sen['sen'].rolling(window=90).mean(), c = colors[i], label=var_name+"_sen PFT="+str(pft), alpha=0.5)
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

        fig, ax     = plt.subplots(figsize=[9,9])
        labels_ctl  = ["lyr1_ctl","lyr2_ctl","lyr3_ctl","lyr4_ctl","lyr5_ctl","lyr6_ctl"]
        labels_sen  = ["lyr1_sen","lyr2_sen","lyr3_sen","lyr4_sen","lyr5_sen","lyr6_sen"]
        labels_sen2 = ["lyr1_sen2","lyr2_sen2","lyr3_sen2","lyr4_sen2","lyr5_sen2","lyr6_sen2"]

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
                ax.plot(np.arange(len(var_ctl[:,0])), var_ctl, ls = "-",  label=labels_ctl, alpha=0.5) #c = "blue",
                ax.plot(np.arange(len(var_sen[:,0])), var_sen, ls = "-.", label=labels_sen, alpha=0.5) #c = "blue",
            elif np.shape(iveg)[0]>1:
                print("Warning: Requiring too many lines in one plot, please only select one vegetation type")
                return
            else:
                var_ctl = np.nanmean(Var_daily_ctl[:,:,:,:],axis=(2,3))
                var_sen = np.nanmean(Var_daily_sen[:,:,:,:],axis=(2,3))
                for i in np.arange(6):
                    var_ctl[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_ctl[:,i,:,:], np.nan),axis=(1,2))
                    var_sen[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_sen[:,i,:,:], np.nan),axis=(1,2))
                print("np.shape(var_ctl) = ",np.shape(var_ctl))
                print("var_ctl[:,0]=",var_ctl[:,0])
                ax.plot(np.arange(len(var_ctl[:,0])), var_ctl, ls = "-",  label=labels_ctl, alpha=0.5) #c = "blue",
                ax.plot(np.arange(len(var_sen[:,0])), var_sen, ls = "-.", label=labels_sen, alpha=0.5) #c = "blue",

            if file_paths_sen_2 != None:
                if iveg == None:
                    var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:,:],axis=(2,3))
                    ax.plot(np.arange(len(var_sen_2[:,0])), var_sen_2, ls = "--", label=var_name+"_sen_2", alpha=0.5)
                else:
                    var_sen_2 = np.nanmean(Var_daily_sen_2[:,:,:,:],axis=(2,3))
                    for i in np.arange(6):
                        var_sen_2[:,i] = np.nanmean(np.where(landcover == iveg, Var_daily_sen_2[:,i,:,:], np.nan),axis=(1,2))
                    ax.plot(np.arange(len(var_sen_2[:,0])), var_sen_2, ls = "-.", label=labels_sen2, alpha=0.5) #c = "blue",

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

    if 1:
        message               = "T_Q_LWdown_detrend_vs_ctl_litter_on_PM"
        var_names             = ["SM35","Evap","Qle","Qh","GPP","Fwsoil","SoilMoist",] #"SoilMoist","Fwsoil","VPD",
        iveg                  = [2,5,6,9,14] #2#[2]#
        time_s                = datetime(2000,1,1,0,0,0,0)
        time_e                = datetime(2020,1,1,0,0,0,0)
        lat_name              = "latitude"
        lon_name              = "longitude"

        file_paths_ctl        = ["/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2000.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2001.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2002.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2003.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2004.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2005.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2006.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2007.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2008.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2009.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2010.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2011.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2012.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2013.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2014.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2015.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2016.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2017.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2018.nc",
                                 "/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/litter_on_PM/outputs/cable_out_2019.nc",]
                                 #"/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/outputs/cable_out_2017.nc",
                                 #"/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/outputs/cable_out_2018.nc",
                                 #"/g/data/w97/mm3972/model/cable/runs/VPD_drought/ctl/outputs/cable_out_2019.nc"]

        file_paths_sen        = [ "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2000.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2001.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2002.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2003.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2004.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2005.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2006.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2007.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2008.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2009.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2010.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2011.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2012.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2013.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2014.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2015.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2016.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2017.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2018.nc",
                                  "/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_LWdown_detrend_2000_2019/litter_on_PM/outputs/cable_out_2019.nc",]
                                 #"/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_detrend_2000_2019/outputs/cable_out_2017.nc",
                                 #"/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_detrend_2000_2019/outputs/cable_out_2018.nc",
                                 #"/g/data/w97/mm3972/model/cable/runs/VPD_drought/T_Q_detrend_2000_2019/outputs/cable_out_2019.nc"]

        for var_name in var_names:
            plot_time_series_diff(file_paths_ctl, file_paths_sen, file_paths_sen_2=None, var_name=var_name, time_s=time_s, time_e=time_e, loc_lat=loc_lat, loc_lon=loc_lon,
                            lat_name=lat_name, lon_name=lon_name, message=message,multi=True, iveg=iveg)
