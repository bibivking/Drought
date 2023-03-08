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
from scipy.signal import detrend
import matplotlib.ticker as mticker
from common_utils import *

def read_AWAP_file(file_path, var_name):

    '''
    Read AWAP data, output time coordinate and variable array
    '''

    print(var_name)

    print("file_path = ", file_path)

    # Initilizing
    time = []

    # Read in data
    var_file = Dataset(file_path, mode='r')
    time_tmp = nc.num2date(var_file.variables['time'][:],var_file.variables['time'].units,
                only_use_cftime_datetimes=False, only_use_python_datetimes=True)

    time = time_tmp
    print("time",time)

    Var_tmp = var_file.variables[var_name][:]
    if hasattr(var_file.variables[var_name], '_FillValue'):
        def_val = var_file.variables[var_name]._FillValue
        Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
    elif hasattr(var_file.variables[var_name], '_fillvalue'):
        def_val = var_file.variables[var_name]._fillvalue
        Var = np.where(Var_tmp == def_val, np.nan, Var_tmp)
    else:
        Var = Var_tmp

    var = Var

    print("=== In read_var_multi_file ===")
    print("time = ", time)
    print("var = ", var)

    return time,var

def detrend_add_intercept(time,var,year_s=1970,scale="monthly"):

    time_sum = np.shape(var)[0]
    lat_sum  = np.shape(var)[1]
    lon_sum  = np.shape(var)[2]
    Var_detr = np.zeros((time_sum,lat_sum,lon_sum))

    print("var dimensions are ", time_sum," x ", lat_sum," x ",lon_sum)

    month  = []
    season = []
    year   = []

    # set
    for t in time:
        print('t = ',t)
        mon = t.month
        yr  = t.year
        print('mon = ',mon)
        print('yr = ',yr)
        print('season = ',(mon+9)%12//3+1)
        month.append(mon+ (yr-year_s)*12)
        season.append((mon+9)%12//3+1 + (yr-year_s)*4) # 1:MAM 2:JJA 3:SON 4:DJF
        year.append(yr)

    year   = np.array(year)
    month  = np.array(month)
    season = np.array(season)

    print("year,month,season",year,month,season)

    if scale == "monthly":
        month_unq    = np.unique(month)
        month_sum    = len(month_unq)
        var_resample = np.zeros((month_sum,lat_sum,lon_sum))
        t = 0
        for m in month_unq:
            var_resample[t,:,:] = np.nanmean(var[month == m,:,:],axis=0)
            t = t + 1

    elif scale == "seasonal":
        season_unq    = np.unique(season)
        season_sum    = len(season_unq)
        var_resample  = np.zeros(season_sum,lat_sum,lon_sum)
        t = 0
        for s in season_unq:
            var_resample[t,:,:] = np.nanmean(var[season == s,:,:],axis=0)
            t = t + 1

    elif scale == "annual":
        year_unq     = np.unique(year)
        year_sum     = len(year_unq)
        var_resample = np.zeros(year_sum,lat_sum,lon_sum)
        t = 0
        for y in year_unq:
            var_resample[t,:,:] = np.nanmean(var[year == y,:,:],axis=0)
            t = t + 1
    else:
        print("Please choose scale!")

    var_detr  = detrend(var_resample, axis=0, type='linear', bp=0, overwrite_data=False)

    var_diff  = var_resample - var_detr

    Var_itc   = var_resample[0,:,:] - var_detr[0,:,:] # intercept
    var_diff  = var_diff - Var_itc                    # remove the intercept from the difference

    if scale == "monthly":
        for t, m in enumerate(month_unq):
            Var_detr[month == m,:,:]  = var[month == m,:,:] - var_diff[t,:,:]
    elif scale == "seasonal":
        for t,s in enumerate(season_unq):
            Var_detr[season == s,:,:] = var[season == s,:,:] - var_diff[t,:,:]
    elif scale == "annual":
        for t,y in enumerate(year_unq):
            Var_detr[year == y,:,:]   = var[year == y,:,:] - var_diff[t,:,:]
    else:
        print("Please choose scale!")

    return Var_detr

def output_detrend_AWAP_daily(file_path, var_name=None, scale="daily", message=None):

    # read dataset
    time, Var = read_AWAP_file(file_path, var_name)

    # remove linear trend
    if scale == "daily":
        Var_detr  = detrend(Var, axis=0, type='linear', bp=0, overwrite_data=False)
        Var_itc   = Var[0,:,:] - Var_detr[0,:,:] # intercept
        Var_detr  = Var_detr + Var_itc # add the intercept back
    else:
        Var_detr  = detrend_add_intercept(time,Var,1970,scale)


    # save the detrend data
    var_out   = np.where(Var_detr == np.nan, 999., Var_detr)
    file      = Dataset(file_path, mode='r+')
    file.variables[var_name][:]=var_out[:]
    file      = None

def output_binary_AWAP_daily(file_path,var_name):

    time, Var = read_AWAP_file(file_path, var_name)

    for i,t in enumerate(time):
        print('t = ',t)
        yr  = t.year
        mon = t.month
        dy  = t.day
        print('yr = ',yr)
        print('mon = ',mon)
        print('dy = ',dy)

        out_file = "/g/data/w97/mm3972/data/AWAP_detrend/"+var_name+"/bom-"+var_name+"-day-"+"%4d%02d%02d-%4d%02d%02d" %(yr,mon,dy,yr,mon,dy)+".flt"
        print("out_file =",out_file)
        flt_file = open(out_file, "wb") # "b":  Binary mode
        flt_file.write(Var[i,:,:])

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

    message               = "AWAP_detrend"
    lat_name              = "latitude"
    lon_name              = "longitude"
    scale                 = "daily"

    #var_name              = "vph09" # rad  rain  tmax  tmean  tmin  vph09  vph15  windspeed
    #file_path             = "/g/data/w97/mm3972/data/AWAP_detrend/"+var_name+"/AWAP_daily_"+var_name+"_1970_2019.nc"
    #output_detrend_AWAP_daily(file_path, var_name, scale=scale, message=message)

    #var_name              = "vph15" # rad  rain  tmax  tmean  tmin  vph09  vph15  windspeed
    #file_path             = "/g/data/w97/mm3972/data/AWAP_detrend/"+var_name+"/AWAP_daily_"+var_name+"_1970_2019.nc"
    #output_detrend_AWAP_daily(file_path, var_name, scale=scale, message=message)

    var_name              = "tmin"
    file_path             = "/g/data/w97/mm3972/data/AWAP_detrend/"+var_name+"/AWAP_daily_"+var_name+"_1970_2019.nc"
    output_binary_AWAP_daily(file_path,var_name)

    var_name              = "vph09"
    file_path             = "/g/data/w97/mm3972/data/AWAP_detrend/"+var_name+"/AWAP_daily_"+var_name+"_1970_2019.nc"
    output_binary_AWAP_daily(file_path,var_name)

    var_name              = "vph15"
    file_path             = "/g/data/w97/mm3972/data/AWAP_detrend/"+var_name+"/AWAP_daily_"+var_name+"_1970_2019.nc"
    output_binary_AWAP_daily(file_path,var_name)
