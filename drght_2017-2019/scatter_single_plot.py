#!/usr/bin/python

'''
Plot spitial map of land diagnosis and parameters from LIS-CABLE
1. per time step
2. time period average
cp /g/data/w97/mm3972/scripts/Groundwater_Atmosphere_Heatwave_Drought/src/Fig8_scatter_deltaT_WTD_PFT.py scatter_single_plot.py
'''

from netCDF4 import Dataset
import numpy as np
import pandas as pd
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import matplotlib.ticker as mticker
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

def mask_by_pft(land_path,case_name):

    print("In mask_by_pft Var")
    
    file_path  = land_path + case_name + '/LIS_output/Landcover_inst/LIS.CABLE.201701-201912.nc'
    f          = Dataset(file_path, mode='r')
    Time       = nc.num2date(f.variables['time'][:],f.variables['time'].units,
                 only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time       = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)
    Var        = f.variables['Landcover_inst']
    var        = time_clip_to_day(time,Var,time_s,time_e)
    # var        = np.where(var == iveg, var, np.nan) # awap_t == current for awap_t in AWAP_time
    
    print(var)

    return var

def read_data(land_path,case_name,var_name,pft,time_s,time_e):

    # ============ Read data ============
    file_path = land_path + case_name + '/LIS_output/' + var_name + '/LIS.CABLE.201701-201912.nc'

    f         = Dataset(file_path, mode='r')

    Time      = nc.num2date(f.variables['time'][:],f.variables['time'].units,
                only_use_cftime_datetimes=False, only_use_python_datetimes=True)
    time      = UTC_to_AEST(Time) - datetime(2000,1,1,0,0,0)
    Var       = f.variables[var_name]
    print("np.shape(Var)",np.shape(Var))
    var       = time_clip_to_day_max(time,Var,time_s,time_e)
    print("np.shape(var)",np.shape(var))

    t_s       = time_s - datetime(2000,1,1,0,0,0)
    t_e       = time_e - datetime(2000,1,1,0,0,0)
    
    time_3D   = np.zeros(np.shape(var))
 
    j=0
    for i in np.arange(t_s.days, t_e.days):
        # print("j=",j,"i=",i)
        time_3D[j,:,:] = i
        j += 1

    var_shrink = np.reshape(var,-1)
    time_shrink= np.reshape(time_3D,-1)
    var_1D     = var_shrink[~ np.isnan(var_shrink)]
    time_1D    = time_shrink[~ np.isnan(var_shrink)]

    df         = pd.DataFrame(var_1D, columns=['var'])
    df['time'] = time_1D

    if pft is not None:
        mask_pft   = mask_by_pft(land_path,case_name)
        pft_shrink = np.reshape(mask_pft,-1)
        pft_1D     = pft_shrink[~ np.isnan(var_shrink)]
        df['pft']  = pft_1D

    print(df)
  
    return df

def plot_spatial_land_days( land_path,case_names,var_name,pft,time_s,time_e, message=None):

    # ============= read data ================
    df_ctl = read_data(land_path,case_names[0],var_name,pft,time_s,time_e)
    df_sen = read_data(land_path,case_names[1],var_name,pft,time_s,time_e)

    # ============ Setting for plotting ============
    cmap     = plt.cm.BrBG #YlOrBr #coolwarm_r
    
    markers  = ["o","o","o","^","s"]
    mrk_sz   = 1.5

    # fig, ax  = plt.subplots(nrows=1, ncols=1, figsize=[12,8],sharex=True, sharey=True, squeeze=True) #
    fig, ax = plt.subplots(figsize=[20, 10])
    # plt.subplots_adjust(wspace=0.0, hspace=0.0)

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 14

    almost_black = '#262626'
    # change the tick colors also to the almost black
    plt.rcParams['ytick.color']     = almost_black
    plt.rcParams['xtick.color']     = almost_black

    # change the text colors also to the almost black
    plt.rcParams['text.color']      = almost_black

    # Change the default axis colors from black to a slightly lighter black,
    # and a little thinner (0.5 instead of 1)
    plt.rcParams['axes.edgecolor']  = almost_black
    plt.rcParams['axes.labelcolor'] = almost_black

    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    mask_1D = df_ctl['pft'] == pft

    sct    = ax.scatter(df_ctl[mask_1D]['time'], df_sen[mask_1D]['var']-df_ctl[mask_1D]['var'],  color='none', edgecolors='red',  s=9, 
                        marker=markers[0], alpha=0.05, cmap=cmap, label='ctl') #edgecolor='none', c='red'
    # sct    = ax.scatter(df_sen[mask_1D]['time']+0.5, df_sen[mask_1D]['var'], edgecolor='none', c='blue', s=7, 
                        # marker=markers[0], alpha=0.05, cmap=cmap, label='alb') 

    if pft is not None:
        fig.savefig("./plots/scatter_"+message+"_pft="+str(pft)+"_Tmax",bbox_inches='tight')
    else:
        fig.savefig("./plots/scatter_"+message+"_Tmax",bbox_inches='tight')

if __name__ == "__main__":

    # =============================== Operation ================================
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

    loc_lat    = [-33,-29]
    loc_lon    = [147,149]
    PFT        = False
    lat_name   = "lat"
    lon_name   = "lon"

    # =========================== Plot =============================
    message = "2017-2019_diff" #"2018_growth_season"

    time_s  = datetime(2017,1,2,0,0,0,0)
    time_e  = datetime(2020,1,1,0,0,0,0)

    if 1:
        land_path = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/'
        case_names= ['drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2',
                     'drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB']
        var_name  = 'Tair_f_inst'
        pft       = 2
        
        # pft      = ["BEF","crop","shrub","grass","barren"]
        # iveg_num = [2, 9, 5, 6, 14]
        plot_spatial_land_days(land_path,case_names,var_name,pft,time_s,time_e,loc_lat=loc_lat,loc_lon=loc_lon,
                                lat_name=lat_name, lon_name=lon_name, message=message)




