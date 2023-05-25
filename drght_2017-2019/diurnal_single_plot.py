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

def read_data(land_path,case_name,var_name,pft,time_s,time_e,loc_lat=None,loc_lon=None,
              lat_name=None,lon_name=None):

    # ============ Read data ============
    file_path = land_path + case_name + '/LIS_output/' + var_name + '/LIS.CABLE.201701-201912.nc'
    time, Var = read_var_multi_file(file_path, var_name, loc_lat, loc_lon, lat_name, lon_name)
    time, Lat = read_var(file_path, lat_name, loc_lat, loc_lon, lat_name, lon_name)
    time, Lon = read_var(file_path, lon_name, loc_lat, loc_lon, lat_name, lon_name)

    plt.contourf(Lat)
    plt.show()
    plt.contourf(Lon)
    plt.show()

    hour_3D   = np.zeros(np.shape(var))
    day_3D    = np.zeros(np.shape(var))
    lat_3D    = np.zeros(np.shape(var))
    lon_3D    = np.zeros(np.shape(var))
    
    for i,t in enumerate(time):
        hour_3D[i,:,:] = t.hours
        day_3D[i,:,:]  = t.days
        lat_3D[i,:,:]  = Lat
        lon_3D[i,:,:]  = Lon


    var_shrink = np.reshape(var,-1)
    hour_shrink= np.reshape(hour_3D,-1)
    day_shrink = np.reshape(day_3D,-1)
    lat_shrink = np.reshape(lat_3D,-1)
    lon_shrink = np.reshape(lon_3D,-1)

    var_1D     = var_shrink[~ np.isnan(var_shrink)]
    hour_1D    = hour_shrink[~ np.isnan(var_shrink)]
    day_1D     = day_shrink[~ np.isnan(var_shrink)]
    lat_1D     = lat_shrink[~ np.isnan(var_shrink)]
    lon_1D     = lon_shrink[~ np.isnan(var_shrink)]

    df         = pd.DataFrame(var_1D, columns=['var'])
    df['hour'] = hour_1D
    df['day']  = day_1D
    df['lat']  = lat_1D
    df['lon']  = lon_1D

    if pft is not None:
        mask_pft   = mask_by_pft(land_path,case_name)
        pft_shrink = np.reshape(mask_pft,-1)
        pft_1D     = pft_shrink[~ np.isnan(var_shrink)]
        df['pft']  = pft_1D

    print(df)
  
    return df

def plot_spatial_land_days(land_path,case_names,var_name,pft,time_s,time_e,loc_lat=None,loc_lon=None,
                           lat_name=None, lon_name=None, message=None):

    # ============= read data ================
    df_ctl = read_data(land_path,case_names[0],var_name,pft,time_s,time_e,loc_lat=loc_lat,loc_lon=loc_lon,
                        lat_name=lat_name, lon_name=lon_name)
    df_sen = read_data(land_path,case_names[1],var_name,pft,time_s,time_e,loc_lat=loc_lat,loc_lon=loc_lon,
                        lat_name=lat_name, lon_name=lon_name)

    t_s       = time_s - datetime(2000,1,1,0,0,0)
    t_e       = time_e - datetime(2000,1,1,0,0,0)
    

    # Prepare dataset 
    df_ctl    = df_ctl[df_ctl["day"].values >= t_s and df_ctl["day"].values < t_e]
    df_sen    = 


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
    if pft is not None:
        mask_1D = df_ctl['pft'] == pft
        sct     = ax.scatter(df_ctl[mask_1D]['time'], df_sen[mask_1D]['var']-df_ctl[mask_1D]['var'],  color='none', edgecolors='red',  s=9, 
                            marker=markers[0], alpha=0.05, cmap=cmap, label='ctl') #edgecolor='none', c='red'
    else: 
        sct     = ax.scatter(df_ctl['time'], df_sen['var']-df_ctl['var'], color='none', edgecolors='red', s=9, marker=markers[0], alpha=0.05, cmap=cmap, label='ctl')

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

    # # small region
    # loc_lat    = [-33,-29]
    # loc_lon    = [147,149]

    # east coast
    loc_lat    = [-33,-27]
    loc_lon    = [152,154]

    PFT        = False
    lat_name   = "lat"
    lon_name   = "lon"

    # =========================== Plot =============================
    message = "2017-2019_diff_east_coast" #"2018_growth_season"

    time_s  = datetime(2017,1,2,0,0,0,0)
    time_e  = datetime(2020,1,1,0,0,0,0)

    if 1:
        land_path = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/'
        case_names= ['drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2',
                     'drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB']

        pft       = None #2
        
        # pft      = ["BEF","crop","shrub","grass","barren"]
        # iveg_num = [2, 9, 5, 6, 14]        
        
        var_name  = 'Tmax' #'Tair_f_inst'
        plot_spatial_land_days(land_path,case_names,var_name,pft,time_s,time_e,loc_lat=loc_lat,loc_lon=loc_lon,
                               lat_name=lat_name, lon_name=lon_name, message=message)

        var_name  = 'Tmin' #'Tair_f_inst'
        plot_spatial_land_days(land_path,case_names,var_name,pft,time_s,time_e,loc_lat=loc_lat,loc_lon=loc_lon,
                               lat_name=lat_name, lon_name=lon_name, message=message)



