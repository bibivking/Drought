#!/usr/bin/env python
__author__    = "MU Mengyuan"

import os
import sys
import glob
import numpy as np
import pandas as pd
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors
from netCDF4 import Dataset
from matplotlib import cm
from matplotlib import ticker
import xarray as xr
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# from convert_units import get_land_var_scale


def plot_map_landinfo_offline_CABLE(file_path, var_names,loc_lat=None,loc_lon=None):
    
    # Loop vars
    for var_name in var_names:
        print(var_name)
        
        # Read gridinfo file 
        var1 = Dataset(file_path[0], mode='r')
        var2 = Dataset(file_path[1], mode='r')    

        Var1 = var1.variables[var_name][:,:]
        Var2 = var2.variables[var_name][:,:]
        Var_d= Var2 - Var1
        
        lons = var1.variables['longitude']
        lats = var1.variables['latitude']
        lon, lat = np.meshgrid(lons, lats)    
        # units = var1.variables[var_name].units
        

        # Plot setting
        fig, ax = plt.subplots(nrows=6, ncols=3, figsize=[10,12],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
        plt.subplots_adjust(wspace=0, hspace=0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

        plt.rcParams['text.usetex']     = False
        plt.rcParams['font.family']     = "sans-serif"
        plt.rcParams['font.serif']      = "Helvetica"
        plt.rcParams['axes.linewidth']  = 1.5
        plt.rcParams['axes.labelsize']  = 14
        plt.rcParams['font.size']       = 14
        plt.rcParams['legend.fontsize'] = 10
        plt.rcParams['xtick.labelsize'] = 14
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

        if var_name in ["iveg","isoil"]:
            clevs   = np.linspace(0.5,15.5, num=16)
            clevs_d = np.linspace(-8,8, num=17)
        elif var_name in ["Albedo","ssat_vec","swilt_vec","sfc_vec","watr"]:
            clevs   = np.linspace(0.0,0.5, num=21)
            clevs_d = np.linspace(-0.2,0.2, num=21)
        elif var_name in ["sand_vec","clay_vec","silt_vec","org_vec"]:
            clevs   = np.linspace(0.0,0.8, num=17)
            clevs_d = np.linspace(-0.4,0.4, num=17)
        elif var_name in ["LAI"]:
            clevs   = np.linspace(0.,8., num=16)
            clevs_d = np.linspace(-8,8., num=17)
        elif var_name in ["bch_vec"]:
            clevs   = np.linspace(0.,13., num=26)
            clevs_d = np.linspace(-5,5., num=11)
        elif var_name in ["hyds_vec"]:
            clevs   = np.linspace(0.,0.02, num=21)
            clevs_d = np.linspace(-0.01,0.01, num=21)
        elif var_name in ["sucs_vec"]:
            clevs   = np.linspace(0.,0.1, num=21)
            clevs_d = np.linspace(-0.1,0.1, num=21)
        elif var_name in ["cnsd_vec"]:
            clevs   = np.linspace(0.2,0.3, num=11)
            clevs_d = np.linspace(-0.1,0.1, num=11)
        elif var_name in ["css_vec"]:
            clevs   = np.linspace(700,900, num=21)
            clevs_d = np.linspace(-200,200, num=21)
        elif var_name in ["org_vec"]:
            clevs   = np.linspace(0.0,0.005, num=11)
            clevs_d = np.linspace(-0.002,0.002, num=11)
        elif var_name in ["rhosoil_vec"]:
            clevs   = np.linspace(500,1500, num=11)
            clevs_d = np.linspace(-1000,1000, num=11)
        else:
            clevs   = np.linspace(np.min(Var1),np.max(Var1), num=21)
            clevs_d = np.linspace(np.max(Var1)*(-1.0),np.max(Var1), num=21)
        
        # Loop column
        for c in np.arange(3):
            if c == 0:
                var_plot = Var1
                title = "AU_NAT"
                cmap    = plt.cm.BrBG
            elif c == 1:
                var_plot = Var2
                title = "SoilGrid"
                cmap    = plt.cm.BrBG
            elif c == 2:
                var_plot = Var_d
                clevs    = clevs_d
                title = "SoilGrid - AU_NAT"
                cmap    = plt.cm.BrBG
                
            # Loop layers
            for l in np.arange(6):
                
                print("start plotting")
                
                if loc_lat == None:
                    ax[l,c].set_extent([135,155,-40,-25])
                else:
                    ax[l,c].set_extent([loc_lon[0],loc_lon[1],loc_lat[0],loc_lat[1]])

                ax[l,c].coastlines(resolution="50m",linewidth=1)
                plot1 = ax[l,c].contourf(lon, lat, var_plot[l,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
            
            ax[0,c].set_title(title)
            cbar = plt.colorbar(plot1, ax=ax[:,c], ticklocation="right",aspect=40, pad=0.01, orientation="vertical") #, shrink=0.06
        # cbar.set_label(var_name,size=10) #,labelpad=15
        plt.title(var_name)#+' '+units)
        plt.savefig('./plots/spatial_map_AU_NAT_vs_Openland_'+var_name+'.png',dpi=1200)

        var1  = None
        var2  = None
        Var1  = None
        Var2  = None
        Var_d = None
        
if __name__ == "__main__":

    #############################################
    #   plot plot_map_landinfo_offline_CABLE    # 
    #############################################
    region    = "Aus"
    file_path = ["/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_CSIRO_AU_NAT_ELEV_DLCM_fix.nc",                 
                 "/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc"]  

    var_names  = ["bch_vec","sucs_vec","ssat_vec","swilt_vec","sfc_vec","css_vec","hyds_vec",
                  "rhosoil_vec","cnsd_vec","sand_vec","clay_vec","silt_vec","org_vec","watr"]
    
    
    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-25]
        loc_lon    = [135,155]
    elif region == "CORDEX":
        loc_lat    = [-52.36,3.87]
        loc_lon    = [89.25,180]
        
    plot_map_landinfo_offline_CABLE(file_path,var_names,loc_lat,loc_lon)






