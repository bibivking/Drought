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


def plot_map_landinfo_offline_CABLE(file_path, var_names):

    # Open the NetCDF4 file (add a directory path if necessary) for reading:
    var = Dataset(file_path, mode='r')

    lons = var.variables['longitude']
    lats = var.variables['latitude']
    lon, lat = np.meshgrid(lons, lats)

    for var_name in var_names:

        print(var_name)
        Var   = var.variables[var_name]
        clevs = np.linspace(0.0,1., num=11)

        # Make plots
        fig, ax = plt.subplots(nrows=6, ncols=1, figsize=[5,10],sharex=True, sharey=True, squeeze=True, subplot_kw={'projection': ccrs.PlateCarree()})
 

        # fig = plt.figure(figsize=(9,5))
        # ax = plt.axes(projection=ccrs.PlateCarree())


        for l in np.arange(6):
            # ax[l].set_extent([110,155,-45,-10])
            ax[l].coastlines(resolution="50m",linewidth=1)
            # Add gridlines
            # # gl = ax[l].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
            # gl.xlabels_top  = False
            # gl.ylabels_right= False
            # gl.xlines       = True
            # gl.xlocator     = mticker.FixedLocator([110,115,120,125,130,135,140,145,150,155])
            # gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15,-10])
            # gl.xformatter   = LONGITUDE_FORMATTER
            # gl.yformatter   = LATITUDE_FORMATTER
            # gl.xlabel_style = {'size':10, 'color':'black'}
            # gl.ylabel_style = {'size':10, 'color':'black'}

        cmap = plt.cm.jet

        plot1 = ax[0].contourf(lon, lat, Var[0,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        plot2 = ax[1].contourf(lon, lat, Var[1,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        plot3 = ax[2].contourf(lon, lat, Var[2,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        plot4 = ax[3].contourf(lon, lat, Var[3,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        plot5 = ax[4].contourf(lon, lat, Var[4,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)
        plot6 = ax[5].contourf(lon, lat, Var[5,:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap)

        cb1 = plt.colorbar(plot1, ax=ax[0], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb2 = plt.colorbar(plot2, ax=ax[1], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb3 = plt.colorbar(plot3, ax=ax[2], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb4 = plt.colorbar(plot4, ax=ax[3], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb5 = plt.colorbar(plot5, ax=ax[4], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb6 = plt.colorbar(plot6, ax=ax[5], orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

        # units_string = var.variables[var_name].units

        # plt.title(var_name, size=16)

        # cb.ax[l].tick_params(labelsize=10)
        # cb.set_label(units_string,size=14,rotation=270,labelpad=15)

        plt.savefig('spatial_map_AU_NAT_C_'+var_name+'.png',dpi=300)

        Var = None

if __name__ == "__main__":

    #############################################
    #   plot plot_map_landinfo_offline_CABLE    # 
    #############################################
    file_path = "/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/gridinfo_AWAP_CSIRO_AU_NAT_new_iveg.nc"
    # var_names  = ["iveg","isoil","sand","clay","silt",
    #               "sfc","ssat","swilt","hyds","bch","sucs",
    #               "LAI","Albedo",'elevation']
    var_names = ["clay_vec","silt_vec","sand_vec"]
    plot_map_landinfo_offline_CABLE(file_path,var_names)






