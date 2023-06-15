#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

'''
Functions:
1. Climdex indices: https://www.climdex.org/learn/indices/
'''

import os
import numpy as np
from osgeo import gdal
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap
from cartopy.feature import NaturalEarthFeature, OCEAN
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from common_utils import *

def tif_to_nc(tif_file_path, output_file, wrf_file):

    # for opening the raster read-only and saving it on f variable.
    f = gdal.Open(tif_file_path, gdal.GA_ReadOnly) 
    print('f',f)

    # Copy the transformation to a variable, it will be useful later.
    GT_input = f.GetGeoTransform()
    print('GT_input',GT_input)

    # Read the bands of your raster using GetRasterBand
    albedo_imp = f.GetRasterBand(1)
    lai_imp    = f.GetRasterBand(2)

    # Transform the read band into an array, to work with it.
    albedo  = albedo_imp.ReadAsArray()
    lai     = lai_imp.ReadAsArray()
    print('albedo',albedo)

    # Read the size of your array
    size1,size2=albedo.shape
    print('size1,size2=',size1,size2)

    Albedo = np.float32(albedo)
    LAI    = np.float32(lai)
    Albedo = np.where(Albedo < 0, -9999. ,Albedo)
    LAI    = np.where(LAI < 0, -9999. ,LAI)
    print('Albedo',Albedo)

    # make nc file
    f                   = Dataset(output_file, 'w', format='NETCDF4')

    ### Create nc file ###
    f.history           = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date     = "%s" % (datetime.now())
    f.description       = 'Transfer '+tif_file_path+' to netcdf file, created by MU Mengyuan'

    # set dimensions
    f.createDimension('north_south', size1)
    f.createDimension('east_west', size2)
    f.Conventions       = "CF-1.0"

    # create variables
    latitude            = f.createVariable('lat', 'f4', ('north_south', 'east_west'))
    latitude.long_name  = "latitude"
    latitude.units      = "degree_north"
    latitude._CoordinateAxisType = "Lat"

    longitude           = f.createVariable('lon', 'f4', ('north_south', 'east_west'))
    longitude.long_name = "longitude"
    longitude.units     = "degree_east"
    longitude._CoordinateAxisType = "Lon"

    lai_out             = f.createVariable('LAI_imp', 'f4', ('north_south', 'east_west'),fill_value=-9999.)
    lai_out.long_name   = "Importance of ΔLAI to ΔTmax"
    lai_out[:]          = LAI[::-1,:]

    alb_out             = f.createVariable('Albedo_imp', 'f4', ('north_south', 'east_west'),fill_value=-9999.)
    alb_out.long_name   = "Importance of Δalbedo to ΔTmax"
    alb_out[:]          = Albedo[::-1,:]
    
    wrf                 = Dataset(wrf_file,'r')
    lat                 = wrf.variables['lat'][:]
    lon                 = wrf.variables['lon'][:]
    latitude[:]         = lat
    longitude[:]        = lon   

    f.close()

def plot_RF_map(output_file,message=None):

    # =================== Plotting spatial map ===================
    f           = Dataset(output_file, 'r')
    LAI_imp     = f.variables['LAI_imp'][:,:] # nseason, nindex, nlat, nlon
    Albedo_imp  = f.variables['Albedo_imp'][:,:]
    lat         = f.variables['lat'][:,:]
    lon         = f.variables['lon'][:,:]
    LAI_imp     = np.where(LAI_imp <0, np.nan, LAI_imp)
    Albedo_imp  = np.where(Albedo_imp <0, np.nan, Albedo_imp)

    # regrid to
    lat_out     = np.arange(-39,-24,0.04)
    lon_out     = np.arange(135,155,0.04)

    lon_out_2D,lat_out_2D = np.meshgrid(lon_out,lat_out)
    nlat        = len(lat_out)
    nlon        = len(lon_out)

    lat_in_1D   = lat.flatten()
    lon_in_1D   = lon.flatten()
    LAI_in_1D_tmp = LAI_imp.flatten()
    ALB_in_1D   = Albedo_imp.flatten()

    LAI_in_1D   = LAI_in_1D_tmp[~np.isnan(LAI_in_1D_tmp)]
    ALB_in_1D   = ALB_in_1D[~np.isnan(LAI_in_1D_tmp)]
    lat_in_1D   = lat_in_1D[~np.isnan(LAI_in_1D_tmp)]    # here I make nan in values as the standard
    lon_in_1D   = lon_in_1D[~np.isnan(LAI_in_1D_tmp)]
    
    LAI_regrid  = griddata((lat_in_1D, lon_in_1D), LAI_in_1D, (lat_out_2D, lon_out_2D), method='nearest')
    ALB_regrid  = griddata((lat_in_1D, lon_in_1D), ALB_in_1D, (lat_out_2D, lon_out_2D), method='nearest')

    var_imp     = np.zeros((nlat,nlon))
    var_imp     = np.where(np.all([ALB_regrid>0.005,      LAI_regrid<0.005],axis=0),1,var_imp) # albedo domainate 
    var_imp     = np.where(np.all([LAI_regrid>0.005,      ALB_regrid<0.005],axis=0),2,var_imp) # LAI domainate 
    var_imp     = np.where(np.all([ALB_regrid>LAI_regrid, LAI_regrid>0.005],axis=0),3,var_imp) # both important by albedo domainate 
    var_imp     = np.where(np.all([LAI_regrid>ALB_regrid, ALB_regrid>0.005],axis=0),4,var_imp) # both important by lai domainate 

    # ____ plotting ____
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=[5,5],sharex=True, sharey=True, squeeze=True,
                            subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=0., hspace=-0.05) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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

    # set the box type of sequence number
    props = dict(boxstyle="round", facecolor='white', alpha=0.0, ec='white')
    # choose colormap

    states= NaturalEarthFeature(category="cultural", scale="50m",
                                        facecolor="none",
                                        name="admin_1_states_provinces_shp")


    # ======================= Set colormap =======================
    cmap    = ListedColormap(["lightgray",     # 0: 
                              "orange",    # 1: albedo domainate 
                              "blue",          # 2: LAI domainate 
                              "lightgreen",    # 3: both important by albedo domainate 
                              "green",         # 4: both important by lai domainate 
                              ]) #plt.cm.tab10 #BrBG

    ax.coastlines(resolution="50m",linewidth=1)
    ax.set_extent([135,155,-39,-24])
    ax.add_feature(states, linewidth=.5, edgecolor="black")

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
    gl.top_labels   = False
    gl.right_labels = False
    gl.bottom_labels= True
    gl.left_labels  = True
    gl.xlines       = False
    gl.ylines       = False
    gl.xlocator     = mticker.FixedLocator(np.arange(125,160,1))
    gl.ylocator     = mticker.FixedLocator(np.arange(-40,-20,1))
    gl.xformatter   = LONGITUDE_FORMATTER
    gl.yformatter   = LATITUDE_FORMATTER
    gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
    gl.ylabel_style = {'size':12, 'color':almost_black}

    extent   = (135, 155, -39, -24)
    plot1    = ax.imshow(var_imp, origin="lower", extent=extent, vmin=-0.5, vmax=4.5, transform=ccrs.PlateCarree(), cmap=cmap)
    ax.add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
    cbar = plt.colorbar(plot1, ax=ax, ticklocation="right", pad=0.01, orientation="vertical",
                        aspect=20, shrink=0.6) # cax=cax,
    cbar.ax.tick_params(labelsize=8, labelrotation=45)

    if message != None:
        plt.savefig('./plots/spatial_map_RF_LAI_ALB_imp_'+message+'.png',dpi=300)
    else:
        plt.savefig('./plots/spatial_map_RF_LAI_ALB_imp.png',dpi=300)

if __name__ == "__main__":

    tif_file_path = '/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/nc_files/VarImp_deltaTmax_deltaLAI_deltaAlbedo_2017_2020_DJF.tif'
    output_file   = '/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/nc_files/VarImp_deltaTmax_deltaLAI_deltaAlbedo_2017_2020_DJF.nc'
    wrf_file      = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2/WRF_output/T2/wrfout_201701-202002.nc"

    # make nc file
    # tif_to_nc(tif_file_path, output_file, wrf_file)

    plot_RF_map(output_file)