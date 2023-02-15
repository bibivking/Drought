import sys
import cartopy
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from scipy.interpolate import griddata
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from common_utils import *

def read_off_wb(file_path,offline_path):
    '''
    read off wb
    '''
    print("read_off_wb")
    #"SE Aus":
    loc_lat         = [-40,-25]
    loc_lon         = [135,155]
     
    time_s          = datetime(2017,1,1,0,0,0,0)
    time_e          = datetime(2017,1,1,23,59,0,0)

    time, wb_tmp    = read_var(offline_path, 'SoilMoist', loc_lat, loc_lon, 'latitude', 'longitude')
    wb              = spital_var(time,wb_tmp,time_s,time_e)
    print(time)
    time, GWwb_tmp  = read_var(offline_path, 'GWMoist', loc_lat, loc_lon, 'latitude', 'longitude')
    GWwb            = spital_var(time,GWwb_tmp,time_s,time_e)    
    
    # read lat lon in
    offline         = Dataset(offline_path, mode='r')
    lat_in          = offline.variables["latitude"][:]  
    lon_in          = offline.variables["longitude"][:]  
    
    # read lat lon out
    file            = Dataset(file_path, mode='r')
    lat_out         = file.variables["lat"][:,:]
    lon_out         = file.variables["lon"][:,:]
    nlat            = len(lat_out[:,0])
    nlon            = len(lon_out[0,:])
    nsoil           = 6
    
    # define
    wb_regrid       = np.zeros((nsoil,nlat,nlon))
    GWwb_regrid     = np.zeros((nlat,nlon))
    
    for l in np.arange(6):
        wb_regrid[l,:,:] = regrid_data(lat_in, lon_in, lat_out, lon_out, wb[l,:,:]) 
    
    GWwb_regrid = regrid_data(lat_in, lon_in, lat_out, lon_out, GWwb)
        
    return (wb_regrid,GWwb_regrid)

def read_rst_wb(file_path, lis_rst_path):
    
    '''
    read rst wb
    '''
    print("read_rst_wb")
    
    den_rat     = 0.921
    
    lis_rst     = Dataset(lis_rst_path, mode='r')
    wbliq       = lis_rst.variables["WB"][:,:]
    wbice       = lis_rst.variables["WBICE"][:,:]
    wb          = wbliq + wbice*den_rat
    GWwb        = lis_rst.variables["GWWB"][:]
    lat_in      = lis_rst.variables["lat"][:]  
    lon_in      = lis_rst.variables["lon"][:]  
    
    file        = Dataset(file_path, mode='r')
    lat_out     = file.variables["lat"][:,:]
    lon_out     = file.variables["lon"][:,:]
    nlat        = len(lat_out[:,0])
    nlon        = len(lon_out[0,:])
    nsoil       = 6
    
    wb_regrid   = np.zeros((nsoil,nlat,nlon))
    
    for l in np.arange(6):
        wb_regrid[l,:,:] = griddata((lon_in, lat_in), wb[l,:], (lon_out, lat_out), method="linear")
        
    GWwb_regrid = griddata((lon_in, lat_in), GWwb, (lon_out, lat_out), method="linear")
        
    return (wb_regrid, GWwb_regrid)

def read_lis_wb(file_path):
    
    '''
    read lis wb
    '''
    print("read_lis_wb")
    
    lis_out     = Dataset(file_path, mode='r')
    wb          = lis_out.variables["SoilMoist_inst"][0,:,:,:]
    GWwb        = lis_out.variables["GWwb_tavg"][0,:,:]
        
    return (wb, GWwb)

def spatial_map_single_plot_diff(file_path, lis_rst_path, offline_path, wrf_path):
    
    wb_rst, GWwb_rst = read_rst_wb(file_path, lis_rst_path)
    wb_off, GWwb_off = read_off_wb(file_path, offline_path)
    wb_lis, GWwb_lis = read_lis_wb(file_path)
    
    # wb_diff          = wb_rst - wb_off
    # GWwb_diff        = GWwb_rst - GWwb_off
    
    # wb_diff          = wb_lis - wb_off
    # GWwb_diff        = GWwb_lis - GWwb_off
    
    wb_diff          = wb_lis - wb_rst
    GWwb_diff        = GWwb_lis - GWwb_rst

    # read lat and lon outs
    wrf              = Dataset(wrf_path,  mode='r')
    lons             = wrf.variables['XLONG'][0,:,:]
    lats             = wrf.variables['XLAT'][0,:,:]

    # =========================== plot wb ===========================
    for l in np.arange(6):
        fig = plt.figure(figsize=(6,5))
        ax  = plt.axes(projection=ccrs.PlateCarree())

        plt.rcParams['text.usetex']     = False
        plt.rcParams['font.family']     = "sans-serif"
        plt.rcParams['font.serif']      = "Helvetica"
        plt.rcParams['axes.linewidth']  = 1.5
        plt.rcParams['axes.labelsize']  = 14
        plt.rcParams['font.size']       = 14
        plt.rcParams['legend.fontsize'] = 14
        plt.rcParams['xtick.labelsize'] = 14
        plt.rcParams['ytick.labelsize'] = 14

        almost_black                    = '#262626'
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

        # color bar
        cmap  = plt.cm.seismic

        # start plotting
        ax.set_extent([135,155,-40,-25])
        ax.coastlines(resolution="50m",linewidth=1)

        # Add gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
        gl.xlabels_top   = False
        gl.ylabels_right = False
        gl.xlines        = True

        gl.xlocator      = mticker.FixedLocator([135,140,145,150,155])
        gl.ylocator      = mticker.FixedLocator([-40,-35,-30,-25])

        gl.xformatter    = LONGITUDE_FORMATTER
        gl.yformatter    = LATITUDE_FORMATTER
        gl.xlabel_style  = {'size':10, 'color':'black'}
        gl.ylabel_style  = {'size':10, 'color':'black'}

        clevs = [-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2,0.25,0.3]

        plt.contourf(lons, lats, wb_diff[l,:,:], clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') # 

        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        cb.ax.tick_params(labelsize=10)
        plt.title("wb_diff_lyr="+str(l), size=16)

        # message = "rst-off_wb_lyr="+str(l)
        message = "lis-rst_wb_lyr="+str(l)
        
        plt.savefig('/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/plots/WTD_sudden_change/spatial_map_'+message+'.png',dpi=300)
        
    # ========================== plot GWwb ==========================
    fig = plt.figure(figsize=(6,5))
    ax  = plt.axes(projection=ccrs.PlateCarree())

    plt.rcParams['text.usetex']     = False
    plt.rcParams['font.family']     = "sans-serif"
    plt.rcParams['font.serif']      = "Helvetica"
    plt.rcParams['axes.linewidth']  = 1.5
    plt.rcParams['axes.labelsize']  = 14
    plt.rcParams['font.size']       = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    almost_black                    = '#262626'
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

    # color bar
    cmap  = plt.cm.seismic

    # start plotting
    ax.set_extent([135,155,-40,-25])
    ax.coastlines(resolution="50m",linewidth=1)

    # Add gridlines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='black', linestyle='--')
    gl.xlabels_top   = False
    gl.ylabels_right = False
    gl.xlines        = True

    gl.xlocator      = mticker.FixedLocator([135,140,145,150,155])
    gl.ylocator      = mticker.FixedLocator([-40,-35,-30,-25])

    gl.xformatter    = LONGITUDE_FORMATTER
    gl.yformatter    = LATITUDE_FORMATTER
    gl.xlabel_style  = {'size':10, 'color':'black'}
    gl.ylabel_style  = {'size':10, 'color':'black'}

    clevs = [-0.05,-0.04,-0.03,-0.02,-0.01,-0.005,0.005,0.01,0.02,0.03,0.04,0.05]

    plt.contourf(lons, lats, GWwb_diff, clevs, transform=ccrs.PlateCarree(), cmap=cmap, extend='both') # 

    cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
    cb.ax.tick_params(labelsize=10)
    plt.title("GWwb_diff", size=16)

    # message = "rst-off_GWwb"
    message = "lis-rst_GWwb"
    
    plt.savefig('/g/data/w97/mm3972/scripts/Drought/drght_2017-2019/plots/WTD_sudden_change/spatial_map_'+message+'.png',dpi=300)        

if __name__ == "__main__":
    
    file_path      = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/LIS_output/LIS.CABLE.201701-201701.d01.nc"
    wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/uniform_soil_param/drght_2017_2019/run_Jan2017/WRF_output/wrfout_d01_2017-01-01_11:00:00"
    lis_rst_path   = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/offline_rst_output/output_1719_drght/LIS_RST_CABLE_201701011100_restart_2017_Day2.d01.nc"
    offline_path   = "/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc"
    
    spatial_map_single_plot_diff(file_path, lis_rst_path, offline_path, wrf_path)
