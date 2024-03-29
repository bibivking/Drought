#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

'''
Plot spitial map of land diagnosis and parameters from LIS-CABLE
1. per time step
2. time period average
'''

import numpy as np
import shapefile as shp
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from netCDF4 import Dataset
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature, OCEAN
from mpl_toolkits.axes_grid1 import make_axes_locatable
from convert_units import get_land_var_scale, get_land_var_range_diff
from common_utils import *

def plot_spatial_fwsoil_qle_qh_tmax_qair(land_20181201_path,land_20191201_path, land_clim_2018_path, land_clim_2019_path,
                                        wrf_path, shape_path = None, loc_lat=None, loc_lon=None, seconds=None, message=None):

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=4, ncols=5, figsize=[16,10],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=-0.6, hspace=0.0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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
    # ocean = NaturalEarthFeature('cultural', 'ocean', scale='50m',
    #                         edgecolor='none', facecolor="lightgray")

    # ======================= Set colormap =======================
    cmap_1     = plt.cm.seismic
    blue2white = truncate_colormap(cmap_1, minval=0., maxval=0.5)
    white2red  = truncate_colormap(cmap_1, minval=0.5, maxval=1.)
    cmap_2     = plt.cm.seismic_r
    red2white  = truncate_colormap(cmap_2, minval=0., maxval=0.5)
    white2blue = truncate_colormap(cmap_2, minval=0.5, maxval=1.)

    # ======================= Read WRF file =======================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    texts   = [ "(a)","(b)","(c)","(d)","(e)",
                "(f)","(g)","(h)","(i)","(j)",
                "(k)","(l)","(m)","(n)","(o)",
                "(p)","(q)","(r)","(s)","(t)"]

    label_x = ["Δ$β$",
               "ΔQ$\mathregular{_e}$",
               "ΔQ$\mathregular{_h}$",
               "ΔT$\mathregular{_{max}}$",
               "Δq"]

    label_y = ["2018/2019","2019/2020","2018/2019","2019/2020"]
    loc_y   = [0.5,0.5,0.5,0.5]
    #[0.62,0.55,0.48,0.42]

    cnt     = 0

    # ==================== Set up files ====================
    for i in np.arange(4):
        if i == 0:
            time_s = datetime(2018,12,1,0,0,0,0)
            time_e = datetime(2019,2,28,23,59,0,0)
            year_s = 2018
            land_1201_path = land_20181201_path
            land_clim_path = land_clim_2018_path
        elif i == 1:
            time_s = datetime(2019,1,14,0,0,0,0)
            time_e = datetime(2019,1,26,23,59,0,0)
            year_s = 2018
            land_1201_path = land_20181201_path
            land_clim_path = land_clim_2018_path
        elif i == 2:
            time_s = datetime(2019,12,1,0,0,0,0)
            time_e = datetime(2020,2,29,23,59,0,0)
            year_s = 2019
            land_1201_path = land_20191201_path
            land_clim_path = land_clim_2019_path
        elif i == 3:
            time_s = datetime(2019,12,16,0,0,0,0)
            time_e = datetime(2020,1,7,23,59,0,0)
            year_s = 2019
            land_1201_path = land_20191201_path
            land_clim_path = land_clim_2019_path

        land_1201_files = [ land_1201_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02.d01.nc",]
        land_clim_files = [ land_clim_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02.d01.nc",]

        time, Dec1_T    = read_var_multi_file(land_1201_files, "Tair_f_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_T    = read_var_multi_file(land_clim_files, "Tair_f_inst", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qle  = read_var_multi_file(land_1201_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qle  = read_var_multi_file(land_clim_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_FWsoil = read_var_multi_file(land_1201_files, "FWsoil_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_FWsoil = read_var_multi_file(land_clim_files, "FWsoil_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qh   = read_var_multi_file(land_1201_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qh   = read_var_multi_file(land_clim_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qair = read_var_multi_file(land_1201_files, "Qair_f_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qair = read_var_multi_file(land_clim_files, "Qair_f_inst", loc_lat, loc_lon, "lat", "lon")

        print("time = ", time)
        dec1_Tmax       = spatial_var_max(time,Dec1_T,time_s,time_e)
        clim_Tmax       = spatial_var_max(time,Clim_T,time_s,time_e)

        dec1_Qle        = spatial_var(time,Dec1_Qle,time_s,time_e)
        clim_Qle        = spatial_var(time,Clim_Qle,time_s,time_e)

        dec1_FWsoil     = spatial_var(time,Dec1_FWsoil,time_s,time_e)
        clim_FWsoil     = spatial_var(time,Clim_FWsoil,time_s,time_e)

        dec1_Qh         = spatial_var(time,Dec1_Qh,time_s,time_e)
        clim_Qh         = spatial_var(time,Clim_Qh,time_s,time_e)

        dec1_Qair       = spatial_var(time,Dec1_Qair,time_s,time_e)
        clim_Qair       = spatial_var(time,Clim_Qair,time_s,time_e)


        # "Tmax":
        tmax_diff       = dec1_Tmax-clim_Tmax

        # "Qle":
        qle_diff        = dec1_Qle-clim_Qle

        # "FWsoil":
        fw_diff         = dec1_FWsoil-clim_FWsoil

        # "Qh":
        qh_diff         = dec1_Qh-clim_Qh

        # "Qair":
        q_diff          = dec1_Qair-clim_Qair


        # ==================== Start to plot ====================
        for j in np.arange(5):

            ax[i,j].coastlines(resolution="50m",linewidth=1)
            ax[i,j].set_extent([135,155,-40,-23.6])
            ax[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            if j == 0:
                gl.ylabels_left   = True
            else:
                gl.ylabels_left   = False
            if i == 3:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            # set y label
            if j == 0:
                ax[i,j].text(-0.255, loc_y[i], label_y[i], va='bottom', ha='center',
                              rotation='vertical', rotation_mode='anchor',
                              transform=ax[i,j].transAxes)
                if i == 0:
                    ax[i,j].text(-0.36, 0, "Summer", va='bottom', ha='center',
                                  rotation='vertical', rotation_mode='anchor',fontsize=14,
                                  transform=ax[i,j].transAxes)
                    ax[i,j].text(-0.36, -2, "Heatwave Period", va='bottom', ha='center',
                                  rotation='vertical', rotation_mode='anchor',fontsize=14,
                                  transform=ax[i,j].transAxes)

        # left - FWsoil
        clevs1  = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5] # -1,-0.9,-0.8,-0.7,-0.6, ,0.6,0.7,0.8,0.9,1
        plot1   = ax[i,0].contourf(lon, lat, fw_diff, levels=clevs1, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        cbar    = ax[i,0].text(0.02, 0.15, texts[cnt], transform=ax[i,0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray") # lightgray

        # middle left - Qle
        clevs2   = [ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100] #-140, -120, , 120, 140
        plot2    = ax[i,1].contourf(lon, lat, qle_diff, levels=clevs2, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        ax[i,1].text(0.02, 0.15, texts[cnt+1], transform=ax[i,1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle - Qh
        clevs3  = [ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100]
        # [ -140, -120, -100, -80, -60, -40, -20, -10]
        plot3   = ax[i,2].contourf(lon, lat, qh_diff, levels=clevs3, transform=ccrs.PlateCarree(),cmap=cmap_1,extend='both') #
        ax[i,2].text(0.02, 0.15, texts[cnt+2], transform=ax[i,2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle right - Tmax
        tmax_diff= np.where(np.isnan(fw_diff), np.nan, tmax_diff)
        clevs4   = [-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.25, 0.25, 0.5, 1., 1.5, 2., 2.5, 3.]
        plot4    = ax[i,3].contourf(lon, lat, tmax_diff, levels=clevs4, transform=ccrs.PlateCarree(),cmap=cmap_1,extend='both') #
        ax[i,3].text(0.02, 0.15, texts[cnt+3], transform=ax[i,3].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,3].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right - Qair
        clevs5  = [-2,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2.]
        plot5   = ax[i,4].contourf(lon, lat, q_diff*1000., levels=clevs5, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        ax[i,4].text(0.02, 0.15, texts[cnt+4], transform=ax[i,4].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,4].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
        clevs   = None

        if shape_path != None:
        # Requires the pyshp package
            sf = shp.Reader(shape_path)

            for shape in sf.shapeRecords():
                x = [n[0] for n in shape.shape.points[:]]
                y = [n[1] for n in shape.shape.points[:]]
                ax[i,0].plot(x,y,c="black")
                ax[i,1].plot(x,y,c="black")
                ax[i,2].plot(x,y,c="black")
                ax[i,3].plot(x,y,c="black")
                ax[i,4].plot(x,y,c="black")

        # set top x label
        if i == 0:
            ax[i,0].set_title(label_x[0])#, fontsize=12)
            ax[i,1].set_title(label_x[1])#, fontsize=12)
            ax[i,2].set_title(label_x[2])#, fontsize=12)
            ax[i,3].set_title(label_x[3])#, fontsize=12)
            ax[i,4].set_title(label_x[4])#, fontsize=12)


        # set bottom colorbar
        if i == 3:

            # left - Fwsoil
            cbar = plt.colorbar(plot1, ax=ax[:,0], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "-"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
            # cbar.ax.set_yticklabels(["-0.5","-0.4","-0.3","-0.2","-0.1","0.1","0.2","0.3","0.4","0.5"])

            # middle left - Qle
            cbar = plt.colorbar(plot2, ax=ax[:,1], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100])
            # cbar.ax.set_yticklabels(["-100","-80","-60","-40","-20","-5","5","20","40","60","80","100"])

            # middle - Qh
            cbar = plt.colorbar(plot3, ax=ax[:,2], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100])
            # cbar.ax.set_yticklabels(["-100","-80","-60","-40","-20","-5","5","20","40","60","80","100"])

            # middle right - Tmax
            cbar = plt.colorbar(plot4, ax=ax[:,3], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "$\mathregular{^o}$C"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.25, 0.25, 0.5, 1., 1.5, 2., 2.5, 3.])
            # cbar.ax.set_yticklabels(["-3","-2.5","-2","-1.5","-1","-0.5","-0.25","0.25","0.5","1.","1.5","2.","2.5","3."])

            # right - q
            cbar = plt.colorbar(plot5, ax=ax[:,4], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "g kg$\mathregular{^{-1}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-2,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2.])
            # cbar.ax.set_yticklabels(["-2"," ","-1.5"," ","-1."," ","-0.5"," "," ","0.5"," ","1"," ","1.5"," ","2"])
        cnt = cnt + 5

    plt.savefig('./plots/spatial_map_fwsoil_Qle_Qh_Tmax_q_2018_2020_summer_HW_version2.pdf',dpi=300)

def plot_spatial_fwsoil_qle_qh_tmax_qrate(land_20181201_path,land_20191201_path, land_clim_2018_path, land_clim_2019_path,
                                        wrf_path, shape_path = None, loc_lat=None, loc_lon=None, seconds=None, message=None):

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=4, ncols=5, figsize=[16,10],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=-0.65, hspace=0.005) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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
    # ocean = NaturalEarthFeature('cultural', 'ocean', scale='50m',
    #                         edgecolor='none', facecolor="lightgray")

    # ======================= Set colormap =======================
    cmap_1     = plt.cm.seismic
    blue2white = truncate_colormap(cmap_1, minval=0., maxval=0.5)
    white2red  = truncate_colormap(cmap_1, minval=0.5, maxval=1.)
    cmap_2     = plt.cm.seismic_r
    red2white  = truncate_colormap(cmap_2, minval=0., maxval=0.5)
    white2blue = truncate_colormap(cmap_2, minval=0.5, maxval=1.)

    # ======================= Read WRF file =======================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    texts   = [ "(a)","(b)","(c)","(d)","(e)",
                "(f)","(g)","(h)","(i)","(j)",
                "(k)","(l)","(m)","(n)","(o)",
                "(p)","(q)","(r)","(s)","(t)"]

    label_x = ["Δ$β$",
               "ΔQ$\mathregular{_e}$",
               "ΔQ$\mathregular{_h}$",
               "ΔT$\mathregular{_{max}}$",
               "Δq"]

    label_y = ["2018/2019","2019/2020","2018/2019","2019/2020"]
    loc_y   = [0.5,0.5,0.5,0.5]
    #[0.62,0.55,0.48,0.42]

    cnt     = 0

    # ==================== Set up files ====================
    for i in np.arange(4):
        if i == 0:
            time_s = datetime(2018,12,1,0,0,0,0)
            time_e = datetime(2019,2,28,23,59,0,0)
            year_s = 2018
            land_1201_path = land_20181201_path
            land_clim_path = land_clim_2018_path
        elif i == 1:
            time_s = datetime(2019,1,14,0,0,0,0)
            time_e = datetime(2019,1,26,23,59,0,0)
            year_s = 2018
            land_1201_path = land_20181201_path
            land_clim_path = land_clim_2018_path
        elif i == 2:
            time_s = datetime(2019,12,1,0,0,0,0)
            time_e = datetime(2020,2,29,23,59,0,0)
            year_s = 2019
            land_1201_path = land_20191201_path
            land_clim_path = land_clim_2019_path
        elif i == 3:
            time_s = datetime(2019,12,16,0,0,0,0)
            time_e = datetime(2020,1,7,23,59,0,0)
            year_s = 2019
            land_1201_path = land_20191201_path
            land_clim_path = land_clim_2019_path

        land_1201_files = [ land_1201_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12_select.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01_select.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02_select.d01.nc",]
        land_clim_files = [ land_clim_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12_select.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01_select.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02_select.d01.nc",]

        time, Dec1_T    = read_var_multi_file(land_1201_files, "Tair_f_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_T    = read_var_multi_file(land_clim_files, "Tair_f_inst", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qle  = read_var_multi_file(land_1201_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qle  = read_var_multi_file(land_clim_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_FWsoil = read_var_multi_file(land_1201_files, "FWsoil_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_FWsoil = read_var_multi_file(land_clim_files, "FWsoil_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qh   = read_var_multi_file(land_1201_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qh   = read_var_multi_file(land_clim_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qair = read_var_multi_file(land_1201_files, "Qair_f_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qair = read_var_multi_file(land_clim_files, "Qair_f_inst", loc_lat, loc_lon, "lat", "lon")

        print("time = ", time)
        dec1_Tmax       = spatial_var_max(time,Dec1_T,time_s,time_e)
        clim_Tmax       = spatial_var_max(time,Clim_T,time_s,time_e)

        dec1_Qle        = spatial_var(time,Dec1_Qle,time_s,time_e)
        clim_Qle        = spatial_var(time,Clim_Qle,time_s,time_e)

        dec1_FWsoil     = spatial_var(time,Dec1_FWsoil,time_s,time_e)
        clim_FWsoil     = spatial_var(time,Clim_FWsoil,time_s,time_e)

        dec1_Qh         = spatial_var(time,Dec1_Qh,time_s,time_e)
        clim_Qh         = spatial_var(time,Clim_Qh,time_s,time_e)

        dec1_Qair       = spatial_var(time,Dec1_Qair,time_s,time_e)
        clim_Qair       = spatial_var(time,Clim_Qair,time_s,time_e)


        # "Tmax":
        tmax_diff       = dec1_Tmax-clim_Tmax

        # "Qle":
        qle_diff        = dec1_Qle-clim_Qle

        # "FWsoil":
        fw_diff         = dec1_FWsoil-clim_FWsoil

        # "Qh":
        qh_diff         = dec1_Qh-clim_Qh

        # "Qair":
        q_diff          = (dec1_Qair-clim_Qair)/clim_Qair*100.


        # ==================== Start to plot ====================
        for j in np.arange(5):

            ax[i,j].coastlines(resolution="50m",linewidth=1)
            ax[i,j].set_extent([135,155,-40,-23.6])
            ax[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            if j == 0:
                gl.ylabels_left   = True
            else:
                gl.ylabels_left   = False
            if i == 3:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            # set y label
            if j == 0:
                ax[i,j].text(-0.255, loc_y[i], label_y[i], va='bottom', ha='center',
                              rotation='vertical', rotation_mode='anchor',
                              transform=ax[i,j].transAxes)
                if i == 0:
                    ax[i,j].text(-0.36, 0, "Summer", va='bottom', ha='center',
                                  rotation='vertical', rotation_mode='anchor',fontsize=14,
                                  transform=ax[i,j].transAxes)
                    ax[i,j].text(-0.36, -2, "Heatwave Period", va='bottom', ha='center',
                                  rotation='vertical', rotation_mode='anchor',fontsize=14,
                                  transform=ax[i,j].transAxes)

        # left - FWsoil
        clevs1  = [-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5] # -1,-0.9,-0.8,-0.7,-0.6, ,0.6,0.7,0.8,0.9,1
        plot1   = ax[i,0].contourf(lon, lat, fw_diff, levels=clevs1, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        cbar    = ax[i,0].text(0.02, 0.15, texts[cnt], transform=ax[i,0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,0].add_feature(OCEAN,edgecolor='none', facecolor="white") #, facecolor="lightgray"# lightgray

        # middle left - Qle
        clevs2   = [ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100] #-140, -120, , 120, 140
        plot2    = ax[i,1].contourf(lon, lat, qle_diff, levels=clevs2, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        ax[i,1].text(0.02, 0.15, texts[cnt+1], transform=ax[i,1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,1].add_feature(OCEAN,edgecolor='none', facecolor="white")

        # middle - Qh
        clevs3  = [ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100]
        # [ -140, -120, -100, -80, -60, -40, -20, -10]
        plot3   = ax[i,2].contourf(lon, lat, qh_diff, levels=clevs3, transform=ccrs.PlateCarree(),cmap=cmap_1,extend='both') #
        ax[i,2].text(0.02, 0.15, texts[cnt+2], transform=ax[i,2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,2].add_feature(OCEAN,edgecolor='none', facecolor="white") #, facecolor="lightgray"

        # middle right - Tmax
        tmax_diff= np.where(np.isnan(fw_diff), np.nan, tmax_diff)
        clevs4   = [-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.25, 0.25, 0.5, 1., 1.5, 2., 2.5, 3.]
        plot4    = ax[i,3].contourf(lon, lat, tmax_diff, levels=clevs4, transform=ccrs.PlateCarree(),cmap=cmap_1,extend='both') #
        ax[i,3].text(0.02, 0.15, texts[cnt+3], transform=ax[i,3].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,3].add_feature(OCEAN,edgecolor='none', facecolor="white") #, facecolor="lightgray"

        # right - Qair rate
        clevs5  = [-16,-14,-12,-10,-8,-6,-4,-2,2,4,6,8,10,12,14,16]
        plot5   = ax[i,4].contourf(lon, lat, q_diff, levels=clevs5, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        ax[i,4].text(0.02, 0.15, texts[cnt+4], transform=ax[i,4].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,4].add_feature(OCEAN,edgecolor='none', facecolor="white") #, facecolor="lightgray"
        clevs   = None

        if shape_path != None:
        # Requires the pyshp package
            sf = shp.Reader(shape_path)

            for shape in sf.shapeRecords():
                x = [n[0] for n in shape.shape.points[:]]
                y = [n[1] for n in shape.shape.points[:]]
                ax[i,0].plot(x,y,c="black")
                ax[i,1].plot(x,y,c="black")
                ax[i,2].plot(x,y,c="black")
                ax[i,3].plot(x,y,c="black")
                ax[i,4].plot(x,y,c="black")

        # set top x label
        if i == 0:
            ax[i,0].set_title(label_x[0])#, fontsize=12)
            ax[i,1].set_title(label_x[1])#, fontsize=12)
            ax[i,2].set_title(label_x[2])#, fontsize=12)
            ax[i,3].set_title(label_x[3])#, fontsize=12)
            ax[i,4].set_title(label_x[4])#, fontsize=12)


        # set bottom colorbar
        if i == 3:

            # left - Fwsoil
            cbar = plt.colorbar(plot1, ax=ax[:,0], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.32)
            color_label= "-"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
            # cbar.ax.set_yticklabels(["-0.5","-0.4","-0.3","-0.2","-0.1","0.1","0.2","0.3","0.4","0.5"])

            # middle left - Qle
            cbar = plt.colorbar(plot2, ax=ax[:,1], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.32)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100])
            # cbar.ax.set_yticklabels(["-100","-80","-60","-40","-20","-5","5","20","40","60","80","100"])

            # middle - Qh
            cbar = plt.colorbar(plot3, ax=ax[:,2], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.32)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100])
            # cbar.ax.set_yticklabels(["-100","-80","-60","-40","-20","-5","5","20","40","60","80","100"])

            # middle right - Tmax
            cbar = plt.colorbar(plot4, ax=ax[:,3], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.32)
            color_label= "$\mathregular{^o}$C"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.25, 0.25, 0.5, 1., 1.5, 2., 2.5, 3.])
            # cbar.ax.set_yticklabels(["-3","-2.5","-2","-1.5","-1","-0.5","-0.25","0.25","0.5","1.","1.5","2.","2.5","3."])

            # right - q
            cbar = plt.colorbar(plot5, ax=ax[:,4], ticklocation="right", pad=0.05, orientation="horizontal",
                                aspect=20, shrink=0.32)
            color_label= "%"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-2,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2.])
            # cbar.ax.set_yticklabels(["-2"," ","-1.5"," ","-1."," ","-0.5"," "," ","0.5"," ","1"," ","1.5"," ","2"])
        cnt = cnt + 5

    plt.savefig('./plots/spatial_map_fwsoil_Qle_Qh_Tmax_qrate_2018_2020_summer_HW_version2.pdf',dpi=300)



def plot_spatial_sm_qle_qh_tmax_qair(land_20181201_path,land_20191201_path, land_clim_2018_path,land_clim_2019_path,
                                     wrf_path, shape_path = None, loc_lat=None, loc_lon=None, seconds=None, message=None):

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=4, ncols=5, figsize=[16,10],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=-0.6, hspace=0.0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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
    # ocean = NaturalEarthFeature('cultural', 'ocean', scale='50m',
    #                         edgecolor='none', facecolor="lightgray")

    # ======================= Set colormap =======================
    cmap_1     = plt.cm.seismic
    blue2white = truncate_colormap(cmap_1, minval=0., maxval=0.5)
    white2red  = truncate_colormap(cmap_1, minval=0.5, maxval=1.)
    cmap_2     = plt.cm.seismic_r
    red2white  = truncate_colormap(cmap_2, minval=0., maxval=0.5)
    white2blue = truncate_colormap(cmap_2, minval=0.5, maxval=1.)

    # ======================= Read WRF file =======================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    texts   = [ "(a)","(b)","(c)","(d)","(e)",
                "(f)","(g)","(h)","(i)","(j)",
                "(k)","(l)","(m)","(n)","(o)",
                "(p)","(q)","(r)","(s)","(t)"]

    label_x = ["ΔSM$\mathregular{_2m}$",
               "ΔQ$\mathregular{_e}$",
               "ΔQ$\mathregular{_h}$",
               "ΔT$\mathregular{_{max}}$",
               "Δq"]

    label_y = ["2018/2019","2019/2020","2018/2019","2019/2020"]
    loc_y   = [0.5,0.5,0.5,0.5]
    #[0.62,0.55,0.48,0.42]

    cnt     = 0

    # ==================== Set up files ====================
    for i in np.arange(4):
        if i == 0:
            time_s = datetime(2018,12,1,0,0,0,0)
            time_e = datetime(2019,2,28,23,59,0,0)
            year_s = 2018
            land_1201_path = land_20181201_path
            land_clim_path = land_clim_2018_path
        elif i == 1:
            time_s = datetime(2019,1,14,0,0,0,0)
            time_e = datetime(2019,1,26,23,59,0,0)
            year_s = 2018
            land_1201_path = land_20181201_path
            land_clim_path = land_clim_2018_path
        elif i == 2:
            time_s = datetime(2019,12,1,0,0,0,0)
            time_e = datetime(2020,2,29,23,59,0,0)
            year_s = 2019
            land_1201_path = land_20191201_path
            land_clim_path = land_clim_2019_path
        elif i == 3:
            time_s = datetime(2019,12,16,0,0,0,0)
            time_e = datetime(2020,1,7,23,59,0,0)
            year_s = 2019
            land_1201_path = land_20191201_path
            land_clim_path = land_clim_2019_path

        land_1201_files = [ land_1201_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02.d01.nc",]
        land_clim_files = [ land_clim_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02.d01.nc",]

        time, Dec1_T    = read_var_multi_file(land_1201_files, "Tair_f_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_T    = read_var_multi_file(land_clim_files, "Tair_f_inst", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qle  = read_var_multi_file(land_1201_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qle  = read_var_multi_file(land_clim_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_SM   = read_var_multi_file(land_1201_files, "SoilMoist_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_SM   = read_var_multi_file(land_clim_files, "SoilMoist_inst", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qh   = read_var_multi_file(land_1201_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qh   = read_var_multi_file(land_clim_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qair = read_var_multi_file(land_1201_files, "Qair_f_inst", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qair = read_var_multi_file(land_clim_files, "Qair_f_inst", loc_lat, loc_lon, "lat", "lon")

        print("time = ", time)
        dec1_Tmax       = spatial_var_max(time,Dec1_T,time_s,time_e)
        clim_Tmax       = spatial_var_max(time,Clim_T,time_s,time_e)

        dec1_Qle        = spatial_var(time,Dec1_Qle,time_s,time_e)
        clim_Qle        = spatial_var(time,Clim_Qle,time_s,time_e)

        dec1_SM_tmp     = spatial_var(time,Dec1_SM,time_s,time_e)
        clim_SM_tmp     = spatial_var(time,Clim_SM,time_s,time_e)

        dec1_SM         = (dec1_SM_tmp[0,:,:]*0.005+
                           dec1_SM_tmp[1,:,:]*0.075+
                           dec1_SM_tmp[2,:,:]*0.154+
                           dec1_SM_tmp[3,:,:]*0.409+
                           dec1_SM_tmp[4,:,:]*1.085+
                           dec1_SM_tmp[5,:,:]*0.272
                           )/2.

        clim_SM         = (clim_SM_tmp[0,:,:]*0.005+
                           clim_SM_tmp[1,:,:]*0.075+
                           clim_SM_tmp[2,:,:]*0.154+
                           clim_SM_tmp[3,:,:]*0.409+
                           clim_SM_tmp[4,:,:]*1.085+
                           clim_SM_tmp[5,:,:]*0.272
                           )/2.

        dec1_Qh         = spatial_var(time,Dec1_Qh,time_s,time_e)
        clim_Qh         = spatial_var(time,Clim_Qh,time_s,time_e)

        dec1_Qair       = spatial_var(time,Dec1_Qair,time_s,time_e)
        clim_Qair       = spatial_var(time,Clim_Qair,time_s,time_e)


        # "Tmax":
        tmax_diff       = dec1_Tmax-clim_Tmax

        # "Qle":
        qle_diff        = dec1_Qle-clim_Qle

        # "FWsoil":
        sm_diff         = dec1_SM-clim_SM

        # "Qh":
        qh_diff         = dec1_Qh-clim_Qh

        # "Qair":
        q_diff          = dec1_Qair-clim_Qair


        # ==================== Start to plot ====================
        for j in np.arange(5):

            ax[i,j].coastlines(resolution="50m",linewidth=1)
            ax[i,j].set_extent([135,155,-40,-23.6])
            ax[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            if j == 0:
                gl.ylabels_left   = True
            else:
                gl.ylabels_left   = False
            if i == 3:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            # set y label
            if j == 0:
                ax[i,j].text(-0.255, loc_y[i], label_y[i], va='bottom', ha='center',
                              rotation='vertical', rotation_mode='anchor',
                              transform=ax[i,j].transAxes)
                if i == 0:
                    ax[i,j].text(-0.36, 0, "Summer", va='bottom', ha='center',
                                  rotation='vertical', rotation_mode='anchor',fontsize=14,
                                  transform=ax[i,j].transAxes)
                    ax[i,j].text(-0.36, -2, "Heatwave Period", va='bottom', ha='center',
                                  rotation='vertical', rotation_mode='anchor',fontsize=14,
                                  transform=ax[i,j].transAxes)

        # left - SM
        clevs1  = [-0.4,-0.3,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.3,0.4] # -1,-0.9,-0.8,-0.7,-0.6, ,0.6,0.7,0.8,0.9,1
        plot1   = ax[i,0].contourf(lon, lat, sm_diff, levels=clevs1, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        cbar    = ax[i,0].text(0.02, 0.15, texts[cnt], transform=ax[i,0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray") # lightgray

        # middle left - Qle
        clevs2   = [ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100] #-140, -120, , 120, 140
        plot2    = ax[i,1].contourf(lon, lat, qle_diff, levels=clevs2, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        ax[i,1].text(0.02, 0.15, texts[cnt+1], transform=ax[i,1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle - Qh
        clevs3  = [ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100]
        # [ -140, -120, -100, -80, -60, -40, -20, -10]
        plot3   = ax[i,2].contourf(lon, lat, qh_diff, levels=clevs3, transform=ccrs.PlateCarree(),cmap=cmap_1,extend='both') #
        ax[i,2].text(0.02, 0.15, texts[cnt+2], transform=ax[i,2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle right - Tmax
        tmax_diff= np.where(np.isnan(qle_diff), np.nan, tmax_diff)
        clevs4   = [-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.25, 0.25, 0.5, 1., 1.5, 2., 2.5, 3.]
        plot4    = ax[i,3].contourf(lon, lat, tmax_diff, levels=clevs4, transform=ccrs.PlateCarree(),cmap=cmap_1,extend='both') #
        ax[i,3].text(0.02, 0.15, texts[cnt+3], transform=ax[i,3].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,3].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right - Qair
        clevs5  = [-2,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2.]
        plot5   = ax[i,4].contourf(lon, lat, q_diff*1000., levels=clevs5, transform=ccrs.PlateCarree(),cmap=cmap_2,extend='both') #
        ax[i,4].text(0.02, 0.15, texts[cnt+4], transform=ax[i,4].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[i,4].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")
        clevs   = None

        if shape_path != None:
        # Requires the pyshp package
            sf = shp.Reader(shape_path)

            for shape in sf.shapeRecords():
                x = [n[0] for n in shape.shape.points[:]]
                y = [n[1] for n in shape.shape.points[:]]
                ax[i,0].plot(x,y,c="black")
                ax[i,1].plot(x,y,c="black")
                ax[i,2].plot(x,y,c="black")
                ax[i,3].plot(x,y,c="black")
                ax[i,4].plot(x,y,c="black")

        # set top x label
        if i == 0:
            ax[i,0].set_title(label_x[0])#, fontsize=12)
            ax[i,1].set_title(label_x[1])#, fontsize=12)
            ax[i,2].set_title(label_x[2])#, fontsize=12)
            ax[i,3].set_title(label_x[3])#, fontsize=12)
            ax[i,4].set_title(label_x[4])#, fontsize=12)


        # set bottom colorbar
        if i == 3:

            # left - SM
            cbar = plt.colorbar(plot1, ax=ax[:,0], ticklocation="right", pad=0.02, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "m$\mathregular{^3}$ m$\mathregular{^{-3}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5])
            # cbar.ax.set_yticklabels(["-0.5","-0.4","-0.3","-0.2","-0.1","0.1","0.2","0.3","0.4","0.5"])

            # middle left - Qle
            cbar = plt.colorbar(plot2, ax=ax[:,1], ticklocation="right", pad=0.02, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100])
            # cbar.ax.set_yticklabels(["-100","-80","-60","-40","-20","-5","5","20","40","60","80","100"])

            # middle - Qh
            cbar = plt.colorbar(plot3, ax=ax[:,2], ticklocation="right", pad=0.02, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "W m$\mathregular{^{-2}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([ -100, -80, -60, -40, -20, -5, 5, 20, 40, 60, 80, 100])
            # cbar.ax.set_yticklabels(["-100","-80","-60","-40","-20","-5","5","20","40","60","80","100"])

            # middle right - Tmax
            cbar = plt.colorbar(plot4, ax=ax[:,3], ticklocation="right", pad=0.02, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "$\mathregular{^o}$C"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-3.,-2.5,-2.,-1.5,-1.,-0.5, -0.25, 0.25, 0.5, 1., 1.5, 2., 2.5, 3.])
            # cbar.ax.set_yticklabels(["-3","-2.5","-2","-1.5","-1","-0.5","-0.25","0.25","0.5","1.","1.5","2.","2.5","3."])

            # right - q
            cbar = plt.colorbar(plot5, ax=ax[:,4], ticklocation="right", pad=0.02, orientation="horizontal",
                                aspect=20, shrink=0.4)
            color_label= "g kg$\mathregular{^{-1}}$"
            cbar.set_label(color_label, loc='center',size=14)
            cbar.ax.tick_params(labelsize=8, rotation=90)
            # cbar.ax.set_yticks([-2,-1.75,-1.5,-1.25,-1.,-0.75,-0.5,-0.25,0.25,0.5,0.75,1,1.25,1.5,1.75,2.])
            # cbar.ax.set_yticklabels(["-2"," ","-1.5"," ","-1."," ","-0.5"," "," ","0.5"," ","1"," ","1.5"," ","2"])
        cnt = cnt + 5

    plt.savefig('./plots/spatial_map_SM_Qle_Qh_Tmax_q_2018_2020_summer_HW.pdf',dpi=300)

def plot_spatial_Rnet_LW_SW(land_1201_path, land_clim_path, atmo_1201_path, atmo_clim_path,
                            wrf_path, shape_path = None, loc_lat=None, loc_lon=None, seconds=None, message=None):

    # ======================= Make plots ========================
    # Three integers (nrows, ncols, index)

    fig, ax = plt.subplots(nrows=4, ncols=5, figsize=[16,14],sharex=True, sharey=True, squeeze=True,
                           subplot_kw={'projection': ccrs.PlateCarree()})
    plt.subplots_adjust(wspace=-0.44, hspace=0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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
    cmap       = plt.cm.seismic

    # ======================= Read WRF file =======================
    # use WRF output's lat & lon, since LIS output has default value
    wrf = Dataset(wrf_path,  mode='r')
    lon = wrf.variables['XLONG'][0,:,:]
    lat = wrf.variables['XLAT'][0,:,:]

    texts   = [ "(a)","(b)","(c)","(d)","(e)",
                "(f)","(g)","(h)","(i)","(j)",
                "(k)","(l)","(m)","(n)","(o)",
                "(p)","(q)","(r)","(s)","(t)"]

    label_x = ["ΔR$\mathregular{_{net}}$",
               "ΔLW$\mathregular{_{net}}$",
               "ΔLW$\mathregular{_{up}}$",
               "ΔLW$\mathregular{_{dn}}$",
               "ΔSW$\mathregular{_{net}}$"]

    label_y = ["2018-2019 Summer","2019-2020 Summer","2018-2019 HW","2019-2020 HW"]
    loc_y   = [0.63,0.55,0.47,0.38]
    cnt     = 0

    # ==================== Set up files ====================
    for i in np.arange(4):
        if i == 0:
            time_s = datetime(2018,12,1,0,0,0,0)
            time_e = datetime(2019,2,28,23,59,0,0)
            year_s = 2018
        elif i == 1:
            time_s = datetime(2019,1,14,0,0,0,0)
            time_e = datetime(2019,1,26,23,59,0,0)
            year_s = 2018
        elif i == 2:
            time_s = datetime(2019,12,1,0,0,0,0)
            time_e = datetime(2020,2,29,23,59,0,0)
            year_s = 2019
        elif i == 3:
            time_s = datetime(2019,12,16,0,0,0,0)
            time_e = datetime(2020,1,7,23,59,0,0)
            year_s = 2019

        land_1201_files = [ land_1201_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01.d01.nc",
                            land_1201_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02.d01.nc",]
        land_clim_files = [ land_clim_path+"LIS.CABLE."+str(year_s)+"12-"+str(year_s)+"12.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"01-"+str(year_s+1)+"01.d01.nc",
                            land_clim_path+"LIS.CABLE."+str(year_s+1)+"02-"+str(year_s+1)+"02.d01.nc",]

        atmo_1201_files = [ atmo_1201_path+"wrfout_d01_"+str(year_s)+"-12-01_01:00:00",
                            atmo_1201_path+"wrfout_d01_"+str(year_s+1)+"-01-01_01:00:00",
                            atmo_1201_path+"wrfout_d01_"+str(year_s+1)+"-02-01_01:00:00",]
        atmo_clim_files = [ atmo_clim_path+"wrfout_d01_"+str(year_s)+"-12-01_01:00:00",
                            atmo_clim_path+"wrfout_d01_"+str(year_s+1)+"-01-01_01:00:00",
                            atmo_clim_path+"wrfout_d01_"+str(year_s+1)+"-02-01_01:00:00",]

        time, Dec1_Qle   = read_var_multi_file(land_1201_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qle   = read_var_multi_file(land_clim_files, "Qle_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qh    = read_var_multi_file(land_1201_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qh    = read_var_multi_file(land_clim_files, "Qh_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Qg    = read_var_multi_file(land_1201_files, "Qg_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Qg    = read_var_multi_file(land_clim_files, "Qg_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Lwnet = read_var_multi_file(land_1201_files, "Lwnet_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Lwnet = read_var_multi_file(land_clim_files, "Lwnet_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_Swnet = read_var_multi_file(land_1201_files, "Swnet_tavg", loc_lat, loc_lon, "lat", "lon")
        time, Clim_Swnet = read_var_multi_file(land_clim_files, "Swnet_tavg", loc_lat, loc_lon, "lat", "lon")

        time, Dec1_LWUPB = read_wrf_surf_var_multi_files(atmo_1201_files, "LWUPB", loc_lat, loc_lon)
        time, Clim_LWUPB = read_wrf_surf_var_multi_files(atmo_clim_files, "LWUPB", loc_lat, loc_lon)

        time, Dec1_LWDNB = read_wrf_surf_var_multi_files(atmo_1201_files, "LWDNB", loc_lat, loc_lon)
        time, Clim_LWDNB = read_wrf_surf_var_multi_files(atmo_clim_files, "LWDNB", loc_lat, loc_lon)

        dec1_Qle         = spatial_var(time,Dec1_Qle,time_s,time_e)
        clim_Qle         = spatial_var(time,Clim_Qle,time_s,time_e)

        dec1_Qh          = spatial_var(time,Dec1_Qh,time_s,time_e)
        clim_Qh          = spatial_var(time,Clim_Qh,time_s,time_e)

        dec1_Qg          = spatial_var(time,Dec1_Qg,time_s,time_e)
        clim_Qg          = spatial_var(time,Clim_Qg,time_s,time_e)

        dec1_Lwnet       = spatial_var(time,Dec1_Lwnet,time_s,time_e)
        clim_Lwnet       = spatial_var(time,Clim_Lwnet,time_s,time_e)

        dec1_Swnet       = spatial_var(time,Dec1_Swnet,time_s,time_e)
        clim_Swnet       = spatial_var(time,Clim_Swnet,time_s,time_e)

        dec1_LWUPB       = spatial_var(time,Dec1_LWUPB,time_s,time_e)
        clim_LWUPB       = spatial_var(time,Clim_LWUPB,time_s,time_e)

        dec1_LWDNB       = spatial_var(time,Dec1_LWDNB,time_s,time_e)
        clim_LWDNB       = spatial_var(time,Clim_LWDNB,time_s,time_e)

        # "Rnet":
        rnet_diff        = (dec1_Qle+dec1_Qh+dec1_Qg)-(clim_Qle+clim_Qh+clim_Qg)

        # "Lwnet_tavg":
        lwnet_diff       = dec1_Lwnet-clim_Lwnet

        # "Swnet_tavg":
        swnet_diff       = dec1_Swnet-clim_Swnet

        # "LWup"
        lwup_diff       = dec1_LWUPB-clim_LWUPB

        # "LWdn"
        lwdn_diff       = dec1_LWDNB-clim_LWDNB

        # ==================== Start to plot ====================
        for j in np.arange(5):

            ax[i,j].coastlines(resolution="50m",linewidth=1)
            ax[i,j].set_extent([130,155,-44,-20])
            ax[i,j].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i,j].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([125,130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-45,-40,-35,-30,-25,-20,-15])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            if j == 0:
                gl.ylabels_left   = True
            else:
                gl.ylabels_left   = False
            if i == 2:
                gl.xlabels_bottom = True
            else:
                gl.xlabels_bottom = False

            # set y label
            if j == 0:
                ax[i,j].text(-0.23, loc_y[i], label_y[i], va='bottom', ha='center',
                                rotation='vertical', rotation_mode='anchor',
                                transform=ax[i,j].transAxes)


            clevs          = [-22,-18,-14,-10,-6,-2,2,6,10,14,18,22]
            clevs_label    = ["-22","-18","-14","-10","-6","-2","2","6","10","14","18","22"]
            # clevs    = [-25,-20,-15,-10,-5,5,10,15,20,25]

            # left - Rnet
            plot1    = ax[i,0].contourf(lon, lat, rnet_diff, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
            ax[i,0].text(0.02, 0.15, texts[cnt], transform=ax[i,0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
            ax[i,0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

            # left middle - LWnet
            plot2    = ax[i,1].contourf(lon, lat, lwnet_diff, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
            ax[i,1].text(0.02, 0.15, texts[cnt+1], transform=ax[i,1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
            ax[i,1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

            # middle - LWup
            lwup_diff = np.where(np.isnan(rnet_diff), np.nan, lwup_diff)
            plot3   = ax[i,2].contourf(lon, lat, lwup_diff, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
            ax[i,2].text(0.02, 0.15, texts[cnt+2], transform=ax[i,2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
            ax[i,2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

            # right middle - LWdn
            lwdn_diff = np.where(np.isnan(rnet_diff), np.nan, lwdn_diff)
            plot4   = ax[i,3].contourf(lon, lat, lwdn_diff, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
            ax[i,3].text(0.02, 0.15, texts[cnt+3], transform=ax[i,3].transAxes, fontsize=14, verticalalignment='top', bbox=props)
            ax[i,3].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

            # right - SWnet
            plot5   = ax[i,4].contourf(lon, lat, swnet_diff, levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both')
            ax[i,4].text(0.02, 0.15, texts[cnt+4], transform=ax[i,4].transAxes, fontsize=14, verticalalignment='top', bbox=props)
            ax[i,4].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

            if shape_path != None:
                # Requires the pyshp package
                sf = shp.Reader(shape_path)

                for shape in sf.shapeRecords():
                    x = [n[0] for n in shape.shape.points[:]]
                    y = [n[1] for n in shape.shape.points[:]]
                    ax[i,0].plot(x,y,c="black")
                    ax[i,1].plot(x,y,c="black")
                    ax[i,2].plot(x,y,c="black")
                    ax[i,3].plot(x,y,c="black")
                    ax[i,4].plot(x,y,c="black")

            if i == 3 and j == 4:
                cbar = plt.colorbar(plot1, ax=ax, ticklocation="right", pad=0.05, orientation="horizontal",
                                    aspect=40, shrink=0.8) # cax=cax,
                cbar.set_ticks(clevs)
                cbar.set_ticklabels(clevs_label) # cax=cax,
                color_label= "W m$\mathregular{^{-2}}$"
                cbar.set_label(color_label, loc='center',size=14)
                cbar.ax.tick_params(labelsize=12)

            # set top x label
            if i == 0:
                ax[i,0].set_title(label_x[0])#,labelpad=-0.1)#, fontsize=12)
                ax[i,1].set_title(label_x[1])#,labelpad=-0.1)#, fontsize=12)
                ax[i,2].set_title(label_x[2])#,labelpad=-0.1)#, fontsize=12)
                ax[i,3].set_title(label_x[3])#,labelpad=-0.1)#, fontsize=12)
                ax[i,4].set_title(label_x[4])#,labelpad=-0.1)#, fontsize=12)

            cnt = cnt + 5
    plt.savefig('./plots/spatial_map_Rnet_LW_SW_2018_2020_summer_HW.png',dpi=300)

if __name__ == "__main__":

    # ======================================
    # Decks for plot_spatial_tmax_qle_fwsoil
    # ======================================

    shape_path = "/g/data/w97/ad9701/drought_2017to2020/drought_focusArea/smooth_polygon_drought_focusArea.shp"

    region = "SE Aus" #"SE Aus" #"CORDEX" #"SE Aus"

    if region == "Aus":
        loc_lat    = [-44,-10]
        loc_lon    = [112,154]
    elif region == "SE Aus":
        loc_lat    = [-40,-23.6]
        loc_lon    = [134,155]
    elif region == "CORDEX":
        loc_lat    = [-52.36,3.87]
        loc_lon    = [89.25,180]

    case_name  = "bl_pbl2_mp4_sf_sfclay2" #"bl_pbl5_mp6_sf_sfclay1" #
    case_20181201  = "drght_2017_2019_"+case_name+"_20181201"
    case_20191201  = "drght_2017_2019_"+case_name+"_20191201"
    case_clim_2018 = "drght_2017_2019_"+case_name+"_climatology_2018"
    case_clim_2019 = "drght_2017_2019_"+case_name+"_climatology_2019"

    wrf_path       = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_clim_2019+"/WRF_output/wrfout_d01_2019-12-01_01:00:00"
    land_20181201_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_20181201+"/LIS_output/"
    land_20191201_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_20191201+"/LIS_output/"
    land_clim_2018_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_clim_2018+"/LIS_output/"
    land_clim_2019_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_clim_2019+"/LIS_output/"
    atmo_20181201_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_20181201+"/WRF_output/"
    atmo_20191201_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_20191201+"/WRF_output/"
    atmo_clim_2018_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_clim_2018+"/WRF_output/"
    atmo_clim_2019_path = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_HW/"+case_clim_2019+"/WRF_output/"
    # plot_spatial_Rnet_LW_SW(land_1201_path, land_clim_path, atmo_1201_path, atmo_clim_path, wrf_path, shape_path)
    plot_spatial_fwsoil_qle_qh_tmax_qrate(land_20181201_path,land_20191201_path, land_clim_2018_path, land_clim_2019_path, wrf_path, shape_path, loc_lat=loc_lat, loc_lon=loc_lon)

    # plot_spatial_sm_qle_qh_tmax_qair(land_20181201_path,land_20191201_path, land_clim_2018_path, land_clim_2019_path,wrf_path, shape_path, loc_lat=loc_lat, loc_lon=loc_lon)
