import sys
import cartopy
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from scipy.interpolate import griddata, interp1d
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature, OCEAN
from common_utils import *

def plot_LAI(LAI_obs_path, LAI_wrf_path, wrf_path, LIS_input, time_s, time_e):

    nlat = 439
    nlon = 529
    PFT  = True
    # read LIS_input LAI
    lis_file             = Dataset(LIS_input, mode='r')
    LAI_lis_tmp          = lis_file.variables['LAI'][:]
    LAI_lis_t            = np.zeros((12*3+2,nlat,nlon))
    LAI_lis_t[0,:,:]     = LAI_lis_tmp[11,:,:]
    LAI_lis_t[1:13,:,:]  = LAI_lis_tmp
    LAI_lis_t[13:25,:,:] = LAI_lis_tmp
    LAI_lis_t[25:37,:,:] = LAI_lis_tmp
    LAI_lis_t[37,:,:]    = LAI_lis_tmp[0,:,:]

    # middle of each month during 2016 Dec -2020 Jan
    time_lis_in = [ 6193,
                    6224,6254,6283,6314,6344,6375,6405,6436,6467,6497,
                    6528,6558,6589,6619,6648,6679,6709,6740,6770,6801,
                    6832,6862,6893,6923,6954,6984,7013,7044,7074,7105,
                    7135,7166,7197,7227,7258,7288,
                    7319]

    time_init    = datetime(2000,1,1,0,0,0)
    deltatime_s  = time_s - time_init
    deltatime_e  = time_e - time_init
    time_out     = np.arange(deltatime_s.days,deltatime_e.days+1,1)
    f            = interp1d(time_lis_in, LAI_lis_t, kind='linear',axis=0)
    LAI_lis      = f(time_out)

    # # read WRF LAI
    # LAI_wrf = np.zeros((365*3,nlat,nlon))
    # dom     = [31,28,31,30,31,30,31,31,30,31,30,31]
    # d_s     = 0
    # for yr in np.arange(2017,2020):
    #     print("year="+str(yr))
    #     for mth in np.arange(1,13):
    #         d_e = d_s + dom[mth-1]
    #         print("month="+str(mth))
    #         if mth < 10:
    #             LAI_wrf_file = LAI_wrf_path + "LIS.CABLE."+str(yr)+"0"+str(mth)+"-"+str(yr)+"0"+str(mth)+".d01.nc"
    #         else:
    #             LAI_wrf_file = LAI_wrf_path + "LIS.CABLE."+str(yr)+str(mth)+"-"+str(yr)+str(mth)+".d01.nc"
    #
    #         time_wrf, LAI_wrf_tmp = read_var(LAI_wrf_file, "LAI_inst")
    #
    #     print("np.shape(LAI_wrf)")
    #
    #     LAI_wrf[d_s:d_e,:,:]  = time_clip_to_day(time_wrf[:-10], LAI_wrf_tmp[:-10,:,:], time_s, time_e)
    #     print(np.shape(LAI_wrf))
    #     # # ================= testing plot ===================
    #     # for i in np.arange(dom[mth-1]):
    #     #     fig, axs = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True,
    #     #                             subplot_kw={'projection': ccrs.PlateCarree()})
    #     #     plot     = axs.contourf( lon, lat, LAI_wrf[i,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic,extend='both') # levels=clevs,
    #     #     cbar     = plt.colorbar(plot, ax=axs, ticklocation="right", pad=0.05, orientation="horizontal",aspect=40, shrink=0.8) # cax=cax,
    #     #     plt.savefig('./plots/spatial_map_LAI_WRF_i='+str(i)+'.png',dpi=300)
    #     #     fig = None
    #     #     axs = None
    #     #     plot = None

    # read obs LAI
    time_obs,LAI_obs_tmp = read_var(LAI_obs_path, "LAI")
    LAI_obs              = time_clip_to_day(time_obs, LAI_obs_tmp, time_s, time_e)

    LAI_obs              = np.where(LAI_obs ==0, np.nan,LAI_obs)
    # LAI_wrf              = np.where(LAI_wrf ==-9999., np.nan,LAI_wrf)
    LAI_lis              = np.where(LAI_lis ==-9999., np.nan,LAI_lis)
    LAI_obs_mean         = np.nanmean(LAI_obs,axis=0)
    # LAI_wrf_mean         = np.nanmean(LAI_wrf,axis=0)
    LAI_lis_mean         = np.nanmean(LAI_lis,axis=0)

    LAI_obs_time_series  = np.nanmean(LAI_obs,axis=(1,2))
    # LAI_wrf_time_series  = np.nanmean(LAI_wrf,axis=(1,2))
    LAI_lis_time_series  = np.nanmean(LAI_lis,axis=(1,2))

    # read WRF PFT
    PFT_wrf_file = LAI_wrf_path + "LIS.CABLE.201702-201702.d01.nc"
    pft_wrf      = Dataset(PFT_wrf_file, mode='r')
    PFT_tmp      = pft_wrf.variables['Landcover_inst'][0,:,:]
    PFT          = [PFT_tmp] * (365*3)
    print("np.shape(PFT)")
    print(np.shape(PFT))
    LAI_obs_pft  = np.zeros((17,365*3,nlat,nlon))
    LAI_lis_pft  = np.zeros((17,365*3,nlat,nlon))

    for i in np.arange(1,18,1):
        LAI_obs_pft[i-1,:,:,:] = np.where(PFT == i, LAI_obs, np.nan)
        LAI_lis_pft[i-1,:,:,:] = np.where(PFT == i, LAI_lis, np.nan)
    LAI_obs_mean_pft           = np.nanmean(LAI_obs_pft,axis=1)
    LAI_lis_mean_pft           = np.nanmean(LAI_lis_pft,axis=1)
    print("np.shape(LAI_obs_mean_pft)")
    print(np.shape(LAI_obs_mean_pft))
    LAI_obs_time_series_pft    = np.nanmean(LAI_obs_pft,axis=(2,3))
    LAI_lis_time_series_pft    = np.nanmean(LAI_lis_pft,axis=(2,3))
    print("np.shape(LAI_obs_time_series_pft)")
    print(np.shape(LAI_obs_time_series_pft))

    # read lat and lon outs
    wrf                  = Dataset(wrf_path,  mode='r')
    lon                  = wrf.variables['XLONG'][0,:,:]
    lat                  = wrf.variables['XLAT'][0,:,:]

    # # ================= testing plot ===================
    if 0:
        for i in np.arange(1,18,1):
            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True,
                                    subplot_kw={'projection': ccrs.PlateCarree()})
            plot     = axs.contourf( lon, lat, LAI_obs_mean_pft[i-1,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic,extend='both') # levels=clevs,
            cbar     = plt.colorbar(plot, ax=axs, ticklocation="right", pad=0.05, orientation="horizontal",aspect=40, shrink=0.8) # cax=cax,
            plt.savefig('./plots/spatial_map_OBS_PFT_i='+str(i)+'.png',dpi=300)
            fig = None
            axs = None
            plot = None

        for i in np.arange(1,18,1):
            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True,
                                    subplot_kw={'projection': ccrs.PlateCarree()})
            plot     = axs.contourf( lon, lat, LAI_lis_mean_pft[i-1,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic,extend='both') # levels=clevs,
            cbar     = plt.colorbar(plot, ax=axs, ticklocation="right", pad=0.05, orientation="horizontal",aspect=40, shrink=0.8) # cax=cax,
            plt.savefig('./plots/spatial_map_WRF_PFT_i='+str(i)+'.png',dpi=300)
            fig = None
            axs = None
            plot = None

    # =================== Plotting time series ===================
    if 0:
        ls_mark       = ['-','--','-.',':','--',',','o','v','^','<','>','1','2','3','4','s','p']
        cleaner_dates = ["2017","2018", "2019", "2020" ]
        xtickslocs    = [0,     365,      730,   1095  ]

        fig, ax       = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True)

        print("np.shape(LAI_obs_time_series)")
        print(np.shape(LAI_obs_time_series))
        cnt = 0
        for i in [2,5,6,9,14]:
            ax.plot(np.arange(len(LAI_obs_time_series[:])), LAI_obs_time_series_pft[i-1,:], c = 'blue', ls=ls_mark[cnt], label='obs'+str(i), alpha=0.5) #
            ax.plot(np.arange(len(LAI_obs_time_series[:])), LAI_lis_time_series_pft[i-1,:], c = 'red', ls=ls_mark[cnt], label='clim'+str(i), alpha=0.5) #
            cnt = cnt + 1
        ax.legend()
        fig.tight_layout()

        plt.savefig('./plots/time_series_LAI.png',dpi=300)


    # =================== Plotting spatial map ===================
    if 1:
        # for j in np.arange(17):
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[10,4],sharex=True, sharey=True, squeeze=True,
                                subplot_kw={'projection': ccrs.PlateCarree()})
        # plt.subplots_adjust(wspace=-0.44, hspace=0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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
        cmap    = plt.cm.BrBG
        texts   = [ "(a)","(b)","(c)","(d)","(e)",
                    "(f)","(g)","(h)","(i)","(j)",
                    "(k)","(l)","(m)","(n)","(o)",
                    "(p)","(q)","(r)","(s)","(t)"]

        label_x = ["LAI$\mathregular{_{obs}}$",
                    "LAI$\mathregular{_{clim}}$",
                    "ΔLAI$\mathregular{_{clim-obs}}$",]

        label_y = ["Annual LAI","Spring LAI","Summer LAI","Autumn LAI","Winter LAI"]
        loc_y   = [0.63,0.55,0.47,0.38]
        cnt     = 0


        for i in np.arange(3):

            ax[i].coastlines(resolution="50m",linewidth=1)
            ax[i].set_extent([135,155,-39,-23])
            ax[i].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25,-20])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            gl.xlabels_bottom = True
            gl.ylabels_left   = True


        clevs          = np.arange(0,6,0.2)
        clevs_diff     = [-1,-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8,1.]

        # left - LAI obs
        # LAI_lis_mean
        plot1    = ax[0].contourf(lon, lat, LAI_lis_mean[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both') #
        ax[0].text(0.02, 0.15, texts[0], transform=ax[0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle - LAI clim
        plot2    = ax[1].contourf(lon, lat, LAI_obs_mean[:,:], levels=clevs,transform=ccrs.PlateCarree(),cmap=cmap,extend='both') # levels=clevs,
        ax[1].text(0.02, 0.15, texts[1], transform=ax[1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right - LAI clim - obs
        plot3   = ax[2].contourf(lon, lat, LAI_lis_mean[:,:]-LAI_obs_mean[:,:],levels=clevs_diff, transform=ccrs.PlateCarree(),cmap=cmap,extend='both') #  levels=clevs_diff,
        ax[2].text(0.02, 0.15, texts[2], transform=ax[2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        cbar = plt.colorbar(plot1, ax=ax[0], ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.8) # cax=cax,
        cbar.ax.tick_params(labelsize=12)

        cbar = plt.colorbar(plot2, ax=ax[1], ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.8) # cax=cax,
        cbar.ax.tick_params(labelsize=12)

        cbar = plt.colorbar(plot3, ax=ax[2], ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.8)
        cbar.ax.tick_params(labelsize=12)

        # set top x label
        ax[0].set_title(label_x[0])#,labelpad=-0.1)#, fontsize=12)
        ax[1].set_title(label_x[1])#,labelpad=-0.1)#, fontsize=12)
        ax[2].set_title(label_x[2])#,labelpad=-0.1)#, fontsize=12)

        plt.savefig('./plots/spatial_map_LAI_2017_2019.png',dpi=300)
        fig = None
        ax  = None
        plot1 = None
        plot2 = None
        plot3 = None

def plot_ALBEDO(ALB_obs_path, ALB_wrf_path, wrf_path, LIS_input, time_s, time_e):
    
    nlat = 439
    nlon = 529
    PFT  = True
    
    time_init    = datetime(2000,1,1,0,0,0)
    deltatime_s  = time_s - time_init 
    deltatime_e  = time_e - time_init
    time_out     = np.arange(deltatime_s.days,deltatime_e.days+1,1)
    
    # =========== read albedo from lis_input.d01.nc ===========
    if 0:
        lis_file             = Dataset(LIS_input, mode='r')
        ALB_lis_tmp          = lis_file.variables['ALBEDO'][:]
        ALB_lis_t            = np.zeros((12*3+2,nlat,nlon))
        ALB_lis_t[0,:,:]     = ALB_lis_tmp[11,:,:]
        ALB_lis_t[1:13,:,:]  = ALB_lis_tmp
        ALB_lis_t[13:25,:,:] = ALB_lis_tmp
        ALB_lis_t[25:37,:,:] = ALB_lis_tmp
        ALB_lis_t[37,:,:]    = ALB_lis_tmp[0,:,:]

        # middle of each month during 2016 Dec -2020 Jan
        time_lis_in = [ 6193,
                        6224,6254,6283,6314,6344,6375,6405,6436,6467,6497,
                        6528,6558,6589,6619,6648,6679,6709,6740,6770,6801,
                        6832,6862,6893,6923,6954,6984,7013,7044,7074,7105,
                        7135,7166,7197,7227,7258,7288,
                        7319]

        f            = interp1d(time_lis_in, ALB_lis_t, kind='linear',axis=0)
        ALB_lis      = f(time_out)

    # =========== read albedo from wrflowinp_d01_2017-12 ===========
    if 1:
        # read WRF ALB
        ALB_lis = np.zeros((365*3,nlat,nlon)) # 20170101-20200101
        
        for yr in np.arange(2017,2020):
            print("year="+str(yr))
            if yr == 2020:
                mth_e = 7
            else:
                mth_e = 13
            for mth in np.arange(1,mth_e):
                # d_e = d_s + dom[mth-1]
                print("month="+str(mth))
                if mth < 10:
                    ALB_wrf_file = ALB_wrf_path + "LIS.CABLE."+str(yr)+"0"+str(mth)+"-"+str(yr)+"0"+str(mth)+".d01.nc"
                else:
                    ALB_wrf_file = ALB_wrf_path + "LIS.CABLE."+str(yr)+str(mth)+"-"+str(yr)+str(mth)+".d01.nc"

                time_wrf, ALB_wrf_tmp = read_var(ALB_wrf_file, "Albedo_inst")

                time_from_2017     = datetime(2017,1,1,0,0,0) - time_init
                time_from_2017_day = time_from_2017.days
                d_s                = time_wrf[0].days - time_from_2017_day
                if yr == 2019 and mth_e == 12:
                    d_e  = time_wrf[-1].days - time_from_2017_day
                else:
                    d_e  = time_wrf[-1].days - time_from_2017_day +1
                print("d_s=",d_s)
                print("d_e=",d_e)
                
                ALB_lis[d_s:d_e,:,:]  = time_clip_to_day(time_wrf[:], ALB_wrf_tmp[:,:,:], time_s, time_e)
                print(np.shape(ALB_lis))
                print("ALB_lis",ALB_lis)
        
    # =========== read obs ALB ===========
    time_obs,ALB_obs_tmp = read_var(ALB_obs_path, "shortwave")
    print(time_obs)
    
    # interpolate to daily 
    Time_days    = np.zeros(len(time_obs),dtype=int)
    
    for i in np.arange(len(time_obs)):
        Time_days[i] = time_obs[i].days
    print(Time_days)
        
    g            = interp1d(Time_days, ALB_obs_tmp, kind='linear',axis=0)
    ALB_obs      = g(time_out)
    
    # =========== read obs ABL ===========
    ALB_obs      = np.where(ALB_obs ==0, np.nan,ALB_obs)
    ALB_lis      = np.where(ALB_lis ==-9999., np.nan,ALB_lis)
    
    print("ALB_obs ", ALB_obs)
    print("ALB_lis ", ALB_lis)
    
    ALB_obs_mean         = np.nanmean(ALB_obs,axis=0)
    ALB_lis_mean         = np.nanmean(ALB_lis,axis=0)

    ALB_obs_time_series  = np.nanmean(ALB_obs,axis=(1,2))
    ALB_lis_time_series  = np.nanmean(ALB_lis,axis=(1,2))

    # read WRF PFT
    PFT_wrf_file = ALB_wrf_path + "LIS.CABLE.201702-201702.d01.nc"
    pft_wrf      = Dataset(PFT_wrf_file, mode='r')
    PFT_tmp      = pft_wrf.variables['Landcover_inst'][0,:,:]
    PFT          = [PFT_tmp] * (365*3)
    print("np.shape(PFT)")
    print(np.shape(PFT))
    ALB_obs_pft  = np.zeros((17,365*3,nlat,nlon))
    ALB_lis_pft  = np.zeros((17,365*3,nlat,nlon))

    for i in np.arange(1,18,1):
        ALB_obs_pft[i-1,:,:,:] = np.where(PFT == i, ALB_obs, np.nan)
        ALB_lis_pft[i-1,:,:,:] = np.where(PFT == i, ALB_lis, np.nan)
    ALB_obs_mean_pft           = np.nanmean(ALB_obs_pft,axis=1)
    ALB_lis_mean_pft           = np.nanmean(ALB_lis_pft,axis=1)
    print("np.shape(ALB_obs_mean_pft)")
    print(np.shape(ALB_obs_mean_pft))
    ALB_obs_time_series_pft    = np.nanmean(ALB_obs_pft,axis=(2,3))
    ALB_lis_time_series_pft    = np.nanmean(ALB_lis_pft,axis=(2,3))
    print("np.shape(ALB_obs_time_series_pft)")
    print(np.shape(ALB_obs_time_series_pft))

    # read lat and lon outs
    wrf                  = Dataset(wrf_path,  mode='r')
    lon                  = wrf.variables['XLONG'][0,:,:]
    lat                  = wrf.variables['XLAT'][0,:,:]

    # # ================= testing plot ===================
    if 0:
        for i in np.arange(1,18,1):
            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True,
                                    subplot_kw={'projection': ccrs.PlateCarree()})
            plot     = axs.contourf( lon, lat, ALB_obs_mean_pft[i-1,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic,extend='both') # levels=clevs,
            cbar     = plt.colorbar(plot, ax=axs, ticklocation="right", pad=0.05, orientation="horizontal",aspect=40, shrink=0.8) # cax=cax,
            plt.savefig('./plots/spatial_map_OBS_PFT_i='+str(i)+'.png',dpi=300)
            fig = None
            axs = None
            plot = None

        for i in np.arange(1,18,1):
            fig, axs = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True,
                                    subplot_kw={'projection': ccrs.PlateCarree()})
            plot     = axs.contourf( lon, lat, ALB_lis_mean_pft[i-1,:,:], transform=ccrs.PlateCarree(),cmap=plt.cm.seismic,extend='both') # levels=clevs,
            cbar     = plt.colorbar(plot, ax=axs, ticklocation="right", pad=0.05, orientation="horizontal",aspect=40, shrink=0.8) # cax=cax,
            plt.savefig('./plots/spatial_map_WRF_PFT_i='+str(i)+'.png',dpi=300)
            fig = None
            axs = None
            plot = None

    # =================== Plotting time series ===================
    if 1:
        ls_mark       = ['-','--','-.',':','--',',','o','v','^','<','>','1','2','3','4','s','p']
        cleaner_dates = ["2017","2018", "2019", "2020" ]
        xtickslocs    = [0,     365,      730,   1095  ]

        fig, ax       = plt.subplots(nrows=1, ncols=1, figsize=[4,4],sharex=True, sharey=True, squeeze=True)

        print("np.shape(ALB_obs_time_series)")
        print(np.shape(ALB_obs_time_series))
        cnt = 0
        for i in [2,5,6,9,14]:
            ax.plot(np.arange(len(ALB_obs_time_series[:])), ALB_obs_time_series_pft[i-1,:], c = 'blue', ls=ls_mark[cnt], label='obs'+str(i), alpha=0.5) #
            ax.plot(np.arange(len(ALB_obs_time_series[:])), ALB_lis_time_series_pft[i-1,:], c = 'red', ls=ls_mark[cnt], label='clim'+str(i), alpha=0.5) #
            cnt = cnt + 1
        ax.legend()
        fig.tight_layout()

        plt.savefig('./plots/time_series_ALB.png',dpi=300)


    # =================== Plotting spatial map ===================
    if 1:
        # for j in np.arange(17):
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[10,4],sharex=True, sharey=True, squeeze=True,
                                subplot_kw={'projection': ccrs.PlateCarree()})
        # plt.subplots_adjust(wspace=-0.44, hspace=0) # left=0.15,right=0.95,top=0.85,bottom=0.05,

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
        cmap    = plt.cm.BrBG
        texts   = [ "(a)","(b)","(c)","(d)","(e)",
                    "(f)","(g)","(h)","(i)","(j)",
                    "(k)","(l)","(m)","(n)","(o)",
                    "(p)","(q)","(r)","(s)","(t)"]

        label_x = ["ALB$\mathregular{_{obs}}$",
                    "ALB$\mathregular{_{clim}}$",
                    "ΔALB$\mathregular{_{obs-clim}}$",]

        label_y = ["Annual ALB","Spring ALB","Summer ALB","Autumn ALB","Winter ALB"]
        loc_y   = [0.63,0.55,0.47,0.38]
        cnt     = 0


        for i in np.arange(3):

            ax[i].coastlines(resolution="50m",linewidth=1)
            ax[i].set_extent([135,155,-39,-23])
            ax[i].add_feature(states, linewidth=.5, edgecolor="black")

            # Add gridlines
            gl = ax[i].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color=almost_black, linestyle='--')
            gl.xlabels_top  = False
            gl.ylabels_right= False
            gl.xlines       = False
            gl.ylines       = False
            gl.xlocator     = mticker.FixedLocator([130,135,140,145,150,155,160])
            gl.ylocator     = mticker.FixedLocator([-40,-35,-30,-25,-20])
            gl.xformatter   = LONGITUDE_FORMATTER
            gl.yformatter   = LATITUDE_FORMATTER
            gl.xlabel_style = {'size':12, 'color':almost_black}#,'rotation': 90}
            gl.ylabel_style = {'size':12, 'color':almost_black}

            gl.xlabels_bottom = True
            gl.ylabels_left   = True


        clevs          = np.arange(0,0.3,0.005)
        clevs_diff     = [-0.035,-0.03,-0.025,-0.02,-0.015,-0.01,-0.005,0.005,0.01,0.015,0.02,0.025,0.03,0.035]

        # left - ALB obs
        # ALB_lis_mean
        plot1    = ax[0].contourf(lon, lat, ALB_obs_mean[:,:], levels=clevs, transform=ccrs.PlateCarree(),cmap=cmap,extend='both') #
        ax[0].text(0.02, 0.15, texts[0], transform=ax[0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[0].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # middle - ALB clim
        plot2    = ax[1].contourf(lon, lat, ALB_lis_mean[:,:], levels=clevs,transform=ccrs.PlateCarree(),cmap=cmap,extend='both') # levels=clevs,
        ax[1].text(0.02, 0.15, texts[1], transform=ax[1].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[1].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        # right - ALB clim - obs
        plot3   = ax[2].contourf(lon, lat, ALB_obs_mean[:,:]-ALB_lis_mean[:,:],levels=clevs_diff, transform=ccrs.PlateCarree(),cmap=cmap,extend='both') #  levels=clevs_diff,
        ax[2].text(0.02, 0.15, texts[2], transform=ax[2].transAxes, fontsize=14, verticalalignment='top', bbox=props)
        ax[2].add_feature(OCEAN,edgecolor='none', facecolor="lightgray")

        cbar = plt.colorbar(plot1, ax=ax[0], ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.8) # cax=cax,
        cbar.ax.tick_params(labelsize=12)

        cbar = plt.colorbar(plot2, ax=ax[1], ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.8) # cax=cax,
        cbar.ax.tick_params(labelsize=12)

        cbar = plt.colorbar(plot3, ax=ax[2], ticklocation="right", pad=0.08, orientation="horizontal",
                            aspect=40, shrink=0.8)
        cbar.ax.tick_params(labelsize=12)

        # set top x label
        ax[0].set_title(label_x[0])#,labelpad=-0.1)#, fontsize=12)
        ax[1].set_title(label_x[1])#,labelpad=-0.1)#, fontsize=12)
        ax[2].set_title(label_x[2])#,labelpad=-0.1)#, fontsize=12)

        plt.savefig('./plots/spatial_map_ALB_2017_2019_obs_LAI.png',dpi=300)
        fig = None
        ax  = None
        plot1 = None
        plot2 = None
        plot3 = None

if __name__ == "__main__":
    
    # LAI
    if 0:
        LAI_obs_path  = "/g/data/w97/mm3972/data/MODIS/whittakerSmoothed_MCD15A3H_c61_LAI_for_WRF_daily_20170101_20200630.nc"
        LAI_wrf_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/LIS_output/"
        wrf_path      = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/WRF_output/wrfout_d01_2017-01-01_11:00:00"
        LIS_input     = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/bdy_data/lis_input.d01.nc"

        time_s        = datetime(2017,1,1,0,0,0,0)
        time_e        = datetime(2019,12,31,23,59,0,0)

        plot_LAI(LAI_obs_path, LAI_wrf_path, wrf_path, LIS_input, time_s, time_e)
    
    # ALBEDO
    if 1:
        ALB_obs_path  = "/g/data/w97/mm3972/data/MODIS/MODIS_Albedo/MCD43C3_bigWRFroi_smoothed-Albedo_BSA_shortwave_5000m_for_WRF_2016_2020.nc"
        ALB_wrf_path  = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI/LIS_output/"
        wrf_path      = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/WRF_output/wrfout_d01_2017-01-01_11:00:00"
        LIS_input     = "/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/bdy_data/lis_input.d01.nc"

        time_s        = datetime(2017,1,1,0,0,0,0)
        time_e        = datetime(2019,12,31,23,59,0,0)

        plot_ALBEDO(ALB_obs_path, ALB_wrf_path, wrf_path, LIS_input, time_s, time_e)
