#!/usr/bin/env python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"

import sys
import cartopy
import numpy as np
from netCDF4 import Dataset,num2date
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature


# '''
# plot a simple spatial map
# '''

# # var_names  = ["SAND_VEC","CLAY_VEC","SILT_VEC","OC_VEC","BULK_DEN_VEC"]
# # var_names  = ["SAND","CLAY","SILT","OC","BULK_DEN"]
# # var_names  = ["sand_vec","clay_vec","silt_vec","org_vec","rhosoil_vec"]
# var_names  = ["CLAY"]
# for var_name in var_names:
#     for layer in np.arange(1,7,1):
#         message   =  "Openlandmap_soilcomposition_CORDEX_180E_depth_varying"#_lyr="+str(layer)
#         path      = "/g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/"
#         file_path = path +"Openlandmap_soilcomposition_CORDEX_180E_depth_varying_lyr"+str(layer)+".nc"

#         # message   =  "Openlandmap_ELEV_DLCM_lyr="+str(layer)
#         # file_path = "/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc"

#         file   = Dataset(file_path, mode='r')
#         var    = file.variables[var_name][:]#[layer-1,:,:]

#         # ================== Start Plotting =================
#         fig, ax = plt.subplots()
#         clevs = np.arange(0.1,1.1,0.1)
#         ax = plt.contourf(var,clevs)
#         cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
#         message = message + "_" + var_name
#         plt.savefig('./plots/spatial_map_'+message+'.png',dpi=300)


'''
plot a simple spatial difference map
'''

file_name1 = "/g/data/w97/mm3972/data/AWAP_detrend/AWAP_detrend_10km/Tair/AWAP.Tair.3hr.2000.nc"
# file_name1 = "/g/data/w97/mm3972/model/cable/runs/VPD_drought/old_runs/100th_check/outputs/cable_out_2017.nc"
file_name2 = "/g/data/w97/mm3972/data/AWAP_3h_v1_10km/Tair/AWAP.Tair.3hr.2000.nc"
# file_name2 =  "/g/data/w97/mm3972/model/cable/runs/runs_4_coupled/gw_after_sp30yrx3/outputs/cable_out_2000-2019.nc"
var_name = "Tair"
message  = "Tair"
scale    = 1
clevs    = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,-0.5,-0.1,0.1,0.5,1,2,3,4,5,6,7,8,9,10]
# [-10,-1,-0.01,0.01,1,10]

f1         = Dataset(file_name1, mode='r')
f2         = Dataset(file_name2, mode='r')

var1       = np.nanmean(f1.variables[var_name][:366,:,:],axis=0)
var2       = np.nanmean(f2.variables[var_name][:366,:,:],axis=0)
var_diff   = (var2-var1)*scale

print(var_diff)

fig, ax = plt.subplots()

plot    = ax.contourf(var_diff,clevs,cmap="BrBG")#"viridis_r") #"BrBG") 
cb      = plt.colorbar(plot)#, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)

plt.show()
plt.savefig('./plots/spatial_map_'+message+'_'+var_name+'_diff.png',dpi=300)
