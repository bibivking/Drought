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


'''
plot a simple spatial map
'''

# var_names  = ["SAND_VEC","CLAY_VEC","SILT_VEC","OC_VEC","BULK_DEN_VEC"]
# var_names  = ["SAND","CLAY","SILT","OC","BULK_DEN"]
# var_names  = ["sand_vec","clay_vec","silt_vec","org_vec","rhosoil_vec"]
var_names  = ["CLAY"]
for var_name in var_names:
    for layer in np.arange(1,7,1):
        message   =  "Openlandmap_soilcomposition_CORDEX_180E_depth_varying"#_lyr="+str(layer)
        path      = "/g/data/w97/mm3972/scripts/wrf_scripts/make_LIS_landinfo/nc_file/"
        file_path = path +"Openlandmap_soilcomposition_CORDEX_180E_depth_varying_lyr"+str(layer)+".nc"

        # message   =  "Openlandmap_ELEV_DLCM_lyr="+str(layer)
        # file_path = "/g/data/w97/mm3972/model/cable/src/CABLE-AUX/offline/mmy_gridinfo_AU/gridinfo_AWAP_OpenLandMap_ELEV_DLCM_fix.nc"

        file   = Dataset(file_path, mode='r')
        var    = file.variables[var_name][:]#[layer-1,:,:]

        # ================== Start Plotting =================
        fig, ax = plt.subplots()
        clevs = np.arange(0.1,1.1,0.1)
        ax = plt.contourf(var,clevs)
        cb = plt.colorbar(ax=ax, orientation="vertical", pad=0.02, aspect=16, shrink=0.8)
        message = message + "_" + var_name
        plt.savefig('./plots/spatial_map_'+message+'.png',dpi=300)
