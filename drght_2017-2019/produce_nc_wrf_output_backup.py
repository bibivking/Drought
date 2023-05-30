#!/usr/bin/python

__author__ = "Mengyuan Mu"
__email__  = "mu.mengyuan815@gmail.com"


import os
import sys
import glob
import argparse
import netCDF4 as nc
from datetime import datetime, timedelta
import numpy as np
from wrf import (getvar, interplevel, get_cartopy, cartopy_xlim,
                 cartopy_ylim, to_np, latlon_coords, ALL_TIMES)
from common_utils import *

def cmd_line_parser():
    '''
    Get var_name and var_unit from command-line arguments
    '''

    p      = argparse.ArgumentParser() # define and parse command-line arguments, options, and subcommands
    p.add_argument("-n", dest="var_name", default='T2', help="Name of the variable")
    p.add_argument("-c", dest="case_name", default='drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2', help="Name of the experiment")
    # p.add_argument("-u", dest="var_unit", default='???', help="Units of the variable")
    args   = p.parse_args()

    return (args.var_name, args.case_name)


def read_wrf_var(file_paths, var_name, var_unit=None):

    # Open the NetCDF file
    encoding = 'utf-8' # Times in WRF output is btype, convert to string

    time     = []
    pre_time = datetime(2000, 1, 1, 0, 0, 1)
    file_num = len(file_paths)

    for num, file_path in enumerate(file_paths):

        ncfile      = nc.Dataset(file_path)

        # Read lat and lon
        if num == 0:
            lat  = ncfile.variables['XLAT'][0,:,:]
            lon  = ncfile.variables['XLONG'][0,:,:]

        # set the next file name
        if num+1 < file_num:
            ncfile_next    = nc.Dataset(file_paths[num+1])
            next_file_time = datetime.strptime(str(ncfile_next.variables['Times'][0,:], encoding),'%Y-%m-%d_%H:%M:%S')
            ncfile_next.close()

        # Read time
        ntime    = len(ncfile.variables['Times'][:,0])
        for i in np.arange(ntime):
            cur_time = datetime.strptime(str(ncfile.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
            if cur_time > pre_time and cur_time < next_file_time:
                time.append(cur_time)
            elif cur_time <= pre_time and cur_time < next_file_time:
                # Raise an exception and end the script
                raise Exception("Error: time doesn't monotonically increase in %s i=%s" % (file_path, str(i)))
            elif cur_time >= next_file_time:
                # this time step shouldn't be excluded
                time_step_stop = i
                # jump out of the loop
                break
            else:
                raise Exception("An error occurred, please check")
            time_step_stop = None
            pre_time       = cur_time

        # Read values and attributes
        Var_tmp = getvar(ncfile, var_name, timeidx=ALL_TIMES)
        var_unit = Var_tmp.units
        var_desc = Var_tmp.description
        var_coor = Var_tmp.coordinates

        if time_step_stop is not None:
            Var_tmp = Var_tmp[:time_step_stop,:,:]
        if num ==0:
            var     = Var_tmp[:,:,:]
        else:
            var     = np.append(var,Var_tmp[:,:,:],axis=0)

    print("var",var)
    print("type(time) ",type(time),"np.shape(time) ",np.shape(time))
    print("type(lat) ", type(lat), "np.shape(lat) ", np.shape(lat))
    print("type(lon) ", type(lon), "np.shape(lon) ", np.shape(lon))
    print("type(var) ", type(var), "np.shape(var) ", np.shape(var))

    return time, lat, lon, var, var_unit, var_desc, var_coor

# def read_wrf_var_4D(file_paths, var_name, var_unit=None):
#
#     # Open the NetCDF file
#     encoding = 'utf-8' # Times in WRF output is btype, convert to string
#
#     time_tmp = []
#     Var_tmp  = []
#     pre_time = datetime(2000, 1, 1, 0, 0, 1)
#
#     for num, file_path in enumerate(file_paths):
#
#         ncfile   = nc.Dataset(file_path)
#
#         # Read lat and lon
#         if num == 0:
#             lat  = ncfile.variables['XLAT'][0,:,:]
#             lon  = ncfile.variables['XLONG'][0,:,:]
#
#         # Read time
#         ntime    = len(ncfile.variables['Times'][:,0])
#         for i in np.arange(ntime):
#             cur_time = datetime.strptime(str(ncfile.variables['Times'][i,:], encoding),'%Y-%m-%d_%H:%M:%S')
#             if cur_time.year > pre_time.year:
#                 time_tmp.append(temp)
#             else:
#                 # Raise an exception and end the script
#                 raise Exception("An error occurred in time: cur_time.year <= pre_time.year in %s i=%s" % (file_path, str(i)))
#
#             pre_time = cur_time
#
#         # Read values
#         if num == 0:
#             # Extract the pressure and variable
#             if p_hgt == "p":
#                 hgt = getvar(ncfile, "pressure",timeidx=ALL_TIMES)
#             elif p_hgt == "hgt":
#                 hgt = getvar(ncfile, "height_agl",timeidx=ALL_TIMES)
#             print(np.shape(hgt))
#
#             Var     = getvar(ncfile, var_name, units=var_unit, timeidx=ALL_TIMES)
#             print("np.shape(Var)",np.shape(Var))
#             raise Exception("stop for checking")
#         else:
#             Var_tmp = getvar(ncfile, var_name, units=var_unit, timeidx=ALL_TIMES)
#             Var     = np.append(Var,Var_tmp,axis=0)
#             print("np.shape(Var)",np.shape(Var))
#             raise Exception("stop for checking")
#
#     # var     = interplevel(Var, p, height=[????])
#
#     print("type(time) ",type(time),"np.shape(time) ",np.shape(time))
#     print("type(lat) ", type(lat), "np.shape(lat) ", np.shape(lat))
#     print("type(lon) ", type(lon), "np.shape(lon) ", np.shape(lon))
#     print("type(var) ", type(var), "np.shape(var) ", np.shape(var))
#
#     return time, lat, lon, var, hgt # ???

def make_wrf_nc(file_paths, file_out_path, var_name):

    Time, Lat, Lon, Var, var_unit, var_desc, var_coor= read_wrf_var(file_paths, var_name)

    nLat = np.shape(Lat)[0]
    nLon = np.shape(Lat)[1]
    print("nLat=",nLat,"nLon=",nLon)

    # create file and write global attributes
    f = nc.Dataset(file_out_path, 'w', format='NETCDF4')

    ### Create nc file ###
    f.history           = "Created by: %s" % (os.path.basename(__file__))
    f.creation_date     = "%s" % (datetime.now())
    f.description       = 'wrf output content '+var_name+', created by MU Mengyuan'

    # set dimensions
    f.createDimension('time', None)
    f.createDimension('north_south', nLat)
    f.createDimension('east_west', nLon)
    f.Conventions       = "CF-1.0"

    # create variables
    time                = f.createVariable('time', 'f4', ('time'))
    time.standard_name  = "time"
    time.units          = "seconds since 2017-01-01 00:00:00"
    time.calendar       = "proleptic_gregorian"
    time.axis           = "T"

    latitude            = f.createVariable('lat', 'f4', ('north_south', 'east_west'))
    latitude.long_name  = "latitude"
    latitude.units      = "degree_north"
    latitude._CoordinateAxisType = "Lat"

    longitude           = f.createVariable('lon', 'f4', ('north_south', 'east_west'))
    longitude.long_name = "longitude"
    longitude.units     = "degree_east"
    longitude._CoordinateAxisType = "Lon"

    var                 = f.createVariable(var_name, 'f4', ('time', 'north_south', 'east_west'))
    var.units           = var_unit
    var.standard_name   = var_name
    var.description     = var_desc
    var.coordinates     = var_coor

    latitude[:,:]       = Lat
    longitude[:,:]      = Lon
    var[:,:,:]          = Var
    for i,t in enumerate(Time):
        time[i]         = (t - datetime(2017, 1, 1, 0, 0, 0)).total_seconds()

    f.close()

if __name__ == "__main__":

    # ======================= Option =======================
    # Get var_name and var_unit from command-line arguments
    var_name, case_name = cmd_line_parser()

    # ======================= Setting path =======================
    atmo_path     = '/g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/'+case_name +'/WRF_output/'
    file_out_path = atmo_path+var_name+'/wrfout_201701-201912.nc'

    # glob.glob: retrieve files and directories that match a specified pattern.
    # sorted: rank the file paths
    file_paths    = sorted(glob.glob(atmo_path + "wrfout_d01_2017-01*00:00"))
    print(file_paths)

    # ==================== Generate nc file ======================
    # make one variable wrf output netcdf file
    make_wrf_nc(file_paths, file_out_path, var_name)
