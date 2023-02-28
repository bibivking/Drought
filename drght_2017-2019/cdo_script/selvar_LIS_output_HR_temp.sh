#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q hugemem
#PBS -l walltime=10:00:00
#PBS -l mem=1470GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module load cdo
#var_name=""
# cdo selvar,lat,lon,time,Swnet_tavg,Lwnet_tavg,Qle_tavg,Qh_tavg,Qg_tavg /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904.nc /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_energy.nc

#cdo selvar,lat,lon,time,Snowf_tavg,Rainf_tavg,Evap_tavg,ECanop_tavg,TVeg_tavg,FWsoil_tavg,ESoil_tavg,WaterTableD_tavg,GWwb_tavg,Qs_tavg,Qsb_tavg,SoilMoist_inst /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904.nc /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_water.nc

cdo selvar,Tair_f_inst,VegT_tavg,AvgSurfT_tavg,SoilTemp_inst /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904.nc /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_temp.nc

#cdo selvar,Psurf_f_inst,Qair_f_inst,Tair_f_inst,Wind_f_inst,LWdown_f_inst,SWdown_f_inst /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904.nc /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_met.nc

# cdo selvar,Albedo_inst,LAI_inst,Landcover_inst /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904.nc /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_ALB_LAI.nc
