#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normal
#PBS -l walltime=15:00:00
#PBS -l mem=190GB
#PBS -l ncpus=8
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module use /g/data/hh5/public/modules 
module load conda/analysis3-22.04
nccompress -r -o -np 8 /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_energy.nc
nccompress -r -o -np 8 /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_temp.nc
nccompress -r -o -np 8 /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_met.nc
nccompress -r -o -np 8 /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB_HR/LIS_output/LIS.CABLE.201701-201904_ALB_LAI.nc
