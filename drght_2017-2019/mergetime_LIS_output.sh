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
cdo mergetime /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/LIS_output/LIS.CABLE.??????-??????.d01.nc /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/LIS_output/LIS.CABLE.201701-202006.nc
