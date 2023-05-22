#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q hugemem
#PBS -l walltime=20:00:00
#PBS -l mem=1470GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module load cdo
cdo mergetime /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/WRF_output/wrfout_20* /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/drght_2017_2019_bl_pbl2_mp4_sf_sfclay2/WRF_output/wrfout_201701-201702_CLDFRA
