#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q express
#PBS -l walltime=5:00:00
#PBS -l mem=190GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+scratch/w97+gdata/w97+gdata/ua8

module use /g/data/hh5/public/modules
module load conda/analysis3-22.01
cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019
python spatial_map_single_plot.py
#time_series_wrf_cable.py
#spatial_map_time_series_LAI.py
# python spatial_map_T_R_metrics.py
