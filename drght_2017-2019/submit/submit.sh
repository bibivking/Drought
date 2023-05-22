#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalbw
#PBS -l walltime=3:30:00
#PBS -l mem=256GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-22.01
cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019
#python spatial_map_time_series_LAI_ALB.py
#python spatial_map_time_series_LAI_ALB.py
#python spatial_map_single_plot.py
python time_series_wrf_cable.py
