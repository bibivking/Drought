#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q megamem
#PBS -l walltime=1:40:00
#PBS -l mem=1024GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+scratch/w97+gdata/w97+gdata/ua8

module use /g/data/hh5/public/modules
module load conda/analysis3-22.01
cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019/
python spatial_map_LAI_Albedo_Drought_relationship.py
#spatial_map_T_R_metrics.py
