#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q expressbw
#PBS -l walltime=1:20:00
#PBS -l mem=190GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/ua8

module use /g/data/hh5/public/modules
module load conda/analysis3-22.01
cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019_VPD
python spatial_map_compare_obs.py
