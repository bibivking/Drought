#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q express
#PBS -l walltime=0:10:00
#PBS -l mem=20GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module use /g/data/hh5/public/modules 
module load conda/analysis3-unstable
cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019

python time_series_offline_spinup.py
