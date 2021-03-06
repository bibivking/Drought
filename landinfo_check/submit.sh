#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q express
#PBS -l walltime=4:00:00
#PBS -l mem=190GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module use /g/data/hh5/public/modules 
module load conda/analysis3-unstable
cd /g/data/w97/mm3972/scripts/Drought/landinfo_check
python compare_soil_dataset.py
