#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=190GB
#PBS -l ncpus=8
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97+scratch/w97+gdata/w97

module use /g/data/hh5/public/modules 
module load conda/analysis3-22.04
nccompress -r -o -np 8 /g/data/w97/mm3972/scripts/Drought/drght_2017-2019/nc_files/deltaTmax_deltaLAI_deltaAlbedo_2017_2020_DJF.nc

