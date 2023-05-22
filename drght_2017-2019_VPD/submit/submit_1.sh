#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q expressbw
#PBS -l walltime=2:20:00
#PBS -l mem=190GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/rt52+gdata/zz93+gdata/hh5+gdata/w97+scratch/w97+gdata/ua8

module load ncl
cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019_VPD
ncl Fig1a_time_serial_txtout_GRACE_vs_CABLE.ncl
