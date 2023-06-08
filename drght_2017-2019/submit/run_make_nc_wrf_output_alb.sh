#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -l walltime=30:00:00
#PBS -l mem=256GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normalbw
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-22.04

var_names=("td2" "T2" "rh2" "slp" "pw" "cloudfrac" "cape_2d" "p" "th" "rh" "tc" "ua" "va")

for var_name in "${var_names[@]}"; do

case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/WRF_output/
mkdir ${var_name}

cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019
python ./produce_nc_wrf_output.py -n ${var_name} -c ${case_name}

done
