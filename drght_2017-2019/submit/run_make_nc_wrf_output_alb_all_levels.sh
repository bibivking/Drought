#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -l walltime=24:00:00
#PBS -l mem=500GB
#PBS -l ncpus=8
#PBS -j oe
#PBS -q normalsr
#PBS -l wd
#PBS -l storage=gdata/hh5+gdata/w97

module use /g/data/hh5/public/modules
module load conda/analysis3-22.04

case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
var_names=("rh2")
#"wa" "z" "ua" "va" "th" "rh" "PBLH" "p" "td2" "T2" "rh2" "slp" "pw" "cloudfrac" "cape_2d" "tc"

for var_name in "${var_names[@]}"; do

# cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/WRF_output/
# mkdir ${var_name}

cd /g/data/w97/mm3972/scripts/Drought/drght_2017-2019
python ./produce_nc_wrf_output_all_levels.py -n ${var_name} -c ${case_name}

cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/WRF_output/${var_name}
rm -r tmp.nc_compress
nccompress -r -o -np 8 wrfout_201912-202002_hourly.nc

done
