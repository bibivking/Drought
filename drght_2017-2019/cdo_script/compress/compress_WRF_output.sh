case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
case_short_name="ctl"
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB" 
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"

var_names=("cape_2d" "cloudfrac" "p" "pw" "rh" "rh2" "slp" "T2" "tc" "td2" "th" "ua" "va")

for var_name in "${var_names[@]}"; do

cat > ./compress_WRF_${case_short_name}_${var_name}.pbs << EOF_compress
#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalbw
#PBS -l walltime=2:30:00
#PBS -l mem=256GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/hh5

module load cdo

cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/WRF_output/${var_name}

module use /g/data/hh5/public/modules 
module load conda/analysis3-22.04
nccompress -r -o -np 8 /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/WRF_output/${var_name}/wrfout_201701-202002.nc

EOF_compress

qsub ./compress_WRF_${case_short_name}_${var_name}.pbs

done

