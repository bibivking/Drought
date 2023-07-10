case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
case_type="ctl"

var_names=("wa" "PBLH" "td2" "T2" "rh2" "slp" "pw" "cloudfrac" "cape_2d" "p" "th" "rh" "tc" "ua" "va")


for var_name in "${var_names[@]}"; do

cat > ./selvar_WRF_output_${case_type}_${var_name}.pbs << EOF_selvar
#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalsr
#PBS -l walltime=2:30:00
#PBS -l mem=256GB
#PBS -l ncpus=8
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97+gdata/hh5

module load cdo

cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/WRF_output/${var_name}

cdo seldate,2017-12-01,2018-02-28 wrfout_201701-202002.nc wrfout_201712-201802.nc
cdo seldate,2018-12-01,2019-02-28 wrfout_201701-202002.nc wrfout_201812-201902.nc

module use /g/data/hh5/public/modules
module load conda/analysis3-22.04
nccompress -r -o -np 8 wrfout_201712-201802.nc
nccompress -r -o -np 8 wrfout_201812-201902.nc

EOF_selvar

qsub ./selvar_WRF_output_${case_type}_${var_name}.pbs

done
