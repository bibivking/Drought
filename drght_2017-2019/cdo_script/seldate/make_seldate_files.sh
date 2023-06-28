case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB" 
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
case_type="sen"

var_names=("Qair_f_inst" "AvgSurfT_tavg" "Wind_f_inst" "Swnet_tavg" "Lwnet_tavg" "GPP_tavg" "NPP_tavg" "VegT_tavg" "Landcover_inst" "Albedo_inst" "LAI_inst" "Qle_tavg" "Qh_tavg" "Qg_tavg" "Evap_tavg" "TVeg_tavg" "ESoil_tavg" "FWsoil_tavg" "WaterTableD_tavg" "GWwb_tavg" "SoilMoist_inst" "Rainf_tavg" "Tair_f_inst") 


for var_name in "${var_names[@]}"; do

cat > ./selvar_LIS_output_${case_type}_${var_name}.pbs << EOF_selvar
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

cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/LIS_output/${var_name}

cdo seldate,2019-12-01,2020-02-29 LIS.CABLE.201701-202002.nc LIS.CABLE.201912-202002.nc

module use /g/data/hh5/public/modules
module load conda/analysis3-22.04
nccompress -r -o -np 8 LIS.CABLE.201912-202002.nc

EOF_selvar

qsub ./selvar_LIS_output_${case_type}_${var_name}.pbs

done
