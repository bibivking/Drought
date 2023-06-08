case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"
case_short_name="ctl"
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB" 
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"

var_names=("Qair_f_inst" "AvgSurfT_tavg" "Wind_f_inst" "Swnet_tavg" "Lwnet_tavg" "GPP_tavg" "NPP_tavg" "VegT_tavg" "Landcover_inst" "Albedo_inst" "LAI_inst" "Qle_tavg" "Qh_tavg" "Qg_tavg" "Evap_tavg" "TVeg_tavg" "ESoil_tavg" "FWsoil_tavg" "WaterTableD_tavg" "GWwb_tavg" "SoilMoist_inst" "Rainf_tavg" "Tair_f_inst")

for var_name in "${var_names[@]}"; do

cat > ./compress_LIS_output_${case_short_name}_${var_name}.pbs << EOF_compress
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

cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/LIS_output/${var_name}
rm *.d01.nc

module use /g/data/hh5/public/modules 
module load conda/analysis3-22.04
nccompress -r -o -np 8 /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/LIS_output/${var_name}/LIS.CABLE.201701-202002.nc

EOF_compress

qsub ./compress_LIS_output_${case_short_name}_${var_name}.pbs

done

