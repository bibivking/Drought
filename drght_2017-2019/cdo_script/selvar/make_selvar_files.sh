case_name="drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB"
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2_obs_LAI_ALB" 
#"drght_2017_2019_bl_pbl2_mp4_ra5_sf_sfclay2"

var_names=("Qair_f_inst" "AvgSurfT_tavg" "Wind_f_inst" "Swnet_tavg" "Lwnet_tavg" "GPP_tavg" "NPP_tavg" "VegT_tavg" "Landcover_inst" "Albedo_inst" "LAI_inst" "Qle_tavg" "Qh_tavg" "Qg_tavg" "Evap_tavg" "TVeg_tavg" "ESoil_tavg" "FWsoil_tavg" "WaterTableD_tavg" "GWwb_tavg" "SoilMoist_inst" "Rainf_tavg" "Tair_f_inst") 


for var_name in "${var_names[@]}"; do

cat > ./selvar_LIS_output_sen_${var_name}.pbs << EOF_selvar
#!/bin/bash

#PBS -m ae
#PBS -P w97
#PBS -q normalbw
#PBS -l walltime=2:30:00
#PBS -l mem=256GB
#PBS -l ncpus=1
#PBS -j oe
#PBS -l wd
#PBS -l storage=gdata/w97

module load cdo

cd /g/data/w97/mm3972/model/wrf/NUWRF/LISWRF_configs/Tinderbox_drght_LAI_ALB/${case_name}/LIS_output

mkdir ${var_name}
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201701-201701.d01.nc ./${var_name}/LIS.CABLE.201701-201701.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201702-201702.d01.nc ./${var_name}/LIS.CABLE.201702-201702.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201703-201703.d01.nc ./${var_name}/LIS.CABLE.201703-201703.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201704-201704.d01.nc ./${var_name}/LIS.CABLE.201704-201704.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201705-201705.d01.nc ./${var_name}/LIS.CABLE.201705-201705.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201706-201706.d01.nc ./${var_name}/LIS.CABLE.201706-201706.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201707-201707.d01.nc ./${var_name}/LIS.CABLE.201707-201707.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201708-201708.d01.nc ./${var_name}/LIS.CABLE.201708-201708.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201709-201709.d01.nc ./${var_name}/LIS.CABLE.201709-201709.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201710-201710.d01.nc ./${var_name}/LIS.CABLE.201710-201710.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201711-201711.d01.nc ./${var_name}/LIS.CABLE.201711-201711.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201712-201712.d01.nc ./${var_name}/LIS.CABLE.201712-201712.d01.nc

cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201801-201801.d01.nc ./${var_name}/LIS.CABLE.201801-201801.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201802-201802.d01.nc ./${var_name}/LIS.CABLE.201802-201802.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201803-201803.d01.nc ./${var_name}/LIS.CABLE.201803-201803.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201804-201804.d01.nc ./${var_name}/LIS.CABLE.201804-201804.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201805-201805.d01.nc ./${var_name}/LIS.CABLE.201805-201805.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201806-201806.d01.nc ./${var_name}/LIS.CABLE.201806-201806.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201807-201807.d01.nc ./${var_name}/LIS.CABLE.201807-201807.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201808-201808.d01.nc ./${var_name}/LIS.CABLE.201808-201808.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201809-201809.d01.nc ./${var_name}/LIS.CABLE.201809-201809.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201810-201810.d01.nc ./${var_name}/LIS.CABLE.201810-201810.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201811-201811.d01.nc ./${var_name}/LIS.CABLE.201811-201811.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201812-201812.d01.nc ./${var_name}/LIS.CABLE.201812-201812.d01.nc

cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201901-201901.d01.nc ./${var_name}/LIS.CABLE.201901-201901.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201902-201902.d01.nc ./${var_name}/LIS.CABLE.201902-201902.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201903-201903.d01.nc ./${var_name}/LIS.CABLE.201903-201903.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201904-201904.d01.nc ./${var_name}/LIS.CABLE.201904-201904.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201905-201905.d01.nc ./${var_name}/LIS.CABLE.201905-201905.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201906-201906.d01.nc ./${var_name}/LIS.CABLE.201906-201906.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201907-201907.d01.nc ./${var_name}/LIS.CABLE.201907-201907.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201908-201908.d01.nc ./${var_name}/LIS.CABLE.201908-201908.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201909-201909.d01.nc ./${var_name}/LIS.CABLE.201909-201909.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201910-201910.d01.nc ./${var_name}/LIS.CABLE.201910-201910.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201911-201911.d01.nc ./${var_name}/LIS.CABLE.201911-201911.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.201912-201912.d01.nc ./${var_name}/LIS.CABLE.201912-201912.d01.nc

cdo selvar,lat,lon,time,${var_name} LIS.CABLE.202001-202001.d01.nc ./${var_name}/LIS.CABLE.202001-202001.d01.nc
cdo selvar,lat,lon,time,${var_name} LIS.CABLE.202002-202002.d01.nc ./${var_name}/LIS.CABLE.202002-202002.d01.nc

cd ${var_name}
cdo mergetime LIS.CABLE.??????-??????.d01.nc LIS.CABLE.201701-202002.nc

EOF_selvar

qsub ./selvar_LIS_output_sen_${var_name}.pbs

done
