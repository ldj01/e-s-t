#!/bin/bash

main_folderpath=$1
cloudmask_path=$2
cloudmask_name=$3
lst_sst_data=$4
results_filename=$5

lsat_scenename=${cloudmask_name:0:21}

IFS=$'\n' read -d'' -r -a lines < "${lst_sst_data}"
count=0

for line in "${lines[@]}"
  do
  words=( $line )
  landsat="${words[0]}"
  modis="${words[1]}"
  which_point="${words[2]}"
  lsat_lst="${words[3]}"
  modis_sst="${words[4]}"
  easting="${words[5]}"
  northing="${words[6]}"
  qual="${words[7]}"
  stdev="${words[8]}"
  numzeros="${words[9]}"
  lobs="${words[10]}"
  trans="${words[11]}"
  upwell="${words[12]}"
  downwell="${words[13]}"
  dist_output_file="${main_folderpath}${lsat_scenename}_pt${which_point}.txt"

  if [[ $count -gt 0 ]]
    then
    if [[ ${lsat_scenename} == ${landsat} ]]
      then
      
      if [[ ${easting:0:1} -eq 0 ]]
        then
        echo "setting distance to 0..."
        distance=0
      else       
        idl -e "DIST2CLOUD_MAIN, '${cloudmask_path}', '${cloudmask_name}', '${easting}','${northing}','${dist_output_file}"
        
        distance=$( cat "${dist_output_file}" )
        rm ${dist_output_file}

      fi
      
     echo "${lsat_scenename}, ${which_point}, ${easting}, ${northing}, ${lsat_lst}, ${modis_sst}, ${qual}, ${stdev}, ${numzeros}, ${lobs}, ${trans}, ${upwell}, ${downwell}, ${distance}" >> "${main_folderpath}${results_filename}"    

    fi

  fi
  ((count++))
done


