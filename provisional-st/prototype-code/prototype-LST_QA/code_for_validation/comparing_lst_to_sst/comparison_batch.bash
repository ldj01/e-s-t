#!/bin/bash

# COMPARISON_BATCH.BASH #
#########################
#
# This is a batch file that will read in a list of Landsat and corresponding
# MODIS scenes, then use IDL scripts to find up to five points for each
# scene pair that would be appropriate for comparing the Landsat Land Surface 
# Temperature (LST) and the MODIS Sea Surface Temperature (SST). The results
# get written out to a text file.


#
# Define the follwing things:
#
# 1) the "path_row" being processed (helps define other paths)
# 2) the path to the main folder where output files should go
# 3) the path and name of the .txt file with Lsat and MODIS scene names
# 4) the path where the Landsat LST results files are
# 5) the path where the georeferenced sst and quality images are
# 6) the path and name of the Landsat LUT
# 7) the path where the temporary output files will go (will be deleted)
# 8) the name of the output file, which will be written to the main folder
#

path_row=$1
year=$2

#path_row="107_19"
#year="2013"
main_folderpath="/dirs/home/phd/kga1099/globalL7"
modis_lsat_lut="${main_folderpath}/scenelists/${path_row}/modis_scenenames_${path_row}_${year}.txt"
lsat_folderpath="/dirs/home/phd/kga1099/globalL7/landsat_files/${path_row}/"
modis_folderpath="/dirs/home/phd/kga1099/globalL7/modis_files/${path_row}/${year}/"
lut_path="${main_folderpath}/LUT7_new.txt"
datadump_path="${main_folderpath}/"
results_filename="last_study_${path_row}_${year}.txt"

#
# Write column headers to output file
#
echo "LSAT_SCENE             MODIS_SCENE    NO.  LSAT_LST  MODIS_SST   Easting    Northing    Qual   StdevBest    Numzeros   lobs    trans    upwell    downwell" >> ${results_filename}


#
# Read in the .txt file with Landsat and MODIS scenenames
#
IFS=$'\n' read -d'' -r -a lines < "${modis_lsat_lut}"
count=0

#
# For each line in the input file, grab the scenenames and run the
# COMPARE_MODIS_LSAT_GLOBAL procedure
#
for line in "${lines[@]}"
  do
  words=( $line )
  lsat_scenename="${words[0]}"
  modis_scenename="${words[1]}"
  if [[ count -gt 0 ]]
    then
    idl -e "COMPARE_MODIS_LSAT_GLOBAL, '${lsat_folderpath}', '${modis_folderpath}', '${lsat_scenename}', '${modis_scenename}', '${lut_path}','${datadump_path}'"
    ifs=$'\n' read -d'' -r -a data < "${datadump_path}${lsat_scenename}_lst_sst_data.txt"
    num_entries=${#data[@]}
    num=0
    range=$(( ${num_entries}/11 - 1 ))
    while [[ num -le $range ]]
      do
      entry1=$(( 11*($num) ))
      entry2=$(( 11*($num) + 1 ))
      entry3=$(( 11*($num) + 2 ))
      entry4=$(( 11*($num) + 3 ))
      entry5=$(( 11*($num) + 4 )) 
      entry6=$(( 11*($num) + 5 ))
      entry7=$(( 11*($num) + 6 ))
      entry8=$(( 11*($num) + 7 ))
      entry9=$(( 11*($num) + 8 ))
      entry10=$(( 11*($num) + 9 ))
      entry11=$(( 11*($num) + 10 ))
      
      echo "${lsat_scenename}  ${modis_scenename} ${num}    ${data[$entry1]}   ${data[$entry2]} ${data[$entry3]} ${data[$entry4]} ${data[$entry5]} ${data[$entry6]} ${data[$entry7]} ${data[$entry8]} ${data[$entry9]} ${data[$entry10]} ${data[$entry11]}" >> ${results_filename} 
      ((num++))
    done
    rm "${datadump_path}${lsat_scenename}_lst_sst_data.txt" 
  fi

  ((count++))
done


