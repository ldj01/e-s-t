#!/bin/bash

directory=$1

year=$2
month=$3
day=$4
date=$year$month$day

if [[ $year -le  1992 ]]; then
   var=100
elif [[ ($year -ge 1993) && ($year -le 2000) ]]; then
   var=200
else
  var=300
fi

echo $var

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies https://goldsmr3.sci.gsfc.nasa.gov:443/opendap/MERRA/MAI3CPASM.5.2.0/${year}/${month}/MERRA${var}.prod.assim.inst3_3d_asm_Cp.${date}.hdf.nc4 -O ${directory}/merra.nc4

#wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies https://goldsmr3.sci.gsfc.nasa.gov/opendap/MERRA/MAI3CPASM.5.2.0/$year/$month/MERRA$var.prod.assim.inst3_3d_asm_Cp.$date.hdf -O $directory/merra.hdf

#wget ftp://goldsmr3.sci.gsfc.nasa.gov/data/s4pa/MERRA/MAI3CPASM.5.2.0/$year/$month/MERRA$var.prod.assim.inst3_3d_asm_Cp.$date.hdf -O $directory/merra.hdf
