#!/bin/bash

directory=$1

year=$2
month=$3
day=$4
date=$year$month$day

if [[ $year -le  1991 ]]; then
   var=100
elif [[ ($year -ge 1992) && ($year -le 2000) ]]; then
   var=200
elif [[ ($year -ge 2001) && ($year -le 2010) ]]; then
  var=300
else
  var=400
fi

echo $var
echo $date

wget --user=<TBD> --password=<TBD> --no-clobber --no-check-certificate https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4/${year}/${month}/MERRA2_${var}.inst3_3d_asm_Np.${date}.nc4 -O ${directory}/merra2_${date}.nc4

