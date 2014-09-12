#!/bin/bash 

home=`pwd`
echo home = $home

directory="temp"
if [ ! -d "$home/$directory" ]; then
mkdir $home/$directory 
fi

year=$1
month=$2
day=$3
hour=$4
min=$5
sec=$6

if [ ${month} -lt 10 ]
then
	month=0${month}
fi
if [ ${day} -lt 10 ]
then
	day=0${day}
fi
echo ${month}
echo ${day}

#  determine three hour increment before and after acquisition time
rem1=$(($hour % 3))
rem2=$((3-$rem1))
hour1=$(($hour-$rem1))
hour2=$(($hour+$rem2))
echo hour1 = $hour1
echo hour2 = $hour2
if [ ${hour1} -lt 10 ]
then
	hour1=0${hour1}
fi
if [ ${hour2} -lt 10 ]
then
	hour2=0${hour2}
fi
echo ${hour1}
echo ${hour2}

#  define variables for scripts to obtain NARR GRIB files
yearmo=$year$month
date=$year$month$day
echo $yearmo
echo $date

#  define files for scripts to obtain NARR GRIB files containing the geopotential height variable
fileHGT1='script_HGT_'$hour1
echo fileHGT1 = $fileHGT1
fileHGT2='script_HGT_'$hour2
echo fileHGT2 = $fileHGT2

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the geopotential height variable
cat script_HGT_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $directory/$fileHGT1
cat script_HGT_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $directory/$fileHGT2

#  define files for scripts to obtain NARR GRIB files containing the specific humidity variable
fileSHUM1='script_SHUM_'$hour1
echo fileSHUM1 = $fileSHUM1
fileSHUM2='script_SHUM_'$hour2
echo fileSHUM2 = $fileSHUM2

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the specific humidity variable
cat script_SHUM_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $directory/$fileSHUM1
cat script_SHUM_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $directory/$fileSHUM2

#  define files for scripts to obtain NARR GRIB files containing the temperature variable
fileTMP1='script_TMP_'$hour1
echo fileTMP1 = $fileTMP1
fileTMP2='script_TMP_'$hour2
echo fileTMP2 = $fileTMP2

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the temperature variable
cat script_TMP_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $directory/$fileTMP1
cat script_TMP_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $directory/$fileTMP2

#  copy folder with GRIB files to directory
copy=${home}/GRIB
destination=${directory}/GRIB
cp -r $copy $destination

#  change permissions on script files
chmod 755 $directory/$fileHGT1
chmod 755 $directory/$fileTMP1
chmod 755 $directory/$fileSHUM1

#  copy script files for time before to appropriate GRIB directory
mv $directory/$fileHGT1 $destination/script_HGT
mv $directory/$fileTMP1 $destination/script_TMP
mv $directory/$fileSHUM1 $destination/script_SHUM

#  change to appropriate GRIB directory and run scripts
cd $destination
if [ ! -d "HGT" ]; then
mkdir HGT 
fi
if [ ! -d "SHUM" ]; then
mkdir SHUM 
fi
if [ ! -d "TMP" ]; then
mkdir TMP 
fi
./script_HGT
./script_SHUM
./script_TMP
./HGT_grb2txt
./SHUM_grb2txt
./TMP_grb2txt

#  change back to home directory
cd $home

# make directories to contain results for this specific landsat image
mkdir $directory/HGT_1
mkdir $directory/TMP_1
mkdir $directory/SHUM_1

#  copy results from GRIB directory to directory for this specific Landsat image
mv $home/$destination/HGT/* $directory/HGT_1
mv $home/$destination/TMP/* $directory/TMP_1
mv $home/$destination/SHUM/* $directory/SHUM_1

#  change permissions on script files
chmod 755 $directory/$fileHGT2
chmod 755 $directory/$fileTMP2
chmod 755 $directory/$fileSHUM2

#  copy script files for time after to appropriate GRIB directory
mv $directory/$fileHGT2 $destination/script_HGT
mv $directory/$fileTMP2 $destination/script_TMP
mv $directory/$fileSHUM2 $destination/script_SHUM

#  change to appropriate GRIB directory and run scripts
cd $destination
./script_HGT
./script_SHUM
./script_TMP
./HGT_grb2txt
./SHUM_grb2txt
./TMP_grb2txt

#  change back to home directory
cd $home

# make directories to contain results for this specific landsat image
mkdir $directory/HGT_2
mkdir $directory/TMP_2
mkdir $directory/SHUM_2

#  copy results from GRIB directory to directory for this specific Landsat image
mv $home/$destination/HGT/* $directory/HGT_2
mv $home/$destination/TMP/* $directory/TMP_2
mv $home/$destination/SHUM/* $directory/SHUM_2

echo DONE DOWNLOADING NARR

date
