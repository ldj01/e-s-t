#!/bin/bash 

home=`pwd`
echo home = $home

#directory="temp"
#if [ ! -d "$home/$directory" ]; then
#mkdir $home/$directory 
#fi

# copy needed files from $BIN directory
cp $BIN/script_HGT_generic .
cp $BIN/script_SHUM_generic .
cp $BIN/script_TMP_generic .
cp $BIN/get_grib.pl .
cp $BIN/get_inv.pl .
cp $BIN/HGT_grb2txt .
cp $BIN/SHUM_grb2txt .
cp $BIN/TMP_grb2txt .
cp $BIN/wgrib .

while read line
do
year=$(echo $line|cut -d ',' -f 1)
month=$(echo $line|cut -d ',' -f 2)
day=$(echo $line|cut -d ',' -f 3)
hour=$(echo $line|cut -d ',' -f 4)
min=$(echo $line|cut -d ',' -f 5)
sec=$(echo $line|cut -d ',' -f 6)
done < "datetime.txt"

echo year = $year
echo sec = $sec

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

# change to the temp directory

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the geopotential height variable
cat script_HGT_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $fileHGT1
cat script_HGT_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $fileHGT2

#  define files for scripts to obtain NARR GRIB files containing the specific humidity variable
fileSHUM1='script_SHUM_'$hour1
echo fileSHUM1 = $fileSHUM1
fileSHUM2='script_SHUM_'$hour2
echo fileSHUM2 = $fileSHUM2

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the specific humidity variable
cat script_SHUM_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $fileSHUM1
cat script_SHUM_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $fileSHUM2

#  define files for scripts to obtain NARR GRIB files containing the temperature variable
fileTMP1='script_TMP_'$hour1
echo fileTMP1 = $fileTMP1
fileTMP2='script_TMP_'$hour2
echo fileTMP2 = $fileTMP2

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the temperature variable
cat script_TMP_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $fileTMP1
cat script_TMP_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $fileTMP2

#  change permissions on script files
chmod 755 $fileHGT1
chmod 755 $fileTMP1
chmod 755 $fileSHUM1

#  copy script files for time before to appropriate GRIB directory

#mv $directory/$fileHGT1 $destination/script_HGT
#mv $directory/$fileTMP1 $destination/script_TMP
#mv $directory/$fileSHUM1 $destination/script_SHUM

#  change to appropriate GRIB directory and run scripts
#cd $destination
if [ ! -d "HGT" ]; then
mkdir HGT 
fi
if [ ! -d "SHUM" ]; then
mkdir SHUM 
fi
if [ ! -d "TMP" ]; then
mkdir TMP 
fi
./$fileHGT1
./$fileSHUM1
./$fileTMP1
./HGT_grb2txt
./SHUM_grb2txt
./TMP_grb2txt

# make directories to contain results for this specific landsat image
mkdir $home/HGT_1
mkdir $home/TMP_1
mkdir $home/SHUM_1

#  copy results from GRIB directory to directory for this specific Landsat image
mv $home/HGT/* $home/HGT_1
mv $home/TMP/* $home/TMP_1
mv $home/SHUM/* $home/SHUM_1

#  change permissions on script files
chmod 755 $fileHGT2
chmod 755 $fileTMP2
chmod 755 $fileSHUM2

#  copy script files for time after to appropriate GRIB directory
#mv $directory/$fileHGT2 $destination/script_HGT
#mv $directory/$fileTMP2 $destination/script_TMP
#mv $directory/$fileSHUM2 $destination/script_SHUM

#  change to appropriate GRIB directory and run scripts
#cd $destination
./$fileHGT2
./$fileSHUM2
./$fileTMP2
./HGT_grb2txt
./SHUM_grb2txt
./TMP_grb2txt

#  change back to home directory
cd $home

# make directories to contain results for this specific landsat image
mkdir $home/HGT_2
mkdir $home/TMP_2
mkdir $home/SHUM_2

#  copy results from GRIB directory to directory for this specific Landsat image
mv $home/HGT/* $home/HGT_2
mv $home/TMP/* $home/TMP_2
mv $home/SHUM/* $home/SHUM_2

echo DONE DOWNLOADING NARR

date
