#!/bin/bash 

##	DEFINING HOME DIRECTORY AND INPUTS

#  This program must be run from directory containing the first input and the necessary supporting files and programs
#  First input is directory containing Landsat metadata and the directory where results will be put
#  Second input in the basename of the Landsat files

date
hostname

home="$(pwd)"
echo "home = $home"

directory="$1"
echo "directory = $directory"

imageBase=$2
echo "imageBase = $imageBase"

#height=$3

# make scene list
#scenelistname='scenelist'
#echo "$imageBase" >> "$scenelistname.txt"

#  define Landsat metadata file
metaFile="${directory}/${imageBase}_MTL.txt"
echo metaFile = "$metaFile"

cat "$metaFile" | grep -a = | sed s/\ =\ /=/g > "${directory}/${imageBase}_meta.txt"
source "${directory}/${imageBase}_meta.txt"

##	PARSE NECESSARY DATA FROM LANDSAT METADATA

# parse path and row from Landsat metadata
path="$(echo ${WRS_PATH} | bc)"
row="$(echo ${WRS_ROW} | bc)"

echo "path = $path"
echo "row = $row"

#  parse year, month, and day of acquisition from Landsat metadata
year="$(echo "${DATE_ACQUIRED}" | awk -F'-' '{print $1}')"
month="$(echo "${DATE_ACQUIRED}" | awk -F'-' '{print $2}')"
day="$(echo "${DATE_ACQUIRED}" | awk -F'-' '{print $3}')" 

echo "year = $year"
echo "month = $month"
echo "day = $day"

#  parse hour, minute, and second of acquisition from Landsat metadata
hour="$(echo "${SCENE_CENTER_TIME}" | awk -F':' '{print $1}')"
min="$(echo "${SCENE_CENTER_TIME}" | awk -F':' '{print $2}')"
sec="$(echo "${SCENE_CENTER_TIME}" |  awk -F':' '{print $3}'| sed 's/Z//')"

echo "hour = $hour"
echo "min = $min"
echo "sec = $sec"

#  parse number of samples and number of lines in the image from Landsat metadata
landsatSamples="${THERMAL_SAMPLES}"
landsatLines="${THERMAL_LINES}"

echo "landsatSamples = $landsatSamples"
echo "landsatLines = $landsatLines"

#  parse pixel size from Landsat metadata
pixelSize="${GRID_CELL_SIZE_THERMAL}"

echo "pixelSize = $pixelSize"

#  parse latitude, longitude coordinates of image corners from Landsat metadata
UL_LAT="${CORNER_UL_LAT_PRODUCT}"
UL_LON="${CORNER_UL_LON_PRODUCT}"
UR_LAT="${CORNER_UR_LAT_PRODUCT}"
UR_LON="${CORNER_UR_LON_PRODUCT}"
LL_LAT="${CORNER_LL_LAT_PRODUCT}"
LL_LON="${CORNER_LL_LON_PRODUCT}"
LR_LAT="${CORNER_LR_LAT_PRODUCT}"
LR_LON="${CORNER_LR_LON_PRODUCT}"

echo "UL_LAT = $UL_LAT"
echo "UL_LON = $UL_LON"
echo "UR_LAT = $UR_LAT"
echo "UR_LON = $UR_LON"
echo "LL_LAT = $LL_LAT"
echo "LL_LON = $LL_LON"
echo "LR_LAT = $LR_LAT"
echo "LR_LON = $LR_LON"

#  parse UTM coordinates of image corners from Landsat metadata
UL_EAST="${CORNER_UL_PROJECTION_X_PRODUCT}"
UL_NORTH="${CORNER_UL_PROJECTION_Y_PRODUCT}"
UR_EAST="${CORNER_UR_PROJECTION_X_PRODUCT}"
UR_NORTH="${CORNER_UR_PROJECTION_Y_PRODUCT}"
LL_EAST="${CORNER_LL_PROJECTION_X_PRODUCT}"
LL_NORTH="${CORNER_LL_PROJECTION_Y_PRODUCT}"
LR_EAST="${CORNER_LR_PROJECTION_X_PRODUCT}"
LR_NORTH="${CORNER_LR_PROJECTION_Y_PRODUCT}"

echo "UL_EAST = $UL_EAST"
echo "UL_NORTH = $UL_NORTH"
echo "UR_EAST = $UR_EAST"
echo "UR_NORTH = $UR_NORTH"
echo "LL_EAST = $LL_EAST"
echo "LL_NORTH = $LL_NORTH"
echo "LR_EAST = $LR_EAST"
echo "LR_NORTH = $LR_NORTH"

#  parse UTM zone from Landsat metadata
zone="${UTM_ZONE}"

echo "zone = $zone"

# determine which Landsat for current scene
whichLandsat="${SPACECRAFT_ID:(-1)}"



if [[ "${whichLandsat}" -eq 8 ]]
then
	whichLandsat=11
fi

echo "whichLandsat = $whichLandsat"

#  obtain values to convert digital to radiance for thermal band

case "${whichLandsat}" in
	5) LMAX6="${RADIANCE_MAXIMUM_BAND_6}"
  	   LMIN6="${RADIANCE_MINIMUM_BAND_6}"
	   QCALMAX6="${QUANTIZE_CAL_MAX_BAND_6}"
	   QCALMIN6="${QUANTIZE_CAL_MIN_BAND_6}"
	   ;;
	7) LMAX6="${RADIANCE_MAXIMUM_BAND_6_VCID_1}"
	   LMIN6="${RADIANCE_MINIMUM_BAND_6_VCID_1}"
	   QCALMAX6="${QUANTIZE_CAL_MAX_BAND_6_VCID_1}"
	   QCALMIN6="${QUANTIZE_CAL_MIN_BAND_6_VCID_1}"
	   ;;
	10) LMAX6="${RADIANCE_MAXIMUM_BAND_10}"
  	    LMIN6="${RADIANCE_MINIMUM_BAND_10}"
	    QCALMAX6="${QUANTIZE_CAL_MAX_BAND_10}"
	    QCALMIN6="${QUANTIZE_CAL_MIN_BAND_10}"
	    ;;
	11) LMAX6="${RADIANCE_MAXIMUM_BAND_11}"
  	    LMIN6="${RADIANCE_MINIMUM_BAND_11}"
	    QCALMAX6="${QUANTIZE_CAL_MAX_BAND_11}"
	    QCALMIN6="${QUANTIZE_CAL_MIN_BAND_11}"
	    ;;
esac

echo "LMAX6 = $LMAX6"
echo "LMIN6 = $LMIN6"
echo "QCALMAX6 = $QCALMAX6"
echo "QCALMIN6 = $QCALMIN6"


: <<'COMMENT'
##	DOWNLOAD NARR DATA

#  determine three hour increment before and after acquisition time
rem1="$(($hour % 3))"
rem2="$((3-$rem1))"
hour1="$(($hour-$rem1))"
hour2="$(($hour+$rem2))"
echo "hour1 = $hour1"
echo "hour2 = $hour2"
if [[ "${hour2}" -eq 24 ]]
then
        hour2=0
fi
if [[ "${hour1}" -lt 10 ]]
then
	hour1="0${hour1}"
fi
if [[ "${hour2}" -lt 10 ]]
then
	hour2="0${hour2}"
fi
echo "${hour1}"
echo "${hour2}"


#  define variables for scripts to obtain NARR GRIB files
yearmo="${year}${month}"
date="${year}${month}${day}"
echo "$yearmo"
echo "$date"



#  define files for scripts to obtain NARR GRIB files containing the geopotential height variable
fileHGT1="script_HGT_$hour1"
echo "fileHGT1 = $fileHGT1"
fileHGT2="script_HGT_$hour2"
echo "fileHGT2 = $fileHGT2"

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the geopotential height variable
cat script_HGT_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $directory/$fileHGT1
cat script_HGT_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $directory/$fileHGT2

#  define files for scripts to obtain NARR GRIB files containing the specific humidity variable
fileSHUM1="script_SHUM_$hour1"
echo "fileSHUM1 = $fileSHUM1"
fileSHUM2="script_SHUM_$hour2"

echo "fileSHUM2 = $fileSHUM2"

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the specific humidity variable
cat script_SHUM_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $directory/$fileSHUM1
cat script_SHUM_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $directory/$fileSHUM2

#  define files for scripts to obtain NARR GRIB files containing the temperature variable
fileTMP1="script_TMP_$hour1"
echo "fileTMP1 = $fileTMP1"
fileTMP2="script_TMP_$hour2"
echo "fileTMP2 = $fileTMP2"

#  modify generic scripts with acquisition dates and hours to obtain NARR GRIB files containing the temperature variable
cat script_TMP_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour1'/' > $directory/$fileTMP1
cat script_TMP_generic | sed 's/year/'$year'/' | sed 's/mo/'$month'/' | sed 's/dy/'$day'/' | sed 's/rh/'$hour2'/' > $directory/$fileTMP2

#  copy folder with GRIB files to directory
copy="${home}/GRIB"
destination="${directory}/GRIB"
cp -r "$copy" "$destination"

#  change permissions on script files
chmod 755 "$directory/$fileHGT1"
chmod 755 "$directory/$fileTMP1"
chmod 755 "$directory/$fileSHUM1"

#  copy script files for time before to appropriate GRIB directory
mv "${directory}/${fileHGT1}" "${destination}/script_HGT"
mv "${directory}/${fileTMP1}" "${destination}/script_TMP"
mv "${directory}/${fileSHUM1}" "$destination/script_SHUM"

#  change to appropriate GRIB directory and run scripts
cd "$destination"
./script_HGT
./script_SHUM
./script_TMP
./HGT_grb2txt
./SHUM_grb2txt
./TMP_grb2txt


#  change back to home directory
cd "$home"

# make directories to contain results for this specific landsat image
mkdir "$directory/HGT_1"
mkdir "$directory/TMP_1"
mkdir "$directory/SHUM_1"

#  copy results from GRIB directory to directory for this specific Landsat image
mv $home/$destination/HGT/* "$directory/HGT_1"
mv $home/$destination/TMP/* "$directory/TMP_1"
mv $home/$destination/SHUM/* "$directory/SHUM_1"

#  change permissions on script files
chmod 755 "$directory/$fileHGT2"
chmod 755 "$directory/$fileTMP2"
chmod 755 "$directory/$fileSHUM2"

#  copy script files for time after to appropriate GRIB directory
mv "$directory/$fileHGT2" "$destination/script_HGT"
mv "$directory/$fileTMP2" "$destination/script_TMP"
mv "$directory/$fileSHUM2" "$destination/script_SHUM"

#  change to appropriate GRIB directory and run scripts
cd "$destination"
./script_HGT
./script_SHUM
./script_TMP
./HGT_grb2txt
./SHUM_grb2txt
./TMP_grb2txt

#  change back to home directory
cd "$home"

# make directories to contain results for this specific landsat image
mkdir "$directory/HGT_2"
mkdir "$directory/TMP_2"
mkdir "$directory/SHUM_2"

#  copy results from GRIB directory to directory for this specific Landsat image
mv $home/$destination/HGT/* "$directory/HGT_2"
mv $home/$destination/TMP/* "$directory/TMP_2"
mv $home/$destination/SHUM/* "$directory/SHUM_2"

echo "DONE DOWNLOADING NARR"

date

#  define albedo at which modtran was run
alb=0.1
#  define number of height modtran was run at at each NARR point
numHeights=9

source /home/kga1099/.bashrc
module load envi 
pwd

#  call IDL program to generate tape5 file and commandList
idl -e "LST_NARR_STEP1, '$home/', '$directory/', '$imageBase', $year, $month, $day, $hour, $min, $sec,\
	              $landsatSamples, $landsatLines, $pixelSize, $zone,\
		      $UL_EAST, $UL_NORTH, $UR_EAST, $UR_NORTH, $LL_EAST, $LL_NORTH, $LR_EAST, $LR_NORTH,\
		      $UL_LAT, $UL_LON, $UR_LAT, $UR_LON, $LL_LAT, $LL_LON, $LR_LAT, $LR_LON"

echo "DONE GENERATING FILES"

date 

cmdLst="${directory}/commandList"
#  change permission on commandList
chmod 755 "$cmdLst"

cseLst="${directory}/caseList"
#  change permission on caseList
chmod 755 "$cseLst"


#  perform modtran runs by calling commandList
$cmdLst

echo "MODTRAN RUNS COMPLETED"

date

echo "PARSING TAPE6 FILES"

#  for each case in caseList (for each modtran run), copy program to delete headers and parse wavelength
#  and total radiance from tape6 file
for CASE in $(cat $cseLst); do

	#Create link.  If create fails, loop and sleep until successful (don't forget to yell!)
	if ! ln ~/elim2.sed "${CASE}"; then
                while [[ ! -e "${CASE}/elim2.sed" ]] ; do 
                         logger "LINK CREATE FAILED, case=${CASE} pwd=`pwd`" 
                         sleep 5 ; ln ~/elim2.sed "${CASE}"
                done
        fi

	./tape6parser.bash ${CASE}
done

echo "TAPE6s PARSED"

date

echo "GENERATING NARR PARAMETERS"

COMMENT

cseLst="${directory}/caseList"
#  determine number of cases in caseList
numCases="$(wc $cseLst | awk '{print $1}')"
#  determine number of NARR points in current Landsat scene
alb=0.1
numHeights=9
numPoints="$((${numCases}/3/${numHeights}))"
module load envi

#  call IDL program to generate parameters for each height and NARR point
idl -e "LST_NARR_STEP2, '$home/', '$directory/', $numPoints, $numHeights, $alb, $whichLandsat"

echo "NARR PARAMETERS GENERATED"

date

echo "GENERATING PIXEL PARAMETERS"

demname=DEM/"${path}_${row}"/*.tif
demFile="$demname"
#demFile='DEM/'$path'_'$row'/'$path'_'$row'_DEM.tif'

#  call IDL program to generate parameters for each pixel
idl -e "LST_NARR_STEP3, '$home/', '$directory/', '$imageBase', $numPoints, $numHeights,\
		            $landsatSamples, $landsatLines, $pixelSize, $zone,\
      		            $UL_EAST, $UL_NORTH, $UR_EAST, $UR_NORTH, $LL_EAST, $LL_NORTH, $LR_EAST, $LR_NORTH,\
                	    $UL_LAT, $UL_LON, $UR_LAT, $UR_LON, $LL_LAT, $LL_LON, $LR_LAT, $LR_LON,\
                            $LMAX6, $LMIN6, $QCALMAX6, $QCALMIN6, '$demFile', $whichLandsat"

echo "PIXEL PARAMETERS GENERATED"

echo "REMOVING INTERMEDIATE FILES"

#if [ -a ${home}/${directory}/${imagebase}_LSTresults.tif ]
#then
#rm -r $directory/HGT*
#rm -r $directory/SHUM*
#rm -r $directory/TMP*
#rm -r $directory/GRIB
#rm -r $directory/[0-9]*
#rm $directory/atmosphericParameters.txt
#rm $directory/new*
#rm $directory/tempLayers.txt
#rm $directory/caseList
#rm $directory/commandList

#else
#B
#   exit 1
#fi

echo "FILES REMOVED"

date
