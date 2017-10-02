#! /bin/bash
set -e
# Kelly G. Laraby
# PhD Candidate in Imaging Science
# Rochester Institute of Technology
# kga1099@rit.edu
#
#
# get_modis_for_L7_scenes.bash
# ----------------------------
# Main file that gets he Landsat corner coordinates and loops through up to 10 attempts to download the MODIS SST scene
# that encapsulates the Landsat image. For each attempt, the IDL program is_right_modis_scene.pro determines if the
# downloaded image is the correct one, and once it has been found all further download attempts cease and the information
# is recorded to text file. The arguments for this script are strings containing paths to the Landsat MTL files, the 
# desired download location, the Landsat scene currently being considered, and the text file that results get recorded in.
# 
#
# Background: Given a certain Landsat 7 scene, it is sometimes useful to get the MODIS SST image that captures
# the same area on the ground. Landsat 7 and MODIS have a very similar orbit, and the image capture times are usually
# about 20 minutes apart. MODIS SST images can be downloaded from https://oceancolor.gsfc.nasa.gov either by browsing
# manually or by using direct data access. Since it is not obvious by the MODIS scene identifier which one overlaps 
# the Landsat scene of interest, one would usually choose the browse method of downloading SST images. But if there
# are a large number of Landsat scenes that we want to find corresponding MODIS SST images for, it would be preferable
# to make a way to download the correct scenes automatically.
#
# The solution that is implemented is to consider the Landsat image capture time, then start by downloading a MODIS
# SST file that was captured soon after that (SST capture times are always between 10-30 minutes behind Landsat).
# There is an IDL program that will check if the Landsat scene falls within the SST image that was downloaded, and 
# if it does not the file is deleted and a new SST image that was captured 5 minutes later is downloaded. Up to 10
# attempts are allowed to find the correct image, and if it is found before then the program stops and keeps the 
# correct SST file. The name of the file is recorded along with the corresponding Landsat scene and Landsat's 
# corner coordinate information, which can be used to georeference the SST images later.

# this script will retrieve the appropriate MODIS scene for each L7 scene
# it will download until it finds the right one
# this runs out of the directory that the file exists in

# Required arguments
landsat_folder=$1   # path to where Landsat MTL file lives
download_folder=$2  # path to where MODIS SST files will be placed
landsat_scene=$3    # Landsat scene name that we will find the right SST file for
result_filename=$4  # path (including file name) where the Landsat and MODIS scenenames, and
                    # Landsat corner coordinates will be stored


# read in Landsat info using MTL
metaFile="${landsat_folder}/${landsat_scene}_MTL.txt"

cat "${metaFile}" | grep -a = | sed s/\ =\ /=/g > "${landsat_folder}/${landsat_scene}_meta.txt"
source "${landsat_folder}/${landsat_scene}_meta.txt"

# get image aquisition date
year="$(echo "${DATE_ACQUIRED}" | awk -F'-' '{print $1}')"
month="$(echo "${DATE_ACQUIRED}" | awk -F'-' '{print $2}')"
day="$(echo "${DATE_ACQUIRED}" | awk -F'-' '{print $3}')"

doy="${landsat_scene:13:3}"

# get image aquisition time
ttime="${SCENE_CENTER_TIME}" 
hour="${ttime:0:2}"
min="${ttime:3:2}"

# get zone and corner lat/lons
utm_zone="${UTM_ZONE}"
UL_LAT="${CORNER_UL_LAT_PRODUCT}"
UL_LON="${CORNER_UL_LON_PRODUCT}"
LR_LAT="${CORNER_LR_LAT_PRODUCT}"
LR_LON="${CORNER_LR_LON_PRODUCT}"


# url for direct data access
website="http://oceandata.sci.gsfc.nasa.gov/cgi/getfile"

# end of each SST file
file_end=".L2_LAC_SST.nc"

# based on L7 acquisition, find first MODIS time guess
# MODIS is always at least 10 minutes behind Landsat so we can skip the first few
rem1=$((10#$min % 5))
rem2=$((5 - rem1))
startmin=$((10#$min + rem2 + 10))
starthour=${hour}

# change minutes and hour if necessary (e.g. if 10 minutes after Landsat goes over 60 minutes,
# or if the hour and/or minutes are less than ten they need to have a zero in front in order
# to get the right file name)
if [[ $startmin -ge 60 ]]
  then
  startmin=$(($startmin - 60))
  starthour=$((10#$starthour + 1))
fi

if [[ ${#startmin} -lt 2 ]]
  then
  startmin="0${startmin}"
fi

if [[ ${#starthour} -lt 2 ]]
  then
  starthour="0${starthour}"
fi

modis_hour="${starthour}"
modis_min="${startmin}"

count=1
isright=0


echo "=============================================================================="
echo ' '
echo "  FINDING MODIS SCENE FOR ${landsat_scene} (acquired at ${hour}:${min})"
echo "___________________________________________________________________"
echo ' '

echo "modis_hour=$modis_hour"
echo "modis_min=$modis_min"


# Begin downloading loop
while isright=0
do

  echo ' '
  echo '-------------'
  echo "  ATTEMPT # ${count}:"
  echo '-------------'
  echo ' '

  # if 10 tries have been made, give up
  if [[ $count -eq 11 ]]
  then
    echo "I've tried 10 times and haven't found the right scene, so I gave up."
    idl_modisname='not found'
    break
  fi
  
  # downoad file
  filename="T${year}${doy}${modis_hour}${modis_min}00${file_end}"
  echo "    trying $filename ....."
  echo ' '
  wget -nv "${website}/${filename}" -O "${filename}"
  echo ' '

  # unzip and move to downloads folder
  idl_modisname=${filename}
  mv "${filename}" "${download_folder}/" 
  
  # run idl program to test if MODIS scene is correct for given Landsat scene
  # this function creates a temporary file that contains a 0 or 1 to indicate if the modis scene was correct
  idl -e "IS_RIGHT_MODIS_SCENE, '${download_folder}','${idl_modisname}','${UL_LAT}','${UL_LON}','${LR_LAT}','${LR_LON}'"
  echo ' '
  
  # sleep to make sure the file can be found (in the past there would sometimes be an error where the file made by the 
  # above IDL program could not be found because it hadn't finished creating it, so waiting 0.1 seconds eliminates 
  # this problem)
  sleep .1
  isright=$(cat "${idl_modisname}.txt")
  
  # if temporary file contains a 1, the downloaded file is kept
  # a zero results in the deletion of the file and increments the modis aquisition time by 5 minutes
  if [[ isright -eq 1 ]]
    then
    break
  else
    rm "${idl_modisname}.txt"
    rm "${download_folder}/${idl_modisname}"
    
    # increment modis aquisition time by 5 minutes (may require changing hour as well)
    if [[ $modis_min -eq 55 ]]
      then
      ((modis_hour++))
      modis_min="00"
     
     if [[ ${#modis_hour} -lt 2 ]]
       then
       modis_hour="0${modis_hour}"
     fi

    elif [[ $modis_min -eq 00 ]]
      then
      modis_min="05"
    else
      modis_min=$((modis_min + 5))
    fi
  fi
  
  ((count++))
done

# clean up
rm "${landsat_folder}/${landsat_scene}_meta.txt"
rm "${idl_modisname}.txt"

echo "The correct MODIS file for ${landsat_scene} is ${idl_modisname}"

# write results to text file
echo "${landsat_scene}     ${idl_modisname:0:14}     ${utm_zone}     ${UL_LAT}    ${UL_LON}    ${LR_LAT}    ${LR_LON}" >> "${result_filename}"
echo '=============================================================================='
