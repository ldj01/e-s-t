#!/bin/bash

directory=$1
file1=$2
file2=$3

# cd to directory
cd "$directory"

# define which grb2 message lines contain geopotential height for desired pressure levels
hgtlines=( 2 10 18 26 34 42 51 60 69 78 87 96 105 114 123 132 141 150 159 168 177 186 195 204 213 222 231 240 249 258 267 276 285 294 303 312 321 )

# grib message lines where temperature at desired pressure levels
tmplines=( 3 11 19 27 35 43 52 61 70 79 88 97 106 115 124 133 142 151 160 169 178 187 196 205 214 223 232 241 250 259 268 277 286 295 304 313 322 )

# grib lines where relative humidity is
# rhlines=( 4 12 20 28 36 44 53 62 71 80 89 98 107 116 125 134 143 152 161 170 179 188 197 206 215 224 233 242 251 260 269 278 287 296 305 314 323 )

# grib lines where specific humidity is
shumlines=( 5 13 21 29 37 45 54 63 72 81 90 99 108 117 126 135 144 153 162 171 180 189 198 207 216 225 234 243 252 261 270 279 288 297 306 315 324 ) 

# define pressure levels
plevels=( 0001 0002 0003 0005 0007 0010 0020 0030 0050 0070 0100 0125 0150 0175 0200 0225 0250 0300 0350 0400 0450 0500 0550 0600 0650 0700 0750 0775 0800 0825 0850 0875 0900 0925 0950 0975 1000 )

# make directories
mkdir $directory/HGT_1 $directory/HGT_2 $directory/TMP_1 $directory/TMP_2 $directory/SHUM_1 $directory/SHUM_2


# convert HGTs to txt
count=-1
for line in ${hgtlines[@]}
do
      ((count++))
      p="${plevels[$count]}"
      echo "pressure= $p"
      filename="HGT1_${p}mbar.txt"
      ./wgrib2 -d $line -text "$filename" "$file1"
      filename="HGT2_${p}mbar.txt"
      ./wgrib2 -d $line -text "$filename" "$file2"
done
mv HGT1*.txt $directory/HGT_1/
mv HGT2*.txt $directory/HGT_2/


# convert TMPs to txt
count=-1
for line in ${tmplines[@]}
do
      ((count++))
      p="${plevels[$count]}"
      echo "pressure= $p"
      filename="TMP1_${p}mbar.txt"
      ./wgrib2 -d $line -text "$filename" "$file1"
      filename="TMP2_${p}mbar.txt"
      ./wgrib2 -d $line -text "$filename" "$file2"
done
mv TMP1_*.txt $directory/TMP_1/
mv TMP2_*.txt $directory/TMP_2/


# convert RHs to txt
#count=-1
#for line in ${rhlines[@]}
#do
#      count=`expr $count + 1`
#      p=${plevels[$count]}
#      echo pressure= $p
#      filename=RH1_${p}mbar.txt
#      ./wgrib2 -d $line -text "$filename" pgbhnl.gdas."$file1".grb2
#      filename=RH2_${p}mbar.txt
#      ./wgrib2 -d $line -text "$filename" pgbhnl.gdas."$file2".grb2
#done
#mv RH1_*.txt ./RH_1/
#mv RH2_*.txt ./RH_2/


# convert SHUMs to txt
count=-1
for line in ${shumlines[@]}
do
      ((count++))
      p="${plevels[$count]}"
      echo "pressure= $p"
      filename="SHUM1_${p}mbar.txt"
      ./wgrib2 -d $line -text "$filename" "$file1"
      filename="SHUM2_${p}mbar.txt"
      ./wgrib2 -d $line -text "$filename" "$file2"
done
mv SHUM1_*.txt $directory/SHUM_1/
mv SHUM2_*.txt $directory/SHUM_2/

