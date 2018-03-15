#!/bin/bash
#
# This test script is based on the RIT prototype script, but is modified as 
# follows:
#
# - It does not use SLURM, since we aren't using SLURM for ST processing.
#   We will only process 1 scene per directory so we don't need a job ID to 
#   distinguish files.
# - It does not use /dev/shm (shared memory), since we aren't using that for ST
# - It scans the MODTRAN output files for different strings, since the MODTRAN
#   version we are using gives somewhat different output
#

#  change to directory containing MODTRAN output of current case
pushd $1

#  determine line where pertinent output ends and delete everything below that
linenum1=`grep -n "MULTIPLE SCATTERING CALCULATION" tape6 |head -1| awk '{print $1}'| awk -F':' '{print $1}'`
sed "$linenum1,$ d" < tape6 > tempfile

chmod 777 tempfile

#  determine line where pertinent output begins and delete everything before that
linenum2=`grep -n "WAVELENGTH" tempfile |head -1| awk '{print $1}' |awk -F':' '{print $1}'`
sed "1,$linenum2 d" < tempfile > tempfile2

chmod 777 tempfile2

#  delete headers that divide data within file so data is continuous columns
sed -f elim2.sed < tempfile2 > tempfile3

chmod 777 tempfile3

#  determine if tape6 contains warning
linenum3=`grep -n "Warning" tempfile3 | head -1 | awk '{print $1}'|awk -F':' '{print $1}'`
linenum4=`grep -n "WARNING" tempfile3 | tail -1 | awk '{print $1}'|awk -F':' '{print $1}'`

#  if tape6 contains warning, parse warning and write to tempfile4
#  then print (wavelength, total radiance) to parsed file
#  else print (wavelength, total radiance) to parsed file

if [ $linenum3 > 0 ]; then
   sed "$linenum3,$linenum4 d" < tempfile3 > tempfile4
   chmod 777 tempfile4
   awk '{print $2, $13}' tempfile4 > parsed
   rm elim2.sed
   rm tempfile*
else
   awk '{print $2, $13}' tempfile3 > parsed
   rm elim2.sed
   rm tempfile*
fi

popd

