#!/bin/bash

#  change to directory containing MODTRAN output of current case
pushd $1

#  determine line where pertinent output ends and delete everything below that
linenum1=`grep -n "MULTIPLE SCATTERING " tape6 |head -1| awk '{print $1}'| awk -F':' '{print $1}'`
sed "$linenum1,$ d" < tape6 > tempfile

echo linenum1 = $linenum1

chmod 777 tempfile

#  determine line where pertinent output begins and delete everything before that
linenum2=`grep -n "WAVELENGTH" tempfile |head -1| awk '{print $1}' |awk -F':' '{print $1}'`
sed "1,$linenum2 d" < tempfile > tempfile2

echo linenum2 = $linenum2

chmod 777 tempfile2

#  delete headers that divide data within file so data is continuous columns
sed -f elim2.sed < tempfile2 > tempfile3

chmod 777 tempfile3

#  determine if tape6 contains warning
linenum3=`grep -n "WARNING" tempfile3 | head -1 | awk '{print $1}'|awk -F':' '{print $1}'`
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

