#!/bin/bash

#  change to directory containing MODTRAN output of current case
pushd $1

#  determine line where pertinent output ends and delete everything below that
linenum1=`grep -n "1 MULTIPLE" tape6 |head -1| awk '{print $1}'| awk -F':' '{print $1}'`
sed "$linenum1,$ d" < tape6 > /dev/shm/{$SLURM_JOB_ID}_tempfile

chmod 777 /dev/shm/{$SLURM_JOB_ID}_tempfile

#  determine line where pertinent output begins and delete everything before that
linenum2=`grep -n "0 FREQ" /dev/shm/{$SLURM_JOB_ID}_tempfile |head -1| awk '{print $1}' |awk -F':' '{print $1}'`
sed "1,$linenum2 d" < /dev/shm/{$SLURM_JOB_ID}_tempfile > /dev/shm/{$SLURM_JOB_ID}_tempfile2

chmod 777 /dev/shm/{$SLURM_JOB_ID}_tempfile2

#  delete headers that divide data within file so data is continuous columns
sed -f elim2.sed < /dev/shm/{$SLURM_JOB_ID}_tempfile2 > /dev/shm/{$SLURM_JOB_ID}_tempfile3

chmod 777 /dev/shm/{$SLURM_JOB_ID}_tempfile3

#  determine if tape6 contains warning
linenum3=`grep -n "WARNING" /dev/shm/{$SLURM_JOB_ID}_tempfile3 | head -1 | awk '{print $1}'|awk -F':' '{print $1}'`
linenum4=`grep -n "WARNING" /dev/shm/{$SLURM_JOB_ID}_tempfile3 | tail -1 | awk '{print $1}'|awk -F':' '{print $1}'`

#  if tape6 contains warning, parse warning and write to tempfile4
#  then print (wavelength, total radiance) to parsed file
#  else print (wavelength, total radiance) to parsed file

if [ $linenum3 > 0 ]; then
   sed "$linenum3,$linenum4 d" < /dev/shm/{$SLURM_JOB_ID}_tempfile3 > /dev/shm/{$SLURM_JOB_ID}_tempfile4
   chmod 777 /dev/shm/{$SLURM_JOB_ID}_tempfile4
   awk '{print $2, $13}' /dev/shm/{$SLURM_JOB_ID}_tempfile4 > parsed
   rm elim2.sed
   rm /dev/shm/{$SLURM_JOB_ID}_*
else
   awk '{print $2, $13}' /dev/shm/{$SLURM_JOB_ID}_tempfile3 > parsed
   rm elim2.sed
   rm /dev/shm/{$SLURM_JOB_ID}_*
fi

popd

