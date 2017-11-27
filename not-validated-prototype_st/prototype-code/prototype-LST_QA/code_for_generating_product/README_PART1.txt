                    Landsat Surface Temperature Generation Code, as of 05/20/2017

                                         Kelly G. Laraby
                                PhD Candidate for Imaging Science
                                Rochester Institute of Technology
                                         kga1099@rit.edu




CONTENTS OF THIS FILE
=====================

   * Introduction

   * Required Files

   * Running LST

   * Changes You May Need to Make

   * Troubleshooting






INTRODUCTION
============

   Welcome to the Landsat Surface Temperature (LST) tool! The purpose of this project is to
be able to accurately report the surface temperature at every pixel within a Landsat scene. This
requires knowledge of the atmosphere of the scene at the time of acquisition, which can be
gained from reanalysis data (relies on radiosonde measurements). The atmospheric data are then
used to perform atmospheric compensation via MODTRAN, and after implementing some fancy tricks we 
can calculate the transmission, upwelled radiance, and downwelled radiance for each reanalysis 
point. Then, simple interpolation is performed to obtain these three parameters on a per-pixel 
level. The final output is an array containing the original Landsat thermal band, the elevation, and
the three parameters of interest for every pixel in the scene. See part 2 to explore the ways you
can use this array of parameters to obtain information such as surface temperature at specific 
points or determine distance to the nearest cloud!

   Here are some specifics on the capabilities of the current LST process:

      * Can be used for Landsat 5, 7, or 8 (bands 10 and 11 for Landsat 8) 
      * Can use either NARR reanalysis for North America or MERRA for global scenes
      * Can be run in "zeroDEM" mode, which lets the user specify a single elevation for the whole scene
      * Uses MODTRAN 4v3 for atmospheric compensation




REQUIRED FILES
==============

   There are many files that are needed in order for LST to run smoothly (or at all). Currently, the
LST code is set up to look for all the files it needs in the same directory that the main code is
initialized. Below is a list of all the files needed and a brief description of what they are:

 INITIAL FILES
   submit-narr.sh                ........ submits LST jobs to a scheduler, calls LST code for NARR
   submit-merra.sh               ........ submits LST jobs to a scheduler, calls LST code for MERRA
   submit-narr-zerodem.sh        ........ submits NARR zeroDEM job, where user specifies elevation
   submit-merra-zerodem.sh       ........ submits MERRA zeroDEM job, where user specifies elevation
   LSTprocess_narr.bash          ........ main program. downloads NARR and calls IDL programs
   LSTprocess_merra.bash         ........ main program. downloads MERRA and calls IDL programs
   LSTprocess_narr_zeroDEM.bash  ........ main program for NARR with user specified elevation
   LSTprocess_merra_zeroDEM.bash ........ main program for MERRA with user specified elevation

 
 IDL PROGRAMS  
   lst_narr_step1.pro          ............ step 1 in LST process using NARR
   lst_narr_step2.pro          ............ step 2 in LST process using NARR
   lst_narr_step3.pro          ............ step 3 in LST process using NARR
   lst_merra_step1.pro         ............ step 1 in LST process using MERRA
   lst_merra_step2.pro         ............ step 2 in LST process using MERRA
   lst_merra_step3.pro         ............ step 3 in LST process using MERRA
   lst_narr_step3_zerodem.pro  ............ step 3 in LST process using NARR and single elevation
   lst_merra_step3_zerodem.pro ............ step 3 in LST process using MERRA and single elevation
   calculate_lobs.pro          ............ calculates observed radiance
   calculate_lt.pro            ............ calculates blackbody radiance from temperature
   convert_geopot2geomet.pro   ............ converts geopotential height to geometric height
   convert_ll_utm.pro          ............ converts lat/lon coordinates to UTM
   convert_rad_temp.pro        ............ converts radiance to temperature using lookup table
   convert_sh_rh.pro           ............ converts specific humity to relative humidity
   distance_in_utm.pro         ............ calculates distance between two UTM coordinates
   interp_to_height.pro        ............ interpolates parameters to the elevation of a pixel 
   interp_to_location.pro      ............ interpolates parameters to spatial location of a pixel
   planck_eq.pro               ............ calculates radiance using wavelength and temperature

 OTHER SCRIPTS AND FILES
   DEM                    .............. folder containing DEM files for various path/rows
   GRIB                   .............. folder containing GRIB manipulation scripts
   LSTscenes2run          .............. folder containing scenes to run LST on
   download_merra.bash    .............. script that downloads MERRA data
   elim2.sed              .............. deletes headers in tape6 files
   coordinates.txt        .............. lat/lon grid information for the NARR GRIB data
   head.txt               .............. file needed to make tape5 files
   tail.txt               .............. file needed to make tape5 files
   L5.rsp                 .............. spectral response data for Landsat 5
   L7.rsp                 .............. spectral response data for Landsat 7
   L8_B10.rsp             .............. spectral response data for Landsat 8, band 10
   L8_B11.rsp             .............. spectral response data for Landsat 8, band 11
   script_HGT_generic     .............. gets geopotential heights from NARR GRIB data 
   script_SHUM_generic    .............. gets specific humidity from NARR GRIB data
   script_TMP_generic     .............. gets temperatures from NARR GRIB data
   stanAtm.txt            .............. standard atmosphere for MODTRAN
   tape6parser.bash       .............. parses tape6 files, the output MODTRAN files





   
RUNNING LST
===========

   Much of the LST process is automatic, but there are a few things that need to happen first. Here
is a general guideline of how to ensure a smooth run:

   1) Put the thermal band (if Landsat 7 use B6_VCID_1) and MTL file for every scene you want to
      run in a folder. Given folder is called 'LSTscenes2run'.

   2) Make sure you have the proper DEM files in the DEM folder for every path/row you plan to use,
      or if you want to use a single elevation for the whole scene (e.g. if you only care about 
      water pixels) you must use the "zeroDEM" version of the code and specify an elevation [km].

   3) If you have a parallel computing machine that can have jobs "submitted" to it, then continue
      on to step 4. If not, or not sure, you will have to write your own bash script that runs the
      LST code for every scene in the folder you made in step 1, and you should skip to step 7.     
   
   4) If you want to run LST using NARR, open submit-narr.sh (or submit-narr-zerodem.sh)
      If you want to run LST using MERRA, open submit-merra.sh (or submit-merra-zerodem.sh)

   5) Change the foldername in ' for pathrow in "foldername"; do ' to the folder where the scenes
      you want to run are. Then, if necessary, change the way the script submits the job, because
      it is most likely different for other scheduler machines.

   6) Run the job submission file (e.g. './submit-narr.sh') to submit the jobs to the scheduler  

   7) Once all the jobs are complete, there should be a folder for every scene you ran. Inside each
      folder should be a file ending in either "LSTparams_NARR.tif" of "LSTparams_MERRA.tif", which
      is the array that contains the parameters for every pixel in the Landsat scene. You will need 
      this file for part 2 of the LST process, so that you can estimate the surface temperature at 
      a certain point and also calculate distance to the nearest cloud.
      



CHANGES YOU MAY NEED TO MAKE
============================

   There are a few instances where the LST uses a hard coded path to find/use something. Here is a
list of such documented cases:

   In submit-narr.sh
      submit-merra.sh

        line 20: if Landsat 8 is being used, you must specify which band you want to use the LST
                 process on.


   In submit-narr-zerodem.sh
      submit-merra-zerodem.sh
        
        line 8: Define elevation [km] that you want to be used for the whole scene.

        line 23: if Landsat 8 is being used, you must specify which band you want to use the LST
                 process on.

        
   In lst_narr_step1.pro

        line 542: command that links DATA and runs MODTRAN is recorded. May need to change path.

   In lst_merra_step1.pro

        line 711: command that links DATA and runs MODTRAN is recorded. May need to change path.




TROUBLESHOOTING
===============

   Unfortunately, there is no thorough troubleshooting guide for the LST process yet. The closest
thing is the terminal.txt file that will show up in every results folder, and it may indicate where
something went wrong. The biggest sign that something went wrong is if no LSTparams.tif file shows
up when the LST process is complete. The first thing to do is to look at the terminal.txt file, and
answer the following questions:

   1) Was the reanalysis data successfully downloaded?
        Both NARR and MERRA may have limitations on how many downloads can be made at once (this
        limit is unknown for NARR but the limit is 10 for MERRA). If the terminal.txt file shows
        an error that says "login incorrect" when trying to download MERRA data, it is because
        more than 10 instances of downloading data exist (so in this case just try again later).
        Sometimes the terminal.txt file will say the MERRA file does not exist, so in that case
        you're out of luck. 

   2) Is the path to MODTRAN correct?
        If the program can't find the path to MODTRAN, the terminal.txt file  should say something
        to that effect.

   3) Are all the LST files in the right spot?
        The terminal.txt file should say something whenever it can't find a file. Also make sure the
        folder and file names that are called by the code haven't been changed or that they match.

   4) Are the IDL programs that are called in the LSTprocess code all lowercase?
        If the program names haven't been changed this should not be an issue, but know that they
        need to be lowercase because Unix is case sensitive and IDL is not.     

   If the problem is not easily found, you may want to look at the temporary files that are normally
deleted at the end of the LST process. Just uncomment the lines at the end of either LSTprocess_narr
or LSTprocess_merra that delete the intermediate files. That way, you can make sure the tape5 files
were generated, that the tape6 files are present and not full of 0's, and that the commandlist has
the correct path for MODTRAN.



Good Luck, and enjoy!
