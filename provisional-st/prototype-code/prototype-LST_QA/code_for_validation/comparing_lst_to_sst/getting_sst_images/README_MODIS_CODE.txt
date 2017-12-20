                               Collection of MODIS-related scripts

                                         Kelly G. Laraby
                                PhD Candidate for Imaging Science
                                Rochester Institute of Technology
                                         kga1099@rit.edu




CONTENTS OF THIS FILE
=====================

   * Introduction & Background

   * Required Files

   * How to Run




INTRODUCTION & BACKGROUND
=========================

  The motivation behind this collection of scripts has to do with the validation of the Landsat
Land Surface Temperature (LST) product. There is an algorithm that predicts surface temperature from
Landsat thermal imagery, but it is important to be able to compare it to "truth data" so the 
accuracy of the algorithm can be quantified. Studies have shown (not presented here) that the MODIS 
Sea Surface Temperature (SST) product is an adequate source of ground truth as long as the best
quality pixels are used. This will only be useful for validating LST predictions using Landsat 7 
imagery because Landsat 7 and MODIS (aboard the TERRA satellite) have a very similar orbit, which
means that images of the same spatial area are captured around the same time. MODIS SST files can be 
downloaded from https://oceancolor.gsfc.nasa.gov either by browsing manually or by using direct data 
access. Unfortunately, it is not obvious by the name of the SST files which one overlaps the Landsat 
scene of interest, making it difficult to automate the downloading process. This is because there 
are SST images captured every 5 minutes, and the one that encapsulates the Landsat image of interest
can have an image acquisition time somewhere from 10-30 minutes after the Landsat acquisition time.
Since it was necessary to compare many hundreds of Landsat LST images to MODIS SST images, effort
was placed towards resolving this issue, which saved a lot of time in the long run.

The solution that was implemented consisted of using the Landsat image capture time to guess the
MODIS acquisition time (the time must be known because it is part of the file name used to download
the SST images). "Guessing" the MODIS acquisition time refers to the fact that the SST image that
encapsulates the Landsat scene may have been captured 10 minutes after Landsat, or maybe 15-30 min
after Landsat. The first "guess," then, would be that the time difference was 10 minutes, and that
SST file gets downloaded. An IDL program checks if the Landsat scene falls within the SST image that
was downloaded, and if it does not the file is deleted and a new SST image that was captured five 
minutes later is downloaded. Up to 10 attempts are allowed to find the correct image, and if it is 
found before then the program stops and keeps the correct SST file. The name of the SST file is 
recorded along with the  Landsat scene name and Landsat's corner coordinate information.

After the appropriate SST files have been downloaded, they must be georeferenced before they can be
compared to Land Surface Temperature predictions. This can easily done using ENVI Classic, but it
would be very tedious to repeat this process hundreds of times. Therefore, an IDL program was
created to be able to georeference a given list of SST files. This program also subsets the image
because MODIS imagery covers a vast spatial area, where the Landsat scene only occupies about a
tenth of the image. Subsetting the image to this smaller area reduces the size of the file that is
saved at the end of the program, which aids computation time for other scripts that perform the
actual comparision between LST and SST images. The comparison code is not included here, because
it would require the user to have all the LST imagery. Instead, the user is equipped to be able to
run the download process and the georeferencing process for a few example scenes. 
 


REQUIRED FILES
==============

   There are many files that are needed in order for LST to run smoothly (or at all). Currently, the
LST code is set up to look for all the files it needs in the same directory that the main code is
initialized. Below is a list of all the files needed and a brief description of what they are:

 DOWNLOADING MODIS SST FILES
   submit-modisjob.sh             ........ submits job for every MTL file in target folder
   get_modis_for_L7scenes.bash    ........ main program for downloading MODIS SST files
   is_right_modis_scene.pro       ........ IDL program that checks if SST file is correct
   Landsat MTL.txt files          ........ metadata files that contain info for each Landsat scene

 
 GEOREFERENCING MODIS SST FILES
   change_modis_windows.pro       ........ georeferences a list of SST files (Windows OS)
   convert_ll_utm.pro             ........ converts lat/lon coordinates to UTM coordinates
   input text file                ........ text file that was produced by the SST download process



   
HOW TO RUN
==========

 DOWNLOADING MODIS SST FILES
   
   1) This process is designed to be run on a Unix OS. If using a parallel processing system, 
      open submit-modisjob.sh and change paths (path to Landsat MTL files, path to where 
      downloads go, path and name of output text file). If not using a parallel processing 
      environment, open modisjob-batch.sh, and change paths just mentioned. This will take a 
      little longer, but only two examples are included.

   2) Either run submit-modisjob.sh or modisjob-batch.sh, and check that a text file appears
      that contains the Landsat scnene name, MODIS SST scene name, UTM zone, and Landsat upper
      left and lower right lat/lons. In other words, there should be 7 columns, and a row of
      information for each case that was run (and also a header line). At this point, the
      download process is done and the downloaded SST files can be used to perform the 
      georeferencing process, which is described in the next section. 


 GEOREFERENCING MODIS SST FILES

   1) This process must be run on Windows OS, and you must have both IDL and ENVI classic open.
      The program change_modis_windows.pro must be given the path to where the original SST files
      are located, the path to where the text file produced from the downloading process is kept,
      and the name of that text file. Make sure the function convert_ll_utm.pro is in the working
      directory.

   2) The georeferenced images will automatically be saved to folder called "georeferenced" that
      will be in the folder that was specified by the first arguement.
 

