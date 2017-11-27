#/bin/bash
#
# Kelly G. Laraby
# PhD Candidate in Imaging Science
# Rochester Institute of Technology
# kga1099@rit.edu
#
#
# submit-modisjob.sh
# ------------------
# Job submission file to run several instances of downloading MODIS Sea Surface Temperature (SST) files 
# in parallel. 
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
 



# Name the slurm submission file that this script is going to submit to slurm.
  jobfile="get_modis_for_L7scenes.bash"


# Decide the following paths/names:
# 1) folder where Landsat 7 MTL's are
# 2) folder where MODIS will be downloaded
# 3) name of file that will be written that says which MODIS scene goes with which L7 scene
  path_row="232_17"
  year="extra"
  landsat_mtl_path="/home/kga1099/MODIS/compare_to_lsat7/landsat/mtls/L7global/${path_row}/${year}/"
  download_path="/home/kga1099/MODIS/compare_to_lsat7/modis/hdfs/${path_row}"
  results_filename="modis_scenenames_${path_row}_${year}.txt"

# Write header line to results file
  echo "Landsat Scene             MODIS Scene        Zone   UL_lat      UL_lon       LR_lat      LR_lon" >> "${results_filename}"
  

# For each MTL file in target folder, submit a job that will download the correct MODIS SST file 
# for each L7 scene, then add the scene names and Landsat corner coordinates to the results file
  for mtl_file in "${landsat_mtl_path}"/*_MTL.txt
    do
    base=$(basename $mtl_file)
    currentscene=${base:0:21}

    # Give our job a meaningful name
    jobname="${currentscene}"
    echo "Submitting job $jobname"
 
    # Setup where we want the output from each job to go
    outfile="./logfile_${jobname}.txt"
    
    # load idl module
    module load envi

    # Submit the job
    sbatch --qos free --partition work -x woodcrest-[01-24]  --mem=2000 -J $jobname -o ${outfile} ${jobfile} ${landsat_mtl_path} ${download_path} ${currentscene} ${results_filename}
  done
