#/bin/bash

# Another constant variable used to name the slurm submission file that
#   this script is going to submit to slurm.
  jobfile="dist2cloud_perpixel.bash"


# decide the following paths/names:
# 1) folder where Landsat 7 MTL's are
# 2) folder where MODIS will be downloaded
# 3) name of file that says which MODIS scene goes with which L7 scene
  cloudmask_folderpath="/home/kga1099/global_clouds/cloud_files/43_36/2015/"
  imagebase='LE70430362015241EDC00'
  output_filepath="./${imagebase}_cloud_distance_image.tif"

  
  cloudmask_full_path="${cloudmask_folderpath}${imagebase}_cfmask.tif"

  # Give our job a meaningful name
  jobname="${imagebase}"
  echo "Submitting job $jobname"
 
  # Setup where we want the output from each job to go
  outfile="./logfile_${jobname}.txt"
    
  # load idl module
  module load envi

  # Submit the job
  sbatch --qos free --partition work -x woodcrest-[01-24]  --mem=2000 -J $jobname -o ${outfile} ${jobfile} ${cloudmask_full_path} ${output_filepath}
