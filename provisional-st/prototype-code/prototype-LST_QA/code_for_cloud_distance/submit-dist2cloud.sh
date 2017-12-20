#/bin/bash

# Another constant variable used to name the slurm submission file that
#   this script is going to submit to slurm.
  jobfile="dist2cloud_global.bash"


# decide the following paths/names:
# 1) folder where Landsat 7 MTL's are
# 2) folder where MODIS will be downloaded
# 3) name of file that says which MODIS scene goes with which L7 scene
  path_row="43_36"
  year="2009"
  main_folderpath='/home/kga1099/global_clouds/'
  cloudmask_path="/home/kga1099/global_clouds/cloud_files/${path_row}/${year}/"
  lst_sst_file="/home/kga1099/global_clouds/input_files/${path_row}/last_study_${path_row}_${year}.txt"
  results_filename="full_last_study_${path_row}_${year}.txt"

  echo "Landsat Scene, No., Easting, Northing, Lsat_lst, Modis_sst, AvgQual, AreaStdev, NumZeros, Lobs, tau, up, down, Dist2cloud" >> "${results_filename}"
  
  for cloud_file in "${cloudmask_path}"/*_cfmask.tif
    do
    cloudmask_name=$(basename ${cloud_file})
    currentscene=${cloudmask_name:0:21}

    # Give our job a meaningful name
    jobname="${currentscene}"
    echo "Submitting job $jobname"
 
    # Setup where we want the output from each job to go
    outfile="./logfile_${jobname}.txt"
    
    # load idl module
    module load envi

    # Submit the job
    sbatch --qos free --partition work -x woodcrest-[01-24]  --mem=2000 -J $jobname -o ${outfile} ${jobfile} ${main_folderpath} ${cloudmask_path} ${cloudmask_name} ${lst_sst_file} ${results_filename}
  done
