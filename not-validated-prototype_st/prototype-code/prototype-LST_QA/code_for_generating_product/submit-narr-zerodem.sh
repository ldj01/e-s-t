#/bin/bash

# Another constant variable used to name the slurm submission file that
#   this script is going to submit to slurm.
jobfile="LSTprocess_narr_zeroDEM.bash"

# Define elevation [km] to be used for the whole scene
height=0.0

for pathrow in "LSTscenes2run"; do

	for file in "$pathrow/"*_MTL.txt; do

                # Extract basename
                base=`basename $file .txt | awk -F'_' '{print $1}'`
                
                # Make directory to contain results
                mkdir $base

                which=${base:2:1}
		if [ ${which} -eq 8 ]
		then
			which=11
		fi

		# Give our job a meaningful name
		jobname=$base
		echo "Submitting job $jobname"

		# Setup where we want the output from each job to go
		outfile=$base/terminal.txt

		case ${which} in
			5)  cp $pathrow/$base'_B6.TIF' $base
                	    cp $pathrow/$base'_MTL.txt' $base
			    ;;
			7)  cp $pathrow/$base'_B6_VCID_1.TIF' $base
			    cp $pathrow/$base'_MTL.txt' $base
			    ;;
			10)  cp $pathrow/$base'_B10.TIF' $base
			     cp $pathrow/$base'_MTL.txt' $base
			     ;;
			11)  cp $pathrow/$base'_B11.TIF' $base
			     cp $pathrow/$base'_MTL.txt' $base
			     ;;
		esac
		
		# Submit the job
		sbatch --qos free --partition work --mem=6000 -J $jobname -o $outfile $jobfile $base $base $height
	done;
done
