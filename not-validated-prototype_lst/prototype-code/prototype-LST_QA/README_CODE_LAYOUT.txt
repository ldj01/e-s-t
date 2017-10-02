LST Code Layout, June 2017
Dr. Kelly G. Laraby
kga1099@rit.edu


This file's purpose is to present an overview of the code layout and where to
find various pieces of code. More specific README files are located within some
of the subfolders. Most if not all IDL scripts have significant comments
throughout them as well.

The Landsat Surface Temperature code is divided into three folders for easier
navigation, with the following brief descriptions:

  code_for_generating_product
      All the programs and scripts needed to generate a 5 band geotiff containing
      atmopsheric parameters, an LST geotiff, and an option to generate a LST
      uncertainty geotiff. Currently, a fake emissivity image is generated that 
      has every pixel set to the emissivity of water. Additionally, the LST
      uncertainty code assumes that a distance to nearest cloud image has been 
      generated for a given Landsat scene.

  code_for_cloud_distance
      Contains code to generate a distance to cloud geotiff for a given Landsat
      cloud mask.

  code_for_validation
      Contains code that was used to perform validation studies (some are for
      use with buoys and some are for use with MODIS SST images).



Misc Notes:

  1) Code that generates 5 band geotiff has been run on Landsat 5, 7, and 8. The
     part that generates the LST band should also work for these instruments
     with no issue. For Landsat 8 the product has to be generated individually
     for each of the two thermal bands.

  2) The code that estimates LST uncertainty is currently only in place for the 
     scripts that use MERRA. A small study must be done for NARR to get some
     coefficients, but the uncertainty code will remain largely the same.

  3) Any file that starts with "submit-" and ends in ".sh" is a script that was
     used to submit jobs to a scheduler on an RIT research server. By removing
     the command from these files that start with "sbatch," you can instead
     kick off the job right from the within the submit script.

  4) The validation code that compares to buoys has been used for Landsat 5 & 7,
     but the GLOBAL validation code can only be used for Landsat 7 because it 
     uses Modis Sea Surface Temperature (SST) as truth.

  5) Again, make sure to read the other README files, especially README_MAIN.txt
     because it gives great detail on all the scripts involved, how to run the
     provided example, and what lines you'll need to change in order for the
     code to work properly or at all.
