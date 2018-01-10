README for cloud distance code

Kelly G. Laraby
PhD Candidate of Imaging Science
Rochester Institute of Technology
kga1099@rit.edu


This README file discusses the code that deals with "distance to nearest cloud."
Distance to nearest cloud is used to estimate the uncertainty in the Landsat
Surface Temperature (LST) for every pixel. In order to do this, a distance to
cloud image must be generated at Landsat resolution. There is a script to
generate these images for any given Landsat cloud mask, but there are also
scripts that were used for the global validation of Landsat 7 that calculates 
distance to nearest cloud for individual locations within a Landsat scene.
These latter files can be found with the validation code. 


GENERATING DISTANCE TO CLOUD IMAGES

Currently, the code that generates the distance to nearest cloud images is
separate from the main LST code. If the option is set in the LST code to
estimate uncertainty, then it expects to find all the necessary distance image
files in a specified folder.
  
The file generate_dist2cloud_image.pro is the script that accepts the path to
the Landsat cloud mask file, the output path to write the result to, and an
indication of which Landsat instrument is being used (SLC gap is dealt with if
it is Landsat 7). It has only been used to test individual examples, so there is
no script to generate images for many Landsat scenes. There is a Landsat cloud
mask included in this folder.


