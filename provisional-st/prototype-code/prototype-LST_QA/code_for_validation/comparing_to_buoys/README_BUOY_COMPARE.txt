Code for Comparing LST to Buoys
Dr. Kelly Laraby
06/20/2017

This file explains how to run the code that records LST values and buoy-measured
surface temperatures, as well as other useful info. This code does not have a 
bash script that calls the IDL code, so either one can be made or it can be 
initiated from the IDL command line.

Below are a list of the scripts used for performing buoy validation:

  atmparams_buoys.pro
      Main script, takes path to file that contains Landsat scene IDs, a path to
      the output file, and an option to use provided lat lons. If this option is
      not set, the script will search through a buoy history file to get the 
      correct lat lon. Otherwise it gets the lat lons from a provided file. Note
      that when atmospheric parameters are retrieved (see next function), the 
      option to calculate distance to cloud is set. So if distance to cloud info
      is not desired, the option needs to be removed manually.

  get_atm_params.pro
      Function called by the main script, reads in 5 band geotiff file that is 
      created by the main LST code, and calculates atmospheric parameters and
      LST at a specific lat lon, which in this case is a buoy's lat lon. This
      function also has an option to calculate distance to nearest cloud and
      return the result with the other parameters. This option utilizes the
      two scripts below.
 
  get_cloud_mask.pro
      Function that reads in Landsat cloud mask, and alters values so that
      there is a binary image indicating cloud/no cloud. It returns the 
      binary image.

  calc_nearest_cloud_from_loc.pro
      Function that takes a binary cloud mask & geotiff, and an easting and 
      northing to calculate distance to nearest cloud at that location. The
      script calculates the distance from the location of interest to every
      cloud pixel in the cloud mask image, then returns the smallest value
      in units of kilometers.

  convert_ll_utm.pro
      Converts lat/lon to UTM coordinates.

  convert_rad_temp.pro
      Converts radiance(s) to temperature(s) using a lookup table.
 
  whichBuoyInScene.txt
      Contains paths, rows, and buoy numbers that fall within them.

  buoyHistory.txt
      Contains history info for many N.A. buoys.

  buoylatlon.txt
      Contains paths, rows, buoy numbers, and specific lat/lons for them. Useful
      for buoys that are not in the history file.
 
  LUT5.txt
      Lookup table to use for calculating LST.

  scenelist.txt
      Text file that contains Landsat IDs to get validation results for.
