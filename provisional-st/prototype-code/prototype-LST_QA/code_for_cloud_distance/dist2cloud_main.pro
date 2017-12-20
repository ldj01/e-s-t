

PRO dist2cloud_main, cloudmask_path, cloudmask_name, easting, northing, results_file

imagebase = STRMID(cloudmask_name, 0, 21)

cfmask = READ_TIFF(cloudmask_path+cloudmask_name, GEOTIFF=cloudmask_geotiff)

cloudmask = GET_CLOUDMASK(cloudmask_path, imagebase)

distance = CALC_NEAREST_CLOUD_GLOBAL(cloudmask, cloudmask_geotiff, easting, northing)

OPENW, lun, results_file, /GET_LUN
  PRINTF, lun,  distance
CLOSE, lun
FREE_LUN, lun

END 
