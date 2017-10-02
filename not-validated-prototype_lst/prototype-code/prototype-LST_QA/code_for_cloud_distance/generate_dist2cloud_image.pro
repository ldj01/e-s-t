; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    GENERATE_DIST2CLOUD_IMAGE
;
; PURPOSE:
;    IDL PROCEDURE
;    This script generates a "distance to nearest cloud" image from a Landsat cloud mask image. This utilizes idl's
;    MORPH_DISTANCE function, and multiplies the result by the Landsat pixel size so that the distance image is in km.
;    
; OUTLINE:
;    1) read in cloud mask image
;    2) set everything that is not labeled as a cloud to 0
;    3) use MORPH_DISTANCE to calculate distances from non-zero pixels (which are now just clouds)
;    4) multiply result by 0.03 to get distances in km (pixel size is 0.03km)
;    5) write dist2cloud image to a geotiff
;
; CALL SEQUENCE: 
;    GENERATE_DIST2CLOUD_IMAGE,  cloudmask_filepath, $             ; full path to cloud mask file
;                                output_filepath, $                ; full path to output file
;                                which_landsat, $                  ; integer indicating which landsat imstrument us used
;
;
;
; REQUIRED PROGRAMS (in working directory):      
;    N/A
;
; REQUIRED FILES (in folderpath directory):      
;    N/A
;
; MODIFICATIONS:
;    May, 2017        Original code
;
;

PRO GENERATE_DIST2CLOUD_IMAGE, cloudmask_filepath, output_filepath, which_landsat
  
  ;
  ; read in cloud mask
  ;
  cmask = READ_TIFF(cloudmask_filepath, GEOTIFF = cmask_geotiff)
  
  ;
  ; set everything but clouds to zero
  ;
  cmask_binary = cmask
  cmask_binary[WHERE(cmask NE 4)] = 0
  
  ;
  ; use morph_distance to calculate distance to nearest cloud
  ;
  dist2cloud_image = MORPH_DISTANCE(cmask_binary, /BACKGROUND, NEIGHBOR_SAMPLING=3)
  dist2cloud_image = dist2cloud_image * 0.03 ; multiply by pixel size in km
 
  ; 
  ; if using Landsat 7, set gap locations to -1
  ;
  IF which_landsat EQ 7 THEN BEGIN
    ; get SLC gap locations
    gap_locations = WHERE(cmask EQ 255,/NULL)
    dist2cloud_image[gap_locations] = -1.0
  ENDIF

  ; 
  ; write image to tiff file
  ;
  WRITE_TIFF, output_filepath, dist2cloud_image, GEOTIFF=cmask_geotiff, /FLOAT


  END


