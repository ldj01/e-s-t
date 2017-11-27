; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    CALC_NEAREST_CLOUD_FROM_LOC
;
; PURPOSE:
;    IDL FUNCTION
;    Calculates the distance from a specified location to each cloud pixel in order to determine the distance to the nearest cloud pixel.
;    
; OUTLINE:
;    1) retreive upper left coordinates and zone for cloudmask
;    2) convert buoy lat/long to UTM
;    3) calculate difference between Landsat UTM and buoy UTM and divide by pixel size, which allows array indexing
;    4) index the cloudmask where there are clouds (value of 255), and obtain cols/rows
;    5) calculate distance from buoy col/row to each cloud pixel col/row
;    6) save distances in an array and output to tiff file
;    7) return smallest distance
;
; CALL SEQUENCE: 
;    CALC_NEAREST_CLOUD_FROM_LOC,   cloudmask, $               ; cloud mask to be altered
;                                 cmask_geotiff, $           ; cloud mask geotiff 
;                                 easting, $                 ; easting for location to do calculations
;                                 northing, $                ; northing for location to do calculations 
;
; RESTRICTIONS:
;   N/A
;
;
; REQUIRED PROGRAMS (in working directory):      
;   N/A 
;
; REQUIRED FILES (in folderpath directory):      
;   N/A
;
; MODIFICATIONS:
;    January, 2015        Original code
;
;

FUNCTION CALC_NEAREST_CLOUD_FROM_LOC, cloudmask, cmask_geotiff, easting, northing
  
  cmask_ULeast= cmask_geotiff.modeltiepointtag[3]
  cmask_ULnorth= cmask_geotiff.modeltiepointtag[4]
  zone = STRMID(STRCOMPRESS(STRING(cmask_geotiff.PROJECTEDCSTYPEGEOKEY)),4,2)
  
  pixsizeX = cmask_geotiff.MODELPIXELSCALETAG[0]
  pixsizeY = cmask_geotiff.MODELPIXELSCALETAG[1]
  
  IF cmask_ULnorth LT 0 THEN BEGIN
    cmask_ULnorth = cmask_ULnorth + 10000000
  ENDIF 
  
  ;
  ; calculate difference between Landsat UTM and buoy UTM
  ; for easting, do buoy-Landsat (b/c buoy has to be to the right of the Landsat easting)
  ; for northing, do landsat-buoy (b/c buoy has to be below the Landsat northing)
  ; put in terms of pixels for indexing
  ;
  utm_col= DOUBLE(ROUND((easting - cmask_ULeast)/pixsizeX))
  utm_row= DOUBLE(ROUND((cmask_ULnorth - northing)/pixsizeY))
  
  ;
  ; find where cloud pixels are
  ;
  s = SIZE(cloudmask)
  ncol = s[1]
  cloudind = WHERE(cloudmask EQ 255)
  cloudCols = cloudind MOD ncol
  cloudRows = cloudind / ncol
  
  distances= DBLARR(1,N_ELEMENTS(cloudCols))
  
  ;
  ; calculate distance to each cloud pixel and build up distance arrays
  ;
  PRINT, 'Calculating distance to all cloud pixels...'
  FOR i= 0, N_ELEMENTS(cloudCols)-1 DO BEGIN
    
    col_diff= ABS(utm_col - cloudCols[i])
    row_diff= ABS(utm_row - cloudRows[i])
    distances[i]= SQRT( col_diff^2 + row_diff^2 ) 
    
    dist_in_km=distances[i]*pixsizeX/1000
      
  ENDFOR
  
  ;
  ; find smallest distance in kilometers
  ;
  minIndex= WHERE( distances EQ MIN(distances) )
  minDistance= distances[minIndex[0]]*pixsizeX/1000
  RETURN, minDistance
  
END


