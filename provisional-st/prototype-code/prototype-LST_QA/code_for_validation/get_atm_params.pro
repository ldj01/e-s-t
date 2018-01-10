; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
;
; NAME:  
;    GET_ATM_PARAMS
;
; PURPOSE:  
;    IDL FUNCTION
;    This code calculates and returns all the atmospheric parameters above the buoy in a particular Landsat scene. There is also
;    the option to calculate the distance to nearest cloud.
;
; OUTLINE:
;    1) read in LST results and get geotiff info
;    2) calculate difference between Landsat UTM and buoy UTM
;    3) select radiance values at and around buoy location, and change size
;       if there are many "zero" pixels due to the broken SLC in Landsat 7.
;    4) if cloudmask keyword is set, get the distance to nearest cloud
;    5) take the average of each atmospheric parameter for the pixels around the buoy (not including gap pixels)
;    6) pick emissivity of water and lookup table based on which Landsat is being used
;    7) use the big equation to get observed radiance, and use LUT to get temperature
;    8) return array of parameters
;
;
; CALL SEQUENCE: 
;    GET_ATM_PARAMS,    folderpath, $                   ; path where results_float.tif and LUTs are located (include '\' after folder)
;                       imagebase, $                    ; scene identifier (e.g. LE70160302008113EDC00)
;                       buoylat, $                      ; latitute of buoy
;                       buoylon, $                      ; longitude of buoy (negative if in western hemisphere)
;                       whichLandsat, $                 ; which Landsat is being used (5,7,10,11, where 10&11 are bands on L8)
;                       GETDIST2CLOUD=getdist2cloud, $  ; if set, checks if clouds are over buoy                               
;
;
; RESTRICTIONS:
;    1) LUT is created with 500 entries
;    2) when checking for clouds, more than 50 "cloud pixels" in the 10x10 surrounding the buoy means the scene should be skipped.
;       For now, this number is arbitrary and may need to be adjusted.
;
; REQUIRED PROGRAMS (in working directory):   
;    CONVERT_LL_UTM.PRO
;    CONVERT_RAD_TEMP.PRO
;    GET_CLOUD_MASK.PRO (optional)
;    CALC_NEAREST_CLOUD.PRO (optional)
;                                            
;  
; REQUIRED FILES (in folderpath directory):   
;    imagebase_LSTresults.tif (for every scene)
;    LUT7.txt (or other LUT)
;    cloud mask files (optional. see get_cloudmask.pro)
;
; MODIFICATIONS:
;    November, 2015       Original code
;
;
FUNCTION GET_ATM_PARAMS, folderpath, imagebase, buoylat, buoylon, whichLandsat, GETDIST2CLOUD=getdist2cloud

  ;
  ; initialize parameters array
  ; 
  IF KEYWORD_SET(getdist2cloud) THEN BEGIN
     params=FLTARR(7,1)
  ENDIF ELSE BEGIN
     params=FLTARR(6,1)
  ENDELSE
  ;
  ; read in LST results and get geotiff info
  ;
  file=folderpath + imagebase + '_LSTresults.tif'
  IF FILE_TEST(file) EQ 0 THEN BEGIN
    PRINT, 'ERROR! File not found: ' + file
    STOP
  ENDIF
  results=READ_TIFF(file,GEOTIFF=geotiff)
  Landsat_UL_East= geotiff.modeltiepointtag[3]
  Landsat_UL_North= geotiff.modeltiepointtag[4]
  zone = STRMID(STRCOMPRESS(STRING(geotiff.PROJECTEDCSTYPEGEOKEY)),4,2) 

  ;
  ; convert buoy lat/long to UTM, and retrieve pixel size. 
  ; should be 30x30 meters for Landsat.
  ;
  buoyUTM= CONVERT_LL_UTM(buoylat, buoylon, zone)
  buoyEast= buoyUTM[1]
  buoyNorth= buoyUTM[2] 
  pixsizeX = geotiff.MODELPIXELSCALETAG[0]
  pixsizeY = geotiff.MODELPIXELSCALETAG[1]

  ;
  ; calculate difference between Landsat UTM and buoy UTM (get answer in pixels)
  ; for easting, do buoy-Landsat (b/c buoy has to be to the right of the Landsat easting)
  ; for northing, do landsat-buoy (b/c buoy has to be below the Landsat northing)
  ;
  Easting_diff= ROUND((buoyEast - Landsat_UL_East)/pixsizeX)
  Northing_diff= ROUND((Landsat_UL_North - buoyNorth)/pixsizeY)

  ;
  ; select radiance values at and around buoy location
  ; change size if there are many "zero" pixels due to the broken SLC in Landsat 7
  ;
  pix_around_buoy= results[*,Easting_diff-2:Easting_diff+2,Northing_diff-2:Northing_diff+2]

  IF N_ELEMENTS(WHERE(pix_around_buoy[0,*,*] EQ 0 )) GT 15 THEN BEGIN
    pix_around_buoy= results[*,Easting_diff-5:Easting_diff+4,Northing_diff-5:Northing_diff+4]
  ENDIF

  ;
  ; if cloudmask keyword is set, get distance to nearest cloud
  ;
  IF KEYWORD_SET(getdist2cloud) THEN BEGIN
    cfmaskpath= folderpath + imagebase +'_cfmask.tif'
    cfmask=READ_TIFF(cfmaskpath,GEOTIFF=cmask_geotiff)
    cloudmask= GET_CLOUDMASK(folderpath,imagebase)
    minDistance= CALC_NEAREST_CLOUD(folderpath, imagebase, buoylat, buoylon, cloudmask, cmask_geotiff)
    params[6,0]= minDistance
  ENDIF
  
  ;
  ; retrieve atmospheric parameters for these pixels
  ; eliminate pixels that have a value of zero
  ; only keep values that the mean-value is less than the standard deviation
  ;
  radiance_therm = pix_around_buoy[0,*,*]
  nonzero= where(radiance_therm NE 0)
  rad_therm_nozeros= radiance_therm[nonzero]
  rad_therm_keep = WHERE(ABS(rad_therm_nozeros-MEAN(rad_therm_nozeros)) LT STDDEV(rad_therm_nozeros))
  rad_therm_avg = MEAN(rad_therm_nozeros[rad_therm_keep])

  transmission = pix_around_buoy[2,*,*]
  trans_nozeros= transmission[nonzero]
  trans_keep = WHERE(ABS(trans_nozeros-MEAN(trans_nozeros)) LT STDDEV(trans_nozeros))
  trans_avg = MEAN(trans_nozeros[trans_keep])

  upwelled = pix_around_buoy[3,*,*]
  up_nozeros= upwelled[nonzero]
  up_keep = WHERE(ABS(up_nozeros-MEAN(up_nozeros)) LT STDDEV(up_nozeros))
  up_avg = MEAN(up_nozeros[up_keep])

  downwelled = pix_around_buoy[4,*,*]
  down_nozeros= downwelled[nonzero]
  down_keep = WHERE(ABS(down_nozeros-MEAN(down_nozeros)) LT STDDEV(down_nozeros))
  down_avg = MEAN(down_nozeros[down_keep])

  ;
  ; pick emissivity of water based on whichLandsat
  ;
  CASE whichLandsat OF
    5:  emissivity = 0.98939488
    7:  emissivity = 0.98996972
    10: emissivity = 0.99243912
    11: emissivity = 0.98595542
  ENDCASE

  ;
  ; get lookup table based on whichLandsat
  ;
  CASE whichLandsat OF
    5:  LUTfile = folderpath+'LUT5.txt'
    7:  LUTfile = folderpath+'LUT7.txt'
    10: LUTfile = folderpath+'LUT8_B10.txt'
    11: LUTfile = folderpath+'LUT8_B11.txt'
  ENDCASE

  ;
  ; import LUT
  ;
  LUT = MAKE_ARRAY(2,500) 
  line = ''
  count=-1
  OPENR, lun, LUTfile, /GET_LUN  
  WHILE NOT EOF(lun) DO BEGIN & $
    count=count+1
    READF, lun, line & $
    LUT[0,count]= FLOAT(STRSPLIT(line, string(' '), /EXTRACT))
  ENDWHILE
  LUT=LUT[*,0:count]
  CLOSE, lun
  FREE_LUN, lun

  ;
  ; THE BIG EQUATION, and use LUT to get temperature
  ;
  rad_observed = (((rad_therm_avg - up_avg)/trans_avg)-(1-emissivity)*down_avg)/emissivity
  temperature = CONVERT_RAD_TEMP(rad_observed,LUT)

  params[0]=rad_therm_avg
  params[1]=trans_avg
  params[2]=up_avg
  params[3]=down_avg
  params[4]=rad_observed
  params[5]=temperature

  ; print results
  ; PRINT, '# of "zero" pixels in 5x5 around buoy= ', N_ELEMENTS(WHERE(pix_around_buoy[0,*,*] EQ 0 ))
  ; PRINT, 'Band6= ', rad_therm_avg
  ; PRINT, 'transmission= ', trans_avg
  ; PRINT, 'upwelled= ', up_avg
  ; PRINT, 'downwelled= ', down_avg
  ; PRINT, 'Radiance Observed= ',rad_observed
  ; PRINT, 'Temp= ', temperature

  RETURN, params
END

