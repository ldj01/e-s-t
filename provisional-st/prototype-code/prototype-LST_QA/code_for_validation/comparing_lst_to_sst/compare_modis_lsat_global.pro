; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:
;    COMPARE_MODIS_LSAT_GLOBAL
;
; PURPOSE:
;    IDL PROCEDURE
;    For a pair of Landsat 7 scenes and corresponding MODIS SST scenes, find up to seven locations that would be
;    appropriate for comparing temperature values (LST for Landsat and SST for MODIS). These are the follwing "rules"
;    for choosing "appropriate" comparison locations:
;       1) There must be a 10x10 area that has no land pixels or pixels outside the Landsat scene
;       2) The center or near center pixel of this area must be best quality
;       3) The surrounding 3x3 SST values must not have a STD > 0.5 K
;       4) The surrounding 3x3 quality values must have at least 5 best quality pixels,
;
;
; OUTLINE:
;    1) Read in Landsat LST file, MODIS SST and qual files
;    2) Find best quality pixels in MODIS scene and pick largest area that is at least 3 pixels in each dimension
;    3) Use 3x3 averaging kernel in the area
;    4) Get lat/lon locations for everypixel and locate them in Landsat scene
;    5) Use 5x5 averaging kernel on Landsat scene and retrieve values at lat/lon locations
;
;
; CALL SEQUENCE:
;    COMPARE_MODIS_LSAT_GLOBAL,    lsat_folderpath, $          ; path where landsat LST files are (include slash)
;                                  modis_folderpath, $         ; path where MODIS SST and QUAL bands are (include slash)
;				   lsat_scenename, $           ; name of Landsat LST scene (no extensions) 
;				   modis_scenename, $          ; name of MODIS scene (no extensions)
;                                  lut_path, $                 ; path and name of Ladnsat LUT file (include extension)
;                                  datadump_path, $            ; path where temporary output files go (include slash)
;                                 
;
;
; RESTRICTIONS:
;    N/A
;
; REQUIRED PROGRAMS (in working directory):
;    
;
; REQUIRED FILES (in folderpath directory):
;    
;
; MODIFICATIONS:
;    December, 2015      Original code
;

PRO COMPARE_MODIS_LSAT_GLOBAL, lsat_folderpath, modis_folderpath, lsat_scenename, modis_scenename, lut_path, datadump_path


lst_sst_data=FLTARR(5,1)

   lsat_file   = lsat_scenename
   modis_file  = modis_scenename

  PRINT, ' '
  PRINT, "lsat_file = ",lsat_file, "     modis_file = ", modis_file
  PRINT, ' '
  
  
  ;
  ; check if LST file is present (LST might not have been generated for all files in list)
  ; else read in LST results and band6
  ;
  lsat_endfile = "_LSTparams_MERRA.tif"
  lsat_lstfile = lsat_folderpath + lsat_file + lsat_endfile
  exists = FILE_TEST(lsat_lstfile)
  IF exists EQ 0 THEN BEGIN
    PRINT, 'file ',lsat_file+lsat_endfile, ' does not exist, skipping...'
    GOTO, SKIP_SCENE
  ENDIF ELSE BEGIN 
    okay = QUERY_TIFF(lsat_lstfile)
    IF (okay) THEN BEGIN
      lst_array = READ_TIFF(lsat_lstfile,GEOTIFF=lst_geotiff)
      band6 = REFORM(lst_array[0,*,*])    
    ENDIF ELSE BEGIN
      PRINT, 'Landsat file is not good. Skipping ...'
      GOTO, SKIP_SCENE
    ENDELSE
  ENDELSE


  ;
  ; check if sst file is present 
  ;
  modis_sstfile = modis_folderpath + modis_file + '_sst.tif'
  exists = FILE_TEST(modis_sstfile)
  IF exists EQ 0 THEN BEGIN
    PRINT, 'file ',modis_file+'_sst.tif', ' does not exist, skipping...'
    GOTO, SKIP_SCENE
  ENDIF ELSE BEGIN
    sst_file  = modis_folderpath + modis_file + '_sst.tif'
    qual_file = modis_folderpath + modis_file + '_qual.tif'
    sst  = READ_TIFF(sst_file,GEOTIFF=sst_geotiff)
    qual = READ_TIFF(qual_file)
  ENDELSE

  
 
  

  ;
  ; Skip scene if there are no good quality pixels
  ;
  find_zeros = WHERE(qual EQ 0D,count)
  IF count GT 0 THEN BEGIN
    best_pix = ARRAY_INDICES(qual,find_zeros)
  ENDIF ELSE BEGIN
    PRINT, ' '
    PRINT, '================================================='
    PRINT, 'No best quality pixels in image, skipping ...'
    PRINT, '================================================='
    PRINT, ' '
    GOTO, SKIP_SCENE
  ENDELSE
  
  ;
  ; check if whole quality image is zeros, which means something is wrong with the image
  ;
  IF N_ELEMENTS(WHERE(qual EQ 0)) EQ N_ELEMENTS(qual) THEN BEGIN
    PRINT, ' ' 
    PRINT, '================================================='
    PRINT, 'SST image is bad, skipping ....'
    PRINT, '================================================='
    PRINT, ' '
    GOTO, SKIP_SCENE
  ENDIF
  
  ;
  ; Obtain SST values and locations for an area of best quality pixels
  ;
  numcols = 10
  numrows = 10
  
  ;
  ; Find where zeros are in Landsat image
  ; 
  zeros_mask = MAKE_L7_ZEROS_MASK(band6);, lst_geotiff)
  
  ;
  ; Make new qual image by turning areas outside the Landsat scene to a -1
  ;
  new_qual_image = MODIFY_QUALVALS_WITH_LSAT(qual, sst_geotiff, lst_geotiff, zeros_mask)
  
  ;
  ; Find indexes within MODIS image that will be used for comparison
  ;
  indexes_and_utms = FIND_COMPARE_PTS_IN_MODIS(sst,new_qual_image,sst_geotiff,numcols,numrows)
  
  
  IF N_ELEMENTS(WHERE(indexes_and_utms EQ 0.0)) EQ N_ELEMENTS(indexes_and_utms) THEN BEGIN
     PRINT, ' '
     PRINT, '================================================='
     PRINT, 'No useable pixels in this image set. Skipping ... '
     PRINT, '================================================='
     PRINT, ' '
     GOTO, SKIP_SCENE  
  ENDIF ELSE BEGIN
    where_valid = WHERE(indexes_and_utms[0,*] NE 0.0,/NULL)
    valid_data = indexes_and_utms[*,where_valid]
    modis_column = valid_data[0,*]
    modis_row = valid_data[1,*]
    eastings = valid_data[2,*]
    northings = valid_data[3,*]
    stdev_of_best = valid_data[4,*]
 ENDELSE
  
  ; get sst vals for each spot and save image of modis scene with selected areas
  modis_with_boxes = new_qual_image
  sst_avg    = FLTARR(1,N_ELEMENTS(valid_data[0,*]))
  qual_avg   = FLTARR(1,N_ELEMENTS(valid_data[0,*]))
  stdev_best = FLTARR(1,N_ELEMENTS(valid_data[0,*]))
  num_best   = INTARR(1,N_ELEMENTS(valid_data[0,*]))
  
  FOR pts = 0, N_ELEMENTS(valid_data[0,*])-1 DO BEGIN
    pt_col = valid_data[0,pts]
    pt_row = valid_data[1,pts]
    pt_col_grow = INDGEN(3,1) + pt_col - 1
    pt_row_grow = INDGEN(3,1) + pt_row - 1
    nx=N_ELEMENTS(pt_col_grow) & ny=N_ELEMENTS(pt_row_grow)
    cols_mesh = REFORM( pt_col_grow#REPLICATE(1,ny),  N_ELEMENTS(pt_row_grow)*N_ELEMENTS(pt_col_grow), 1 )
    rows_mesh = REFORM( pt_row_grow##REPLICATE(1,nx), N_ELEMENTS(pt_row_grow)*N_ELEMENTS(pt_col_grow), 1 )
    modis_with_boxes[cols_mesh, rows_mesh] = 5
    ;
    ; get sst values
    ;
    sstvals = sst[cols_mesh,rows_mesh]*0.005 +273.15
    sstqual = new_qual_image[cols_mesh,rows_mesh]
    bestqual = WHERE(sstqual EQ 0,pixcount,/NULL)
    sst_avg[pts] = MEAN(REFORM(sstvals[bestqual],pixcount,1))
    qual_avg[pts] = MEAN(REFORM(sstqual,N_ELEMENTS(sstqual),1))
    stdev_best[pts] = stdev_of_best[pts]
    num_best[pts] = pixcount
    
  ENDFOR

  
  ;
  ; Find params for UTM locations
  ;
  ;lut_path = main_folderpath + 'LUT7.txt'
  params = GET_LST_AT_UTMS(lst_array, lst_geotiff, eastings, northings, lut_path)
  lsurf = params[4,*]
  trans = params[1,*]
  upwell = params[2,*]
  downwell = params[3,*]
  temps = params[5,*] 
  

  
  
  lst_sst_data = FLTARR(11,N_ELEMENTS(sst_avg))
  lst_sst_data[0,*] = FLOAT(temps)
  lst_sst_data[1,*] = FLOAT(sst_avg)
  lst_sst_data[2,*] = FLOAT(ROUND(eastings))
  lst_sst_data[3,*] = FLOAT(ROUND(northings))
  lst_sst_data[4,*] = FLOAT(qual_avg)
  lst_sst_data[5,*] = FLOAT(stdev_best)
  lst_sst_data[6,*] = FIX(num_best)
  lst_sst_data[7,*] = FLOAT(lsurf)
  lst_sst_data[8,*] = FLOAT(trans)
  lst_sst_data[9,*] = FLOAT(upwell)
  lst_sst_data[10,*] = FLOAT(downwell)
  
  
SKIP_SCENE:

data_filename = lsat_file + '_lst_sst_data.txt'
OPENW, lun, datadump_path + data_filename, /GET_LUN
  PRINTF, lun, FORMAT = '(F12," ",F12," ",F16," ",F16, " ",F12," ",F12," ",I1, " ", F12, " ", F12, " ", F12, " ", F12)', lst_sst_data
CLOSE, lun
FREE_LUN, lun

  
END
  
 
 
