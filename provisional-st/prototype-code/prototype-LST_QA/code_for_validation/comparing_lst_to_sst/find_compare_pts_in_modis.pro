; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
;
; NAME:
;    GET_SST_FOR_AREA
;
; PURPOSE:
;    IDL FUNCTION
;    This code finds up to 7 good comparison points in the SST scene and returns info about the locations and the SST values.
;    
;
; OUTLINE:
;    1) read in MODIS SST and qual bands with geotiff info
;    2) find where best quality pixels are
;    3) calculate largest area of best pixels
;    4) use 3x3 averaging kernel on area
;    5) get sst values for line of pixels within area
;    6) get UTM locations for each pixel
;    7) return sst and qual values, return lat lons
;
;
; CALL SEQUENCE:
;    GET_SST_AREA,      folderpath, $                   ; path where modis files are (include '\' after folder)
;                       modis_scene                     ; MODIS scene name (no extension)
;
; RESTRICTIONS:
;    N/A
;    
;    
; REQUIRED PROGRAMS (in working directory):

;
;
; REQUIRED FILES (in folderpath directory):

;
; MODIFICATIONS:
;    December, 2015       Original code
;
;
FUNCTION FIND_COMPARE_PTS_IN_MODIS, sst,new_qual_image, sst_geotiff, numcols, numrows

  ;indexes_and_utms = 0


  max_qual_allowed = 0D
  min_goodpix_allowed=5
  total_points=7

  modis_ul_east  = sst_geotiff.modeltiepointtag[3]
  modis_ul_north = sst_geotiff.modeltiepointtag[4]
  zone = STRMID(STRCOMPRESS(STRING(sst_geotiff.PROJECTEDCSTYPEGEOKEY)),4,2)
  pixsize = sst_geotiff.MODELPIXELSCALETAG[0]
  image_size = SIZE(sst)

  IF modis_ul_north LT 0 THEN BEGIN
  modis_ul_north = modis_ul_north + 10000000
  ENDIF

  
  ;
  ; Choose size of subset to contain all zeros (should be multiples of 3, maybe)
  ;
  area_numcols = numcols
  area_numrows = numrows
  area = MAKE_ARRAY(area_numcols,area_numrows,/FLOAT)
  
  ;
  ; Figure out how many subsets can fit across and down the image
  ; Use pad to center the subsets
  ;
  col_numfits = FLOOR(DOUBLE(image_size[1]) / DOUBLE(area_numcols))
  col_pad = FLOOR( (image_size[1] - col_numfits*area_numcols) / 2 )
  row_numfits = FLOOR(DOUBLE(image_size[2]) / DOUBLE(area_numrows)) 
  row_pad = FLOOR( (image_size[2] - row_numfits*area_numrows) / 2 )
  
  ;
  ; Subtract one (unless already zero) for indexing purposes. Basically, 
  ; if the subset doesn't fit evenly I pad on either side so that the 
  ; group of subsets are centered ( I don't really to find a group of zeros
  ; on the edge anyway because it could be outside the Landsat scene.
  ;
  IF col_pad GT 0 THEN col_pad = col_pad - 1
  IF row_pad GT 0 THEN row_pad = row_pad - 1 
  
  
  ;
  ; Loop through rows and columns looking for a subset-sized area of zeros
  ;
  indexes_and_utms = FLTARR(5,total_points)
  numsamples = 0
  FOR row = row_pad, row_pad + (row_numfits-1)*area_numrows, area_numrows DO BEGIN
    
    FOR column = col_pad, col_pad + (col_numfits-1)*(area_numcols) DO BEGIN
      
      ; retrieve current subset and calculate average quality
      check_area = new_qual_image[column:column+area_numcols-1,row:row+area_numrows-1]     
      sst_area = sst[column:column+area_numcols-1,row:row+area_numrows-1]
      ;avg_qual_area = MEAN(REFORM(check_area,area_numrows*area_numcols,1))
      
      ; make qual image showing selected area as black box (will bw saved later)
      img=new_qual_image
      img[column:column+area_numcols-1,row:row+area_numrows-1]=5
     
      ; if any of the qualities are -1 or 4 then move to next area. Also make sure no NANs in SST
      IF N_ELEMENTS(WHERE( REFORM(check_area,area_numrows*area_numcols,1) EQ -1,/NULL)) GT 0 THEN CONTINUE
      IF N_ELEMENTS(WHERE( REFORM(check_area,area_numrows*area_numcols,1) EQ  4,/NULL)) GT 0 THEN CONTINUE
      IF N_ELEMENTS(WHERE( FINITE(sst_area, /NAN) EQ 1, /NULL)) GT 0 THEN CONTINUE
      
   
        ; figure out which columns and rows make up the sample within subset (find center pixel)
        numpix = 1
        sample_rows = INTARR(numpix,1) + row + FLOOR(area_numrows/2D)
        sample_cols = INDGEN(numpix,1) + column + FLOOR((area_numcols - numpix)/2D)
        
        ; check if enough pixels nearby are good quality ( 5 or more )
        FOR check = 0, 4 DO BEGIN
          ; find 3x3 area around sample to check for sst variations
          grow_rows = INDGEN(3,1) + sample_rows[0] - 1
          grow_cols = INTARR(numpix+2,1)
          grow_cols[0] = MIN(sample_cols) - 1
          grow_cols[numpix+1] = MAX(sample_cols) + 1
          grow_cols[1:numpix] = sample_cols

          nx=N_ELEMENTS(grow_rows) & ny=N_ELEMENTS(grow_cols)
          grow_rows_mesh = REFORM( grow_rows#REPLICATE(1,ny),  N_ELEMENTS(grow_cols)*N_ELEMENTS(grow_rows), 1 )
          grow_cols_mesh = REFORM( grow_cols##REPLICATE(1,nx), N_ELEMENTS(grow_cols)*N_ELEMENTS(grow_rows), 1 )
         
          current_sst_3by3  = sst[grow_cols_mesh, grow_rows_mesh]
          current_qual_3by3 = new_qual_image[grow_cols_mesh, grow_rows_mesh]
          best_pixels_index = WHERE(current_qual_3by3 LE max_qual_allowed AND current_qual_3by3 GE 0)
          IF N_ELEMENTS(best_pixels_index) GE min_goodpix_allowed THEN BREAK
          ;else
          sample_cols = sample_cols + 1
        ENDFOR
        
        ; If none of the five options were good, continue to next iteration
        IF N_ELEMENTS(best_pixels_index) LT min_goodpix_allowed THEN CONTINUE
        
        ; check if surrounding area  has highly varying temperatures
        grow_sst = REFORM(sst[grow_cols_mesh,grow_rows_mesh],N_ELEMENTS(grow_cols)*N_ELEMENTS(grow_rows),1)
        best_pixels_vector = REFORM(current_sst_3by3, N_ELEMENTS(current_sst_3by3), 1)
        
        stdev_of_best_pix = STDEV(best_pixels_vector*.005)
	stdev_of_area = STDEV(grow_sst*0.005)
        ;IF stdev_of_best_pix GT 0.5D THEN BEGIN 
        IF stdev_of_area GT 0.5D THEN BEGIN
          PRINT,'High Temperature gradiant. Moving on ...'
          CONTINUE
        ENDIF
        
             
        numsamples = numsamples + 1
        IF numsamples GT total_points THEN CONTINUE
        
        
        ; turn the indices into meshgrids so we can index the whole image (only if dims are different sizes)
;        IF N_ELEMENTS(sample_rows) NE N_ELEMENTS(sample_cols) THEN BEGIN
;          nx=n_elements(sample_rows) & ny=n_elements(sample_cols)
;          sample_rows_mesh = REFORM( sample_rows#REPLICATE(1,ny),  N_ELEMENTS(sample_cols)*N_ELEMENTS(sample_rows), 1 ) 
;          sample_cols_mesh = REFORM( sample_cols##REPLICATE(1,nx), N_ELEMENTS(sample_cols)*N_ELEMENTS(sample_rows), 1 )
;          sample_rows = sample_rows_mesh
;          sample_cols = sample_cols_mesh
;        ENDIF
        
        indexes_and_utms[0,numsamples-1] = sample_cols
        indexes_and_utms[1,numsamples-1] = sample_rows

        ;
        ; Find UTM points for pixel in sample
        ;
        pixel_eastings = modis_ul_east + pixsize*sample_cols
        pixel_northings = modis_ul_north - pixsize*sample_rows

        indexes_and_utms[2,numsamples-1] = pixel_eastings
        indexes_and_utms[3,numsamples-1] = pixel_northings
        indexes_and_utms[4,numsamples-1] = stdev_of_area
        

        column = column + area_numcols

    ENDFOR
  ENDFOR
  
  IF numsamples EQ 0 THEN BEGIN
    indexes_and_utms = 0
  ENDIF
  RETURN, indexes_and_utms
  
  
END

