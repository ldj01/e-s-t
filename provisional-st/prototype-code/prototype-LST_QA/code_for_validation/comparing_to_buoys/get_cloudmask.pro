; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    GET_CLOUDMASK
;
; PURPOSE:  
;    IDL FUNCTION
;    Using the Landsat reflectance product, a binary mask of cloud locations is created in order to determine if there are clouds
;    above a buoy. 
;
; OUTLINE:
;    1) import cfmask.tif
;    
;       If using Landsat 7:                                               If NOT using Landsat 7:
;         2) import sr_fill_qa.tif, the binary image of slc gaps             2) get rid of cloud shadow so that the cloudmask only has 
;         3) invert values in gap mask, and give gap locations a                values 0 (for no clouds), 156 (gaps), and 255 (clouds).
;            value of 100 to tell the difference between gap pixels          3) return cloudmask
;            and no-cloud pixels, which have a value of 0.
;         4) multiply gapmask and cfmask
;         5) get rid of cloud shadow so that the cloudmask only has
;            values 0 (for no clouds), 156 (gaps), and 255 (clouds).
;         6) return cloudmask
;    
;    
; CALL SEQUENCE: 
;    GET_CLOUDMASK,   folderpath, $               ; path where input files are (include '\' after folder)
;                     imagebase, $                ; scene identifier (e.g. LE70160302008113EDC00)
;                                     
; RESTRICTIONS:
;    1) hard coded paths
;    2) The gaps are given a byte value of 100 so the difference between no clouds and the gap can be known.
;
; REQUIRED PROGRAMS (in working directory):
;    N/A
;
; REQUIRED FILES (in folderpath directory):      
;    scene_cfmask.tif
;    scene_sr_fill_qa.tif (if using L7)
;    
; MODIFICATIONS:
;    December, 2014       Original code
;

FUNCTION GET_CLOUDMASK, folderpath, imagebase

  ;              
  ; import cfmask 
  ;               
  cfmaskpath= folderpath + imagebase +'_cfmask.tif'
  cfmask=READ_TIFF(cfmaskpath)
  
  HELP, cfmask
  PRINT, min(cfmask)
  PRINT, max(cfmask)
  
  ;                                                
  ; if using Landsat 7, must import sr_fill_qa.tiff
  ;                                                 
  whichLandsat=STRMID(imagebase,2,1)

  IF whichLandsat EQ '7' THEN BEGIN
  
    s = SIZE(cfmask)
    ncol = s[1]

    ;
    ; Find where there are clouds (values of 4), fill value is 255
    ;
    cloudind = WHERE(cfmask EQ 4)
    cols = cloudind MOD ncol
    rows = cloudind / ncol

    gapind = WHERE(cfmask EQ 255)
    gapcols = gapind MOD ncol
    gaprows = gapind / ncol

    cloudmask=BYTARR(s[1],s[2])
    cloudmask[*]=0
    cloudmask[gapcols,gaprows]=156
    IF cloudind[0] EQ -1 THEN BEGIN
      RETURN, cloudmask
    ENDIF
    cloudmask[cols,rows]=255

    RETURN, cloudmask
  ENDIF ELSE BEGIN
  
    ;
    ; if NOT using Landsat 7, do the following.
    ;
    s = SIZE(cfmask)
    ncol = s[1]
    cloudind = WHERE(cfmask EQ 4)
    IF cloudind[0] EQ -1 THEN BEGIN
      cloudmask[*]=0
      RETURN, cloudmask
    ENDIF
    cols = cloudind MOD ncol
    rows = cloudind / ncol
  
    cloudmask=BYTARR(s[1],s[2])
    cloudmask[*]=0
    cloudmask[cols,rows]=255
  
    RETURN, cloudmask
  ENDELSE

END
