; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    ATMPARAMS_BUOYS
;
; PURPOSE:  
;    IDL PROCEDURE
;    This code is designed to calculate the atmospheric parameters for specific buoys within Landsat scenes. This allows a list of scenes to be
;    run automatically and writes a text file with the atmospheric parameters and other info for each scene in the list. 
;    
;   
; OUTLINE:   
;    1) Read in scenelist.txt, whichBuoyInScene.txt, and buoyHistory.txt
;    2) For each scene in scene list, extract the Landsat number, path, row, and the acquisition date
;    3) Find which buoy station(s) are in each scene and use check_buoy_history to obtain the appropriate lat/long based on date
;    4) Run get_atm_params for each buoy to obtain atmospheric parameters
;    5) Write .txt file that shows scene name, buoy #, buoy lat/long, the 5 atmo. params, and the predicted temperature
;
; CALL SEQUENCE: 
;    atmparams_buoys,    folderpath, $                    ; path where read-in files and result.tifs are (include '\' after folder)
;                        results_filename, $              ; desired name of output .txt file (include .txt)
;                        MYLATLONS=mylatlons, $           ; if set, program will skip checking buoy history and will use provided lat/lons
;                            
;                                        
; RESTRICTIONS:
;    1) Code will only work for a maximum of 2 buoys per scene. Will be increased as necessary.
;    2) read-in files are currently hard-coded in, so their names should not be altered.
;
; REQUIRED PROGRAMS (in working directory):      
;    GET_ATM_PARAMS.PRO   
;    CONVERT_LL_UTM.PRO
;    CONVERT_RAD_TEMP.PRO
;    CHECK_BUOY_HISTORY.PRO
;    GET_CLOUD_MASK.PRO (optional)
;                                                 
; REQUIRED FILES (in folderpath directory):      
;    scenelist.txt
;    whichBuoyInScene.txt
;    buoyHistory.txt
;    LUT7.txt
;    cloud mask files (optional. see create_cloud_mask.pro)   
;    buoylatlons.txt (optional)                                             
;
; MODIFICATIONS:
;    November, 2014       Original code
;    February, 2015       Added capibility to get distance to nearest cloud and write it to file
;    February, 2015       Added /mylatlons keyword so buoys that aren't in the history file can be used
;
PRO atmparams_buoys, mainfolder, results_filename, MYLATLONS=mylatlons

  which_anal='_merra/'
  ;
  ; import scenelist
  ;
  scenelistpath= mainfolder + 'scenelist.txt'
  nlines = FILE_LINES(scenelistpath)
  scenelist = MAKE_ARRAY(1,nlines, /STRING)
  line = ''
  count=-1
  OPENR, lun, scenelistpath, /GET_LUN
  WHILE NOT EOF(lun) DO BEGIN & $
    count=count+1
    READF, lun, line & $
    scenelist[0,count]=line
  ENDWHILE
  CLOSE, lun
  FREE_LUN, lun


 
  IF KEYWORD_SET(mylatlons) THEN BEGIN
    
    ;
    ; import file that lists buoy number and corresponding lat/lon
    ;
    buoylatlonpath= mainfolder + 'buoylatlon.txt'
    nlines = FILE_LINES(buoylatlonpath)
    buoylatlons = FLTARR(5,nlines)
    line = ''
    count=-1
    OPENR, lun, buoylatlonpath, /GET_LUN
    WHILE NOT EOF(lun) DO BEGIN & $
      count=count+1
      READF, lun, line & $
      buoylatlons[0,count]=FLOAT(STRSPLIT(line, STRING(9B), /EXTRACT))
    ENDWHILE
    CLOSE, lun
    FREE_LUN, lun
    
    ; initiate variables that will be inputs for get_atm_params.pro
    filename=STRARR(1,500)
    folderpath=STRARR(1,500)
    buoylat=DBLARR(1,500)
    buoylon=DBLARR(1,500)
    buoynum=INTARR(1,500)
    whichLandsat=INTARR(1,500)


    k=0
    FOR i= 0,N_ELEMENTS(scenelist)-1 DO BEGIN
      filename_temp=scenelist[i]
      whichLandsat_temp=STRMID(scenelist[0,i],2,1)
      path_temp=STRMID(scenelist[0,i],4,2)
      row_temp=STRMID(scenelist[0,i],7,2)
      
      ; check which buoys are in the scene, and whether they should be used
      ind=WHERE(path_temp EQ FIX(buoylatlons[0,*]) AND row_temp EQ FIX(buoylatlons[1,*]))
      IF ind[0] EQ -1 THEN BEGIN
        PRINT, 'ERROR: BUOY FOR PATH/ROW NOT IN FILE!'
        STOP
      ENDIF
      ; if only one buoy
      IF N_ELEMENTS(ind) EQ 1 THEN BEGIN
        buoy=buoylatlons[2,ind]
        buoynum[0,k]=buoy
        buoylat[0,k]=buoylatlons[3,ind]
        buoylon[0,k]=buoylatlons[4,ind]
        filename[0,k]=filename_temp
        folderpath[0,k]=mainfolder + path_temp + '_' + row_temp + '/' + filename[0,k] + which_anal
        whichLandsat[0,k]= whichLandsat_temp
        k=k+1
        CONTINUE
      ENDIF ELSE BEGIN
        
        ; case where 2 buoys are present
        ; first buoy
        buoy=buoylatlons[2,ind[0]]
        buoynum[0,k]=buoy
        buoylat[0,k]=buoylatlons[3,ind[0]]
        buoylon[0,k]=buoylatlons[4,ind[0]]
        filename[0,k]=filename_temp
        folderpath[0,k]=mainfolder + path_temp + '_' + row_temp + '/' + filename[0,k] + which_anal
        whichLandsat[0,k]= whichLandsat_temp

        ; second buoy
        k=k+1
        buoy=buoylatlons[2,ind[1]]
        buoynum[0,k]=buoy
        buoylat[0,k]=buoylatlons[3,ind[1]]
        buoylon[0,k]=buoylatlons[4,ind[1]]
        filename[0,k]=filename_temp
        folderpath[0,k]=mainfolder + path_temp + '_' + row_temp + '/' + filename[0,k] + which_anal
        whichLandsat[0,k]= whichLandsat_temp
        k=k+1
        CONTINUE
        ENDELSE
    ENDFOR
    
  ENDIF ELSE BEGIN
  
  
  ;
  ; IF MYLATLONS NOT SET, DO THE FOLLOWING
  ;
  
  ;
  ; import file that labels which buoy station #'s correspond to which scenes
  ;
  whichbuoypath= mainfolder + 'whichBuoyInScene.txt'
  nlines = FILE_LINES(whichbuoypath) 
  buoysinscene = MAKE_ARRAY(5,nlines)
  line = ''
  count=-1
  OPENR, lun, whichbuoypath, /GET_LUN
  WHILE NOT EOF(lun) DO BEGIN & $
    count=count+1
    READF, lun, line & $
    buoysinscene[0,count]=FLOAT(STRSPLIT(line, string(9B), /EXTRACT))
  ENDWHILE
  CLOSE, lun
  FREE_LUN, lun
  
  
  ;
  ; import buoy history file
  ;
  buoyhistpath= mainfolder + 'BuoyHistory.txt'
  nlines = FILE_LINES(buoyhistpath)
  buoyhistfull = STRARR(11,nlines)
  line = ''
  count=-1
  OPENR, lun, buoyhistpath, /GET_LUN
  WHILE NOT EOF(lun) DO BEGIN & $
    count=count+1
    READF, lun, line & $
    buoyhistfull[0,count]=STRSPLIT(line, ',', /EXTRACT)
  ENDWHILE
  CLOSE, lun
  FREE_LUN, lun

  ;
  ; select only columns 0,1,2,6,7 from buoy history array
  ;
  s=SIZE(buoyhistfull)
  buoyhist=DBLARR(5,s[2]-1)
  buoyhist[0:2,*]=DOUBLE(buoyhistfull[0:2,1:*])
  buoyhist[3:4,*]=DOUBLE(buoyhistfull[6:7,1:*])

  ;
  ; initiate variables that will be inputs for get_buoy_params.pro
  ;
  filename=STRARR(1,500)
  whichLandsat=INTARR(1,500)
  buoylat=FLTARR(1,500)
  buoylon=FLTARR(1,500)
  buoynum=INTARR(1,500)
  yearday=INTARR(1,500)
  usebuoy=FLTARR(3,500)

  ;
  ; retrieve info from scene names (from scenelist), which will be input for get_buoy_params.pro
  ; k acts as a counter, so that indexing works for scenes with multiple buoys
  ;
  k=0
  FOR i= 0,N_ELEMENTS(scenelist)-1 DO BEGIN
    
    filename_temp=scenelist[i]
    whichLandsat_temp=STRMID(scenelist[0,i],2,1)
    path_temp=STRMID(scenelist[0,i],4,2)
    row_temp=STRMID(scenelist[0,i],7,2)
    yearday_temp=STRMID(scenelist[0,i],9,7)
    
    ;
    ; find which buoy(s) are in the scene
    ;
    ind=WHERE(path_temp EQ FIX(buoysinscene[0,*]) AND row_temp EQ FIX(buoysinscene[1,*]))
 
    ;
    ; if only one buoy
    ;
    IF N_ELEMENTS(ind) EQ 1 THEN BEGIN
      buoy=buoysinscene[2,ind]
      usebuoy[*,k]= CHECK_BUOY_HISTORY(buoyhist,buoy[0],yearday_temp)
      buoynum[0,k]=buoy
      buoylat[0,k]= usebuoy[1,k]
      buoylon[0,k]=usebuoy[2,k]
      filename[0,k]=filename_temp  
      whichLandsat[0,k]= whichLandsat_temp
      k=k+1
      CONTINUE
    ENDIF ELSE BEGIN
      
      ;
      ; case where 2 buoys are present
      ; first buoy
      ;
      buoy=buoysinscene[2,ind[0]]
      usebuoy[*,k]= CHECK_BUOY_HISTORY(buoyhist,buoy[0],yearday_temp)
      buoynum[0,k]=buoy
      buoylat[0,k]= usebuoy[1,k]
      buoylon[0,k]=usebuoy[2,k]
      filename[0,k]=filename_temp
      whichLandsat[0,k]= whichLandsat_temp
    
      ;
      ; second buoy
      ;
      k=k+1
      buoy=buoysinscene[2,ind[1]]
      usebuoy[*,k]= CHECK_BUOY_HISTORY(buoyhist,buoy[0],yearday_temp)
      buoynum[0,k]=buoy
      buoylat[0,k]= usebuoy[1,k]
      buoylon[0,k]=usebuoy[2,k]
      filename[0,k]=filename_temp
      whichLandsat[0,k]= whichLandsat_temp
      k=k+1
      CONTINUE
    ENDELSE
  
  ENDFOR

  ENDELSE
  ;
  ; truncate arrays
  ;
  filename=filename[*,0:k-1]
  buoynum=buoynum[*,0:k-1]
  buoylat=buoylat[*,0:k-1]
  buoylon=buoylon[*,0:k-1]
  whichLandsat=whichLandsat[*,0:k-1]
  
  IF KEYWORD_SET(mylatlons) EQ 0 THEN BEGIN
    usebuoy=usebuoy[*,0:k-1]
  ENDIF 
  
  IF KEYWORD_SET(getdist2cloud) THEN BEGIN
    params=FLTARR(7,k)
  ENDIF ELSE BEGIN
    params=FLTARR(6,k)
  ENDELSE

  ; this is for my bash stuff. Move LUT to each results folder
  lutpath='/home/kga1099/LUT5.txt'

  ;
  ; GET THE PARAMETERS
  ;
  FOR i=0,N_ELEMENTS(buoylat)-1 DO BEGIN

    command= 'cp ' + lutpath + ' ' + folderpath[i] + '/'
    SPAWN, command
  
    IF KEYWORD_SET(getdist2cloud) THEN BEGIN
      
      PRINT, 'Obtaining parameters for ',filename[i],'.....'
      params[*,i]= GET_ATM_PARAMS(folderpath,filename[i],buoylat[i],buoylon[i], whichLandsat[i], /getdist2cloud)
        
    ENDIF  ELSE BEGIN
      
      IF KEYWORD_SET(mylatlons) EQ 0 THEN BEGIN
        IF usebuoy[0,i] EQ 0 THEN BEGIN
          PRINT, STRJOIN(['Buoy',STRCOMPRESS(UINT(buoynum[i])),'not valid for scene',filename[i],'. Params given values of 9.999999.'])
          params[*,i]=9.999999
          CONTINUE
        ENDIF
      ENDIF
      
          PRINT, 'Obtaining parameters for ',filename[i],'.....'
          params[*,i]= GET_ATM_PARAMS(folderpath[i],filename[i],buoylat[i],buoylon[i], whichLandsat[i])
    
    ENDELSE
        
  ENDFOR
  


  ;
  ; write results to file
  ; 
  resultfile=mainfolder + results_filename
  OPENW,unit,resultfile, /GET_LUN
  FOR i=-1,k-1 DO BEGIN
    IF i EQ -1 THEN BEGIN
      
      ;
      ; Write header
      ;
      IF KEYWORD_SET(getdist2cloud) THEN BEGIN
        PRINTF, unit, FORMAT = '("Scene",T27,"Buoy #",T38,"Buoylat",T50,"Buoylon",T64,"Band6",T78,"Transmission",T95,"Upwelled",T108,"Downwelled",T123,"Obs Radiance",T140,"Temperature [K]",T160,"Dist to Nearest Cloud [km]")'
        CONTINUE
      ENDIF ELSE BEGIN
        PRINTF, unit, FORMAT = '("Scene",T27,"Buoy #",T38,"Buoylat",T50,"Buoylon",T64,"Band6",T78,"Transmission",T95,"Upwelled",T108,"Downwelled",T123,"Obs Radiance",T140,"Temperature [K]")'
        CONTINUE
      ENDELSE
      
    ENDIF
    
    ;
    ; write parameters
    ;
    IF KEYWORD_SET(getdist2cloud) THEN BEGIN
      PRINTF, unit, FORMAT = '(A21,T27,I5,T35,F10.4,T49,F10.4,T62,F10.6,T76,F10.6,T93,F10.6,T106,F10.6,T121,F10.6,T140,F10.6,T160,F15.4)', $
                            filename[i],UINT(buoynum[i]),buoylat[i],buoylon[i],params[0,i],params[1,i],params[2,i],params[3,i],params[4,i],params[5,i],params[6,i]
    ENDIF ELSE BEGIN
      PRINTF, unit, FORMAT = '(A21,T27,I5,T35,F10.4,T49,F10.4,T62,F10.6,T76,F10.6,T93,F10.6,T106,F10.6,T121,F10.6,T140,F10.6)', $
        filename[i],UINT(buoynum[i]),buoylat[i],buoylon[i],params[0,i],params[1,i],params[2,i],params[3,i],params[4,i],params[5,i]
    ENDELSE
        
  ENDFOR
  CLOSE, unit
  FREE_LUN, unit
  

end
