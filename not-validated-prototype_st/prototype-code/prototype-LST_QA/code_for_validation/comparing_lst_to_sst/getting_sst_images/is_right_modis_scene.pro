; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    IS_RIGHT_MODIS_SCENE
;
; PURPOSE:  
;    IDL PROCEDURE
;    This program reads in a MODIS SST file, and compares its corner coordinates to given Landsat corner coordinates,
;    and determines if the Landsat scene is encapsulated by the MODIS scene. If it is, a "1" is written to a temporary
;    file, otherwise a "0" is written there. Ideal for use with get_modis_for_L7_scenes.bash    
;   
; OUTLINE:   
;    1) Read in MODIS SST file
;    2) Get upper left and lower right lat/lons for MODIS using NCDF routines
;    3) Compare MODIS corner lat/lons to given Landsat corner lat/lons
;    4) If Landsat is within MODIS scene, write "1" to file; otherwise write a "0"
;
; CALL SEQUENCE: 
;    IS_RIGHT_MODIS_SCENE,   modis_folder, $                ; path where MODIS SST file is
;                            modis_filename, $              ; scene name for MODIS SST file
;                            ul_lat, $                      ; upper left lattitude for Landsat scene
;                            ul_lon, $                      ; upper left longitude for Landsat scene
;                            lr_lat, $                      ; lower right lattitude for Landsat scene
;                            lr_lon                         ; lower right longitude for Landsat scene
;                            
;                                        
; RESTRICTIONS:
;    1) This program was designed to be called by get_modis_for_L7_scenes.bash. It can certainly be used 
;       independently, but the fact that the result of the program is written to a file may not be ideal.
;
; REQUIRED PROGRAMS (in working directory):      
;    NONE                  
;                       
; REQUIRED FILES (in modis_folder directory):      
;    MODIS SST file
;
; MODIFICATIONS:
;    May,  2015       Original code
;    June, 2016       MODIS SST file format changed, so HDF routines were changed to NCDF routines



PRO is_right_modis_scene, modis_folder, modis_filename, ul_lat, ul_lon, lr_lat, lr_lon

; read in modis file
id_modis = NCDF_OPEN(modis_folder + '/' + modis_filename)

;
; get start lat/lons and end lat/lons (start is first line, end is last line) from MODIS image,
; which are located within certain attributes
;
id_beglon = NCDF_ATTNAME(id_modis, /GLOBAL, 40) ; westlon  is attribute # 40
id_endlon = NCDF_ATTNAME(id_modis, /GLOBAL, 39) ; eastlon  is attribute # 39
id_beglat = NCDF_ATTNAME(id_modis, /GLOBAL, 37) ; northlat is attribute # 37
id_endlat = NCDF_ATTNAME(id_modis, /GLOBAL, 38) ; southlat is attribute # 38

NCDF_ATTGET, id_modis, /GLOBAL, id_beglon, beglon  
NCDF_ATTGET, id_modis, /GLOBAL, id_endlon, endlon 
NCDF_ATTGET, id_modis, /GLOBAL, id_beglat, beglat 
NCDF_ATTGET, id_modis, /GLOBAL, id_endlat, endlat


; give more intuitive names
UL_LAT_MO = beglat
UL_LON_MO = beglon 
LR_LAT_MO = endlat 
LR_LON_MO = endlon 


; landsat latlons
UL_LAT_LA = FLOAT(ul_lat)
UL_LON_LA = FLOAT(ul_lon)
LR_LAT_LA = FLOAT(lr_lat)
LR_LON_LA = FLOAT(lr_lon)


; Check if Landsat corner coordinates are encapsulated my MODIS corner coordinates.
; If they are, a 1 gets written to a temporary text file, otherwise a 0 gets written.
answer=0
IF ( UL_LAT_MO gt UL_LAT_LA ) AND ( UL_LON_MO lt UL_LON_LA) AND ( LR_LAT_MO lt LR_LAT_LA ) AND ( LR_LON_MO gt LR_LON_LA) THEN BEGIN
  answer=1
ENDIF


; write answer to temporary file
ans=STRING(answer)
OPENW, unit, modis_filename+'.txt', /GET_LUN
PRINTF, unit, ans
CLOSE, unit
FREE_LUN, unit

END
