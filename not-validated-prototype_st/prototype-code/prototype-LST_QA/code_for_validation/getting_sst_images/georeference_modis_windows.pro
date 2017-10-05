; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    GEOREFERENCE_MODIS_WINDOWS
;
; PURPOSE:  
;    IDL PROCEDURE
;    This program wil georeference any given number of Ocean Color MODIS Sea Surface Temperature (SST) files using ENVI. 
;    Both the SST and the quality band are georeferenced and then subset to the area in the image where the corresponding
;    Landsat image overlaps. Images are saved as geotiffs, which can then be used to compare points between the SST images
;    and Landsat Land Surface Temperature predictions. Meant for use on Windows OS.
;    
;   
; OUTLINE:   
;    1) Read in text file that contains Landsat scene names, MODIS scene names, UTM zone, and the corner coordinates for 
;       the Landsat image
;    2) For each MODIS scene, georeference SST and quality bands so they are in UTM coordinates
;    3) Find upperleft and lowerright points of Landsat scene in MODIS scene
;    4) Subset MODIS image and save sst and quality bands as geotiffs
;
; CALL SEQUENCE: 
;    GEOREFERENCE_MODIS_WINDOWS,   modis_folderpath, $             ; path where MODIS SST files are (include '\' after folder)
;                                  sceneinfo_folderpath, $         ; path where the input text file that has 7 columns of info:
;                                                                    Landsat scenename, MODIS scenename, zone, Landsat UL corner lat,
;                                                                    Landsat UL corner lon, Landsat LR corner lat, Landsat LR corner lon
;                                  sceneinfo_filename              ; name of input text file (include '.txt')
;                            
;                            
;                                        
; RESTRICTIONS/LIMITATIONS:
;    1) The program expects the input text file to have 7 columns of information (listed in call sequence), and 1 header line.
;    2) The output files are automatically placed within a folder called "georeferenced" in the folder specified by the first arguement.
;    3) The output files are automatically named by adding "_sst.tif" or "_qual.tif" to the MODIS scene name.
;    4) User must already have ENVI Classic open.
;
; REQUIRED PROGRAMS (in working directory):        
;    CONVERT_LL_UTM.PRO
;    ENVI CLASSIC
; 
;                                                 
; REQUIRED FILES:      
;    sceneinfo_file     (text file found in "sceneinfo_folderpath")
;    MODIS SST file     (found in "modis_folderpath")                                           
;
;
; MODIFICATIONS:
;    May,  2015        Original code
;    June, 2016        Altered code to accomodate for new MODIS SST file format

PRO georeference_modis_windows, modis_folderpath, sceneinfo_folderpath, sceneinfo_filename

  
  ; ==========================================================
  ; STEP 1: IMPORT TEXT FILE CONTAINING LANDSAT AND MODIS INFO
  ; ==========================================================
  
  sceneinfo_complete_path = sceneinfo_folderpath + sceneinfo_filename
  nlines = FILE_LINES(sceneinfo_complete_path)
  
  ; Assumes there are 7 columns in the text file and one header line, which will not be included in this array 
  sceneinfo_entries = MAKE_ARRAY(7,nlines-1,/STRING)
  
  sceneinfo_lines = ''
  count = -1
  OPENR, lun, sceneinfo_complete_path, /GET_LUN
    
    WHILE NOT EOF(lun) DO BEGIN & $

      READF, lun, sceneinfo_lines & $
        
        ; If count is -1 then we are at the first line of the file, AKA the header line.
        ; We will skip this iteration so the header doesn't get included in scene_info_entries.
        IF count EQ -1 THEN BEGIN
          count = count + 1
          CONTINUE
        ENDIF ELSE BEGIN
          ; Split the current line into individual entries using a space as a seperator, and add them to the array scene_info_entries.
          sceneinfo_entries[*,count] = STRSPLIT(sceneinfo_lines, ' ', /EXTRACT)
          count = count + 1
        ENDELSE
        
    ENDWHILE
    
  CLOSE, lun
  
  ; 
  ; For each modis scene, get the sst and quality band, georeference them, subset them, and save them as geotiffs.
  ;
  FOR i=0,N_ELEMENTS(sceneinfo_entries[0,*])-1 DO BEGIN
    
    ;
    ; Create results folder if it does not exist already
    ;
    output_folder  = 'georeferenced\'
    exists = FILE_TEST(modis_folderpath + output_folder)
    IF exists EQ 0 THEN BEGIN
      FILE_MKDIR, modis_folderpath + output_folder
    ENDIF

    ;
    ; Identify current scene and Landsat corner lat/lons (add or substract .1 to make sure whole Landsat scene 
    ; gets included when the MODIS images get subset later on)
    ;
    lsat_file   = sceneinfo_entries[0,i]
    modis_file  = sceneinfo_entries[1,i]
    utm_zone    = sceneinfo_entries[2,i]
    lsat_ul_lat = FLOAT(sceneinfo_entries[3,i]) + .1
    lsat_ul_lon = FLOAT(sceneinfo_entries[4,i]) - .1
    lsat_lr_lat = FLOAT(sceneinfo_entries[5,i]) - .1
    lsat_lr_lon = FLOAT(sceneinfo_entries[6,i]) + .1

    PRINT, 'lsat_scene = ', lsat_file
    PRINT, 'modis_scene = ', modis_file
    
    ;
    ; Define input and output paths (Currently named by adding _sst.tif or _qual.tif to the MODIS scene name)
    ;
    fname = modis_folderpath + modis_file + '.L2_LAC_SST.nc'
    output_path   = modis_folderpath + output_folder
    out_sst_tiff  = output_path + STRMID(modis_file,0,14) + '_sst.tif'
    out_qual_tiff = output_path + STRMID(modis_file,0,14) + '_qual.tif'
    
    ; Convert Landsat corners to UTM (output gives zone, Easting, and Northing for each lat/lon pair)
    upperleftUTM  = CONVERT_LL_UTM(lsat_ul_lat, lsat_ul_lon, utm_zone)
    lowerrightUTM = CONVERT_LL_UTM(lsat_lr_lat, lsat_lr_lon, utm_zone)
    
    ;
    ; Use NCDF functions to acquire SST and quality data 
    ;
    IF FILE_TEST(fname) EQ 0 THEN BEGIN
      PRINT, "SST file does not exist! Skipping ..."
      CONTINUE
    ENDIF ELSE BEGIN
      file_id = NCDF_OPEN(fname)
      group_ids = NCDF_GROUPSINQ(file_id)
      data_id = group_ids[2]
      nav_id = group_ids[3]
      datavar_ids = NCDF_VARIDSINQ( data_id )
      varinq = NCDF_VARINQ(nav_id, datavar_ids[1])
      ; SST and quality data lives in the data group, lat/lon data lives in the navigation group
      NCDF_VARGET, data_id, 0, sst_data
      NCDF_VARGET, data_id, 1, qual_data
      NCDF_VARGET, nav_id,  1, lat_data
      NCDF_VARGET, nav_id,  0, lon_data
    ENDELSE
    
    ;
    ; Stack SST, quality, latitude, and longitude arrays and send to ENVI
    ;
    s = SIZE(sst_data)
    super_data = MAKE_ARRAY(s[1],s[2],4)
    super_data[*,*,0] = sst_data
    super_data[*,*,1] = qual_data
    super_data[*,*,2] = lat_data
    super_data[*,*,3] = lon_data
    
    ENVI_ENTER_DATA, super_data*1.0, R_FID=super_id
    
    
    ; ==========================
    ; STEP 2: GEOREFERNCE IMAGES
    ; ==========================
  
    ;
    ; Define input and output projections, and check if landsat scene is in north or south hemisphere  
    ;
    iproj = envi_proj_create(/geographic, datum='WGS-84')
    IF lsat_ul_lat GT 0 THEN BEGIN
      oproj = envi_proj_create(/utm, zone=utm_zone, datum='WGS-84')
    ENDIF ELSE BEGIN
      oproj = envi_proj_create(/utm, zone=utm_zone, /south, datum='WGS-84')
    ENDELSE
    
    ; Build GLT file from latitude (Y_POS) and longitude (X_POS) arrays
    ENVI_DOIT, 'ENVI_GLT_DOIT', DIMS=dims, I_PROJ=iproj, /IN_MEMORY, O_PROJ=oproj, $ 
      OUT_NAME='test', PIXEL_SIZE=1000.0, R_FID=glt_id, ROTATION=0, X_FID=super_id, $ 
      X_POS=3, Y_FID=super_id, Y_POS=2
      
    ; Georeference SST image using GLT file (band 1 --> POS=0)
    ENVI_DOIT, 'ENVI_GEOREF_FROM_GLT_DOIT', BACKGROUND=-1.0, FID=super_id, $
      GLT_FID=glt_id, /IN_MEMORY, OUT_NAME='georef_sst', POS=0, R_FID=georef_sst_id
        
    ; Georeference Quality image using GLT file (band 2 --> POS=1)
    ENVI_DOIT, 'ENVI_GEOREF_FROM_GLT_DOIT', BACKGROUND=-1.0, FID=super_id, $
      GLT_FID=glt_id, /IN_MEMORY, OUT_NAME='georef_qual', POS=1, R_FID=georef_qual_id
    


    ; ======================================================================
    ; STEP 3: FIND UPPER LEFT AND LOWER RIGHT LANDSAT CORNERS IN MODIS IMAGE
    ; ======================================================================
    
    ;
    ; Get image size info of the georeferenced images
    ;
    ENVI_FILE_QUERY, georef_qual_id, dims=qualdims
    new_image = ENVI_GET_DATA(DIMS=qualdims, FID=georef_qual_id, POS=0)
    new_image_size = SIZE(new_image)
    image_num_cols = new_image_size[1]
    image_num_rows = new_image_size[2]

    ;
    ; Convert Landsat corner UTMS to pixel locations within the georeferenced MODIS scene
    ;
    xutm=[upperleftUTM[1], lowerrightUTM[1]]
    yutm=[upperleftUTM[2], lowerrightUTM[2]]    
    ENVI_CONVERT_FILE_COORDINATES, georef_sst_id, xpix, ypix, xutm, yutm
    IF xpix[0] LT 0 THEN xpix[0] = 0
    IF ypix[0] LT 0 THEN ypix[0] = 0
    ; If pixel location is outside image, then set location to last row or column
    IF xpix[1] GT image_num_cols-1 THEN xpix[1] = image_num_cols-1
    IF ypix[1] GT image_num_rows-1 THEN ypix[1] = image_num_rows-1

    ;
    ; Define subset array, which must be a 5 element array of long type
    ;
    subset_array=MAKE_ARRAY(5,1,/LONG)
    subset_array[0]= -1L
    subset_array[1]= FLOOR(xpix[0])
    subset_array[2]= FLOOR(xpix[1])
    subset_array[3]= FLOOR(ypix[0])
    subset_array[4]= FLOOR(ypix[1])
    
    ; ==================================================
    ; STEP 4: SUBSET ARRAYS AND WRITE IMAGES TO GEOTIFFS 
    ; ==================================================
    
    ; Subsetting the image and saving it is done in the same step
    ENVI_OUTPUT_TO_EXTERNAL_FORMAT, FID=georef_sst_id, DIMS=subset_array, POS=0, OUT_NAME=out_sst_tiff, /TIFF
    ENVI_OUTPUT_TO_EXTERNAL_FORMAT, FID=georef_qual_id, DIMS=subset_array, POS=0, OUT_NAME=out_qual_tiff, /TIFF

  ENDFOR
END