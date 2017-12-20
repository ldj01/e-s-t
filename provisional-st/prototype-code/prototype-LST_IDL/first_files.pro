; Monica Cook
; 6 March 2012
;
;
; NAME:  FIRST_FILES
;  
; PURPOSE:  IDL PROCEDURE
;   Creates directories and writes tape5 file, caseList, and commandList
;
; CALL SEQUENCE: FIRST_FILES, home, $            ;directory containing programs and supporting files
;                             directory, $       ;directory containing landsat metadata and location for results
;                             imageBase, $       ;basename of landsat files
;                             year, $            ;year of aquisition of landsat scene
;                             month, $           ;month of acquisition of landsat scene
;                             day, $             ;day of acquisition of landsat scene 
;                             hour, $            ;hour of acquisition of landsat scene
;                             min, $             ;minute of acquisition of landsat scene
;                             sec, $             ;second of acquisition of landsat scene
;                             landsatSamples, $  ;number of columns in thermal landsat scene
;                             landsatLines, $    ;number of rows in thermal landsat scene
;                             pixelSize, $       ;size of pixel in thermal landsat scene (in meters)
;                             zone, $            ;utm zone of landsat scene (for utm corrdinates)
;                             UL_EAST, $         ;easting coordinate of upper left corner
;                             UL_NORTH, $        ;northing coordiante of upper left corner
;                             UR_EAST, $         ;easting coordinate of upper right corner
;                             UR_NORTH, $        ;northing coordinate of upper right corner
;                             LL_EAST, $         ;easting coordinate of lower left corner
;                             LL_NORTH, $        ;northing coordinate of lower left corner
;                             LR_EAST, $         ;easting coordinate of lower right corner
;                             LR_NORTH,$         ;northing coordinate of lower right corner
;                             UL_LAT, $          ;latitude of upper left corner
;                             UL_LON, $          ;longitude of upper left corner
;                             UR_LAT, $          ;latitude of upper right corner
;                             UR_LON, $          ;longitude of upper right corner
;                             LL_LAT, $          ;latitude of lower left corner
;                             LL_LON, $          ;longitude of lower left corner
;                             LR_LAT, $          ;latitude of lower right corner
;                             LR_LON             ;longitude of lower right corner
;    
;  RESTRICTIONS:
;    Called from faster.bash
;  
;  REQUIRED PROGRAMS AND FILES (in home directory):   coordinates.txt
;                                                     CONVERT_GEOMETRIC_GEOPOTENTIAL.pro
;                                                     CONVERT_SH_RH.pro
;                                                     stanAtm.txt
;                                                     head.txt
;                                                     tail.txt
;  

PRO FIRST_FILES, home, $
                 directory, $
                 imageBase, $
                 year, $
                 month, $
                 day, $
                 hour, $
                 min, $
                 sec, $
                 landsatSamples, $
                 landsatLines, $
                 pixelSize, $
                 zone, $
                 UL_EAST, $
                 UL_NORTH, $
                 UR_EAST, $
                 UR_NORTH, $
                 LL_EAST, $
                 LL_NORTH, $
                 LR_EAST, $
                 LR_NORTH, $
                 UL_LAT, $
                 UL_LON, $
                 UR_LAT, $
                 UR_LON, $
                 LL_LAT, $
                 LL_LON, $
                 LR_LAT, $
                 LR_LON

   ;convert inputs
   home = STRING(home)
   directory = STRING(directory)
   imageBase = STRING(imageBase)
   year = DOUBLE(year)
   month = DOUBLE(month)
   day = DOUBLE(day)
   hour = DOUBLE(hour)
   min = DOUBLE(min)
   sec = DOUBLE(sec)
   landsatSamples = DOUBLE(landsatSamples)
   landsatLines = DOUBLE(landsatLines)
   pixelSize = DOUBLE(pixelSize)
   zone = DOUBLE(zone)
   UL_EAST = DOUBLE(UL_EAST)
   UL_NORTH = DOUBLE(UL_NORTH)
   UR_EAST = DOUBLE(UR_EAST)
   UR_NORTH = DOUBLE(UR_NORTH)
   LL_EAST = DOUBLE(LL_EAST)
   LL_NORTH = DOUBLE(LL_NORTH)
   LR_EAST = DOUBLE(LR_EAST)
   LR_NORTH = DOUBLE(LR_NORTH)
   UL_LAT = DOUBLE(UL_LAT)
   UL_LON = DOUBLE(UL_LON)
   UR_LAT = DOUBLE(UR_LAT)
   UR_LON = DOUBLE(UR_LON)
   LL_LAT = DOUBLE(LL_LAT)
   LL_LON = DOUBLE(LL_LON)
   LR_LAT = DOUBLE(LR_LAT)
   LR_LON = DOUBLE(LR_LON)
   
   
;##### Read in GRIB data and pull out points pertinent to current Landsat scene
   
   ;read in file containing latitude and longitude grid information for the grib data
   OPENR, 20, home+'coordinates.txt'
   coordinates = MAKE_ARRAY(4, 96673, /DOUBLE)
   READF, 20, coordinates
   CLOSE, 20
   FREE_LUN, 20
   
   ;pull out latitude and reform to 349x277 grid
   narrLat = coordinates[2,*]
   lat = REFORM(narrLat, 349, 277)
   
   ;pull out longitude, manipulate range to [-180,180] and reform to 349x277 grid
   narrLon = coordinates[3,*]
   east = WHERE(narrLon GT 180.0)
   narrLon[east] = 360.0 - narrLon[east]
   west = WHERE(narrLon LE 180.0)
   narrLon[west] = (-1)*narrLon[west]   
   lon = REFORM(narrLon, 349, 277)
   
   ;pull out i and j values and reform to 349 x 277 grid
   i = coordinates[0,*]
   eye = REFORM(i, 349, 277)
   j = coordinates[1,*]
   jay = REFORM(j, 349, 277)
   
   ;define 29 pressure levels in grib data
   p = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, $
        600, 550, 500, 450, 400, 350, 300, 275, 250, 225, 200, 175, 150, 125, 100]
        
   ;create arrays for height, specific humidity and temperature for time before Landsat acquisition
   hgt_1 = MAKE_ARRAY(96673, 29, /DOUBLE)
   shum_1 = MAKE_ARRAY(96673, 29, /DOUBLE)
   tmp_1 = MAKE_ARRAY(96673, 29, /DOUBLE)
   
   ;read in height, specific humidity, and temperature from text files for time before Landsat acquisition     
   FOR i = 0, N_ELEMENTS(p)-1 DO BEGIN
   
      file = STRMID(STRCOMPRESS(STRING(p[i])), 1)+'.txt'
      hgtFile = directory+'HGT_1/'+file
      shumFile = directory+'SHUM_1/'+file
      tmpFile = directory+'TMP_1/'+file
      
      OPENR, 20, hgtFile
      tempHGT = MAKE_ARRAY(96675,/DOUBLE)
      READF, 20, tempHGT
      CLOSE, 20
      FREE_LUN, 20
      hgt_1[*,i] = tempHGT[2:96674]
      
      OPENR, 20, shumFile
      tempSHUM = MAKE_ARRAY(96675,/DOUBLE)      
      READF, 20, tempSHUM
      CLOSE, 20
      FREE_LUN, 20
      shum_1[*,i] = tempSHUM[2:96674]
      
      OPENR, 20, tmpFile
      tempTMP = MAKE_ARRAY(96675,/DOUBLE)            
      READF, 20, tempTMP
      CLOSE, 20
      FREE_LUN, 20
      tmp_1[*,i] = tempTMP[2:96674]
               
   ENDFOR
   
   ;create arrays for height, specific humidity and temperature for time after Landsat acquisition
   hgt_2 = MAKE_ARRAY(96673, 29, /DOUBLE)
   shum_2 = MAKE_ARRAY(96673, 29, /DOUBLE)
   tmp_2 = MAKE_ARRAY(96673, 29, /DOUBLE)
   
   ;read in height, specific humidity, and temperature from text files for time after Landsat acquisition     
   FOR i = 0, N_ELEMENTS(p)-1 DO BEGIN
   
      file = STRMID(STRCOMPRESS(STRING(p[i])), 1)+'.txt'
      hgtFile = directory+'HGT_2/'+file
      shumFile = directory+'SHUM_2/'+file
      tmpFile = directory+'TMP_2/'+file
      
      OPENR, 20, hgtFile
      tempHGT = MAKE_ARRAY(96675,/DOUBLE)
      READF, 20, tempHGT
      CLOSE, 20
      FREE_LUN, 20
      hgt_2[*,i] = tempHGT[2:96674]
      
      OPENR, 20, shumFile
      tempSHUM = MAKE_ARRAY(96675,/DOUBLE)      
      READF, 20, tempSHUM
      CLOSE, 20
      FREE_LUN, 20
      shum_2[*,i] = tempSHUM[2:96674]
      
      OPENR, 20, tmpFile
      tempTMP = MAKE_ARRAY(96675,/DOUBLE)            
      READF, 20, tempTMP
      CLOSE, 20
      FREE_LUN, 20
      tmp_2[*,i] = tempTMP[2:96674]
               
   ENDFOR
   
   ;determine if landsat is in the northern or southern hemisphere. '6' = northern hemisphere, '7' = southern hermisphere.
   IF UL_LAT GE 0 THEN landsatHemi = 6 ELSE landsatHemi = 7
   
   ;expand range to include NARR points outside image for edge pixels
   UL_LAT = UL_LAT + 0.5
   UL_LON = UL_LON - 0.5
   LR_LAT = LR_LAT - 0.5
   LR_LON = LR_LON + 0.5
      
   ;determine what points in the NARR dataset fall within the Landsat image using logical operators
   ;lessThanLat and greaterThanLat are values where the NARR values are less than or greater than the edges of the Landsat
   ;corners values respectively
   ;pixels that are true in both fall within the Landsat scene
   ;the same thing is done with longitude values
   lessThanLat = DOUBLE(lat LT UL_LAT)
   greaterThanLat = DOUBLE(lat GT LR_LAT)
   keepLat = lessThanLat*greaterThanLat
   greaterThanLon = DOUBLE(lon GT UL_LON)
   lessThanLon = DOUBLE(lon LT LR_LON)
   keepLon = greaterThanLon*lessThanLon
   
   ;values that are true in both keepLat and keepLon fall within the Landsat image
   ;convert indices into (x,y) values in the NARR dataset
   ;Because the Landsat is in easting/northing originally and the NARR data is Lambert conformal, this is not a perfect square
   keep = keepLat * keepLon
   inLandsat = WHERE(keep NE 0)
   
   ;determine indices to pull out rectangle of NARR points
   iindices = [MIN(eye[inLandsat])-1,MAX(eye[inLandsat])-1]
   jindices = [MIN(jay[inLandsat])-1,MAX(jay[inLandsat])-1]
   
   ;extract coordinates within this rectangle
   ivalues = eye[iindices[0]:iindices[1],jindices[0]:jindices[1]]
   jvalues = jay[iindices[0]:iindices[1],jindices[0]:jindices[1]]
   latvalues = lat[iindices[0]:iindices[1],jindices[0]:jindices[1]]
   lonvalues = lon[iindices[0]:iindices[1],jindices[0]:jindices[1]]

   ;determine indices of included NARR points within full array
   NARRindices = MAKE_ARRAY(N_ELEMENTS(ivalues), /DOUBLE)
   ivector = REFORM(ivalues,N_ELEMENTS(ivalues))
   jvector = REFORM(jvalues,N_ELEMENTS(jvalues))

   FOR i = 0, N_ELEMENTS(NARRindices)-1 DO BEGIN
      NARRindices[i] = WHERE(eye EQ ivector[i] AND jay EQ jvector[i])
   ENDFOR
      
   ;determine how many NARR points fall within the landsat scene
   numPoints = N_ELEMENTS(ivalues)  
 
   ;create pressure array [numPoints x 29 pressure levels]
   pressure = MAKE_ARRAY(numPoints, 29, /DOUBLE)
   FOR i = 0, numPoints-1 DO pressure[i,*] = p
   
   ;define arrays of height, temperature and specific humidity values that are [numPoints x 29 pressure levels]
   ;for times before and after Landsat acquisition
   ;these are the NARR values pertinent to the current landsat scene
   landsatHGT_1 = hgt_1[NARRindices,*]
   landsatTMP_1 = tmp_1[NARRindices,*]
   landsatSHUM_1 = shum_1[NARRindices,*]
   landsatHGT_2 = hgt_2[NARRindices,*]
   landsatTMP_2 = tmp_2[NARRindices,*]
   landsatSHUM_2 = shum_2[NARRindices,*]
   
   ;convert grib data to variables to be input to MODTRAN
   geometricHeight_1 = CONVERT_GEOPOTENTIAL_GEOMETRIC(landsatHGT_1, lat[NARRindices])
   geometricHeight_1 = geometricHeight_1/1000D
   geometricHeight_2 = CONVERT_GEOPOTENTIAL_GEOMETRIC(landsatHGT_2, lat[NARRindices])
   geometricHeight_2 = geometricHeight_2/1000D
   relativeHumidity_1 = CONVERT_SH_RH(landsatSHUM_1, landsatTMP_1, pressure)
   relativeHumidity_2 = CONVERT_SH_RH(landsatSHUM_2, landsatTMP_2, pressure)
   temperature_1 = landsatTMP_1
   temperature_2 = landsatTMP_2
   
;##### Temporal interpolation   
   
   ;determine three hour-increment before and after scene center scan time
   rem1 = hour MOD 3
   rem2 = 3 - rem1
   hour1 = DOUBLE(hour - rem1)
   hour2 = DOUBLE(hour + rem2)
      
   ;round to nearest minute
   IF sec GE 30 THEN min = min + 1
   
   ;convert hour-min acquisition time to decimal time
   time = hour + (min/60D)
   
   ;linearly interpolate geometric height, relative humidity, and temperature for NARR points within
   ;Landsat scene
   ;this is the NARR data corresponding to the acquisition time of the Landsat image converted to appropriated
   ;variable for MODTRAN input
   geometricHeight = geometricHeight_1 + (time-hour1)*((geometricHeight_2 - geometricHeight_1)/(hour2 - hour1))
   relativeHumidity = relativeHumidity_1 + (time-hour1)*((relativeHumidity_2 - relativeHumidity_1)/(hour2 - hour1))
   temperature = temperature_1 + (time-hour1)*((temperature_2 - temperature_1)/(hour2 - hour1))
      
;##### Building tape5 files

   ;read in file containing standard mid lat summer atmosphere information to be used for upper layers
   OPENR, 20, home+'stanAtm.txt'
   stanAtm = MAKE_ARRAY(4, 30, /DOUBLE)
   READF, 20, stanAtm
   CLOSE, 20
   FREE_LUN, 20
   
   ;separate variable in standard atmosphere
   stanGeoHeight = stanAtm[0,*]
   stanPress = stanAtm[1,*]
   stanTemp = stanAtm[2,*]
   stanRelHum = stanAtm[3,*]
   
   ;determine index of last layer in standard atmosphere
   stanNum = N_ELEMENTS(stanGeoHeight)
   stanLast = stanNum - 1
   
   ;define an array of ground altitudes at which to run modtran for each NARR point
   gndalt = [0.0, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 3.6, 4.05]
   numElevations = N_ELEMENTS(gndalt)
   
   ;determine number of MODTRAN runs
   numCases = numPoints*numElevations*3
   ;initialize arrays
   caseList = MAKE_ARRAY(numCases, /STRING)
   commandList = MAKE_ARRAY(numCases, /STRING)
   entry = 0

   ;iterate through all NARR points within the landsat scene
   FOR i = 0, numPoints-1 DO BEGIN
   PRINT, i

      ;create a directory for the current NARR point
      latString = STRMID(STRCOMPRESS(STRING(narrLat[NARRindices[i]])),1,6)
      IF narrLon[NARRindices[i]] LT 0 THEN BEGIN
         lonString = STRMID(STRCOMPRESS(STRING(narrLon[NARRindices[i]])),2,6)
      ENDIF ELSE BEGIN
         tempLon = 360D - narrLon[NARRindices[i]]
         lonString = STRMID(STRCOMPRESS(STRING(tempLon)),1,6)
      ENDELSE
      currentPoint = directory+latString+'_'+lonString
      command = "mkdir "+currentPoint
      SPAWN, command      
         
      ;set lowest altitude is the first geometric height at that NARR point (if positive)
      ; (if negative set to zero)
      IF geometricHeight[i,0] LT 0 THEN gndalt[0] = 0.000 ELSE gndalt[0] = geometricHeight[i,0]
   
      ;pull out arrays for current NARR point
      p = REFORM(pressure[i,*])
      t = REFORM(temperature[i,*])
      hgt = REFORM(geometricHeight[i,*])
      rh = REFORM(relativeHumidity[i,*])
      
      ;define number of levels and index of maximum level in NARR data
      numLevels = N_ELEMENTS(p)
      maxLevel = numLevels-1
      
      ;determine latitude and longitude of current NARR point and insert into tail file
      command = "cat "+home+"tail.txt | sed 's/latitu/"+latString+"/' > "+home+directory+"newTail.txt"
      SPAWN, command
      command = "cat "+home+directory+"newTail.txt | sed 's/longit/"+lonString+"/' > "+home+directory+"newTail2.txt"
      SPAWN, command

      ;determine current julian day
      ;first determine if current year is a leap year
      ;A year is a leap year if it is divisible by 4 but not if its divisible by 100 except when its divisible by 400
      leap = 0
      IF (year MOD 4) EQ 0 THEN BEGIN
         leap = 1
         IF (year MOD 100) EQ 0 THEN BEGIN
            IF (year MOD 400) EQ 0 THEN leap = 1 ELSE leap = 0
         ENDIF
      ENDIF
      CASE month OF
         1: JDAY = (day)
         2: JDAY = 31 + (day)
         3: IF leap THEN JDAY = 31 + 29 + (day) ELSE JDAY = 31 + 28 + (day)
         4: IF leap THEN JDAY = 91 + (day) ELSE JDAY = 90 + (day)
         5: IF leap THEN JDAY = 121 + (day) ELSE JDAY = 120 + (day)
         6: IF leap THEN JDAY = 152 + (day) ELSE JDAY = 151 + (day)
         7: IF leap THEN JDAY = 182 + (day) ELSE JDAY = 181 + (day)
         8: IF leap THEN JDAY = 213 + (day) ELSE JDAY = 212 + (day)
         9: IF leap THEN JDAY = 244 + (day) ELSE JDAY = 243 + (day)
         10: IF leap THEN JDAY = 274 + (day) ELSE JDAY = 273 + (day)
         11: IF leap THEN JDAY = 305 + (day) ELSE JDAY = 304 + (day)
         12: IF leap THEN JDAY = 335 + (day) ELSE JDAY = 334 + (day)
      ENDCASE

      ;insert current julian day into tail file
      JDAY = FIX(JDAY)         
      IF JDAY GE 100 THEN BEGIN
         jay = STRMID(STRCOMPRESS(STRING(JDAY)),1)
      ENDIF ELSE BEGIN
         jay = STRCOMPRESS(STRING(JDAY))
      ENDELSE
      command = "cat "+home+directory+"newTail2.txt | sed 's/jay/"+jay+"/' > "+home+directory+"newTail3.txt"
      SPAWN, command      
      
      ;iterature through all ground altitudes at which modtran is run
      FOR j = 0, numElevations-1 DO BEGIN
            
         ;create a directory for the current height
         gdalt = STRING(gndalt[j], FORMAT = '(F5.3)')                  
         currentGdalt = currentPoint+'/'+gdalt
         command = "mkdir "+currentGdalt
         SPAWN, command

         ;print current location and height to command line         
         PRINT, latString, ' ', lonString, ' ', gdalt
          
         ;determine layers below current gndalt and closest index above and below
         delete = WHERE(hgt LT gndalt[j])
         indexBelow = N_ELEMENTS(delete)-1
         indexAbove = N_ELEMENTS(delete)

print,"j,indexBelow,indexAbove=",j,indexBelow,indexAbove
         
         ;linearly interpolate pressure, temperature, and relative humidity to gndalt for lowest layer
         newPressure = p[indexBelow]+(gndalt[j]-hgt[indexBelow])*((p[indexAbove]-p[indexBelow])/(hgt[indexAbove]-hgt[indexBelow]))        
         newTemperature = t[indexbelow]+(gndalt[j]-hgt[indexBelow])*((t[indexAbove]-t[indexBelow])/(hgt[indexAbove]-hgt[indexBelow]))        
         newRelativeHumidity = rh[indexBelow]+(gndalt[j]-hgt[indexBelow])*((rh[indexAbove]-rh[indexBelow])/(hgt[indexAbove]-hgt[indexBelow]))
         
         ;create arrays containing only layers to be included in current tape5 file
         tempGeoHeight = [gndalt[j], [hgt[indexAbove:maxLevel]]]         
         tempPress = [newPressure, [p[indexAbove:maxLevel]]]
         tempTemp = [newTemperature, [t[indexAbove:maxLevel]]]
         tempRelHum = [newRelativeHumidity, [rh[indexAbove:maxLevel]]]
         
         ;modtran throws as error when there are two identical layers in the tape5 file
         ;if the current ground altitude and the next highest layer are close enough, eliminate interpolated layer
         IF ABS(gndalt[j]-hgt[indexAbove]) LT 0.001 THEN BEGIN
            tempGeoHeight = hgt[indexAbove:maxLevel]
            tempPress = p[indexAbove:maxLevel]
            tempTemp = t[indexAbove:maxLevel]
            tempRelHum = rh[indexAbove:maxLevel]
         ENDIF
         
         ;determine index of last layer
         last = N_ELEMENTS(tempGeoHeight)-1
         
         ;determine maximum height of NARR layers and where the standard atmosphere is greater than this
         maxHeight = hgt[maxLevel]
         above = WHERE(stanGeoHeight GT maxHeight)
         ;if there are at least three layers above to highest NARR layer, add standard atmosphere layers
         IF N_ELEMENTS(above) GE 3 THEN BEGIN
         
            ;interpolate to 2 layers above highest NARR layer to create a smooth transition
            interpolateTo = above[2]
                  
            ;linearly interpolate height, pressure, temp, and rel hum to create a smooth transition between
            ;the NARR layers and the standard upper atmosphere
            newHeight = (stanGeoHeight[interpolateTo]+tempGeoHeight[last])/2D 
            newPressure2 = tempPress[last]+(newHeight-tempGeoHeight[last])*((stanPress[interpolateTo]-tempPress[last])/(stanGeoHeight[interpolateTo]-tempGeoHeight[last]))
            newTemperature2 = tempTemp[last]+(newHeight-tempGeoHeight[last])*((stanTemp[interpolateTo]-tempTemp[last])/(stanGeoHeight[interpolateTo]-tempGeoHeight[last]))            
            newRelativeHumidity2 = tempRelHum[last]+(newHeight-tempGeoHeight[last])*((stanRelHum[interpolateTo]-tempRelHum[last])/(stanGeoHeight[interpolateTo]-tempGeoHeight[last]))

            ;concatenate NARR layers, new layer, and standard atmosphere layers
            tempGeoHeight = [tempGeoHeight, newHeight, stanGeoHeight[interpolateTo:stanLast]]         
            tempPress = [tempPress, newPressure2, stanPress[interpolateTo:stanLast]]
            tempTemp = [tempTemp, newTemperature2, stanTemp[interpolateTo:stanLast]]
            tempRelHum = [tempRelHum, newRelativeHumidity2, stanRelHum[interpolateTo:stanLast]]
            
            ;determine index of last layer
            last = N_ELEMENTS(tempGeoHeight)-1 

         ENDIF
         
         ;write atmospheric layers to a text file in format proper for tape5 file
         OPENW, unit, home+directory+'tempLayers.txt', /GET_LUN
         FOR k = 0, last DO BEGIN
            PRINTF, unit, tempGeoHeight[k], tempPress[k], tempTemp[k], tempRelHum[k], 0, 0, 'AAH             ', $
                          FORMAT = '(F10.3,E10.3,E10.3,E10.3,E10.3,E10.3,A16)'
         ENDFOR
         CLOSE, unit
         FREE_LUN, unit
                           
         ;determine number of layers for current ground altitude and insert into head file
         numLayers = N_ELEMENTS(tempGeoHeight)
         IF numLayers GE 100 THEN BEGIN
            nml = STRMID(STRCOMPRESS(STRING(numLayers)),1)
         ENDIF ELSE BEGIN
            nml = STRCOMPRESS(STRING(numLayers))
         ENDELSE
         command = "cat "+home+"head.txt | sed 's/nml/"+nml+"/' > "+home+directory+"newHead.txt"
         SPAWN, command
         
         ;insert current ground altitude into head file
         command = "cat "+home+directory+"newHead.txt | sed 's/gdalt/"+gdalt+"/' > "+home+directory+"newHead2.txt"
         SPAWN, command
         
         ;define arrays containing [temperature,albedo]  pairs at which to run modtran
         tmp = ['273','310','000']
         alb = ['0.0','0.0','0.1']         
                  
         ;iterate through pairs
         FOR k = 0, 2 DO BEGIN
            
            ;create directory for the current temperature
            currentTemp = currentGdalt+'/'+tmp[k]
            command = "mkdir "+currentTemp
            SPAWN, command
            
            ;insert current temperature into head file
            command = "cat "+home+directory+"newHead2.txt | sed 's/tmp/"+tmp[k]+"/' > "+home+directory+"newHead3.txt"
            SPAWN, command

            ;create directory for the current albedo
            currentAlb = currentTemp+'/'+alb[k]
            command = "mkdir "+currentAlb
            SPAWN, command            
            
            ;insert current albedo into head file
            command = "cat "+home+directory+"newHead3.txt | sed 's/alb/"+alb[k]+"/' > "+home+directory+"newHead4.txt"
            SPAWN, command            
                  
            ;concatenate head file, atmospheric layers, and tail file to create a tape5 file for modtran
            ;specific to this location and ground altitude with variables for temperature and albedo
            headFile = home+directory+'newHead4.txt '
            tailFile = home+directory+'newTail3.txt'
            tempLayers = home+directory+'tempLayers.txt '
            newFile = currentAlb+'/tape5'
            command = 'cat '+headFile+tempLayers+tailFile+' > '+newFile
            SPAWN, command
               
            ;create string for case list containing location of current tape5 file
            ;create string for command list containing commands for modtran run
            ;iterate entry count
            currentCase = currentAlb
            caseList[entry] = home+currentCase
            commandList[entry] = 'pushd ' + home + currentCase + '; ln -s /home/sguo/MODTRAN/DATA; /home/sguo/MODTRAN/Mod90_5.2.2.exe; popd'
            entry = entry+1

         ;end temperatures loop   
         ENDFOR      
      ;end ground altitudes loop         
      ENDFOR  
   ;end narr points loop
   ENDFOR

   ;write caseList to a file
   OPENW, unit, home+directory+'caseList', /GET_LUN
   FOR k = 0, N_ELEMENTS(WHERE(caseList NE ''))-1 DO BEGIN
      PRINTF, unit, caseList[k]
   ENDFOR
   CLOSE, unit
   FREE_LUN, unit
   
   ;write commandList to a file
   OPENW, unit, home+directory+'commandList', /GET_LUN
   FOR k = 0, N_ELEMENTS(WHERE(commandList NE ''))-1 DO BEGIN
      PRINTF, unit, commandList[k]
   ENDFOR
   CLOSE, unit
   FREE_LUN, unit   

END
