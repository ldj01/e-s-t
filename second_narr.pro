; Monica Cook
; 2 February 2012
;
;
; NAME:  SECOND_NARR
;  
; PURPOSE:  IDL PROCEDURE
;   Generate transmission, upwelled radiance, and downwelled radiance at each height at each NARR point
;
; CALL SEQUENCE: SECOND_NARR, home, $            ;directory containing programs and supporting files
;                             directory, $       ;directory containing landsat metadata and location for results
;                             numPoints, $       ;number of NARR points within Landsat scene
;                             numHeights, $      ;number of height at each NARR point
;                             alb, $             ;albedo of second MODTRAN run at each height at each NARR point
;			      whichLandsat	 ;which Landsat sensor this image is captured with
;    
;  RESTRICTIONS:
;    Called from faster.bash
;    Already called GET_FILES_FASTER.pro, commandList, and parsed tape6 files
;  
;  REQUIRED PROGRAMS AND FILES (in home directory):   L5v2.rsp
;                                                     caseList
;                                                     CALCULATE_LT.pro
;                                                     CALCULATE_LOBS.pro
;
;
;
;

PRO SECOND_NARR, home, $
                 directory, $
                 numPoints, $
                 numHeights, $
                 alb, $
	         whichLandsat

   ;convert inputs
   home = STRING(home)
   directory = STRING(directory)
   numPoints = FIX(numPoints)
   numHeights = FIX(numHeights)
   alb = DOUBLE(alb)
   whichLandsat = FIX(whichLandsat)
   
   ;define emissivity from albedo
   ems = 1-alb

   CASE whichLandsat OF
      5: BEGIN
         OPENR, 20, home+'L5v2.rsp'
         spectralResponse = MAKE_ARRAY(2, 171, /DOUBLE)
         READF, 20, spectralResponse
         CLOSE, 20
         FREE_LUN, 20
      END
      7: BEGIN
         OPENR, 20, home+'L7.rsp'
         spectralResponse = MAKE_ARRAY(2, 47, /DOUBLE)
         READF, 20, spectralResponse
         CLOSE, 20
         FREE_LUN, 20
      END
      10: BEGIN
         OPENR, 20, home+'L8_B10.rsp'
         spectralResponse = MAKE_ARRAY(2, 101, /DOUBLE)
         READF, 20, spectralResponse
         CLOSE, 20
         FREE_LUN, 20
      END
      11: BEGIN
         OPENR, 20, home+'L8_B11.rsp'
         spectralResponse = MAKE_ARRAY(2, 101, /DOUBLE)
         READF, 20, spectralResponse
         CLOSE, 20
         FREE_LUN, 20
      END
   ENDCASE

   ;read caseList into string array
   caseListFile = home+directory+'caseList'
   ;determine number of entries in caseList file
   command = "wc "+caseListFile+" | awk '{print $1}'"
   SPAWN, command, numEntries
   ;numEntries = FIX(numEntries)   
   numEntries = FIX(numEntries[0])   
    
  OPENR, 20, caseListFile
   caseList = MAKE_ARRAY(1, numEntries, /STRING)
   READF, 20, caseList
   CLOSE, 20
   FREE_LUN, 20
   
   ;calculate Lt for each temperature
   tempRadiance273 = CALCULATE_LT(273, spectralResponse)
   tempRadiance310 = CALCULATE_LT(310, spectralResponse)
 
   ;initialize counter
   ;counter is where to extract from caseList
   counter = 0      
   
   ;initialize place
   ;place is where to place entries in results
   place = 0
   
   ;create array to contain calculated output
   results = MAKE_ARRAY(6, numPoints*numHeights, /STRING)
   
   ;iterate through all points in the scene
   FOR i = 0, numPoints-1 DO BEGIN
     
      ;iterate through all heights at each points
      FOR j = 0, numHeights-1 DO BEGIN
      
         ;determine current latlon and height
         ;depends on number of steps in path
         extract = caseList[counter]
         command = "echo "+extract+" | tr '/' '\n'"
         SPAWN, command, specs
         lat_lon = specs[7]
         height = specs[8]
         command = "echo "+lat_lon+" | tr '_' '\n'"
         SPAWN, command, coordinates
         lat = coordinates[0]
         lon = coordinates[1]
         
         ;determine number of entries in current file
         command = "wc "+extract+"/parsed | awk '{print $1}'"
         SPAWN, command, numEntries
         numEntries = FIX(numEntries[0])

print,"i,j,numEntries =",i,j,numEntries
               
         ;for each height, read in radiance inforomation for three modtran runs
         ;columns of array are organized:
         ;wavelength | 273,0.0 | 310,0.0 | 000,0.1
         currentData = MAKE_ARRAY(4, numEntries, /DOUBLE)
         
         ;initialize index
         ;index is where to put data from parsed file into current data
         index = 0
      
         ;iterature through three pairs of parameters
         FOR k = 0, 2 DO BEGIN
            
            ;define current file
            currentFile = caseList[counter]+'/parsed'

            OPENR, 20, currentFile
            temp = MAKE_ARRAY(2, numEntries, /DOUBLE)
            READF, 20, temp
            CLOSE, 20
            FREE_LUN, 20
               
            ;put arrays into data array for current point at current height
            IF index EQ 0 THEN BEGIN
               currentData[0,*] = temp[0,*]
               index = index + 1
            ENDIF
               
            currentData[index,*] = temp[1,*]
            
            IF k EQ 2 THEN BEGIN
               ;determine temperature at lowest atmospheric layer (when MODTRAN is run at 0K)
               zeroTape6 = caseList[counter]+'/tape6'
               command = 'grep "TARGET-PIXEL (H2) SURFACE TEMPERATURE" '+ zeroTape6 + " | awk '{print $7}'"
               SPAWN, command, zeroTemp
               zeroTemp = DOUBLE(zeroTemp[0])
            ENDIF               
            
            ;iterate counter and index
            counter = counter+1
            index = index+1 
                                                   
         ENDFOR
                                       
         ;extract wavelengths from array of current data
         wavelengths = currentData[0,*]
         
         ;***parameters from 3 modtran runs
                                   
         ; Lobs = Lt*tau + Lu
         ; m = tau
         ; b = Lu
         x = MAKE_ARRAY(2,2)
         x[0,*] = 1
         x[1,0] = tempRadiance273
         x[1,1] = tempRadiance310
         y = MAKE_ARRAY(1,2)
         y[0,0] = CALCULATE_LOBS(wavelengths, currentData[1,*], spectralResponse)
         y[0,1] = CALCULATE_LOBS(wavelengths, currentData[2,*], spectralResponse)
         a = INVERT(TRANSPOSE(x)##x)##TRANSPOSE(x)##y
         tau = a[1]
         Lu = a[0]
      
         ;determine Lobs and Lt when modtran was run a 0 K - calculate downwelled
         tempRadianceZero = CALCULATE_LT(zeroTemp, spectralResponse)
         obsRadianceZero = CALCULATE_LOBS(wavelengths, currentData[3,*], spectralResponse)
         ;Ld = (((Lobs-Lu)/tau) - (Lt*ems))/(1-ems)
         Ld = (((obsRadianceZero-Lu)/tau)-(tempRadianceZero*ems))/(1-ems)         
            
         ;convert results to strings
         tau = STRING(tau)
         Lu = STRING(Lu)
         Ld = STRING(Ld)
         
         ;put results into results array
         results[*,place] = [lat, lon, height, tau, Lu, Ld]

         ;iterate place
         place = place+1
         
         print, place
              
      ;end for loop iterating heights                                       
      ENDFOR
   ;end for loop iterating narr points   
   ENDFOR
   
   ;determine file to output results to
   file = home+directory+'/atmosphericParameters.txt'
   
   ;write results to a file
   OPENW, unit, file, /GET_LUN
   PRINTF, unit, results
   CLOSE, unit
   FREE_LUN, unit


END
