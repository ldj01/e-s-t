; Monica Cook
; 9 November 2011

; NAME:  CONVERT_RAD_TEMP
;  
; PURPOSE:  IDL FUNCTION
;   Convert from radiance to temperature
;
; CALL SEQUENCE: CONVERT_RAD_TEMP, radiance, $    ;radiance due to temperature
;                                  LUT            ;look up table relating radiance and temperature
;    
;  RESTRICTIONS:
;    Already generated results for current Landsat Scene
;
;  REQUIRED PROGRAMS AND FILES (in home directory):   RAD_TO_TEMP.pro
;       

FUNCTION CONVERT_RAD_TEMP, L, LUT

   ;input is given radiance, want to find corresponding temperature
   targetRadiance = L
   ;extract temperatures and radiances from look up table
   temp = LUT[0,*]
   radiance = LUT[1,*]
   
   ;determine size of target radiance and temperature arrays
   Radsize = SIZE(targetRadiance,/DIMENSIONS)
   IF N_ELEMENTS(Radsize) EQ 1 THEN BEGIN
      IF Radsize[0] EQ 1 THEN BEGIN
        numSamples=1
        numLines=1
        targetRadiance=MAKE_ARRAY(1,1,VALUE=targetRadiance)
      ENDIF ELSE BEGIN
        numSamples=Radsize[0]
        numLines=1
      ENDELSE
   ENDIF ELSE BEGIN
      numSamples = Radsize[0]
      numLines = Radsize[1]
   ENDELSE   
   
   targetTemp = MAKE_ARRAY(numSamples, numLines, /DOUBLE)
   
   FOR i = 0, numSamples-1 DO BEGIN
      FOR j = 0, numLines - 1 DO BEGIN
      
      IF targetRadiance[i,j] NE 0 THEN BEGIN
      
         ;find index in table of closest value below input radiance
         below = WHERE(radiance LE targetRadiance[i,j])
         closestBelow = below[N_ELEMENTS(below)-1]
   
         ;find index in table of closest value above input radiance
         above = WHERE(radiance GE targetRadiance[i,j])
         closestAbove = above[0]
   
         ;extract radiances from table just below and just above input radiance for linear interpolation
         underRadiance = radiance[closestBelow]
         overRadiance = radiance[closestAbove]
   
         ;extract temperatures from table just below and just above input radiance for linear interpolation 
         underTemp = temp[closestBelow]
         overTemp = temp[closestAbove]
   
         ;calculate slope and intercept from two extracted points (temp, radiance) above 
         slope = (overRadiance - underRadiance)/(overTemp - underTemp)
         intercept = overRadiance - slope*overTemp
   
         ;calcualte target temperature from linear interpolation of two nearest points
         targetTemp[i,j] = (targetRadiance[i,j] - intercept)/(slope)
                           
      ENDIF
         
      ENDFOR
   ENDFOR
 
   ;return target temperature
   RETURN, targetTemp 
   
END
