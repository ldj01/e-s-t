; Kelly G. Laraby
; PhD Candidate in Imaging Science
; Rochester Institute of Technology
; kga1099@rit.edu
;
; NAME:  
;    CONVERT_LL_UTM
;
; PURPOSE:  
;    IDL FUNCTION
;    This function converts latitude/longitude coordinates to UTM coordinates. The latitude and longitude information
;    can either be in scalar or vector form, and the ouput will follow suit with the zone, Easting and Northing values.
;    This function expects longitude to range from -180 to 180 degrees.
;    
;   
; OUTLINE:   
;    1) Convert degrees to radians
;    2) Define constants
;    3) Calculate meridonial arc length
;    4) Calculate UTM coefficients
;    5) Calculate Easting and Northing, and add 10 million to northing for southern hemisphere cases
;    6) Return zone, easting, and northing values
;
; CALL SEQUENCE: 
;    CONVERT_LL_UTM,        lat, $         ; latitude  (either single value or vector)
;                           lon, $         ; longitude (either single value or vector)
;                           utm            ; utm zone for lat/lon points (either single value or vector)
;                            
;                            
;                                        
; RESTRICTIONS/LIMITATIONS:
;    1) If vectors are given for latitude and longitude, but a scalar value is provided for the zone, it is assumed
;       that the zone should be used for all the lat/lon pairs.
;    2) Longitude values must range from -180 to 180 degrees rather than 0 to 360 degrees.
;    
;
; REQUIRED PROGRAMS (in working directory):        
;    NONE
; 
;                                                 
; REQUIRED FILES:      
;    NONE                                      
;
;
; MODIFICATIONS:
;    May,  2013        Original code
;    July, 2014        Added checks to make sure inputs are valid (e.g. they match in size)


FUNCTION CONVERT_LL_UTM, lat, lon, zone
    
   ;
   ; Make sure input arguements are valid (scalar or 1-D vectors, and same size)
   ;
   lat_size  = SIZE(lat)
   lon_size  = SIZE(lon)
   zone_size = SIZE(zone)
   
   IF lat_size(0) GT 1 OR lon_size(0) GT 1 OR zone_size(0) GT 1 THEN BEGIN
     MESSAGE, "Latitude, longitude, and zone values must be scalar or 1-D vectors."
   ENDIF ELSE IF lat_size(1) NE lon_size(1) THEN BEGIN
     MESSAGE, "Vectors for latitude and longitude must be the same length."
   ENDIF ELSE BEGIN
     numPoints = N_ELEMENTS(lat)
   ENDELSE
   
   ; If zone argument is scalar, the same value is used for all lat/lon pairs given.
   ; If the zone is a vector, it must have the same length as the lat/lon vectors.
   IF zone_size(0) EQ 1 THEN BEGIN
     IF zone_size(1) NE lon_size(1) THEN BEGIN
       MESSAGE, "UTM vector must be the same size as the lat/lon vectors."
     ENDIF ELSE BEGIN
       output_zone = zone
     ENDELSE
   ENDIF ELSE BEGIN
    output_zone = MAKE_ARRAY(1, numPoints, VALUE = zone)
   ENDELSE
    
   
   ; Calculate zone central meridian in degrees and radians
   zoneCM    = 6 * zone - 183
   zoneCMrad = (zoneCM * !DPI) / 180D
      
   ; Calculate delta longitude in degrees and radians
   deltaLon = lon - zoneCM
   p = (deltaLon * !DPI) / 180D
   
   ; Convert lat and lon points from decimal degrees to radians
   latrad = ( lat * !DPI ) / 180D
   lonrad = ( lon * !DPI ) / 180D
   
   ;
   ; Define datum constants
   ;
   a = 6378137D             ;equatorial radius
   b = 6356752.3            ;polar radius
   f = 0.0033528            ;flattening
   invf = 298.25722         ;inverse flattening
   mr = SQRT(a * b)         ;mean radius
   k0 = 0.9996              ;scale factor
   e = SQRT(1 - (b / a)^2)  ;eccentricity
   
   eprimesqrd = (e * e) / (1 - e * e)
   n = (a - b) / (a + b)
   rho = a * (1 - e * e) / ( (1 - (e * SIN(latrad))^2)^(3D/2) )
   nu = a / ( (1 - (e * SIN(latrad))^2)^(1D/2) )
  
   ;
   ; Calculate meridional arc length
   ; 
   A0 = a * (1 - n + (5 * n * n / 4) * (1 - n) + (81 * n^(4D/64)) * (1 - n))
   B0 = (3 * a * n / 2) * (1 - n - (7 * n * n / 8) * (1 - n) + 55 * n^(4D/64))
   C0 = (15 * a * n * n /16) * (1 - n + (3 * n * n / 4) * (1 - n))
   D0 = (35 * a * n^(3D/48)) * (1 - n + 11 * n * n / 16)
   E0 = (315 * a * n^(4D/51)) * (1 - n)
   
   S = A0*latrad - B0*SIN(2*latrad) + C0*SIN(4*latrad) - D0*SIN(6*latrad) + E0*SIN(8*latrad)
   
   ;
   ; Coefficients for UTM coordinates
   ;
   Ki = S * k0
   Kii = nu * SIN(latrad) * COS(latrad) * K0/2
   Kiii = ((nu * SIN(latrad) * COS(latrad)^3)/24) * (5 - TAN(latrad)^2 + 9*eprimesqrd * COS(latrad)^2 + 4*eprimesqrd^2 * COS(latrad)^4) * k0
   Kiv = nu * COS(latrad) * k0
   Kv = (cos(latrad))^3 * (nu/6) * (1 - TAN(latrad)^2 + eprimesqrd * COS(latrad)^2) * k0
   
   ;
   ; Calculate UTM coordinates
   ;
   northing = (Ki + Kii*p*p + Kiii * p^4)
   easting = 500000D + (Kiv * p + Kv * p^3)
   
   
   
   ; If the lat/lon location is in the southern hemisphere, you must add 10 million to the Northing result
   which_lats_are_south = WHERE(lat LT 0, /NULL)
   IF which_lats_are_south NE !NULL THEN northing(which_lats_are_south) = northing(which_lats_are_south) + 10000000D
   
   easting     = REFORM(easting,  1, numPoints)
   northing    = REFORM(northing, 1, numPoints)
   output_zone = REFORM(output_zone, 1, numPoints)
      
   ; Return pertinent values
   RETURN, [output_zone, easting, northing]

END