;Monica Cook
;26 September 2011

;
; Convert array of geopotential heights to array of geometric heights given latitude
;

FUNCTION CONVERT_GEOPOTENTIAL_GEOMETRIC, geoHeight, lat

   ;convert input variables
   geopotential = geoHeight
   radLat = lat*!DPI/180

   ;determine number of pressure levels
   size = SIZE(geopotential, /DIMENSIONS)
   numPressures = size[1]

   ;define constants
   g_0 = 9.80665 ;ms^-2
   R_max = 6378.137 ;km
   R_min = 6356.752 ;km
   
   ;define variable based on latitude
   g = 9.80616*(1-0.002637*COS(2*radlat)+0.0000059*(COS(2*radlat))^2)
   radius = SQRT(1/(((COS(radlat))^2/R_max^2)+((SIN(radlat))^2/R_min^2))) ;km
   radius = radius * 1000 ;m
   gravityRatio = g/g_0
   
   geometric = MAKE_ARRAY(N_ELEMENTS(lat), numPressures, /DOUBLE)
   
   ;convert geopotential height to geometric height
   
   FOR i = 0, N_ELEMENTS(lat)-1 DO BEGIN
      geometric[i,*] = (geopotential[i,*]*radius[i])/(gravityRatio[i]*radius[i] - geopotential[i,*])
   ENDFOR

   RETURN, geometric


END
