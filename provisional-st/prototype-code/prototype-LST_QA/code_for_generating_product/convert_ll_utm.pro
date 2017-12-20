;Monica Cook
;17 October 2011

;
;Convert Lat,Lon coordinates to UTM coordinates in specified zone
;


FUNCTION CONVERT_LL_UTM, lat, lon, zone

   numPoints = N_ELEMENTS(lon)
   
   ;initialize array containing zone
   zone = MAKE_ARRAY(1, numPoints, VALUE = zone)
   
   ;calculate zone central meridian in degrees and radians
   zoneCM = 6*zone-183
   zoneCMrad = zoneCM*!DPI/180D
      
   ;calculate delta longitude in degrees and radians
   deltaLon = lon - zoneCM
   p = deltaLon*!DPI/180D
   
   ;convert lat and lon points from decimal degrees to radians
   latrad = lat*!DPI/180D
   lonrad = lon*!DPI/180D
   
   ;define datum constants
   a = 6378137D  ;equatorial radius
   b = 6356752.3  ;polar radius
   f = 0.0033528  ;flattening
   invf = 298.25722  ;inverse flattening
   mr = SQRT(a*b)  ;mean radius
   
   k0 = 0.9996  ;scale factor
   e = SQRT(1-(b/a)^2)  ;eccentricity
   eprimesqrd = e*e/(1-e*e)
   n = (a-b)/(a+b)
   
   rho = a*(1-e*e)/((1-(e*SIN(latrad))^2)^(3D/2))
   nu = a/((1-(e*SIN(latrad))^2)^(1D/2))
   
   ;calculate meridional arc length
   
   A0 = a*(1-n+(5*n*n/4)*(1-n)+(81*n^4D/64)*(1-n))
   B0 = (3*a*n/2)*(1-n-(7*n*n/8)*(1-n)+55*n^4D/64)
   C0 = (15*a*n*n/16)*(1-n+(3*n*n/4)*(1-n))
   D0 = (35*a*n^3D/48)*(1-n+11*n*n/16)
   E0 = (315*a*n^4D/51)*(1-n)
   
   S = A0*latrad - B0*SIN(2*latrad) + C0*SIN(4*latrad) - D0*SIN(6*latrad) + E0*SIN(8*latrad)
   
   ;coefficients for UTM coordinates
   
   Ki = S*k0
   Kii = nu*SIN(latrad)*COS(latrad)*K0/2
   Kiii = ((nu*SIN(latrad)*COS(latrad)^3)/24)*(5-TAN(latrad)^2+9*eprimesqrd*COS(latrad)^2+4*eprimesqrd^2*COS(latrad)^4)*k0
   Kiv = nu*COS(latrad)*k0
   Kv = (cos(latrad))^3*(nu/6)*(1-TAN(latrad)^2+eprimesqrd*COS(latrad)^2)*k0
   
   ;calculate UTM coordinates
   
   northing = (Ki+Kii*p*p+Kiii*p^4)
   easting = 500000 + (Kiv*p+Kv*p^3)
   
   easting = REFORM(easting, 1, numPoints)
   northing = REFORM(northing, 1, numPoints)
      
   ;return pertinent values
   RETURN, [zone, easting, northing]

END