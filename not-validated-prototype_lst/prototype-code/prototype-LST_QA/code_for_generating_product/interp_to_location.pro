;Monica Cook
;19 October 2011

;
;Interpolate to location of current pixel
;


FUNCTION INTERP_TO_LOCATION, coordinates, data, interpolateToEasting, interpolateToNorthing

   coordinates = DOUBLE(coordinates)
   data = DOUBLE(data)
   enew = DOUBLE(interpolateToEasting)
   nnew = DOUBLE(interpolateToNorthing)
   
   ;shepard's method
   
   p = 1
   
   h = MAKE_ARRAY(N_ELEMENTS(coordinates[0,*]))
   w = MAKE_ARRAY(N_ELEMENTS(coordinates[0,*]))
      
   FOR i = 0, N_ELEMENTS(h)-1 DO BEGIN
   
      h[i] = SQRT( (coordinates[0,i]-enew)^2 + (coordinates[1,i]-nnew)^2 )
   
   ENDFOR
   
   denominator = TOTAL(h^(-p))
   
   FOR i = 0, N_ELEMENTS(w)-1 DO BEGIN
      
      w[i] = (h[i]^(-p))/denominator
      
   ENDFOR
   
   parameters = MAKE_ARRAY(3, /DOUBLE)
   FOR i = 0, 2 DO BEGIN
      
      parameters[i] = TOTAL(w*data[i,*])
      
   ENDFOR
            
   RETURN, parameters
   
END
