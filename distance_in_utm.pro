;Monica Cook
;17 October 2011

;
;Calculate distances between UTM coordiantes
;

FUNCTION DISTANCE_IN_UTM, e1, n1, e2, n2

   s = 0.9996 ;scale factor
   r = 6378137 ;Earth radius
   
   SR1 = s/(COS(e1/r))
   SR2 = s/(COS(((e2-e1)/6)/r))
   SR3 = s/(COS(e2/r))
   
   Edist = ((e2-e1)/6)*(SR1+4*SR2+SR3)
   
   d = SQRT(Edist^2+(n2-n1)^2)
   
   RETURN, d
   
END
