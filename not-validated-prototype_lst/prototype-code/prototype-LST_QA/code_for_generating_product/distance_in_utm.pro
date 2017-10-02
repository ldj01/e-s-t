;Monica Cook
;17 October 2011

;
;Calculate distances between UTM coordiantes
;

; MODIFICATIONS
;
; 9/28/15 -- subtract eastings by 500,000 m
;            correct simpson's rule

FUNCTION DISTANCE_IN_UTM, e1, n1, e2, n2
   
   east1 = e1 - 500000D
   east2 = e2 - 500000D

   s = 0.9996D ;scale factor
   r = 6378137D ;Earth radius
   
   SR1 = s/(COS(east1/r))
   SR2 = s/(COS(((east2+east1)/2)/r))
;  SR2 = s/(COS(((e2-e1)/6)/r))
   SR3 = s/(COS(east2/r))
   
   Edist = ((e2-e1)/6)*(SR1+4*SR2+SR3)
   
   d = SQRT(Edist^2+(n2-n1)^2)
    
   RETURN, d
   
END
