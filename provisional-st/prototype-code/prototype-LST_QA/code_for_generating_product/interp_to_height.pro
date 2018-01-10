;Monica Cook
;19 October 2011

;
;Inteprolate to height of current pixel
;

FUNCTION INTERP_TO_HEIGHT, height, data, interpolateTo

   ;Determine points below interpolateTo
   below = WHERE(height LT interpolateTo)
   IF below[0] EQ -1 THEN closestBelow = 0 ELSE closestBelow = below[N_ELEMENTS(below)-1]
   underHeight = height[closestBelow]
   underVariable = [data[0,closestBelow], data[1,closestBelow], data[2,closestBelow]]
   
   ;Determine points above interpolateTo
   above = WHERE(height GE interpolateTo)
   IF above[0] Eq -1 THEN closestAbove = N_ELEMENTS(height)-1 ELSE closestAbove = above[0]
   overHeight = height[closestAbove]
   overVariable = [data[0,closestAbove], data[1,closestAbove], data[2,closestAbove]]
   
   ;check that height falls within range
   ;if it does, interpolate parameters
   ;if it doesn't, use first or last value of each variable
      
   ;Calculate slope and intercept
   m = (overVariable-underVariable)/(overHeight-underHeight)
   b = overVariable - m*overHeight
   
   ;Calculate and return value of new variable
   new = m*interpolateTo + b
   IF closestAbove EQ closestBelow THEN new = underVariable
         
   RETURN, new
   
END
