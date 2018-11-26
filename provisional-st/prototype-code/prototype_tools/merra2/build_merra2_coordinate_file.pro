PRO build_merra2_coordinate_file

   ;Build the MERRA-2 coordinate text file.  This is based on some code
   ;extracted from the RIT prototype, but modified for MERRA-2.  It uses 
   ;the same format at the NARR coordinate file.

   ;define sampling for latitude and longitude of merra-2 grid
   deltalambda = 0.625
   deltaphi = 0.5

   ;define lon grid for MERRA-2 data
   i = INDGEN(576)+1
   lon = -180 + deltalambda*(i-1)
   longrid = MAKE_ARRAY(576,361, /DOUBLE)
   eye = MAKE_ARRAY(576,361, /DOUBLE)
   FOR g = 0, 360 DO BEGIN
      longrid[*,g] = lon
      eye[*,g] = i
   ENDFOR

   ;define lat grid for MERRA-2 data
   j = INDGEN(361)+1
   lat = -90 + deltaphi*(j-1)
   latgrid = MAKE_ARRAY(576,361, /DOUBLE)
   jay = MAKE_ARRAY(576,361, /DOUBLE)
   FOR g = 0, 575 DO BEGIN
      latgrid[g,*] = lat
      jay[g,*] = j
   ENDFOR

   FOR x = 1, 361 DO BEGIN
       FOR y = 1, 576 DO BEGIN
           PRINT, y, ' ', x, ' ', lat[x-1], ' ', lon[y-1]
       ENDFOR
   ENDFOR

END
