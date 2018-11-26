PRO build_merra_coordinate_file

   ;Build the MERRA coordinate text file.  This is based on some code 
   ;extracted from the RIT prototype.  It uses the same format at the 
   ;NARR coordinate file.

   ;define sampling for latitude and longitude of merra grid
   deltalambda = 1.25
   deltaphi = 1.25

   ;define lon grid for MERRA data
   i = INDGEN(288)+1
   lon = -180 + deltalambda*(i-0.5)
   longrid = MAKE_ARRAY(288,144, /DOUBLE)
   eye = MAKE_ARRAY(288,144, /DOUBLE)
   FOR g = 0, 143 DO BEGIN
      longrid[*,g] = lon
      eye[*,g] = i
   ENDFOR

   ;define lat grid for MERRA data
   j = INDGEN(144)+1
   lat = -90 + deltaphi*(j-0.5)
   latgrid = MAKE_ARRAY(288,144, /DOUBLE)
   jay = MAKE_ARRAY(288,144, /DOUBLE)
   FOR g = 0, 287 DO BEGIN
      latgrid[g,*] = lat
      jay[g,*] = j
   ENDFOR

   FOR x = 1, 144 DO BEGIN
       FOR y = 1, 288 DO BEGIN
           PRINT, y, ' ', x, ' ', lat[x-1], ' ', lon[y-1]
       ENDFOR
   ENDFOR

END
