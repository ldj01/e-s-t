;
; This is based on the original RIT prototype, with the following differences:
;
; - Instead of hardcoding emissivity, we read it from the emissivity band that
;   we have.
;

PRO generate_lst_image, home, directory, imagebase, which_landsat

; read in atmosphere layers
atmosphere_layers = READ_TIFF(home + directory + imagebase + '_LSTparams_MERRA.tif',GEOTIFF=lst_geotiff)

nonzeros = WHERE(atmosphere_layers[0,*,*] NE 0)

Lobs_array = atmosphere_layers[0,*,*]
tau_array  = atmosphere_layers[2,*,*]
Lu_array   = atmosphere_layers[3,*,*]
Ld_array   = atmosphere_layers[4,*,*]

Lobs = Lobs_array(nonzeros)
tau = tau_array(nonzeros)
Lu = Lu_array(nonzeros)
Ld = Ld_array(nonzeros)

array_size = SIZE(Lobs_array)
slice_size = SIZE(Lobs)

; psuedo emiss array
; emiss = DBLARR(slice_size[1])
; emiss[*] = 0.98996972

emiss_array = READ_TIFF(imagebase + '_emis.tif')
emiss = emiss_array(nonzeros)


;emiss_stdev = DBLARR(slice_size[1])
;emiss_stdev[*] = 0.0

;
; import LUT
;
CASE which_landsat OF
  5:  LUTfile = home + 'LUT5.txt'
  7:  LUTfile = home + 'LUT7_new.txt'
  10: LUTfile = home + 'LUT8_B10.txt'
  11: LUTfile = home + 'LUT8_B11.txt'
ENDCASE


nlines = FILE_LINES(LUTfile)
line = ''
count=-1
OPENR, lun, LUTfile, /GET_LUN
LUT = MAKE_ARRAY(2, nlines, /DOUBLE)
READF, lun, LUT
CLOSE, lun
FREE_LUN, lun


surface_radiance = (( (Lobs - Lu) / tau ) - (1 - emiss) * Ld) / emiss
HELP, surface_radiance
PRINT, SIZE(surface_radiance, /DIMENSIONS)
LST = CONVERT_RAD_TEMP(surface_radiance,LUT)

LST_image = DBLARR(array_size[1], array_size[2], array_size[3])
LST_image[nonzeros] = LST

WRITE_TIFF, home + directory + imagebase + '_lst.tif', LST_image, GEOTIFF=lst_geotiff, /FLOAT
END
