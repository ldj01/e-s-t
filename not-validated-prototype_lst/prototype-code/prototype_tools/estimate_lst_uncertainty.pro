; This tool is similar to the prototype estimate_lst_uncertainty script.  The
; differences include:
;
; - It operates on Collection 1 input. 
; - It expects all input in the current directory.
; - It requires the input scene to be hardcoded instead of using a parameter.
; - Instead of building an emissivity band hardcoded to water emissivity, and
;   an emissivity standard deviation band hardcoded to 0, it reads these bands.
; - It takes as input the production atmospheric parameteter bands.  These are
;   separate band files, instead of a single file as in the prototype.
; - Before running it, the production ST procedure must be run so the pixel_qa,
;   atmospheric parameter, and emissivity bands are available.
; - Before running it, convert_espa_to_gtif must be run, since it (like the
;   prototype version) expects input to be in geotiff format.
; - It produces a number of intermediate bands used for comparison/debugging.
; - It expects the get_unknown_uncertainty script to be available with similar
;   updates.
; - It expects the other "generate QA" scripts to be available.  These are
;   unchanged from the original prototype.
; - It excludes processing on input thermal radiance fill locations instead of 
;   zero locations.
; - It fixes what looks like a bug in the prototype when writing the final
;   output file.
; - It eliminates some steps in the prototype that are not used.
;
;   The purpose of this script is to assist in development/debugging.  It is 
;   similar to the original prototype, but with the above changes.  It allows 
;   comparison with the production code results at the various steps where 
;   debug outputs are created.
;

; input arguments should probably be pointers
PRO estimate_lst_uncertainty

; Update this with whatever scene you are testing.
imagebase = 'LE07_L1TP_014035_20150615_20160902_01_T1' ; L7 scene

; read in atmosphere layers
Lobs_array = READ_TIFF(imagebase + '_lst_thermal_radiance.tif',GEOTIFF=lst_geotiff)
emiss_array = READ_TIFF(imagebase + '_landsat_emis.tif')
emiss_stdev_array = READ_TIFF(imagebase + '_landsat_emis_stdev.tif')
tau_array = READ_TIFF(imagebase + '_lst_atmospheric_transmittance.tif')
Lu_array = READ_TIFF(imagebase + '_lst_upwelled_radiance.tif')
Ld_array = READ_TIFF(imagebase + '_lst_downwelled_radiance.tif')

nonfill = WHERE(Lobs_array[*,*] NE -9999)

Lobs = Lobs_array(nonfill)
tau = tau_array(nonfill)
Lu = Lu_array(nonfill)
Ld = Ld_array(nonfill)
emiss = emiss_array(nonfill)
emiss_stdev = emiss_stdev_array(nonfill)

array_size = SIZE(Lobs_array)
slice_size = SIZE(Lobs)

; read in distance image
dist2cloud_image = READ_TIFF(imagebase + '_cloud_distance.tif')
distance = dist2cloud_image(nonfill)

;
; Calculate partials
;
dLT_dTAU = (Lu - Lobs) / ( emiss * tau^2)
dLT_dLU  = -1 / ( tau * emiss )
dLT_dLD  = ( emiss - 1) / emiss
dLT_dLOBS = 1 / ( tau * emiss )

S_TAU = get_transmission_uncertainty(tau)
S_LU  = get_upwelled_uncertainty(Lu)
S_LD  = get_downwelled_uncertainty(Ld)

; Write some intermediate files for debugging purposes
dLT_dTAU_out = DBLARR(array_size[1], array_size[2])
dLT_dLU_out = DBLARR(array_size[1], array_size[2])
dLT_dLD_out = DBLARR(array_size[1], array_size[2])
dLT_dLOBS_out = DBLARR(array_size[1], array_size[2])
S_TAU_out = DBLARR(array_size[1], array_size[2])
S_LU_out = DBLARR(array_size[1], array_size[2])
S_LD_out = DBLARR(array_size[1], array_size[2])
dist_out = DBLARR(array_size[1], array_size[2])
dLT_dTAU_out[nonfill] = dLT_dTAU
dLT_dLU_out[nonfill] = dLT_dLU
dLT_dLD_out[nonfill] = dLT_dLD
dLT_dLOBS_out[nonfill] = dLT_dLOBS
S_TAU_out[nonfill] = S_TAU
S_LU_out[nonfill] = S_LU
S_LD_out[nonfill] = S_LD
dist_out[nonfill] = distance
WRITE_TIFF, imagebase + '_dLT_dTAU.tif', dLT_dTAU_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_dLT_dLU.tif', dLT_dLU_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_dLT_dLD.tif', dLT_dLD_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_dLT_dLOBS.tif', dLT_dLOBS_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_S_TAU.tif', S_TAU_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_S_LU.tif', S_LU_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_S_LD.tif', S_LD_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_dist.tif', dist_out, GEOTIFF=lst_geotiff, /FLOAT

; Atmosphere Uncertainty
S_A = ( dLT_dTAU * S_TAU )^2 + ( dLT_dLU * S_LU )^2 + ( dLT_dLD * S_LD )^2 

;
; Instrument Uncertainty
;
landsat_uncertainty = 0.032 ; [K], for Landsat 7 
S_I = ( dLT_dLOBS * landsat_uncertainty )^2 

;
; Emissivity Uncertainty
;
S_E = ( (Lu - Lobs + Ld * tau) / ( tau * emiss^2 ) * emiss_stdev )^2 

;
; Cross correlation terms
;
S_P = get_cross_correlation(dLT_dTAU, dLT_dLU, dLT_dLD, S_TAU, S_LU, S_LD)


; Unknown Uncertainty

;unknown_uncertainty = get_unknown_uncertainty(distance, tau)
; The last 3 parameters are not needed in the original prototype, but they
; are used here to help build debug files. 
unknown_uncertainty = get_unknown_uncertainty(distance, tau, lst_geotiff, array_size, nonfill)
S_U = (unknown_uncertainty)^2

; Write some intermediate files for debugging purposes
S_ALL = S_A + S_I + S_E + S_P + S_U
S_A_out = DBLARR(array_size[1], array_size[2])
S_I_out = DBLARR(array_size[1], array_size[2])
S_E_out = DBLARR(array_size[1], array_size[2])
S_P_out = DBLARR(array_size[1], array_size[2])
S_U_out = DBLARR(array_size[1], array_size[2])
S_ALL_out = DBLARR(array_size[1], array_size[2])
S_A_out[nonfill] = S_A
S_I_out[nonfill] = S_I
S_E_out[nonfill] = S_E
S_P_out[nonfill] = S_P
S_U_out[nonfill] = S_U
S_ALL_out[nonfill] = S_ALL
WRITE_TIFF, imagebase + '_SA.tif', S_A_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SI.tif', S_I_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SE.tif', S_E_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SP.tif', S_P_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SU.tif', S_U_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SALL.tif', S_ALL_out, GEOTIFF=lst_geotiff, /FLOAT

lst_uncertainty = sqrt( S_A + S_I + S_E + S_P + S_U)

lst_uncertainty_array = DBLARR(array_size[1], array_size[2])

lst_uncertainty_array[nonfill] = lst_uncertainty

WRITE_TIFF, imagebase + '_uncertainty.tif', lst_uncertainty_array, GEOTIFF=lst_geotiff, /FLOAT
END
