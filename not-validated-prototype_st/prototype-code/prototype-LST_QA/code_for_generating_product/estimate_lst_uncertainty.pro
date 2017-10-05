
; input arguments should probably be pointers
PRO estimate_lst_uncertainty, home, directory, imagebase, folderpath_to_distance_image

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

; read in distance image
dist2cloud_image = READ_TIFF(folderpath_to_distance_image + imagebase + '_dist2cloud.tif')
distance = dist2cloud_image(nonzeros)

; psuedo emiss array
emiss = DBLARR(slice_size[1])
emiss[*] = 0.98996972

emiss_stdev = DBLARR(slice_size[1])
emiss_stdev[*] = 0.0



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

; Atmosphere Uncertainty
S_A = ( dLT_dTAU * S_TAU )^2 + ( dLT_dLU * S_LU )^2 + ( dLT_dLD * S_LD )^2 
SA_tau_part = dLT_dTAU * S_TAU
SA_lu_part  = dLT_dLU * S_LU
SA_ld_part  = dLT_dLD * S_LD



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

unknown_uncertainty = get_unknown_uncertainty(distance, tau)
S_U = (unknown_uncertainty)^2

lst_uncertainty = sqrt( S_A + S_I + S_E + S_P + S_U)
lst_uncertainty_array = DBLARR(array_size[1], array_size[2], array_size[3])
lst_uncertainty_array[nonzeros] = lst_uncertainty

WRITE_TIFF, home + directory + imagebase + '_uncertainty.tif', lst_uncertainty_array, GEOTIFF=lst_geotiff, /FLOAT
END
