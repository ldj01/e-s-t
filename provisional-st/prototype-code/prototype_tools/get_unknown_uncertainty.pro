; This is similar to the original prototype version, with the following
; changes:
;
; - It is expected to be run with the modified estimate_lst_uncertainty script,
;   with all conditions that implies.
; - Like the modified estimate_lst_uncertainty script, this builds some 
;   intermediate debug files.
; - Like the modified estimate_lst_uncertainty script, this eliminates some 
;   steps that aren't used. 

;FUNCTION get_unknown_uncertainty, cloud_distances, transmission_values
FUNCTION get_unknown_uncertainty, cloud_distances, transmission_values, lst_geotiff, array_size, nonfill


  ;
  ; matrix of "unknown errors," which was calculated from observed and 
  ; predicted LST errors from the L7 global validation study.  
  ;
  unknown_error_matrix =  [ [2.3749, 2.6962, 2.5620, 2.1131], $
                            [1.9912, 1.6789, 1.4471, 1.2739], $
                            [1.7925, 1.0067, 0.9143, 0.6366], $
                            [1.9416, 1.3558, 0.7604, 0.6682], $
                            [1.3861, 0.8269, 0.7404, 0.3125] ]

  ;
  ; tau bins are 0.3 - 0.55, 0.55 - 0.7, 0.7 - 0.85, 0.85 - 1.0
  ; cloud bins are 0 - 1 km, 1 - 5 km, 5 - 10 km, 10 - 40 km, 40 - inf
  ; 
  ; so tau_interp and cloud_interp should be a vector of the center 
  ; values for each bin, but we also want anything outside the entire range
  ; to be equal to the nearest value from the unknown matrix.
  ;
  tau_interp = [0.0,  0.425, 0.625, 0.775, 0.925, 1.0]
  cld_interp = [0, .5, 3, 7.5, 25, 82.5, 200]
  
  tau_interp_length = N_ELEMENTS(tau_interp)
  cld_interp_length = N_ELEMENTS(cld_interp)
  
  ;
  ; from input transmission values, find closest indices from tau_interp vector, calculate
  ; step in vector, then calculate fractional tau index.
  ;
  tau_close_index = VALUE_LOCATE(tau_interp, transmission_values)
  tau_step = tau_interp[tau_close_index+1] - tau_interp[tau_close_index]
  tau_frac_index = tau_close_index + ((transmission_values - tau_interp[tau_close_index]) / tau_step)
  tau_frac_index[WHERE(transmission_values EQ 1.0,/NULL)] = tau_interp_length - 1

  ;
  ; from input cloud distance values, find closest indices from cld_interp vector, calculate
  ; step in vector, then calculate fractional cloud index.
  ;
  cld_close_index = VALUE_LOCATE(cld_interp, cloud_distances)
  cld_step = cld_interp[cld_close_index+1] - cld_interp[cld_close_index]
  cld_frac_index = cld_close_index + ((cloud_distances - cld_interp[cld_close_index]) / cld_step)
  cld_frac_index[WHERE(cloud_distances EQ 200.0,/NULL)] = cld_interp_length - 1

  ; Write some intermediate files for debugging purposes. 
  tau_close_index_out = DBLARR(array_size[1], array_size[2])
  tau_step_out = DBLARR(array_size[1], array_size[2])
  tau_frac_index_out = DBLARR(array_size[1], array_size[2])
  tau_close_index_out[nonfill] = tau_close_index 
  tau_step_out[nonfill] = tau_step
  tau_frac_index_out[nonfill] = tau_frac_index 
  WRITE_TIFF, 'tau_close.tif', tau_close_index_out, GEOTIFF=lst_geotiff, /FLOAT
  WRITE_TIFF, 'tau_step.tif', tau_step_out, GEOTIFF=lst_geotiff, /FLOAT
  WRITE_TIFF, 'tau_frac.tif', tau_frac_index_out, GEOTIFF=lst_geotiff, /FLOAT
  cld_close_index_out = DBLARR(array_size[1], array_size[2])
  cld_step_out = DBLARR(array_size[1], array_size[2])
  cld_frac_index_out = DBLARR(array_size[1], array_size[2])
  cld_close_index_out[nonfill] = cld_close_index 
  cld_step_out[nonfill] = cld_step
  cld_frac_index_out[nonfill] = cld_frac_index 
  WRITE_TIFF, 'cld_close.tif', cld_close_index_out, GEOTIFF=lst_geotiff, /FLOAT
  WRITE_TIFF, 'cld_step.tif', cld_step_out, GEOTIFF=lst_geotiff, /FLOAT
  WRITE_TIFF, 'cld_frac.tif', cld_frac_index_out, GEOTIFF=lst_geotiff, /FLOAT
 
 ;
 ; calculate interpolated unknown uncertainty
 ;
  unknown_uncertainty_interpolated = INTERPOLATE(unknown_error_matrix, tau_frac_index, cld_frac_index)

  ; Write some intermediate files for debugging purposes. 
  unknown_uncertainty_interpolated_out = DBLARR(array_size[1], array_size[2])
  unknown_uncertainty_interpolated_out[nonfill] = unknown_uncertainty_interpolated
  WRITE_TIFF, 'unknown_uncertainty.tif', unknown_uncertainty_interpolated_out, GEOTIFF=lst_geotiff, /FLOAT

RETURN, unknown_uncertainty_interpolated
END
