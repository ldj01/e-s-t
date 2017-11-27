

FUNCTION get_unknown_uncertainty, cloud_distances, transmission_values


  ;
  ; matrix of "unknown errors," which was calculated from observed and 
  ; predicted LST errors from the L7 global validation study.  
  ;
  unknown_error_matrix =  [ [2.3905, 2.7150, 2.5762, 2.1302], $
                            [2.0158, 1.7028, 1.4872, 1.3053], $
                            [1.8156, 1.0619, 0.9760, 0.7264], $
                            [1.9715, 1.3853, 0.8110, 0.7295], $
                            [1.4160, 0.8752, 0.7948, 0.4269] ]
                            

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
  
  
  ;
  ; expand the matrix (add column at position 0 that matches the original first
  ; column, add column at the end that matches the original last column, and 
  ; the same concept for the rows)
  ;
  extra_beginning = unknown_error_matrix[0,*]
  extra_ending = unknown_error_matrix[3,*]
  expanded1 = [extra_beginning, unknown_error_matrix, extra_ending]
  extra_beg = expanded1[*,0]
  extra_end = expanded1[*,4]
  expanded = [ [extra_beg], [expanded1], [extra_end] ]
  
  tau_interp_length = N_ELEMENTS(tau_interp)
  cld_interp_length = N_ELEMENTS(cld_interp)
  
  tau_full_vector = INTERPOLATE(tau_interp,0.001*INDGEN((tau_interp_length/.001),1))
  cld_full_vector = INTERPOLATE(cld_interp,0.001*INDGEN((cld_interp_length/.001),1))
  
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
  
 
 ;
 ; calculate interpolated unknown uncertainty
 ;
  unknown_uncertainty_interpolated = INTERPOLATE(unknown_error_matrix, tau_frac_index, cld_frac_index)

RETURN, unknown_uncertainty_interpolated
END