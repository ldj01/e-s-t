; Dr. Kelly G. Laraby
; 04/15/2017
;
; This function calcualtes the transmission uncertainty term, which is part of the lst uncertainty estimation.
;

FUNCTION  get_transmission_uncertainty, tau_values

  tau_uncertainty = MAKE_ARRAY(N_ELEMENTS(tau_values))

  ;
  ; Set a lower bound, above which a quadratic fit will be made,
  ; and below which last value of the fit line is extended as a constant. 
  ; The lower bound was set as simply the smallest transmission value from 
  ; the validation set that was used.
  ;
  tau_value_where_real_data_ends = 0.3005
  
  ;
  ; Transmission coefficients calculated from MODTRAN simulations using MERRA
  ;
  tau_poly_coeff1 = -0.050527295343549
  tau_poly_coeff2 = 0.029930689405143
  tau_poly_coeff3 = 0.019127148003052

  ;
  ; Calculate uncertainty values for the quadratic region
  ;  
  quadratic_region = WHERE(tau_values GE tau_value_where_real_data_ends)
  tau_uncertainty[quadratic_region] = tau_poly_coeff1 * tau_values[quadratic_region]^2 + tau_poly_coeff2 * tau_values[quadratic_region] + tau_poly_coeff3

  ;
  ; Calculate uncertainty values for the constant region
  ;
  constant_region = WHERE(tau_values LT tau_value_where_real_data_ends)
  tau_uncertainty[constant_region] = tau_poly_coeff1 * tau_value_where_real_data_ends^2 + tau_poly_coeff2 * tau_value_where_real_data_ends + tau_poly_coeff3

RETURN, tau_uncertainty
END
