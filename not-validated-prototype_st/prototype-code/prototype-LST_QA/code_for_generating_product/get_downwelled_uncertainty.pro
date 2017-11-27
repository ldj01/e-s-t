; Dr. Kelly G. Laraby
; 04/15/2017
;
; This function calcualtes the downwelled radiance uncertainty term, which is part of the lst uncertainty estimation.
;




FUNCTION  get_downwelled_uncertainty, ld_values


  downwelled_uncertainty = MAKE_ARRAY(N_ELEMENTS(ld_values))

  ;
  ; Set a lower bound, above which a quadratic fit will be made,
  ; and below which last value of the fit line is extended as a constant.
  ; The lower bound was set as simply the smallest transmission value from
  ; the validation set that was used.
  ;
  ld_value_where_real_data_ends = 7.2307

  ;
  ; Transmission coefficients calculated from MODTRAN simulations using MERRA
  ;  
  ld_poly_coeff1 = -0.005291631300309
  ld_poly_coeff2 = 0.073763835328557
  ld_poly_coeff3 = -0.007066004902229

  ;
  ; Calculate uncertainty values for the quadratic region
  ;
  quadratic_region = WHERE(ld_values LE ld_value_where_real_data_ends)
  downwelled_uncertainty[quadratic_region] = ld_poly_coeff1 * ld_values[quadratic_region]^2 + ld_poly_coeff2 * ld_values[quadratic_region] + ld_poly_coeff3

  ;
  ; Calculate uncertainty values for the constant region
  ;
  constant_region = WHERE(ld_values GT ld_value_where_real_data_ends)
  downwelled_uncertainty[constant_region] = ld_poly_coeff1 * ld_value_where_real_data_ends^2 + ld_poly_coeff2 * ld_value_where_real_data_ends + ld_poly_coeff3

RETURN, downwelled_uncertainty
END
