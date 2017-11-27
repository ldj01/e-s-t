; Dr. Kelly G. Laraby
; 04/15/2017
;
; This function calcualtes the upwelled radiance uncertainty term, which is part of the lst uncertainty estimation.
;


FUNCTION  get_upwelled_uncertainty, lu_values

  upwelled_uncertainty = MAKE_ARRAY(N_ELEMENTS(lu_values))

  ;
  ; Set a lower bound, above which a quadratic fit will be made,
  ; and below which last value of the fit line is extended as a constant.
  ; The lower bound was set as simply the smallest transmission value from
  ; the validation set that was used.
  ;
  lu_value_where_real_data_ends = 5.6788
 
  ;
  ; Transmission coefficients calculated from MODTRAN simulations using MERRA
  ; 
  lu_poly_coeff1 = -0.007147307279714
  lu_poly_coeff2 = 0.082335806862813
  lu_poly_coeff3 = -0.006782188536986

  ;
  ; Calculate uncertainty values for the quadratic region
  ;
  quadratic_region = WHERE(lu_values LE lu_value_where_real_data_ends)
  upwelled_uncertainty[quadratic_region] = lu_poly_coeff1 * lu_values[quadratic_region]^2 + lu_poly_coeff2 * lu_values[quadratic_region] + lu_poly_coeff3

  ;
  ; Calculate uncertainty values for the constant region
  ;
  constant_region = WHERE(lu_values GT lu_value_where_real_data_ends)
  upwelled_uncertainty[constant_region] = lu_poly_coeff1 * lu_value_where_real_data_ends^2 + lu_poly_coeff2 * lu_value_where_real_data_ends + lu_poly_coeff3

RETURN, upwelled_uncertainty
END
