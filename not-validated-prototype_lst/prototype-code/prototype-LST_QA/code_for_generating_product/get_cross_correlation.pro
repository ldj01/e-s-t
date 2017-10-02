; Dr. Kelly G. Laraby
; 04/15/2017
;
; This function calcualtes the cross correlation term, which is part of the lst uncertainty estimation.
;


FUNCTION  get_cross_correlation, dLT_dTAU, dLT_dLU, dLT_dLD, S_TAU, S_LU, S_LD

  ;
  ; Correlation coefficients from MODTRAN simulations using MERRA.
  ;
  corr_tau_lu = -0.9899;
  corr_tau_ld = -0.9857;
  corr_lu_ld  =  0.9965;

  ;
  ; Calculate cross correlation terms
  ;
  part_tau_lu = 2 * corr_tau_lu * dLT_dTAU * dLT_dLU * S_TAU * S_LU;
  part_tau_ld = 2 * corr_tau_ld * dLT_dTAU * dLT_dLD * S_TAU * S_LD;
  part_lu_ld =  2 * corr_lu_ld  * dLT_dLU  * dLT_dLD * S_LU  * S_LD;

  ;
  ; Calculate cross correlation
  ;
  cross_correlation = part_tau_lu + part_tau_ld + part_lu_ld;

RETURN, cross_correlation
END
