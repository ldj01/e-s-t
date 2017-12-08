
; input arguments should probably be pointers
PRO estimate_lst_uncertainty_Test2, home, directory, imagebase, folderpath_to_distance_image

satelliteNumber = 8
;; read in atmosphere layers
;atmosphere_layers = READ_TIFF(home + directory + imagebase + '_LSTparams_MERRA.tif',GEOTIFF=lst_geotiff)

;              Lobs_array = read_binary('/Users/aarongerace/Desktop/LSTpart1/LC80380382016201LGN00/LE07_L1TP_016038_20111109_20160912_01_T1_st_thermal_radiance.img', DATA_TYPE=4, DATA_DIMS=[8031,7051,1])
;              tau_array  = read_binary('/Users/aarongerace/Desktop/LSTpart1/LC80380382016201LGN00/LE07_L1TP_016038_20111109_20160912_01_T1_st_atmospheric_transmittance.img', DATA_TYPE=4, DATA_DIMS=[8031,7051,1])
;              Lu_array   = read_binary('/Users/aarongerace/Desktop/LSTpart1/LC80380382016201LGN00/LE07_L1TP_016038_20111109_20160912_01_T1_st_upwelled_radiance.img', DATA_TYPE=4, DATA_DIMS=[8031,7051,1])
;              Ld_array   = read_binary('/Users/aarongerace/Desktop/LSTpart1/LC80380382016201LGN00/LE07_L1TP_016038_20111109_20160912_01_T1_st_downwelled_radiance.img', DATA_TYPE=4, DATA_DIMS=[8031,7051,1])
              Lobs_array = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayLobs.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
              tau_array  = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayTau.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
              Lu_array   = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayLu.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
              Ld_array   = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayLd.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
            
            
              nonzeros = WHERE(Lobs_array NE 0)
;nonzeros = WHERE(atmosphere_layers[0,*,*] NE 0)

;Lobs_array = atmosphere_layers[0,*,*]
;tau_array  = atmosphere_layers[2,*,*]
;Lu_array   = atmosphere_layers[3,*,*]
;Ld_array   = atmosphere_layers[4,*,*]

Lobs = Lobs_array(nonzeros)
tau = tau_array(nonzeros)
Lu = Lu_array(nonzeros)
Ld = Ld_array(nonzeros)

array_size = SIZE(Lobs_array)
slice_size = SIZE(Lobs)

;; read in distance image
;dist2cloud_image = READ_TIFF(folderpath_to_distance_image + imagebase + '_dist2cloud.tif')
              dist2cloud_image = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayCloudDistance.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
distance = dist2cloud_image(nonzeros)

    RayEmisUnc_array   = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayEmisUnc.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
    RayEmis_array   = read_binary('/Users/aarongerace/Desktop/Test_13_32/RayEmis.img', DATA_TYPE=4, DATA_DIMS=[7361,7321,1])
    
    emiss = RayEmis_array(nonzeros)
    emiss_stdev = RayEmisUnc_array(nonzeros)
;; psuedo emiss array
;emiss = DBLARR(slice_size[1])
;emiss[*] = 0.98996972
;
;emiss_stdev = DBLARR(slice_size[1])
;emiss_stdev[*] = 0.0



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
CASE satelliteNumber OF

  4: begin
    PRINT, 'four'
    landsat_uncertainty = 0.08959421
    end
  5: begin
    PRINT, 'five'
    landsat_uncertainty = 0.08959421
    K1 = 607.76
    K2 = 1260.56
    end
  7: begin
    PRINT, 'seven'
    landsat_uncertainty = 0.04471454
    K1 = 666.09
    K2 = 1282.71
;    K1_CONSTANT_BAND_6_VCID_2 = 666.09
;    K2_CONSTANT_BAND_6_VCID_2 = 1282.71
    end

  ELSE: begin
    PRINT, 'eight'
    landsat_uncertainty = 0.03577072
    K1 = 774.8853
    K2 = 1321.0789
;    K1_CONSTANT_BAND_11 = 480.8883
;    K2_CONSTANT_BAND_11 = 1201.1442
    end

ENDCASE

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

unknown_uncertainty = get_unknown_uncertainty_NewMatrix(distance, tau)
S_U = (unknown_uncertainty)^2


;  ###########  BEGIN CODE CHANGES  #########
      
;            ##############  THIS BLOCK JUST GENERATES THE LST IMAGE AND THE SURFACE EMITTED RADIANCE IMAGE
             
            Le_Image=(Lobs_array-Lu_array-(1-emiss)*Ld_array*tau_array)/(emiss*tau_array)
            Le_Image(where(abs(Le_Image) gt 30))=0
            LST_Image=K2/(alog((K1/Le_Image)+1))

;;            ##############  
      
            LST4Calculations = LST_Image(nonzeros)  ;  THIS PUT IN THE FORM THAT KELLY USES FOR THE UNCERTAINTY CALCULATIONS
            Le4Calculations = Le_Image(nonzeros)    ;  THIS PUT IN THE FORM THAT KELLY USES FOR THE UNCERTAINTY CALCULATIONS
      
      
;           THIS GIVES RESIDUAL TEMPERATURE ABOUT THE NOMINAL TEMPERATURE
            Delta_Temps=S_U/2.0
            LST_Minus_Delta=LST4Calculations-Delta_Temps
            LST_Plus_Delta=LST4Calculations+Delta_Temps

;;            CONVERT TO RADIANCE
;            S_U_Radiance_High=0.0005*LST_Plus_Delta^2-0.1798*LST_Plus_Delta+15.969
;            S_U_Radiance_Low=0.0005*LST_Minus_Delta^2-0.1798*LST_Minus_Delta+15.969
;            S_U_Radiance=S_U_Radiance_High-S_U_Radiance_Low
;            S_U_Radiance(where(S_U_Radiance lt 0))=0    ;I'M GETTING ODD VALUES IN THE SURROUND SO I ZERO THEM OUT
            S_U_Radiance_High=K1/((exp(K2/LST_Plus_Delta)) - 1)
            S_U_Radiance_Low=K1/((exp(K2/LST_Minus_Delta)) - 1)
            S_U_Radiance_Low(where(S_U_Radiance_Low lt 0))=0
            S_U_Radiance=S_U_Radiance_High-S_U_Radiance_Low
            S_U_Radiance(where(S_U_Radiance lt 0))=0    ;I'M GETTING ODD VALUES IN THE SURROUND SO I ZERO THEM OUT



;            CALCULATE THE UNCERTAINTY WITH COMMON UNITS (OF RADIANCE)
            lst_uncertainty_radiance = sqrt( S_A + S_I + S_P + S_U_Radiance)


;           THIS GIVES RESIDUAL RADIANCE 
            Delta_Radiance=lst_uncertainty_radiance/2.0
            Radiance_Minus_Delta=Le4Calculations-Delta_Radiance
            Radiance_Plus_Delta=Le4Calculations+Delta_Radiance
            
;;            CONVERT BACK TO TEMPERATURE
;            Temp_uncertainty_High=-0.0035*Radiance_Plus_Delta^4 + 0.1653*Radiance_Plus_Delta^3 - 2.9369*Radiance_Plus_Delta^2 + 29.627*Radiance_Plus_Delta + 170.1
;            Temp_uncertainty_Low=-0.0035*Radiance_Minus_Delta^4 + 0.1653*Radiance_Minus_Delta^3 - 2.9369*Radiance_Minus_Delta^2 + 29.627*Radiance_Minus_Delta + 170.1
;            Temp_Uncertainty=Temp_uncertainty_High-Temp_uncertainty_Low
;            Temp_Uncertainty(where(Temp_Uncertainty gt 100))=0  ;I'M GETTING ODD VALUES IN THE SURROUND SO I ZERO THEM OUT
            Temp_uncertainty_High=K2/(alog((K1/Radiance_Plus_Delta)+1))
            Temp_uncertainty_Low=K2/(alog((K1/Radiance_Minus_Delta)+1))
            Temp_Uncertainty=Temp_uncertainty_High-Temp_uncertainty_Low
            Temp_Uncertainty(where(Temp_Uncertainty gt 100))=0  ;I'M GETTING ODD VALUES IN THE SURROUND SO I ZERO THEM OUT

;            PUT BACK IN IMAGE FORM
            Temp_Uncertainty_Image = DBLARR(array_size[1], array_size[2]);, array_size[3])
            Temp_Uncertainty_Image[nonzeros] = Temp_Uncertainty
            envi_enter_data,Temp_Uncertainty_Image
            envi_enter_data,LST_Image
;  ###########  END CODE CHANGES  #########


;lst_uncertainty = sqrt( S_A + S_I + S_E + S_P + S_U)
;lst_uncertainty_array = DBLARR(array_size[1], array_size[2], array_size[3])
;lst_uncertainty_array[nonzeros] = lst_uncertainty

;WRITE_TIFF, home + directory + imagebase + '_uncertainty.tif', lst_uncertainty_array, GEOTIFF=lst_geotiff, /FLOAT
END

