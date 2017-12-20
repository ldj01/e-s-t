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
;   10/18/2017 update:
;
; Several changes were made to the production ST procedure that require
; corresponding changes to this procedure for comparison purposes.
;
; - Fill locations for emissivity, emissivity standard deviation, and cloud 
;   distance are set to fill in the output product.  This is typically most
;   significant for emissivity.  Tau, Ld, and Lu need no be processed since
;   their fill values are based on thermal radiance.
; - The output ST QA product is scaled by 100.  The production version is int16.
;   WRITE_TIFF only offers unsigned int (/SHORT) so this is kept as FLOAT
;   (because the nodata value is -9999) but rounded.
; - The emissivity uncertainty algorithm from Glynn Hulley used used to 
;   compute S_E.
; - Landsat uncertainty values for L4, L5, and L8 are used.  Also the L7 value
;   is updated.  The new values are from Pat Scaramuzza.
; - Filenames are updated (for example, "st" is used instead of "lst").
;
;   11/30/2017 update:
;
;   Several changes were made to the prototye and production ST procedures
;   that require corresponding changes to this procedure for comparison 
;   purposes.
;
; - A problem was found with the instrument uncertainty values.  These were
;   in K, but should be in radiance to be consistent with other units.  The
;   values were converted for the supported satellites using the wavelength
;   range of the instruments used in the ST procedure.
; - A problem was found with the S_U unknown uncertainty value.  This was in
;   K, but was being used with other values in radiance units.  This needed
;   to be converted to radiance, and the final result needed to be converted
;   back to K.  The values are uncertainty in K/radiance, not K/radiance 
;   values, so for the conversions to work they need to be done in the nominal
;   range (e.g.: an uncertainty of 4K really is something like 280K +/- 2, and
;   the conversion needs to be done on 278K/282K, not 4K).
; - New debug steps were added to write bands for the new processing steps.
; - New test scenes were introduced. 
; - Also, in get_unknown_uncertainty.pro, the uncertainty matrix was updated.
;
;   Note that the planckInv function and inline radiance <--> temperature
;   conversions used here match the original prototype, but in the production
;   code, conversions based on k1/k2 metadata parameters for the instrument
;   are used instead.  Those conversions will be needed here if this will be
;   used to compare to the production QA script.
;


; input arguments should probably be pointers
PRO estimate_lst_uncertainty

MULT_FACTOR = 100.0

; Update this with whatever scene you are testing.
;imagebase = 'LE07_L1TP_014035_20150615_20160902_01_T1' ; L7 scene
; imagebase = 'LE07_L1TP_013033_20130330_20160908_01_T2' ; L7 scene
imagebase = 'LC08_L1TP_013032_20130330_20170310_01_T1' ; L8 scene
whichLandsat = strmid(imagebase, 3, 1)
satNum = fix(whichLandsat) 

; read in atmosphere layers
Lobs_array = READ_TIFF(imagebase + '_st_thermal_radiance.tif',GEOTIFF=lst_geotiff)
emis_array = READ_TIFF(imagebase + '_emis.tif')
emis_stdev_array = READ_TIFF(imagebase + '_emis_stdev.tif')
tau_array = READ_TIFF(imagebase + '_st_atmospheric_transmittance.tif')
Lu_array = READ_TIFF(imagebase + '_st_upwelled_radiance.tif')
Ld_array = READ_TIFF(imagebase + '_st_downwelled_radiance.tif')

nonfill = WHERE(Lobs_array[*,*] NE -9999)

Lobs = Lobs_array(nonfill)
tau = tau_array(nonfill)
Lu = Lu_array(nonfill)
Ld = Ld_array(nonfill)
emis = emis_array(nonfill)
emis_stdev = emis_stdev_array(nonfill)

array_size = SIZE(Lobs_array)
slice_size = SIZE(Lobs)

; read in distance image
dist2cloud_image = READ_TIFF(imagebase + '_st_cloud_distance.tif')
distance = dist2cloud_image(nonfill)

; read fill locations
emis_fill_locations = WHERE(emis_array[*,*] EQ -9999)
emis_stdev_fill_locations = WHERE(emis_stdev_array[*,*] EQ -9999)
distance_fill_locations = WHERE(distance[*,*] EQ -9999)

;
; Calculate partials
;
dLT_dTAU = (Lu - Lobs) / ( emis * tau^2)
dLT_dLU  = -1 / ( tau * emis )
dLT_dLD  = ( emis - 1) / emis
dLT_dLOBS = 1 / ( tau * emis )

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
; Instrument Uncertainty in radiance
;
CASE satNum OF
   4: landsat_uncertainty = 0.08959421
   5: landsat_uncertainty = 0.08959421
   7: landsat_uncertainty = 0.04471454
   8: landsat_uncertainty = 0.03577072
ENDCASE

CASE satNum OF
   4: K1 = 671.62
   5: K1 = 607.76
   7: K1 = 666.09
   8: K1 = 774.8853
ENDCASE
CASE satNum OF
   4: K2 = 1284.3
   5: K2 = 1260.56
   7: K2 = 1282.71
   8: K2 = 1321.0789
ENDCASE


S_I = ( dLT_dLOBS * landsat_uncertainty )^2 

;
; Emissivity Uncertainty
; This S_E algorithm is from Glynn Hulley.
;
; S_E = ( (Lu - Lobs + Ld * tau) / ( tau * emis^2 ) * emis_stdev )^2  
; 
CASE satNum OF
   4: emis_regfit = 0.00085135
   5: emis_regfit = 0.0013
   7: emis_regfit = 0.0011
   8: emis_regfit = 0.00093909
ENDCASE

eret13 = 0.0164
eret14 = 0.0174
eret = sqrt((eret13^2 + eret14^2) / 2)

S_E = sqrt((emis_stdev^2 + emis_regfit^2 + eret^2) / 3)

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

Le_Image=(Lobs_array-Lu_array-(1-emis_array)*Ld_array*tau_array)/(emis_array*tau_array)
Le_Image(where(Le_Image gt 30))=0
a=Le_Image
LST_Image=K2/(alog((K1/Le_Image)+1))

LST4Calculations = LST_Image(nonfill)
Le4Calculations = Le_Image(nonfill)

; THIS GIVES RESIDUAL TEMPERATURE ABOUT THE NOMINAL TEMPERATURE
Delta_Temps=S_U/2.0
LST_Minus_Delta=LST4Calculations-Delta_Temps
LST_Plus_Delta=LST4Calculations+Delta_Temps

; CONVERT TO RADIANCE
S_U_Radiance_High=K1/((exp(K2/LST_Plus_Delta)) - 1)
S_U_Radiance_Low=K1/((exp(K2/LST_Minus_Delta)) - 1)
S_U_Radiance=S_U_Radiance_High-S_U_Radiance_Low
;I'M GETTING ODD VALUES IN THE SURROUND SO I ZERO THEM OUT
S_U_Radiance(where(S_U_Radiance lt 0))=0

; CALCULATE THE UNCERTAINTY WITH COMMON UNITS (OF RADIANCE)
lst_uncertainty_radiance = sqrt(S_A + S_E + S_I + S_P + S_U_Radiance)

; THIS GIVES RESIDUAL RADIANCE
Delta_Radiance=lst_uncertainty_radiance/2.0
Radiance_Minus_Delta=Le4Calculations-Delta_Radiance
Radiance_Plus_Delta=Le4Calculations+Delta_Radiance

; CONVERT BACK TO TEMPERATURE
Temp_uncertainty_High=K2/(alog((K1/Radiance_Plus_Delta)+1))
Temp_uncertainty_Low=K2/(alog((K1/Radiance_Minus_Delta)+1))
Temp_Uncertainty=Temp_uncertainty_High-Temp_uncertainty_Low
;I'M GETTING ODD VALUES IN THE SURROUND SO I ZERO THEM OUT
Temp_Uncertainty(where(Temp_Uncertainty gt 100))=0

; PUT BACK IN IMAGE FORM
Temp_Uncertainty_Image = DBLARR(array_size[1], array_size[2])
Temp_Uncertainty_Image[nonfill] = Temp_Uncertainty


; Write some intermediate files for debugging purposes
S_ALL = S_A + S_I + S_E + S_P + S_U
S_A_out = DBLARR(array_size[1], array_size[2])
S_I_out = DBLARR(array_size[1], array_size[2])
S_E_out = DBLARR(array_size[1], array_size[2])
S_P_out = DBLARR(array_size[1], array_size[2])
S_U_out = DBLARR(array_size[1], array_size[2])
S_ALL_out = DBLARR(array_size[1], array_size[2])
Le4_out = DBLARR(array_size[1], array_size[2])
LST4_out = DBLARR(array_size[1], array_size[2])
Delta_Temps_out = DBLARR(array_size[1], array_size[2])
Delta_Radiance_out = DBLARR(array_size[1], array_size[2])
LST_Minus_Delta_out = DBLARR(array_size[1], array_size[2])
LST_Plus_Delta_out = DBLARR(array_size[1], array_size[2])
Radiance_Minus_Delta_out = DBLARR(array_size[1], array_size[2])
Radiance_Plus_Delta_out = DBLARR(array_size[1], array_size[2])
S_U_Radiance_High_out = DBLARR(array_size[1], array_size[2])
S_U_Radiance_Low_out = DBLARR(array_size[1], array_size[2])
S_U_Radiance_out = DBLARR(array_size[1], array_size[2])
Temp_Uncertainty_High_out = DBLARR(array_size[1], array_size[2])
Temp_Uncertainty_Low_out = DBLARR(array_size[1], array_size[2])
lst_uncertainty_radiance_out = DBLARR(array_size[1], array_size[2])
S_A_out[nonfill] = S_A
S_I_out[nonfill] = S_I
S_E_out[nonfill] = S_E
S_P_out[nonfill] = S_P
S_U_out[nonfill] = S_U
S_ALL_out[nonfill] = S_ALL
Le4_out[nonfill] = Le4Calculations
LST4_out[nonfill] = LST4Calculations
Delta_Temps_out[nonfill] = Delta_Temps
Delta_Radiance_out[nonfill] = Delta_Radiance
LST_Minus_Delta_out[nonfill] = LST_Minus_Delta
LST_Plus_Delta_out[nonfill] = LST_Plus_Delta 
Radiance_Minus_Delta_out[nonfill] = Radiance_Minus_Delta
Radiance_Plus_Delta_out[nonfill] = Radiance_Plus_Delta 
S_U_Radiance_High_out[nonfill] = S_U_Radiance_High 
S_U_Radiance_Low_out[nonfill] = S_U_Radiance_Low
S_U_Radiance_out[nonfill] = S_U_Radiance 
Temp_Uncertainty_High_out[nonfill] = Temp_Uncertainty_High 
Temp_Uncertainty_Low_out[nonfill] = Temp_Uncertainty_Low
lst_uncertainty_radiance_out[nonfill] = lst_uncertainty_radiance 
WRITE_TIFF, imagebase + '_SA.tif', S_A_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SI.tif', S_I_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SE.tif', S_E_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SP.tif', S_P_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SU.tif', S_U_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_SALL.tif', S_ALL_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Le4.tif', Le4_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_LST4.tif', LST4_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Delta_Temps.tif', Delta_Temps_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Delta_Radiance.tif', Delta_Radiance_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_LST_Minus_Delta.tif', LST_Minus_Delta_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_LST_Plus_Delta.tif', LST_Plus_Delta_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Radiance_Minus_Delta.tif', Radiance_Minus_Delta_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Radiance_Plus_Delta.tif', Radiance_Plus_Delta_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_S_U_Radiance_High.tif', S_U_Radiance_High_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_S_U_Radiance_Low.tif', S_U_Radiance_Low_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_S_U_Radiance.tif', S_U_Radiance_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Temp_Uncertainty_High.tif', Temp_Uncertainty_High_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_Temp_Uncertainty_Low.tif', Temp_Uncertainty_Low_out, GEOTIFF=lst_geotiff, /FLOAT
WRITE_TIFF, imagebase + '_lst_uncertainty_radiance.tif', lst_uncertainty_radiance_out, GEOTIFF=lst_geotiff, /FLOAT

lst_uncertainty_array = DBLARR(array_size[1], array_size[2])

lst_uncertainty_array = round(Temp_Uncertainty_Image * MULT_FACTOR)
lst_uncertainty_array[emis_fill_locations] = -9999 
lst_uncertainty_array[emis_stdev_fill_locations] = -9999 
lst_uncertainty_array[distance_fill_locations] = -9999 

WRITE_TIFF, imagebase + '_uncertainty.tif', lst_uncertainty_array, GEOTIFF=lst_geotiff, /FLOAT
END
