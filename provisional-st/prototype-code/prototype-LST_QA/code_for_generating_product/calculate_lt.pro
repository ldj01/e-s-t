;Monica Cook
;5 October 2011

;
;Calculate blackbody radiance from temperature using spectral repsonse function
;

FUNCTION CALCULATE_LT, temperature, spectralResponse

   ;split spectral reponse into wavelength and repsonse
   landsatWavelength = spectralResponse[0,*]
   landsatResponse = spectralResponse[1,*]
   
   ;integrate spectral response over wavelength
   RSRintegral = INT_TABULATED(landsatWavelength, landsatResponse)
   
   ;using planck's equaiton to calculate radiance at each wavelength for current temp
   blackbodyRadiance = PLANCK_EQ(landsatWavelength, temperature)

   ;convert to W/cm^2 sr micron to match modtran units    
   blackbodyRadiance = blackbodyRadiance/(100^2)
      
   ;multiply the caluclated planck radiance by the spectral reponse and integrate over wavelength
   ;to get one number for current temp
   product = blackbodyRadiance*landsatResponse
   tempIntegral = INT_TABULATED(landsatWavelength, product)
      
   ;divide above result by integral of spectral response function
   radianceResult = tempIntegral/RSRintegral
   
   RETURN, radianceResult

END


