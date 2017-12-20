;Monica Cook
;5 October

;
;Calculate observed radiance from modtran results and the spectral reponse function
;

FUNCTION CALCULATE_LOBS, wavelengths, modtran, spectralResponse

   ;split modtran results into wavelength and radiance
   modtranWavelength = wavelengths
   modtranRadiance = modtran[0,*]
   
   ;split spectral reponse into wavelength and repsonse
   landsatWavelength = spectralResponse[0,*]
   landsatResponse = spectralResponse[1,*]

   ;integrate spectral response over wavelength
   RSRintegral = INT_TABULATED(landsatWavelength, landsatResponse)

   ;inteprolate modtran radiance to landsat wavelengths
   tempRad = INTERPOL(modtranRadiance, modtranWavelength, landsatWavelength)
   
   ;multiply modtran radiance and spectral response and integrate over wavelength
   product = tempRad*landsatResponse
   tempIntegral = INT_TABULATED(landsatWavelength, product)
   
   ;divide above result by spectral response integral to generate observed radiance
   result = tempIntegral/RSRintegral
   
   RETURN, result
   
END
