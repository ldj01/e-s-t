;Monica Cook
;26 September 2011

;
;Given array of specific humidities, temperature, and pressure, generate array of relative humidities
;

FUNCTION CONVERT_SH_RH, specHum, tempK, pressure

   ;convert input variables
   tempC = tempK - 273.15
   TC = tempC
   TK = tempK
   q = specHum
   p = pressure

   ;given constants
   NL = 6.0221415e23 ;mil^-1
   R = 8.31447215 ;Jmol^-1K^-1
   Mh20 = 18.01534 ;gmol^-1
   Mdry = 28.9644 ;gmol^-1
   
   a0w = 6.107799961
   a1w = 4.436518521e-1 
   a2w = 1.428945805e-2
   a3w = 2.650648471e-4
   a4w = 3.0312403963e-6
   a5w = 2.034080948e-8
   a6w = 6.136820929e-11
   
   ;calculate vapor pressure at given temperature
   ewater = a0w + TC*(a1w +TC*(a2w + TC*(a3w + TC*(a4w + TC*(a5w + TC*(a6w*TC)))))) ;hPa
   
   loge2 = -0.58002206e4/TK + 0.13914993-0.48640239e-1*TK+0.41764768e-4*(TK^2)-0.14452093e-7*(TK^3)+0.65459673*ALOG(TK)
   e2 = EXP(loge2) ;Pa
   
   goff = -7.90298D*(373.16/TK-1)+5.02808*ALOG10(373.16/TK)-1.3816e-7*(10^(11.344*(1-(TK/373.16)))-1)+8.1328e-3*(10^(-3.49149*(373.16/TK-1))-1)+ALOG10(1013.246) ;hPa
   ew = 10^goff
   
   ;calculate volume mixing ratio
   Xh20 = q*Mdry/(Mh20-q*Mh20+q*mdry)
   
   ;calculate partial pressure
   Ph20 = Xh20*p 
   
   ;calculate relative humidity
   RH = (Ph20/ew)*100   

   ;two alternative equations for calculating dew point temperature if later desired
;   TD = (RH/100D)^(1/8D)*(112+0.9*TC)+0.1*TC-112
;   TD2 = (243.5*ALOG(RH/100)+17.67*(TC/(243.5+TC)))/(17.67-ALOG(RH/100)+17.67*(TC/(243.5+TC)))
;   
   ;converting dew point temperature to kelvin   
;   dewPoint = TD + 273.15
;   dewPoint2 = TD2 + 273.15
       
   ;return array of corresponding relative humidities
   RETURN, RH

END
