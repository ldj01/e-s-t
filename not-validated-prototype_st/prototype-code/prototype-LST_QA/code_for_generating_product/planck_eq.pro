FUNCTION PLANCK_EQ, wavelengths, temp

    ;##Lambda intervals of Landsat5 spectral response locations microns ## units: m
    lambda = wavelengths*10^(-6.0)

    ;## Planck Const hecht pg, 585 ## units: Js
    h = DOUBLE(6.6260755*10.0^(-34.0))

    ;## Boltzmann Gas Const halliday et 2001 ## units: J/K
    K = DOUBLE(1.3806503*10.0^(-23.0 ))

    ;## Speed of Light ## units: m/s
    c = 299792458.0

    ;## Compute the Planck Blackbody Eq [W/m^2 sr um] ##  
    L = 2.0*h*c^(2.0) * ((10^(-6.0))*lambda^(-5.0)) * (1.0/(EXP((h*c)/(lambda * K * Temp))-1.0))
    RETURN, L
    
END
