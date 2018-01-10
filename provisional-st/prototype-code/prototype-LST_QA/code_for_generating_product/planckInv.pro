function planckInv, L, lambda1
;% Calculate the temperature from the inverse Planck function 
;%  given the spectral radiance and wavelength.
;% 
;%   lambda1 (input):  wavelength at which to evaluate [microns]
;%   L      (input):  spectral radiance at given wavelength 
;%                     [watts per meter squared per steradian per micron] 
;%   T     (output):  apparent temperature of radiance [Kelvin]
;
;% NOTE: L must be in terms of [W/m2/sr/um]

lambda1 = lambda1 * 1e-6;    % convert wavelength to meters (from microns)
L = L * 1e6;               % convert radiance to W/m2/sr/m (from W/m2/sr/um)

;% physical constants:
c = 2.99792458e8;     % speed of light [m/s]
h = 6.6260755e-34;    % Planck's constant [J s]
k = 1.3806503e-23;     % Boltzmann constant [J/K]0658

T = 2 * h * c * c / (L * (lambda1^5));
T = alog(T + 1);
T = (h * c / (k * lambda1)) / T;
return,T

end
