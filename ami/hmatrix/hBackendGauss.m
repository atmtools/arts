%------------------------------------------------------------------------
% NAME:    hBackendGauss
%
%          Includes into H a backend with Gaussian channel responses.
%          All channels are assumed to have the same response.
%          The response of the channels is normalised and the
%          response values do not need to be normalised.
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hBackendGauss(H,f_sensor,za_sensor,
%                                                f_centre,fwhm,width,npoints)
%
% RETURN:  H           H matrix after backend
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequency grid
% IN:      H           H matrix before the backend
%          f_sensor    input frequency grid
%          za_sensor   zenith angles
%          f_centre    centre frequencies of the channels
%          fwhm        full width at half mean of the channel response [Hz]
%          width       the width of the response to consider [Hz] 
%          npoints     number of values to use for definition of the response.
%                      NPOINTS must be an odd number > 2
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,f_sensor] = ...
              hBackendGauss(H,f_sensor,za_sensor,f_centre,fwhm,width,npoints)


%=== Check input
if npoints < 3
  error('The number of points must be >= 3');
end
if ~isodd(npoints)
  error('The number of points must be an odd number');
end


%=== Convert FWHM and WIDTH to a standard deviation values
si     = fwhm/sqrt(2*log(2))/2;
nsi    = width/sqrt(2*log(2))/2;


%=== Set up grid for channel response
f_back = linspace(-nsi,nsi,npoints);


%=== Calculate channel response at F_BACK
back   = exp(-1*f_back.^2/(2*si*si))/si/sqrt(2*pi);


%=== Get H for the backend
[Hback,f_sensor] = h_backend(f_sensor,f_centre,za_sensor,f_back,back);


%=== Include Hback in H
H      = h_x_h(Hback,H);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor)
