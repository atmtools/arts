%------------------------------------------------------------------------
% NAME:    hBackendGauss
%
%          Includes into H a backend with Gaussian channel responses.
%          All channels are assumed to have the same response.
%          The response of the channels is normalised and the
%          response values do not need to be normalised.
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hBackendGauss(H,f_sensor,za_sensor,
%                                        f_centre,fwhm,width,spacing,o_ch,o_y)
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
%          width       the (total) width of the response to consider [Hz] 
%          spacing     maximum spacing of the abscissa to use for the 
%                      channel response [Hz]
%          o_ch        linear (=1) or cubic (=3) treatment of the channel
%                      response
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 
%          00.11.16  Included linear/cubic flags (PE) 


function [H,f_y,za_y,f_sensor] = ...
      hBackendGauss(H,f_sensor,za_sensor,f_centre,fwhm,width,spacing,o_ch,o_y)


%=== Convert FWHM and WIDTH to a standard deviation values
si     = fwhm/sqrt(2*log(2))/2;
nsi    = width/sqrt(2*log(2))/2;


%=== Set up grid for channel response
npoints = 2*ceil(width/spacing) + 1;
f_back  = linspace(-nsi,nsi,npoints);


%=== Calculate channel response at F_BACK
back   = exp(-1*f_back.^2/(2*si*si))/si/sqrt(2*pi);


%=== Get H for the backend
[Hback,f_sensor] = h_backend(f_sensor,f_centre,za_sensor,f_back,back,o_ch,o_y);


%=== Include Hback in H
H      = h_x_h(Hback,H);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor);
