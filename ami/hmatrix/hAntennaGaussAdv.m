%------------------------------------------------------------------------
% NAME:    hAntennaGaussAdv
%
%          Includes a Gaussian antenna pattern into H with options to
%          scale the pattern with frequency and to consider a moving antenna.
%
% FORMAT:  [H,f_y,za_y,za_sensor] = hAntennaGaussAdv(H,f_sensor,za_sensor,
%                               za_obs,fwhm,width,npoints,fscale,f0,move,dza)
%
% RETURN:  H           H matrix after antenna
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          za_sensor   new zenith angles, set to ZA_OBS
% IN:      H           H matrix before the antenna
%          f_sensor    frequencies
%          za_sensor   input zenith angles
%          za_obs      zenith angles observed by the sensor
%          fwhm        full width at half mean of the antenna pattern [degs]
%          width       width of the antenna pattern to consider [degs] 
%          npoints     number of values to use for definition of the pattern.
%                      NPOINTS must be an odd number
%          fscale      flag to scale the pattern with frequency
%          f0          reference frequency for frequency scaling, i.e. for 
%                      which frequency FWHM is valid
%          move        flag to consider a moving antenna with a constant
%                      scanning velocity during the integration
%          dza         total movement during the integration [deg]
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,za_sensor] = hAntennaGaussAdv(...
           H,f_sensor,za_sensor,za_obs,fwhm,width,npoints,fscale,f0,move,dza)


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


%=== Set up zenith angle grid for antenna pattern
za_ant = linspace(-nsi,nsi,npoints);


%=== Calculate antenna pattern at ZA_ANT
ant    = exp(-1*za_ant.^2/(2*si*si))/si/sqrt(2*pi);


%=== Get H for the antenna pattern
[Hant,za_new] = ...
         h_antenna(f_sensor,za_sensor,za_obs,za_ant,ant,fscale,f0,move,dza);


%=== Include Hant in H
H = h_x_h(Hant,H);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor)
