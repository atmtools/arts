%------------------------------------------------------------------------
% NAME:    hAntennaGauss
%
%          Includes a Gaussian antenna pattern into H.
%          
% FORMAT:  [H,f_y,za_y,za_sensor] = hAntennaGauss(H,f_sensor,za_sensor,
%                                                  za_obs,fwhm,width,npoints)
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
%          npoints     number of values to use for definition of the pattern
%                      NPOINTS must be an odd number
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,za_sensor] = ...
               hAntennaGauss(H,f_sensor,za_sensor,za_obs,fwhm,width,npoints)


[H,f_y,za_y,za_sensor] = ...
     hAntennaGaussAdv(H,f_sensor,za_sensor,za_obs,fwhm,width,npoints,0,0,0,0);





