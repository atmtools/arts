%------------------------------------------------------------------------
% NAME:    hAntennaGauss
%
%          Includes a Gaussian antenna pattern into H.
%          
% FORMAT:  [H,f_y,za_y,za_sensor] = hAntennaGauss(H,f_sensor,za_sensor,
%                                         za_obs,fwhm,width,spacing,o_ant,o_y)
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
%          width       total width of the antenna pattern to consider [degs] 
%          spacing     maximum spacing of the abscissa to use for the 
%                      antenna pattern [deg]
%          o_ant       linear (=1) or cubic (=3) treatment of the antenna 
%                      pattern
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 
%          00.11.16  Included linear/cubic flags (PE) 


function [H,f_y,za_y,za_sensor] = ...
       hAntennaGauss(H,f_sensor,za_sensor,za_obs,fwhm,width,spacing,o_ant,o_y)


[H,f_y,za_y,za_sensor] = hAntennaGaussAdv(H,f_sensor,za_sensor,za_obs,...
                                         fwhm,width,spacing,o_ant,o_y,0,0,0,0);





