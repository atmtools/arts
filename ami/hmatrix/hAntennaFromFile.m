%------------------------------------------------------------------------
% NAME:    hAntennaFromFile
%
%          Includes an antenna pattern, specified in a file, into H.
%          The antenna file shall have ARTS format, be a 2 column matrix
%          where column 1 is (relative) zenith angles and column 2
%          antenna pattern values.
%          The response of the antenna pattern is normalised and the
%          antenna values do not need to be normalised.          
%
% FORMAT:  [H,f_y,za_y,za_sensor] = hAntennaFromFile(H,f_sensor,za_sensor,
%                                                             za_obs,filename)
%
% RETURN:  H           H matrix after antenna
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          za_sensor   new zenith angles, set to ZA_OBS
% IN:      H           H matrix before the antenna
%          f_sensor    frequencies
%          za_sensor   input zenith angles
%          za_obs      zenith angles observed by the sensor
%          filename    name on file with antenna specification
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,za_sensor] = ...
                      hAntennaFromFile(H,f_sensor,za_sensor,za_obs,filename)


[H,f_y,za_y,za_sensor] = ...
            hAntennaFromFileAdv(H,f_sensor,za_sensor,za_obs,filename,0,0,0,0);


