%------------------------------------------------------------------------
% NAME:    hAntennaFromFileAdv
%
%          Includes an antenna pattern, specified in a file, into H with 
%          the options to scale the pattern with frequency and to consider 
%          a moving antenna.
%          The antenna file shall have ARTS format, be a 2 column matrix
%          where column 1 is (relative) zenith angles and column 2
%          antenna pattern values.
%          The response of the antenna pattern is normalised and the
%          antenna values do not need to be normalised.
%
% FORMAT:  [H,f_y,za_y,za_sensor] = hAntennaFromFileAdv(H,f_sensor,za_sensor,
%                                za_obs,filename,o_ant,o_y,fscale,f0,move,dza)
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
%          o_ant       linear (=1) or cubic (=3) treatment of the antenna 
%                      pattern
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%          fscale      flag to scale the pattern with frequency
%          f0          reference frequency for frequency scaling, i.e. for 
%                      which frequency FWHM is valid
%          move        flag to consider a moving antenna with a constant
%                      scanning velocity during the integration
%          dza         total movement during the integration [deg]
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 
%          00.11.16  Included linear/cubic flags (PE) 


function [H,f_y,za_y,za_sensor] = hAntennaFromFileAdv(H,f_sensor,za_sensor,...
                                 za_obs,filename,o_ant,o_y,fscale,f0,move,dza)


%=== Read the antenna file
A = read_datafile( filename, 'MATRIX' );
 

%=== Get H for the antenna pattern
[Hant,za_sensor] = h_antenna(f_sensor,za_sensor,za_obs,A(:,1),A(:,2),...
                                                o_ant,o_y,fscale,f0,move,dza);


%=== Include Hant in H
H = h_x_h( Hant, H );


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys( f_sensor, za_sensor );
