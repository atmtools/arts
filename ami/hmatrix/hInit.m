%------------------------------------------------------------------------
% NAME:    hInit
%
%          Init the H matrix and associated variables.
%          The H and Hd matrices are set to 1.
%          F_SENSOR is set to F_MONO
%          ZA_SENSOR is set to ZA_PENCIL
%          The vectors F_Y and ZA_Y have the same length as Y, the vector
%          containing the appended spectra. The elements of F_Y and ZA_Y
%          are the frequency and zenith angle of the corresponding value
%          of Y. These vectors are accordingle set to
%            F_Y = [f(1),f(2),...,f(nf),f(1),f(2),...,f(nf),...]
%            ZA_Y = [za(1),za(1),...,za(1),za(2),za(2),...,za(2),...]
%          where f(1) is the first element of F_MONO etc., nf the length
%          of F_MONO and za(1) the first element of ZA_PENCIL etc.
%
% FORMAT:  [H,f_y,za_y,Hd,f_sensor,za_sensor] = hInit(f_mono,za_pencil)
%
% RETURN:  H            total H matrix
%          f_y          frequency vector        
%          za_y         zenith angle vector
%          Hd           H matrix for data reduction
%          f_sensor     sensor frequencies
%          za_sensor    sensor zenith angles
% IN:      f_mono       monochromatic frequencies
%          za_pencil    pencil beam zenith angles
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,Hd,f_sensor,za_sensor] = hInit(f_mono,za_pencil)


H          = 1;
Hd         = 1;
f_sensor   = f_mono;
za_sensor  = za_pencil;

[f_y,za_y] = h_fix_ys(f_sensor,za_sensor);
