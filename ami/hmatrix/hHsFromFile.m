%------------------------------------------------------------------------
% NAME:    hHsFromFile
%
%          Loads a sensor H matrix and includes this in the total H.
%          With other words, the loaded H matrix shall model one, or
%          several sensor features.
%          The H matrix shall be saved with hSave (as H, the Hd matrix is
%          ignored). The saved F_Y, ZA_Y, F_SENSOR and ZA_SENSOR are
%          returned.
%
% FORMAT:  [H,f_y,za_y,f_sensor,za_sensor] = hHsFromFile(H,basename)
%
% RETURN:  H           new total H matrix
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequencies
%          za_sensor   new zenith angles
% IN:      H           old total H matrix
%          basename    basename for H file 
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,f_sensor,za_sensor] = hHsFromFile(H,basename)


%=== Load the H matrix and associated grids
[Hnew,f_y,za_y,U,f_sensor,za_sensor] = hLoad(basename);


%=== Perform multiplication
H    = h_x_h(Hnew,H);



