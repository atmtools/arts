%------------------------------------------------------------------------
% NAME:    hLoad
%
%          Loads a H matrix and associated varaibles stored in Matlab format
%          The file must have the extension .H.mat
%          The data is preferably stored by hSave.
%
% FORMAT:  [H,f_y,za_y,Hd,f_sensor,za_sensor,basename] = hLoad(basename)
%
% RETURN:  H           total H matrix
%          f_y         frequency vector
%          za_y        zenith angle vector 
%          Hd          data reduction H matrix
%          f_sensor    sensor frequencies
%          za_sensor   sensor zenith angles
% IN:      basename    ARTS basename 
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 

function [H,f_y,za_y,Hd,f_sensor,za_sensor,basename] = hLoad(basename)


%=== Create full file name
name = sprintf('%s.H.mat',basename);


load(name);
