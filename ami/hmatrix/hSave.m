%------------------------------------------------------------------------
% NAME:    hSave
%
%          Saves a H matrix and associated varaibles in Matlab format
%          The file gets the extension .H.mat
%
% FORMAT:  hSave(H,f_y,za_y,Hd,f_sensor,za_sensor,basename)
%
% RETURN:  -
% IN:      H           total H matrix
%          f_y         frequency vector
%          za_y        zenith angle vector 
%          Hd          data reduction H matrix
%          f_sensor    sensor frequencies
%          za_sensor   sensor zenith angles 
%          basename    ARTS basename
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson. 


function hSave(H,f_y,za_y,Hd,f_sensor,za_sensor,basename)


name = sprintf('%s.H.mat',basename);


save(name,'H','f_y','za_y','Hd','f_sensor','za_sensor');
