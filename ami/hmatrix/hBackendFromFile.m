%------------------------------------------------------------------------
% NAME:    hBackendGauss
%
%          Includes into H a backend where the channel response is defined 
%          in a file. All channels are assumed to have the same response.
%          The channel file shall have ARTS format, be a 2 column matrix
%          where column 1 is (relative) frequencies and column 2
%          resonse values.
%          The response of the channels is normalised and the
%          response values do not need to be normalised.
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hBackendFromFile(H,f_sensor,za_sensor,
%                                                          f_centre,filename)
%
% RETURN:  H           H matrix after backend
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequency grid
% IN:      H           H matrix before the backend
%          f_sensor    input frequency grid
%          za_sensor   zenith angles
%          f_centre    centre frequencies of the channels
%          filename    name on file with channel response specification
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,f_sensor] = ...
                     hBackendFromFile(H,f_sensor,za_sensor,f_centre,filename)


%=== Read the file defining the channel response
A = read_datafile(filename);


%=== Get H for the backend
[Hback,f_sensor] = h_backend(f_sensor,f_centre,za_sensor,A(:,1),A(:,2));


%=== Include Hback in H
H = h_x_h(Hback,H);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor);
