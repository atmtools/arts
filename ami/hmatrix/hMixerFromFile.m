%------------------------------------------------------------------------
% NAME:    hMixerFromFile
%
%          Includes into H a mixer/sideband filtering, specified in a file.
%          The sideband filter file shall have ARTS format, be a 2 
%          column matrix where column 1 is frequencies and column 2
%          sideband filter values.
%          The response of the sideband filtering is normalised and the
%          response values do not need to be normalised.
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hMixerFromFileAdv(H,f_sensor,za_sensor,
%                                                      lo,fprimary,filename)
%
% RETURN:  H           H matrix after antenna
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequency grid 
% IN:      H           H matrix after the mixer
%          f_sensor    input frequency grid
%          za_sensor   zenith angles
%          lo          LO frequency
%          fprimary    a frequency inside the primary band (!=LO)
%          filename    name on file with sideband filtering specification
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson 


function [H,f_y,za_y,f_sensor] = ...
                hMixerFromFileAdv(H,f_sensor,za_sensor,lo,fprimary,filename)


%=== Read the sideband-mixer file
A = read_datafile(filename);
 

%=== Get H
[Hmix,f_sensor] = h_mixer(f_sensor,za_sensor,lo,fprimary,A(:,1),A(:,2));


%=== Include Hant in H
H = h_x_h(Hmix,H);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor);
