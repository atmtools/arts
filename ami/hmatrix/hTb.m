%------------------------------------------------------------------------
% NAME:    hMixerConst
%
%          Includes into H conversion to brightness temperatures.
%
% FORMAT:  H = hTb(H,f_sensor,za_sensor)
%
% RETURN:  H           H matrix after intensity conversion
% IN:      H           H matrix before intensity conversion
%          f_sensor    frequencies
%          za_sensor   zenith angles
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson 


function H = hTb(H,f_sensor,za_sensor)


global  SPEED_OF_LIGHT  BOLTZMAN_CONST 


%=== Main sizes
nf   = length(f_sensor);
nza  = length(za_sensor);


%=== Conversion factor for each frequency
r    = (SPEED_OF_LIGHT*SPEED_OF_LIGHT/2/BOLTZMAN_CONST) ./ f_sensor.^2;


%=== Create conversion matrix
Htb  = spdiags(repmat(vec2col(r),nza,1),0,nf*nza,nf*nza);


%=== Include Hant in H
H    = h_x_h(Htb,H);
