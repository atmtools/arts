%------------------------------------------------------------------------
% NAME:    hMixerConst
%
%          Includes into H a mixer/sideband filtering, with a constant
%          sideband ratio.
%          The sideband ratio is treated as follows:
%            Response of image band = RATIO
%            Response of primary bad = 1 - RATIO
%
% FORMAT:  [H,f_y,za_y,f_sensor] = hMixerConst(H,f_sensor,za_sensor,
%                                                       lo,fprimary,ratio,o_y )
%
% RETURN:  H           H matrix after antenna
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          f_sensor    new frequency grid 
% IN:      H           H matrix before the mixer
%          f_sensor    input frequency grid
%          za_sensor   zenith angles
%          lo          LO frequency
%          fprimary    a frequency inside the primary band (!=LO)
%          ratio       relative sideband response
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%------------------------------------------------------------------------

% HISTORY: 00.08.25  Created by Patrick Eriksson 


function [H,f_y,za_y,f_sensor] = ...
                      hMixerConst(H,f_sensor,za_sensor,lo,fprimary,ratio,o_y)


%=== Get H
nf = length(f_sensor);
if fprimary > lo
  [Hmix,f_sensor] = h_mixer(f_sensor,za_sensor,lo,fprimary,...
     [f_sensor(1) lo lo+0.1 f_sensor(nf)],[ratio ratio 1-ratio 1-ratio],1,o_y);
else
  [Hmix,f_sensor] = h_mixer(f_sensor,za_sensor,lo,fprimary,...
     [f_sensor(1) lo-0.1 lo f_sensor(nf)],[1-ratio 1-ratio ratio ratio],1,o_y);
end


%=== Include Hmix in H
H = h_x_h(Hmix,H);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f_sensor,za_sensor);
