%------------------------------------------------------------------------
% NAME:    read_arts
%
%          Reads a ARTS spectrum and converts it to brightness
%          temperatures.
%          The mean of the lowest and highest absorption frequencies
%          is used for the conversion to brightness temperatures.
%
% FORMAT:  [y,f] = read_arts(basename)
%
% RETURN:  y         The spectrum [K]
%          f         The absorption frequencies [Hz]
% IN:      basename  The ARTS basename 
%------------------------------------------------------------------------

% HISTORY: 24.05.00  Created by Patrick Eriksson. 

function [y,f] = read_tb(basename)

global BOLTZMAN_CONST SPEED_OF_LIGHT


%=== Read frequencies
f  = read_arts(basename,'f_abs');
f0 = (f(1)+f(length(f)))/2;


%=== Read spectrum
y = read_arts(basename,'y');


%=== Convert to brightness temperatures
y = y * (SPEED_OF_LIGHT^2/f0^2/2/BOLTZMAN_CONST);
