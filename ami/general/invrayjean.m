%------------------------------------------------------------------------
% NAME:    invrayjean(v,i)
% 
%          Calculates the "Rayleigh-Jean" temperature using the
%          Rayleigh-Jean approximation of Planck function.
%
%          The frequency and intensity vectors/matrices must have the 
%          same size.
%
% FORMAT:  tb = invrayjean(v,i)
%
% RETURN:  tb      brightness temperatures [K]
% IN:      v       frequencies [Hz]
%          i       intensities [W/m^2Hz]
%=====================================================

% HISTORY: 16.11.00  Copied and modified from Skuld by Patrick Eriksson. 


function tb = invrayjean(v,i)

global SPEED_OF_LIGHT BOLTZMAN_CONST

tb = (SPEED_OF_LIGHT*SPEED_OF_LIGHT*i)./(2*BOLTZMAN_CONST*v.^2);
