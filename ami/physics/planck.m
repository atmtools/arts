%------------------------------------------------------------------------
% NAME:    planck
%
%          Calculates values of the Planck function.
%
%          The frequency and vectors/matrices must either have the same
%          size, or one of them can be a scalar.           
%
% FORMAT:  i = planck(v,t)
%
% RETURN:  i      intensity [W/m2Hzsr] 
% IN:      v      frequencies [Hz]
%          t      temperatures [K]
%------------------------------------------------------------------------

% HISTORY: 16.11.00  Copied and modified from Skuld by Patrick Eriksson. 


function i = planck(v,t)

global SPEED_OF_LIGHT PLANCK_CONST BOLTZMAN_CONST


i	= PLANCK_CONST*v./(BOLTZMAN_CONST*t);
i	= 2*PLANCK_CONST*v.^3./(SPEED_OF_LIGHT*SPEED_OF_LIGHT*(exp(i)-1));
