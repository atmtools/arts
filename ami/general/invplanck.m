%------------------------------------------------------------------------
% NAME:    invplanck
% 
%          Calculates brigthness temperatures by "inverting" the
%          Planck function.
%
%          The frequency and intensity vectors/matrices must have the 
%          same size.
%
% FORMAT:  tb = invplanck(v,i)

% RETURN:  tb      brightness temperatures [K]
% IN:      v       frequencies [Hz]
%          i       intensities [W/m^2Hz]
%=====================================================

% HISTORY: 16.11.00  Created by Patrick Eriksson. 


function tb = invplanck(v,i)

global SPEED_OF_LIGHT BOLTZMAN_CONST PLANCK_CONST

tb = PLANCK_CONST * v / BOLTZMAN_CONST ./ ...
            log( 2* PLANCK_CONST * v.^3 / SPEED_OF_LIGHT^2 ./ i + 1 );

