%------------------------------------------------------------------------
% NAME:     alpha_between_r1_r2
%
%           Calculates the angular distance between two radii in polar
%           cordinate system for a given zenith angle.
%
%           The zenith angle is the angle between the zenith direction
%           and the observation direction. No distinction is made between
%           positive and negative angles, ang the returned angle is 
%           always positive. 
%
%           Maximum allowed absolute value for PSI1 is 180.
%
% FORMAT:   alpha = alpha_between_r1_r2( r1, psi1, r2 )
%
% OUT:      alpha   The angular distance (always positive).
% IN:       r1      The radius at the starting point.
%           psi     The zenith angle at the r1 point (-180 - 180).
%           r2      The radius at the end point.
%------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function alpha = alpha_between_r1_r2( r1, psi1, r2 )


if abs(psi1) > 180 
  error('The 
end

global DEG2RAD RAD2DEG

if( abs(psi1) <= 90 )
  psi  = abs(psi1);
  psi0 = psi;
  as   = -1; 
else
  psi  = 180 - abs(psi1);
  psi0 = -psi;
  as   = 1; 
end

alpha = psi0 + as * RAD2DEG * asin( sin(DEG2RAD*psi) * r1 / r2 );
