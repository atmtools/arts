%------------------------------------------------------------------------
% NAME:     isinteger
%
%           Determines if a number is an integer.
%
% FORMAT:   bool = isinteger(n [,neps])
%
% RETURN:   bool     1 if n integer, else 0.
% IN:       n        A number.
% OPTIONAL: neps     The deviation from an integer to tolerate in number
%                    of eps (see help eps). Default is 1.   
%------------------------------------------------------------------------

% HISTORY: 2000.12.18  Created by Patrick Eriksson. 


function bool = isinteger(n,neps)

if ~exist('neps'), neps = 1; end

bool = abs( (n-round(n)) ) <= neps*eps;
