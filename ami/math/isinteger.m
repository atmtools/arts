%------------------------------------------------------------------------
% NAME:     isinteger
%
%           Determines if the elements of a matrix are integers.
%
% FORMAT:   bool = isinteger(n [,neps])
%
% RETURN:   bool     1 if n integer, else 0.
% IN:       n        A matrix.
% OPTIONAL: neps     The deviation from an integer to tolerate in number
%                    of eps (see help eps). Default is 1.   
%------------------------------------------------------------------------

% HISTORY: 2000.12.18  Created by Patrick Eriksson. 


function bool = isinteger(n,neps)

if ~exist('neps'), neps = 1; end

if ~isnumeric( n )
  error('Only numeric input is allowed.');
end

bool = abs( (n-round(n)) ) <= neps*eps;
