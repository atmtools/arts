%------------------------------------------------------------------------
% NAME:    pbasis_x_pol
%
%          Gives the polynomial coefficient basis after multiplication
%          with a known polynomial.
%
%          The coefficients in B1 are of the type returned by PBASIS, while
%          the coefficients in P are of the type returned by POLYFIT. 
%
% FORMAT:  b2 = pbasis_x_pol(b1,p)
%
% RETURN:  b2          New polynomial coefficient basis
% IN:      b1          Old polynomial coefficient basis
%          p           Polynomial coefficients (from e.g. POLYFIT)
%------------------------------------------------------------------------

% HISTORY: 00.11.15  Created by Patrick Eriksson (PE).

% This function does not use the full potential of Matlab for consistency
% reasons when it will be ported to ARTS.


function b2 = pbasis_x_pol(b1,p)

nx = size(b1,1);
n1 = size(b1,2);
np = length(p);
n2 = n1 + np - 1;

b2 = zeros(nx,n2);

for i1 = 1:n1
  for ip = 1:np
    for ix = 1:nx
      b2(ix,i1+ip-1) = b2(ix,i1+ip-1) + b1(ix,i1)*p(ip);
    end
  end
end
