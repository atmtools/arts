%------------------------------------------------------------------------
% NAME:    pbasis_integrate
%
%          Performs an integration of a polynomial basis (see PBASIS).
%
%          If the function G has the same abscissa as the polynomial
%          basis, the integral of G between XA and XB can be calculated
%          as 
%             W*G
%          where W is the output of this function and G is assumed to be
%          a column vector.
%          
% FORMAT:  w = pbasis_integrate(b,xa,xb)
%
% RETURN:  w           Integration weights for each point of the basis.
% IN:      b           Polynomial coefficient basis.
%          xa          Lower integration limit.
%          xb          Upper integration limit.
%------------------------------------------------------------------------

% HISTORY: 00.11.15  Created by Patrick Eriksson (PE).

% This function does not use the full potential of Matlab for consistency
% reasons when it will be ported to ARTS.


function w = pbasis_integrate(b,xa,xb)

nx    = size(b,1);
order = size(b,2);

w = zeros(1,nx);

for i = order:-1:1
  a   = (xb^i-xa^i)/i;
  ib  = order - i + 1;
  for ix = 1:nx
    w(ix) = w(ix) + a * b(ix,ib);
  end
end
