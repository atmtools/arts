%------------------------------------------------------------------------
% NAME:     join_grids
%
%           Joins and sorts grids. Points that are closer than the specified
%           min spacing are removed. The grids can be logarithmical ans/or
%           in reversed order.
%
% FORMAT:   g = join_grids( grids, minspacing [, islog, reverse] )
%
% OUT:      g            Final grid.
% IN:       grids        Cell array with grids.
%           minspacing   Minimium spacing between points. If the grids
%                        are logarithmic (ISLOG=1), MINSPACING is in
%                        the unit of decades (base 10).
% OPTIONAL: islog        Flag to indicate log grids. 
%           reverse      Flag to indicate that the grids are in reversed
%                        order. Should be the normal choice with ISLOG=1.
%------------------------------------------------------------------------

% HISTORY: 2000.12.22  Created by Patrick Eriksson.


function g = join_grids( grids, minspacing, islog, reverse )


%=== Set default values
if ~exist('islog'),   islog = 0;    end
if ~exist('reverse'), reverse = 0;  end


%=== Join the grids
g = [];
for i = 1:length(grids)
  g = [g;vec2col(grids{i})];
end


%=== Sort
if ~reverse
  g = sort(g);
else
  g = sort(g);
  g = g(length(g):-1:1);
end


%=== Remove points that are too close
i  = 2;
while i <= length(g)

  if ~islog
    dg = abs( g(i) - g(i-1) );
  else
    dg = abs( log10(g(i)) - log10(g(i-1)) );
  end

  if dg < (minspacing-10*eps)
    g = g([1:(i-1),(i+1):length(g)]);
  else
    i = i + 1;
  end

end





