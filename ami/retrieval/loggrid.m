%------------------------------------------------------------------------
% NAME:    loggrid
%
%          Creates a vector with picewise logarithmically equally spaced 
%          points. For example, LOGGRID([100 10 1],[3,2]) gives 
%          [100.0000 46.4159 21.5443 10.0000 3.1623 1.0000] 
%
% FORMAT:  g = loggrid( ppoints, n_pr_decade )
%
% OUT:     g             Logarithmically spaced vector.
% IN:      ppoints       Brake points for the vector.
%          n_per_decade  Number of points per decade. See the example above.
%                        The length of this vector must be the length of
%                        PPOINTS-1.
%------------------------------------------------------------------------

% HISTORY: 2000.12.21  Created by Patrick Eriksson.


function g = loggrid( ppoints, n_pr_decade )


np = length( ppoints );


%=== Check input
if np ~= ( length(n_pr_decade) + 1 )
  error('The lengths of the two given vectors do not match.')
end
%
%for i = 1:(np-1)
%  if ~isinteger( ...
%         (log10(ppoints(i+1))-log10(ppoints(i))) / (1/n_pr_decade(i)), 5 )
%    error('A consistent grid cannot be created (the ppoints are not seperated as needed for n_scaleheight).');
%  end
%end


%=== Create the grid
g = [];
for i = 1:(np-1)
  n = round(abs(log10(ppoints(i))-log10(ppoints(i+1))) * n_pr_decade(i) + 1);
  g = [g;logspace(log10(ppoints(i)),log10(ppoints(i+1)),n)'];
end


%=== Remove points duplicated
ind = find(diff(g)~=0);
g   = g([ind;length(g)]);

