%------------------------------------------------------------------------
% NAME:     fill_grid
%
%           Ensures that a grid is not too coarse. Where grid points are
%           seperated more than MAXSPACING, values are added to the grid.
%
% FORMAT:   g = fill_grid( grid, maxspacing [, islog] )
%
% OUT:      g            Final grid.
% IN:       grid         Original grid.
%           maxspacing   Maximium spacing between points. If the grid
%                        is logarithmic (ISLOG=1), MAXSPACING is in
%                        the unit of decades (base 10).
% OPTIONAL: islog        Flag to indicate log grid. 
%------------------------------------------------------------------------

% HISTORY: 2000.12.22  Created by Patrick Eriksson.


function g = fill_grid( grid, maxspacing, islog )


%=== Set default values
if ~exist('islog'),   islog = 0;    end


%=== Create vector with new values
gnew = [];
for i = 1:(length(grid)-1)

  if ~islog
    dg = abs( grid(i) - grid(i+1) );
  else
    dg = abs( log10(grid(i)) - log10(grid(i+1)) );
  end

  if dg > maxspacing
    if ~islog
      gtmp = linspace( grid(i), grid(i+1),...
                           ceil( abs(grid(i+1)-grid(i)) /maxspacing ) + 1 );
      gnew = [gnew;gtmp(2:(length(gtmp)-1))'];
    else
      gtmp = logspace(log10(grid(i)), log10(grid(i+1)),...   
             ceil( abs(log10(grid(i+1))-log10(grid(i))) /maxspacing ) + 1 );
      gnew = [gnew;gtmp(2:(length(gtmp)-1))'];
    end
  end

end


%=== Join the two grids
if ~isempty( gnew )

  grids{1} = grid;
  grids{2} = gnew;

  g = join_grids( grids, 0, islog, grid(2) < grid(1) ); 

end




