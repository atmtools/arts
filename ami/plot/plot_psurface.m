%-----------------------------------------------------------------------------
% NAME:     plot_psurface
%
%           Plots a pressure surface.
%
% FORMAT:   plot_psurface( alpha1, alpha2, r1, r2, nps, ps )
%
% IN:       alpha1   Start latitude.
%           alpha2   End latitude.
%           r1       Radius at the start latitude.
%           r2       Radius at the end latitude.
%           nps      The surface is plotted as a number of straight lines.
%                    This argument is the number of lines to use. 
%           ps       Plotting symbol to use, for exmple, '--'.
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function plot_psurface( alpha1, alpha2, r1, r2, nps, ps )


dalpha  = alpha2 - alpha1;
c       = (r2-r1) / dalpha; 
[x1,y1] = cyl2cart( r1, alpha1 );

for k = 1:nps
  da = k * dalpha / nps;
  [x2,y2] = cyl2cart( r1+c*da, alpha1+da );   
  plot( [x1,x2], [y1,y2], ps );
  x1 = x2;
  y1 = y2;
end

