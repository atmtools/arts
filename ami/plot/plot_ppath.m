%-----------------------------------------------------------------------------
% NAME:     plot_ppath
%
%           Plots a propagation path.
%
%           The function handles so far only 2D and is suitable to use with 
%           plot_atm_2d.
%
%           The length unit of the figure is assumed to be kilometers.
%
% FORMAT:   plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid )
%
% OUT:      h             Handle to each used plotting symbols. 
%           ltext         Text suitable for a legend. A legend can be added
%                         to the figure by: legend(h,ltext);
% IN:       r_s           Radius of sensor position.
%           alpha_s       Latitude of sensor position.
%           ppath         Propagation path structure.
%           alpha_grid    Latitude grid.
%           r_geoid       Geoid radius.
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function [h,ltext] = plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid )


% Length unit is set to km
lscale = 1e3;
lunit  = 'km';


% Plot sensor
[x1,y1] = cyl2cart( r_s, alpha_s );
h(1) = plot( x1/lscale, y1/lscale, 'r*' );  
ltext{1} = 'sensor position';


% Plot ppath
for i = 1:ppath.np
  r = interp1( alpha_grid, r_geoid, ppath.pos(i,2) );
  [x2,y2] = cyl2cart( r+ppath.z(i), ppath.pos(i,2) );
  h(2) = plot( x2/lscale, y2/lscale, 'ro' );  
  if i > 1
    plot( [x1,x2]/lscale, [y1,y2]/lscale, 'r-' );
  end
  x1 = x2;
  y1 = y2;
end
ltext{2} = 'path points';



% Tangent point?
if( ~isempty(ppath.tan_pos) )
  r = interp1( alpha_grid, r_geoid, ppath.tan_pos(2) );
  [x1,y1] = cyl2cart( r+ppath.tan_pos(1), ppath.tan_pos(2) );
  h(3) = plot( x1/lscale, y1/lscale, 'r+' ); 
  ltext{3} = 'tangent point';
end
