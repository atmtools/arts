%-----------------------------------------------------------------------------
% NAME:     plot_1gridbox_2d
%
%           Plots a 2D grid box.
%
%           The present figure is cleared and the figure is left with hold on.
%
%           The length unit of the figure is set to kilometers.
%
% FORMAT:   plot_1gridbox_2d( r1, r2, r3, r4, alpha1, alpha2, rground1, 
%                                                                   rground2 )
%
% IN:       r1         Radius of the box at the lower left corner.
%           r2         Radius of the box at the lower right corner.
%           r3         Radius of the box at the upper right corner.
%           r4         Radius of the box at the upper left corner.
%           alpha1     Lower latitude boundary of the grid box.
%           alpha2     Upper latitude boundary of the grid box.
%           rground1   Ground radius at ALPHA1.
%           rground2   Ground radius at ALPHA2.
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function plot_1gridbox_2d( r1, r2, r3, r4, alpha1, alpha2, rground1, rground2 )


% Number of points for plotting of a pressure surface
nps = max([ 3, ceil(alpha2-alpha1)*2 ]);


% Length unit is set to km
lscale = 1e3;
lunit  = 'km';


% Clear figure
clf;


% Plot vertical limits
[x1,y1] = cyl2cart( r1, alpha1 );
[x2,y2] = cyl2cart( r2, alpha2 );
[x3,y3] = cyl2cart( r3, alpha2 );
[x4,y4] = cyl2cart( r4, alpha1 );
%
plot( [x1,x4]/lscale, [y1,y4]/lscale, 'k:' );
hold on
plot( [x2,x3]/lscale, [y2,y3]/lscale, 'k:' );


% Plot pressure surfaces
plot_psurface( alpha1, alpha2, r1/lscale, r2/lscale, nps, 'k:' )
plot_psurface( alpha1, alpha2, r4/lscale, r3/lscale, nps, 'k:' )


% Plot the ground
if( exist('rground1','var') ) 
  plot_psurface( alpha1, alpha2, rground1/lscale, rground2/lscale, nps, ...
                                                                    'b-', 1.5 )
end


xlabel(['x [',lunit,']'])
ylabel(['y [',lunit,']'])



