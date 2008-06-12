% plot_result( [do_1d] )
%
% A simple function to plot the simulation results. Set do_1d to plot
% results from the special 1D test case.

function [f,Y] = plot_result(do_1d)

if ~nargin  |  ~do_1d
  f      = xmlLoad( 'TestOdinSMR.y_f.xml.generated' ); 
  f_grid = xmlLoad( 'TestOdinSMR.sensor_response_f_grid.xml.generated' ); 
  ztan   = xmlLoad( 'TestOdinSMR.ztan.xml.generated' );
  y      = xmlLoad( 'TestOdinSMR.y.xml.generated' );
else
  f      = xmlLoad( 'TestOdinSMR_1D.y_f.xml.generated' ); 
  ztan   = xmlLoad( 'TestOdinSMR_1D.ztan.xml.generated' );
  y      = xmlLoad( 'TestOdinSMR_1D.y.xml.generated' );
end

Y = reshape( y, length(f_grid), length(ztan) );

h = plot( f/1e9, y, '.', f_grid/1e9, Y(:,end:-1:1), '-' );


xlabel( 'Frequency [GHz]' );
ylabel( 'RJ Tb [K]' );


L = {};
%
for i = 1 :length(ztan)
  L{i} = sprintf( '%.1f km', ztan(i)/1e3 );
end
%
legend( h(2:end), L{end:-1:1}, 'Location', 'Best' );

whos