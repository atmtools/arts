% plot_result( [do_1d] )
%
% A simple function to plot the simulation results. Set do_1d to plot
% results from the special 1D test case.

function plot_result(do_1d)

if ~nargin  |  ~do_1d
  f    = xmlLoad( 'TestOdinSMR.sensor_response_f.xml' ); 
  ztan = xmlLoad( 'TestOdinSMR.ztan.xml' );
  y    = xmlLoad( 'TestOdinSMR.y.xml' );
else
  f    = xmlLoad( 'TestOdinSMR_1D.sensor_response_f.xml' ); 
  ztan = xmlLoad( 'TestOdinSMR_1D.ztan.xml' );
  y    = xmlLoad( 'TestOdinSMR_1D.y.xml' );
end

Y = reshape( y, length(f), length(ztan) );

h = plot( f/1e9, Y(:,end:-1:1), '-' );

xlabel( 'Frequency [GHz]' );
ylabel( 'RJ Tb [K]' );


L = {};
%
for i = 1 :length(ztan)
  L{i} = sprintf( '%.1f km', ztan(i)/1e3 );
end
%
legend( L{end:-1:1} );
