f    = xmlLoad( 'TestOdinSMR.sensor_response_f.xml' ); 
ztan = xmlLoad( 'TestOdinSMR.ztan.xml' );

y = xmlLoad( 'TestOdinSMR.y.xml' );
Y = reshape( y, length(f), length(ztan) );

h = plot( f/1e9, Y(:,end:-1:1), '-' );

xlabel( 'Frequency [GHz]' );
ylabel( 'RJ Tb [K]' );


for i = 1 :length(ztan)
  L{i} = sprintf( '%.1f km', ztan(i)/1e3 );
end
%
legend( L{end:-1:1} );