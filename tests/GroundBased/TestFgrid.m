% Plots results from TestFgrid.arts

function TestFgrid
  
arts TestFgrid.arts;

x1 = xmlLoad('TestFgrid.x1.xml.generated');
y1 = xmlLoad('TestFgrid.y1.xml.generated');

x2 = xmlLoad('TestFgrid.x2.xml.generated');
ye = xmlLoad('TestFgrid.ye.xml.generated');

y3 = xmlLoad('TestFgrid.y3.xml.generated');
y5 = xmlLoad('TestFgrid.y5.xml.generated');

yl = interp1( x1, y1, x2 );


f = (x2 - 1.108360400e+11)/1e6;

plot( f, yl-ye, f, y3-ye, f, y5-ye );

xlabel( 'Frequency [MHz]' );
ylabel( 'Error [K]' );

legend( 'Linear', 'poly 3', 'poly 5' );


