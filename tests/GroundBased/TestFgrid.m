% Plots results from TestFgrid.arts

function TestFgrid
  
arts TestFgrid.arts;

x1 = xmlLoad('TestFgrid.x1.xml.generated');
y1 = xmlLoad('TestFgrid.y1.xml.generated');

x2 = xmlLoad('TestFgrid.x2.xml.generated');
ye = xmlLoad('TestFgrid.ye.xml.generated');

x3 = xmlLoad('TestFgrid.x3.xml.generated');
y3 = xmlLoad('TestFgrid.y3.xml.generated');

x5 = xmlLoad('TestFgrid.x5.xml.generated');
y5 = xmlLoad('TestFgrid.y5.xml.generated');


dyl = interp1( x1, y1, x2 ) - ye;
dy3 = interp1( x3, y3, x2 ) - ye;
dy5 = interp1( x5, y5, x2 ) - ye;

f = (x2 - 1.108360400e+11)/1e6;

plot( f, dyl, f, dy3, f, dy5 );

xlabel( 'Frequency [MHz]' );
ylabel( 'Error [K]' );

legend( 'Linear', 'P3/N2', 'P5/N4' );


