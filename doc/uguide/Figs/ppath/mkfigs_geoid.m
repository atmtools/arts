function mkfigs_geoid


figure(1)
%
lat = 0:90;
%
r1 = wgs84( 1, lat, 0 );
r2 = wgs84( 2, lat );
%
plot( lat, r2/1e3, '-', lat, r1/1e3, '--' );
axis([0 90 6330 6400])
legend('ellipsiod radius','curvature radius',4);
xlabel('Latitide [deg]')
ylabel('Radius [km]')

figure(2)
%
lat0 = 0.5:89.5;
lat = 0:90;
%
r2 = wgs84( 2, lat );
%
plot( lat0, diff(r2) );
axis([0 90 -400 0])
xlabel('Latitide [deg]')
ylabel('Radius difference [m/deg]')



%- Define a 2D atmosphere
%
figure(3)
clf
%
nz       = 81;
nlat     = 21;
%
dim      = 2;
lat_grid = linspace(25,65,nlat);
lon_grid = [];
p_grid   = logspace( 5, 0, nz );
z_field  = linspace( 0, 80e3, nz )' * ones(1,nlat);
r_geoid  = wgs84(2,lat_grid');
z_ground = 0 * ones(nlat,1);
%
a_los = 99.8;
a_pos = [ (6378+90)*1e3 45-(a_los-90)];
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, -1, z_ground, [] );
plot( P.pos(:,2), P.z/1e3, 'b-' );
%
hold on
grid
axis( [ 40 50 5 30 ] )
%
xlabel('Latitude [degrees]');
ylabel('Altitude [km]');
set_labels( gca, 'FontSize', 12, 'FontWeight', 'bold' );
%
r_geoid  = wgs84( 1, P.tan_pos(2), 0 ) * ones(nlat,1);
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                           a_pos, a_los, -1, z_ground, [], 0 );
plot( P.pos(:,2), P.z/1e3, 'r--' );


return

if yes_or_no('Print figures')
  print wgs84_dz.eps -depsc
  ! epstopdf wgs84_dz.eps
  ! rm wgs84_dz.eps
  figure(2)
  print wgs84_latdiff.eps -depsc
  ! epstopdf wgs84_latdiff.eps
  ! rm wgs84_latdiff.eps
  figure(1)
  print wgs84_radii.eps -depsc
  ! epstopdf wgs84_radii.eps
  ! rm wgs84_radii.eps
end


