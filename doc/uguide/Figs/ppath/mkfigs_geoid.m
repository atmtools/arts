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


if answer_is_yes('Print figures')
  print wgs84_latdiff.eps -depsc
  ! epstopdf wgs84_latdiff.eps
  figure(1)
  print wgs84_radii.eps -depsc
  ! epstopdf wgs84_radii.eps
end


