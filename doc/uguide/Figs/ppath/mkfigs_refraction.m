function mkfigs_refraction


%- Define a 2D atmosphere
%
figure(1)
clf
%
nz       = 81;
nlat     = 21;
%
dim      = 2;
lat_grid = linspace(70,110,nlat);
lon_grid = [];
p_grid   = logspace( 5, 0, nz );
z_field  = linspace( 0, 80e3, nz )' * ones(1,nlat);
r_geoid  = 6378e3 * ones(nlat,1);
z_ground = 0 * ones(nlat,1);
%
a_pos = [ (6378+90)*1e3 70];
a_los = 99.42;
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, -1, z_ground, [] );
plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, '--' );
%
hold on
grid
xlabel('Latitude distance [degrees]');
ylabel('Altitude [km]');
set_labels( gca, 'FontSize', 12, 'FontWeight', 'bold' );
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
             a_pos, a_los, -1, z_ground, [], 2e3, 200+z_field*0, 0+z_field*0 );
plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, '-' );
%
a_los = 99
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, -1, z_ground, [] );
plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, 'r--' );
%
hold on
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
             a_pos, a_los, -1, z_ground, [], 2e3, 200+z_field*0, 0+z_field*0 );
plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, 'r-' );





if yes_or_no('Print figures')
  print ppath_refr1.eps -depsc
  ! epstopdf ppath_refr1.eps
  ! rm ppath_refr1.eps
end

 



