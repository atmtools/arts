function mkfigs_ppath


%- Define a 2D atmosphere
%
nz       = 9;
nlat     = 9;
%
dim      = 2;
lat_grid = linspace(70,110,nlat);
lon_grid = [];
p_grid   = logspace( 5, 4, nz );
z_field  = linspace( 1, 17, nz )' * ones(1,nlat);
r_geoid  = 50 * ones(nlat,1);
z_ground = 1.5;
cb_lims  = [];
%
figure(1)
clf
%
a_pos = [ 51.5 108 ];
a_los = -75;
%
n = 8;
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, -1, z_ground, cb_lims );
atmplot_pol2cart( P.pos(1,1), P.pos(1,2), 'k*' );
hold on
atmplot_pol2cart( P.pos(1:n,1), P.pos(1:n,2), 'm-o' );
%
for i = 2 : n
  ip   = 1 + P.gp_p{i}.idx;
  ilat = 1 + P.gp_lat{i}.idx;
  
  arts_plot_atmgrids( dim, lat_grid(ilat+(0:1)), lon_grid, ...
                    z_field(ip+(0:1),ilat+(0:1)), r_geoid(ilat+(0:1)), [] );
end
%
axis equal
axis tight
axes_frame( gca, 'off' )




figure(2)
clf
%
a_pos = [ 65.5 108 ];
a_los = -102.1;
%
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, 3, z_ground, cb_lims );
[x,y] = atmplot_pol2cart( P.pos(1,1), P.pos(1,2) );
h(1) = plot( x, y, 'k*' );
hold on
[x,y] = atmplot_pol2cart( P.pos(:,1), P.pos(:,2) );
h(2) = plot( x, y, 'm-o' );
[x,y] = atmplot_pol2cart( P.tan_pos(1), P.tan_pos(2) );
h(3) = plot( x, y, 'k+' );
%
for i = 2 : P.np
  ip   = 1 + P.gp_p{i}.idx;
  ilat = 1 + P.gp_lat{i}.idx;
  
  hg = arts_plot_atmgrids( dim, lat_grid(ilat+(0:1)), lon_grid, ...
                    z_field(ip+(0:1),ilat+(0:1)), r_geoid(ilat+(0:1)), [] );
end
%
axis equal
axis tight
axes_frame( gca, 'off' )
%
legend( [hg(1), h], 'Grid cells','Sensor position','Propagation path','Tangent point');



figure(3)
clf
%
lat_grid = [85 95];
z_field  = [10 10.2;11 10.9];
%
arts_plot_atmgrids( 2, lat_grid, [], z_field, [0 0]', [] );
hold on
axis equal
axis tight
axes_frame( gca, 'off' )
%
r0   = mean(z_field(:,1));
lat0 = lat_grid(1);
%
r1   = r0+0.1;
lat1 = lat_grid(2);
atmplot_pol2cart( [r0 r1], [lat0 lat1], 'mo-' );
%
lat1 = mean(lat_grid);
r1   = interp1( lat_grid, z_field(1,:), lat1 );
atmplot_pol2cart( [r0 r1], [lat0 lat1], 'mo-' );
%
lat1 = mean(lat_grid)+2;
r1   = interp1( lat_grid, z_field(2,:), lat1 );
atmplot_pol2cart( [r0 r1], [lat0 lat1], 'mo-' );
%
lat0 = lat_grid(2) - 3;%
r0   = interp1( lat_grid, z_field(1,:), lat0 );
%
r1   = z_field(1,2)+0.3;
lat1 = lat_grid(2);
atmplot_pol2cart( [r0 r1], [lat0 lat1], 'b--o' );
%
lat1 = lat0+1;
r1   = interp1( lat_grid, z_field(2,:), lat1 );
atmplot_pol2cart( [r0 r1], [lat0 lat1], 'b--o' );
%
r1   = z_field(2,1)-0.3;
lat1 = lat_grid(1);
atmplot_pol2cart( [r0 r1], [lat0 lat1], 'b--o' );



figure(4)
clf
%
lat_grid = [-2 2];
lon_grid = [-2 2];
z_field  = repmat( [10;11], [1 2 2] );
%
arts_plot_atmgrids( 3, lat_grid, lon_grid, z_field, [0 0;0 0], [] );
hold on
%
r0   = 10;
lat0 = 0;
lon0 = 0;
%
atmplot_sph2cart( r0, lat0, lon0, 'k*' );
%
r1   = 10.8;
lat1 = 0;
lon1 = 0.2;
%
atmplot_sph2cart( [r0 r1], [lat0 lat1], [lon0 lon1], 'k+' );
[x,y,z] = atmplot_sph2cart( r1, lat1, lon1 );
h=text( z, x, y+.05, 'l_{in}' );
%
r1   = 12.4;
lat1 = 0;
lon1 = 0.6;
%
atmplot_sph2cart( [r0 r1], [lat0 lat1], [lon0 lon1], 'k+' );
atmplot_sph2cart( [r0 r1], [lat0 lat1], [lon0 lon1], 'm-' );
[x,y,z] = atmplot_sph2cart( r1, lat1, lon1 );
h=text( z, x, y+.05, 'l_{out}^i' );
%
r1   = 11.6;
lat1 = 0;
lon1 = 0.4;
%
atmplot_sph2cart( [r0 r1], [lat0 lat1], [lon0 lon1], 'k+' );
[x,y,z] = atmplot_sph2cart( r1, lat1, lon1 );
h=text( z, x, y+.05, 'l_{out}^{i+1}' );
%
r1   = 11.2;
lat1 = 0;
lon1 = 0.3;
%
atmplot_sph2cart( [r0 r1], [lat0 lat1], [lon0 lon1], 'k+' );
[x,y,z] = atmplot_sph2cart( r1, lat1, lon1 );
h=text( z, x, y+.05, 'l_{out}^{i+2}' );
%
r1   = 11;
lat1 = 0;
lon1 = 0.25;
%
atmplot_sph2cart( [r0 r1], [lat0 lat1], [lon0 lon1], 'mo' );
[x,y,z] = atmplot_sph2cart( r1, lat1, lon1 );
h=text( z, x, y+.05, 'l_{end}' );
%
axis equal
axis tight
axes_frame( gca, 'off' )
view([80 30]);





if yes_or_no('Print figures')
  print ppath_3Dsearch.eps -depsc
  ! epstopdf ppath_3Dsearch.eps
  ! rm ppath_3Dsearch.eps
  figure(3)
  print ppath_ex3.eps -depsc
  ! epstopdf ppath_ex3.eps
  ! rm ppath_ex3.eps
  figure(2)
  print ppath_ex2.eps -depsc
  ! epstopdf ppath_ex2.eps
  ! rm ppath_ex2.eps
  figure(1)
  print ppath_ex1.eps -depsc
  ! epstopdf ppath_ex1.eps
  ! rm ppath_ex1.eps
end

 



