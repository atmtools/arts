function mkfigs_ppath_cases



%--- Allowed paths without cloud box
%
%- Define a 2D atmosphere
%
nz       = 9;
nlat     = 9;
%
dim      = 2;
lat_grid = linspace(45,135,nlat);
lon_grid = [];
p_grid   = logspace( 5, 4, nz );
z_field  = linspace( 1, 17, nz )' * ones(1,nlat);
r_geoid  = 50 * ones(nlat,1);
z_ground = 1.5;
cb_lims  = [];
%
figure(1)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 0 );
axes_frame( gca, 'off' )
hold on
%
a_pos = [ 58 130 ];
a_los = -60;
%
plot_path( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                      a_pos, a_los, -1, z_ground, cb_lims, 1 );
%
a_los = -130;
%
plot_path( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                      a_pos, a_los, -1, z_ground, cb_lims, 1 );
%
a_pos = [ 75 55 ];
a_los = 130;
%
P = plot_path( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                      a_pos, a_los, -1, z_ground, cb_lims, 0 );
atmplot_pol2cart( a_pos(1), a_pos(2), 'k*' );
atmplot_pol2cart( [a_pos(1) P.pos(1,1)], [a_pos(2) P.pos(1,2)], 'm--' );
%
a_los = 170;
%
P = plot_path( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                      a_pos, a_los, -1, z_ground, cb_lims, 0 );
atmplot_pol2cart( [a_pos(1) P.pos(1,1)], [a_pos(2) P.pos(1,2)], 'm--' );
%
atmplot_pol2cart( [a_pos(1) a_pos(1)], [a_pos(2) a_pos(2)+45], 'm--' );
%
axis tight
%
a_pos = [ 58 130 ];
a_los = -60;
%
plot_path( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                      a_pos, a_los, -1, z_ground, cb_lims, 1 );



%--- Allowed paths with cloud box
%
%- Define a 1D atmosphere
%
nz       = 9;
%
dim      = 1;
lat_grid = [45 135];
lon_grid = [];
p_grid   = logspace( 5, 4, nz );
z_field  = linspace( 1, 17, nz )';
r_geoid  = 50;
z_ground = 1.5;
cb_lims  = {4,6};
%
figure(2)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 0 );
axes_frame( gca, 'off' )
hold on
%
a_pos = 56;
lat0  = 60; 
a_los = 90;
%
atmplot_pol2cart( a_pos(1), lat0, 'k*' );
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                         a_pos, a_los, -1, z_ground, cb_lims );
atmplot_pol2cart( P.pos(:,1), lat0+P.pos(:,2), 'm-o' );
%
a_los = 150;
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                         a_pos, a_los, -1, z_ground, cb_lims );
atmplot_pol2cart( P.pos(:,1), lat0+P.pos(:,2), 'm-o' );
%
a_pos = 66;
a_los = 110;
%
atmplot_pol2cart( a_pos(1), lat0, 'k*' );
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                         a_pos, a_los, -1, z_ground, cb_lims );
atmplot_pol2cart( P.pos(:,1), lat0+P.pos(:,2), 'm-o' );
%
a_pos = 60;
lat0  = 111; 
%
atmplot_pol2cart( a_pos(1), lat0, 'k*' );
atmplot_pol2cart( a_pos(1), lat0, 'mo' );
%
a_pos = 75;
lat0  = 90;
a_los = 175;
%
atmplot_pol2cart( a_pos(1), lat0, 'k*' );
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                         a_pos, a_los, -1, z_ground, cb_lims );
atmplot_pol2cart( P.pos(:,1), lat0+P.pos(:,2), 'm-o' );
atmplot_pol2cart( [a_pos(1) P.pos(1,1)], lat0+[0 P.pos(1,2)], 'm--' );
%
axis tight



%--- Not allowed paths
%
%- Define a 2D atmosphere
%
nz       = 9;
nlat     = 9;
%
dim      = 2;
lat_grid = linspace(60,120,nlat);
lon_grid = [];
p_grid   = logspace( 5, 4, nz );
z_field  = linspace( 1, 17, nz )' * ones(1,nlat);
r_geoid  = 50 * ones(nlat,1);
z_ground = 1.5;
cb_lims  = [];
%
figure(3)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 0 );
axes_frame( gca, 'off' )
hold on
%
atmplot_pol2cart( 70, 90, 'k*' );
atmplot_pol2cart( [70 75], [90 100], 'm--' );
%
atmplot_pol2cart( 60, 55, 'k*' );
atmplot_pol2cart( [60 60], [55 60], 'm--' );
%
atmplot_pol2cart( 46, 100, 'k*' );
atmplot_pol2cart( [46 51.5], [100 115], 'm--' );
%
atmplot_pol2cart( 74, 125, 'k*' );
atmplot_pol2cart( [74 66], [125 120], 'm--' );
%
atmplot_pol2cart( 64, 115, 'k*' );
atmplot_pol2cart( [64 60], [115 120], 'm--' );
%
atmplot_pol2cart( 53, 70, 'k*' );
atmplot_pol2cart( [53 65], [70 60], 'm--' );
%
atmplot_pol2cart( 74, 70, 'k*' );
atmplot_pol2cart( [74 70], [70 55], 'm--' );
%
axis tight





if yes_or_no('Print figures')
  print ppath_badcases.eps -depsc
  ! epstopdf ppath_badcases.eps
  ! rm ppath_badcases.eps
  figure(2)
  print ppath_cases1.eps -depsc
  ! epstopdf ppath_cases1.eps
  ! rm ppath_cases1.eps
  figure(1)
  print ppath_cases2.eps -depsc
  ! epstopdf ppath_cases2.eps
  ! rm ppath_cases2.eps
end




function P = plot_path( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                             a_pos, a_los, lmax, z_ground, cb_lims, do_sensor )
  P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, lmax, z_ground, cb_lims );
  atmplot_pol2cart( P.pos(:,1), P.pos(:,2), 'm-o' );
  if do_sensor
    atmplot_pol2cart( P.pos(1,1), P.pos(1,2), 'k*' );
  end
  if ~isempty( P.tan_pos )
    atmplot_pol2cart( P.tan_pos(1), P.tan_pos(2), 'k+' );
  end
return