function mkfigs_atm_dims


%--- 1D
%
dim      = 1;
lat_grid = [45 135];
lon_grid = [];
z_field  = (1:18)';
r_geoid  = 50;
z_ground = 1.5;
cb_lims  = {5 8};
%
figure(1)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 1 );
axes_frame( gca, 'off' )
axis tight
hl = legend( h, hs );
disp('Adjust the legend');
pause;



%--- 2D
%
nlat     = 16;
%
dim      = 2;
lat_grid = linspace(45,135,nlat);
lon_grid = [];
z_field  = (1:18)'*ones(1,nlat) - ones(18,1)*linspace(0,0.8,nlat);
r_geoid  = 50 * ones(nlat,1);
z_ground = 1.5 + 0.5*randn(nlat,1);
cb_lims  = [];
%
figure(2)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 1 );
axes_frame( gca, 'off' )
axis tight
hl = legend( h, hs );
disp('Adjust the legend');
pause;



%--- 3D cross section
%
nlat     = 16;
%
dim      = 2;
lat_grid = linspace(45,135,nlat);
lon_grid = [];
z_field  = (1:18)'*ones(1,nlat);
r_geoid  = 50 * ones(nlat,1);
z_ground = 1.5 + 0.5*randn(nlat,1);
cb_lims  = {5 8 7 12};
%
figure(3)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 1 );
axes_frame( gca, 'off' )
axis tight



%--- 3D
%
nlat     = 5;
nlon     = 5;
%
dim      = 3;
lat_grid = linspace(-30,30,nlat);
lon_grid = linspace(-30,30,nlon);
z_field  = repmat( (1:5:41)'*ones(1,nlat), [1 1 nlon] );
r_geoid  = 50 * ones(nlat,nlon);
z_ground = 1.5 + 0.5*randn(nlat,nlon);
cb_lims  = [];
%
figure(4)
clf
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 1 );
axes_frame( gca, 'off' )
axis tight
%hl = legend( h, hs );
%disp('Adjust the legend');
%pause;



if yes_or_no('Print figures')
  print atm_dim_3d.eps -depsc
  ! epstopdf atm_dim_3d.eps
  ! rm  atm_dim_3d.eps
  figure(3)
  print atm_dim_3dcross.eps -depsc
  ! epstopdf atm_dim_3dcross.eps
  ! rm  atm_dim_3dcross.eps
  figure(2)
  print atm_dim_2d.eps -depsc
  ! epstopdf atm_dim_2d.eps
  ! rm  atm_dim_2d.eps
  figure(1)
  print atm_dim_1d.eps -depsc
  ! epstopdf atm_dim_1d.eps
  ! rm  atm_dim_1d.eps
end
