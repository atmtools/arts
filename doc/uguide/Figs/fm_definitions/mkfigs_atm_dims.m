function mkfigs_atm_dims


%--- 1D
%
dim      = 1;
lat_grid = [45 135];
lon_grid = [];
z_field  = (1:18)';
r_geoid  = 50;
z_ground = 1.5;
cb_lims  = [6 9];
%
figure(1)
clf
fit_to_paper( 'a5l' );
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 1 );
axes_frame( gca, 'off' )
axis tight
hl = legend( h, hs );
disp('Adjust the legend');
%pause;



%--- 2D
%
nlat     = 16;
%
dim      = 2;
lat_grid = linspace(45,135,nlat);
lon_grid = [];
z_field  = (1:18)'*ones(1,nlat);% - ones(18,1)*linspace(0,0.8,nlat);
r_geoid  = 50 * ones(nlat,1);
z_ground = 1.5 + 2*randn(nlat,1)
cb_lims  = [6 9 8 13];
%
figure(2)
clf
fit_to_paper( 'a5l' );
[h,hs] = arts_plot_atmgrids( dim, lat_grid, lon_grid, z_field,...
                                               r_geoid, z_ground, cb_lims, 1 );
axes_frame( gca, 'off' )
axis tight
hl = legend( h, hs );
disp('Adjust the legend');
%pause;

return

if answer_is_yes('Print figures')
  print atm_dim_2d.eps -depsc
  ! epstopdf atm_dim_2d.eps
  figure(1)
  print atm_dim_1d.eps -depsc
  ! epstopdf atm_dim_1d.eps
end
