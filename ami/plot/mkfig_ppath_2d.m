%-----------------------------------------------------------------------------
% NAME:     mkfig_ppath_2d
%
%           This is just a temporary function that will be removed from AMI
%           at a later stage.
%
% FORMAT:   -
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


dim        = 2;

p_grid     = 0:12;
alpha_grid = -10:2:10;
beta_grid  = [];

z_field        = zeros( 1, length(alpha_grid), length(p_grid) );
z_field(1,:,:) = repmat( 0:1e3:1e3*(length(p_grid)-1), length(alpha_grid), 1 );
for i = 1:length(alpha_grid)
  z_field(1,i,:) = z_field(1,i,:) + (i-1)*50;
end

global EARTH_RADIUS

r0 = 50e3;

r_geoid    = r0 + zeros( 1, length(alpha_grid) );
%r_geoid    = linspace( r0, r0+2e3, length(alpha_grid) );
z_ground   = 500*ones( 1, length(alpha_grid) );
%z_ground   = linspace( 0, 2e3, length(alpha_grid) );
 
refr_on          = 0;
blackbody_ground = 0;
cloudbox_on      = 1;
cloudbox_limits  = [8 4 10];
ppath_lmax       = 200e3;

r_s       = r0 + 3.25e3;
alpha_s   = -4;
beta_s    = [];
psi_s     = -170;
omega_psi = [];


ppath = ppath_2d_3d(dim,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,refr_on,blackbody_ground,cloudbox_on,cloudbox_limits,ppath_lmax,r_s,alpha_s,beta_s,psi_s,omega_psi);

close all

plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits);

plot_ppath(r_s,alpha_s,ppath,alpha_grid, r_geoid);

zoom on
