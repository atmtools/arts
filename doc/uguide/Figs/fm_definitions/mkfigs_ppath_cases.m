function mkfigs_ppath_cases

figure(1)

%--- Define a 2D atmosphere

dalpha      = 90; 
nalpha      = 15;
alpha_grid  = linspace( -dalpha/2, dalpha/2, nalpha );
beta_grid   = [];
r_geoid     = 50e3 + zeros( 1, nalpha );
z_ground    = 2e3*ones( 1, nalpha );
z_ground(3) = 1500;
z_ground(5) = 2500;
z_ground(10)= 1e3;
z_ground(11)= 1e3;

z_field(1,:,:) = repmat( 1e3:2e3:18e3, nalpha, 1 );
for i = 1:nalpha
  z_field(1,i,:) = z_field(1,i,:) - (i-1)*1e3/nalpha;
end
p_grid = 1:size(z_field,3); 

cloudbox_on      = 1;
cloudbox_limits  = [1 5 4 10];

plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,1);
%
r_s     = 56e3;
alpha_s = -dalpha/2 + 1;
beta_s  = [];
omega_s = 0;
%
psi_s   = 30;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,0,cloudbox_on,cloudbox_limits,100e3,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
psi_s   = 95;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,0,cloudbox_on,cloudbox_limits,100e3,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
psi_s   = 120;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,0,cloudbox_on,cloudbox_limits,100e3,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
r_s     = 70e3;
alpha_s = dalpha/2 - 1;
%
psi_s   = -119;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,1,cloudbox_on,cloudbox_limits,100e3,r_s,alpha_s,beta_s,psi_s,omega_s);
[h,ltext] = plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
psi_s   = -130;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,1,cloudbox_on,cloudbox_limits,100e3,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
psi_s   = -170;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,1,cloudbox_on,cloudbox_limits,100e3,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
hf = gca;
hl = legend2(8,h,ltext); 
set(hl,'Visible','Off'); 
axis equal 
hide_axes(hf);
disp('Adjust the legend');
pause;
axis tight



figure(2)

%--- Define a new 2D atmosphere

dalpha      = 4; 
nalpha      = 5;
alpha_grid  = linspace( -dalpha/2, dalpha/2, nalpha );
beta_grid   = [];
r_geoid     = 50e3 + zeros( 1, nalpha );
z_ground    = 0e3*ones( 1, nalpha );
z_ground(4:5) = 500;

clear z_field
z_field(1,:,:) = repmat( 0e3:1e3:1e3, nalpha, 1 );
p_grid = 1:size(z_field,3); 

cloudbox_on      = 0;
cloudbox_limits  = [1 5 4 10];

plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,1);
%
r_s     = 51e3;
alpha_s = 0.25;
beta_s  = [];
omega_s = 0;
%
psi_s   = 138;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,0,cloudbox_on,cloudbox_limits,500,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
psi_s   = -138;
ppath = ppath_2d_3d(2,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,0,0,cloudbox_on,cloudbox_limits,500,r_s,alpha_s,beta_s,psi_s,omega_s);
plot_ppath( r_s, alpha_s, ppath, alpha_grid, r_geoid );
%
hf = gca;
axis equal 
hide_axes(hf);
axis tight


if answer_is_yes('Print figures')
  print ppath_cases2.eps -depsc
  ! epstopdf ppath_cases2.eps
  figure(1)
  print ppath_cases1.eps -depsc
  ! epstopdf ppath_cases1.eps
end
