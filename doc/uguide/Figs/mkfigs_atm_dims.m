function mkfigs_atm_dims


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
 
cloudbox_on      = 1;
cloudbox_limits  = [5 4 10];


figure(1)
[h,ltext]=plot_atm_1d(dalpha,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,1);
hf = gca;
hl = legend2(8,h,ltext); 
set(hl,'Visible','Off'); 
axis equal 
hide_axes(hf);
disp('Adjust the legend');
pause;
axis tight


figure(2)
[h,ltext]=plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,1);
hf = gca;
%hl = legend(h,ltext); 
%set(hl,'Visible','Off'); 
axis equal 
hide_axes(hf);
%disp('Adjust the legend');
%pause;
axis tight


if if_yes('Print figures')
  print atm_dim_2d.eps -depsc
  ! ps2pdf atm_dim_2d.eps atm_dim_2d.pdf
  figure(1)
  print atm_dim_1d.eps -depsc
  ! ps2pdf atm_dim_1d.eps atm_dim_1d.pdf
end
