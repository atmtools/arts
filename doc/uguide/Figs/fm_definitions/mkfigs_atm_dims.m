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
cloudbox_limits  = [0 5];


figure(1)
[h,ltext]=plot_atm_1d(dalpha,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,1);
for i = 1:size(z_field,3)
  hp = plot( 0, (r_geoid(1,i)+z_field(1,1,i))/1e3, 'k.' );
end
nh = length(h)+1;
h(nh) = hp;
ltext{nh} = 'atmospheric field';
hf = gca;
hl = legend2(8,h,ltext); 
legend boxoff 
axis equal 
hide_axes(hf);
disp('Adjust the legend');
pause;
axis tight


figure(2)
%
cloudbox_limits  = [3 5 4 10];
%
[h,ltext]=plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,1);
for i = 1:size(z_field,3)
  for j = 1:nalpha
    [x1,y1] = cyl2cart( r_geoid(1,j)+z_field(1,j,i), alpha_grid(j) );
    hp = plot( x1/1e3, y1/1e3, 'k.' );
  end
end
nh = length(h)+1;
h(nh) = hp;
ltext{nh} = 'atmospheric field';
hf = gca;
hl = legend2(8,h,ltext); 
legend boxoff 
axis equal 
hide_axes(hf);
disp('Adjust the legend');
pause;
axis tight


if answer_is_yes('Print figures')
  print atm_dim_2d.eps -depsc
  ! epstopdf atm_dim_2d.eps
  figure(1)
  print atm_dim_1d.eps -depsc
  ! epstopdf atm_dim_1d.eps
end
