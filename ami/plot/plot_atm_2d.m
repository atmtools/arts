%-----------------------------------------------------------------------------
% NAME:     plot_atm_2d
%
%           Plots the grid mesh, the ground and cloud box of a 2D atmosphere.
%
%           The present figure is cleared and the figure is left with hold on.
%
%           The length unit of the figure is set to kilometers.
%
%           The variable dimensions as in ARTS. For example, Z_FIELD is a
%           Tensor3.
%
% FORMAT:   [ h, ltext ] = plot_atm_2d( alpha_grid, z_field, r_geoid, z_ground,
%                                   cloudbox_on, cloudbox_limits, [do_geoid ] )
%
% OUT:      h             Handle to each used linestyle. 
%           ltext         Text suitable for a legend. A legend can be added
% IN:       alpha_grid    Latitude grid.
%           z_field       Altitudes of the pressure surfaces at the crossings
%                         of the pressure and latitude grids.
%           r_geoid       Geoid radius.
%           z_ground      Ground altitude.
%           cloudbox_on   Flag for activating the cloud box.
%           cloudbox_limits   Limits for the cloud box.
% OPTIONAL: do_geoid      Flag to also plot the geoid. Default is 0. 
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function [h,ltext] = plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,do_geoid)


% Number of points per degree for plotting of a pressure surface
%
npd = 0.25;


% Length unit is set to km
lscale = 1e3;
lunit  = 'km';


% Sizes
%
n_alpha = length(alpha_grid);
n_z     = size(z_field,3);


% Clear figure
clf

ltext{1} = 'pressure surfaces';
ltext{2} = 'ground';


for i = 1:n_alpha

  % Plot vertical grid face
  [x1,y1] = cyl2cart( r_geoid(1,i)+z_field(1,i,1), alpha_grid(i) );
  [x2,y2] = cyl2cart( r_geoid(1,i)+z_field(1,i,n_z), alpha_grid(i) );
  h(1) = plot( [x1,x2]/lscale, [y1,y2]/lscale, 'k:' );
  hold on

  if i > 1

    np = max([ 2, ceil(npd*abs(alpha_grid(i)-alpha_grid(i-1))) ]);

    % Plot horisontal grid faces
    for j = 1:n_z
      plot_psurface( alpha_grid(i-1), alpha_grid(i), ...
               (r_geoid(1,i-1)+z_field(1,i-1,j))/lscale, ...
                            (r_geoid(1,i)+z_field(1,i,j))/lscale, np, 'k:' );
    end

    % The ground
    h(2) = plot_psurface( alpha_grid(i-1), alpha_grid(i), ...
               (r_geoid(1,i-1)+z_ground(1,i-1))/lscale, ...
                         (r_geoid(1,i)+z_ground(1,i))/lscale, np, 'b-', 1.5 );

    % The geoid
    if exist('do_geoid','var') & do_geoid
      h(3) = plot_psurface( alpha_grid(i-1), alpha_grid(i), ...
                      r_geoid(1,i-1)/lscale, r_geoid(1,i)/lscale, np, 'g-.' );
      ltext{3} = 'geoid';
    end
  end
end


% Plot cloud box
if exist('cloudbox_on','var') & cloudbox_on
  i = cloudbox_limits(2);
  [x1,y1] = cyl2cart( r_geoid(1,i)+z_field(1,i,1), alpha_grid(i) );
  [x2,y2] = cyl2cart( r_geoid(1,i)+z_field(1,i,cloudbox_limits(1)), ...
                                                               alpha_grid(i) );
  plot( [x1,x2]/lscale, [y1,y2]/lscale, 'm-' );
  i = cloudbox_limits(3);
  [x1,y1] = cyl2cart( r_geoid(1,i)+z_field(1,i,1), alpha_grid(i) );
  [x2,y2] = cyl2cart( r_geoid(1,i)+z_field(1,i,cloudbox_limits(1)), ...
                                                               alpha_grid(i) );
  h2 = plot( [x1,x2]/lscale, [y1,y2]/lscale, 'm-' );
  h = [ h, h2 ];
  ltext{length(ltext)+1} = 'cloud box';
  for i = cloudbox_limits(2):(cloudbox_limits(3)-1)
    np = max([ 2, ceil(npd*abs(alpha_grid(i)-alpha_grid(i+1))) ]);
    plot_psurface( alpha_grid(i), alpha_grid(i+1), ...
               (r_geoid(1,i)+z_field(1,i,cloudbox_limits(1)))/lscale, ...
        (r_geoid(1,i+1)+z_field(1,i+1,cloudbox_limits(1)))/lscale, np, 'm-' );
  end
  
end


xlabel(['x [',lunit,']'])
ylabel(['y [',lunit,']'])


