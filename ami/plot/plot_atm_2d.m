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
% FORMAT:   plot_atm_2d( alpha_grid, z_field, r_geoid, z_ground,
%                                               cloudbox_on, cloudbox_limits )
%
% IN:       alpha_grid    Latitude grid.
%           z_field       Altitudes of the pressure surfaces at the crossings
%                         of the pressure and latitude grids.
%           r_geoid       Geoid radius.
%           z_ground      Ground altitude.
%           cloudbox_on   Flag for activating the cloud box.
%           cloudbox_limits   Limits for the cloud box.
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function plot_atm_2d(alpha_grid,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits)


% Number of points per degree for plotting of a pressure surface
%
npd = 5;


% Length unit is set to km
lscale = 1e3;
lunit  = 'km';


% Sizes
%
n_alpha = length(alpha_grid);
n_z     = size(z_field,3);


% Clear figure
clf


for i = 1:n_alpha

  % Plot vertical grid face
  [x1,y1] = cyl2cart( r_geoid(1,i)+z_field(1,i,1), alpha_grid(i) );
  [x2,y2] = cyl2cart( r_geoid(1,i)+z_field(1,i,n_z), alpha_grid(i) );
  plot( [x1,x2]/lscale, [y1,y2]/lscale, 'k:' );
  hold on

  if i > 1

    % Plot horisontal grid faces
    for j = 1:n_z
      plot_psurface( alpha_grid(i-1), alpha_grid(i), ...
               (r_geoid(1,i-1)+z_field(1,i-1,j))/lscale, ...
                            (r_geoid(1,i)+z_field(1,i,j))/lscale, npd, 'k:' );
    end

    % The ground
    plot_psurface( alpha_grid(i-1), alpha_grid(i), ...
               (r_geoid(1,i-1)+z_ground(1,i-1))/lscale, ...
                              (r_geoid(1,i)+z_ground(1,i))/lscale, npd, 'b-' );
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
  plot( [x1,x2]/lscale, [y1,y2]/lscale, 'm-' );
  for i = cloudbox_limits(2):(cloudbox_limits(3)-1)
    plot_psurface( alpha_grid(i), alpha_grid(i+1), ...
               (r_geoid(1,i)+z_field(1,i,cloudbox_limits(1)))/lscale, ...
        (r_geoid(1,i+1)+z_field(1,i+1,cloudbox_limits(1)))/lscale, npd, 'm-' );
  end
  
end


xlabel(['x [',lunit,']'])
ylabel(['y [',lunit,']'])


