%-----------------------------------------------------------------------------
% NAME:     plot_atm_1d
%
%           Plots the vertical grid, the ground and cloud box of a 1D 
%           atmosphere.
%
%           The present figure is cleared and the figure is left with hold on.
%
%           The length unit of the figure is set to kilometers.
%
%           The variable dimensions as in ARTS. For example, Z_FIELD is a
%           Tensor3.
%
% FORMAT:   [ h, ltext ] = plot_atm_1d( dalpha, z_field, r_geoid, z_ground,
%                                   cloudbox_on, cloudbox_limits, [do_geoid ] )
% 
% OUT:      h             Handle to each used linestyle. 
%           ltext         Text suitable for a legend. A legend can be added
%                         to the figure by: legend(h,ltext);
% IN:       dalpha        The angular size of the plotted atmosphere. That is,
%                         the value 180 gives a half cylinder.
%           z_field       Altitudes of the pressure surfaces at the crossings
%                         of the pressure and latitude grids.
%           r_geoid       Geoid radius.
%           z_ground      Ground altitude.
%           cloudbox_on   Flag for activating the cloud box.
%           cloudbox_limits   Limits for the cloud box.
% OPTIONAL: do_geoid      Flag to also plot the geoid. Default is 0. 
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function [h,ltext] = plot_atm_1d(dalpha,z_field,r_geoid,z_ground,cloudbox_on,cloudbox_limits,do_geoid)


% Number of points per degree for plotting of a pressure surface
%
npd = 0.25;


% Length unit is set to km
lscale = 1e3;
lunit  = 'km';


% Sizes
%
n_z     = size(z_field,3);


% Clear figure
clf
hold on

np = max([ 2, ceil(npd*abs(dalpha)) ]);


% Plot horisontal grid faces
for j = 1:n_z
  h =  plot_psurface( -dalpha/2, dalpha/2, ....
                            (r_geoid(1,1)+z_field(1,1,j))/lscale, ...
                            (r_geoid(1,1)+z_field(1,1,j))/lscale, np, 'k:' );
end
ltext{1} = 'pressure surfaces';


% The ground
h(2) = plot_psurface( -dalpha/2, dalpha/2, ...
                         (r_geoid(1,1)+z_ground(1,1))/lscale, ...
                         (r_geoid(1,1)+z_ground(1,1))/lscale, np, 'b-', 2 );
ltext{2} = 'ground';


% The geoid
if exist('do_geoid','var') & do_geoid
  h(3) = plot_psurface( -dalpha/2, dalpha/2, ...
                        r_geoid(1,1)/lscale, r_geoid(1,1)/lscale, np, 'g-.' );
  ltext{3} = 'geoid';
end


% Plot cloud box
if exist('cloudbox_on','var') & cloudbox_on
  h2 = plot_psurface( -dalpha/2, dalpha/2, ...
            (r_geoid(1,1)+z_field(1,1,cloudbox_limits(1)))/lscale, ...
            (r_geoid(1,1)+z_field(1,1,cloudbox_limits(1)))/lscale, np, 'm-' );
  h = [ h, h2 ];
  ltext{length(ltext)+1} = 'cloud box';
end


xlabel(['x [',lunit,']'])
ylabel(['y [',lunit,']'])


