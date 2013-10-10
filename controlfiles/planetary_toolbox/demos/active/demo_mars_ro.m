function demo_mars_ro

% Here we let Eart be the "GPS"
O.leo_altitude = 600e3;
O.gps_altitude = 60e9;  % Roughly closest distance Earth - MARS   
O.gps_movement = 'disappearing';
O.frequency    = 6e9;
O.lmax         = 2e3;
O.lraytrace    = 200;
O.z_surface    = 10;

O.z_impact_max = 50e3;
O.z_impact_dz  = 500;
O.z_impact4t0  = O.z_impact_max;
O.f_sampling   = 4;


A.planet       = 'mars';
A.atmfunc      = @qarts_add_mars_planettbox;

A.Ls           = 2;
A.daytime      = 1;
A.dust         = 0;
A.solar        = 1;
%
A.interp_order = 1;
%
A.pmin         = 1e-3;
%
A.basespecies  = [ 0, 1, 6, 9 ];
A.ch4species   = 1;
A.h2ospecies   = 1;
A.Necase       = 4;


[R,T] = arts_radioocc_1D( [], O, A );

plot( R.bangle, R.z_impact/1e3 )
%
xlabel( 'Bending angle [deg]' );
ylabel( 'Imapct height [km]' );