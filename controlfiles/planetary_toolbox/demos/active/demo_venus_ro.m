function demo_venus_ro

% Here we let Earth be the "GPS"
O.leo_altitude = 600e3;
O.gps_altitude = 38e12;  % Roughly closest distance Earth - Venus   
O.gps_movement = 'none';
O.frequency    = 6e9;
O.lmax         = 1e3;
O.lraytrace    = 200;
O.z_surface    = 10;

O.slta_max     = 200e3;
O.slta_min     = 0e3;
O.slta_n       = 20;
O.z_impact4t0  = 100e3;
O.f_sampling   = 4;


A.planet       = 'venus';
A.atmfunc      = @qarts_add_venus_planettbox;

A.atmo         = 3;               % Atmospheric scenario
%
A.basespecies  = [ 1, 5 ];        % This is CO and N2
A.h2ospecies   = 1;               % Level of water vapour
A.hdospecies   = 3;               % Level of HDO
A.Necase       = 2;               % Free electron case. Note that 5-6 needs
%                                 % atm=0-2, while 0-4 to needs atmo=3-4
A.interp_order = 1;               % Linear interpolation of fields (higher
%                                   values risky
A.pmin         = 1e-6;            % Min pressure to consider. This value
                                  % crops around 200 km

O.gps_altitude = 600e3;    
A.pmin         = 1e-99;


%- Perform calculation
%
[R,T] = arts_radioocc_1D_power( [], O, A );


%-Plot result
%
plot( R.bangle, R.z_impact/1e3 )
%
xlabel( 'Bending angle [deg]' );
ylabel( 'Imapct height [km]' );
