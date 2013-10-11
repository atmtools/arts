function demo_mars_ro

% Here we let Earth be the "GPS"
O.leo_altitude = 600e3;
O.gps_altitude = 60e12;     % Roughly closest distance Earth - Mars 
O.gps_movement = 'none';
O.frequency    = 2.3e9;     % These are the two frequencies of Mars-Express
%O.frequency    = 8.4e9;    % Select one.   
O.lmax         = 2e3;
O.lraytrace    = 200;
O.z_surface    = 10;

O.z_impact_max = 250e3;
O.z_impact_dz  = 3e3;
O.z_impact4t0  = O.z_impact_max;
O.f_sampling   = 4;


A.planet       = 'mars';
A.atmfunc      = @qarts_add_mars_planettbox;

A.Ls           = 2;               % Season
A.daytime      = 0;               % Day / night
A.dust         = 0;               % Dust level
A.solar        = 2;               % Solar activity
%
A.basespecies  = [ 0, 1, 6, 9 ];  % This is CO, CO2, N2 and N2O
A.ch4species   = 1;               % Standard CH4
A.h2ospecies   = 1;               % Include all water as one species
A.Necase       = 0;               % Free electron case. Note that 4 needs
%                                 % daytime=1, while 0-3 to needs daytime=0
A.interp_order = 1;               % Linear interpolation of fields (higher
%                                   values risky)
A.pmin         = 1e-3;            % Min pressure to consider. This value
                                  % crops a bit above 100 km

O.gps_altitude = 600e3;    
A.pmin         = 1e-99;


%- Perform calculation
%
[R,T] = arts_radioocc_1D( [], O, A );


%-Plot result
%
plot( R.bangle, R.z_impact/1e3 )
%
xlabel( 'Bending angle [deg]' );
ylabel( 'Imapct height [km]' );