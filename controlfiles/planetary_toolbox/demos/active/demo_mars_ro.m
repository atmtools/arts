% DEMO_MARS_RO   Simulates a Mars radio experiment
%
%    The demo covers a satellite in a orbit around Mars. 
%
%    Bending angles and geometrical/length variables are calculated. See
%    demo_earth_ro for a calculation involving attenuation quantities.
%
% FORMAT   [R,T,O,A] = demo_mars_ro
%
% See *arts_radioocc_1D* for definition of the output arguemnts.

% Patrick Eriksson 2013-10-16

function [R,T,O,A] = demo_mars_ro

%-------
% O part
%-------

% Satellite orbit altitude. Note that the satellite is assumed to be in a
% circular orbit around the planet.
%
O.tra_altitude = 1000e3;
O.tra_movement = 'disappearing';


% The closest distance Earth - Mars is about 60e12.
%
% Using the actual distance seems to work here. However, please, double-check
% the results by testing that the results are not changed by using e.g. 1e10.
%
O.rec_altitude = 60e12;
O.rec_movement = 'none';

% These are the two frequencies of Mars-Express. Select one.
%
O.frequency    = 2.3e9;     
%O.frequency    = 8.4e9; % This one gives much lower impact of free electrons

% For higher accuracy, decrease these values, particularly lraytrace.
% lmax mainly affects the accuracy for the tangent point determination
%
O.lmax         = 10e3;
O.lraytrace    = 1e3;

% Surface altitude (the input data start around 5m, and 0m does not work)
%
O.z_surface    = 5;

% These settings determine end point of occultation, how dense calculations
% that are performed etc. See *arts_radioocc_1D* for details.
%
% These actual settings gives a rough overview
%
O.z_impact_max = 250e3;
O.z_impact_dz  = 2.5e3;
O.z_impact4t0  = O.z_impact_max;  % Sets reference point for time
O.f_sampling   = 4;

  
%-------
% A part
%-------

% Don't change these (if you don't switch planet)!
%
A.planet       = 'mars';
A.atmfunc      = @qarts_add_mars_planettbox;


% See DemoMarsAtmo1D.arts for a specification of these settings. These
% fields of A are named exactly as DemoMarsAtmo1D.arts.
%
A.Ls           = 2;               % Season
A.daytime      = 0;               % Day / night
A.dust         = 0;               % Dust level
A.solar        = 2;               % Solar activity
%
A.basespecies  = [ 1, 6 ];        % This is CO2 and N2
A.ch4species   = 1;               % Standard CH4
A.h2ospecies   = 1;               % Include all water as one species
A.Necase       = 1;               % Free electron case. Note that 4 needs
%                                 % daytime=1, while 0-3 to needs daytime=0
A.interp_order = 1;               % Linear interpolation of fields (higher
%                                   values risky)
A.pmin         = 1e-99;           % Min pressure to consider. This value
                                  % gives no cropping. An example:
%---                              % 1e-3 crops around 100 km

% Regarding species selected above, remember that we here only calculate
% refractive index, that is done by *refr_index_airMWgeneral*. No need to
% include species not considered by this WSM (on the other hand, but no big
% overhead to load more species, such as CH4).


%- Perform calculation
%
[R,T] = arts_radioocc_1D( [], O, A );


%- Plot results

tstring = sprintf( ...
    'Mars: Season %d, Daytime %d, Solar %d, Necase %d, %.2f GHz', ...
               A.Ls, A.daytime, A.solar, A.Necase, O.frequency/1e9 );
fs = 12;

figure(1)
clf  
plot( R.bangle, R.z_impact/1e3 )
%
grid
xlabel( 'Bending angle [deg]', 'FontSize', fs );
ylabel( 'Impact height [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

figure(2)
clf  
plot( R.l_optpath-R.l_geometric, R.slta/1e3 )
%
grid
xlabel( 'Excess range [m]', 'FontSize', fs );
ylabel( 'Straight-line tangent altitude [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

figure(3)
clf  
plot( T.t, T.bangle )
%
grid
xlabel( 'Relative time [s]', 'FontSize', fs );
ylabel( 'Bending angle [deg]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )
