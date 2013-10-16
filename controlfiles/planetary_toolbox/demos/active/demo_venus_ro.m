% DEMO_VENUS_RO   Simulates a Venus radio experiment
%
%    The demo covers a satellite in a orbit around Venus. 
%
%    Bending angles and geometrical/length variables are calculated. See
%    demo_earth_ro for a calculation involving attenuation quantities.
%
% FORMAT   demo_venus_ro

% Patrick Eriksson 2013-10-11

function demo_venus_ro

%-------
% O part
%-------

% Satellite orbit altitude. Note that the satellite is assumed to be in a
% circular orbit around the planet.
%
O.gps_altitude = 600e3;
O.gps_movement = 'disappearing';

% Here we let Earth be the "LEO"
%
% The closest distance Earth - Venus is about 38e12, but using this value
% was found to give some occasionaly bad results.
% The quantities affcted by the distance selected are absolut path lengths and
% free space loss. But excess range is OK if just a high value is selected.
% For a distance of >= 38e12, the free space loss is basically constant
% during the occulation, as the relative change in distance is very small.
%  
%O.leo_altitude = 38e12;
O.leo_altitude = 1e10;    % 1e10 seems to give stable results
O.leo_movement = 'none';

% These are the two frequencies of Venus-Express. Select one.
%
%O.frequency    = 2.3e9;     
O.frequency    = 8.4e9;   % This one gives much lower impact of free
                          % electrons, but higher gas absorption


% For higher accuracy, decrease these values, particularly lraytrace.
% lmax mainly affects the accuracy for the tangent point determination
%
O.lmax         = 10e3;
O.lraytrace    = 100;

% Surface altitude 
%
% Here of not important as the rays don't get below about 45 km
%
O.z_surface    = 1e3;

% These settings determine end point of occultation, how dense calculations
% that are performed etc. See *arts_radioocc_1D_power* for details.
%
% These actual settings gives a rough overview
%
O.slta_max     = 200e3;
O.slta_min     = -200e3;
O.slta_n       = 50;
O.z_impact4t0  = 100e3;  % Sets reference point for time
O.f_sampling   = 4;

  
%-------
% A part
%-------

% Don't change these (if you don't switch planet)!
%
A.planet       = 'venus';
A.atmfunc      = @qarts_add_venus_planettbox;

A.atmo         = 3;               % Atmospheric scenario
%
A.basespecies  = [1,3,6,14];      % This is CO2, H2SO4, N2 and SO2
A.h2ospecies   = 1;               % Level of water vapour
A.hdospecies   = 3;               % Level of HDO
A.Necase       = 2;               % Free electron case. Note that 5-6 needs
%                                 % atm=0-2, while 0-4 to needs atmo=3-4
A.interp_order = 1;               % Linear interpolation of fields (higher
%                                   values risky
A.pmin         = 1e-6;            % Min pressure to consider. This value
                                  % crops around 200 km


% Here we need to a special solution, to avoid the need to manually generate
% "linefiles" 
%
Q = qarts;
%
% Increasing df gives better accuracy (more transitions are included in the
% calculation, but gives slower calculations). 10 GHz should be sufficient if
% not max accuracy is demanded. However, the absorption is totally dominated by
% continuum terms, df can in principle be zero.
df = 10e9;
%
Q.ABS_WSMS{end+1} = sprintf( ['abs_linesReadFromSplitArtscat(abs_lines,',...
    'abs_species,"spectroscopy/Perrin/",%.2e,%.2e)'], ...
                                 max([O.frequency-df,0]), O.frequency+df );
Q.ABS_WSMS{end+1} = 'abs_lines_per_speciesCreateFromLines';
Q.ABS_WSMS{end+1} = 'abs_lines_per_speciesAddMirrorLines';


%- Perform calculation
%
[R,T] = arts_radioocc_1D_power( Q, O, A );



%- Plot results

tstring = sprintf( 'Venus: Atmosphere %d, Necase %d, %.2f GHz', ...
                                           A.atmo, A.Necase, O.frequency/1e9 );
fs = 14;

figure(1)
clf  
plot( R.bangle, R.z_impact/1e3 )
%
grid
xlabel( 'Bending angle [deg]', 'FontSize', fs );
ylabel( 'Impact height [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )
axis([-1e-3 1e-3 80 200])


figure(2)
clf  
clonefig(1,2)
axis([0 5.5 45 95])

figure(3)
clf  
plot( R.l_optpath-R.l_geometric, R.slta/1e3 )
%
grid
xlabel( 'Excess range [m]', 'FontSize', fs );
ylabel( 'Straight-line tangent altitude [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

% Don't show total or free space as we are not using the exact length. But
% report what the free spae loss should be.
%
lmin = 38e12;
db0 = 10*log10( 4 * pi * lmin^2);  
%
figure(4)
clf  
plot( -10*log10(R.tr_atmos), R.z_tan/1e3, 'r-', ...
      -10*log10(R.tr_defoc), R.z_tan/1e3, 'b-', 'LineWidth', 2 );
%
grid
xlabel( 'Attenuation [dB]', 'FontSize', fs );
ylabel( 'Tangent height [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )
legend( ' Absorption', ' Defocusing' );
axis([0 22 35 105])
ax = axis;
text( ax(2)/4, ax(3)+0.7*diff(ax(3:4)), sprintf(['Free space loss is at ',...
'least %.0f db, and is\nbasically constant during the occultation'], db0 ));

figure(5)
clf  
plot( T.t, T.bangle )
%
grid
xlabel( 'Relative time [s]', 'FontSize', fs );
ylabel( 'Bending angle [deg]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

