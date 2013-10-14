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
O.leo_altitude = 600e3;

% Here we let Earth be the "GPS", and thus calculates backwards but this has
% no impact on the result.
%
% The closest distance Earth - Venus is about 38e12, but using this value
% was found to give numerical problems, and a lower (but still high)
% value is used. 
% The quantities affcted by the distance selected are absolut path lengths and
% free space loss. But excess range and the variation of free space loss during
% the occultation are OK if a just a high value is selected.
%  
O.gps_altitude = 1e8;  
O.gps_movement = 'none';

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
O.slta_max     = 100e3;
O.slta_min     = -20e3;
O.slta_n       = 21;
O.z_impact4t0  = O.slta_max;  % Sets reference point for time
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

tstring = sprintf( 'Venus: Atmosphere %d,Necase %d, %.2f GHz', ...
                                           A.atmo, A.Necase, O.frequency/1e9 );
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

% Attenuation dominated by free space loss and remove this start value of this
% dominatating part
%
db0 = 10*log10( 4 * pi * O.gps_altitude^2);  %-10*log10(R.tr_space(1));
%
figure(3)
clf  
plot( -10*log10(R.tr_space)-db0, R.z_tan/1e3, 'm-', ...
      -10*log10(R.tr_atmos), R.z_tan/1e3, 'r-', ...
      -10*log10(R.tr_defoc), R.z_tan/1e3, 'b-', ...
      -10*log10(R.tr_total)-db0, R.z_tan/1e3, 'k-', 'LineWidth', 2 );
%
grid
xlabel( 'Attenuation [dB]', 'FontSize', fs );
ylabel( 'Tangent height [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )
legend( sprintf(' Free space - %.2f dB',db0), ' Absorption', ...
        ' Defocusing', sprintf(' Total - %.2f dB',db0) );
ax = axis;
text( ax(2)/4, ax(3)+0.7*diff(ax(3:4)), sprintf(['%.2f dB is the free ',...
'space loss for\n the used planet distance of %.1e m'], db0, O.gps_altitude ));

figure(4)
clf  
plot( T.t, T.bangle )
%
grid
xlabel( 'Relative time [s]', 'FontSize', fs );
ylabel( 'Bending angle [deg]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

