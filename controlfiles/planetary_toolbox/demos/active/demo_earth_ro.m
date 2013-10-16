% DEMO_EARTH_RO  Earth radio occultation demo (troposphere + attenuation)
%
%    A demonstration of simulating a GPS occultation, using the function
%    *arts_radioocc_1D* combined with *qarts_add_fascode*. See these
%    functions for details around the settings inside the script.
%
%    The function focuses on the part passing the lowest approx 20 km.
%    Both bending angles, attenuation quantities and Faraday rotation
%    are determined.
%
%    These calculations are slow as exact transmitter - reciver paths are
%    determined. If you just need bending angles, use *arts_radioocc_1D* for
%    (much) faster calculations!
%
% FORMAT   demo_earth_ro

% 2013-10-16   Created by Patrick Eriksson.

function demo_earth_ro
  
% When we select tropical, we obtain free electrons for a low latitude 
% (10N to be precise). This latitude is shifted if you switch Fascode
% atmosphere. The longitude is selected by Q.LON_TRUE.
  
% Define observation
%
O.tra_altitude = 20200e3;
O.tra_movement = 'disappearing';
O.rec_altitude = 820e3;
O.tra_movement = 'disappearing';
O.frequency    = 1575.42e6;
O.lmax         = 2e3;
O.lraytrace    = 100;
O.z_surface    = 10;

O.z_impact_max = 20e3;
O.z_impact_dz  = 100;
O.z_impact4t0  = O.z_impact_max;
O.f_sampling   = 4;

O.do_atten     = true;
O.do_faraday   = true;


% Init Q and set a longitude
%
Q              = qarts;
Q.LON_TRUE     = 5;


% Define the atmosphere
%
A.planet       = 'earth';
A.atmfunc      = @qarts_add_fascode;
A.fascode_atm  = 'tropical';
%
% The species that contribute to refraction and absorption
% (Here we use fiull absorption models for the gases)
A.ABS_SPECIES(1).TAG{1} = 'N2-SelfContStandardType';
A.ABS_SPECIES(2).TAG{1} = 'O2-PWR98';
A.ABS_SPECIES(3).TAG{1} = 'H2O-PWR98';
A.ABS_SPECIES(4).TAG{1} = 'free_electrons';
%
A.igrf_year    = 2010;  % Year for magnetic field
%
A.ionosphere   = true;
A.iri_sun      = 'max'; 
A.iri_season   = 'winter';
A.iri_utc      = 0;
%
% IRI gives a high Ne all the way up to 2000 km, close to constant from 1000
% km. Seems too high and we use the option to crop!
A.p_min        = z2p_simple( 1200e3 );


[R,T] = arts_radioocc_1D( Q, O, A );


%- Plot results

tstring = sprintf( ...
    'Earth: Season %s, Solar activity %s, %.2f GHz', ...
               A.iri_season, A.iri_sun, O.frequency/1e9 );
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
db0 = -10*log10(R.tr_space(1));
%
figure(3)
clf  
plot( -10*log10(R.tr_space)-db0, R.z_impact/1e3, 'm-', ...
      -10*log10(R.tr_atmos), R.z_impact/1e3, 'r-', ...
      -10*log10(R.tr_defoc), R.z_impact/1e3, 'b-', ...
      -10*log10(R.tr_total)-db0, R.z_impact/1e3, 'k-', 'LineWidth', 2 );
%
grid
xlabel( 'Attenuation [dB]', 'FontSize', fs );
ylabel( 'Impact height [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )
legend( sprintf(' Free space - %.1f dB',db0), ' Absorption', ...
        ' Defocusing', sprintf(' Total - %.1f dB',db0) );
ax = axis;
text( ax(2)/2.5, mean(ax(3:4)), ...
       sprintf('%.1f dB is the free space loss at start',db0) );

figure(4)
clf  
plot( R.faraday, R.z_impact/1e3 )
%
grid
xlabel( 'Faraday rotation [deg]', 'FontSize', fs );
ylabel( 'Impact height [km]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

figure(5)
clf  
plot( T.t, T.bangle )
%
grid
xlabel( 'Relative time [s]', 'FontSize', fs );
ylabel( 'Bending angle [deg]', 'FontSize', fs );
title( tstring, 'FontSize', fs+2 )

