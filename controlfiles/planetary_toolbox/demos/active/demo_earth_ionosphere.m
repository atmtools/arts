% DEMO_EARTH_IONOSPHERE  Earth radio occultation demo (ionospheric focus)
%
%    A demonstration of simulating a GPS occultation, using the function
%    *arts_radioocc_1D* combined with *qarts_add_fascode*. See these
%    functions for details around the settings inside the script.
%
%    Focus is put on high altitudes, with the aim of generating a plot
%    similar to Fig 3 in "Analysis and validation of GPS/MET radio
%    occultation data in the ionosphere", by Schreiner et al., Radio
%    Science, 1999. The figure here is extended down to surface level, and
%    includes results for the case of no ionosphere.
%
%    Note that there is no attempt to recreate the exact profile of free
%    electrons used by Schreiner et al., IRI is instead used. Highly
%    varying results are obtained when changing the IRI settings (such as
%    longitude and UTC).
%
% FORMAT   [R,T,O,A] = demo_earth_ionosphere
%
% See *arts_radioocc_1D* for definition of the output arguemnts.

% 2013-10-16   Created by Patrick Eriksson.

function demo_earth_ionosphere

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

O.z_impact_max = 700e3;
O.z_impact_dz  = 20e3;
O.z_impact4t0  = O.z_impact_max;
O.f_sampling   = 0.5;


% Init Q and set a longitude
%
Q          = qarts;
Q.LON_TRUE = 5;


% Define a planet, without ionosphere, and calculate
%
% Note that no data are generated above top-of-the-atmosphere, 
% here about at 95 km 
%
A.planet       = 'earth';
A.atmfunc      = @qarts_add_fascode;
A.fascode_atm  = 'tropical';
%
% The species that contribute to refraction
A.ABS_SPECIES(1).TAG{1} = 'N2';
A.ABS_SPECIES(2).TAG{1} = 'O2';
A.ABS_SPECIES(3).TAG{1} = 'H2O';
A.ABS_SPECIES(4).TAG{1} = 'free_electrons';
%
R1 = arts_radioocc_1D( Q, O, A );


% Add an ionosphere, for low solar activity and night
%
A.igrf_year    = 2010;  % Year for magnetic field
%
A.ionosphere   = true;
A.iri_sun      = 'min'; 
A.iri_season   = 'winter';
A.iri_utc      = 0;
%
% IRI gives a high Ne all the way up to 2000 km, close to constant from 1000
% km. Seems too high and we use the option to crop!
A.p_min        = z2p_simple( 1200e3 );
%
R2 = arts_radioocc_1D( Q, O, A );


% Switch to high solar activity and daytime
%
A.iri_sun      = 'max'; 
A.iri_season   = 'winter';
A.iri_utc      = 12;
%
R3 = arts_radioocc_1D( Q, O, A );


fac = 1000 * pi/180;


% Set-up figure
%
figure(1),clf
set_figsize( 250, 100 );
h = add_plot_row( 0.8, 0.68*ones(3,1), 0.07, 0.12, 'r', 0.22 );

fs = 12;

% Bending angles
axes( h(1) ); 
plot( [0 0],[0 O.z_impact_max/1e3], 'k:', ...
        fac*R1.bangle, R1.z_impact/1e3, '-', ... 
        fac*R2.bangle, R2.z_impact/1e3, '-.', ... 
        fac*R3.bangle, R3.z_impact/1e3, '--' );
axis( [ -0.5 0.5 0 O.z_impact_max/1e3] );
xlabel( 'Bending angle [mrad]', 'FontSize', fs )
ylabel( 'Impact height [km]', 'FontSize', fs )


% Excess phase
axes( h(2) ); 
plot( [0 0],[0 O.z_impact_max/1e3], 'k:', ...
      R1.l_optpath-R1.l_geometric, R1.z_impact/1e3, '-', ... 
      R2.l_optpath-R2.l_geometric, R2.z_impact/1e3, '-.', ... 
      R3.l_optpath-R3.l_geometric, R3.z_impact/1e3, '--' );
axis( [ -150 0 0 O.z_impact_max/1e3] );
xlabel( 'Excess range [m]', 'FontSize', fs )


% Deviation to SLTA
axes( h(3) ); 
l=plot( [0 0],[0 O.z_impact_max/1e3], 'k:', ...
      (R1.z_tan-R1.slta)/1e3, R1.z_impact/1e3, '-', ... 
      (R2.z_tan-R2.slta)/1e3, R2.z_impact/1e3, '-.', ... 
      (R3.z_tan-R3.slta)/1e3, R3.z_impact/1e3, '--' );
axis( [ -1.5 1.5 0 O.z_impact_max/1e3] );
xlabel( 'Deviation to SLTA [km]', 'FontSize', fs )

legend( l(2:end), 'No ionosphere', 'Low solar activity', ...
                                    'High solar activity' );


