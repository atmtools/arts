% DEMO_PLANETS_LIMBCOMP   Demo comparing limb spectra for the planets
%
%   This function make use of the interface functions to the planet
%   atmosphere set-up templates. Atmospheric cases used are hard-coded.
%   Included planets are Venus, Earth, Mars and Jupiter.
%
%   Limb spectra are simulated, with the same tangent altitude for all four
%   planets. 
%
%   This is just a demonstration, and the results give just a rough
%   overview. For example, species continua can be a lacking in the results.
%   Some hard-coded settings assume a relatively high tangent altitude (for
%   Earth should be OK-ish down to about 20 km). Continuum absorption can be
%   lacking, or use large extrapolation from defined temperature range.
%
% FORMAT [f,Y] = demo_planets_limbcomp( [ztan,f0,df,nf] )
%
% OUT   f      Frequency vector for spectra.
%       Y      simulated spectra. Col1 = Venus, Col2 = Earth, Col3 = Mars
%              and Col4 = Jupiter.
% OPT   ztan   Tangent altitude. Default is 50 km.
%       f0     Centre frequency of the calculation. Default is 552.2 GHz.
%       df     Frequency width of the calculation. Default is 10 GHz.
%       nf     Number of frequency points. Default is 10e3.
%       workfolder   Default is to use a temporary folder for the
%                    calculations. If you specify a folder with this
%                    argument, you will find the control and data files 
%                    (for Jupiter) in this folder after the calculations 
%                    are finished.

% 2013-10-22   Created by Patrick Eriksson.


function [f,Y] = demo_planets_limbcomp(varargin)
%
[ztan,f0,df,nf,workfolder] = optargs( varargin, ...
                                          { 50e3, 552.2e9, 10e9, 10e3, [] } );
  
  
% Settings for the different planets:
 
%--- Venus -------------------------------------------------------------------
%
Av.atmo         = 3;               % Atmospheric scenario
%
% Species included: Everything (beside duplicating CO2-cont.)
Av.basespecies  = [ 0:2,4:15];      
Av.h2ospecies   = 1;               % Level of water vapour
Av.hdospecies   = 3;               % Level of HDO
Av.Necase       = [];              % Free elctrons anyhow ignored
Av.interp_order = 1;               % Linear interpolation of fields (higher
%                                    values risky
Av.pmin         = 1e-6;            % Min pressure to consider. This value
                                   % crops around 200 km

%--- Earth --------------------------------------------------------------------
%  
Ae.fascode_atm             = 'midlatitude-winter';
Ae.ABS_SPECIES(1).TAG{1}   = 'H2O';
Ae.ABS_SPECIES(1).TAG{2}   = 'H2O-ForeignContStandardType';
Ae.ABS_SPECIES(2).TAG{1}   = 'O2';
Ae.ABS_SPECIES(3).TAG{1}   = 'N2-SelfContStandardType';
Ae.ABS_SPECIES(4).TAG{1}   = 'O3';
Ae.ABS_SPECIES(5).TAG{1}   = 'CO';
Ae.ABS_SPECIES(6).TAG{1}   = 'N2O';
Ae.ABS_SPECIES(7).TAG{1}   = 'HNO3';
Ae.ABS_SPECIES(8).TAG{1}   = 'HOCl';
Ae.ABS_SPECIES(9).TAG{1}   = 'ClO';
Ae.ABS_SPECIES(10).TAG{1}  = 'H2O2';
Ae.ABS_SPECIES(11).TAG{1}  = 'NO';
% Some species not included as not in Fascode: OCS, HO2  

%--- Mars --------------------------------------------------------------------
%
Am.Ls           = 2;               % Season
Am.daytime      = 0;               % Day / night
Am.dust         = 0;               % Dust level
Am.solar        = 2;               % Solar activity
%
% Species included: Everything (beside duplicating CO2-cont.)
Am.basespecies  = [ 0:2 4:15 ];
Am.ch4species   = 1;               % Standard CH4
Am.h2ospecies   = 1;               % Include all water as one species
Am.Necase       = [];              % Free electrons anyhow ignored
Am.interp_order = 1;               % Linear interpolation of fields (higher
%                                    values risky)
Am.pmin         = 1e-3;            % This crops around 100 km

%--- Jupiter ------------------------------------------------------------------
%
Aj.atmo         = 0;               % Atmospheric scenario
%
% Species included: Everything
Aj.basespecies  = [ 0:13 ];      
Aj.h2ospecies   = 1;               % Level of water vapour
Aj.nh3species   = 1;               % Level of NH3
Aj.ch4species   = 1;               % All of CH4
Aj.h2species    = 0;               % Included HD
Aj.Necase       = [];              % Free elctrons anyhow ignored
Aj.interp_order = 1;               % Linear interpolation of fields (higher
%                                    values risky
Aj.pmin         = 5e-3;            % Min pressure to consider, about 500 km
Aj.pmax         = 1e5;             % Max pressure to consider.

%------------------------------------------------------------------------------ 
%------------------------------------------------------------------------------ 
%------------------------------------------------------------------------------ 


%- Set-up Q
%
Q0                       = qarts;
%
Q0.F_GRID                = linspace( f0-df/2, f0+df/2, nf )';
%
Q0.INCLUDES{end+1}       = 'Position for inserting planet include file';  
Q0.INCLUDES{end+1}       = fullfile( 'ARTS_INCLUDES', 'agendas.arts' ); 
Q0.INCLUDES{end+1}       = fullfile( 'ARTS_INCLUDES', 'continua.arts' );
%
Q0.ATMOSPHERE_DIM        = 1;
Q0.STOKES_DIM            = 1;
%
Q0.CLOUDBOX_DO           = false;
Q0.J_DO                  = false;
Q0.SENSOR_DO             = false;
%
Q0.YCALC_WSMS            = { 'yCalc' };
Q0.PPATH_LMAX            = 25e3;
Q0.IY_UNIT               = 'PlanckBT';
%
Q0.PPATH_AGENDA               = { 'ppath_agenda__FollowSensorLosPath'   };
Q0.PPATH_STEP_AGENDA          = { 'ppath_step_agenda__GeometricPath'    };
Q0.BLACKBODY_RADIATION_AGENDA = { 'blackbody_radiation_agenda__Planck'  };
Q0.IY_SPACE_AGENDA            = { 'iy_space_agenda__CosmicBackground'   };
Q0.IY_SURFACE_AGENDA          = { 'iy_surface_agenda__UseSurfaceRtprop' };
Q0.IY_MAIN_AGENDA             = { 'iy_main_agenda__Emission'            };
%
Q0.SENSOR_POS            = 1000e3;
%
Q0.Z_SURFACE             = 1e3;  % Here not important, as limb sounding
%
Q0.ABSORPTION            = 'OnTheFly';  % A simple, but slow, option!
Q0.ABS_LINESHAPE         = 'Voigt_Kuntz6';
Q0.ABS_LINESHAPE_CUTOFF  = -1;          % Turn off these to save some time
Q0.ABS_LINESHAPE_FACTOR  = 'no_norm';   % Not generally recommended!
%
f_extra = 3e9;   
Q0.ABS_WSMS{end+1} = sprintf( ['abs_linesReadFromSplitArtscat(abs_lines,',...
    'abs_species,"spectroscopy/Perrin/",%.2e,%.2e)'], ...
                                Q0.F_GRID(1)-f_extra, Q0.F_GRID(end)+f_extra );
Q0.ABS_WSMS{end+1} = 'abs_lines_per_speciesCreateFromLines';
%
Q0.ABS_WSMS{end+1} = [ 'ReadXML(abs_cia_data, "spectroscopy/cia/',...
                       'hitran2011/hitran_cia2012_adapted.xml.gz")' ];
% Qarts default for abs_xsec_agenda is not including CIA:
Q0.ABS_XSEC_AGENDA = { 'abs_xsec_per_speciesInit', ...
                       'abs_xsec_per_speciesAddLines', ... 
                       'abs_xsec_per_speciesAddConts', ...
                       'abs_xsec_per_speciesAddCIA(T_extrapolfac=3)' };



%- Set-up output arguments
%
f = Q0.F_GRID;
Y = repmat( NaN, length(f), 4 );


% Run Venus
%
Q = Q0;
%
Q.INCLUDES{1}  = fullfile( 'ARTS_INCLUDES', 'planet_venus.arts' );
Q.REFELLIPSOID = ellipsoidmodels( 'SphericalVenus' );
%
Q.SENSOR_LOS   = geomztan2za( Q.REFELLIPSOID(1), Q.SENSOR_POS, ztan );
%
Q = qarts_add_venus_planettbox( Av, Q, [] );
%
Y(:,1) = arts_y( Q, workfolder );
%
do_plot(f,Y)


% Run Earth
%
Q = Q0;
%
Q.INCLUDES{1}  = fullfile( 'ARTS_INCLUDES', 'planet_earth.arts' );
Q.REFELLIPSOID = ellipsoidmodels( 'SphericalEarth' );
%
Q.SENSOR_LOS   = geomztan2za( Q.REFELLIPSOID(1), Q.SENSOR_POS, ztan );
%
Q = qarts_add_fascode( Ae, Q, [] );
%
Y(:,2) = arts_y( Q, workfolder );
%
do_plot(f,Y)


% Run Mars
%
Q = Q0;
%
Q.INCLUDES{1}  = fullfile( 'ARTS_INCLUDES', 'planet_mars.arts' );
Q.REFELLIPSOID = ellipsoidmodels( 'SphericalMars' );
%
Q.SENSOR_LOS   = geomztan2za( Q.REFELLIPSOID(1), Q.SENSOR_POS, ztan );
%
Q = qarts_add_mars_planettbox( Am, Q, [] );
%
Y(:,3) = arts_y( Q, workfolder );
%
do_plot(f,Y)


% Run Jupiter
%
Q = Q0;
%
Q.INCLUDES{1}  = fullfile( 'ARTS_INCLUDES', 'planet_jupiter.arts' );
Q.REFELLIPSOID = ellipsoidmodels( 'SphericalJupiter' );
%
Q.SENSOR_LOS   = geomztan2za( Q.REFELLIPSOID(1), Q.SENSOR_POS, ztan );
%
Q = qarts_add_jupiter_planettbox( Aj, Q, [] );
%
Y(:,4) = arts_y( Q, workfolder );
%
do_plot(f,Y)


% Finish plot
%  
if ~nargout 
  axis( [ [min(f),max(f)]/1e9 0.8*min(Y(:)) 1.25*max(Y(:)) ] )
  title( sprintf('%.1f km tangent altitude (%d frequency points)',...
                                                              ztan/1e3, nf ) );
  grid
end

return



function do_plot(f,Y)
  if ~nargout 
    semilogy( f/1e9, Y );
    xlabel( 'Frequency [GHz]' );
    ylabel( 'Brightness temperature [Tb]' );
    legend( 'Venus', 'Earth', 'Mars', 'Jupiter', 'Location', 'North' );
    drawnow;
  end
return
