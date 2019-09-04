function setup_wfun_test(do_mie)
%
if nargin < 1, do_mie = false; end


f_radar = 94e9;


%- Calculate scattering properties
%
if do_mie
  f         = f_radar+[-1 1]*1e9;
  t         = [193 233 273]';
  d         = [50e-6:100e-6:2e-3];
  %
  nf        = length( f );
  nt        = length( t );
  rfr_index = zeros( nf, 1 );
  %
  for i = 1 : nf
    for j = 1 : nt
      rfr_index(i,j) = sqrt( eps_ice_matzler06( f(i), t(j) ) ); 
    end
  end
  %
  for i = 1 : length(d)
    ssp = mie_arts_scat_data( vec2col(f), t, rfr_index, [0:1:180]', d(i)/2 );
    scat_data{1}{i} = ssp;
    scat_meta{1}{i}.version = 3;
    scat_meta{1}{i}.description = 'Ice sphere';    
    scat_meta{1}{i}.source      = 'Created by Atmlab''s Mie code';
    scat_meta{1}{i}.refr_index  = 'Matzler 2006';
    scat_meta{1}{i}.mass        = 917 * pi * d(i)^3 / 6;
    scat_meta{1}{i}.diameter_max                    = d(i);
    scat_meta{1}{i}.diameter_volume_equ             = d(i);
    scat_meta{1}{i}.diameter_area_equ_aerodynamical = d(i);
    
  end
  %
  xmlStore( fullfile(pwd,'testdata','scat_data_wfuntest.xml'), scat_data, ...
            'ArrayOfArrayOfSingleScatteringData', 'binary' );
  xmlStore( fullfile(pwd,'testdata','scat_meta_wfuntest.xml'), scat_meta, ...
            'ArrayOfArrayOfScatteringMetaData', 'binary' );
end


% Approx 0 to 16km, in steps of 250 m
p = logspace(5,4,65)';

% IWC profile is a "spike" at position 41
iwc     = zeros(size(p));
iwc(41) = 1e-3;


% Store to files
xmlStore( fullfile(pwd,'testdata','f_grid_wfuntest.xml'), f_radar, 'Vector' );
xmlStore( fullfile(pwd,'testdata','p_grid_wfuntest.xml'), p, 'Vector' );

xmlStore( fullfile(pwd,'testdata','particle_bulkprop_names.xml'), {'IWC'}, ...
                                                               'ArrayOfString' );
xmlStore( fullfile(pwd,'testdata','particle_bulkprop_field.xml'), ...
                                                   [iwc'], 'Tensor4' );
