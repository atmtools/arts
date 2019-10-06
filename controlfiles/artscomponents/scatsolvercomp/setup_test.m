function setup_test


workfolder = 'testdata';


%- Clear-sky atmosphere
%
z  = [0:250:12e3]';
nz = length( z );
%
[p,t] = z2p_cira86( z, 10, 23 );
%    
xmlStore( fullfile(workfolder,'p_grid.xml'), p, 'Vector', 'binary' );
xmlStore( fullfile(workfolder,'t_field.xml'), t, 'Tensor3', 'binary' );
xmlStore( fullfile(workfolder,'z_field.xml'), z, 'Tensor3', 'binary' );
%
%- H2O 
%
rh            = 0.8;
h2o           = zeros( size(t) );
ind1          = find( t>=273.15 );
h2o(ind1)     = rh * e_eq_water(t(ind1))./p(ind1);
ind2          = find( t<273.15 );
h2o(ind2)     = rh * e_eq_ice(t(ind2))./p(ind2);
%
% VMR:
% N2 and O2 at zero H2O
n2 = 0.781;
o2 = 1 - n2;
%
vmr_field = zeros( 3, nz );   % Order is N2, O2, H2O
% 
vmr_field(3,:) = h2o';
fac            = 1- h2o;
vmr_field(1,:) = fac*n2;
vmr_field(2,:) = fac*o2;
%
xmlStore( fullfile(workfolder,'vmr_field.xml'), vmr_field, 'Tensor4', 'binary' );


%- RWC and IWC
%
[iwc,rwc] = deal( zeros( size( t ) ) );
%
rwc(ind1)        = 0.1e-3;
iwc(ind2(1:15))  = 0.1e-3;
iwc(ind2(16:24)) = 0.05e-3;
iwc(ind2(25))    = 0.01e-3;
%
xmlStore( fullfile(workfolder,'particle_bulkprop_names.xml'), {'RWC','IWC'}, ...
                                                          'ArrayOfString', 'binary' );
xmlStore( fullfile(workfolder,'particle_bulkprop_field.xml'), ...
                                                   [rwc';iwc'], 'Tensor4', 'binary' );


%- Scattering data
%
scat_data = [];
scat_meta = [];
%
% RWC
%
load( '~/Data/SSDB/StandardHabits/FullSet/LiquidSphere.mat' )
%
[S] = crop_s(S);
%
j = 0;
for i = unique( [1 : 10 : length(S), length(S)] );
  j = j + 1;
  scat_data{1}{j} = S(i);
  scat_meta{1}{j} = M(i);
end
%
% IWC
%
load( '~/Data/SSDB/StandardHabits/FullSet/8-ColumnAggregate.mat' )
%
[S] = crop_s(S);
%
j = 0;
for i = unique( [1 : 5 : length(S), length(S)] );
  j = j + 1;
  scat_data{2}{j} = S(i);
  scat_meta{2}{j} = M(i);
end
%
xmlStore( fullfile(workfolder,'scat_data.xml'), scat_data, ...
          'ArrayOfArrayOfSingleScatteringData', 'binary' );
xmlStore( fullfile(workfolder,'scat_meta.xml'), scat_meta, ...
          'ArrayOfArrayOfScatteringMetaData', 'binary' );



function [S] = crop_s(S)
  %
  it = 1:2:5;
  iv = zeros(3,1);
  %
  for i = 1 : length(S)
    iv(1)     = max( find( S(i).f_grid <= 31.5e9 ) ); 
    [~,iv(2)] = min( abs( S(i).f_grid - 165e9 ) ); 
    iv(3)     = min( find( S(i).f_grid>=666e9 ) ); 
    nz = length(S(i).za_grid); 
    if  nz == 181
      iz = 1 : 181;
    elseif  nz == 901
      iz = 1 : 2 : 901;
    else
      iz = 1:nz;
    end
    %
    S(i).f_grid       = S(i).f_grid(iv);
    S(i).f_grid/1e9
    S(i).T_grid       = S(i).T_grid(it); 
    S(i).za_grid      = S(i).za_grid(iz); 
    S(i).pha_mat_data = S(i).pha_mat_data(iv,it,iz,:,:,:,:); 
    S(i).ext_mat_data = S(i).ext_mat_data(iv,it);; 
    S(i).abs_vec_data = S(i).abs_vec_data(iv,it);;     
  end
return