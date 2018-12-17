function nd = setup_radar(R,workfolder)


%- Calculate scattering properties
%
f         = R.f_radar;
t         = R.t_ref   + [-5 0 5];
%
nf        = length( f );
nt        = length( t );
rfr_index = zeros( nf, 1 );
%
for i = 1 : nf
  for j = 1 : nt
    rfr_index(i,j) = sqrt( eps_water_liebe93( f(i), t(j) ) ); 
  end
end
%
ssp = mie_arts_scat_data( vec2col(f), t, rfr_index, [0:1:180]', R.d_droplet/2 );
%
scat_data{1}{1} = ssp;


%- Create pnd_field
%
nd           = 10^(R.dbz0/10) * (1e-3/R.d_droplet)^6;
%
pnd            = zeros( 1, diff(R.cbox_limits)+1 );
pnd(:,1:end-1) = nd;


%- Store to files
%
xmlStore( fullfile(workfolder,'sensor_pos.xml'), R.sensor_z, 'Matrix' );
xmlStore( fullfile(workfolder,'f_grid.xml'), R.f_radar, 'Vector' );
xmlStore( fullfile(workfolder,'range_bins.xml'), R.range_bins, 'Vector' );
xmlStore( fullfile(workfolder,'t_ref.xml'), R.t_ref, 'Numeric' );
xmlStore( fullfile(workfolder,'pnd_field.xml'), pnd, 'Tensor4' );
xmlStore( fullfile(workfolder,'scat_data.xml'), scat_data, ...
          'ArrayOfArrayOfSingleScatteringData', 'binary' );
xmlStore( fullfile(workfolder,'cbox_limits.xml'), ...
          {R.cbox_limits(1)-1,R.cbox_limits(2)-1}, 'ArrayOfIndex' );

