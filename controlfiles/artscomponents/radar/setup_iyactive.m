function setup_iyactive


R.sensor_z    = 100e3;
R.f_radar     = 94e9;
R.range_bins  = 0:500:10e3;
R.t_ref       = 273.15;
R.d_droplet   = 50e-6;
R.dbz0        = -30;
R.cbox_limits = [ 1 101 ];  % 1-based indexing


%- Create ARTS input
%
setup_radar( R, fullfile(pwd,'testdata') );


%- Store reference value
%
xmlStore( fullfile(pwd,'testdata','dbz_ref.xml'), R.dbz0, 'Numeric' );
