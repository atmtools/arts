function mkfigs_refraction


%- Geometrical vs. refraction path
%
figure(1)
clf
%
nz       = 81;
nlat     = 21;
%
dim      = 2;
lat_grid = linspace(70,110,nlat);
lon_grid = [];
p_grid   = logspace( 5, 0, nz );
z_field  = linspace( 0, 80e3, nz )' * ones(1,nlat);
r_geoid  = 6378e3 * ones(nlat,1);
z_ground = 0 * ones(nlat,1);
%
a_pos = [ (6378+90)*1e3 70];
a_los = 99.42;
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, -1, z_ground, [] );
plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, '--' );
%
hold on
grid
xlabel('Latitude distance [degrees]');
ylabel('Altitude [km]');
set_labels( gca, 'FontSize', 12, 'FontWeight', 'bold' );
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
             a_pos, a_los, -1, z_ground, [], 2e3, 200+z_field*0, 0+z_field*0 );
plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, '-' );
%
a_los = 99
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
                                       a_pos, a_los, -1, z_ground, [] );
h(1) = plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, 'r--' );
%
hold on
%
P = arts_ppath( dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, ...
             a_pos, a_los, -1, z_ground, [], 2e3, 200+z_field*0, 0+z_field*0 );
h(2) = plot( P.pos(:,2)-P.pos(1,2), P.z/1e3, 'r-' );
%
hl = legend(h,'Geometric calculations','With refraction');
scale_axes( hl, 1.2 );




%- Gradients of refraction
%
figure(1)
clf
%
nz       = 33;
nlat     = 91;
%
dim      = 2;
lat_grid = linspace( 0, 90, nlat );
lon_grid = [];
p_grid   = logspace( 5, 3, nz );
[t,z,h2o]  = atmdata( p_grid );
t_field  = t' * ones(1,nlat);
z_field  = z' * ones(1,nlat);
vmr_field1 = zeros(1,nz,nlat,1);
vmr_field2 = zeros(1,nz,nlat,1);
vmr_field2(1,:,:,1) = h2o' * ones(1,nlat);
%
pv   = p_grid;
latv = lat_grid;
%
tmpfolder = create_tmpfolder;
cfile = fullfile( tmpfolder, 'cfile.arts' );
Q.TMPFOLDER = tmpfolder;
%
xmlStore( fullfile( tmpfolder, 'p_grid.xml' ), p_grid, 'Vector');
xmlStore( fullfile( tmpfolder, 'lat_grid.xml' ), lat_grid, 'Vector');
xmlStore( fullfile( tmpfolder, 't_field.xml' ), t_field, 'Tensor3');
xmlStore( fullfile( tmpfolder, 'z_field.xml' ), z_field, 'Tensor3');
xmlStore( fullfile( tmpfolder, 'vmr_field.xml' ), vmr_field1, 'Tensor4');
xmlStore( fullfile( tmpfolder, 'pv.xml' ), pv, 'Vector');
xmlStore( fullfile( tmpfolder, 'latv.xml' ), latv, 'Vector');
%
qtool( gradtemplate, cfile, Q );
call_fmodel( ['-r20 ',cfile] );
G1 = xmlLoad( fullfile( tmpfolder, 'G.xml' ) );
%
xmlStore( fullfile( tmpfolder, 'vmr_field.xml' ), vmr_field2, 'Tensor4');
call_fmodel( ['-r20 ',cfile] );
G2 = xmlLoad( fullfile( tmpfolder, 'G.xml' ) );
%
delete_tmpfolder( tmpfolder );
%
figure(2)
clf
semilogy( (squeeze(G1(1,:,nz))-1)*1e6, pv/100, 'r--', ...
          (squeeze(G2(1,:,nz))-1)*1e6, pv/100, 'b-' )
grid
set( gca, 'Ydir', 'rev' );
xlabel('Refractivity [ppm]');
ylabel('Pressure [hPa]');
set_labels( gca, 'FontSize', 12, 'FontWeight', 'bold' );
%
figure(3)
clf
semilogy( squeeze(G1(2,:,nz)), pv/100, 'r--', ...
          squeeze(G2(2,:,nz)), pv/100, 'b-' )
grid
set( gca, 'Ydir', 'rev' );
xlabel('\partialn/\partialr [1/m]');
ylabel('Pressure [hPa]');
set_labels( gca, 'FontSize', 12, 'FontWeight', 'bold' );
axis tight
ax = axis;
axis( [ax(1) 0 ax(3:4)] );
%
figure(4)
clf
[c,h] = contour( latv, pv/100, -squeeze(G2(3,:,:)), [1e-14 2e-13 5e-13, ...
             1e-12 2e-12 5e-12 1e-12 2e-11 5e-11 1e-11 2e-10 5e-10 1e-10]  );
set( gca, 'YScale', 'log' );
set( gca, 'Ydir', 'rev' );
clabel( c, h, 'manual' );
%clabel( c, h, 'FontSize', 10, 'Rotation', 0, 'LabelSpacing', 250 );
xlabel('Latitude [degree]');
ylabel('Pressure [hPa]');
set_labels( gca, 'FontSize', 12, 'FontWeight', 'bold' );

keyboard

return

if yes_or_no('Print figures')
  print ppath_dndlat.eps -depsc
  ! epstopdf ppath_dndlat.eps
  ! rm ppath_dndlat.eps
  figure(3)
  print ppath_dndr.eps -depsc
  ! epstopdf ppath_dndr.eps
  ! rm ppath_dndr.eps
  figure(2)
  print ppath_N.eps -depsc
  ! epstopdf ppath_N.eps
  ! rm ppath_N.eps
  figure(1)
  print ppath_refr1.eps -depsc
  ! epstopdf ppath_refr1.eps
  ! rm ppath_refr1.eps
end

 


function S = gradtemplate

S = { ...
'Main {', ...
'output_file_formatSetBinary{}', ...
'AtmosphereSet2D {}', ...
'ReadXML(p_grid){"$Q.TMPFOLDER$$Q.FILESEP$p_grid.xml"}', ...
'ReadXML(lat_grid){"$Q.TMPFOLDER$$Q.FILESEP$lat_grid.xml"}', ...
'ReadXML(z_field){"$Q.TMPFOLDER$$Q.FILESEP$z_field.xml"}', ...
'ReadXML(t_field){"$Q.TMPFOLDER$$Q.FILESEP$t_field.xml"}', ...
'ReadXML(vmr_field){"$Q.TMPFOLDER$$Q.FILESEP$vmr_field.xml"}', ...
'gas_speciesSet{species = [ "H2O" ]}', ...
'r_geoidWGS84{}', ...
'MatrixSetTakingSizeFromMatrix(z_ground,r_geoid){value = 0}', ...
'ReadXML(vector_1){"$Q.TMPFOLDER$$Q.FILESEP$pv.xml"}', ...
'ReadXML(vector_2){"$Q.TMPFOLDER$$Q.FILESEP$latv.xml"}', ...
'AgendaSet( refr_index_agenda ) {', ...
'    refr_indexThayer{}', ...
'}', ...
'RefrIndexFieldAndGradients(tensor4_1,vector_1,vector_2,vector_1){}', ...
'WriteXML(tensor4_1){"$Q.TMPFOLDER$$Q.FILESEP$G.xml"}', ...
'}'
};



function [t,z,h20] = atmdata( p )

% Profiles copied from FASCODE midlatitude winter

T = [
1.013000e+05 2.942000e+02 0.000000e+00
9.020000e+04 2.897000e+02 1.000000e+03
8.020000e+04 2.852000e+02 2.000000e+03
7.100000e+04 2.792000e+02 3.000000e+03
6.280000e+04 2.732000e+02 4.000000e+03
5.540000e+04 2.672000e+02 5.000000e+03
4.870000e+04 2.612000e+02 6.000000e+03
4.260000e+04 2.547000e+02 7.000000e+03
3.720000e+04 2.482000e+02 8.000000e+03
3.240000e+04 2.417000e+02 9.000000e+03
2.810000e+04 2.353000e+02 1.000000e+04
2.430000e+04 2.288000e+02 1.100000e+04
2.090000e+04 2.223000e+02 1.200000e+04
1.790000e+04 2.158000e+02 1.300000e+04
1.530000e+04 2.157000e+02 1.400000e+04
1.300000e+04 2.157000e+02 1.500000e+04
1.110000e+04 2.157000e+02 1.600000e+04
9.500000e+03 2.157000e+02 1.700000e+04
8.120000e+03 2.168000e+02 1.800000e+04
6.950000e+03 2.179000e+02 1.900000e+04
5.950000e+03 2.192000e+02 2.000000e+04
5.100000e+03 2.204000e+02 2.100000e+04
4.370000e+03 2.216000e+02 2.200000e+04
3.760000e+03 2.228000e+02 2.300000e+04
3.220000e+03 2.239000e+02 2.400000e+04
2.770000e+03 2.251000e+02 2.500000e+04
1.907000e+03 2.284500e+02 2.750000e+04
1.320000e+03 2.337000e+02 3.000000e+04
9.300000e+02 2.390000e+02 3.250000e+04
6.520000e+02 2.452000e+02 3.500000e+04
4.640000e+02 2.513000e+02 3.750000e+04
3.330000e+02 2.575000e+02 4.000000e+04
2.410000e+02 2.637000e+02 4.250000e+04
1.760000e+02 2.699000e+02 4.500000e+04
1.290000e+02 2.752000e+02 4.750000e+04
9.510000e+01 2.757000e+02 5.000000e+04
5.150000e+01 2.693000e+02 5.500000e+04
2.720000e+01 2.571000e+02 6.000000e+04
1.390000e+01 2.401000e+02 6.500000e+04
6.700000e+00 2.181000e+02 7.000000e+04
3.000000e+00 1.961000e+02 7.500000e+04
1.200000e+00 1.741000e+02 8.000000e+04
4.480000e-01 1.651000e+02 8.500000e+04
1.640000e-01 1.650000e+02 9.000000e+04
6.200000e-02 1.783000e+02 9.500000e+04];


H2O = [
1.013000e+05 1.877431e-02
9.020000e+04 1.379119e-02
8.020000e+04 9.687274e-03
7.100000e+04 5.988691e-03
6.280000e+04 3.815319e-03
5.540000e+04 2.226857e-03
4.870000e+04 1.510684e-03
4.260000e+04 1.020322e-03
3.720000e+04 6.466887e-04
3.240000e+04 4.132300e-04
2.810000e+04 2.474165e-04
2.430000e+04 9.562926e-05
2.090000e+04 2.945937e-05
1.790000e+04 8.006502e-06
1.530000e+04 5.004493e-06
1.300000e+04 3.401986e-06
1.110000e+04 3.302817e-06
9.500000e+03 3.200749e-06
8.120000e+03 3.152610e-06
6.950000e+03 3.202477e-06
5.950000e+03 3.301676e-06
5.100000e+03 3.452385e-06
4.370000e+03 3.601543e-06
3.760000e+03 3.852613e-06
3.220000e+03 4.001522e-06
2.770000e+03 4.203033e-06
1.907000e+03 4.452612e-06
1.320000e+03 4.703154e-06
9.300000e+02 4.854007e-06
6.520000e+02 4.953077e-06
4.640000e+02 5.002648e-06
3.330000e+02 5.103407e-06
2.410000e+02 5.304241e-06
1.760000e+02 5.454255e-06
1.290000e+02 5.505127e-06
9.510000e+01 5.503745e-06
5.150000e+01 5.353522e-06
2.720000e+01 5.003623e-06
1.390000e+01 4.402591e-06
6.700000e+00 3.703458e-06
3.000000e+00 2.953037e-06
1.200000e+00 2.101321e-06
4.480000e-01 1.331086e-06
1.640000e-01 8.505575e-07
6.200000e-02 5.447699e-07];


t   = interp1(log10(T(:,1)),T(:,2),log10(p));
z   = interp1(log10(T(:,1)),T(:,3),log10(p));
h20 = interp1(log10(H2O(:,1)),H2O(:,2),log10(p));
