function fig_refrtest


%= The ray tracing step to use
%
l = 1e3;

%------------------------------------------------------------------------------


%=== Create a temporary directory
%
tmpdir    = temporary_directory( '/tmp' );


%=== Set atmosphere to use
if strcmp( whoami, 'patrick' )
  Q.ATM = '/u/patrick/ARTS/arts-data/atmosphere/fascod/midlatitude-summer';
else
  error('Unknown user')
end


%=== Set Q fields to run ARTS
Q.ARTS       = 'arts';
Q.ARTS_LEVEL = 2;


%== Set the ray tracing step length to use
Q.L_STEP = l; 


%=== Create control file from template and run ARTS
[cfile,basename] = qtool( Q, tmpdir, 'fig_refrtest.tmplt' );
qp_arts( Q, cfile ); 


%=== Load data
z_tan = read_artsvar( basename, 'z_tan' );
za    = read_artsvar( basename, 'za_pencil' );
L     = read_artsvar( basename, 'los' );


%=== Pick out psi angel for tangent points
m   = length( L.psi );
psi = zeros( m, 1 );
for j = 1:m
  psi(j) = L.psi{j}(1);
end


%=== Calculate the alpha angel
a = 2 * ( 90 + psi - za );


%=== Plot
semilogx( a, z_tan/1e3 );
xlabel('Bending angle [deg]')
ylabel('Tangent altitude [km]')
axis([1e-3 2 0 40])
grid


%=== Delete the temporary directory
%
delete_tmp_dir( '/tmp', tmpdir );



