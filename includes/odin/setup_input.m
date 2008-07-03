function setup_input
  
% Both Qsmr and Atmlab are needed to run this script.
  
% All storing de-activated (by %) to not cause any unnecessary svn commits.
% Active the storing of the files you want to change.
  
  
%- Input is taken from operational software: Qsmr
%
inpath = '/home/patrick/Projects/SMR/Qsmr/Input';


%- Nominal LO settings
%
lo = [ 497.88e9 548.502e9 ];
%
%xmlStore( 'lo.SM_AC2ab.xml', lo(1), 'Numeric', 'DOUBLE' );
%xmlStore( 'lo.SM_AC1e.xml', lo(2), 'Numeric', 'DOUBLE' );

  

%- Frequency grids
%
f = read_datafile( fullfile(inpath,'Fmono','200mK','f_mono.SM_AC2ab.aa'), ...
                                                                    'Vector' );
%xmlStore( 'f_mono.SM_AC2ab.xml', f, 'Vector', 'DOUBLE' );
%
f = read_datafile( fullfile(inpath,'Fmono','200mK','f_mono.SM_AC1e.aa'), ...
                                                                    'Vector' );
%xmlStore( 'f_mono.SM_AC1e.xml', f, 'Vector', 'DOUBLE' );



%- Antenna
%
% Odin-SMR is making a continous scanning, with changing integration times.
% The broadening is here considered by sub-function.
%
t = 0.875;   %Shortest integration time
%t = 1.875;   %Intermediate integration time
%t = 3.875;   %Longest integration time
%
A = read_datafile( fullfile(inpath,'Antenna','Lab','antenna.SM_AC2ab.aa'),...
                                                                    'Matrix' ); 
A = broaden_apattern(A,t);
%
G.name      = 'Antenna response function';
G.gridnames = { 'Polarisation', 'Frequency', 'Zenith angle', 'Azimuth angle' };
G.grids     = { {'1'}, 501.3e9, A(:,1), 0 };
G.dataname  = 'Response';
G.data(1,1,:,1) = A(:,2);
%
name = sprintf( 'antenna.SM_AC2ab.%.0fms.xml', 1e3*t );
xmlStore( name, G, 'GriddedField4' );
%
A = read_datafile( fullfile(inpath,'Antenna','Lab','antenna.SM_AC1e.aa'),...
                                                                    'Matrix' ); 
A = broaden_apattern(A,t);
%
G.name      = 'Antenna response function';
G.gridnames = { 'Polarisation', 'Frequency', 'Zenith angle', 'Azimuth angle' };
G.grids     = { {'1'}, 544.5e9, A(:,1), 0 };
G.dataname  = 'Response';
G.data(1,1,:,1) = A(:,2);
%
name = sprintf( 'antenna.SM_AC1e.%.0fms.xml', 1e3*t );
xmlStore( name, G, 'GriddedField4' );

return



%- Sideband response
%
% Mimics the model in qsmr_ssbfilter
%
t       = 290;   % Assumed satellite temperature [K]
%
r_min   = 10^-2.3;
t0      = 292;
l0      = 9.472e-3;
ct      = 1.26e-6;
sb_path = -410e-6;
%
fif     = 3.9e9;    % Centre IF
df      = 1.25e9;   % Width of IF band
f       = [ -fif + linspace(-df/2,df/2,30)'; 
             fif + linspace(-df/2,df/2,30)' ];
%
l = l0 + ct*(t-t0) + sb_path/2;
%
% Note that RF shall be used here (lo+f)
r = qsmr_r_diplexer2( lo(1), lo(1)+3.9e9, lo(1)+f, 0 ) .* ...     % LO part
    qsmr_r_diplexer( l, lo(1)+f, r_min );                         % Filter part
%
G.name      = 'Sideband response function';
G.gridnames = { 'Frequency' };
G.grids     = { f };
G.dataname  = 'Response';
G.data      = r;
%
%xmlStore( 'sideband.SM_AC2ab.xml', G, 'GriddedField1' );
%
%
r_min   = 10^-1.44;
t0      = 291;
l0      = 9.552e-3;
ct      = 0;
sb_path = -383.5e-6;
%
fif     = 4.0e9;    % Centre IF
df      = 0.85e9;   % Width of IF band
f       = [ -fif + linspace(-df/2,df/2,30)'; 
             fif + linspace(-df/2,df/2,30)' ];
%
l = l0 + ct*(t-t0) + sb_path/2;
r = qsmr_r_diplexer2( lo(2), lo(2)+3.9e9, lo(2)+f, 0 ) .* ...     % LO part
    qsmr_r_diplexer( l, lo(2)+f, r_min );                         % Filter part
%
G.grids     = { f };
G.data      = r;
%
%xmlStore( 'sideband.SM_AC1e.xml', G, 'GriddedField1' );  



%- Backend respons
%
A = read_datafile( fullfile(inpath,'Backend','finalresponse.2001khz.aa'),...
                                                                    'Matrix' ); 
%
G.name      = 'Backend channel response function';
G.gridnames = { 'Frequency' };
G.grids     = { A(:,1) };
G.dataname  = 'Response';
G.data      = A(:,2);
%
%xmlStore( 'backend_channel_response.xml', {G}, 'ArrayOfGriddedField1' );



%- Nominal position of spectrometer channels
%
f = 1e9 * [501.170:0.001:501.592,501.970:0.001:502.392];
%
%xmlStore( 'f_backend.SM_AC2ab.xml', f, 'Vector' );
%
f = 1e9 * [544.1:0.001:544.9];
%
%xmlStore( 'f_backend.SM_AC1e.xml', f, 'Vector' );







function A2 = broaden_apattern(A,t)
  
  %- Scanning speed is 750 m/s. Convert to an approximate angle.
  %
  R  = constants( 'EARTH_RADIUS' );
  zp = 600e3;
  %
  v = geomztan2za( R, zp, 20e3 ) - geomztan2za( R, zp, 20e3+750 ); 
  
  %- Select n "middle points" over relative zenith angle scanning range
  %
  n = 11;
  %
  dza = (v*t/2) * linspace( -1+1/n, 1-1/n, n ); 
  
  %- Include scanning by simply adding antenna pattern for each dza
  %
  p = 0;
  %
  for i = 1:n
    p = p + interp1( A(:,1)+dza(i), A(:,2), A(:,1) );
  end

  %- The interpolation generates some NaNs
  %- Normalise to max value
  %
  ind = find( ~isnan(p) );
  %
  A2 = [ A(ind,1) p(ind)/max(p(ind)) ];  
  
return