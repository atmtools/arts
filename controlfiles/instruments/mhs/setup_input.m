function setup_input

% These variables end with _multi in arts, but the names are here shortened
% to the names used for single mixer cases
  
% MHS-B is a true DSB instrument, but one sideband must be assigned as 
% main band and the upper band is here selected.


% Channels to include:
%
ind = 16:20;


%--- Index of "channel 0"
%
i0 = 15;

%--- LO frequencies
%
lo = [ 89e9, 157e9, 183.311e9, 183.311e9, 190.311e9 ];
%
%xmlStore( 'mhs.lo.xml', lo(ind-i0), 'Vector' );


%--- Main sideband (see above)
%
sideband_mode = { 'upper', 'upper', 'upper', 'upper', 'upper' };
%
%xmlStore( 'mhs.sideband_mode.xml', {sideband_mode{ind-i0}}, 'ArrayOfString' );


%--- Sideband response
%
% Is equal everywhere but frequency grids must be adjusted to channel
% width and spacing.
%
Gtmplt.name      = 'Sideband response function';
Gtmplt.gridnames = { 'Frequency' };
Gtmplt.grids     = [];
Gtmplt.dataname  = 'Response';
Gtmplt.data      = [0.5 0.5];
%
[G{1:5}]         = deal( Gtmplt );
%
G{1}.grids = { 1400e6*[-1 1] };
G{2}.grids = G{1}.grids;
G{3}.grids = { 3.500e6*[-1 1] };
G{4}.grids = G{3}.grids;
G{5}.grids = {1.100e6*[-1 1] };
%
%xmlStore( 'mhs.sideband_response.xml', {G{ind-i0}}, 'ArrayOfGriddedField1' );
xmlStore( 'MHS.SIDEBAND_RESPONSE.XML', {G{ind-i0}}, 'ArrayOfGriddedField1' );


%--- Center position of backend channels
%
% Frequencies must match selection of "main band"
%
G = { [89.7e9], [157.7e9], [184.311e9], [186.311e9], [190.861e9] };
%
%xmlStore( 'mhs.f_backend.xml', {G{ind-i0}}, 'ArrayOfVector' );


%--- Backend channel response functions
%
Gtmplt.name   = 'Backend channel response function';
Gtmplt.data   = [1 1];
%
G             = [];
%
G{1}{1}       = Gtmplt;
G{1}{1}.grids = { [-700 700]*1e6 };
%
G{2}{1}       = G{1}{1};
%
G{3}{1}       = Gtmplt;
G{3}{1}.grids = { [-250 250]*1e6 };
%
G{4}{1}       = Gtmplt;
G{4}{1}.grids = { [-500 500]*1e6 };
%
G{5}{1}       = Gtmplt;
G{5}{1}.grids = { [-550 550]*1e6 };
%
%xmlStore( 'mhs.backend_channel_response.xml', {G{ind-i0}}, 'ArrayOfArrayOfGriddedField1' );


