%------------------------------------------------------------------------
% NAME:    init
%
%          Initiliaze the ARTS Matlab interface (AMI).
%
% FORMAT:  init( [ warn ] )
%
% RETURN:  -
% IN:      -
% OPTIONAL warn   Boolean for display of warnings. Default is 1.
%                 0 suppresses warnings.
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function init( warn )


if nargin == 0
  warn = 1;
end


%=== Extend the Matlab search path to include the AMI dirs.
%=== It is assumed that this script is executed in the top AMI dir
addpath([pwd,'/arts']);
addpath([pwd,'/file']);
addpath([pwd,'/general']);
addpath([pwd,'/hmatrix']);
addpath([pwd,'/math']);
addpath([pwd,'/path']);
addpath([pwd,'/physics']);
addpath([pwd,'/plot']);
addpath([pwd,'/retrieval']);


%=== Copy the file src/wsv.txt to ami/arts
%
filein  = fullfile( fileparts( pwd ), fullfile( 'src', 'auto_wsv.txt' ) );
fileout = fullfile( pwd, fullfile( 'arts', 'auto_wsv.txt' ) );
%
if ~exist(fileout,'file')  |  ( filedate(fileout,1) < filedate(filein,1) )
  copyfile( filein, fileout );
end


%=== Define global physical constants
%
% To use the constants, type in the script for example:
% global EARTH_RADIUS
%
global EARTH_RADIUS RAD2DEG DEG2RAD PLANCK_CONST SPEED_OF_LIGHT
global BOLTZMAN_CONST AVOGADROS_NUMB COSMIC_BG_TEMP SUN_TEMP
global NAT_LOG_2 ATM2HPA
%
EARTH_RADIUS   = 6.378e6;
RAD2DEG        = 57.29577951308232;
DEG2RAD        = 0.01745329251994;
PLANCK_CONST   = 6.626180e-34;
SPEED_OF_LIGHT = 2.99792458e8;
BOLTZMAN_CONST = 1.380662e-23;
AVOGADROS_NUMB = 6.0220450e26;
COSMIC_BG_TEMP = 2.735;
SUN_TEMP       = 6000.0;
NAT_LOG_2      = 0.69314718055994;
ATM2HPA        = 1.01325e3;                                   



%=== Try to add the arts-data top folder to the search path
%
c = get_artsconfig('ARTS_DATA_PATH');
%
if isstr( c )  &  ~strcmp( c, '"no"' )
  % Remove "-chars
  c = c( 2:(length(c)-1) );
  addpath( c );
else
  if warn
    fprintf('WARNING, could not determine path of arts-data.\n');
  end
end 