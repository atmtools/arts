%------------------------------------------------------------------------
% NAME:    init
%
%          Initiliaze the ARTS Matlab interface (AMI).
%
% FORMAT:  init
%
% RETURN:  -
% IN:      -
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function init


%=== Extend the Matlab search path to include the AMI dirs.
%=== It is assumed that this script is executed in the top AMI dir
%addpath([pwd,'/retrieval']);


%=== Copy the file src/wsv.txt to ami
%eval('!cp ../src/auto_wsv.txt arts');


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
