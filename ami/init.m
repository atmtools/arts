%------------------------------------------------------------------------
% NAME:    init
%
%          Initiliaze the ARTS Matlab interface (AMI).
%          For the moment only the Matlab search path is extended.
%          The script must be exucuted from the AMI top directory.
%
% FORMAT:  init
%
% RETURN:  -
% IN:      -
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function init


%=== Extend the Matlab search path to include the AMI dirs.
%= As a start, it is assumed that this script is executed in the top AMI dir
dir0 = pwd;
path(dir0 ,path);

