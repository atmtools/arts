%-----------------------------------------------------------------------------
% NAME:     arts_path
%
%           Returns the search path to the top directory of arts.
%
%           This function can be used to generate full names of files inside
%	    arts. The top folder of AMI must be in the Matlab search path.
%
% FORMAT:   p = arts_path
%
% OUT:      p   Path to the top directory of arts.
%-----------------------------------------------------------------------------

% HISTORY: 2002-08-15  Created by Patrick Eriksson


function p = arts_path

% We are taking two steps upwards in the folder ladder by calling
% fileparts two extra times. 

p = fileparts( fileparts( fileparts( which('arts_path') ) ) );

