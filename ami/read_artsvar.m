%------------------------------------------------------------------------
% NAME:    read_artsvar
%
%          Reads a ARTS variable.
%
%          The data is read from the file
%             basename.varname.am
%
%          See further READ_DATAFILE
%
% FORMAT:  x = read_artsvar(basename,varname)
%
% RETURN:  x           the data
% IN:      basename    the ARTS basename 
%          varname     variable name
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function x = read_artsvar(basename,varname)


%=== Create full file name
name = sprintf('%s.%s.am',basename,varname);


%=== Read the data by using READ_DATAFILE
x    = read_datafile(name);
