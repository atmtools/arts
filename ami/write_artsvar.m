%------------------------------------------------------------------------
% NAME:    write_artsvar
%
%          Writes a ARTS variable to a file in ARTS format.
%          The data is written to a file called
%             basename.varname.a
%          See further WRITE_DATAFILE
%
% FORMAT:  write_artsvar(basename,varname,x [,prec])
%
% RETURN:  -
% IN:      basename    the ARTS basename
%          varname     variable name
%          x           the data to store
% OPTIONAL prec        number of decimals to use, default 6
%                      If PREC=0, integer values are assumed          
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function write_artsvar(basename,varname,x,prec)


%=== Create full file name
name = sprintf('%s.%s.am',basename,varname);


%=== Create heading text
heading = sprintf('This file contains the ARTS variable %s.',varname);

if nargin == 3
  write_datafile(name,x,heading);
else
  write_datafile(name,x,heading,prec);
end
