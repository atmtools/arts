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
% OPTIONAL prec        number of digits to use, default 6          
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function write_artsvar(basename,varname,x,prec)


%=== Create full file name
name = sprintf('%s.%s.am',basename,varname);


if nargin == 3
  write_datafile(name,x);
else
  write_datafile(name,x,prec);
end
