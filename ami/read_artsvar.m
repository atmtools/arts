%------------------------------------------------------------------------
% NAME:    read_artsvar
%
%          Reads a ARTS variable from binary or ASCII data file.
%
%          Allowed file names are
%             basename.varname.aa   (ASCII)
%             basename.varname.ab   (binary)
%
%          The function picks the latest file for the variable of interest.
%          To force reading of a specific file format use FORCE.            
%
% FORMAT:  x = read_artsvar(basename,varname [,force])
%
% RETURN:  x           Read data.
% IN:      basename    ARTS basename 
%          varname     Variable name.
% OPTIONAL force       'a' or 'b' to force reading of ASCII or binary,
%                      respectively.
%------------------------------------------------------------------------


% HISTORY: 00.04.12  First version by Patrick Eriksson (PE).
%          00.11.10  Included binary files (PE)


function x = read_artsvar(basename,varname,force)


%=== Get ARTS data type
artstype = get_artstype(varname);


%== Create file names
aname = sprintf('%s.%s.aa',basename,varname);;
bname = sprintf('%s.%s.ab',basename,varname);;


%=== Select file type
if exist('force')
  if force(1) == 'a'
    name = aname;
  else
    name = bname;
  end
else
  aexist = exist(aname);
  bexist = exist(bname);
  if ~aexist & ~bexist
    error(['No file for ',varname,' is found.'])
  elseif aexist & ~bexist
    name = aname;
  elseif ~aexist & bexist
    name = bname;
  else
    if filedate(aname,1) < filedate(bname,1)
      name = bname;
    else
      name = aname;
    end
  end
end


%=== Read the data by using READ_DATAFILE
x    = read_datafile(name,artstype);
