%------------------------------------------------------------------------
% NAME:    write_artsvar
%
%          Writes data to ASCII or binary ARTS file.
%
%          Default file format for MATRIX, ARRAYofVECTOR and ARRAYofMATRIX
%          is binary. For other data type, the default is ASCII.
%
%          See further READ_DATAFILE
%
% FORMAT:  write_artsvar(basename,varname,x,force,heading,prec)
%
% RETURN:  -
% IN:      basename    ARTS basename.
%          varname     Variable name.
%          x           The data to store.
% OPTIONAL force       'a' or 'b' to force writing to ASCII or binary,
%                      respectively.
%          heading     Heading text as character matrix (see STR2MAT).
%                      See further WRITE_DATAFILE ([] is OK).
%          prec        Precision for numeric type.
%                      See further WRITE_DATAFILE.
%------------------------------------------------------------------------

% HISTORY: 00.04.12  First version by Patrick Eriksson (PE).
%          00.11.09  Included binary files (PE).


function write_artsvar(basename,varname,x,force,heading,prec)


%=== Get ARTS data type
artstype = get_artstype(varname);


%== Create file names
aname = sprintf('%s.%s.aa',basename,varname);;
bname = sprintf('%s.%s.ab',basename,varname);;


%=== Convert ARRAYof to AO and make artstype uppercase
if strncmp(artstype,'ARRAYof',7)
  shorttype = ['AO',artstype(8:length(artstype))];
end
shorttype = upper(artstype);


%=== Select file type
if exist('force')
  if force(1) == 'a'
    name = aname;
  else
    name = bname;
  end
else
  if strcmp(shorttype,'MATRIX') | strcmp(shorttype,'AOMATRIX') | ...
     strcmp(shorttype,'AOVECTOR')
    name = bname;
  else
    name = aname;
  end
end


%=== Create heading text
if ~exist('heading','var')
  heading = sprintf('This file contains the ARTS variable %s.',varname);
end

  
%=== Write data
if exist('prec')
  write_datafile(name,x,artstype,prec,heading);
else
  write_datafile(name,x,artstype,[],heading);
end

