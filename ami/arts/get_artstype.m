%------------------------------------------------------------------------
% NAME:    get_artstype
%
%          Returns the ARTS data type for a workspace variable.
%          The file wsv.txt is used as data base. This file is
%          copied automatically to AMI by init.m.
%
% FORMAT:  typename = get_artstype(varname)
%
% RETURN:  typename    Type name
% IN:      varname     Variable name.
%------------------------------------------------------------------------


% HISTORY: 00.11.10  Created by Patrick Eriksson (PE).


function typename = get_artstype(varname)


%=== Open the file wsv.txt
fid = fopen('auto_wsv.txt','r');
if fid < 0
  error('Error while opening wsv.txt. Does the file exist?');
end


%=== Read until the given variable is found
%=== No empty lines are allowed between the "VARIABLE:" and "DATA TYPE:" lines
ok = 0;
s = fgets(fid);
while ~ok & isstr(s)
  s = noblanks(s);
  if strncmp(s,'VARIABLE:',9)
    s = s(10:(length(s)-1));
    if strcmp(s,varname)
      ok = 1;
      s = fgets(fid);
      s = noblanks(s);
      typename = s(10:(length(s)-1));
    end
  end
  s = fgets(fid);
end

if ~ok
  error(['The variable (',varname,') is not an ARTS workspace variable']);
end


%=== Close the file
fclose(fid);



function s = noblanks(s)
  ind = find(s~=' ');
  s = s(ind);
return
