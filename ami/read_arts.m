%------------------------------------------------------------------------
% NAME:    read_arts
%
%          Reads data from a ARTS output file.
%
%          The data is read from the file
%             basename.varname.a
%
%          The data is returned as a matrix or an array.
%          The arrays are indexed by using curly braces.
%          For example, to get matrix 2, type
%             a = x{2}
%
% FORMAT:  x = read_arts(basename,varname)
%
% RETURN:  x         The data.
% IN:      basename  The ARTS basename 
%          varname   Variable name.
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function x = read_arts(basename,varname)


%=== Create full file name
name = sprintf('%s.%s.a',basename,varname);


%=== Open file for reading
fid = fopen(name,'r');
if fid < 0
  error('Error while opening file for writing');
end


%== Read until line not begins with #
c1 = fscanf(fid,'%c',1);
while ( c1=='#' )
  while ( c1~=10 )
    c1 = fscanf(fid,'%c',1);
  end
  c1 = fscanf(fid,'%c',1);
end
fseek(fid,-1,'cof');


%=== Read number of matrices and size of first matrix
s    = fgets(fid);
nmat = sscanf(s,'%f ');
if isempty(nmat) | (nmat<1)
  error('Could not read number of matrices');
end


%=== Read the data, return matrix or cell
if nmat == 1
  x = read_mat(fid);
else
  x = cell(nmat,1);
  for i = 1:nmat
    x{i} = read_mat(fid);
  end
end


%=== Close the file
fclose(fid);




%=== Sub-functions =================================================

function x = read_mat(fid)
  %= Read size
  s   = fgets(fid);
  dim = sscanf(s,'%f %f');
  %= Check size
  if length(dim) ~= 2
    error('Could not read matrix size');
  end
  %= Read the data
  x = fscanf(fid,'%f',[dim(2) dim(1)]);
  x = x';
  %= Read row brake
  fgets(fid);

return
