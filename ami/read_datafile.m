%------------------------------------------------------------------------
% NAME:    read_datafile
%
%          Reads data from a file in ARTS data format.
%
%          The data is returned as a matrix or an array.
%          The arrays are indexed by using curly braces.
%          For example, to get matrix 2, type
%             a = x{2}
%
% FORMAT:  x = read_datafile(filename)
%
% RETURN:  x             the data.
% IN:      filename      full file name
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function x = read_datafile(filename)


%=== Open file for reading
fid = fopen(filename,'r');
if fid < 0
  error('Error while opening file for writing');
end


%== Read until line not begins with #
s    = fgets(fid);
while s(1) == '#'
  s    = fgets(fid);
end


%=== Read number of matrices
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
  s    = fgets(fid);
  while s(1) == '#'
    s    = fgets(fid);
  end
  dim = sscanf(s,'%f %f');
  %= Check size
  if length(dim) ~= 2
    error('Could not read matrix size');
  end
  %= Read the data
  c1 = fscanf(fid,'%c',1);
  while ( c1=='#' )
    while ( c1~=10 )
      c1 = fscanf(fid,'%c',1);
    end
    c1 = fscanf(fid,'%c',1);
  end
  fseek(fid,-1,'cof');
  x = fscanf(fid,'%f',[dim(2) dim(1)]);
  if size(x,1)*size(x,2) < dim(1)*dim(2)
    error('The given size and actual data do not match'); 
  end
  x = x';
  %= Read row brake
  fgets(fid);

return
