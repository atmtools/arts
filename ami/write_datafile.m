%------------------------------------------------------------------------
% NAME:    write_datafile
%
%          Writes data to a file in ARTS format.
%
%          An ARTS ArrayOfMatrix is in Matlab given as a cell. See also
%          READ_DATAFILE.
%
% FORMAT:  write_datafile(filename,x, [heading,prec])
%
% RETURN:  -
% IN:      filename    full file name
%          x           the data to store
% OPTIONAL heading     heading text as character matrix (see STR2MAT)
%                      The function puts in '# ' at the start of each line. 
%                      Heading can be empty ([])
%          prec        number of decimals to use, default 6
%                      If PREC=0, integer values are assumed          
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 


function write_datafile(filename,x,heading,prec)


%=== Check input
if ~exist('prec')
  prec = 6;
end
if prec < 0
  error('The precision must be >= 0');
end


%=== Open file for writing
fid = fopen(filename,'w');
if fid < 0
  error('Error while opening file for writing');
end


%=== Print heading
if exist('heading') & ~isempty(heading)
  for i = 1:size(heading,1)
    fprintf(fid,'# %s\n',heading(i,:));
  end
end
fprintf(fid,'#\n');
fprintf(fid,'# (this file is created by the Matlab script write_datafile)\n');
fprintf(fid,'#\n');


%=== Write the data

%= Temporary solution to handle sparse
if strcmp(class(x),'sparse')
  x = full(x);
end

%= Character array
if strcmp(class(x),'char')
  error('Character arrays are not handled');

%= Character array
elseif strcmp(class(x),'struct')
  error('Structure arrays are not handled');

%= Vector or matrix
elseif strcmp(class(x),'double')
  fprintf(fid,'1\n');
  write_mat(fid,x,prec);

%= Array of matrices
elseif strcmp(class(x),'cell')
  fprintf(fid,'%d\n',length(x));
  for i = 1:length(x)
    write_mat(fid,x{i},prec);
  end

else
  error('Unknown data type');
end  


%=== Close the file
fclose(fid);




%=== Sub-functions =================================================

function write_mat(fid,x,prec)
  fprintf(fid,'%d  %d\n',size(x,1),size(x,2));
  if prec > 0
    for i = 1:size(x,1)
      eval(['fprintf(fid,''%.',int2str(prec),'e '',x(i,:));'])
      fprintf(fid,'\n');
    end
  else
    for i = 1:size(x,1)
      eval(['fprintf(fid,''%d '',x(i,:));'])
      fprintf(fid,'\n');
    end
  end
return
