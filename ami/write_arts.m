%------------------------------------------------------------------------
% NAME:    write_arts
%
%          Writes data to a file in ARTS format.
%          The data is written to a file called
%             basename.varname.a
%          For data formats, see READ_ARTS
%
% FORMAT:  write_arts(basename,varname,x)
%
% RETURN:  -
% IN:      basename  The ARTS basename.
%          varname   Variable name.
%          x         The data to store.
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function write_arts(basename,varname,x,xtype,prec)


%=== Create full file name
name = sprintf('%s.%s.am',basename,varname);


%=== Open file for writing
fid = fopen(name,'w');
if fid < 0
  error('Error while opening file for writing');
end


%=== Print some comments and the data
fprintf(fid,'# This file is created by the Matlab script write_to_file\n');


%=== Write the data

%= Vector or matrix
if strcmp(class(x),'double')
  fprintf(fid,'1\n');
  write_mat(fid,x);

%= Array of matrices
else
  fprintf(fid,'%d\n',length(x));
  for i = 1:length(x)
    write_mat(fid,x{i});
  end

end  


%=== Close the file
fclose(fid);




%=== Sub-functions =================================================

function write_mat(fid,x)
  fprintf(fid,'%d  %d\n',size(x,1),size(x,2));
  for i = 1:size(x,1)
    fprintf(fid,'%.6e ',x(i,:));
    fprintf(fid,'\n');
  end
return
