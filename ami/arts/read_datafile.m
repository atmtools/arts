%------------------------------------------------------------------------
% NAME:    read_datafile
%
%          Reads data from a file in ARTS binary or ASCII data format.
%
%          The data is returned in the Matlab format closest to the
%          the ARTS format for the variable to read.
%
%          Note that arrays are indexed by using curly braces.
%          For example, to get matrix 2 of a matrix array, type
%             a = x{2}
%
% FORMAT:  x = read_datafile(filename,artstype [,isascii])
%
% RETURN:  x           The read data.
% IN:      filename    File name.
%          artstype    String with the ARTS data type, e.g. MATRIX.
%                      It is allowed to shorten "ARRAYof" to "AO", e.g.
%                      it is OK to type AOMATRIX.
%                      No distinction between upper and lower case letters.
% OPTIONAL isascii     Flag to indicate if the file holds ASCII or bin.
%                      If not given and the file has extension ".aa",
%                      ISASCII is set to 1. With the extension ".ab",
%                      ISASCII is set to 0. Default is ISASCII=0.
%------------------------------------------------------------------------

% HISTORY: 00.04.12  First version by Patrick Eriksson (PE).
%          00.11.09  Included binary files (PE)


function x = read_datafile(filename,artstype,isascii)


%== Check input
if ~ischar(filename)
  error('The file name must be given as a string');
end
if ~ischar(artstype)
  error('The ARTS data type name must be given as a string');
end


%=== Set default value for ISASCII
if ~exist('isascii'), 
  %= Default value is 0
  isascii = 0;
  l   = length(filename);
  if l > 3
    if strncmp(filename((-2:0)+l),'.aa',3)
      isascii = 1;
    end
  end
end 


%=== Convert ARRAYof to AO and make artstype uppercase
if strncmp(upper(artstype),'ARRAYOF',7)
  artstype = ['AO',artstype(8:length(artstype))];
end
artstype = upper(artstype);



%============================================================================
%============================================================================
%===
%=== ASCII
%===
%============================================================================
%============================================================================
 
if isascii 

  %=== Open file for reading
  fid = fopen(filename,'r');
  if fid < 0
    error(sprintf(...
  'Error while opening file for reading:\n%s\nDoes the file exist?',filename));
  end

  %== Read until line not begins with #
  s    = fgets(fid);
  while s(1) == '#'
    s    = fgets(fid);
  end

  %=== Read number of matrices or strings
  nmat = sscanf(s,'%f ');
  if isempty(nmat) | (nmat<1)
    error(['Could not read number of matrices or strings. ',...
           'Is it a valid ARTS ASCII file?']);
  end

  %=== INDEX and NUMERIC
  if strcmp(artstype,'INDEX') | strcmp(artstype,'NUMERIC') 
    if nmat > 1
      error('You wanted a scalar, but the file contains an array.')
    end
    x = read_matrix(fid);
    if any( size(x) ~= 1 )
      error('You wanted a scalar, but the file contains a vector or matrix.')
    end

  %=== VECTOR, MATRIX and AOINDEX
  elseif strcmp(artstype,'VECTOR') | strcmp(artstype,'MATRIX') | ...
         strcmp(artstype,'AOINDEX')
    if nmat > 1
      if ~strcmp(artstype,'AOINDEX')
        error('You wanted a vector/matrix, but the file contains an array.')
      else
        error(...
     'You wanted an index array, but the file contains another type of array.')
      end
    end
    x = read_matrix(fid);
    if strcmp(artstype,'VECTOR')
      if min(size(x)) > 1
        error('You wanted a vector, but the file contains a matrix.')
      end
      x = vec2col(x);
    elseif strcmp(artstype,'AOINDEX') 
      if min(size(x)) > 1
        error('You wanted an index array, but the file contains a matrix.')
      end
      a = vec2col(x);
      n = length(a);
      x = cell(n,1);
      for i = 1:n
        x{i} = a(i);
      end
    end

  %=== AOVECTOR and AOMATRIX
  elseif strcmp(artstype,'AOVECTOR') | strcmp(artstype,'AOMATRIX') 
    x = cell(nmat,1);
    for i = 1:nmat
      x{i} = read_matrix(fid);
      if strcmp(artstype,'AOVECTOR')
        if min(size(x{i})) > 1
          error(...
            'You wanted a vector array, but the file contains a matrix array.')
        end
        x{i} = vec2col(x{i});
      end
    end

  %=== STRING
  elseif strcmp(artstype,'STRING')
    if nmat > 1
      error('You wanted a string, but the file contains an array.')
    end
    x = read_string(fid);

  %=== AOSTRING
  elseif strcmp(artstype,'AOSTRING')
    x = cell(nmat,1);
    for i = 1:nmat
      x{i} = read_string(fid);
    end

  %=== ??????????????
  else
    error(['Unknown ARTS data type (',artstype,')'])
  end

  %=== Close the file
  fclose(fid);



%============================================================================
%============================================================================
%===
%=== Binary
%===
%============================================================================
%============================================================================

else

  %=== Make sure that all HDF files are closed
  %=== (these files are not closed when errors occur)
  hdfml('closeall');

  %=== Open file for reading
  fid = binfile_open_in(filename);

  %=== INDEX
  if strcmp(artstype,'INDEX')
    x = binfile_read_numeric(filename,fid,'INDEX','SCALAR',1,1);

  %=== AOINDEX
  elseif strcmp(artstype,'AOINDEX') 
    x = binfile_read_numeric(filename,fid,'INDEXARRAY','ARRAY',0,1);
    %a = binfile_read_numeric(filename,fid,'INDEXARRAY','ARRAY',0,1);
    %n = length(a);
    %x = cell(n,1);
    %for i = 1:n
    %  x{i} = a(i);
    %end

  %=== NUMERIC
  elseif strcmp(artstype,'NUMERIC')
    x = binfile_read_numeric(filename,fid,'NUMERIC','SCALAR',1,1); 

  %=== VECTOR
  elseif strcmp(artstype,'VECTOR')
    x = binfile_read_numeric(filename,fid,'VECTOR','VECTOR',0,1);

  %=== MATRIX
  elseif strcmp(artstype,'MATRIX') | strcmp(artstype,'SYMMETRIC') 
    x = binfile_read_numeric(filename,fid,'MATRIX','MATRIX',0,0);

  %=== AOVECTOR
  elseif strcmp(artstype,'AOVECTOR') 
    n = binfile_read_numeric(filename,fid,'N_VECTOR','SCALAR',1,1);
    x = cell(n,1);
    for i = 1:n
      x{i} = binfile_read_numeric(filename,fid,['VECTOR',int2str(i-1)],...
                                                                'VECTOR',0,0);
    end

  %=== AOMATRIX
  elseif strcmp(artstype,'AOMATRIX') 
    n = binfile_read_numeric(filename,fid,'N_MATRIX','SCALAR',1,1);
    x = cell(n,1);
    for i = 1:n
      x{i} = binfile_read_numeric(filename,fid,['MATRIX',int2str(i-1)],...
                                                                'MATRIX',0,0);
    end

  %=== STRING
  elseif strcmp(artstype,'STRING')
    x = binfile_read_string(filename,fid,'STRING','STRING',0,0);

  %=== AOSTRING
  elseif strcmp(artstype,'AOSTRING')
    n = binfile_read_numeric(filename,fid,'N_STRING','SCALAR',1,1);
    x = cell(n,1);
    for i = 1:n
      x{i} = binfile_read_string(filename,fid,['STRING',int2str(i-1)],...
                                                                'STRING',0,1);
    end

  %=== LOS
  elseif strcmp(artstype,'LOS')
    n = binfile_read_numeric(filename,fid,'N_LOS.P','SCALAR',1,1);
    x.p = cell(n,1);
    for i = 1:n
      x.p{i} = binfile_read_numeric(filename,fid,['LOS.P',int2str(i-1)],...
                                                                'VECTOR',0,0);
    end
    %
    n = binfile_read_numeric(filename,fid,'N_LOS.PSI','SCALAR',1,1);
    x.psi = cell(n,1);
    for i = 1:n
      x.psi{i} = binfile_read_numeric(filename,fid,...
                                       ['LOS.PSI',int2str(i-1)],'VECTOR',0,0);
    end
    %
    n = binfile_read_numeric(filename,fid,'N_LOS.Z','SCALAR',1,1);
    x.z = cell(n,1);
    for i = 1:n
      x.z{i} = binfile_read_numeric(filename,fid,['LOS.Z',int2str(i-1)],...
                                                                'VECTOR',0,0);
    end
    %
    x.l_step = binfile_read_numeric(filename,fid,'LOS.L_STEP','VECTOR',0,1);
    %
    x.ground = binfile_read_numeric(filename,fid,'LOS.GROUND','ARRAY',0,1);
    x.start  = binfile_read_numeric(filename,fid,'LOS.START','ARRAY',0,1);
    x.stop   = binfile_read_numeric(filename,fid,'LOS.STOP','ARRAY',0,1);


  %=== ??????????????
  else
    error(['Unknown ARTS data type (',artstype,')'])
  end

  %=== Close the file
  binfile_close(fid,filename);

end

%=== Some checks of read data
if strcmp(artstype,'INDEX') 
  if x < 0
    error('You wanted an INDEX, but the file value is negative.')
  end     
  if (x-floor(x)) ~= 0
    error('You wanted an INDEX, but the file value is not an integer.')
  end     
end
if strcmp(artstype,'AOSIZET') 
  for i = 1:length(x)
    if x(i) < 0
      error('You wanted INDEX values, but the file contains negative values.')
    end     
    if (x(i)-floor(x(i))) ~= 0
      error('You wanted an INDEX values, but the file contains non-integer.')
    end     
  end
end



%============================================================================
%============================================================================
%===
%=== ASCII help functions
%===
%============================================================================
%============================================================================

%=== Reads a matrix
function x = read_matrix(fid)
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


%=== Reads a string
function x = read_string(fid)
  x    = fgets(fid);
  %= Remove line break
  x    = x(1:(length(x)-1));
return




%============================================================================
%============================================================================
%===
%=== Binary help functions
%===
%============================================================================
%============================================================================

% The functions below are copies of the corresponding function in 
% src/file.cc. For description on input and output of the functions and 
% more comments, see file.cc.
%

function fid = binfile_open_in( filename )

  % Check if the file is HDF
  if hdfh('ishdf',filename) < 0 
    error(['The file (',filename,') is not a HDF file.'])
  end

  % Open the file for reading
  fid = hdfh('open',filename, 'DFACC_READ', 0 );
  if  fid < 0 
    error(['Cannot open ',filename,' for reading.'])
  end

  % Initialize the VS interface
  if hdfv('start',fid) < 0
    error(['Cannot initialize the VS interafce in file: ',filename])
  end
return


function binfile_close(fid,filename)

  % Terminate access to the VS interface
  if hdfv('end',fid) < 0 
    error(['Cannot terminate access to the VS interface in: ',filename]);
  end

  % Close the file
  if hdfh('close',fid) < 0
    error(['Cannot close ',filename,'.'])
  end
return


function [vdata_id,nrows,ncols] = ...
            binfile_read_init(fid,filename,dataname,storagetype,nrows0,ncols0)

  % Find the Vdata in the file
  vdata_ref = hdfvs('find',fid,dataname);
  if vdata_ref <= 0
    error(['Cannot find the data ',dataname,' in file ',filename,'. Maybe the file contains data of other type.']);
  end

  % Attach the Vdata  
  vdata_id = hdfvs('attach',fid,vdata_ref,'r');
  if vdata_id <= 0
    error(['Cannot attach the data ',dataname,' in file ',filename]);
  end

  % Get number of rows and columns
  v = hdfvs('getattr', vdata_id,'vdata',0);
  nrows = double(v(1));
  ncols = double(v(2));

  % Check if number of rows and columns are as expected
  if ((nrows0>0)&(nrows~=nrows0)) | ((ncols0>0)&(ncols~=ncols0))
    error('The data have not the expected size.');
  end

  % Set fields to read (if not empty)
  if (nrows > 0) & (ncols > 0)
    if hdfvs('setfields',vdata_id,storagetype) < 0
      error(['Cannot find the field ',storagetype,' in file ',filename,'. Maybe the file contains data of other type.'])
    end
  end
return


function binfile_read_end(vdata_id,filename,dataname)
  if hdfvs('detach',vdata_id) < 0
    error(['Cannot detach the field ',dataname,' in file ',filename])
  end
return


function type_in_file = binfile_get_datatype(vdata_id)
  [type_in_file,status] = hdfvs('getclass',vdata_id);
return


%=== Reads numeric data of any type and size
%=== This function corresponds roughly to binfile_read2 in file.cc.
function x = ...
         binfile_read_numeric(filename,fid,dataname,storagetype,nrows0,ncols0)

  [vdata_id,nrows,ncols] = ...
            binfile_read_init(fid,filename,dataname,storagetype,nrows0,ncols0);

  if (nrows > 0) & (ncols > 0)
    [a,count] = hdfvs('read',vdata_id,nrows*ncols);
    x = reshape(double(a{1}),ncols,nrows)';
  else
    x = [];
  end

  binfile_read_end(vdata_id,filename,dataname);
return


%=== Reads a string
%=== This function corresponds to binfile_read3 in file.cc.
function x = ...
         binfile_read_string(filename,fid,dataname,storagetype,nrows0,ncols0)

  [vdata_id,nrows,ncols] = ...
            binfile_read_init(fid,filename,dataname,storagetype,nrows0,ncols0);

  if (nrows > 0) & (ncols > 0)
    [a,count] = hdfvs('read',vdata_id,nrows*ncols);
    x = reshape(char(a{1}),ncols,nrows);
  else
    x = [];
  end

  binfile_read_end(vdata_id,filename,dataname);
return
