%------------------------------------------------------------------------
% NAME:    write_datafile
%
%          Writes data to a ASCII or binary file in ARTS format.
%
%          An ARTS ArrayOfMatrix is in Matlab given as a cell. See
%          READ_DATAFILE.
%
% FORMAT:  write_datafile(filename,x,artstype [,prec,heading,toascii])
%
% RETURN:  -
% IN:      filename    Full file name.
%          x           The data to store.
%          artstype    String with the ARTS data type, e.g. MATRIX.
%                      It is allowed to shorten "ARRAYof" to "AO", e.g.
%                      it is OK to type AOMATRIX
%                      No distinction between upper and lower case letters.
% OPTIONAL prec        Precision for numeric type.
%                      For binary files, PREC is 1-2 where
%                        1: store as float (4 bytes)
%                        2: store as double(8 bytes)
%                      For ASCII files, PREC gives the number of decimals 
%                      to use.
%                      PREC can be empty ([]).
%                      Default: 2 for binary, 6 for ASCII.
%                      The precision is neglected for strings.
%          heading     Heading text as character matrix (see STR2MAT).
%                      The text is only included in ASCII files.
%                      The function puts in '# ' at the start of each line. 
%                      The time for writing is put in automatically.
%                      HEADING can be empty ([], default).
%          toascii     Flag to indicate  ASCII or binary output.
%                      If not given and the file has extension ".aa",
%                      TOASCII is set to 1. With the extension ".ab",
%                      ISASCII is set to 0. Default is ISASCII=0.
%------------------------------------------------------------------------

% HISTORY: 00.04.12  First version by Patrick Eriksson (PE).
%          00.11.09  Included binary files (PE)


function write_datafile(filename,x,artstype,prec,heading,toascii)



%=== Check input and set default values
if ~ischar(filename)
  error('The file name must be given as a string');
end
if ~exist('toascii')
  %= Default value is 0
  toascii = 0;
  l   = length(filename);
  if l > 3
    if strncmp(filename((-2:0)+l),'.aa',3)
      toascii = 1;
    end
  end
end 
if ~exist('prec') | isempty(prec)
  if toascii
    prec = 6;
  else
    prec = 2;
  end 
end
if prec < 0
  error('The precision must be >= 0');
end
if ~exist('heading')
  heading = []; 
end


%=== Convert ARRAYof to AO and make artstype uppercase
if strncmp(artstype,'ARRAYof',7)
  artstype = ['AO',artstype(8:length(artstype))];
end
artstype = upper(artstype);


%=== Temporary solution to handle sparse
if strcmp(class(x),'sparse')
  x = full(x);
end


%=== Make some checks of the data to store
if strncmp(artstype,'AO',2)
  if ~iscell(x)
    error('Input is not an array.')
  end
end
if iscell(x)
  for i = 1:length(x)
    if (isstr(x{1})~=isstr(x{i})) | (isnumeric(x{1})~=isnumeric(x{i}))
      error('All objects of an array must either be numeric or strings');
    end
  end
end
if strcmp(artstype,'INDEX') 
  if x < 0
    error('Input is not an INDEX (negative).')
  end     
  if (x-floor(x)) ~= 0
    error('Input is not an INDEX (non-integer).')
  end     
end
if strcmp(artstype,'AOINDEX') 
  for i = 1:length(x)
    if x{i} < 0
      error('Input is not an ARRAYofINDEX (negative value(s)).')
    end     
    if (x{i}-floor(x{i})) ~= 0
      error('Input is not an ARRAYofINDEX (non-integer value(s)).')
    end     
  end
end
if strcmp(artstype,'VECTOR') 
  if min(size(x)) > 1
    error('Input is not a vector (it is a matrix).')
  end
  x = vec2col(x);     
end
if strcmp(artstype,'AOVECTOR') 
  for i = 1:length(x)
    if min(size(x{i})) > 1
      error('ARRAYofVECTOR cannot contain matrices.')
    end     
    x{i} = vec2col(x{i});     
  end
end



%============================================================================
%============================================================================
%===
%=== ASCII
%===
%============================================================================
%============================================================================
 
if toascii 
  
  %=== Open file for writing
  fid = fopen(filename,'w');
  if fid < 0
    error(sprintf('Error while opening file for writing:\n%s',filename));
  end
  
  %=== Print heading
  if ~isempty(heading)
    for i = 1:size(heading,1)
      fprintf(fid,'# %s\n',heading(i,:));
    end
  end
  fprintf(fid,'#\n');
  fprintf(fid,'# This file is created by the Matlab script WRITE_DATAFILE.\n');
  cl = clock;
  if cl(5)>9
    fprintf(fid,sprintf('# %d-%d-%d, %d.%d\n#\n',cl(1:5)));  
  else
    fprintf(fid,sprintf('# %d-%d-%d, %d.0%d\n#\n',cl(1:5)));  
  end

  %=== Print number of matrices or strings
  if strncmp(artstype,'AO',2)
    fprintf(fid,'%d\n',length(x));
  else
    fprintf(fid,'1\n');
  end
    
  %=== INDEX
  if strcmp(artstype,'INDEX')
    write_matrix(fid,x,0);

  %=== AOINDEX
  elseif strcmp(artstype,'AOINDEX') 
    n = length(x);
    a = zeros(n,1);
    for i = 1:n
      a(i) = x{i};
    end
    write_matrix(fid,a,0);

  %=== NUMERIC
  elseif strcmp(artstype,'NUMERIC')
    write_matrix(fid,x,prec);

  %=== VECTOR
  elseif strcmp(artstype,'VECTOR')
    write_matrix(fid,x,prec);

  %=== MATRIX
  elseif strcmp(artstype,'MATRIX') 
    write_matrix(fid,x,prec);

  %=== AOVECTOR and AOMATRIX
  elseif strcmp(artstype,'AOVECTOR') | strcmp(artstype,'AOMATRIX') 
    for i = 1:length(x)
      write_matrix(fid,x{i},prec);
    end

  %=== STRING
  elseif strcmp(artstype,'STRING')
    write_string(fid,x);

  %=== AOSTRING
  elseif strcmp(artstype,'AOSTRING')
    for i = 1:length(x)
      write_string(fid,x{i});
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

  %=== Open file for writing
  fid = binfile_open_out(filename);

  %=== INDEX
  if strcmp(artstype,'INDEX')
    binfile_write(fid,filename,'INDEX','SCALAR','INDEX',x,prec);

  %=== AOINDEX
  elseif strcmp(artstype,'AOINDEX') 
    binfile_write(fid,filename,'INDEXARRAY','ARRAY','INDEX',x,prec);

  %=== NUMERIC
  elseif strcmp(artstype,'NUMERIC')
    binfile_write(fid,filename,'NUMERIC','SCALAR','NUMERIC',x,prec);

  %=== VECTOR
  elseif strcmp(artstype,'VECTOR')
    binfile_write(fid,filename,'VECTOR','VECTOR','NUMERIC',x,prec);

  %=== MATRIX
  elseif strcmp(artstype,'MATRIX') 
    binfile_write(fid,filename,'MATRIX','MATRIX','NUMERIC',x,prec);

  %=== AOVECTOR
  elseif strcmp(artstype,'AOVECTOR')
    n = length(x);
    binfile_write(fid,filename,'N_VECTOR','SCALAR','INDEX',n,prec);
    for i = 1:n
     binfile_write(fid,filename,['VECTOR',int2str(i-1)],'VECTOR','NUMERIC',...
                                                                   x{i},prec);
    end

  %=== AOMATRIX
  elseif strcmp(artstype,'AOMATRIX') 
    n = length(x);
    binfile_write(fid,filename,'N_MATRIX','SCALAR','INDEX',n,prec);
    for i = 1:n
     binfile_write(fid,filename,['MATRIX',int2str(i-1)],'MATRIX','NUMERIC',...
                                                                   x{i},prec);
    end

  %=== STRING
  elseif strcmp(artstype,'STRING')
    binfile_write(fid,filename,'STRING','STRING','CHAR',x,prec);
    x = binfile_read_string(filename,fid,'STRING','STRING',0,0);

  %=== AOSTRING
  elseif strcmp(artstype,'AOSTRING')
    n = length(x);
    binfile_write(fid,filename,'N_STRING','SCALAR','INDEX',n,prec);
    for i = 1:n
      binfile_write(fid,filename,['STRING',int2str(i-1)],'STRING','CHAR',...
                                                                   x{i},prec);
    end

  %=== ??????????????
  else
    error(['Unknown ARTS data type (',artstype,')'])
  end

  %=== Close the file
  binfile_close(fid,filename);

end



%============================================================================
%============================================================================
%===
%=== ASCII help functions
%===
%============================================================================
%============================================================================

%=== Writes a matrix
function write_matrix(fid,x,prec)
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

%=== Writes a string
function write_string(fid,x)
  fprintf(fid,'%s\n',x);
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

function fid = binfile_open_out( filename )

  % Open the file for writing
  fid = hdfh('open',filename, 'DFACC_CREATE', 0 );
  if  fid < 0 
    error(['Cannot open ',filename,' for writing.'])
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


function binfile_write_size(filename,dataname,vdata_id,nrows,ncols)
  if  hdfvs('setattr',vdata_id,'vdata','SIZE',int32([nrows ncols])) < 0
    error(['Cannot write size data for ',dataname,' in file ',filename])
  end
return


function binfile_write(fid,filename,dataname,storagetype,atomictype,x,prec)

  %=== Get size
  nrows = size(x,1);
  ncols = size(x,2);

  %=== Create a new vdata
  vdata_id = hdfvs('attach',fid,-1,'w');
  if vdata_id < 0 
    error(['Cannot create a new vdata in file ',filename]);
  end
 
  %=== Set name of the vdata
  if hdfvs('setname',vdata_id,dataname) < 0
    error(['Cannot name the vdata ',dataname,' in file ',filename]);
  end

  %=== Write data size
  binfile_write_size( filename, dataname, vdata_id, nrows, ncols );

  
  %= Create the field
  if strcmp(atomictype,'INDEX')
    status1 = hdfvs('setclass',vdata_id,'UINT');
    status2 = hdfvs('fdefine',vdata_id,storagetype,'uint32',1);
    a = unit32(x);
    x = cell(1,1);
    x{1} = a;

  elseif strcmp(atomictype,'NUMERIC')
    if prec == 2
      status1 = hdfvs('setclass',vdata_id,'DOUBLE');
      status2 = hdfvs('fdefine',vdata_id,storagetype,'double',1);
    elseif prec == 1
      status1 = hdfvs('setclass',vdata_id,'FLOAT');
      status2 = hdfvs('fdefine',vdata_id,storagetype,'float',1);
      x = single(x);
    else
      error('When writing binary, PREC must be 1 or 2.');
    end
    a = x';
    x = cell(1,1);
    x{1} = a(:)';

  elseif strcmp(atomictype,'CHAR')
    status1 = hdfvs('setclass',vdata_id,'CHAR');
    status2 = hdfvs('fdefine',vdata_id,storagetype,'char',1);
    a = x;
    x = cell(1,1);
    x{1} = a;

  else
    error(['The atomic data type ',atomictype',' is not handled']);
  end

  %= Handle error
  if  status1 < 0
    error(['Cannot set class on ',dataname,' in file ',filename]);
  end
  if status2 < 0
    error(['Cannot create the field ',storagetype,' in file ',filename])
  end

  %= Finalize the definition of the field
  if hdfvs('setfields',vdata_id,storagetype) < 0
    error(['Cannot set the field ',storagetype,' in file ',filename])
  end

  %= Do actual writing
  %=== Write data (if not empty)
  if (nrows>0) & (ncols>0)
    nout = nrows*ncols;
    ndone = hdfvs('write',vdata_id, x );
    if ndone ~= nout 
      error(['Could not write all data to field ',storagetype,' in file ',filename,'. Out of memory?']);
    end
  end

  %=== Detach the vdata
  if hdfvs('detach',vdata_id) < 0
    error(['Cannot detach the vdata ',dataname,' in file ',filename]);
  end

return
