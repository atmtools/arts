%------------------------------------------------------------------------
% NAME:    temporary_directory
%
%          Creates an empty temporary directory with unique name.
%
% FORMAT:  tmpdir = temporary_directory(tmparea)
%
% RETURN:  tmpdir      Full path of temporary directory.
% IN:      tmparea     Full path of directory where to place the temporary
%                      directory.
%------------------------------------------------------------------------

% HISTORY: 2000.12.14  Created by Patrick Eriksson.


function tmpdir = temporary_directory(tmparea)


startdir = pwd;

if ~exist(tmparea,'dir')
  error(sprintf('The directory "%s" does not exist.',tmparea));  
end

cd( tmparea );

ready = 0;

while ~ready
  
  tmpdir = sprintf('Tmp%d',floor(1e6*rand(1)));

  if ~exist( tmpdir, 'dir' )
    mkdir(tmpdir);
    ready = 1;
  end

end

cd( startdir );

tmpdir = sprintf('%s/%s',tmparea,tmpdir);
