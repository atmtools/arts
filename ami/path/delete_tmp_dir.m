%------------------------------------------------------------------------
% NAME:    delete_tmp_dir(tmparea,dirname)
%
%          Deletes a temporary directory and its content.
%          The tmparea is used to check that the given directory is in the
%          temporary directory.
%
% FORMAT:  delete_tmp_dir(tmparea,dirname)
%
% RETURN:  -
% IN:      tmparea     Full path of directory where to place the temporary
%                      directory.
%          dirname     Full path of the directory to be removed.
%------------------------------------------------------------------------

% HISTORY: 2000.12.14  Created by Patrick Eriksson.


function delete_tmp_dir(tmparea,dirname)


if ~strncmp( tmparea, dirname, length(tmparea) )
  error('The given directory is not in the area for temporary directories');
end


%=== The Matlab delete function was tested but no good and secury solution 
%=== could be found and system commands are used instead.
if isunix
  eval([ '!rm -r ',dirname ])
else
  error('Unknown computer type.');
end

