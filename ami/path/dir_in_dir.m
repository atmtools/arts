%------------------------------------------------------------------------
% NAME:   dir_in_dir
%
%         Returns the full path for a directory in another directory.
%         The sub-directory is created if not exists.
%
% FORMAT: s = dir_in_dir( maindir, subdir )
%
% OUT:    s         Full path for SUBDIR
% IN:     maindir   A main folder.
%         subdir    A directory in MAINDIR.
%------------------------------------------------------------------------

% HISTORY: 2000.12.22  Created by Patrick Eriksson.


function s = dir_in_dir( maindir, subdir )

s = fullfile( maindir, subdir, '' );

if ~exist( s , 'dir' )

  if exist( s, 'file' )
    error(sprintf(...
              'Cannot create directory %s as file with that name exists.',s));
  end

  mkdir( maindir, subdir );

end

