%-----------------------------------------------------------------------------
% NAME:     get_artsconfig
%
%           Determines the configuration of arts.
%
%	    The function reads config.h in the top arts folder and returns
%	    the configuration setting as given in that file.
%
% FORMAT:   c = get_artsconfig( setting )
%
% OUT:      c         The confugure variable as given in config.h. For 
%                     example, '"/u/patrick/ARTS/arts-data"'.
% IN:       setting   String with name of setting, e.g. 'ARTS_DATA_PATH'.
%-----------------------------------------------------------------------------

% HISTORY: 2002-08-15  Created by Patrick Eriksson


function c = get_artsconfig( setting )

c = -1;

fid = fopen( fullfile( arts_path, 'config.h' ) );

if fid < 0
  error('The arts configure file (config.h) could not be found/opened.');
end

s = fgetl( fid );
ready = 0;

while ~ready & isstr(s)

  if strncmp( s, '#define ', 8 )

    % Remove first part
    s = s(9:length(s));

    if strncmp( s, setting, length(setting) )
    
      s = s( (length(setting)+2):length(s) );
      
      c = deblank(s);

      ready = 1;

    end

  end 

  s = fgetl( fid );
  
end


fclose(fid);

