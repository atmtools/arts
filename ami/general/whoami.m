%------------------------------------------------------------------------
% NAME:    whoami
%
%          Returns the used identity as the Unix function with the same 
%          name.
%
%          This function uses system specific commands and so far the
%          only system handled is Unix/Linux.
%
% FORMAT:  uid = whoami
%
% RETURN:  uid   The user identity as a string.
%------------------------------------------------------------------------

% HISTORY: 2001.08.07  Created by Patrick Eriksson. 


function uid = whoami


%=== Unix/Linux
%
if strcmp(computer,'LNX86') | strcmp(computer,'GLNX86') | strcmp(computer,'SOL2')
  %
  [s,uid] = unix('whoami');
  %
  % Remove line feed
  uid = uid( 1 : (length(uid)-1) );


%=== Unknown system
%
else
  error('Unknown computer type.');
end
