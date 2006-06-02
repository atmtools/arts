%------------------------------------------------------------------------
% NAME:     arts_version
%
%           Returns the ARTS version string.
%
%	    If there is no output argument, the string is printed to the 
%	    screen.
%
% FORMAT:   s = arts_version( [ artsname ] )
%
% RETURN:   s          String describing the ARTS version, e.g. 'arts-1.0.51'.
% OPTIONAL: artsname   Name of ARTS executable. Default is 'arts'.      
%------------------------------------------------------------------------

% HISTORY: 2002-12-03   Created by Patrick Eriksson. 


function s = arts_version( artsname )


if nargin == 0
   artsname = 'arts';
end


%= Unix
%
if isunix
  [r,s] = unix( [artsname,' -v'] );

%= Other computer types
%
else
   error('Only Unix type computers are handled.');
end


if r
  error('No ARTS executable with the given name could be find.');
end


if nargout
  s = s( 9 : length(s) );
  s = s( 1 : (min(find(s==' '))-1) );
else
  disp(s);
  clear s;
end
