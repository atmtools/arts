%------------------------------------------------------------------------
% NAME:    arts
%
%          Runs ARTS from Matlab
%
% FORMAT:  arts cfile1 cfile2 cfile3 cfile4 cfile5
%             or
%          arts(cfile1,cfile2,cfile3,cfile4,cfile5)
%
% RETURN:  -
% IN:      cfile
%------------------------------------------------------------------------

% HISTORY: 12.04.00  Created by Patrick Eriksson. 

function arts(cfile1,cfile2,cfile3,cfile4,cfile5)


if nargin == 0
  eval(['!arts '])

elseif nargin == 1
  eval(['!arts ',cfile1])

elseif nargin == 2
  eval(['!arts ',cfile1,' ',cfile2])

elseif nargin == 3
  eval(['!arts ',cfile1,' ',cfile2,' ',cfile3])

elseif nargin == 4
  eval(['!arts ',cfile1,' ',cfile2,' ',cfile3,' ',cfile4])

elseif nargin == 5
  eval(['!arts ',cfile1,' ',cfile2,' ',cfile3,' ',cfile4,' ',cfile5])

else
  error('Sorry, but more than 5 text strings cannot be handled')
end

