%------------------------------------------------------------------------
% NAME:    arts
%
%          Runs ARTS from Matlab
%
% FORMAT:  arts cfile1 cfile2 cfile3 cfile4 cfile5 cfile6 cfile7 cfile8
%             or
%          arts(cfile1,cfile2,cfile3,cfile4,cfile5,cfile6,cfile7,cfile8)
%
% RETURN:  -
% IN:      cfile
%------------------------------------------------------------------------

% HISTORY: 2000-04-12   Created by Patrick Eriksson. 
           2002-07-12   PE: Included computer specific calls (unix etc.)
	                and error handling after a hint from Dietrich Feist.


function arts(cfile1,cfile2,cfile3,cfile4,cfile5,cfile6,cfile7,cfile8)


if nargin == 0
  error('ARTS requieres some input argument (such as a control file).')
end

if nargin < 2, cfile2 = [];   end
if nargin < 3, cfile3 = [];   end
if nargin < 4, cfile4 = [];   end
if nargin < 5, cfile5 = [];   end
if nargin < 6, cfile6 = [];   end
if nargin < 7, cfile7 = [];   end
if nargin < 8, cfile8 = [];   end
if nargin > 8
  error('Sorry, but more than 8 text strings are not handled')
end


%= Unix
%
if isunix
  s = unix(['arts ',cfile1,' ',cfile2,' ',cfile3,' ',cfile4,' ',cfile5,...
                                            ' ',cfile6,' ',cfile7,' ',cfile8]);

%= DOS/Windows
%
elseif strcmp( computer, 'PCWIN' );
  s = dos(['arts ',cfile1,' ',cfile2,' ',cfile3,' ',cfile4,' ',cfile5,...
                                            ' ',cfile6,' ',cfile7,' ',cfile8]);

%= Other computer types
%
% No error handling here
%
else
  eval(['!arts ',cfile1,' ',cfile2,' ',cfile3,' ',cfile4,' ',cfile5,...
                                            ' ',cfile6,' ',cfile7,' ',cfile8]);
  s = 1;
end


if s 
  fprintf('\n\n');
  error('An ARTS error has ocurred (see above).');
end