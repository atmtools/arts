%------------------------------------------------------------------------
% NAME:    title2
%
%          Writes the title to a figure, giving the possibility to
%          easily set the fontsize and the fontweight of the text.
%
% RETURN:  h       handle to the text
% IN:      s       the text (a string)
%          --------
%          fsize   fontsize, default 12
%          bold    flag for bold text (0=normal,1=bold), default 1
% USING:   -
% HISTORY: 990922  Created by Patrick Eriksson.
%------------------------------------------------------------------------
    

function h = title2(s,fsize,bold)

if ~exist('fsize'), fsize = 12; end
if ~exist('bold'),  bold = 1; end


h = title(s,'FontSize',fsize);


if bold
  set(h,'FontWeight','bold');
else
  set(h,'FontWeight','normal');
end
