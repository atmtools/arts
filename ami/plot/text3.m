%------------------------------------------------------------------------
% NAME:    text3
%
%          Writes a text to a mesh figure, i.e. 3-D (not a plot figure, 
%          2-D), giving the possibility to easily set the fontsize and 
%          the fontweight of the text.
%          Use TEXT2 for plot figures.
%
% RETURN:  h       handle to the text
% IN:      x       the x coordinate
%          y       the y coordinate
%          z       the z coordinate
%          s       the text (a string)
%          --------(Below optional parameters)
%          fsize   fontsize, default 12
%          bold    flag for bold text (0=normal,1=bold), default 1
% USING:   -
% HISTORY: 990922  Created by Patrick Eriksson.
%------------------------------------------------------------------------
    

function h = text3(x,y,z,s,fsize,bold)

if ~exist('fsize'), fsize = 12; end
if ~exist('bold'),  bold = 1; end


h = text(x,y,z,s,'FontSize',fsize);


if bold
  set(h,'FontWeight','bold');
else
  set(h,'FontWeight','normal');
end
