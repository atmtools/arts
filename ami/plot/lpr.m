%------------------------------------------------------------------------
% NAME:    lpr
%
%          Alternative printing function.
%          If the figure size has been set by FIG_SIZE, the figure
%          will be printed by the given size.
%          The figure is centred on the page.
%          If the width of the figure is larger than the height, and
%          the width is close or larger than the paper width, the 
%          figure is printed in landscape mode instaed of portrait
%          mode.
%
%          Additional printing arguments can be specified as for the
%          the ordinary printing command (see PRINT). For example,
%          lpr figname -deps
%
% RETURN:  -
% IN:      -------- (Below optional parameters)
%          args   additional arguments to pass on to PRINT
%       
% USING:   -
%------------------------------------------------------------------------

% HISTORY: 990923  Created by Patrick Eriksson.
    

function lpr(args1,args2,args3,args4,args5)


punit  = get(gcf,'PaperUnit');
funit  = get(gcf,'Unit');
orient = get(gcf,'PaperOrient');

set(gcf,'PaperUnit','centi');
set(gcf,'Unit','centi');

psize = get(gcf,'PaperSize');
fsize = get(gcf,'PaperPosition');
fsize = fsize(3:4);


%=== Check if the figure shall be rotated
if (fsize(1)>fsize(2)) & (fsize(1)>0.8*psize(1))
% print in landscape mode
  set(gcf,'PaperOrient','Landscape');
  set(gcf,'PaperPos',[psize(1)/2-fsize(2)/2 psize(1)/2-fsize(2)/2 fsize(1) fsize(2)]);
else
% print in portrait mode
  set(gcf,'PaperOrient','Portrait');
  set(gcf,'PaperPos',[psize(1)/2-fsize(1)/2 psize(2)/2-fsize(2)/2 fsize(1) fsize(2)]);
end

set(gcf,'PaperUnit',punit);
set(gcf,'Unit',funit);


%= Do the printing
if nargin == 0
  print
else
  if ~exist('args1'), args1=[]; end
  if ~exist('args2'), args2=[]; end
  if ~exist('args3'), args3=[]; end
  if ~exist('args4'), args4=[]; end
  if ~exist('args5'), args5=[]; end
  eval(['print ',args1,' ',args2,' ',args3,' ',args4,' ',args5])
end


set(gcf,'PaperOrient',orient);
