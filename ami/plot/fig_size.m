%------------------------------------------------------------------------
% NAME:    fig_size
%
%          Gives the current figure the given size, both on the screen
%          and for later printing. The printing is best done with the 
%          accompanying function LPR.
%          The global parameter SCREENFAC, set in STARTUP, is used for 
%          a correct scaling for the screen.
%          The figure on the screen can be scaled further by the optional 
%          variable SFAC (the paper figure is not changed).
%          Default unit is centimeters.
%
% RETURN:  xsize   horisontal size
% IN:      ysize   vertical size
% IN:      -------- (Below optional parameters)
%          unit    unit for windowsize, default centimeters
%                  UNIT = 'cm' or 'mm' are OK arguments
%          sfac    scaling factor for the screen figure
% USING:   SCREENFAC
%------------------------------------------------------------------------

% HISTORY: 990927  Created by Patrick Eriksson.
    

function fig_size(xsize,ysize,unit,sfac)


if ~exist('unit'), unit = 'centimeters'; end
if ~exist('sfac'), sfac = 1.0; end


if strcmp(unit,'cm'), unit = 'centimeters'; end

if strcmp(unit,'mm'), 
  unit  = 'centimeters'; 
  xsize = xsize / 10;
  ysize = ysize / 10;
end


global SCREENFAC


h     = gcf;

unit1 = get(h,'Unit');
unit2 = get(h,'PaperUnit');

set(h,'Unit',unit);
set(h,'PaperUnit',unit);

pos   = get(h,'Position');
set(h,'Position',[pos(1:2) xsize*SCREENFAC*sfac ysize*SCREENFAC*sfac])

pos   = get(h,'PaperPosition');
set(h,'PaperPosition',[pos(1:2) xsize ysize])

set(h,'Unit',unit1);
set(h,'PaperUnit',unit2);
