%-----------------------------------------------------------------------------
% NAME:     hide_axes
%
%           Hides the axes of a figure and sets the figure to cover the whole
%           of the plot window. 
%
% FORMAT:   hide_axes
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-10  Created by Patrick Eriksson


function hide_axis( h )


if nargin == 0
  h = gca;
end

set(h,'Position',[0 0 1 1]);
set(h,'Visible','Off');

