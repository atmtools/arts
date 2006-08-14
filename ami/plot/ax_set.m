%------------------------------------------------------------------------
% NAME:    ax_set
%
%          Sets objects of the current axes and childrens.
%
%          The setting structure A can have the following fields:
%
%            POSITION
%            FONTSIZE
%            FONTWEIGHT
%            GRID
%            BOX
%            TITLE
%            TITLE_SIZE
%            TITLE_WEIGHT
%            XLABEL
%            XLABEL_SIZE
%            XLABEL_WEIGHT
%            XGRID
%            XTICK
%            XTICKLABEL
%            YLABEL
%            YLABEL_SIZE
%            YLABEL_WEIGHT
%            YGRID
%            YTICK
%            YTICKLABEL
%
%          Fields with a '_' (underscore) refers to properties of childrens.
%          For example, TITLE_SIZE is the size of the title, which is a
%          children to the axes.
%
% FORMAT:  h = ax_set( A )
%
% OUT:     h     Handle to the current axes.
% IN:      A     Setting structure.
%------------------------------------------------------------------------

% HISTORY: 2001.02.25  Created by Patrick Eriksson.
%          2005.05.26  Adding more fields XTICK and YTICK (CJ) 

function h = ax_set( A )


h = gca;


%= General stuff
%
if isfield( A, 'POSITION' )
  set( h, 'Position', A.POSITION );
end
%
if isfield( A, 'FONTSIZE' )
  set( h, 'FontSize', A.FONTSIZE );
end
%
if isfield( A, 'FONTWEIGHT' )
  set( h, 'FontWeight', A.FONTWEIGHT );
end
%
if isfield( A, 'GRID' )
  set( h, 'XGrid', A.GRID );
  set( h, 'YGrid', A.GRID );
end
%
if isfield( A, 'BOX' )
  set( h, 'Box', A.BOX );
end



%= Title
%
if isfield( A, 'TITLE' )
  set( get(h,'Title'), 'String', A.TITLE );
end
%
if isfield( A, 'TITLE_SIZE' )
  set( get(h,'Title'), 'FontSize', A.TITLE_SIZE );
end
%
if isfield( A, 'TITLE_WEIGHT' )
  set( get(h,'Title'), 'FontWeight', A.TITLE_WEIGHT );
end



%= X axis
%
if isfield( A, 'XLABEL' )
  set( get(h,'Xlabel'), 'String', A.XLABEL );
end
%
if isfield( A, 'XLABEL_SIZE' ) 
  set( get(h,'Xlabel'), 'FontSize', A.XLABEL_SIZE );
end
%
if isfield( A, 'XLABEL_WEIGHT' ) 
  set( get(h,'Xlabel'), 'FontWeight', A.XLABEL_WEIGHT );
end
%
if isfield( A, 'XGRID' )
  set( h, 'XGrid', A.XGRID );
end
%
if isfield( A, 'XTICK' )
  set( h, 'XTick', A.XTICK );
end
%
if isfield( A, 'XTICKLABEL' ) 
  set( h, 'XTickLabel', A.XTICKLABEL );
end



%= Y axis
%
if isfield( A, 'YLABEL' )
  set( get(h,'Ylabel'), 'String', A.YLABEL );
end
%
if isfield( A, 'YLABEL_SIZE' )
  set( get(h,'Ylabel'), 'FontSize', A.YLABEL_SIZE );
end
%
if isfield( A, 'YLABEL_WEIGHT' )
  set( get(h,'Ylabel'), 'FontWeight', A.YLABEL_WEIGHT );
end
%
if isfield( A, 'YGRID' )
  set( h, 'YGrid', A.YGRID );
end
%
if isfield( A, 'YTICK' )
  set( h, 'YTick', A.YTICK );
end

%
if isfield( A, 'YTICKLABEL' )
  set( h, 'YTickLabel', A.YTICKLABEL );
end

