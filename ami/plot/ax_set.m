%------------------------------------------------------------------------
% NAME:    ax_set
%
%          Sets objects of the current axes and childrens.
%
% FORMAT:  h = ax_set( A )
%
% OUT:     h     Handle to the current axes.
% IN:      A     Setting structure.
%------------------------------------------------------------------------

% HISTORY: 2001.02.25  Created by Patrick Eriksson.


function h = ax_set( A )


h = gca;


%= General stuff
%
if isfield( A, 'POSITION' ) & ~isfield( A.POSITION )
  set( h, 'Position', A.POSITION );
end
%
if isfield( A, 'FONTSIZE' ) & ~isfield( A.FONTSIZE )
  set( h, 'FontSize', A.FONTSIZE );
end
%
if isfield( A, 'FONTWEIGHT' ) & ~isfield( A.FONTWEIGHT )
  set( h, 'FontWeight', A.FONTWEIGHT );
end
%
if isfield( A, 'GRID' ) & ~isfield( A.GRID )
  set( h, 'XGrid', A.GRID );
  set( h, 'YGrid', A.GRID );
end
%
if isfield( A, 'BOX' ) & ~isfield( A.BOX )
  set( h, 'Box', A.BOX );
end



%= Title
%
if isfield( A, 'TITLE' ) & ~isfield( A.TITLE )
  set( get(h,'Title'), 'String', A.TITLE );
end
%
if isfield( A, 'TITLE_SIZE' ) & ~isfield( A.TITLE_SIZE )
  set( get(h,'Title'), 'FontSize', A.TITLE_SIZE );
end
%
if isfield( A, 'TITLE_WEIGHT' ) & ~isfield( A.TITLE_WEIGHT )
  set( get(h,'Title'), 'FontWeight', A.TITLE_WEIGHT );
end



%= X axis
%
if isfield( A, 'XLABEL' ) & ~isfield( A.XLABEL )
  set( get(h,'Xlabel'), 'String', A.XLABEL );
end
%
if isfield( A, 'XLABEL_SIZE' ) & ~isfield( A.XLABEL_SIZE )
  set( get(h,'Xlabel'), 'FontSize', A.XLABEL_SIZE );
end
%
if isfield( A, 'XLABEL_WEIGHT' ) & ~isfield( A.XLABEL_WEIGHT )
  set( get(h,'Xlabel'), 'FontWeight', A.XLABEL_WEIGHT );
end
%
if isfield( A, 'XGRID' ) & ~isfield( A.XGRID )
  set( h, 'XGrid', A.XGRID );
end
%
if isfield( A, 'XTICKLABEL' ) & ~isfield( A.XTICKLABEL )
  set( h, 'XTickLabel', A.XTICKLABEL );
end



%= Y axis
%
if isfield( A, 'YLABEL' ) & ~isfield( A.YLABEL )
  set( get(h,'Ylabel'), 'String', A.YLABEL );
end
%
if isfield( A, 'YLABEL_SIZE' ) & ~isfield( A.YLABEL_SIZE )
  set( get(h,'Ylabel'), 'FontSize', A.YLABEL_SIZE );
end
%
if isfield( A, 'YLABEL_WEIGHT' ) & ~isfield( A.YLABEL_WEIGHT )
  set( get(h,'Ylabel'), 'FontWeight', A.YLABEL_WEIGHT );
end
%
if isfield( A, 'YGRID' ) & ~isfield( A.YGRID )
  set( h, 'YGrid', A.YGRID );
end
%
if isfield( A, 'YTICKLABEL' ) & ~isfield( A.YTICKLABEL )
  set( h, 'YTickLabel', A.YTICKLABEL );
end

