%------------------------------------------------------------------------
% NAME:   tag2mo
%
%         Return the molecule name for a given tag string.
%
%         For example: tag2mo('O3-666') gives 'O3'.
%
% FORMAT: mo = tag2mo( tag )
%
% OUT:    mo    String with molecule name.
% IN:     tag   Tag string, e.g. H2O-161.
%------------------------------------------------------------------------

% HISTORY: 2002.01.08  Created by Patrick Eriksson.


function mo = tag2mo( tag )


i = find( tag == '-' );

if isempty( i )

  mo = tag;

else

  mo = tag(1:(i(1)-1));

end