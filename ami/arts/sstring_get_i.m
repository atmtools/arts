%------------------------------------------------------------------------
% NAME:   sstring_get_i
%
%         Returns one string of a "string-string".
%
%         See sstring_length for format of a string-string.
%
% FORMAT: s = sstring_get_i(sstring,i)
%
% OUT:    s         String nr i.
% IN:     sstring   A string-string.
%         i         Index in SSTRING.
%------------------------------------------------------------------------

% HISTORY: 2000.12.22  Created by Patrick Eriksson.


function s = sstring_get_i(sstring,i)

if i > sstring_length(sstring)
  error(sprintf('The given string-string (%s)contains not %d strings.',sstring,i));
end

ind = find( sstring == '"' );

s = sstring( [(ind((i-1)*2+1)+1):(ind((i-1)*2+2)-1)] );
