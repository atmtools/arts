%------------------------------------------------------------------------
% NAME:   sstring_length 
%
%         Returns the length of a "string-string". A string-string is a 
%         Matlab string holding a ARTS string array. With other words,
%         a string-string has the format:
%
%            '"string1","string2"'
%
%         The function gives an error messages if the input:
%            not is a Matlab string 
%            0 or an odd number of " are found 
%
% FORMAT: n = sstring_length(sstring)
%
% OUT:    n         Number of string-strings.
% IN:     sstring   A string-string.
%------------------------------------------------------------------------

% HISTORY: 2000.12.22  Created by Patrick Eriksson.


function n = sstring_length(sstring)

if ~ischar(sstring)
  error('A non-string variable passed to sstring_length.');
end

ind = find( sstring == '"' );

n = length(ind);

if (n==0) | isodd(n)
  error(sprintf('Not a valid sstring (%s)',sstring));
end

n = round(n/2);

