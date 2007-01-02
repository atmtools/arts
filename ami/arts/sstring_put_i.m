%------------------------------------------------------------------------
% NAME:   sstring_put_i
%
%         Replaces one string entry of a "string-string".
%
%         See sstring_length for format of a string-string.
%
% FORMAT: sstring_out=sstring_put_i(sstring_in,i,s)
%
% OUT:    sstring_out   modified string-string.
% IN:     sstring_in    a string-string.
%         i             Index in sstring_in.
%         s             string to replace string i.
%------------------------------------------------------------------------

% HISTORY: 2006-08-09  Created by Jo Urban.


function sstring_out = sstring_put_i(sstring_in,i,s)

if i > sstring_length(sstring_in)
  error(sprintf('The given string-string (%s) contains not %d strings.',sstring_in,i));
end

ind = find( sstring_in == '"' );

start_pos = (ind((i-1)*2+1)+1);
stop_pos  = (ind((i-1)*2+2)-1);

s_start = sstring_in( [1:start_pos-1] );
s_stop  = sstring_in( [stop_pos+1:end] );

sstring_out = [s_start,s,s_stop];

return
