%---------------------------------------------------------------
% h = legend2(tsize,h,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
%
% This functions does the same as the normal legend function but
% writes the text with size defined by tsize
%---------------------------------------------------------------


function h = legend2(tsize,h,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)

ha      = gca;
asize   = get(ha,'FontSize');
aweight = get(ha,'FontWeight');

set(ha,'FontSize',tsize,'FontWeight','normal');

n	= nargin-2;
if n == 1	%Fix this better if you can !!
  h = legend(h,s1,0);
elseif n == 2
  h = legend(h,s1,s2,0);
elseif n == 3
  h = legend(h,s1,s2,s3,0);
elseif n == 4
  h = legend(h,s1,s2,s3,s4,0);
elseif n == 5
  h = legend(h,s1,s2,s3,s4,s5,0);
elseif n == 6
  h = legend(h,s1,s2,s3,s4,s5,s6,0);
elseif n == 7
  h = legend(h,s1,s2,s3,s4,s5,s6,s7,0);
elseif n == 8
  h = legend(h,s1,s2,s3,s4,s5,s6,s7,s8,0);
elseif n == 9
  h = legend(h,s1,s2,s3,s4,s5,s6,s7,s8,s9,0);
elseif n == 10
  h = legend(h,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,0);
end

set(ha,'FontSize',asize);
set(ha,'FontWeight',aweight);
