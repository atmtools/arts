%------------------------------------------------------------------------
%
%
% NAME:    xmlReadAttributes
%
%          Reads tag attributes from xml file.
%
% FORMAT:  attrlist = xmlReadAttributes(fid)
%
% RETURN:  attrlist    Attribute list
% IN:      fid         File descriptor
%------------------------------------------------------------------------
%
%

% HISTORY: 2002-09-25  Created by Oliver Lemke.

function attrlist = xmlReadAttributes(fid)

na = 0;
s = ' ';
while s(size(s)) ~= '>'
  s = fscanf (fid, '%s', 1);
  attr = strtok (s, '=');
  pos = strfind (s, '"');
  value = s (pos(1)+1:pos(2)-1);

  na = na + 1;
  attrlist{na,1} = attr;
  attrlist{na,2} = value;
end

