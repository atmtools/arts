%------------------------------------------------------------------------
%
% NAME:    xmlReadArray
%
%          Reads Array from XML file.
%
% FORMAT:  result = xmlReadArray(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-09-25  Created by Oliver Lemke.

function result = xmlReadArray(fid, attrlist)

ne = str2num (xmlGetAttrValue (attrlist, 'nelem'));

e = 0;
while e ~= ne
  e = e + 1;
  result{e} = xmlReadTag(fid);
end
