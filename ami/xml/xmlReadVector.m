%------------------------------------------------------------------------
%
% NAME:    xmlReadVector
%
%          Reads Vector data from XML file.
%
% FORMAT:  result = xmlReadVector(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-11-14  Created by Oliver Lemke.

function result = xmlReadVector(fid, attrlist)

nr = str2num (xmlGetAttrValue (attrlist, 'nrows'));

result = fscanf (fid, '%f', nr);
xmlCheckSize (nr, size (result));

