%------------------------------------------------------------------------
%
% NAME:    xmlReadMatrix
%
%          Reads Matrix data from XML file.
%
% FORMAT:  result = xmlReadMatrix(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-09-25  Created by Oliver Lemke.

function result = xmlReadMatrix(fid, attrlist)

nr = str2num (xmlGetAttrValue (attrlist, 'nrows'));
nc = str2num (xmlGetAttrValue (attrlist, 'ncols'));

result = fscanf (fid, '%f', nr * nc);
result = reshape (result, [nr nc]);
