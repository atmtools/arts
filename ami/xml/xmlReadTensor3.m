%------------------------------------------------------------------------
%
% NAME:    xmlReadTensor3
%
%          Reads Tensor3 data from XML file.
%
% FORMAT:  result = xmlReadTensor3(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-10-18  Created by Oliver Lemke.

function result = xmlReadTensor3(fid, attrlist)

np = str2num (xmlGetAttrValue (attrlist, 'npages'));
nr = str2num (xmlGetAttrValue (attrlist, 'nrows'));
nc = str2num (xmlGetAttrValue (attrlist, 'ncols'));
nelem =  np * nr * nc;

result = fscanf (fid, '%f', nelem);
xmlCheckSize (nelem, size (result));

result = permute (reshape (result, [nc nr np]), [2 1 3]);
