%------------------------------------------------------------------------
%
% NAME:    xmlReadTensor5
%
%          Reads Tensor5 data from XML file.
%
% FORMAT:  result = xmlReadTensor5(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-10-18  Created by Oliver Lemke.

function result = xmlReadTensor5(fid, attrlist)

nv = str2num (xmlGetAttrValue (attrlist, 'nvitrines'));
nb = str2num (xmlGetAttrValue (attrlist, 'nbooks'));
np = str2num (xmlGetAttrValue (attrlist, 'npages'));
nr = str2num (xmlGetAttrValue (attrlist, 'nrows'));
nc = str2num (xmlGetAttrValue (attrlist, 'ncols'));
nelem =  nv * nb * np * nr * nc;

result = fscanf (fid, '%f', nelem);
xmlCheckSize (nelem, size (result));

result = permute (reshape (result, [nc nr np nb nv]), [2 1 3 4 5]);
