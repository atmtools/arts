%------------------------------------------------------------------------
%
% NAME:    xmlReadTensor4
%
%          Reads Tensor4 data from XML file.
%
% FORMAT:  result = xmlReadTensor4(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-10-18  Created by Oliver Lemke.

function result = xmlReadTensor4(fid, attrlist)

nb = str2num (xmlGetAttrValue (attrlist, 'nbooks'));
np = str2num (xmlGetAttrValue (attrlist, 'npages'));
nr = str2num (xmlGetAttrValue (attrlist, 'nrows'));
nc = str2num (xmlGetAttrValue (attrlist, 'ncols'));
nelem =  nb * np * nr * nc;

result = fscanf (fid, '%f', nelem);
xmlCheckSize (nelem, size (result));

result = permute (reshape (result, [nc nr np nb]), [4 3 2 1]);

