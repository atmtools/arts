%------------------------------------------------------------------------
%
% NAME:    xmlReadGasAbsLookup
%
%          Reads GasAbsLookup table from XML file.
%
% FORMAT:  result = xmlReadGasAbsLookup(fid, attrlist)
%
% RETURN:  result      Data read from file
% IN:      fid         File descriptor
% IN:      attrlist    List of attributes
%------------------------------------------------------------------------
%

% HISTORY: 2002-11-28  Created by Oliver Lemke.

function result = xmlReadGasAbsLookup(fid, attrlist)

  result.species  = xmlReadTag(fid);
  result.f_grid   = xmlReadTag(fid);
  result.p_grid   = xmlReadTag(fid);
  result.vmrs_ref = xmlReadTag(fid);
  result.t_ref    = xmlReadTag(fid);
  result.t_pert   = xmlReadTag(fid);
  result.nls_pert = xmlReadTag(fid);
  result.abs      = xmlReadTag(fid);

