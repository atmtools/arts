%------------------------------------------------------------------------
%
% NAME:    xmlCheckSize
%
%          Checks whether enough data field were read from file.
%          Exits with error message if not.
%
% FORMAT:  xmlCheckSize(ne, nr)
%
% IN:      ne         Expected number of items
% IN:      nr         Read number of items
%------------------------------------------------------------------------
%

% HISTORY: 2002-10-18 Created by Oliver Lemke.

function xmlCheckSize(ne, nr)

if nr(1) ~= ne
  error ('Not enough input data found in file');
end

