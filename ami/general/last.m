%------------------------------------------------------------------------
% NAME:    last
%
%          Returns the last element of a vector.
%
% FORMAT:  v = last(x)
%
% RETURN:  v    Last element of x ( v = x(length(x)); ).
% IN:      x    A vector.
%------------------------------------------------------------------------

% HISTORY: 2000.12.18  Created by Patrick Eriksson. 


function v = last(x)

v = x(length(x));
