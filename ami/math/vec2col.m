%------------------------------------------------------------------------
% NAME:    vec2col
%
%          Ensures that a matrix not has less rows than columns.
%          The most common application of this function is to
%          ensure that a vector is a column vector.
%
% FORMAT:  v = vec2col(v)
%
% RETURN:  v           column vector
% IN:      v           row or column vector
%------------------------------------------------------------------------

% HISTORY: 1993  Created by Patrick Eriksson. 


function v = vec2col(v)

[rows,cols] = size(v);

if (cols > rows)
 v = v';
end
