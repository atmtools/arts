%------------------------------------------------------------------------
% NAME:    vec2row
%
%          Ensures that a matrix not has less columnss than rows.
%          The most common application of this function is to
%          ensure that a vector is a row vector.
%
% FORMAT:  v = vec2row(v)
%
% RETURN:  v           row vector
% IN:      v           row or column vector
%------------------------------------------------------------------------

% HISTORY: 1993  Created by Patrick Eriksson. 


function v = vec2row(v)

[rows,cols] = size(v);

if (rows > cols)
 v = v';
end
