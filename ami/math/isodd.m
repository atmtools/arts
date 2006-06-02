%------------------------------------------------------------------------
% NAME:    isodd
%
%          Determines if a number is odd.
%
% FORMAT:  bool = isodd(n)
%
% RETURN:  bool        1 if n odd, else 0
% IN:      n           a number
%------------------------------------------------------------------------

% HISTORY: 00.08.16  Created by Patrick Eriksson. 


function bool = isodd(n)

if ~iswhole(n)
  error('The given value is not an integer.');
end

bool = rem(n,2);
