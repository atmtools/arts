%------------------------------------------------------------------------
% NAME:     isscalar
%
%           Determines if a variable is a scalar, that is, a numeric
%           of size 1x1.
%
% FORMAT:   bool = isscalar( a )
%
% RETURN:   bool     1 if a scalar, else 0.
% IN:       a        A variable.
%------------------------------------------------------------------------

% HISTORY: 2002.01.03  Created by Patrick Eriksson. 


function bool = isscalar( a )


bool = 0;


if isnumeric( a )
  
  if max(size(a)) == 1

    bool = 1;

  end

end

