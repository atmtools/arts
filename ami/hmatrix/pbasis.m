%------------------------------------------------------------------------
% NAME:    pbasis
%
%          Gives a polynomial basis using the Lagrange formula.
%
%          Returns the coefficients of a polynomial sum to express a 
%          function with known abscissa, but where the functional
%          values are not yet known.
%          If we have the function y1 = f(x1), y2 = f(x2), ..., yN = f(xN),
%          (i.e. the length of X is N) the returned polynomial coefficients 
%          approximating f(x) are
%
%          y = (b(1,1)*x^n+b(1,2)*x^(n-1)+...+b(1,n))*y(1) +
%              (b(2,1)*x^n+b(2,2)*x^(n-1)+...+b(2,n))*y(2) +
%              ...
%              (b(n,1)*x^n+b(n,2)*x^(n-1)+...+b(n,n))*y(n) +
%
%          The size of B is accordingly N x N and the column order matches
%          the order of POLYFIT.
%            
%          Valid values of N are 1-4.
%
%          The Lagrange formula is described in Numeric Recipies (Sec. 3.1).
%
% FORMAT:  b = pbasis(x)
%
% RETURN:  b           Polynomial coefficients.
% IN:      x           Function abscissa.
%------------------------------------------------------------------------

% HISTORY: 00.11.15  Created by Patrick Eriksson (PE).

% If someone can find a general solution (i.e. valid for any N),
% please go ahead and implemented it.

% This function does not use the full potential of Matlab for consistency
% reasons when it will be ported to ARTS.


function b = pbasis(x)


n = length(x);
b = zeros(n,n);


switch n

  case 1
    b = 1;

  case 2
    a      = x(1) - x(2);
    b(1,1) = 1/a;
    b(1,2) = -x(2)/a;
    a      = x(2) - x(1);
    b(2,1) = 1/a;
    b(2,2) = -x(1)/a;

  case 3
    a      = (x(1)-x(2))*(x(1)-x(3));
    b(1,1) = 1/a;
    b(1,2) = (-x(2)-x(3))/a;
    b(1,3) = x(2)*x(3)/a;
    a      = (x(2)-x(1))*(x(2)-x(3));
    b(2,1) = 1/a;
    b(2,2) = (-x(1)-x(3))/a;
    b(2,3) = x(1)*x(3)/a;
    a      = (x(3)-x(1))*(x(3)-x(2));
    b(3,1) = 1/a;
    b(3,2) = (-x(1)-x(2))/a;
    b(3,3) = x(1)*x(2)/a;

  case 4
    a      = (x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4));
    b(1,1) = 1/a;
    b(1,2) = (-x(2)-x(3)-x(4))/a;
    b(1,3) = (x(2)*x(3)+x(2)*x(4)+x(3)*x(4))/a;
    b(1,4) = -x(2)*x(3)*x(4)/a;
    a      = (x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4));
    b(2,1) = 1/a;
    b(2,2) = (-x(1)-x(3)-x(4))/a;
    b(2,3) = (x(1)*x(3)+x(1)*x(4)+x(3)*x(4))/a;
    b(2,4) = -x(1)*x(3)*x(4)/a;
    a      = (x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4));
    b(3,1) = 1/a;
    b(3,2) = (-x(1)-x(2)-x(4))/a;
    b(3,3) = (x(1)*x(2)+x(1)*x(4)+x(2)*x(4))/a;
    b(3,4) = -x(1)*x(2)*x(4)/a;
    a      = (x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3));
    b(4,1) = 1/a;
    b(4,2) = (-x(1)-x(2)-x(3))/a;
    b(4,3) = (x(1)*x(2)+x(1)*x(3)+x(2)*x(3))/a;
    b(4,4) = -x(1)*x(2)*x(3)/a;

  otherwise
    error('The length of the input vector must be 1-4.')
end
