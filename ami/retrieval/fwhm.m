%==============================================================================
% function res = fwhm(z,A)
%
% Calculates the full width at half mean of the averaging kernels.
%
% In:	z	vector of heights [m]
%	A	matrix of averaging kernels
%
% Out:	res	FWHM vertical resolution [m]
%==============================================================================

% Patrick Eriksson 1995


function res = fwhm(z,A)


[m,n] = size(A);
res = zeros(m,1);


for i = 1:m
  [mv,i0] = max(A(i,:));

  if mv == 0
    res(i) = NaN;
  else

  i2 = i0-1;
  if (i2 == 0), i2 = 1; end
  while A(i,i2)>(mv/2) 
    i2 = i2-1;
    if (i2 < 1),break; end
  end
  if (i2 >= 1)
    p1 = interp1(A(i,[i2:i2+1]),z([i2:i2+1]),mv/2);
  else
    p1 = z(1);
  end
  i2 = i0+1;
  if (i2 > n), i2 = n; end
  while A(i,i2)>(mv/2)
    i2 = i2+1;
    if (i2 > n), break; end
  end
  if (i2 <= n)
    p2 = interp1(A(i,[i2-1:i2]),z([i2-1:i2]),mv/2);
  else
    p2 = z(n);
  end
  res(i) = p2-p1;

  end
end
