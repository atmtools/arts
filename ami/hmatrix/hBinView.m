%------------------------------------------------------------------------
% NAME:    hBinView
%
%          Includes binning of viewing angles into H and Hd. That is, 
%          averages of subsequent zenith angles are used.
%          This function assumes (and checks) that all zenith angles
%          have the same frequencies. Accordingly, this function should
%          normally be used before any other data reduction.
%          The binning is specified by the vector BINS. The elements of BINS
%          give the number of subsequent zenith angles to combine. For example,
%            BINS = [1 2 1];
%          means that angles 2 and 3 shall be combined while angles 1 and 4
%          shall be untouched.
%
% FORMAT:  [H,f_y,za_y,Hd] = hBinView(H,f_y,za_y,Hd,bins)
%
% RETURN:  H           total H matrix after angle binning
%          f_y         new frequency vector
%          za_y        new zenith angle vector 
%          Hd          data reduction H matrix after angle binning
% IN:      H           total H matrix before binning
%          f_y         input frequency vector 
%          za_y        input zenith angle vector
%          Hd          data reduction H matrix before angle binning
%          bins        binning pattern (see above)
%------------------------------------------------------------------------

% HISTORY: 25.08.00  Created by Patrick Eriksson. 


function [H,f_y,za_y,Hd] = hBinView(H,f_y,za_y,Hd,bins)


%=== Pick out frequency grid
ny  = length(f_y);     % total length
nf  = 0;               % length of frequency grid
while (nf<ny) & (za_y(nf+1)==za_y(1))
  nf = nf + 1;
end
f   = f_y(1:nf);


%=== Pick out zenith angle grid and perform checks
if ny ~= length(za_y)
  error('F_Y and ZA_Y must have the same length');
end
if rem(ny,nf) ~= 0
  error('The length of F_Y and ZA_Y does not allow a single frequency grid for all zenith angles');
end
nza = ny/nf;
za  = zeros(nza,1);
for i = 1:nza
  za(i) = za_y((i-1)*nf+1);
  for j = 1:f
    ind = (i-1)*nf+j;
    if (f_y(ind)~=f(j)) | (za_y(ind)~=za(i)) 
      error('Incorrect F_Y and/or ZA_Y, all zenith angles must have the same frequencies');
    end
  end
end


%=== Check binning pattern
nout  = length(bins);
if any(bins < 1)
  error('All the binning values must be >= 1');
end
if sum(bins) ~= nza
  error(sprintf('The binning pattern does not match the data (number of ZA: %d, sum of BINS: %d)',nza,sum(bins)));
end


%=== Check if there is something to do
if nout == nza
  return
end


%=== Set up ZA_NEW and some vectors defining the local Hd
za_new = zeros(nout,1);         % new zenith angles
row    = ones(nf*nza,1);        % row index in Hd
col    = (1:nf*nza)';           % col index in Hd
w      = ones(nf*nza,1);        % the values in Hd
i1     = 0;                     % I1 and I2, indeces for ZA
for i = 1:nout
  i1 = i1 + 1;
  i2 = i1 + bins(i) - 1;  
  wv = 1/bins(i);

  za_new(i) = mean(za(i1:i2));

  for j = i1:i2
    ind      = (1:nf) + (j-1)*nf;
    row(ind) = (1:nf) + (i-1)*nf;
    w(ind)   = wv;
  end
  i1 = i2;
end


%=== Include the viewing angle binning in H and Hd
Hb = sparse(row,col,w);
H  = h_x_h(Hb,H);
Hd = h_x_h(Hb,Hd);


%=== Create new F_Y and ZA_Y
[f_y,za_y] = h_fix_ys(f,za_new);
