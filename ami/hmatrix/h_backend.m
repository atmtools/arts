%------------------------------------------------------------------------
% NAME:    h_backend
%
%          Returns the H matrix for a backend (spectrometer) where all
%          channels have identical response.
%          The frequency grid valid for the backend is also returned.
%          The response of the channels is normalised and the
%          response values (W_BACK) do not need to be normalised.
%
% RETURN:  H           H for the backend
%          f_new       new frequency grid (=F2)
% IN:      f1          monochromatic frequencies
%          f2          centre frequencies for the channels
%          za          viewing angles
%          f_back      grid points for the channel response
%          w_back      the channel response
%------------------------------------------------------------------------

% HISTORY: 14.11.00  Created for Skuld by Patrick Eriksson. 
%          25.08.00  Adapted to AMI by Patrick Eriksson


function [H,f_new] = h_backend(f1,f2,za,f_back,w_back)


%=== Main sizes
nf1   = length(f1);
nf2   = length(f2);
nza   = length(za);
nback = length(f_back);


%=== Allocate H and create F_NEW
H      = sparse(nza*nf2,nza*nf1);
f_new  = vec2col(f2);


%=== Check if backend channels totally inside frequency grid
if f1(1) > f2(1)+f_back(1)
  error(sprintf('You must increase your frequency grid downwards with %.3e Hz',f1(1)-(f2(1)+f_back(1))));
end
if f1(nf1) < f2(nf2)+f_back(nback)
  error(sprintf('You must increase your frequency grid upwards with %.3e Hz',(f2(nf2)+f_back(nback))-f1(nf1)));
end


%=== Fill H
%
ind1      = 1:nf1;
%
for i = 1:nf2
     
  f   = f1 - f2(i);

  w    = h_weights_integr(f,f_back,w_back);
  w    = w/sum(w);
  ind3 = find(w~=0);

  for j = 1:nza
    ind2 = ind1 + (j-1)*nf1;
    H(i+(j-1)*nf2,ind2(ind3)) = w(ind3);
  end
end

