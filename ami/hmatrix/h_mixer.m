%------------------------------------------------------------------------
% NAME:    h_mixer
%
%          Returns the H matrix for sideband filtering and mixer.
%          The frequency grid valid for the mixer is also returned, but
%          the original frequency values are returned, and not IF.
%          The new frequency grid consists of all frequncies in F
%          inside the primary band.
%          The response of the sideband filter is normalised and the
%          filter values (W_FILTER) do not need to be normalised.
%
% FORMAT:  [H,f_new] = h_mixer(f,za,lo,fprimary,f_filter,w_filter,o_filter)
%
% RETURN:  H           H matrix for the mixer
%          f_new       new frequency grid
% IN:      f           frequencies
%          za          zenith angles
%          lo          LO frequency
%          fprimary    a frequency inside the primary band (!=LO)
%          f_filter    frequencies for W_FILTER
%          w_filter    sideband filter values
%          o_filter    linear (=1) or cubic (=3) treatment of the
%                      sideband filter
%------------------------------------------------------------------------

% HISTORY: 99.11.18  Created for Skuld by Patrick Eriksson. 
%          00.08.25  Adapted to AMI by Patrick Eriksson
%          00.11.16  Included linear/cubic flag (PE) 


function [H,f_new] = h_mixer(f,za,lo,fprimary,f_filter,w_filter,o_filter)


%=== Main sizes
nf    = length(f);
nza   = length(za);
nfilt = size(f_filter,1);


%=== Check if orders are 1 or 3
if ~( (o_filter==1) | (o_filter==3) ) 
  error('The polynomial order must be 1 or 3.');
end


%=== Determine if upper or lower side is primary band
if fprimary > lo
  upper = 1;
elseif fprimary < lo
  upper = 0;
else
  error('The primary band frequency must deviate from the LO frequency');
end


%=== Determine first index above LO
ilo  = min(find(f>lo));
if isempty(ilo) | (ilo==1)
  error('You must have frequencies on both side of the LO frequency');
end


%=== Check that image frequncies cover primary freqiencies
if upper
  if ilo <= 2
    error('There must be at least 2 frequencies in the image band')
  end
  if 2*lo-f(nf) < f(1)
    error(sprintf('You must increase the frequency grid at the bottom end with %.3e Hz',f(1)-(2*lo-f(nf))));
  end
  if 2*lo-f(ilo) > f(ilo-1)
    error(sprintf('You must include frequencies in the lower band closer to the LO, %.3e Hz is missing',(2*lo-f(ilo))-f(ilo-1)));
  end
else
  if ilo == nf
    error('There must be at least 2 frequencies in the image band')
  end
  if 2*lo-f(1) > f(nf)
    error(sprintf('You must increase the frequency grid at the top end with %.3e Hz',(2*lo-f(1))-f(nf)));
  end
  if 2*lo-f(ilo-1) < f(ilo)
    error(sprintf('You must include frequencies in the upper band closer to the LO, %.3e Hz is missing',f(ilo)-(2*lo-f(ilo-1))));
  end
end


%=== Allocate H and create F_NEW
if upper
  H = sparse((nf-ilo+1)*nza,nf*nza);
  f_new = f(ilo:nf);
else
  H = sparse((ilo-1)*nza,nf*nza);
  f_new = f(1:(ilo-1));
end


%=== Fill H
if upper
  iprim  = ilo;
  iimage = ilo-1;
  istep  = -1;
  istop  = nf;
else
  iprim  = 1;
  iimage = ilo;
  istep  = 1;
  istop  = ilo-1;
end
%
nout = istop - iprim + 1;
%
for i = iprim:istop

  %=== Some variables
  fimage = 2*lo - f(i);
  iout   = i - iprim + 1;

  %=== Get sideband filter values
  if ( o_filter == 1 )
    w = interp1(f_filter,w_filter,[f(i),fimage],'linear');
  else
    w = interp1(f_filter,w_filter,[f(i),fimage],'cubic');
  end

  %=== Determine index/indeces of image frequency
  while istep*f(iimage+istep) < istep*fimage
    iimage = iimage + istep;
  end
  
  for j = 1:nza
    %=== Primary frequency
    H(iout+(j-1)*nout,i+(j-1)*nf) = w(1)/(w(1)+w(2));

    %=== Image frequency
    H(iout+(j-1)*nout,iimage+(j-1)*nf) = w(2)/(w(1)+w(2)) * ...
                  abs((fimage-f(iimage+istep))/(f(iimage)-f(iimage+istep)));
    H(iout+(j-1)*nout,iimage+istep+(j-1)*nf) = w(2)/(w(1)+w(2)) * ...
                  abs((fimage-f(iimage))/(f(iimage)-f(iimage+istep)));
  end
end
