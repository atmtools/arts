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
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%------------------------------------------------------------------------

% HISTORY: 99.11.18  Created for Skuld by Patrick Eriksson. 
%          00.08.25  Adapted to AMI by Patrick Eriksson
%          00.11.16  Included linear/cubic flag (PE) 


function [H,f_new] = h_mixer(f,za,lo,fprimary,f_filter,w_filter,o_filter,o_y)


%=== Main sizes
nf    = length(f);
nza   = length(za);
nfilt = size(f_filter,1);


%=== Check if orders are 1 or 3
if ~( (o_filter==1) | (o_filter==3) ) 
  error('The polynomial order must be 1 or 3.');
end
if ~( (o_y==1) | (o_y==3) ) 
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
    error('There must be at least frequencies in the image band');
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


out(1,1);
out(1,'Setting up mixer and sidenband filter transfer matrix (H).');


%=== Some help variables
if upper
  iprim  = ilo;
  iimage = ilo-1;
  istop  = nf;
else
  iprim  = 1;
  iimage = nf;
  istop  = ilo-1;
end
%
di = ceil( o_y/2 );
%
nout = istop - iprim + 1;


%=== Allocate vectors for row and col index and vector for weights
lrow = (o_y+2) * nout * nza;% This variable is used as length for these vectors
nrow = 0;                   % Number of values in the vectors
rows = zeros( 1, lrow );
cols = zeros( 1, lrow );
wgts = zeros( 1, lrow );



%=== Fill H
%
for i = iprim:istop

  out(2,sprintf('Doing frequency %d of %d',i-iprim+1,nout));

  %=== Some variables
  fimage = 2*lo - f(i);
  iout   = i - iprim + 1;

  %=== Get sideband filter values
  if ( o_filter == 1 )
    w = interp1(f_filter,w_filter,[f(i),fimage],'linear');
  else
    w = interp1(f_filter,w_filter,[f(i),fimage],'cubic');
  end

  %=== Determine index of frequency point just above image frequency
  while f(iimage-1) > fimage
    iimage = iimage - 1;
  end

  if upper
    indb = [ max([1,iimage-di]) : min([iimage+di-1,ilo-1]) ];
  else
    indb = [ max([ilo,iimage-di]) : min([iimage+di-1,nf]) ];
  end
  nb = length(indb);
  if nb == 2
    wimage = pbasis( f(indb) ) * [fimage;1];
  elseif nb == 3
    wimage = pbasis( f(indb) ) * [fimage^2;fimage;1];
  else
    wimage = pbasis( f(indb) ) * [fimage^3;fimage^2;fimage;1];
  end
  wimage = wimage';

  for j = 1:nza

    %=== Primary frequency
    irow       = nrow + 1;
    rows(irow) = iout+(j-1)*nout;
    cols(irow) = i+(j-1)*nf;
    wgts(irow) = w(1)/(w(1)+w(2));
    nrow       = nrow + 1;

    %=== Image frequency
    irow       = nrow + (1:nb);
    rows(irow) = iout+(j-1)*nout;
    cols(irow) = indb+(j-1)*nf;
    wgts(irow) = wimage * w(2)/(w(1)+w(2));
    nrow       = nrow + nb;

  end
end


%=== Allocate H and create F_NEW
if upper
  H = sparse( rows(1:nrow), cols(1:nrow), wgts(1:nrow), (nf-ilo+1)*nza,nf*nza);
  f_new = f(ilo:nf);
else
  H = sparse( rows(1:nrow), cols(1:nrow), wgts(1:nrow), (ilo-1)*nza, nf*nza );
  f_new = f(1:(ilo-1));
end


out(1,-1);
