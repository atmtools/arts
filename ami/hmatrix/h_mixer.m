%------------------------------------------------------------------------
% NAME:    h_mixer
%
%          Returns the H matrix for sideband filtering and mixer.
%
%          The frequency grid valid after the mixer is also returned.
%          The output frequencies are observed frequencies (not IF).
%
%          The new frequency grid includes:
%             1. all frequencies of the primary band
%             2. the frequencies of the image band which projection
%                is inside the limits of the primary band
%
%          The image frequencies must cover the projection of all 
%          frequencies of the primary band. 
%
%          The response of the sideband filter is normalised inside the
%          function and the input filter values (W_FILTER) do not need 
%          to be normalised.
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
%          02.01.05  Major revision to include image frequencies in f_new


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


%=== Create new frequency grid
%
if upper
  fprim = f(ilo:nf);
else
  fprim = f(1:(ilo-1));
end
%
fimband = 2*lo - [ last(fprim) fprim(1) ];
%
fimage  = f( find( f>=fimband(1) & f<=fimband(2) ) );
%
f_new = unique( [ fprim; 2*lo-fimage ] );
%
clear fprim fimage fimband
%
n_new = length( f_new );


%=== Allocate vectors for row and col index and vector for weights
%
lrow = 2 * (o_y+1) * n_new * nza;
%                   
rows = zeros( 1, lrow );
cols = zeros( 1, lrow );
wgts = zeros( 1, lrow );
%
nrow = 0;


%=== Half of the polynomial order
%
di = ceil( o_y/2 );


%=== Loop frequencies in F_NEW
%
for i = 1 : n_new


  out(2,sprintf('Doing frequency %d of %d',i,n_new));


  %= Primary and image frequency
  %
  fprim  = f_new(i);
  fimage = 2*lo - fprim;


  %= Determine index in F of frequency just below F_NEW(I) and its image
  % 
  iprim  = max( find( f <= fprim ) ); 
  iimage = max( find( f <= fimage ) ); 


  %= Index for involved frequency points
  %
  if fprim == f(iprim)
    ind_prim = iprim;
  else
    if upper
      ind_prim  = [ max([iprim-di+1,ilo]) : min([iprim+di,nf]) ];
    else
      ind_prim  = [ max([iprim-di+1,1]) : min([iprim+di,ilo-1]) ];
    end
  end
  %
  if fimage == f(iimage)
    ind_image = iimage;
  else
    if upper
      ind_image = [ max([iimage-di+1,1]) : min([iimage+di,ilo-1]) ];
    else
      ind_image = [ max([iimage-di+1,ilo]) : min([iimage+di,nf]) ];
    end
  end


  %=== Get sideband filter values
  %
  if ( o_filter == 1 )
    w = interp1(f_filter,w_filter,[fprim,fimage],'linear');
  else
    w = interp1(f_filter,w_filter,[fprim,fimage],'cubic');
  end


  %=== Weights for primary band
  %
  np = length( ind_prim );
  %
  switch np
    case 1 
      wprim = 1;
    case 2
      wprim = pbasis( f(ind_prim) ) * [fprim;1];
    case 3
      wprim = pbasis( f(ind_prim) ) * [fprim^2;fprim;1];
    case 4
      wprim = pbasis( f(ind_prim) ) * [fprim^3;fprim^2;fprim;1];
    otherwise
      error('Index error for primary band (problem for PE!).');
  end
  wprim = wprim';


  %=== Weights for image band
  %
  ni = length( ind_image );
  %
  switch ni
    case 1 
      wimage = 1;
    case 2
      wimage = pbasis( f(ind_image) ) * [fimage;1];
    case 3
      wimage = pbasis( f(ind_image) ) * [fimage^2;fimage;1];
    case 4
      wimage = pbasis( f(ind_image) ) * [fimage^3;fimage^2;fimage;1];
    otherwise
      error('Index error for image band (problem for PE!).');
  end
  wimage = wimage';


  %= Fill rows, cols, wghts ...
  %
  for j = 1:nza

    %=== Primary frequency
    irow       = nrow + (1:np);
    rows(irow) = i + (j-1)*n_new;
    cols(irow) = ind_prim + (j-1)*nf;
    wgts(irow) = wprim * w(1)/(w(1)+w(2));
    nrow       = nrow + np;

    %=== Image frequency
    irow       = nrow + (1:ni);
    rows(irow) = i + (j-1)*n_new;
    cols(irow) = ind_image + (j-1)*nf;
    wgts(irow) = wimage * w(2)/(w(1)+w(2));
    nrow       = nrow + ni;

  end

end


%=== Allocate H
%
H = sparse( rows(1:nrow), cols(1:nrow), wgts(1:nrow), n_new*nza, nf*nza );


out(1,-1);


return

%===
%=== Old code 
%===


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
