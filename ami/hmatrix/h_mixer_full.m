%------------------------------------------------------------------------
% NAME:    h_mixer_full
%
%          Returns the H matrix for sideband filtering and mixer.
%
%          The frequency grid valid after the mixer is also returned.
%          The output frequencies are observed frequencies (not IF).
%
%          The new frequency grid is composed of both primary and
%          mirrored image frequencies. Only frequencies that are
%          surrounded by grid points, i.e. that can be interpolated,
%          are stored.
%
%          The response of the sideband filter is normalised inside the
%          function and the input filter values (W_FILTER) do not need
%          to be normalised.
%
%          The last two inputs, o_f and o_y, must be set to 1. They are
%          still here to make a transition between h_matrix och
%          h_matrix_full easier.
%
% FORMAT:  [H,f_new] = h_mixer_full(f,za,lo,fprimary,f_filter,w_filter,o_f)
%
% RETURN:  H           H matrix for the mixer
%          f_new       new frequency grid
% IN:      f           frequencies
%          za          zenith angles
%          lo          LO frequency
%          fprimary    a frequency inside the primary band (!=LO)
%          f_filter    frequencies for W_FILTER
%          w_filter    sideband filter values
%          o_f         linear (=1)
%          o_y         linear (=1)
%------------------------------------------------------------------------

% HISTORY: 04.03.04  Adapted from h_matrix and altered by Mattias Ekström.


function [H,f_new] = h_mixer_full(f,za,lo,fprimary,f_filter,w_filter,o_f,o_y)


if( o_f ~= 1  |  o_y ~= 1 )
  error('Only linear treatment of the variables is now allowed.');
end


%=== Main sizes
nf    = length(f);
nza   = length(za);
nfilt = size(f_filter,1);


%=== Determine if upper or lower side is primary band and create frequency
%=== grids for primary and image band. NB: image band is mapped to
%=== frequencies in the primary range.
if fprimary > lo
  upper = 1;
  fprim = f(find(f>lo));
  fimband = flipud(2*lo-f(find(f<lo)));
elseif fprimary < lo
  upper = 0;
  fprim = f(find(f<lo));
  fimband = flipud(2*lo-f(find(f>lo)));
else
  error('The primary band frequency must deviate from the LO frequency');
end


%=== Determine lowest and highest frequencies in the final primary band
%=== and save only frequencies within this range.
f_min = max(fprim(1),fimband(1));
f_max = min(fprim(end),fimband(end));
fprim = fprim(find(fprim>=f_min & fprim<=f_max));
fimband = fimband(find(fimband>=f_min & fimband<=f_max));


%=== Check that the whole primary band is preserved
%=== FIXME: This could be done in a better way, should be errors?
%% if f_min>fprim(1)
%%   warning(sprintf('The primary band is cut off by %.e3 Hz in the upper end due to insufficient image band grid points', f_min-fprim(1)));
%% end
%% if f_max<fprim(end)
%%   warning(sprintf('The primary band is cut off by %.e3 Hz in the upper end due to insufficient image band grid points', fprim(end)-f_max));
%% end


out(1,1);
out(1,'Setting up mixer and sidenband filter transfer matrix (H).');


%=== Create new frequency grid
%
f_new = unique( [ fprim; fimband ] );
%
clear fprim fimband f_min f_max
%
n_new = length( f_new );


%=== Allocate vectors for column index and weights
%
lrow = 3 * n_new;
%
cols = zeros( 1, lrow );
wgts = zeros( 1, lrow );
%
%% nrow = 0;


%=== Construct vector for row index, this can be done independently
rows = reshape(repmat([1:n_new*nza], 3,1), 1, lrow*nza);

%=== Loop frequencies in F_NEW
%
for i = 1 : n_new


  %out(3,sprintf('Doing frequency %d of %d',i,n_new));


  %= Primary and image frequency
  %
  fprim  = f_new(i);
  fimage = 2*lo - fprim;


  %= Determine indices just below and above fprim and fimage
  %
  iprim  = unique([max(find(f<=fprim)) min(find(f>=fprim))]);
  iimage = unique([max(find(f<=fimage)) min(find(f>=fimage))]);


  %=== Get sideband filter values
  %
  w = interp1(f_filter, w_filter, [fprim, fimage], 'linear');


  %= Fill rows, cols and wgts with row- and col indices and weights
  icol = [1:3]+(i-1)*3;
%%   rows(icol) = n_new;
  if length(iprim)==1 && length(iimage)==2
    cols(icol) = [iprim iimage];
    itw = flipud(abs(f(iimage)-fimage))./diff(f(iimage));
    wgts(icol) = [w(1) w(2)*itw']/sum(w);
  elseif length(iprim)==2 && length(iimage)==1
    cols(icol) = [iprim iimage];
    itw = flipud(abs(f(iprim)-fprim))./diff(f(iprim));
    wgts(icol) = [w(1)*itw' w(2)]/sum(w);
  else
    %=== A zero is added at iimage for the missing weight, this is
    %=== handled by sparse by adding zero to the already existing
    %=== element.
    cols(icol) = [iprim iimage iimage];
    wgts(icol) = [w(1) w(2) 0]/sum(w);
  %=== FIXME: Does this need to be here?
%%     error('Mirrored frequency grids does not match original grid');
  end

end


%=== Repeat rows, cols and wgts for all zenith angles
cshift = reshape(repmat([0:nza-1]*nf, [3*n_new 1]), 1, 3*nza*n_new);
cols = repmat(cols, [1 nza])+cshift;
wgts = repmat(wgts, [1 nza]);


%=== Allocate H
%
H = sparse( rows, cols, wgts, n_new*nza, nf*nza );


out(1,-1);


