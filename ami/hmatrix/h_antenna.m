%------------------------------------------------------------------------
% NAME:    h_antenna
%
%          Returns the H matrix for an antenna pattern with the options to
%          scale the pattern with frequency and to consider a moving antenna.
%          The zenith angle grid valid for the antenna is also returned.
%          The response of the antenna pattern is normalised and the
%          antenna values (W_ANT) do not need to be normalised.
%
% FORMAT:  [H,za_new] = h_antenna(f,za1,za2,za_ant,w_ant,o_ant,o_y,...
%                                                   fscale,f0,move,dza,report)
%
% RETURN:  H           antenna H matrix
%          za_new      new zenith angle grid (=ZA2)
% IN:      f           frequencies
%          za1         pencil beam zenith angles
%          za2         zenith angles observed by the sensor
%          za_ant      grid points for the antenna pattern
%          w_ant       antenna pattern
%          o_ant       linear (=1) or cubic (=3) treatment of the antenna 
%                      pattern
%          o_y         linear (=1) or cubic (=3) treatment of spectra
%          fscale      flag to scale the pattern with frequency
%          f0          reference frequency for frequency scaling, i.e. for 
%                      which frequency FWHM is valid
%          move        flag to consider a moving antenna with a constant
%                      scanning velocity during the integration
%          dza         total movement during the integration [deg]
%------------------------------------------------------------------------

% HISTORY: 99.11.12  Created for Skuld by Patrick Eriksson. 
%          00.08.25  Adapted to AMI by Patrick Eriksson
%          00.11.16  Included linear/cubic flags (PE) 


function [H,za_new] = ...
         h_antenna(f,za1,za2,za_ant,w_ant,o_ant,o_y,fscale,f0,move,dza,report)


za1    = vec2col(za1);
za_ant = vec2col(za_ant);
w_ant  = vec2col(w_ant);


%=== Main sizes
n1    = length(za1);
n2    = length(za2);
nf    = length(f);
nant  = length(za_ant);


%=== Check some lenghts
if ( n1 <= o_y )
  error('The number of pencil beam spectra must be > the selected order.')
end
if ( nant <= o_ant )
  error('The number of antenna points must be > the selected order.')
end


%=== Include possible effect of moving antenna
if move
  ant = moving_ant(za_ant,w_ant,dza);
end


%=== Check if antenna pattern totally inside pencil beam grid
if za1(1) > za2(1)+za_ant(1)
  error(sprintf('You must increase your pencil beam grid downwards with %.2f degs.',za1(1)-(za2(1)+za_ant(1))));
end
if za1(n1) < za2(n2)+za_ant(nant)
  error(sprintf('You must increase your pencil beam grid upwards with %.2f degs.',(za2(n2)+za_ant(nant))-za1(n1)));
end


out(1,1);
out(1,'Setting up antenna transfer matrix (H).');
if move
  out(2,sprintf('Antenna movement of %.3f deg. is included.',dza));
end
if fscale
  out(2,sprintf('Frequency scaling is performed (f0=%.3fGHz).',f0/1e9));
end


%=== Allocate vectors for row and col index and vector for weights
%=== Assume that the antenna pattern is described with, on average, 12 values.
lrow = 12 * n2 * nf;       % This variable is used as length for these vectors
nrow = 0;                  % Number of values in the vectors
rows = zeros( 1, lrow );
cols = zeros( 1, lrow );
wgts = zeros( 1, lrow );


%=== Fill H
%
%= With  frequency scaling
if fscale
  for i = 1:n2
    out(2,sprintf('Doing angle %d of %d',i,n2));
    za   = za1 - za2(i);
    ind1 = (i-1)*nf;
    ind2 = 0:nf:(nf*(n1-1));
    for j = 1:nf
      w    = h_weights_integr(za,o_y,za_ant*f0/f(j),o_ant,w_ant);
      s    = sum(w);
      if s > 0
        w    = w/s;
        ind3 = find(w~=0);
        nw   = length(ind3);
        if (nrow+nw) > lrow
          [rows,cols,wgts,lrow] = reallocate(rows,cols,wgts,lrow,i,n2);
        end
        irow = nrow + (1:nw);
        rows(irow) = ind1+j;
        cols(irow) = ind2(ind3)+j;
        wgts(irow) = w(ind3);
        nrow       = nrow + nw;
        %H(ind1+j,ind2(ind3)+j) = w(ind3);
      end
    end
  end

%= Without  frequency scaling
else
  for i = 1:n2
    out(2,sprintf('Doing angle %d of %d.',i,n2));
    za   = za1 - za2(i);
    ind1 = (i-1)*nf;
    ind2 = 0:nf:(nf*(n1-1));
    w    = h_weights_integr(za,o_y,za_ant,o_ant,w_ant);
    s    = sum(w);
    if s > 0
      w    = w/s;
      ind3 = find(w~=0);
      for j = 1:nf
        nw   = length(ind3);
        if (nrow+nw) > lrow
          [rows,cols,wgts,lrow] = reallocate(rows,cols,wgts,lrow,i,n2);
        end
        irow = nrow + (1:nw);
        rows(irow) = ind1+j;
        cols(irow) = ind2(ind3)+j;
        wgts(irow) = w(ind3);
        nrow       = nrow + nw;
        %H(ind1+j,ind2(ind3)+j) = w(ind3);
      end
    end
  end
end


%=== Allocate H and create ZA_NEW
H      = sparse( rows(1:nrow), cols(1:nrow), wgts(1:nrow), nf*n2, nf*n1 );
za_new = vec2col(za2);


out(1,-1);




%=========================================================================

function w_ant = moving_ant(za,w_ant,dza)

  zastep = min(diff(za))/10;

  n      = max([ceil(dza/zastep),2]);
  dzap   = linspace(-dza/2,dza/2,n);
  nza    = length(za);
  A      = zeros(n,nza);
  for i = 1:n
    A(i,:) = interp1([za(1)*2;za+dzap(i);za(nza)*2],[0;w_ant;0],za)';
  end
  w_ant  = trapz(dzap,A)'/dza;

return



function [rows,cols,wgts,lrow] = reallocate(rows,cols,wgts,lrow,i,n2);

  out(2,'Reallocates the vectors to set-up H (can take some time).');

  %=== Estimate the extra space needed. Overestimate with 25% to be safe.
  nextra = round( lrow * (n2/i-1) * 1.25 );
  zvec   = zeros( 1, nextra );
  rows   = [ rows, zvec ];        
  cols   = [ cols, zvec ];        
  wgts   = [ wgts, zvec ];        
  lrow   = lrow + nextra ;

return
