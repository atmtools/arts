%------------------------------------------------------------------------
% NAME:    h_weights_integr
%
%          Determines weights for performing an integration as a matrix
%          multiplication. With other words:
%          We want to calculate the integral of the product of G and H 
%          by a vector multiplication as W*G
%
%          The function H is known, but only the abscissa of G is known.
%
%          The two functions (G and H) can be treated to be piecewise 
%          linear or cubic (ORDERG/H = 1 or 3).
%
% FORMAT:  w = h_weights_integr(fg,orderg,fh,orderh,h)
%
% RETURN:  w          matrix weights
% IN:      fg         position (e.g. frequency) of values of G
%          orderg     assumed polynomial order for G (1 or 3)
%          fh         position of values of H
%          orderh     assumed polynomial order for H (1 or 3)
%          h          values of the H function (at fh)
%------------------------------------------------------------------------

% HISTORY: 99.11.12  Created for Skuld by Patrick Eriksson. 
%          00.08.25  Adapted to AMI by Patrick Eriksson
%          00.11.16  Adapted to handle linear/cubic (PE) 


function w = h_weights_integr(fg,orderg,fh,orderh,h)


%=== Check if orders are 1 or 3
if ~( (orderg==1) | (orderg==3) ) | ~( (orderh==1) | (orderh==3) ) 
  error('The polynomial order must be 1 or 3.');
end


%=== Main sizes
ng   = length(fg);
nh   = length(fh);
w    = zeros(1,ng);


%=== Check lengths
if ( (ng<=orderg) | (nh<=orderh) )
  error('The vector lengths must be larger then the corresponding order.');
end


%=== Check input
if ( length(h) ~= nh )
  error('Lengths of integration grid and function do not match');
end


%=== Init index and loop until freqs. of G are inside freqs. of H
ig   = 1;
ih   = 1;    
while ( (ig<ng) & (fg(ig+1)<=fh(ih)) )
  ig = ig + 1;
end
if ( ig == ng )    % G-range totally below H-range
  return
end


%=== Loop until freqs. of H are inside freqs. of G
while ( (ih<nh) & (fh(ih+1)<=fg(ig)) )
  ih = ih + 1;
end
if ( ih == nh )    % H-range totally below G-range
  return
end


%=== Determine end frequencies for first common frequency range
f1    = max([fg(ig),fh(ih)]);       % min freq.
f2    = min([fg(ig+1),fh(ih+1)]);   % max freq.



%=== Loop until end of G or H frequencies
while ( 1 )

  % Determine the lowest and highest index for the present fit 
  % range of G. If ORDERG==3, the outermost ranges are fitted by
  % second order polynomials.
  if ( orderg == 1 )
    ig1 = ig;
    ig2 = ig + 1;
  elseif ( ig == 1 )
    ig1 = ig;
    ig2 = ig + 2;
  elseif ( ig == (ng-1) )
    ig1 = ig - 1;
    ig2 = ig + 1;
  else
    ig1 = ig - 1;
    ig2 = ig + 2;
  end

  % Do the same for H.
  if ( orderh == 1 )
    ih1 = ih;
    ih2 = ih + 1;
  elseif ( ih == 1 )
    ih1 = ih;
    ih2 = ih + 2;
  elseif ( ih == (nh-1) )
    ih1 = ih - 1;
    ih2 = ih + 1;
  else
    ih1 = ih - 1;
    ih2 = ih + 2;
  end

  % Get coefficients of the polynomial basis for G
  pg = pbasis(fg(ig1:ig2));

  % Get the polynomial coefficients for H
  ph = polyfit(fh(ih1:ih2),h(ih1:ih2),ih2-ih1);

  % Multiplicate G and H
  pgh = pbasis_x_pol(pg,ph);

  % Integrate the product between F1 and F2 and add to WG
  ws   = pbasis_integrate(pgh,f1,f2);
  for j = ig1:ig2
    w(j) = w(j) + ws(j-ig1+1);
  end

  % Determine indeces for next step
  if ( fg(ig+1) == fh(ih+1) )
    ig = ig + 1;
    ih = ih + 1;
  elseif ( fg(ig+1) < fh(ih+1) )
    ig = ig + 1;
  else
    ih = ih + 1;
  end

  % Check if end of G or H has been reached
  if ( ig == ng )
    return
  end
  if ( ih == nh )
    return
  end

  % Update frequence ranges
  f1  = f2 ;
  f2  = min([fg(ig+1),fh(ih+1)]);   % max freq.

end
