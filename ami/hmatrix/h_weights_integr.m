%------------------------------------------------------------------------
% NAME:    h_weights_integr
%
%          Determines weights for performing an integration as a matrix
%          multiplication. With other words:
%          We want to calculate the integral of the product of G and H 
%          by a vector multiplication as W*G
%          The function H is known, but G is unknown.
%          Both G and H are assumed to be linear between the grid points
%
% FORMAT:  w = h_weights_integr(fg,fh,h)
%
% RETURN:  w          matrix weights
% IN:      fg         position (e.g. frequency) of values of G
%          fh         position of values of H
%          h          values of the H function (at fh)
%------------------------------------------------------------------------

% HISTORY: 12.11.99  Created for Skuld by Patrick Eriksson. 
%          25.08.00  Adapted to AMI by Patrick Eriksson


function w = h_weights_integr(fg,fh,h)

%=== Main sizes
ng   = length(fg);
nh  = length(fh);
w      = zeros(1,ng);


%=== Check input
if ( length(h) ~= nh )
  error('Lengths integration grid and function do not match');
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

  % Express G between fg(ig) and fg(ig+1) as g = g(i)(a+bf) + g(i+1)(c-bf) 
  b  = -1 / ( fg(ig+1) - fg(ig) ) ;
  a  = -fg(ig+1) * b ;
  c  = fg(ig) * b ;
  
  % Express H between fh(ih) and fh(ih+1) as h = d + ef 
  e  = ( h(ih+1)-h(ih) ) / ( fh(ih+1)-fh(ih) ) ; 
  d  = h(ih) - e*fh(ih) ;

  % Add to the weights at fg(ig) and fg(ig+1)
  w(ig)   = w(ig)   + f2 * ( a*d + f2*(b*d+a*e)/2 + f2*f2*b*e/3 ) ...
                    - f1 * ( a*d + f1*(b*d+a*e)/2 + f1*f1*b*e/3 ) ;
  w(ig+1) = w(ig+1) + f2 * ( c*d + f2*(c*e-b*d)/2 - f2*f2*b*e/3) ...
                    - f1 * ( c*d + f1*(c*e-b*d)/2 - f1*f1*b*e/3) ;

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
