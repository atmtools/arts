%-----------------------------------------------------------------------------
% NAME:     psurface_crossing
%
%           Calculates the crossing point between a geometrical propagation
%           path and a pressure surface for 2D atmospheres.
%
%           The algorithm to find the crossing point is described in AUG.
%
%           If no crossing point is found, ALPHA is set to 720 degrees.
%
%           Note that RP, PSI and R1 shall be given for the same latitude.
%           The returned angle (ALPHA) is the distance from that latitude.
%
% FORMAT:   alpha = psurface_crossing( rp, psi, r1, c, the_same_surface )
%
% OUT:      alpha   The latitude distance between the starting and crossing
%                   points.
% IN:       rp      The radius of the path the latitude of r1.
%           psi     The zenith angle of the path at the same point.
%           r1      The radius for the pressure surface at the latitude of rp. 
%           c       The slope of the pressure surface. Unit is m/degree.
%                   Increasing radius with incresing latitides is defined as
%                   a positive slope.
%           the_same_surface
%                   If this argument is set to 1, it is assumed that the path
%                   starts at the pressure surface and the value of 0 is
%                   rejected for ALPHA.   
%-----------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function alpha = psurface_crossing( rp, psi, r1, c, the_same_surface )

  global DEG2RAD RAD2DEG
  
  beta = DEG2RAD * ( 180 - abs(psi) );

  c = RAD2DEG * c;

  if( psi < 0 )
    c = -c;
  end
  
  p = zeros(5,1);
  
  cv = cos( beta );
  sv = sin( beta );
  
  p(1) = -c * cv / 6;
  p(2) = -r1 * cv / 6 - c * sv / 2;
  p(3) = -r1 * sv / 2 + c * cv;
  p(4) = r1 * cv + c * sv;
  p(5) = ( r1 - rp ) * sv;

  r = roots( p );

  if the_same_surface
    alpha = RAD2DEG * min( r( find( imag(r)==0 & real(r)>0 ) ) );
  else  
    alpha = RAD2DEG * min( r( find( imag(r)==0 & real(r)>=0 ) ) );
  end

  if isempty( alpha )
    alpha = 720;
  end

  if( psi < 0 )
    alpha = -alpha;
  end
return
