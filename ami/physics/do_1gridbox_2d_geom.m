%------------------------------------------------------------------------
% NAME:     do_1gridbox_2d_geom
%
%           This is a function to test out 2-D the propagation paths 
%           calculations for arts-2. This function will be ported to arts
%           and no further documentation is given here.
%------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


function [r,alpha,l_step,psi,r_istep,alpha_istep,r_tan,alpha_tan,ground,i_ground] = do_1gridbox_2d_geom( r1, r2, r3, r4, alpha1, alpha3, rground1, rground2, blackbody_ground, r_p, alpha_p, psi_p, l_max )

  global DEG2RAD RAD2DEG

  % Asserts 
  % That r_p is inside the box is asserted below.
  assert( r1<=r3, 'r1<=r3' );
  assert( r2<=r4, 'r2<=r4' );
  assert( alpha1<=alpha3, 'alpha1<=alpha3' );
  assert( alpha_p>=alpha1, 'alpha_p>=alpha1' );
  assert( alpha_p<=alpha3, 'alpha_p<=alpha3' );
  assert( abs(psi_p)<=180, 'abs(psi_p)<=180' );

  % Set tangent point variables to NaN
  r_tan     = NaN;
  alpha_tan = NaN;

  % Angular distance from face 1 and 3, and from start point to alpha1/3.
  dalpha13 = alpha3 - alpha1;
  dalpha_l = alpha_p - alpha1;
  dalpha_r = alpha3 - alpha_p;

  % Slope of the pressure surfaces, corresponding to face 2 and 4, and the 
  % ground. The slope is calculated from alpha1 to alpha3. Unit is m/degree.
  c2      = ( r2 - r1 ) / dalpha13;
  c4      = ( r3 - r4 ) / dalpha13;
  cground = ( rground2 - rground1 ) / dalpha13;
  %
  assert( r_p <= r4+c4*dalpha_l, 'r_p<=r4+c4*dalpha_l' );
  assert( r_p >= r1+c2*dalpha_l, 'r_p>=r2+c2*dalpha_l' );

  % Cosine and sine of the psi_p angle
  cv  = cos( DEG2RAD * psi_p );
  sv  = sin( DEG2RAD * psi_p );

  % The value 1 in can_be_end means that the path can end at that face.
  % Face 1 and 3 can easily be sorted out for some cases.
  can_be_end = ones(4,1);
  if( psi_p==0 | abs(psi_p)==180 )
    can_be_end(1) = 0;
    can_be_end(3) = 0;
  elseif( psi_p > 0 )
    can_be_end(1) = 0;
  else
    can_be_end(3) = 0;
  end    

  % Change of box index in the r and alpha directions.
  r_istep     = 0;
  alpha_istep = 0;

  % Check if the start point is on one of the faces, and if this is the case,
  % that the zenith angle is inside the box and remove impossible end faces
  % from can_be_end.
  % Corner points are considered to belong to face 2 and 4.
  
  % Starting at face 2? (Cannot start at 1 or 3 as corner points belongs to 2)
  if( r_p == r1 + c2*dalpha_l )
    tilt = tilt_of_psurface( r_p, c2 );
    if( psi_p > 0 )
      assert( psi_p<=90-tilt, 'psi_p<=90-tilt' );
      can_be_end(1) = 0;
    else
      assert( psi_p>=-90-tilt, 'psi_p>=-90-tilt' );
      can_be_end(3) = 0;
    end
    can_be_end(2) = 0;
    % Handle the corner points
    if( alpha_p == alpha1 )
      assert( psi_p>=0, 'psi_p>=0' );
      can_be_end(1) = 0;
    elseif( alpha_p == alpha3 )
      assert( psi_p<=0, 'psi_p<=0' );
      can_be_end(3) = 0;
    end

  % Starting at face 4? 
  % This face is special as the end face can equal the start face.
  elseif( r_p == r4 + c4*dalpha_l )
    tilt = tilt_of_psurface( r_p, c2 );
    if( psi_p > 0 )
      assert( psi_p>90-tilt, 'psi_p>90-tilt' );
      can_be_end(1) = 0;
    else
      assert( psi_p<-90-tilt, 'psi_p<-90-tilt' );
      can_be_end(3) = 0;
    end
    % Handle the corner points
    if( alpha_p == alpha1 )
      assert( psi_p>=0, 'psi_p>=0' );
      can_be_end(1) = 0;
    elseif( alpha_p == alpha3 )
      assert( psi_p<=0, 'psi_p<=0' );
      can_be_end(3) = 0;
    end

  % Starting at face 1? 
  elseif( alpha_p == alpha1 )
    assert( psi_p>=0, 'psi_p>=0' );
    can_be_end(1) = 0;

  % Starting at face 3? 
  elseif( alpha_p == alpha3 )
    assert( psi_p<=0, 'psi_p<=0' );
    can_be_end(3) = 0;
  end

  % Determine the end face.
  % Calculations are first done ignoring a possible ground intersection.

  % Zenith angle to the corner points
  psi_c    = zeros(4,1); 
  psi_c(1) = za_geom2other_point( r_p, r1, -dalpha_l );
  if( psi_c(1 ) == 180 )
    psi_c(1 ) = -180;
  end    
  psi_c(2) = za_geom2other_point( r_p, r2, dalpha_r );
  psi_c(3) = za_geom2other_point( r_p, r3, dalpha_r );
  psi_c(4) = za_geom2other_point( r_p, r4, -dalpha_l );

  % Face 2 needs special care as it bended inwards and can intersect with
  % the path to face 1 or 3. Calculate angular distance to a possible
  % intersection with face 2.
  dalpha2 = 720;
  if( can_be_end(2) )
    dalpha2 = abs( psurface_crossing( r_p, psi_p, r1+c2*dalpha_l, c2, 0 ) );
  end

  % End face is 1?
  if( can_be_end(1) & psi_p>psi_c(1) & psi_p<psi_c(4) & dalpha2>dalpha_l )
    end_face    = 1;
    dalpha      = -dalpha_l;
    alpha_istep = -1;

  % End face is 3?
  elseif( can_be_end(3) & psi_p<psi_c(2) & psi_p>psi_c(3) & dalpha2>dalpha_r )
    end_face    = 3;
    dalpha      = dalpha_r;
    alpha_istep = 1;
  
  % End face is 4?
  % Again this face is special as start and end face can be the same.
  elseif( can_be_end(4) & psi_p>=psi_c(4) & psi_p<=psi_c(3) )
    end_face   = 4;
    from4to4   = 0;
    if( r_p == r4 + c4*dalpha_l )
      from4to4 = 1;
    end
    dalpha = psurface_crossing( r_p, psi_p, r4+c4*dalpha_l, c4, from4to4);
    r_istep = 1;

  % Only remaining choice is face 2!
  else             %if( can_be_end(2) & ( psi_p<=psi_c(1) | psi_p>=psi_c(2) ) )
    end_face   = 2;
    dalpha = psurface_crossing( r_p, psi_p, r1+c2*dalpha_l, c2, 0 );
    r_istep = -1;
  end

  % Avoid the case of a tangent point exactly at a pressure surface.
  % This is very unlikely, but if it would happen there would be a crash.
  % This can only happen for when the end face is 2. 
  % Ignore then the touching of the pressure surface (a point will be placed 
  % at the tangent point anyhow).
  if( end_face==2 & abs(psi_p-dalpha)-90 == 0 )
    if( psi_p > 0 )
      end_face = 3;
      dalpha   = alpha3 - alpha_p;
    else
      end_face = 1;
      dalpha   = alpha1 - alpha_p;
    end
  end

  % Check if there is a intersection with the ground. If yes, we calculate
  % the PPATH to the ground and continues after the reflection later.
  ground = 0;
  if( rground1>=r1 | rground2>=r2 )
    rground = rground1 + cground * ( alpha_p - alpha1 );
    assert( r_p>rground, 'r_p>rground' );
    dag = psurface_crossing( r_p, psi_p, rground, cground, 0 );
    if( abs(dag) <= abs(dalpha) )
      end_face = 5;
      ground   = 1;
      dalpha   = dag;
    end
  end

  % Handle possible numerical inaccaracy for psi_p=0
  if( psi_p==0 | psi_p==180 )
    dalpha = 0;
  end 

  % Calculate position and angle of end point.
  % On the same time check if we have not only moved vertically or 
  % horisontally to a new grid box. This can happen if the end point is 
  % exactly a corner point, or due to numerical inaccuracy.
  alpha_e = alpha_p + dalpha;
  psi_e   = psi_p   - dalpha;
  if( end_face == 1 )
    alphe_e = alpha1;                     % To be as exact as possible
    r_e     = r_p * sv / sin( DEG2RAD*psi_e );
    if( r_e < r1 )
      r_istep = r_istep - 1;
    elseif( r_e > r4 )
      r_istep = r_istep + 1;
    end
  elseif( end_face == 2 )
    r_e     = r1 + c2 * ( alpha_e - alpha1 );
    if( alpha_e < alpha1 )
      alpha_istep = alpha_istep - 1;
    elseif( alpha_e > alpha3 )
      alpha_istep = alpha_istep + 1;
    end
  elseif( end_face == 3 )
    alphe_e = alpha3;                     % To be as exact as possible
    r_e     = r_p * sv / sin( DEG2RAD*psi_e );
    if( r_e < r2 )
      r_istep = r_istep - 1;
    elseif( r_e > r3 )
      r_istep = r_istep + 1;
    end
  elseif( end_face == 4 )
    r_e     = r4 + c4 * ( alpha_e - alpha1 );
    if( alpha_e < alpha1 )
      alpha_istep = alpha_istep - 1;
    elseif( alpha_e > alpha3 )
      alpha_istep = alpha_istep + 1;
    end
  elseif( end_face == 5 )   % That is, the ground
    r_e     = rground1 + cground * ( alpha_e - alpha1 );
  end 
  if dalpha ~= 0
    l_e  = r_e * sin( DEG2RAD*dalpha ) / sv;
  else
    l_e = abs( r_p - r_e );
  end

  % If l_e is zero something has gone wrong. Probably starting from a corner
  % going out from the box.
  assert( l_e>0, 'l_e>0' );

  % Create a vector with lengths along the PPATH, including starting point.
  % The PPATH is divided into parts if l_e > l_max.
  np = ceil( l_e / l_max );
  l  = linspace( 0, l_e, np+1 );

  % Check if we have passed a tangent point. 
  % If yes, put in the tangent point in l.
  if( abs(psi_p)>=90 & abs(psi_e)<90 )
    l_tan = r_p * cos( DEG2RAD*( 180 - psi_p ) );
    if( all( abs(l-l_tan)/l_max > 1e-5 ) )
      l     = sort( [l,l_tan] );
      np    = np + 1;
    end
    r_tan     = abs( r_p * sv );
    if( psi_p > 0 )
      alpha_tan = alpha_p + psi_p - 90;
    else
      alpha_tan = alpha_p + psi_p + 90;
    end
  end

  % Allocate return vectors and fill 
  r      = zeros(np,1);
  alpha  = zeros(np,1);
  l_step = zeros(np,1);
  psi    = zeros(np,1);
  %
  r(np)      = r_e;
  alpha(np)  = alpha_e;
  l_step(np) = l(np+1) - l(np);
  psi(np)    = psi_e;
  %
  if( np > 1 )
    a  = r_p * r_p;
    for i = 1:np-1
      r(i)      = sqrt( a + l(i+1)*l(i+1) + 2 * r_p * l(i+1) * cv );
      dalpha    = RAD2DEG*asin( l(i+1)*sv / r(i) );
      alpha(i)  = alpha_p + dalpha;
      l_step(i) = l(i+1) - l(i);  
      psi(i)    = psi_p - dalpha;
    end
  end

  % Continue to follow the PPATH if ground intersection and not a blackbody 
  % ground (blackbody_ground). If blackbody_ground, we are ready as the PPATH 
  % then starts there.
  if( ground )
    i_ground = np;
    if( ~blackbody_ground )
      i_ground = np;
      if( psi_e > 0 )
	psi_r  = 180 - psi_e - 2 * tilt_of_psurface( rground1, cground );
      else
	psi_r  = -180 - psi_e - 2 * tilt_of_psurface( rground1, cground );
      end
      % If we are exactly at the alpha boundaries and there is a ground slope,
      % we can bounce out from the box. Fix this with a small disturbance of
      % the angle.
      if( alpha_e == alpha1 )
        alpha_e = alpha1 + (alpha3-alpha1)/1e6;
      elseif( alpha_e == alpha3 )
        alpha_e = alpha3 - (alpha3-alpha1)/1e6;
      end
      [r2,alpha2,l_step2,psi2,r_istep,alpha_istep] = do_1gridbox_2d_geom( r1, r2, r3, r4, alpha1, alpha3, 0, 0, 0, r_e, alpha_e, psi_r, l_max );
      r      = [ r; r2 ];
      alpha  = [ alpha; alpha2 ];
      l_step = [ l_step; l_step2 ];
      psi    = [ psi; psi2 ];
    end
  else
    i_ground = 0;
  end
return



function assert( bool, m )
  if ~bool
    error(m);
  end 
return






