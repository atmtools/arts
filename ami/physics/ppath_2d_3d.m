%------------------------------------------------------------------------
% NAME:     ppath_2d_3d
%
%           This is a function to test out 2-D the propagation paths 
%           calculations for arts-2. This function will be ported to arts
%           and no further documentation is given here.
%
%           The function will only be part of AMI as long this Matlab
%           version identical to the C++ version.
%------------------------------------------------------------------------

% HISTORY: 2002-03-08  Created by Patrick Eriksson


% Limitations for OK propagation paths:
% 1. The propagation path is not allowed to cross the latitude (or longitude
% for 3D) end faces. That is, the path must leave the defined atmopshere 
% at the top. 
% 2. Observations in an upward direction (abs(psi_s) <= 90 ) are only allowed
% from points inside the atmosphere.
% 3. For downlooking observations (abs(psi_s) > 90) the (geometrical) 
% tangent point must be inside the covered latitude and longitude ranges 
% (but can be above the top of the atmosphere).

function ppath = ppath_2d_3d(dim,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,refr_on,blackbody_ground,cloudbox_on,cloudbox_limits,ppath_lmax,r_s,alpha_s,beta_s,psi_s,omega_psi)

  if( dim == 3 )
    error('3D not yet implemented');
  end

  global DEG2RAD RAD2DEG

  % Asserts ...
  % dim <= 3

  % Some sizes
  n_p     = length( p_grid );
  n_alpha = length( alpha_grid );
  if( dim >= 3 )
    n_beta = length( beta_grid );
  end

  % Flag for the path is inside the atmosphere (yes is default).
  do_path = 1;

  % Start by checking if the sensor is inside the atmosphere, but not 
  % below the ground.
  is_inside = is_inside_atm_box( dim, p_grid, alpha_grid, beta_grid, ...
          z_field, r_geoid, z_ground, r_s, alpha_s, beta_s, psi_s, omega_psi );


  %--- Start point inside the atmosphere --------------------------------------
  if( is_inside )

    %- LOS start values
    r1     = r_s;
    alpha1 = alpha_s;  % For dim>=3, beta1 and omega1 are set below
    psi1   = psi_s;    

    %- Find in what grid box the sensor is placed
    if( dim == 2 )
      %
      if( psi_s >= 0 )
	ialpha1 = find_grid_range( alpha_grid, alpha1, 1 ); 
      else
	ialpha1 = find_grid_range( alpha_grid, alpha1, -1 ); 
      end
      %
      r_g    = interp1( alpha_grid, r_geoid(1,:), alpha1 );
      r_grid = r_g + interp1( alpha_grid, squeeze(z_field(1,:,:)), alpha1 );
      if( abs(psi_s) <= 90 )
	iz1 = find_grid_range( r_grid, r1, 1 ); 
      else
	iz1 = find_grid_range( r_grid, r1, -1 ); 
      end

    else
      %
      beta1  = beta_s;
      omega1 = omega_s;    
      %

    end


  %--- Start point outside the atmosphere -------------------------------------
  else

    % Handle cases when the sensor appears to look in the wrong direction.
    if( dim == 2 )
      if( (alpha_s<=alpha_grid(1) & psi_s<0 & psi_s>=-180 ) | ...    
                        (alpha_s>=alpha_grid(n_alpha) & psi_s>0 & psi_s<=180) )
        s = sprintf('The sensor is outside (or at the limit) of the defined atmosphere\nbut looks in the wrong direction (wrong sign for psi?).');
        s = sprintf('%s\nThis case includes nadir looking exactly at the latitude end points.',s);
        error(s);
      end
    else
      %
    end

    %--- If upward, error -----------------------------------------------------
    if( abs(psi_s) <= 90 ) 
      s = sprintf('Upward observations are only allowed from inside the defined atmosphere,\nbut the sensor is outside (or exactly at the top of) this atmosphere.');
      error(s);


    %--- Downward -------------------------------------------------------------
    else

      % The geometrical tangent position
      if( dim == 2 )
        alpha_tan = alpha_s + sign(psi_s) * ( abs(psi_s) - 90 );
      else
        %
      end

      % Is geometrical tangent altitude above the top of the defined 
      % atmosphere? If yes, this OK and the observed spectrum equals the
      % cosmic background raditiaton.
      if( alpha_tan>=alpha_grid(1) & alpha_tan<=alpha_grid(n_alpha) )
      	%
      	if( dim == 2 )
      	  r_g   = interp1( alpha_grid, r_geoid, alpha_tan );
      	  z_top = interp1( alpha_grid, squeeze(z_field(1,:,n_p)), alpha_tan );
      	else
      	  %
      	end
      	%
        z_tan = r_s * sin(DEG2RAD*abs(psi_s)) - r_g;
        %
      	if( z_tan > z_top )
      	  do_path = 0;
        end
      end     % do_path=0 ?


      %- Do something only if do_path=1 ---------------------------------------
      if( do_path )

	%--- Bad paths for 2D -------------------------------------------------
	if( dim == 2 )
  
	  % Define some help variables, also used to calculate the path
	  if( psi_s >= 0 )       % Looking against higher latitudes
	    ialpha1 = 1;          % Index of first corner point in the view dir
	    istep   = 1;          % Index increament in viewing direction
            s_end   = 'lower';    % First lat end in viewing direction
	  else         	         % Looking against lower latitudes
	    ialpha1 = n_alpha;
	    istep   = -1;
            s_end   = 'upper';
	  end
	  %
	  psi0 = abs( psi_s );
	  %
	  a = r_s * sin( DEG2RAD*psi0 ); 
  
          % Handle cases the sensor and the tangent point both are outside (on
	  % the same side) of the covered latitude range.
	  if( istep*alpha_s<istep*alpha_grid(ialpha1) & ...
                                    istep*alpha_tan<istep*alpha_grid(ialpha1) )
	    s = 'The sensor is too far away from the covered latitude range.';
	    s = sprintf('%s\nThe latitude grid must be increased at the %s end.',s,s_end);
	    error(s);
	  end
  
	  % If the sensor is outside the latitude range, check that the path
	  % goes over the closest upper corner of the defined atmosphere.
	  if( istep*alpha_s < istep*alpha_grid(ialpha1) )
	    dalpha = abs( alpha_s - alpha_grid(ialpha1) );
	    rppath  = a / sin( DEG2RAD*(psi0-dalpha) );
	    if( rppath <  r_geoid(1,ialpha1)+z_field(1,ialpha1,n_p) )
	      if( r_s > r_geoid(1,ialpha1)+z_field(1,ialpha1,n_p) )
		  dalpha = dalpha - ( psi0 - 180 + RAD2DEG * ...
                     asin( a / (r_geoid(1,ialpha1)+z_field(1,ialpha1,n_p)) ) );
	      end
	      s = sprintf('The propagation path enters the atmosphere at the %s latitude end face.',s_end);
	      s = sprintf('%s\nYou need to increase the latitude grid with, at least, %.3f degrees at\nthe %s end.',s,dalpha,s_end);
	      error( s );
	    end
	  end

	  % Check that the path does not go over the end (in the viewing
	  % direction) corner point. This can only happen without prevoius
          % error if alpha_tan is outside the covered latitude range.
          % OK paths going over this corner point are already handled above 
          % by setting do_path=0.
          ic2    = n_alpha + 1 - ialpha1;         % Index of the other corner
          if( istep*alpha_tan > istep*alpha_grid(ic2) )
  	    dalpha = abs( alpha_s - alpha_grid(ic2) );
	    rppath = a / sin( DEG2RAD*(psi0-dalpha) );
            rv1    = r_geoid(1,ic2) + z_field(1,ic2,n_p);
            if( rppath >= rv1 )
              if( istep > 0 )
                s_end   = 'upper';
              else
                s_end   = 'lower';
              end
	      s = sprintf('The propagation path passes over the atmosphere but the\ntangent point is not above the covered latitude range.');
	      s = sprintf('%s\nThe path goes %.3f km above the %s end corner and the',s,(rppath-rv1)/1e3,s_end);
	      s = sprintf('%s\ntangent point is %.3f degrees outside the %s latitude limit.',s,abs(alpha_tan-alpha_grid(ic2)),s_end);
	      error( s );
	    end
          end
  
  
	%--- Bad paths for 3D -----------------------------------------------
	else    
	  % 
	end   % Bad paths


        % Find where the path enters the defined atmosphere. Geometrical 
        % calculations can be used as this part takes place in space.

        %- 2D path ------------------------------------------------------------
        if( dim == 2 )

	  % Find the latitude index just before the sensor, in the viewing 
          % direction, or end index if the sensor is outside the grid.
          while( istep*alpha_s > istep*alpha_grid(ialpha1+istep) )
	    ialpha1 = ialpha1 + istep;
          end

          % Step forward until the face where the path enters the atmosphere
          % is found. The enter face extends between ialpha and ialpha+istep,
          % thus ialpha is the corner point closest to the sensor.
          % The rare case that the path both enters and leaves the atmosphere
	  % inside 1 face is covered by check if the tangent point is passed.
	  rv1 = -1;
          while( rv1 < 0 )
            if( istep*alpha_tan > istep*alpha_grid(ialpha1+istep) )
    	      dalpha = abs( alpha_s - alpha_grid(ialpha1+istep) );
	      rppath = a / sin( DEG2RAD*(psi0-dalpha) );
              rv1    = r_geoid(1,ialpha1+istep) + z_field(1,ialpha1+istep,n_p);
              if( rppath > rv1 )
	        ialpha1 = ialpha1 + istep;
	        rv1     = -1;
              end
            else
	      rv1 = 1;
	    end
	  end

	  % Calculate the angular distance between alpha_grid(ialpha) and
	  % the entering point. Note that dalpha is used here and above for
	  % different angular distances.

     	  dalpha = abs( alpha_s - alpha_grid(ialpha1) );
          % If the sensor is inside the latitude range, reverse sign of dalpha.
          if( istep*alpha_s >= istep*alpha_grid(ialpha1) )
            dalpha = -dalpha;
	  end
	  rppath = a / sin( DEG2RAD*(psi0-dalpha) );
          rv1    = r_geoid(1,ialpha1) + z_field(1,ialpha1,n_p);
          rv2    = r_geoid(1,ialpha1+istep) + z_field(1,ialpha1+istep,n_p);
          c      = ( rv2 - rv1 ) / ...
                        abs( alpha_grid(ialpha1+istep) - alpha_grid(ialpha1) );
	  dalpha = psurface_crossing( rppath, psi0-dalpha, rv1, c, 0 );
          % As psurface_crossing is not 100% exact, we better check for 
	  % obvoiusly wrong results.
	  if( dalpha >= abs(alpha_grid(ialpha1+istep)-alpha_grid(ialpha1)) )
	    dalpha = 0.9999*abs(alpha_grid(ialpha1+istep)-alpha_grid(ialpha1));
          end

	  % Now we can calculate the coordinates for the effective start point
	  % of the propagation path.
          r1     = rv1 + c * dalpha;
	  alpha1 = alpha_grid(ialpha1) + istep * dalpha;
	  psi1   = psi_s + alpha_s - alpha1;
	  iz1    = n_p - 1;
          if( istep < 0 )
	    ialpha1 = ialpha1 - 1;
	  end


        %- 3D path ------------------------------------------------------------
        else
          %
        end   % 2/3D paths

      end     % do_path = 1;
    end       % Upward/downward
  end         % Start inside/outside the atmopshere


  % Create  ppath
  %
  ppath = empty_path( dim );
  %
  if( ~do_path )
    %
    ppath.tan_pos    = zeros( dim, 1 );
    ppath.tan_pos(1) = z_tan;
    ppath.tan_pos(2) = alpha_tan;
    if( dim >= 3 )
      ppath.tan_pos(3) = beta_tan;
    end

  else

    % Allocate vector arrays that can hold the path inside each grid box.
    % The starting point is put in as box 1.
    %
    if( dim == 2 )
      nmax = 2*n_p + n_alpha + 1;
    else
      %
      % beta  =
    end
    %
    r        = cell( nmax, 1 );     r{1}      = r1;
    alpha    = cell( nmax, 1 );     alpha{1}  = alpha1; 
    l_step   = cell( nmax-1, 1 );
    psi      = cell( nmax, 1 );     psi{1}    = psi1;
    iz       = zeros( nmax, 1 );    iz(1)     = iz1;
    ialpha   = zeros( nmax, 1 );    ialpha(1) = ialpha1;
    if( dim == 2 )
      ibeta1 = 1;
    end
    ibeta = zeros( nmax, 1 );   ibeta(1) = ibeta1;

    % Counting variables
    nb        = 1;   % Number of boxes
    np        = 1;   % Number of points in ppath
    not_ready = 1;   % End of path reached?

    % If the start point is inside an active cloud box there is no path 
    % to follow. 
    if( cloudbox_on )
      l1=cloudbox_limits(1); l2=cloudbox_limits(2); l3=cloudbox_limits(3);
      if( dim >= 3 )
        l4=cloudbox_limits(4); l5=cloudbox_limits(5);
        bgrid = beta_grid( l4:l5 );
      else
        l4=1; l5=1;
        bgrid = beta_grid;
      end
      if( is_inside_atm_box( dim, p_grid(1:l1), alpha_grid(l2:l3), ...
            bgrid, z_field(l4:l5,l2:l3,1:l1), r_geoid(l4:l5,l2:l3), ...
              z_ground(l4:l5,l2:l3), r_s, alpha_s, beta_s, psi_s, omega_psi ) )
        ppath.background = 'Inside cloud box';
        not_ready = 0;
      end
    end

    % Go to next grid box until end of path
    %
    while( not_ready )

      nb  = nb + 1;

      iz(nb)     = iz1;
      ialpha(nb) = ialpha1;
      ibeta(nb)  = ibeta1;

      rll = r_geoid(ibeta1,ialpha1)   + z_field(ibeta1,ialpha1,iz1);
      rlr = r_geoid(ibeta1,ialpha1+1) + z_field(ibeta1,ialpha1+1,iz1);
      rur = r_geoid(ibeta1,ialpha1+1) + z_field(ibeta1,ialpha1+1,iz1+1);
      rul = r_geoid(ibeta1,ialpha1)   + z_field(ibeta1,ialpha1,iz1+1);
      rg1 = r_geoid(ibeta1,ialpha1)   + z_ground(ibeta1,ialpha1);
      rg2 = r_geoid(ibeta1,ialpha1+1) + z_ground(ibeta1,ialpha1+1);

      if( dim == 2 )
        if( refr_on )
          error('2D PPATH with refraction not yet implemented.');
        else
          [r{nb},alpha{nb},l_step{nb-1},psi{nb},d_iz1,d_ialpha1,r_tan,alpha_tan,ground,i_ground] = do_1gridbox_2d_geom( rll, rlr, rur, rul, alpha_grid(ialpha1), alpha_grid(ialpha1+1), rg1, rg2, blackbody_ground, r1, alpha1, psi1, ppath_lmax );
        end
      else
        error('3D PPATH calculations not yet implemented.');
      end

      n  = length( r{nb} );
      np = np + n;

      % The xxx1 variables are here used to keep track on the presently last
      % point of the ppath.
      r1     = r{nb}(n);
      alpha1 = alpha{nb}(n);
      psi1   = psi{nb}(n);

      % Check position of new grid box?
      %
      iz1     = iz1 + d_iz1;
      ialpha1 = ialpha1 + d_ialpha1; 
      %
      if( iz1 >= n_p )
        not_ready = 0;
      end
      %
      if( ialpha1<1 | ialpha1>=n_alpha | ...
                                   ( dim>=3 & ( ibeta1<1 | ibeta1>=n_beta ) ) )
        s = 'The exit point of the propagation path is not at the top of the atmosphere.';
        if( dim==2 )
          r_g = interp1( alpha_grid, r_geoid, alpha1 );
        else
          %
        end
        s = sprintf('%s\nThe exit point is at:\naltitude = %.3f km\nlatitude = %.3f',s,(r1-r_g)/1e3,alpha1);
        if( dim >= 3 )
          s = sprintf('%s\nlongitude = %.3f',s,beta1);
        end       
        error(s);
      end

      % Intersection with ground?
      if( ground )
	if( blackbody_ground )
          not_ready        = 0;
          ppath.background = 'Blackbody ground';
        else
          ppath.ground     = 1;
          ppath.i_ground   = np - n + i_ground;
        end
      end

      % Intersection with cloud box?
      if( cloudbox_on & iz1<cloudbox_limits(1) )
        if( ialpha1>=cloudbox_limits(2) & ialpha1<cloudbox_limits(3) )
          if( dim==2 | ( dim>=3 & ibeta1>=cloudbox_limits(4) & ...
                                                  ibeta1<cloudbox_limits(5) ) )
            not_ready        = 0;
            ppath.background = 'Surface of cloud box';
          end
        end
      end

      % Tangent point?
      if( ~isnan( r_tan ) )
        ppath.tan_pos    = zeros( dim, 1 );
        ppath.tan_pos(1) = r_tan - interp1( alpha_grid, r_geoid, alpha_tan );
        ppath.tan_pos(2) = alpha_tan;
        if( dim >= 3 )
          ppath.tan_pos(3) = beta_tan;
        end
      end
    end   % while

    ppath.np         = np;
    ppath.i_start    = np - 1;
    ppath.i_stop     = 0;
  
    ppath.pos        = zeros(np,dim);
    ppath.ip_pos     = zeros(np,dim);
    ppath.z          = zeros(np,1);
    ppath.l_step     = zeros(np-1,1);
    ppath.los        = zeros(np,dim-1);
  
    np = 0;
    for i = 1:nb
      for j = 1:length(r{i})
  
	np = np + 1;
  
	% Special stuff for dim = 3.
	if( dim >= 3 )
  
	  % beta
  
	  % los
	  ppath.los(np,2) = omega(i);
	end
  
	% alpha
	ip0                = ialpha(i);
	ppath.ip_pos(np,2) = ip0 + (alpha{i}(j)-alpha_grid(ip0)) / ... 
				       ( alpha_grid(ip0+1) - alpha_grid(ip0) );
	ppath.pos(np,2)    = alpha{i}(j);
  
	% p/z
	ip0                = iz(i);
	if( dim == 2 )
	  ppath.z(np) = r{i}(j) - interp1( alpha_grid, r_geoid, alpha{i}(j) );
	  z_bot       = interp1( alpha_grid, z_field(1,:,ip0), alpha{i}(j) );
	  z_top       = interp1( alpha_grid, z_field(1,:,ip0+1), alpha{i}(j) );
	else
	  %
	end
	ppath.ip_pos(np,1) = (ppath.z(np)-z_bot) / (z_top-z_bot); 
	ppath.pos(np,1)    = p_grid(ip0) + ppath.ip_pos(np,1) * ...
					       ( p_grid(ip0+1) - p_grid(ip0) );
	ppath.ip_pos(np,1) = ppath.ip_pos(np,1) + ip0;
  
	% los
	ppath.los(np,1)    = psi{i}(j);
  
	% l_step
	if( i > 1 )
	  ppath.l_step(np-1) = l_step{i-1}(j);
	end
      end   % for j
    end     % for i

  end       % if do_path


if 0
if do_path
  fprintf('r1      = %.3f km\n', r1/1e3 );
  fprintf('iz1     = %d\n', iz1 );
  fprintf('alpha1  = %.4f\n', alpha1 );
  fprintf('ialpha1 = %d\n', ialpha1 );
  fprintf('psi1    = %.3f\n', psi1 );
else
  fprintf('Path above the atmosphere.\n');
end
end
return


function i1 = find_grid_range( grid, x, direction )

  n = length( grid );

  if( direction > 0 )
    i1    = 1;
    istep = 1;
  else
    i1    = n;
    istep = -1;
  end

  while( i1+istep<n & istep*x >= istep*grid(i1+istep) )
    i1 = i1 + istep;
  end

  if( direction < 0 )
    i1 = i1 - 1;
  end

return




function ppath = empty_path( dim )
  ppath.dim        = dim;
  ppath.np         = 0;
  ppath.i_start    = 0;
  ppath.i_stop     = 0;
  ppath.pos        = [];
  ppath.ip_pos     = [];
  ppath.z          = [];
  ppath.l_step     = [];
  ppath.los        = [];
  ppath.background = 'Cosmic background radiation';
  ppath.ground     = 0;
  ppath.i_ground   = 0;
  ppath.tan_pos    = [];
return



function is_inside = is_inside_atm_box(dim,p_grid,alpha_grid,beta_grid,z_field,r_geoid,z_ground,r_s,alpha_s,beta_s,psi_s,omega_psi)

  is_inside      = 1;

  % Check first lat and lon position.
  % If on the box boundary, viewing direction must be inward.
  if( dim>=2 & ( alpha_s<alpha_grid(1) | alpha_s>last(alpha_grid) ) )
    is_inside = 0;
  end 
  if( dim>=3 & ( beta_s<beta_grid(1) | beta_s>last(beta_grid) ) )
    is_inside = 0;
  end 
  if( dim==2 & ( (alpha_s==alpha_grid(1) & psi_s<0) | ...
                                      (alpha_s==last(alpha_grid) & psi_s>0) ) )
    is_inside = 0;
  end 
  % Put in boundary check for dim=3.

  if( is_inside )

    n_p = size( z_field, 3 );

    if( dim == 1 )
      r_g        = r_geoid(1,1);
      z_bot      = z_ground(1,1);
      z_top      = z_field(1,1,n_p);
    elseif( dim == 2 )
      r_g        = interp1( alpha_grid, r_geoid(1,:), alpha_s );
      z_bot      = interp1( alpha_grid, z_ground(1,:), alpha_s );
      z_top      = interp1( alpha_grid, z_field(1,:,n_p), alpha_s );
    else
      %r_g   =
      %z_bot = 
      %z_top = 
    end

    if( r_g + z_bot >= r_s )
      s = sprintf('The sensor is below or exactly at the ground altitude.\nThe sensor is at %.3f km and the ground altitude is %.3f km.',(r_s-r_g)/1e3,z_bot/1e3);
      error(s);
    end      
 
    % Above the atmosphere?
    if( r_s>r_g+z_top )
      is_inside = 0;
    end

    % On the boundary looking above?
    if( r_s == r_g + z_top )
      if( dim == 1 )
        z_top_tilt = 0;
        znth_angle = psi_s;
      elseif( dim == 2 )
        i1         = find_grid_range( alpha_grid, alpha_s, 1 ); 
        r          = r_geoid(1,i1) + z_field(1,i1,n_p);
        c          = ( r_geoid(1,i1+1) + z_field(1,i1+1,n_p) - r ) / ...
                                         ( alpha_grid(i1+1) - alpha_grid(i1) );
        z_top_tilt = tilt_of_psurface( r, c );
        znth_angle = psi_s;
      else
        %
      end
      if( znth_angle<=90-z_top_tilt & znth_angle>=-90-z_top_tilt )
        is_inside = 0;
      end
    end

  end
