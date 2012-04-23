//! cart2poslos
/*! 
   The inverse of *poslos2cart*.

   The azimuth angle is set to: <br> 
      0 when the zenith angle is 0 or 180.
      atan2(dy,dx) at the poles (lat = +- 90).

   The longitude is set to 0 at the poles (lat = +- 90).

   \param   r     Out: Radius of observation position.
   \param   lat   Out: Latitude of observation position.
   \param   lon   Out: Longitude of observation position.
   \param   za    Out: LOS zenith angle at observation position.
   \param   aa    Out: LOS azimuth angle at observation position.
   \param   x     x-coordinate of observation position.
   \param   y     y-coordinate of observation position.
   \param   z     z-coordinate of observation position.
   \param   dx    x-part of LOS unit vector.
   \param   dy    y-part of LOS unit vector.
   \param   dz    z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void cart2poslos(
             double&   r,
             double&   lat,
             double&   lon,
             double&   za,
             double&   aa,
       const double&   x,
       const double&   y,
       const double&   z,
       const double&   dx,
       const double&   dy,
       const double&   dz )
{
  // Assert that LOS vector is normalised
  assert( abs( sqrt( dx*dx + dy*dy + dz*dz ) - 1 ) < 1e-6 );

  // Spherical coordinates
  cart2sph( r, lat, lon, x, y, z );

  // Spherical derivatives
  const double   coslat = cos( DEG2RAD * lat );
  const double   sinlat = sin( DEG2RAD * lat );
  const double   coslon = cos( DEG2RAD * lon );
  const double   sinlon = sin( DEG2RAD * lon );
  const double   dr   = coslat*coslon*dx    + coslat*sinlon*dy   + sinlat*dz;
  const double   dlat = -sinlat*coslon/r*dx - sinlat*sinlon/r*dy + coslat/r*dz;
  const double   dlon = -sinlon/coslat/r*dx + coslon/coslat/r*dy;

  // LOS angles
  za = RAD2DEG * acos( dr );
  //
  if( za < ANGTOL  ||  za > 180-ANGTOL  )
    { aa = 0; }

  else if( abs( lat ) <= POLELAT )
    {
      aa = RAD2DEG * acos( r * dlat / sin( DEG2RAD * za ) );

      if( isnan( aa ) )
        {
          if( dlat >= 0 )
            { aa = 0; }
          else
            { aa = 180; }
        }
      else if( dlon < 0 )
        { aa = -aa; }
    }

  // For lat = +- 90 the azimuth angle gives the longitude along which the
  // LOS goes
  else
    { aa = RAD2DEG * atan2( dy, dx ); }
}


//! cart2sph
/*! 
   The inverse of *sph2cart*.

   \param   r     Out: Radius of observation position.
   \param   lat   Out: Latitude of observation position.
   \param   lon   Out: Longitude of observation position.
   \param   x     x-coordinate of observation position.
   \param   y     y-coordinate of observation position.
   \param   z     z-coordinate of observation position.

   \author Patrick Eriksson
   \date   2002-12-30
*/
void cart2sph(
             double&    r,
             double&    lat,
             double&    lon,
       const double&    x,
       const double&    y,
       const double&    z )
{
  r   = sqrt( x*x + y*y + z*z );
  lat = RAD2DEG * asin( z / r );
  lon = RAD2DEG * atan2( y, x ); 
}



//! za_geom2other_point
/*!
   Calculates the zenith angle for the geometrical propagation path between
   two specified points.

   The returned zenith angle is valid for point 1. That is, the propagation
   path goes from point 1 to point 2.

   \return         Zenith angle.
   \param   r1     Radius for point 1.
   \param   lat1   Latiytude for point 1.
   \param   r2     Radius for point 2.
   \param   lat2   Latitude for point 2.

   \author Patrick Eriksson
   \date   2002-07-03
*/
Numeric za_geom2other_point(
       const Numeric&  r1,
       const Numeric&  lat1,
       const Numeric&  r2,
       const Numeric&  lat2 )
{
  if( lat2 == lat1 )
    {
      if( r1 <= r2 )
        { return 0; }
      else
        { return 180; }
    }
  else
    {
      // Absolute value of the latitude difference
      const Numeric dlat = abs( lat2 - lat1 );

      // The zenith angle is obtained by a combination of the lawes of sine
      // and cosine.
      Numeric za = dlat + RAD2DEG * asin( r1 * sin( DEG2RAD * dlat ) / 
                 sqrt( r1*r1 + r2*r2 - 2 * r1 * r2 * cos( DEG2RAD * dlat ) ) );

      // Consider viewing direction
      if( lat2 < lat1 )
        { za = -za; }

      return za;
    }
}



//! lat_crossing_3d
/*!
   Calculates where a 3D LOS crosses the specified latitude

   The solution algorithm is described in ATD. See the
   chapter on propagation paths.

   The function only looks for crossings in the forward direction of
   the given zenith angle (neglecting all solutions giving *l* <= 0).
   Note that the tangent point can be passed.
 
   R_NOT_FOUND, LON_NOT_FOUND and L_NOT_FOUND are returned if no solution 
   is found.

   \param   r         Out: Radius of found crossing.
   \param   lon       Out: Longitude of found crossing.
   \param   l         Out: Length along the path to the crossing.
   \param   lat_hit   Target latitude.
   \param   lat_start Latitude of start point.
   \param   aa_start  Azimuth angle at start point.
   \param   x         x-coordinate of start position.
   \param   y         y-coordinate of start position.
   \param   z         z-coordinate of start position.
   \param   dx        x-part of LOS unit vector.
   \param   dy        y-part of LOS unit vector.
   \param   dz        z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2012-02-29
*/
void lat_crossing_3d(
             Numeric&   r,
             Numeric&   lon,
             Numeric&   l,
       const Numeric&   lat_hit,
       const Numeric&   lat_start,
       const Numeric&   za_start,
       const Numeric&   x,
       const Numeric&   y,
       const Numeric&   z,
       const Numeric&   dx,
       const Numeric&   dy,
       const Numeric&   dz )
{
  assert( lat_start >= -90 );
  assert( lat_start <= 90 );
  assert( lat_hit >= -90 );
  assert( lat_hit <= 90 );

  // For za=0/180 and lat+-90 there is no solution
  if( za_start == 0  ||  za_start == 180  ||  abs(lat_hit) == 90 )
    { l = -1; }

  // The expressions below can not be used for lat=0
  else if( abs( lat_hit ) < 1e-7 )
    { l = -z / dz; }

  else
    {
      const Numeric   t2     = pow( tan( DEG2RAD * lat_hit ), 2.0 );
      const Numeric   a      = t2 * ( dx*dx + dy*dy ) - dz*dz;
      const Numeric   b      = 2 * ( t2 * ( x*dx + y*dy ) - z*dz );
      const Numeric   c      = t2 * ( x*x + y*y ) - z*z;
      const Numeric   bb     = b * b;
      const Numeric   ac4    = 4 * a * c;

      // Check if a real solution is possible
      if( ac4 > bb )
        { l = -1; }
      else
        {
          const Numeric   d      = -0.5*b/a;      
          const Numeric   e      = -0.5*sqrt(b*b-4*a*c)/a;      
                Numeric   l1     = d + e;
                Numeric   l2     = d - e;

          // Both lat and -lat can end up as a solution (the sign is lost as
          // tan(lat) is squared). A correct solution requires that l>=0 and
          // that z+l*dz has the same sigh as lat. Set l to -1 if this not
          // fulfilled.
          const Numeric zsign = sign( lat_hit ); 
          if( l1 > 0   &&   abs(sign(z+dz*l1)-zsign)>0.01 ) 
            { l1 = -1;}
          if( l2 > 0   &&   abs(sign(z+dz*l2)-zsign)>0.01 ) 
            { l2 = -1;}

          // If both l1 and l2 are > 0, we want theoretically the smallest
          // value. However, with lat=lat0 the "zero solution" can deviate
          // slightly from zero due to numerical issues, and it is not just to
          // pick the smallest positive value. As a solution, don't except a
          // final l below 1e-6 if not both l1 and l2 are inside [0,1e-6].
          const Numeric lmin = min( l1, l2 );
          const Numeric lmax = max( l1, l2 );
          if( lmin >= 0  && lmax < 1e-6 )
            { l = lmax; }
          else
            {
              if( lmin > 1e-6 ) 
                { l = lmin; }
              else if( lmax > 1e-6 )
                { l = lmax; }
              else
                { l = -1; }
            }
        }
    }

  if( l > 0 )
    {
      const Numeric xp = x+dx*l;
      const Numeric yp = y+dy*l;
      r   = sqrt( pow(xp,2.0) + pow(yp,2.0) + pow(z+dz*l,2.0) );
      lon = RAD2DEG * atan2( yp, xp );
    }
  else  
    { r = R_NOT_FOUND;   lon = LON_NOT_FOUND;   l   = L_NOT_FOUND; }
}



//! lon_crossing_3d
/*!
   Calculates where a 3D LOS crosses the specified longitude.

   The solution algorithm is described in ATD. See the
   chapter on propagation paths.

   The function only looks for crossings in the forward direction of
   the given zenith angle (neglecting all solutions giving *l* <= 0).
   Note that the tangent point can be passed.
 
   R_NOT_FOUND, LAT_NOT_FOUND and L_NOT_FOUND are returned if no solution 
   is found.

   \param   r         Out: Radius of found crossing.
   \param   lat       Out: Latitude of found crossing.
   \param   l         Out: Length along the path to the crossing.
   \param   lon_hit   Target longitude.
   \param   x         x-coordinate of start position.
   \param   y         y-coordinate of start position.
   \param   z         z-coordinate of start position.
   \param   dx        x-part of LOS unit vector.
   \param   dy        y-part of LOS unit vector.
   \param   dz        z-part of LOS unit vector.

   \author Patrick Eriksson
   \date   2012-02-29
*/
void lon_crossing_3d(
             Numeric&   r,
             Numeric&   lat,
             Numeric&   l,
       const Numeric&   lon_hit,
       const Numeric&   lon_start,
       const Numeric&   za_start,
       const Numeric&   aa_start,
       const Numeric&   x,
       const Numeric&   y,
       const Numeric&   z,
       const Numeric&   dx,
       const Numeric&   dy,
       const Numeric&   dz )
{

  if( lon_hit == lon_start  ||  za_start == 0  ||  za_start == 180  ||
                                aa_start == 0  ||  abs(aa_start) == 180 )
    { l = -1; }

  else
    {
      const double   tanlon = tan( DEG2RAD * lon_hit );
      l = ( y - x*tanlon ) / ( dx*tanlon - dy );
    }

  if( l <= 0 )
    { r = R_NOT_FOUND;   lat = LAT_NOT_FOUND;   l   = L_NOT_FOUND; }

  else
    {
      const Numeric zp = z + dz*l;
      r   = sqrt( pow(x+dx*l,2.0) + pow(y+dy*l,2.0) + pow(zp,2.0) );
      lat = RAD2DEG * asin( zp / r );
    }
}



//! plevel_crossing_3d
/*!
   As plevel_crossing_2d but for 3D
 

   \param   r           Out: Radius at crossing.
   \param   lat         Out: Latitude at crossing.
   \param   lon         Out: Longitude at crossing.
   \param   l           Out: Distance between start and crossing points.
   \param   r_start0    In: Radius of start point.
   \param   lat_start   In: Latitude of start point.
   \param   lon_start   In: Longitude of start point.
   \param   za_start    In: LOS zenith angle at start point.
   \param   aa_start    In: LOS azimuth angle at start point.
   \param   ppc         In: Propagation path constant.
   \param   lat1        In: Latitude of lower end.
   \param   lat3        In: Latitude of upper end.
   \param   lon5        In: Longitude of lower end.
   \param   lon6        In: Longitude of upper end.
   \param   r15         In: Radius at lat1/lon5.
   \param   r35         In: Radius at lat3/lon5.
   \param   r36         In: Radius at lat3/lon6.
   \param   r16         In: Radius at lat1/lon6.
   \param   x           In: x-coordinate of start position.
   \param   y           In: y-coordinate of start position.
   \param   z           In: z-coordinate of start position.
   \param   dx          In: x-part of LOS unit vector.
   \param   dy          In: y-part of LOS unit vector.
   \param   dz          In: z-part of LOS unit vector.
   \param   above       In: True if ppath start point is above level. 
                        In: Otherwise false.

   \author Patrick Eriksson
   \date   2012-02-19
*/
void plevel_crossing_3d(
              Numeric&  r,
              Numeric&  lat,
              Numeric&  lon,
              Numeric&  l,
        const Numeric&  r_start0,
        const Numeric&  lat_start,
        const Numeric&  lon_start,
        const Numeric&  za_start,
        const Numeric&  aa_start,
        const Numeric&  x,
        const Numeric&  y,
        const Numeric&  z,
        const Numeric&  dx,
        const Numeric&  dy,
        const Numeric&  dz,
        const Numeric&  ppc,
        const Numeric&  lat1,
        const Numeric&  lat3,
        const Numeric&  lon5,
        const Numeric&  lon6,
        const Numeric&  r15,
        const Numeric&  r35,
        const Numeric&  r36,
        const Numeric&  r16,
        const bool&     above )
{
  assert( za_start <= 180 );
  assert( lat_start >=lat1  &&  lat_start <= lat3 );
  assert( lon_start >=lon5  &&  lon_start <= lon6 );

  const Numeric rmin = min( r15, min( r35, min( r36, r16 ) ) );
  const Numeric rmax = max( r15, max( r35, max( r36, r16 ) ) );

  // The case of negligible slope
  if( rmax-rmin < RTOL/10 )
    {
      // Set r_start, considering impact of numerical problems
      Numeric r_start = r_start0;
              r       = r15;
      if( above )
        { if( r_start < rmax ) { r_start = r = rmax; } }
      else
        { if( r_start > rmin ) { r_start = r = rmin; } }

      r_crossing_3d( lat, lon, l, r, r_start, za_start, ppc, 
                                                         x, y, z, dx, dy, dz );

      // Check if inside [lat1,lat3]
      if( lat > lat3  ||  lat < lat1  || lon > lon6  ||  lon < lon5 )
        { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND;   lon = LON_NOT_FOUND; }  
    }

  // With slope
  else
    {
      // Set r_start, considering impact of numerical problems
      Numeric r_start = r_start0;
      if( above )
        { if( r_start < rmin ) { r_start = rmin; } }
      else
        { if( r_start > rmax ) { r_start = rmax; } }

      // Calculate crossing with closest radius
      if( r_start > rmax )
        {
          r = rmax;
          r_crossing_3d( lat, lon, l, r, r_start, za_start, ppc,
                                                         x, y, z, dx, dy, dz );
        }
      else if( r_start < rmin )
        {
          r = rmin;      
          r_crossing_3d( lat, lon, l, r, r_start, za_start, ppc,
                                                         x, y, z, dx, dy, dz );
        }
      else
        { r = r_start; lat = lat_start; lon = lon_start; l = 0; }
  
      // lat/lon must be inside [lat1,lat3]/[lon5,lon6]] if relevant to continue
      if( lat < lat1  ||  lat > lat3  ||  lon < lon5  ||  lon > lon6 )
        { r = R_NOT_FOUND; }   // lat and lon already set by r_crossing_3d

      // Otherwise continue from found point, considering the level slope 
      else
        {
          // Level radius at lat/lon
          const Numeric  rpl = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                                r15, r35, r36, r16, lat, lon );
          
          // Make adjustment if numerical problems
          if( above )
            { if( r < rpl ) { r = rpl; } }
          else
            { if( r > rpl ) { r = rpl; } }

          // za_start = 0
          if( za_start < ANGTOL )
            {
              if( r >= rpl )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND;   
                  l = L_NOT_FOUND;  lon = LON_NOT_FOUND; }
              else
                { r = rpl; lat = lat_start; lon = lon_start; l = rpl-r_start; }
            }

          // za_start = 180
          else if( za_start > 180-ANGTOL )
            {
              if( r <= rpl )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND;   
                  l = L_NOT_FOUND;  lon = LON_NOT_FOUND; }
              else
                { r = rpl; lat = lat_start; lon = lon_start; l = r_start-rpl; }
            }

          else
            {
              // Azimuth angle:
              Numeric d1, d2, d3, za, aa;
              if( l == 0 )
                { za = za_start; aa = aa_start; }
              else
                { cart2poslos( d1, d2, d3, za, aa, x+dx*l, y+dy*l, z+dz*l, 
                                                                dx, dy, dz ); 
                  assert( abs(d1-r) < 1e-3 );
                  assert( abs(d2-lat) < 1e-8 );
                  assert( abs(d3-lon) < 1e-8 );
                }

              // Level slope at lat/lon
              const Numeric  cpl = plevel_slope_3d( lat1, lat3, lon5, lon6, 
                                            r15, r35, r36, r16, lat, lon, aa );
              // Angular distance from present point to actual crossing
              const Numeric dang = rslope_crossing( r, za, rpl, cpl );

              // Lat and lon at dang
              const Numeric danrad = DEG2RAD * dang;
              const Numeric latrad = DEG2RAD * lat;
              const Numeric aarad  = DEG2RAD * aa;
              const Numeric cosdan = cos( danrad );
              const Numeric sindan = sin( danrad );
              const Numeric coslat = cos( latrad );
              const Numeric sinlat = sin( latrad );
              //
              lat  = RAD2DEG*asin( sinlat*cosdan + coslat*sindan*cos(aarad) );
              lon += RAD2DEG*atan2( sin(aarad)*sindan*coslat,
                                    cosdan-sinlat*sin(DEG2RAD*lat) );

              // lat/lon still inside gridbox? If yes, update r and l
              if( lat < lat1  ||  lat > lat3  ||  lon < lon5  ||  lon > lon6 )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND;   
                  l = L_NOT_FOUND;  lon = LON_NOT_FOUND; }
              else
                {
                  r = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                                r15, r35, r36, r16, lat, lon );
                  l = abs( geompath_l_at_r( ppc, r_start ) -
                           geompath_l_at_r( ppc, r ) );
                }
            }          
        }
    }
}







//! refr_gradients_1d
/*! 
   Determines the refractive index, and the radial gradient.

   The gradient is calculated in pure numerical way. That is, the
   refractive index is calculated for slightly shifted radius and the
   difference to the refractive index at the given point determines
   the gradient.

   The atmosphere is given by its 1D view. That is, the latitude and
   longitude dimensions are removed from the atmospheric fields. For
   example, the temperature is given as a vector.

   \param   ws                  Current Workspace
   \param   refr_index          Output: As the WSV with the same name.
   \param   dndr                Output: Radial gradient of refractive index.
   \param   a_pressure          Pressure in hPa.
   \param   a_temperature       Temperature in K.
   \param   a_vmr_list          Vector with VMR values.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   p_grid              As the WSV with the same name.
   \param   refellipsoid        As the WSV with the same name.
   \param   z_field             The geometric altitude of each pressure level
                                at each latitude.
   \param   t_field             The temperature 2D field.
   \param   vmr_field           The VMR 2D field for each species.
   \param   r                   The radius of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-20
*/
void refr_gradients_1d(
              Workspace&  ws,
              Numeric&    refr_index,
              Numeric&    dndr,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   refellipsoid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        const Numeric&    r )
{ 
   get_refr_index_1d( ws, refr_index, a_pressure,  a_temperature, a_vmr_list, 
                      refr_index_agenda, p_grid, 
                      refellipsoid, z_field, t_field, vmr_field, r );

   const Numeric   n0 = refr_index;

   get_refr_index_1d( ws, refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, p_grid, 
                      refellipsoid, z_field, t_field, vmr_field, r+1 );

   dndr = refr_index - n0;

   refr_index = n0;
}



void ppath_geom_updown_1d(
              Ppath&      ppath,
        ConstVectorView   z_field,
        ConstVectorView   refellipsoid,
        const Numeric&    z_surface,
        const Index&      cloudbox_on, 
     const ArrayOfIndex&  cloudbox_limits )
{
  // Starting radius, zenith angle and latitude
  Numeric r_start, lat_start, za_start;

  // Index of the pressure level being the lower limit for the
  // grid range of interest.
  Index ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( r_start, lat_start, za_start, ip, ppath );

  if( za_start > 85  &&  za_start < 120 )
    {
      throw runtime_error( "This method can not be used for initial zenith "
                           "angles between 85 and 120 deg.!!");
    }

  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  Numeric ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }

  // Upward
  if( za_start < 90 )
    { 
      // Determine number of ppath points
      Index ilastp1;  // Last index + 1
      if( cloudbox_on  &&  cloudbox_limits[0] > ip )
        { ilastp1 = cloudbox_limits[0] + 1; }  // Points inside cloudbox
      else                                     // are handled by 
        { ilastp1 = z_field.nelem(); }         // ppath_start_stepping
      const Index np = ilastp1 - ip;

      ppath_init_structure( ppath, 1, np );
      //
      ppath.constant = ppc;
      //
      // Start point
      ppath.r[0]          = r_start;
      ppath.pos(0,0)      = r_start - refellipsoid[0];
      ppath.los(0,0)      = za_start;
      ppath.pos(0,1)      = lat_start;
      Numeric llast        = geompath_l_at_r( ppc, ppath.r[0] );
      ppath.gp_p[0].idx   = ip;
      ppath.gp_p[0].fd[0] = ( ppath.pos(0,0) - z_field[ip] ) / 
                            ( z_field[ip+1]  - z_field[ip] );
      ppath.gp_p[0].fd[1] = 1 - ppath.gp_p[0].fd[0];
      gridpos_check_fd( ppath.gp_p[0] );
      // Later points
      for( Index i=1; i<np; i++ )
        {
          ppath.pos(i,0)      = z_field[ip+i];
          ppath.r[i]          = refellipsoid[0] + ppath.pos(i,0);
          ppath.los(i,0)      = geompath_za_at_r( ppc, za_start, ppath.r[i] );
          ppath.pos(i,1)      = geompath_lat_at_za( za_start, lat_start,
                                                              ppath.los(i,0) );
          const Numeric lthis  = geompath_l_at_r( ppc, ppath.r[i] );
          ppath.lstep[i-1]   = lthis - llast;
          llast               = lthis;
          ppath.gp_p[i].idx   = ip + i;
          ppath.gp_p[i].fd[0] = 0;
          ppath.gp_p[i].fd[1] = 1;
        }
      // Special treatment of last point
      ppath.gp_p[np-1].idx  -= 1;
      ppath.gp_p[np-1].fd[0] = 1;
      ppath.gp_p[np-1].fd[1] = 0;
    }

  // Downward
  else
    {
      if( ppc > refellipsoid[0] + z_surface )
        {
          ostringstream os;
          os << "This function can not be used for propgation paths\n"
             << "including tangent points. Such a point occurs in this case.";
          throw runtime_error(os.str());
        }

      // Find grid position of surface altitude
      GridPos   gp;
      gridpos( gp, z_field, z_surface );
      
      // Determine number of ppath points. Start assumption is hit with surface
      Index  ilast   = gp.idx+1;  // Index of last pressure level to include
      Index  surface = 1;
      if( cloudbox_on  &&  cloudbox_limits[1] <= ip )
        {                                  // Points inside cloudbox are
          ilast   = cloudbox_limits[1];    // handled by ppath_start_stepping
          surface = 0;
        }  
      const Index np = ip - ilast + 2 + surface;

      ppath_init_structure( ppath, 1, np );
      //
      ppath_set_background( ppath, 2 );
      ppath.constant = ppc;
      //
      // Start point
      ppath.r[0]          = r_start;
      ppath.pos(0,0)      = r_start - refellipsoid[0];
      ppath.los(0,0)      = za_start;
      ppath.pos(0,1)      = lat_start;
      Numeric llast        = geompath_l_at_r( ppc, ppath.r[0] );
      ppath.gp_p[0].idx   = ip;
      ppath.gp_p[0].fd[0] = ( ppath.pos(0,0) - z_field[ip] ) / 
                            ( z_field[ip+1]  - z_field[ip] );
      ppath.gp_p[0].fd[1] = 1 - ppath.gp_p[0].fd[0];
      gridpos_check_fd( ppath.gp_p[0] );
      // Later points
      for( Index i=1; i<np; i++ )
        {
          if( i < np-1  ||  !surface )
            { ppath.pos(i,0)  = z_field[ip-i+1]; }
          else
            { ppath.pos(i,0)  = z_surface; }
          ppath.r[i]          = refellipsoid[0] + ppath.pos(i,0);
          ppath.los(i,0)      = geompath_za_at_r( ppc, za_start, ppath.r[i] );
          ppath.pos(i,1)      = geompath_lat_at_za( za_start, lat_start,
                                                              ppath.los(i,0) );
          const Numeric lthis  = geompath_l_at_r( ppc, ppath.r[i] );
          ppath.lstep[i-1]   = llast - lthis;
          llast               = lthis;
          ppath.gp_p[i].idx   = ip - i + 1;
          ppath.gp_p[i].fd[0] = 0;
          ppath.gp_p[i].fd[1] = 1;
        }
      // Special treatment of last point
      if( surface )
        { gridpos_copy( ppath.gp_p[np-1], gp ); }
    }
}



