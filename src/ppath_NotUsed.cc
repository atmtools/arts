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

      if( std::isnan( aa ) )
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
          // pick the smallest positive value. As a solution, don't accept a
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
   \param   lon_start Longitude of start position.
   \param   za_start  Zenith angle at start position.
   \param   aa_start  Azimuth angle at start position.
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

  // Zenith case
  if( za_start < ANGTOL )
    {
      if( above )
        { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND; 
          l = L_NOT_FOUND;  lon = LON_NOT_FOUND; } 
      else
        {
          lat = lat_start;
          lon = lon_start;
          r   = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15, r35, r36, r16,
                                                                    lat, lon );
          l   = max( 1e-9, r - r_start0 ); // Max to ensure a small positive
        }                                  // step, to handle numerical issues
    }

  // Nadir case
  else if( za_start > 180-ANGTOL )
    {
      if( above )
        {
          lat = lat_start;
          lon = lon_start;
          r   = rsurf_at_latlon( lat1, lat3, lon5, lon6, r15, r35, r36, r16,
                                                                    lat, lon );
          l   = max( 1e-9, r_start0 - r ); // As above
        }  
      else
        { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND; 
          l = L_NOT_FOUND;  lon = LON_NOT_FOUND; } 
    }

  // The general case
  else
    {
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

          r_crossing_3d( lat, lon, l, r, r_start, lat_start, lon_start,
                         za_start, ppc, x, y, z, dx, dy, dz );

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
              r_crossing_3d( lat, lon, l, r, r_start, lat_start, lon_start,
                             za_start, ppc, x, y, z, dx, dy, dz );
            }
          else if( r_start < rmin )
            {
              r = rmin;      
              r_crossing_3d( lat, lon, l, r, r_start, lat_start, lon_start,
                             za_start, ppc, x, y, z, dx, dy, dz );
            }
          else
            { r = r_start; lat = lat_start; lon = lon_start; l = 0; }
  
          // lat/lon must be inside if relevant to continue
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

              // Azimuth angle:
              Numeric d1, d2, d3, za, aa;
              if( l == 0 )
                { za = za_start; aa = aa_start; }
              else
                { cart2poslos( d1, d2, d3, za, aa, x+dx*l, y+dy*l, z+dz*l, 
                               dx, dy, dz, ppc, lat_start, lon_start, 
                               za_start, aa_start ); 
                  assert( abs(d1-r) < 1e-3 );
                  assert( abs(d2-lat) < 1e-8 );
                  assert( abs(d3-lon) < 1e-8 );
                }

              // Level slope at lat/lon
              Numeric  c1, c2; 
              plevel_slope_3d( c1, c2, lat1, lat3, lon5, lon6, 
                               r15, r35, r36, r16, lat, lon, aa );

              // Angular distance from present point to actual crossing
              const Numeric dang = rslope_crossing3d( r, za, rpl, c1, c2 );

              // Lat and lon at dang
              Numeric lat2, lon2;
              //
              latlon_at_aa( lat2, lon2, lat, lon, aa, dang ); 
              //
              lat = lat2;
              lon = lon2;   resolve_lon( lon, lon5, lon6 );

              // lat/lon still inside gridbox? If yes, update r and l
              if( lat < lat1  ||  lat > lat3  ||  lon < lon5  ||  lon > lon6 )
                { r = R_NOT_FOUND;  lat = LAT_NOT_FOUND;   
                  l = L_NOT_FOUND;  lon = LON_NOT_FOUND; }
              else
                {
                  r = rsurf_at_latlon( lat1, lat3, lon5, lon6, 
                                                r15, r35, r36, r16, lat, lon );
                  distance3D( l, r_start, lat_start, lon_start, r, lat, lon );
                }
            }          
        }
    }
}



//! do_gridcell_3d
/*!
   Calculates the geometrical path through a 3D grid cell.

   The function determines the geometrical path from the given start
   point to the boundary of the grid cell. The face where the path
   exits the grid cell is denoted as the end face. The same number
   coding as in *do_gridcell_2d* is used, where the additional longitude
   end faces are numbered as: <br>
   5: The face at the lower longitude point. <br>
   6: The face at the upper longitude point.

   The corner points are numbered as *do_gridcell_2d*, but 5 or 6 is added
   after the latitude number to indicate the longitude. This means that
   r16a, is the corner at lat1, lon6 and pressure level a.

   See further *do_gridcell_2d*.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   lon_v       Out: Vector with longitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   aa_v        Out: Vector with LOS azimuth angle at found path points.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   r_start0    Radius of start point.
   \param   lat_start0  Latitude of start point.
   \param   lon_start0  Longitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   aa_start    LOS azimuth angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   lat1        Latitude of left end face (face 1) of the grid cell.
   \param   lat3        Latitude of right end face (face 3) of the grid cell.
   \param   lon5        Lower longitude limit of the grid cell.
   \param   lon6        Upper longitude limit of the grid cell.
   \param   r15a        Radius of corner: lower p-level,*lat1* and *lon5*.
   \param   r35a        Radius of corner: lower p-level,*lat3* and *lon5*.
   \param   r36a        Radius of corner: lower p-level,*lat3* and *lon6*.
   \param   r16a        Radius of corner: lower p-level,*lat1* and *lon6*.
   \param   r15b        Radius of corner: upper p-level,*lat1* and *lon5*.
   \param   r35b        Radius of corner: upper p-level,*lat3* and *lon5*.
   \param   r36b        Radius of corner: upper p-level,*lat3* and *lon6*.
   \param   r16b        Radius of corner: upper p-level,*lat1* and *lon6*.
   \param   rsurface15  Radius for the surface at *lat1* and *lon5*.
   \param   rsurface35  Radius for the surface at *lat3* and *lon5*.
   \param   rsurface36  Radius for the surface at *lat3* and *lon6*.
   \param   rsurface16  Radius for the surface at *lat1* and *lon6*.

   \author Patrick Eriksson
   \date   2002-11-28
*/
void do_gridcell_3d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   lon_v,
              Vector&   za_v,
              Vector&   aa_v,
              Numeric&  lstep,
              Index&    endface,
        const Numeric&  r_start, 
        const Numeric&  lat_start,
        const Numeric&  lon_start,
        const Numeric&  za_start,
        const Numeric&  aa_start,
        const Numeric&  ppc,
        const Numeric&  lmax,
        const Numeric&  lat1,
        const Numeric&  lat3,
        const Numeric&  lon5,
        const Numeric&  lon6,
        const Numeric&  r15a,
        const Numeric&  r35a,
        const Numeric&  r36a,
        const Numeric&  r16a,
        const Numeric&  r15b,
        const Numeric&  r35b,
        const Numeric&  r36b,
        const Numeric&  r16b,
        const Numeric&  rsurface15,
        const Numeric&  rsurface35,
        const Numeric&  rsurface36,
        const Numeric&  rsurface16 )
{
  // Radius end latitude of end point
  Numeric r, lat, lon, l= L_NOT_FOUND;  // l not always calculated/needed

  endface = 0;

  // Sensor pos and LOS in cartesian coordinates
  Numeric   x, y, z, dx, dy, dz;
  poslos2cart( x, y, z, dx, dy, dz, r_start, lat_start, lon_start, 
                                             za_start,  aa_start ); 

  // Check if crossing with lower pressure level
  plevel_crossing_3d( r, lat, lon, l, r_start, lat_start, lon_start,
                      za_start, aa_start, x, y, z, dx, dy, dz, ppc, 
                      lat1, lat3, lon5, lon6, r15a, r35a, r36a, r16a, true );
  if( r > 0 )
    { endface = 2; }

  // Check if crossing with surface
  if( rsurface15 >= r15a  ||  rsurface35 >= r35a  ||
      rsurface36 >= r36a  ||  rsurface16 >= r16a )
    {
      Numeric rt, latt, lont, lt; 
      plevel_crossing_3d( rt, latt, lont, lt, r_start, lat_start, lon_start,
                          za_start, aa_start, x, y, z, dx, dy, dz, ppc, 
                          lat1, lat3, lon5, lon6, rsurface15, rsurface35, 
                          rsurface36, rsurface16, true );

      if( rt > 0  &&  lt <= l )  // lt<=l to resolve the closest crossing
        { endface = 7;   r = rt;   lat = latt;   lon = lont;   l = lt; }
    }


  // Upper pressure level
  {
    Numeric rt, latt, lont, lt;     
    plevel_crossing_3d( rt, latt, lont, lt, r_start, lat_start, lon_start, 
                          za_start, aa_start, x, y, z, dx, dy, dz, ppc, 
                          lat1, lat3, lon5, lon6, r15b, r35b, r36b, r16b, 
                          false ); 
    if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
      { endface = 4;   r = rt;   lat = latt;   lon = lont;   l = lt; }
  }

  // Latitude 1
  {
    Numeric rt, lont, lt;     
    lat_crossing_3d( rt, lont, lt, lat1, lat_start, za_start,
                                                         x, y, z, dx, dy, dz );
    if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
      { endface = 1;   r = rt;   lat = lat1;   lon = lont;   l = lt; }
  }

  // Latitude 3
  {
    Numeric rt, lont, lt;     
    lat_crossing_3d( rt, lont, lt, lat3, lat_start, za_start,
                                                         x, y, z, dx, dy, dz );
    if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
      { endface = 3;   r = rt;   lat = lat3;   lon = lont;   l = lt; }
  }

  // Longitude 5 (only done if solution is lacking)
  if( !( r>0 && lat>=lat1 && lat<=lat3 && lon>=lon5 && lon<=lon6 ) )
    {
      Numeric rt, latt, lt;     
      lon_crossing_3d( rt, latt, lt, lon5, lon_start, za_start, aa_start,
                                                         x, y, z, dx, dy, dz );
      if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
        { endface = 5;   r = rt;   lat = latt;   lon = lon5;   l = lt; }
    }

  // Longitude 6 (only done if solution is lacking)
  if( !( r>0 && lat>=lat1 && lat<=lat3 && lon>=lon5 && lon<=lon6 ) )
    {
      Numeric rt, latt, lt;     
      lon_crossing_3d( rt, latt, lt, lon6, lon_start, za_start, aa_start,
                                                         x, y, z, dx, dy, dz );
      if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
        { endface = 6;   r = rt;   lat = latt;   lon = lon6;   l = lt; }
    }

  assert( endface );

  // Check if there is a tangent point inside the grid cell. 
  if( za_start > 90  )
    {
      Numeric ltan = geompath_l_at_r( ppc, r_start );
      if( l-ltan > LACC ) 
        { 
          endface = 8; 
          geompath_tanpos_3d( r, lat, lon, l, r_start, lat_start, lon_start, 
                                                     za_start, aa_start, ppc );
        }
    }

  resolve_lon( lon, lon5, lon6 );              


  //--- Create return vectors
  //
  Index n = 1;
  //
  if( lmax > 0 )
    {
      n = Index( ceil( abs( l / lmax ) ) );
      if( n < 1 )
        { n = 1; }
    }
  //
  r_v.resize( n+1 );
  lat_v.resize( n+1 );
  lon_v.resize( n+1 );
  za_v.resize( n+1 );
  aa_v.resize( n+1 );
  //
  r_v[0]   = r_start;
  lat_v[0] = lat_start;
  lon_v[0] = lon_start;
  za_v[0]  = za_start;
  aa_v[0]  = aa_start;
  //
  lstep = l / (Numeric)n;
  // 
  for( Index j=1; j<=n; j++ )
    {
      const Numeric lj  = lstep * (Numeric)j;

      cart2poslos( r_v[j], lat_v[j], lon_v[j], za_v[j], aa_v[j],
                   x+dx*lj, y+dy*lj, z+dz*lj, dx, dy, dz, ppc,
                   lat_start, lon_start, za_start, aa_start );
      resolve_lon( lon_v[j], lon5, lon6 );
   }

  //--- Set last point especially, which should improve the accuracy
  r_v[n]   = r; 
  lat_v[n] = lat;
  lon_v[n] = lon;
  if( endface == 8 )
    { za_v[n] = 90; }
  else
    { za_v[n] = geompath_za_at_r( ppc, za_start, r_v[n] ); }
}






//! do_gridcell_2d
/*!
   Calculates the geometrical path through a 2D grid cell.

   The function determines the geometrical path from the given start
   point to the boundary of the grid cell. The face where the path
   exits the grid cell is denoted as the end face. The following
   number coding is used for the variable *endface*: <br>
   1: The face at the lower latitude point. <br>
   2: The face at the lower (geometrically) pressure level. <br>
   3: The face at the upper latitude point. <br>
   4: The face at the upper (geometrically) pressure level. <br>
   7: The end point is an intersection with the surface. 

   The corner points are names r[lat][a,b]. For example: r3b.
   The latitudes are numbered to match the end faces. This means that
   the lower latitude has number 1, and the upper number 3. The pressure
   levels are named as a and b: <br>
   a: Lower pressure level (highest pressure). <br>
   b: Upper pressure level (lowest pressure).

   Path points are included if *lmax*>0 and the distance to the end
   point is > than *lmax*.

   The return vectors (*r_v* etc.) can have any length when handed to
   the function.

   \param   r_v         Out: Vector with radius of found path points.
   \param   lat_v       Out: Vector with latitude of found path points.
   \param   za_v        Out: Vector with LOS zenith angle at found path points.
   \param   lstep       Out: Vector with length along the path between points.
   \param   endface     Out: Number coding for exit face. See above.
   \param   r_start0    Radius of start point.
   \param   lat_start   Latitude of start point.
   \param   za_start    LOS zenith angle at start point.
   \param   ppc         Propagation path constant.
   \param   lmax        Maximum allowed length along the path. -1 = no limit.
   \param   lat1        Latitude of left end face (face 1) of the grid cell.
   \param   lat3        Latitude of right end face (face 3) of the grid cell.
   \param   r1a         Radius of lower-left corner of the grid cell.
   \param   r3a         Radius of lower-right corner of the grid cell.
   \param   r3b         Radius of upper-right corner of the grid cell (r3b>r3a)
   \param   r1b         Radius of upper-left corner of the grid cell (r1b>r1a).
   \param   rsurface1   Radius for the surface at *lat1*.
   \param   rsurface3   Radius for the surface at *lat3*.

   \author Patrick Eriksson
   \date   2002-11-28
*/
void do_gridcell_2d(
              Vector&   r_v,
              Vector&   lat_v,
              Vector&   za_v,
              Numeric&  lstep,
              Index&    endface,
        const Numeric&  r_start,
        const Numeric&  lat_start,
        const Numeric&  za_start,
        const Numeric&  ppc,
        const Numeric&  lmax,
        const Numeric&  lat1,
        const Numeric&  lat3,
        const Numeric&  r1a,
        const Numeric&  r3a,
        const Numeric&  r3b,
        const Numeric&  r1b,
        const Numeric&  rsurface1,
        const Numeric&  rsurface3 )
{
  // Radius and latitude of end point, and the length to it
  Numeric r, lat, l= L_NOT_FOUND;  // l not always calculated/needed

  endface = 0;

  // Check if crossing with lower pressure level
  plevel_crossing_2d( r, lat, l, r_start, lat_start, za_start, ppc, lat1, lat3, 
                                                              r1a, r3a, true );
  if( r > 0 )
    { endface = 2; }

  // Check if crossing with surface
  if( rsurface1 >= r1a  ||  rsurface3 >= r3a )
    {
      Numeric rt, latt, lt; 
      plevel_crossing_2d( rt, latt, lt, r_start, lat_start, za_start, ppc, 
                                    lat1, lat3, rsurface1, rsurface3, true );

      if( rt > 0  &&  lt <= l )  // lt<=l to resolve the closest crossing
        { endface = 7;   r = rt;   lat = latt;   l = lt; }
    }

  // Upper pressure level
  {
    Numeric rt, latt, lt; 
    plevel_crossing_2d( rt, latt, lt, r_start, lat_start, za_start, ppc, 
                                                 lat1, lat3, r1b, r3b, false );
    if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
      { endface = 4;   r = rt;   lat = latt;  /* l = lt; */ }
  }
  
  // Latitude endfaces
  if( r <= 0 )
    {
      if( za_start < 0 )
        { endface = 1;  lat = lat1; }
      else
        { endface = 3;  lat = lat3; }
      r = geompath_r_at_lat( ppc, lat_start, za_start, lat ); 
    }

  assert( endface );

  //Check if there is a tangent point inside the grid cell. 
  bool tanpoint;
  const Numeric absza = abs( za_start );
  if( absza > 90  &&  ( absza - abs(lat_start-lat) ) < 90 ) 
    { tanpoint = true; }
  else
    { tanpoint = false; }

  geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                          za_start, r, tanpoint, lmax );

  // Force exact values for end point when they are known
  if( endface == 1  ||  endface == 3 )
    { lat_v[lat_v.nelem()-1] = lat; }
}



//! raytrace_1d_linear_basic
/*! 
   Performs ray tracing for 1D with linear steps.

   A geometrical step with length of *lraytrace* is taken from each
   point. The zenith angle for the end point of that step is
   calculated exactly by the expression c = r*n*sin(theta), and a new
   step is taken. The length of the last ray tracing step to reach the
   end radius is adopted to the distance to the end radius.

   The refractive index is assumed to vary linearly between the pressure
   levels.

   As the ray tracing is performed from the last end point, the found path
   will not be symmetric around the tangent point.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of ATD.

   \param   ws           Current Workspace
   \param   r_array      Out: Radius of ray tracing points.
   \param   lat_array    Out: Latitude of ray tracing points.
   \param   za_array     Out: LOS zenith angle at ray tracing points.
   \param   l_array      Out: Distance along the path between ray tracing 
                         points.
   \param   n_array      Out: Refractive index at ray tracing points.
   \param   ng_array     Out: Group refractive index at ray tracing points.
   \param   endface      See do_gridrange_1d.
   \param   refellipsoid As the WSV with the same name.
   \param   p_grid       Pressure grid.
   \param   z_field      Geometrical altitudes (1D).
   \param   t_field      As the WSV with the same name.
   \param   vmr_field    As the WSV with the same name.
   \param   f_grid       As the WSV with the same name.
   \param   lmax         As the WSV ppath_lmax
   \param   refr_index_air_agenda   As the WSV with the same name.
   \param   lraytrace    Maximum allowed length for ray tracing steps.
   \param   ppc          Propagation path constant.
   \param   r_surface    Radius of the surface.
   \param   r1           Radius of lower pressure level.
   \param   r3           Radius of upper pressure level (r3 > r1).
   \param   r            Start radius for ray tracing.
   \param   lat          Start latitude for ray tracing.
   \param   za           Start zenith angle for ray tracing.

   \author Patrick Eriksson
   \date   2002-12-02
*/
void raytrace_1d_linear_basic(
              Workspace&       ws,
              Array<Numeric>&  r_array,
              Array<Numeric>&  lat_array,
              Array<Numeric>&  za_array,
              Array<Numeric>&  l_array,
              Array<Numeric>&  n_array,
              Array<Numeric>&  ng_array,
              Index&           endface,
        ConstVectorView        refellipsoid,
        ConstVectorView        p_grid,
        ConstVectorView        z_field,
        ConstTensor3View       t_field,
        ConstTensor4View       vmr_field,
        ConstVectorView        f_grid,
        const Numeric&         lmax,
        const Agenda&          refr_index_air_agenda,
        const Numeric&         lraytrace,
        const Numeric&         ppc,
        const Numeric&         r_surface,
        const Numeric&         r1,
        const Numeric&         r3,
              Numeric&         r,
              Numeric&         lat,
              Numeric&         za )
{
//  /*
  // We currently do not handle bending of rays back into the medium due to
  // refraction. Following check test, whether this is expected to happen for
  // the given grid cell and ray zenith angle.
  // Generally constant path constant rule applies:
  // n1 * r1 * sin(th1) == n3 * r3 * sin(th3),
  // where largest path constant values that can occur are (n1*r1) and (n3*r3)
  // No solution exists if
  // n1 * r1 * sin(th1) > n3 * r3, hence
  // n3 < n1 * sin(th1) * r1/r3 (similar for n1 < n3 * sin(th3) * r3/r1)
  // Assuming that za is given at r1, we check for n3 >= n1 * sin(za) * r1/r3
  Numeric refr1, refr3, refrg;
  get_refr_index_1d( ws, refr1, refrg, 
                     refr_index_air_agenda, p_grid, refellipsoid, z_field,
                     t_field, vmr_field, f_grid, r1 );
  get_refr_index_1d( ws, refr3, refrg, 
                     refr_index_air_agenda, p_grid, refellipsoid, z_field,
                     t_field, vmr_field, f_grid, r3 );
  if( refr3 < refr1 * (r1/r3) * sin( DEG2RAD * abs(za) ) )
    {
      ostringstream os;
      os << "For path between r1=" << r1 << "(n-1=" << refr1-1. << ") and r2="
         << r3 << "(n-1=" << refr3-1. << "),\n"
         << "path calculation will run into refractive back-bending issues.\n"
         << "We are currently NOT ABLE to handle them. Consider repeating\n"
         << "your calculation without taking refraction into account.";
      throw runtime_error(os.str());
    }
//  */

  // Loop boolean
  bool ready = false;

  // Store first point
  Numeric refr_index_air, refr_index_air_group;
  get_refr_index_1d( ws, refr_index_air, refr_index_air_group, 
                     refr_index_air_agenda, p_grid, refellipsoid, z_field,
                     t_field, vmr_field, f_grid, r );
  r_array.push_back( r );
  lat_array.push_back( lat );
  za_array.push_back( za );
  n_array.push_back( refr_index_air );
  ng_array.push_back( refr_index_air_group );

  // Variables for output from do_gridrange_1d
  Vector    r_v, lat_v, za_v;
  Numeric   lstep, lcum = 0;

  while( !ready )
    {
      // Constant for the geometrical step to make
      const Numeric   ppc_step = geometrical_ppc( r, za );

      // Explicitly check here, that predicted path point is still within
      // gridcell and did not undetectedly slip out. This can happen in case of
      // close-to-lateral angles and high refraction gradient, which causes a
      // bending of the ray back into the medium (kind of reflection). There is
      // no proper handling of this back-bending case yet.
      if( ( r < r1 - RTOL ) || ( r > r3 + RTOL ) )
        {
          throw runtime_error(
            "Ooops. Path undetectedly left the grid cell.\n"
            "This should not happen. But there are issues with cases of high\n"
            "refractivity gradients. Seems you unfortunately encountered such\n"
            "a case. Little to be done about this now (if this is an option\n"
            "for you, run the case without considering refraction). For further\n"
            "details, contact Patrick Eriksson.");
        }

      // Where will a geometric path exit the grid cell?
      do_gridrange_1d( r_v, lat_v, za_v, lstep, endface, r, lat, za, ppc_step, 
                                                       -1, r1, r3, r_surface );
      assert( r_v.nelem() == 2 );

      Numeric za_flagside = za;

      if( lstep <= lraytrace )
        {
          r           = r_v[1];
          lat         = lat_v[1];
          lcum       += lstep;
          za_flagside = za_v[1];
          ready       = true;
        }
      else
        {
          Numeric l;
          if( za <= 90 ) 
            { l = geompath_l_at_r(ppc_step,r) + lraytrace; }
          else
            { 
              l = geompath_l_at_r(ppc_step,r) - lraytrace; 
              if( l < 0 )
                { za_flagside = 80; }     // Tangent point passed!
            }

          r = geompath_r_at_l( ppc_step, l );  // Works als0 for l<0

          lat = geompath_lat_at_za( za, lat, 
                                geompath_za_at_r( ppc_step, za_flagside, r ) );
          lcum += lraytrace;
        }

      // Refractive index at *r*
      get_refr_index_1d( ws, refr_index_air, refr_index_air_group, 
                         refr_index_air_agenda, p_grid, refellipsoid, 
                         z_field, t_field, vmr_field, f_grid, r );

      // Calculate LOS zenith angle at found point.

      const Numeric ppc_local = ppc / refr_index_air; 

      if( r >= ppc_local )
        { za = geompath_za_at_r( ppc_local, za_flagside, r ); }
      else  // If moved below tangent point:
        { 
          r = ppc_local;
          za = 90;
        }

      // Store found point?
      if( ready  ||  ( lmax > 0  &&  lcum + lraytrace > lmax ) )
        {
          r_array.push_back( r );
          lat_array.push_back( lat );
          za_array.push_back( za );
          n_array.push_back( refr_index_air );
          ng_array.push_back( refr_index_air_group );
          l_array.push_back( lcum );
          lcum = 0;
        }
    }  
}



