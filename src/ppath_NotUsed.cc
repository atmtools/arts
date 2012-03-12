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
   \param   rsurface15   Radius for the surface at *lat1* and *lon5*.
   \param   rsurface35   Radius for the surface at *lat3* and *lon5*.
   \param   rsurface36   Radius for the surface at *lat3* and *lon6*.
   \param   rsurface16   Radius for the surface at *lat1* and *lon6*.

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

      if( rt > 0  &&  lt < l )  // lt<l to resolve the closest crossing
        { endface = 7;   r = rt;   lat = latt;   lon = lont;   l = lt; }
    }

  // If crossing found (r>0) we are ready!
  // (plevel_crossing_3d checks if crossing is inside grid box)
  
  // Upper pressure level
  if( r <= 0 )
    {
      plevel_crossing_3d( r, lat, lon, l, r_start, lat_start, lon_start, 
                          za_start, aa_start, x, y, z, dx, dy, dz, ppc, 
                          lat1, lat3, lon5, lon6, r15b, r35b, r36b, r16b, 
                          false );
      if( r > 0 )
        { endface = 4; }
    }
  
  // Latitude and longitude endfaces
  if( r <= 0 )
    {
      // Here we test both sides blindly and takes shortest l as solution:

      // Latitude
      Numeric rlat, latlat = lat1, lonlat, llat;
      Index   eflat = 1;
      //
      lat_crossing_3d( rlat, lonlat, llat, lat1, lat_start, za_start,
                                                         x, y, z, dx, dy, dz );
      {
        Numeric rlat3, lonlat3, llat3;
        lat_crossing_3d( rlat3, lonlat3, llat3, lat3, lat_start, za_start,
                                                         x, y, z, dx, dy, dz );
        if( rlat3 > 0  &&  llat3 < llat )
          { 
            eflat = 3; 
            rlat = rlat3; latlat = lat3; lonlat = lonlat3; llat = llat3; 
          }
      }
 
      // Longitude
      Numeric rlon, latlon, lonlon = lon5, llon;
      Index   eflon = 5;
      //
      lon_crossing_3d( rlon, latlon, llon, lon5, lon_start, za_start, aa_start, 
                                                         x, y, z, dx, dy, dz );
      //      
      {
        Numeric rlon6, latlon6, llon6;
        lon_crossing_3d( rlon6, latlon6, llon6, lon6, lon_start, za_start, 
                                               aa_start, x, y, z, dx, dy, dz );
        if( rlon6 > 0  &&  llon6 < llon )
          { 
            eflon = 6; 
            rlon = rlon6; latlon = latlon6; lonlon = lon6; llon = llon6; 
          }        
      }

      // Pick out solution with shortest l
      if( llat < llon )
        { r = rlat; lat = latlat; lon = lonlat; l = llat; endface = eflat; }
      else
        { r = rlon; lat = latlon; lon = lonlon; l = llon; endface = eflon; }
      assert( lat >= lat1 );
      assert( lat <= lat3 );
      assert( lon >= lon5 );
      assert( lon <= lon6 );
    }

  assert( endface );

  // Check if there is a tangent point inside the grid cell. 
  if( za_start > 90  )
    {
      Numeric ltan = geompath_l_at_r( ppc, r_start);
      if( l > ltan ) 
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
                                       x+dx*lj, y+dy*lj, z+dz*lj, dx, dy, dz );
   }

  //--- Set last point especially, which should improve the accuracy
  r_v[n]   = r; 
  lat_v[n] = lat;
  lon_v[n] = lon;

  //--- Set last zenith angle to be as accurate as possible
  if( za_start < ANGTOL  ||  za_start > 180-ANGTOL )
    { za_v[n] = za_start; }
  else if( endface == 8 )
    { za_v[n] = 90; }
  else
    { za_v[n] = geompath_za_at_r( ppc, za_start, r_v[n] ); }

  //--- Set last azimuth angle to be as accurate as possible for
  //    zenith and nadir observations
  if( abs( lat_start ) < POLELAT  &&  
          ( abs(aa_start) < ANGTOL  ||  abs( aa_start) > 180-ANGTOL ) )
    {  
      aa_v[n]  = aa_start; 
      lon_v[n] = lon_start;
    }

  // Shall lon values be shifted (value 0 and n+1 are already OK)?
  for( Index j=1; j<n; j++ )
    { resolve_lon( lon_v[j], lon5, lon6 ); }
}



