/* Copyright (C) 2000, 2001 Stefan Buehler  <sbuehler@uni-bremen.de>
                            Axel von Engeln <engeln@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////

/*!
  \file   ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   Sun May  5  2002
  
  \brief  Functions releated to calculation of propagation paths.
  
  Functions to determine propagation paths for different atmospheric
  dimensionalities.

  The term propagation path is here shortened to ppath.
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "ppath.h"
#include "logic.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



////////////////////////////////////////////////////////////////////////////
//   Functions related to geometrical propagation paths.
////////////////////////////////////////////////////////////////////////////

//// geomppath_constant ///////////////////////////////////////////////////////
/**
   Calculates the propagation path constant for cases where refraction
   is neglected.

   The constant is r*sin(psi). 

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \return                The propagation path constant.
   \param    r_sensor     The radius of the sensor position.
   \param    psi_sensor   The zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-04-15
*/
Numeric geomppath_constant (
        const Numeric&   r_sensor,
        const Numeric&   psi_sensor )
{
  assert( r_sensor > 0 );
  assert( fabs(psi_sensor) <= 180 );

  return r_sensor * sin( DEG2RAD * psi_sensor );
}



//// geomppath_r2l ////////////////////////////////////////////////////////////
/**
   Calculates the length along the propagation path from the tangent point
   to the given radius.

   Note that the length has no sign, it is not clear from the length at what
   side of the tangent point the radius is placed.

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \return         The length along the ppath to the tangent point.
   \param    ppc   The ppath constant.
   \param    r     The radius of the point of concern.

   \author Patrick Eriksson
   \date   2002-04-15
*/
Numeric geomppath_r2l (
	const Numeric&   ppc, 
        const Numeric&   r )  
{
  assert( ppc >= 0 );
  assert( r >= ppc );

  // Double is hard-coded here to avoid numerical problems
  double a=r*r, b=ppc*ppc;

  return sqrt( a - b );
}



//// geomppath_l2r ////////////////////////////////////////////////////////////
/**
   Calculates the radius for a given length from the tangent point.

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \return                  Tha latitude of the tangent point.
   \param    ppc   The ppath constant.
   \param    psi_sensor     The zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-04-15
*/
Numeric geomppath_l2r (
	const Numeric&   ppc, 
        const Numeric&   l )
{
  assert( ppc >= 0 );

  // Double is hard-coded here to avoid numerical problems
  double a=l*l, b=ppc*ppc;

  return sqrt( a + b );
}



//// geomppath_r2psi //////////////////////////////////////////////////////////
/**
   Calculates the zenith angle for a given radius.

   The propagation direction is here not considered and the returned zenith 
   angle is always between 0 and 90 degrees.

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \return                The zenith angle of the path.
   \param    ppc          The ppath constant.
   \param    psi_sensor   The radius of the sensor position.

   \author Patrick Eriksson
   \date   2002-04-15
*/
Numeric geomppath_r2psi (
	const Numeric&   ppc, 
        const Numeric&   r )
{
  assert( ppc >= 0 );
  assert( r >= ppc );

  return RAD2DEG * asin( ppc / r );
}




////////////////////////////////////////////////////////////////////////////
//   Functions operating on the Ppath structure
////////////////////////////////////////////////////////////////////////////

//// empty_Ppath //////////////////////////////////////////////////////////////
/**
   Sets a Ppath structure to match the case when the propagation path is 
   totally outside the atmosphere. 

   The structure is then basicfally empty. The field *background* is set to
   ""Cosmic background radiation".

   \param   ppath Output:   A Ppath structure.
   \param    dim     The atmospheric dimensionality.

   \author Patrick Eriksson
   \date   2002-04-15
*/
void empty_Ppath( 
	      Ppath&      ppath,
	const Index&      dim )
{
  assert( dim >= 1 );
  assert( dim <= 3 );

  ppath.dim        = dim;
  ppath.np         = 0;
  ppath.pos.resize(0,dim);
  ppath.z.resize(0);
  ppath.l_step.resize(0);
  ppath.gridindex.resize(0,dim);
  ppath.los.resize(0,0);
  ppath.background = "Cosmic background radiation";
  ppath.ground     = 0;
  ppath.i_ground   = 0;
  ppath.tan_pos.resize(0);
  ppath.symmetry   = 0;
  ppath.i_symmetry = 0;
}



////////////////////////////////////////////////////////////////////////////
//   Main functions for 1D paths
////////////////////////////////////////////////////////////////////////////

//// ppath_1d_geom_upward /////////////////////////////////////////////////////
/**
   Calculates the path for a 1D atmosphere, no refraction, an upward
   observation and no cloud box.

   The lowest point of the path (r_start) and the sensor radius (r_sensor),
   if inside the atmosphere, are included as points of the path.

   If the path starts below an active cloud box, the z_field shall be cut at
   the lower limit of the box.

   The returned latitudes give the angular distance to the lowest point of the
   path.

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \param   r Output:            Radius for the points along the path.
   \param   l_step Output:       The distance along the path between the points.
   \param   psi Output:          Zenith angle for the points along the path.
   \param   alpha Output:        The angular distance for the points along the path.
                          Calculated with respect to the lowest point.
   \param   gridindex Output:    Grid index for the points along the path.
   \param    ppc          The ppath constant.
   \param    r_start      Radius for the starting point of the path.
   \param    psi_start    Zenith angle at the starting point of the path.
   \param    r_greoid     Geoid radius.
   \param    z_field      The geometricalk altitude of the pressure surfaces.
   \param    ppath_lmax   Maximum allowed length between path points.
   \param    r_sensor     Radius of the sensor position.

   \author Patrick Eriksson
   \date   2002-04-17
*/
void ppath_1d_geom_upward ( 
	      Vector&   r, 
	      Vector&   l_step, 
	      Vector&   psi, 
	      Vector&   alpha, 
	      Vector&   gridindex,
	const Numeric&  ppc,
	const Numeric&  r_start, 
	const Numeric&  psi_start, 
	const Numeric&  r_geoid, 
	const Vector&   z_field, 
	const Numeric&  ppath_lmax, 
	const Numeric&  r_sensor )
{
  // Number of z_field levels
  const Index nz = z_field.nelem();

  assert( ppc >= 0 );
  assert( r_start >= r_geoid+z_field[0] );
  assert( r_start < r_geoid+z_field[nz-1] );
  assert( psi_start >= 0 );
  assert( psi_start <= 180 );
  assert( ppath_lmax > 0 );
  assert( r_sensor >= r_geoid+z_field[0] );

  // Create a vector holding r_start and the z_field levels passed by the path.
  // If the sensor is inside the atmosphere, the sensor position is also 
  // included. Duplicates of radii are avoided. Store on the same time the 
  // integer part of gridindex.
  //
  // Find first z_field level above r_start (it cannot be 0).
  Index i_first;
  for( i_first=1; r_geoid+z_field[i_first]<=r_start; i_first++ ) {}
  //
  // Shall the sensor be put in? Setting i_sensor=nz, means no.
  Index i_sensor = nz;
  if( r_sensor < r_geoid + z_field[nz-1] )
    {
      for( i_sensor=i_first-1; r_geoid+z_field[i_sensor]<r_sensor; i_sensor++ )
        {}
      if( r_geoid+z_field[i_sensor] == r_sensor )
        i_sensor = nz;
    }
  //
  Vector       r1( nz - i_first + 1 + (i_sensor<nz) );
  ArrayOfIndex i1( r1.nelem() );
  //
  r1[0] = r_start;
  i1[0] = i_first - 1;
  for( Index i=i_first; i<i_sensor; i++ )
    {
      r1[i-i_first+1] = r_geoid + z_field[i];
      i1[i-i_first+1] = i;
    }
  if( i_sensor < nz )
    {
      r1[i_sensor-i_first+1] = r_sensor;
      i1[i_sensor-i_first+1] = i_sensor - 1;
      for( Index i=i_sensor; i<nz; i++ )
	{
          r1[i-i_first+2] = r_geoid + z_field[i];
	  i1[i-i_first+2] = i;
	}
    }


  // Determine how many points that are needed between the points in r1 to 
  // fulfill the length criterium set by ppath_lmax. The distance along the
  // path from the (imaginary) tangent point is stored in l1.
  //
  ArrayOfIndex n1( r1.nelem() - 1 );
  Vector       l1( r1.nelem() );
  Index        ntot=0;   // The total number of points
  //
  l1[0] = geomppath_r2l( ppc, r1[0] );
  //
  for( Index i=0; i<n1.nelem(); i++ )
    {
      l1[i+1] = geomppath_r2l( ppc, r1[i+1] );
      n1[i]   = Index( ceil( (l1[i+1]-l1[i]) / ppath_lmax ) );
      ntot   += n1[i];
    }
  ntot += 1;   // To include the point at the uppermost z_field level

  // We are now ready to create the output vectors.
  //
  r.resize( ntot );
  l_step.resize( ntot - 1 );
  psi.resize( ntot );
  alpha.resize( ntot );
  gridindex.resize( ntot );
  //
  ntot = 0;   // Now used as a counter for nummer of points done
  //
  for( Index i=0; i<n1.nelem(); i++ )
    {
      l_step[ntot] = ( l1[i+1] - l1[i] ) / n1[i];
      for( Index j=0; j<n1[i]; j++ )
	{
	  l_step[ntot]    = l_step[ntot-j];
	  r[ntot]         = geomppath_l2r( ppc, l1[i]+j*l_step[ntot] );
	  psi[ntot]       = geomppath_r2psi( ppc, r[ntot] );
	  alpha[ntot]     = psi[0] - psi[ntot];
	  gridindex[ntot] = i1[i] + ( r[ntot] - (r_geoid+z_field[i1[i]]) ) /
	                                 ( z_field[i1[i]+1] - z_field[i1[i]] );
	  ntot++;
	}
    }
  // The uppermost z_field level must be handled seperately
  r[ntot]         = last( r1 );
  psi[ntot]       = geomppath_r2psi( ppc, r[ntot] );
  alpha[ntot]     = psi[0] - psi[ntot];
  gridindex[ntot] = Numeric( last( i1 ) );
  ntot++;

  // All values fild?
  assert( r.nelem() == ntot );
}



//// ppath_1d /////////////////////////////////////////////////////////////////
/**
   Performs the logistics for determining the path for 1D cases, with or
   without refraction.

   The actual calculations are mainly performed by ppath_1d_geom_upward and
   ppath_1d_refr_upward.

   \param   ppath Output:          Propagation path structure with all fields set.
   \param    p_grid         Pressure grid.
   \param    z_field        The geometricalk altitude of the pressure surfaces.
   \param    r_greoid       Geoid radius.
   \param    refr_on        Refraction flag.
   \param    blackbody_ground Flag to indicate a blackbody ground.
   \param    cloudbox_on    Cloud box flag.
   \param    cloudbox_limits Vertical limits of the cloud box.
   \param    ppath_lmax     Maximum allowed length between path points.
   \param    r_sensor       Radius of the sensor position.
   \param    psi_sensor     Zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-04-17
*/
void ppath_1d (
	      Ppath&          ppath,
	ConstVectorView       p_grid,
	ConstVectorView       z_field,
	const Numeric&        r_geoid,
	const Numeric&        z_ground,
        const Index&          refr_on,
        const Index&          blackbody_ground,
	const Index&          cloudbox_on,
	const ArrayOfIndex&   cloudbox_limits,
	const Numeric&        ppath_lmax,
	const Numeric&        r_sensor,
        const Numeric&        psi_sensor )
{
  // Asserts.
  assert( p_grid.nelem() == z_field.nelem() );
  assert( z_field.nelem() > 1 );
  assert( z_ground >= z_field[0] );
  assert( z_ground < last(z_field) );
  assert( is_bool( refr_on ) );
  assert( is_bool( blackbody_ground ) );
  assert( is_bool( cloudbox_on ) );
  if( cloudbox_on )
    {
      assert( cloudbox_limits.nelem() == 2 );
      assert( cloudbox_limits[0] < cloudbox_limits[1] );
      assert( cloudbox_limits[0] >= 0 );
      if( !blackbody_ground )
        assert( cloudbox_limits[0] == 0 );
      assert( cloudbox_limits[1] < p_grid.nelem() );
    }
  assert( ppath_lmax > 0 );
  assert( psi_sensor >= 0 );
  assert( psi_sensor <= 180 );

  // Check that the sensor not is below the ground.
  if( r_sensor < r_geoid + z_ground )
    {
      ostringstream os;
      os << "The sensor position cannot be below the ground altitude.\n"
         << "The ground altitude is here " << z_ground/1e3 << "km and the " 
         << "sensor altitude is " << (r_sensor-r_geoid)/1e3 << " km."; 
      throw runtime_error( os.str() );
    }


  // Start by setting PPATH to be empty 
  empty_Ppath( ppath, 1 );

  // Calculate the path constant assuming no refraction. For the moment, ppc
  // will only be used outside the atmosphere where the refractive index is 1.
  Numeric ppc = geomppath_constant ( r_sensor, psi_sensor );

  // Get length of z_field and the radius for the top of the atmosphere
  const Index     n_zfield = z_field.nelem();
  const Numeric   r_top    = r_geoid + z_field[n_zfield-1];

  // Radius of upper and lower limit of the cloud box
  Numeric r_cb_low, r_cb_upp;
  if( cloudbox_on )
    {
      r_cb_low = r_geoid + z_field[cloudbox_limits[0]];
      r_cb_upp = r_geoid + z_field[cloudbox_limits[1]];
    }


  // Handle the cases where there is no path to follow. For these cases, the
  // the spectrum equals the background. If the sensor is inside the atmosphere
  // the ppath fields *pos* and *los* must be set to make it possible to 
  // extract the background. 
  //
  Index do_path = 1, put_in_sensor = 0;
  //
  // Path totally outside the atmosphere?. If yes, then ppath is already OK.
  if( r_sensor>=r_top && ( psi_sensor<=90 || ( psi_sensor>90 && ppc>=r_top ) ))
    {
      do_path = 0;
    }
  // Standing on the ground looking and looking down into a blackbody ground?
  if( r_sensor==(r_geoid+z_ground) && psi_sensor>90 && blackbody_ground )
    {
      do_path          = 0;
      put_in_sensor    = 1;
      ppath.background = "Blackbody ground";
    }
  // The sensor is inside an active cloud box?
  if( cloudbox_on && r_sensor>r_cb_low && r_sensor<r_cb_upp )
    {
      do_path          = 0;
      put_in_sensor    = 1;
      ppath.background = "Inside cloud box";
    }
  // The sensor is on the boundary of an active cloud box and looks into 
  // the box?
  if( cloudbox_on && ( ( r_sensor==r_cb_low && psi_sensor<=90 )
                               || ( r_sensor==r_cb_upp && psi_sensor>90 ) ) )
    {
      do_path          = 0;
      put_in_sensor    = 1;
      ppath.background = "Surface of cloud box";
    }


  // Put in position and zenith angle of sensor if there is no path to follow.
  // If there is no path to follow, there is nothing more to do.
  if( !do_path )
    {
      if( put_in_sensor )
	{
          ppath.np = 1;
          ppath.z.resize(1);
          ppath.z[0] = r_sensor - r_geoid;
	  ppath.pos.resize(1,2);
	  ppath.pos(0,0) = ppath.z[0];
	  ppath.pos(0,1) = 0;   // The latitudes start here at the sensor.
	  ppath.los.resize(1,1);
	  ppath.los(0,0) = psi_sensor;
	}
    }
  

  // There is a path to follow.
  else
    {
      // A path from the lowest point and upward is first determined. This path
      // is then mirrored if necessary. The mirror part can be made common for
      // the geometrical and refraction options.

      // Recalculate the path constant if refrection is considered.
      if( refr_on )
	{} //ppc = ???; FIX THIS.

      // Find the lowest point of the path (r_low) and the zenith angle at this
      // point (psi_low). Determine on the same time if there is a ground
      // intersection, if symmetry applies and the radiative background.
      // The case when the sensor is below a cloud box is handled below.
      //
      Numeric r_low, psi_low;
      Numeric r_tan;
      Index stops_at_rlow = 0;        // True for downward obs. but no symmetry
      //
      if( psi_sensor <= 90 )
	{
	  r_low      = r_sensor;
	  psi_low    = psi_sensor;
	  ppath.background = "Cosmic background radiation";
	}
      else
	{
	  if( refr_on )
	    {} // r_tan = ??? FIX THIS.
	  else
	    r_tan = ppc;
	  //
	  // Intersection with cloud box from above?
	  if( cloudbox_on && r_sensor>r_cb_upp && r_tan<r_cb_upp )
	    {
	      r_low   = r_cb_upp;
	      if( refr_on )
		{} // FIX THIS
	      else
		psi_low = geomppath_r2psi ( ppc, r_low );
	      ppath.background = "Surface of cloud box";
	      stops_at_rlow = 1;
	    }
	  //
	  // Intersection with the ground?
	  else if( r_tan < (r_geoid+z_ground) )
	    {
	      r_low = r_geoid + z_ground;
	      if( refr_on )
		{} // FIX THIS
	      else
		psi_low = geomppath_r2psi ( ppc, r_low );
	      if( blackbody_ground )
		{
		  ppath.background = "Blackbody ground";
		  stops_at_rlow = 1;
		}
	      else
		{
		  ppath.ground   = 1;
		  ppath.symmetry = 1;
		  ppath.background     = "Cosmic background radiation";
		}
	    }
	  //
	  // If not any of the cases above, path through clear atmosphere
	  else
	    {
	      r_low          = r_tan;
	      psi_low        = 90;
	      ppath.symmetry = 1;
  	      ppath.background     = "Cosmic background radiation";
	    }
	}

      // Check if the practical upper limit of the atmosphere is the cloud box.
      // The upper limit is specified by the altitude index (i_max).
      //
      Index i_max = n_zfield - 1;
      //
      if( cloudbox_on && r_sensor<=r_cb_low )
	{
	  i_max      = cloudbox_limits[0];
	  if( ppath.background == "Cosmic background radiation" )
	    ppath.background = "Surface of cloud box";
	}


      // Calculate the path from the lowest point to the upper limit.
      // If the sensor is inside the atmosphere and looks down, the sensor
      // position is included as a point of the path.
      // This path is described by the radii (r_path1), the zenith angle 
      // (psi_path1) and lengths along the path(l_path1). Grid indexes are
      // calculated on the same time.
      //
      Vector r_path1, psi_path1, alpha_path1, l_path1, gi_path1;
      //
      if( refr_on )
	{} // FIX THIS
      else
 	ppath_1d_geom_upward( r_path1, l_path1, psi_path1, alpha_path1, 
                    gi_path1, ppc, r_low, psi_low,  r_geoid, 
                             z_field[Range(0,i_max+1)], ppath_lmax, r_sensor );

      
      // Determine at what index of r_path1 the sensor is placed. If the
      // sensor is above the atmosphere, the upper end index is used.
      // The index is then used to mirror the path if a symmetry point exist.
      // The latitude difference between the sensor and r_low is alpha0.
      //
      Index   i_sensor = 0;      // This is correct if an upward observation.
      Numeric alpha0   = 0;      // Again correct for upward observation.
      //
      if( psi_sensor > 90 )
	{
	  if( r_sensor > r_top )   
	    {
	      i_sensor = r_path1.nelem() - 1;
	      alpha0  = last(alpha_path1) + psi_sensor - 180 + last(psi_path1);
	    }
	  else   // The sensor is somewhere in there
	    {
	      for( i_sensor=0; r_sensor!=r_path1[i_sensor]; i_sensor++ ) {}
	      alpha0 = alpha_path1[i_sensor];
	    }
	  ppath.i_symmetry = i_sensor;	    
	}
      //
      // Resize fields z, los and l_step. The field z is used as temporary
      // storage for the radii.
      Index n_path1  = r_path1.nelem();
      if( stops_at_rlow )
	ppath.np = i_sensor + 1;
      else
	ppath.np = n_path1 + i_sensor;
      ppath.pos.resize(ppath.np,2);
      ppath.z.resize( ppath.np );
      ppath.l_step.resize( ppath.np - 1 );
      ppath.los.resize( ppath.np, 1 );
      ppath.gridindex.resize( ppath.np, 1 );
      //
      // If i_sensor>0, put in part between sensor and the lowest point.
      // The zenith angle is here 180 - psi_path1.
      // Latitudes give here the angular distance to the lowest point.
      if( i_sensor > 0 )
	{
	  ppath.z[ Range(0,i_sensor) ] = 
                                        r_path1[ Range(i_sensor,i_sensor,-1) ];
	  ppath.l_step[ Range(0,i_sensor) ] = 
                                      l_path1[ Range(i_sensor-1,i_sensor,-1) ];
	  for( Index i=0; i<i_sensor; i++ )
	    {
	      ppath.pos( i, 1 ) = -alpha_path1[ i_sensor-i ];
	      ppath.los( i, 0 ) = 180 - psi_path1[ i_sensor-i ];
	    }
	  ppath.gridindex( Range(0,i_sensor), 0 ) = 
                                       gi_path1[ Range(i_sensor,i_sensor,-1) ];
	}
      //
      if( stops_at_rlow )
	{
	  ppath.z[ i_sensor ]            = r_path1[0];
	  ppath.pos( i_sensor, 1 )       = alpha_path1[0];
	  ppath.los( i_sensor, 0 )       = psi_path1[0];
	  ppath.gridindex( i_sensor, 0 ) = gi_path1[0];
	}
      else
	{
	  ppath.z[ Range(i_sensor,n_path1) ]            = r_path1;
	  ppath.l_step[ Range(i_sensor,n_path1-1) ]     = l_path1;
	  ppath.pos( Range(i_sensor,n_path1), 1 )       = alpha_path1;
	  ppath.los( Range(i_sensor,n_path1), 0 )       = psi_path1;
	  ppath.gridindex( Range(i_sensor,n_path1), 0 ) = gi_path1;
	}

      // Fill remaining fields of ppath and adjust latitudes with alpha0.
      for( Index i=0; i<ppath.np; i++ )
	{
	  ppath.z[i] -= r_geoid;
	  ppath.pos( i, 0 )  = ppath.z[i];
	  ppath.pos( i, 1 ) += alpha0;
	}
      if( ppath.ground )
	ppath.i_ground = ppath.i_symmetry;
      if( psi_sensor > 90 )
	{
	  ppath.tan_pos.resize(2);
	  ppath.tan_pos[0] = r_tan - r_geoid;
	  ppath.tan_pos[1] = ppath.pos(i_sensor,1) + 90 -ppath.los(i_sensor,0);
	}

    } // do_path=1.
}




////////////////////////////////////////////////////////////////////////////
//   ppath_calculate
////////////////////////////////////////////////////////////////////////////

//// ppath_calculate //////////////////////////////////////////////////////////
/**
   The main function for determining the propagation path for a single
   combination observation position and direction.

   The main task of this function is to check the input.

   \param   ppath Output:          Propagation path structure with all fields set.
   \param    dim            Atmospheric dimensionality.
   \param    p_grid         Pressure grid.
   \param    alpha_grid     Latitude grid.
   \param    beta_grid      Longitude grid.
   \param    z_field        The geometrical altitude of the pressure surfaces.
   \param    r_greoid       Geoid radius.
   \param    z_ground       Ground altitude.
   \param    e_ground       Emissivity of the ground.
   \param    refr_on        Refraction flag.
   \param    cloudbox_on    Cloud box flag.
   \param    cloudbox_limits Vertical limits of the cloud box.
   \param    ppath_lmax     Maximum allowed length between path points.
   \param    sensor_pos     Sensor position.
   \param    sensor_los     Sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-04-17
*/
void ppath_calculate(
	      Ppath&          ppath,
	const Index&          dim,
	ConstVectorView       p_grid,
	ConstVectorView       alpha_grid,
	ConstVectorView       beta_grid,
	const Tensor3&        z_field,
 	ConstMatrixView       r_geoid,
 	ConstMatrixView       z_ground,
	const Tensor3&        e_ground,
        const Index&          refr_on,
	const Index&          cloudbox_on,
	const ArrayOfIndex&   cloudbox_limits,
	const Numeric&        ppath_lmax,
	ConstVectorView       sensor_pos,
	ConstVectorView       sensor_los )
{
  // Check input:
  // Atmospheric grids are checked inside chk_atm_grids.
  // Consistency between the grids and sensor position checked in sub-funs.
  chk_if_in_range( "dim", dim, 1, 3 );
  chk_atm_grids( dim, p_grid, alpha_grid, beta_grid );
  chk_atm_field( "z_field", z_field, dim, p_grid, alpha_grid, beta_grid );
  for( Index row=0; row<z_field.nrows(); row++ )
    {
      for( Index col=0; col<z_field.ncols(); col++ )
	{
	  // I could not get the compliler to accept a solution without dummy!!
	  Vector dummy(z_field.npages());
	  dummy = z_field(Range(joker),row,col);
	  ostringstream os;
	  os << "z_field for latitude nr " << row << " and longitude nr " 
             << col;
	  chk_if_increasing( os.str(), dummy ); 
	}
    }
  chk_atm_surface( "r_geoid", r_geoid, dim, alpha_grid, beta_grid );
  chk_atm_surface( "z_ground", z_ground, dim, alpha_grid, beta_grid );
  for( Index row=0; row<z_ground.nrows(); row++ )
    {
      for( Index col=0; col<z_ground.ncols(); col++ )
	{
	  if( z_ground(row,col)<z_field(0,row,col) || 
                       z_ground(row,col)>=z_field(z_field.npages()-1,row,col) )
	    {
	      ostringstream os;
	      os << "The ground altitude (*z_ground*) cannot be outside of the"
		 << " altitudes in *z_field*.";
		if( dim > 1 )
		  os << "\nThis was found to be the case for:\n" 
		     << "latitude " << alpha_grid[row];
		if( dim > 2 )
		  os << "\nlongitude " << beta_grid[col];
	      throw runtime_error( os.str() );
	    }
	}
    }
  chk_atm_surface( "e_ground (for one frequency)", 
           e_ground(0,Range(joker),Range(joker)), dim, alpha_grid, beta_grid );
  chk_if_bool(  "refr_on", refr_on );
  // The cloud box is checked below.
  chk_if_over_0( "ppath_lmax", ppath_lmax );
  chk_vector_length( "sensor_pos", sensor_pos, dim );
  if( dim<3 )
    chk_vector_length( "sensor_los", sensor_los, 1 );
  else
    {
      chk_vector_length( "sensor_los", sensor_los, 2 );
      chk_if_in_range( "sensor_los azimuth angle", sensor_los[1], -180, 180 );
    }
  chk_if_in_range( "sensor_los zenith angle", sensor_los[0], -180, 180 );


  // Is the ground is perfect blackbody everywhere?
  Index blackbody_ground = 1;
  for( Index i=0; blackbody_ground && i<e_ground.ncols(); i++ )
    {
      if( min( e_ground(Range(joker),Range(joker),i) ) < 0.99999 )
	blackbody_ground = 0;
    }


  // Check of the cloud box
  chk_if_bool(  "cloudbox_on", cloudbox_on );
  if( cloudbox_on )
    {
      if( cloudbox_limits.nelem() != dim*2 )
	{
	  ostringstream os;
	  os << "The array *cloudbox_limits* has incorrect length.\n"
	     << "For dim = " << dim << " the length shall be " << dim*2
	     << " but it is " << cloudbox_limits.nelem() << ".";
	  throw runtime_error( os.str() );
	}
      if( !blackbody_ground && cloudbox_limits[0]!=0 )
	{
	  ostringstream os;
	  os << "The lower pressure limit for cloud box must be 0 when the"
             << "ground\nis not treated to be a blackbody, but the limit is"
	     << "set to be " << cloudbox_limits[0] << ".";
	  throw runtime_error( os.str() );
	}
       if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
	{
	  ostringstream os;
	  os << "Incorrect value(s) for cloud box pressure limit(s) found."
	     << "\nValues are either out of range or upper limit is not "
	     << "greater than lower limit.\nWith present length of "
             << "*p_grid*, OK values are 0 - " << p_grid.nelem()-1
             << ".\nThe pressure index limits are set to " 
	     << cloudbox_limits[0] << " - " << cloudbox_limits[1] << ".";
	  throw runtime_error( os.str() );
	}
      if( dim >= 2 )
	{
	  if( cloudbox_limits[3]<=cloudbox_limits[2] || cloudbox_limits[2]<1 ||
                                cloudbox_limits[3]>=alpha_grid.nelem()-1 )
	    {
	      ostringstream os;
	      os << "Incorrect value(s) for cloud box latitude limit(s) found."
		 << "\nValues are either out of range or upper limit is not "
		 << "greater than lower limit.\nWith present length of "
                 << "*alpha_grid*, OK values are 1 - " << alpha_grid.nelem()-2
                 << ".\nThe latitude index limits are set to " 
		 << cloudbox_limits[2] << " - " << cloudbox_limits[3] << ".";
	      throw runtime_error( os.str() );
	    }
	}
      if( dim >= 3 )
	{
	  if( cloudbox_limits[5]<=cloudbox_limits[4] || cloudbox_limits[4]<1 ||
                                cloudbox_limits[5]>=beta_grid.nelem()-1 )
	    {
	      ostringstream os;
	      os << "Incorrect value(s) for cloud box longitude limit(s) found"
		 << ".\nValues are either out of range or upper limit is not "
		 << "greater than lower limit.\nWith present length of "
                 << "*beta_grid*, OK values are 1 - " << beta_grid.nelem()-2
                 << ".\nThe longitude limits are set to " 
		 << cloudbox_limits[4] << " - " << cloudbox_limits[5] << ".";
	      throw runtime_error( os.str() );
	    }
	}
    }


  // 1D cases are treated seperately to make use of the symmetry of the
  // calculations, while 2D and 3D are treated with the same function.
  if( dim == 1 )
    ppath_1d( ppath, p_grid, z_field(Range(joker),0,0), r_geoid(0,0),
                 z_ground(0,0), refr_on, blackbody_ground, cloudbox_on, 
             cloudbox_limits, ppath_lmax, sensor_pos[0], fabs(sensor_los[0]) );
  else
    throw runtime_error(
                          "2D and 3D PPATH calculations not yet implemented" );
}


void test_new_ppath()
{
  Ppath ppath;
 
  // Dimension
  Index dim = 1;

  // Grids
  const Index n_p = 11;
  Vector p_grid( 10, n_p, -1 );
  Vector alpha_grid(0), beta_grid(0);

  // z_field
  Tensor3 z_field(n_p,1,1);
  for( Index i=0; i<n_p; i++ )
    z_field(i,0,0) = i*1e3;

  // Geoid and ground
  Tensor3 e_ground(5,1,1);
  Matrix r_geoid(1,1), z_ground(1,1);
  Numeric r0 = 6.378e6;
  r_geoid(0,0)  = r0;
  z_ground(0,0) = 0;
  e_ground      = 0.1;

  // Booleans
  Index refr_on = 0;

  // Cloud box
  Index cloudbox_on;
  ArrayOfIndex   cloudbox_limits(2);
  cloudbox_on = 1;
  cloudbox_limits[0] = 0;
  cloudbox_limits[1] = 4;

  // Path step length
  Numeric ppath_lmax;
  ppath_lmax = 20e3;

  // Sensor  
  Vector sensor_pos(1), sensor_los(1);
  sensor_pos[0] = r0 + 0e3;
  sensor_los[0] = 90.1;

  ppath_calculate( ppath, dim, p_grid, alpha_grid, beta_grid, z_field, r_geoid,
 	         z_ground, e_ground, refr_on, cloudbox_on, cloudbox_limits,
                                          ppath_lmax, sensor_pos, sensor_los );

  //  cout << "z = \n" << ppath.z << "\n";
  //  cout << "gridindex = \n" << ppath.gridindex << "\n";
  cout << "pos = \n" << ppath.pos << "\n";
  cout << "los = \n" << ppath.los << "\n";
  cout << "l_step = \n" << ppath.l_step << "\n";
  cout << "np      = " << ppath.np << "\n";
  cout << "bground = " << ppath.background << "\n";
  cout << "ground  = " << ppath.ground << "\n";
  if( ppath.ground )
    cout << "i_ground= " << ppath.i_ground << "\n";
  if( ppath.tan_pos.nelem() > 0 )
  {
    cout << "z_tan   = " << ppath.tan_pos[0]/1e3 << " km\n";
    cout << "a_tan   = " << ppath.tan_pos[1] << " degs\n";
  }  
  cout << "symmet = " << ppath.symmetry << "\n";
  if( ppath.symmetry )
    cout << "i_symmet= " << ppath.i_symmetry << "\n";
}
