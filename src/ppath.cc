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

/**
  \file   ppath.cc

  Functions to determine propagation paths for different atmospheric
  dimensionalities.

  The term propagation path is here shortened to ppath.

  The main function for these calculations is ppath_calculate, that is
  found at the end of the file.

  \author Patrick Eriksson
  \date 2002-04-15 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "ppath.h"
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



//// geomppath_alpha_tan //////////////////////////////////////////////////////
/**
   Calculates the angular position in the observation plane of the tangent
   point. The angular distance in the observation plane is denoted as 
   the latitude for 1D and 2D cases.

   For upward observations (abs(psi)>=90), the tangent point is found along
   an imaginary path behind the sensor. 

   Negative zenith angles means that the sensor is directed towards lower
   latitudes.

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \return                  The latitude of the tangent point.
   \param    alpha_sensor   The latitude of the sensor position.
   \param    psi_sensor     The zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-04-15
*/
Numeric geomppath_alpha_tan (
	const Numeric&   alpha_sensor,
        const Numeric&   psi_sensor )
{
  assert( fabs(psi_sensor) <= 180 );

  Numeric dalpha = psi_sensor - 90; 

  if( fabs( psi_sensor ) > 90 )
    dalpha = -dalpha;

  return alpha_sensor + dalpha;
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
	const double&   ppc,     // double is hard-coded here to avoid
        const double&   r )      // possible numerical problems with float.
{
  assert( r >= ppc );

  return sqrt( r*r - ppc*ppc );
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
	const double&   ppc,     // double is hard-coded here to avoid
        const double&   l )      // possible numerical problems with float.
{
  return sqrt( ppc*ppc + l*l );
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

   See further the ARTS user guide (AUG). Use the index to find where
   this function is discussed. The function is listed as a subentry to
   "internal ARTS functions".

   \retval   ppath   A Ppath structure.
   \param    dim     The atmospheric dimensionality.

   \author Patrick Eriksson
   \date   2002-04-15
*/
void empty_Ppath( 
	      Ppath&      ppath,
	const Index&      dim )
{
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
  // Asserts. Variables are also asserted in ???.
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
      if( ~blackbody_ground )
        assert( cloudbox_limits[0] == 0 );
      assert( cloudbox_limits[1] < p_grid.nelem() );
    }
  assert( ppath_lmax > 0 );
  assert( r_sensor >= r_geoid + z_ground );
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
  if( r_sensor>=r_top && ppc>=r_top )
    {
      do_path = 0;
    }
  // Standing on the ground looking and looking down into a blackbody ground?
  if( r_sensor>=(r_geoid+z_ground) && psi_sensor>90 && blackbody_ground )
    {
      do_path          = 0;
      put_in_sensor    = 0;
      ppath.background = "Blackbody ground";
    }
  // The sensor is inside an active cloud box?
  if( cloudbox_on && r_sensor>r_cb_low && r_sensor<r_cb_upp )
    {
      do_path          = 0;
      put_in_sensor    = 0;
      ppath.background = "Inside cloud box";
    }
  // The sensor is on the boundary of an active cloud box and looks into 
  // the box?
  if( cloudbox_on && ( ( r_sensor==r_cb_low && psi_sensor<=90 )
                               || ( r_sensor==r_cb_upp && psi_sensor>90 ) ) )
    {
      do_path          = 0;
      put_in_sensor    = 0;
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
      String  background;
      Numeric r_tan = -1;      // Radius of tangent point. -1 means upward obs.
      if( psi_sensor >= 90 )
	{
	  r_low      = r_sensor;
	  psi_low    = psi_sensor;
	  background = "Cosmic background radiation";
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
	      background = "Surface of cloud box";
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
		background = "Blackbody ground";
	      else
		{
		  ppath.ground   = 1;
		  ppath.symmetry = 1;
		  background     = "Cosmic background radiation";
		}
	    }
	  //
	  // If not any of the cases above, path through clear atmosphere
	  else
	    {
	      r_low          = r_tan;
	      psi_low        = 90;
	      ppath.symmetry = 1;
  	      background     = "Cosmic background radiation";
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
	  if( background == "Cosmic background radiation" )
	    background = "Surface of cloud box";
	}


      // Calculate the path from the lowest point to the upper limit.
      // If the sensor is inside the atmosphere and looks down, the sensor
      // position is included as a point of the path.
      // This path is described by the radii (r_path1), the zenith angle 
      // (psi_path1) and lengths along the path(l_path1).
      //
      Vector r_path1, psi_path1, l_path1;
      //
      if( refr_on )
	{} // FIX THIS
      else
	{}

      
      // Determine at what index of r_path1 the sensor is placed. If the
      // sensor is above the atmosphere, the upper end index is used.
      // The index is then used to mirror the path if a symmetry point exist.
      //
      Index i_sensor = 0;   // This is correct if an upward observation.
      //
      if( ppath.symmetry )
	{
	  if( r_sensor > r_top )   
	    i_sensor = r_path1.nelem() - 1;
	  else   // The sensor is somewhere in there
	    for( i_sensor=0; r_sensor!=r_path1[i_sensor]; i_sensor++ ) {}
	  ppath.i_symmetry = i_sensor;	    
	}
      //
      // Resize fields z, los and l_step. The field z is used as temporary
      // storage for the radii.
      Index n_path1 = r_path1.nelem();
      ppath.np = n_path1 + i_sensor;
      ppath.z.resize(ppath.np);
      ppath.l_step.resize(ppath.np);
      ppath.los.resize(ppath.np,1);
      //
      // If i_sensor>0, put in part between sensor and the lowest point.
      // The zenith angle is here 180 - psi_path1.
      if( i_sensor > 0 )
	{
	  ppath.z[ Range(0,i_sensor) ] = 
                                        r_path1[ Range(i_sensor,i_sensor,-1) ];
	  ppath.l_step[ Range(0,i_sensor) ] = 
                                        l_path1[ Range(i_sensor,i_sensor,-1) ];
	  for( Index i; i<i_sensor; i++ )
	    ppath.los( i, 0 ) = 180 - psi_path1[ i ];
	}
      //
      ppath.z[ Range(i_sensor,n_path1) ]      = r_path1;
      ppath.l_step[ Range(i_sensor,n_path1) ] = l_path1;
      ppath.los( Range(i_sensor,n_path1), 0 ) = psi_path1;

      // Fill remaining fields of ppath
      ppath.pos.resize(ppath.np,2);
      for( Index i; i<ppath.np; i++ )
	{
	  ppath.pos( i, 0 ) = ppath.z[i] - r_geoid;
	  ppath.pos( i, 1 ) = ppath.los( i, 0 ) - 90;
	}
      if( ppath.ground )
	ppath.i_ground = ppath.symmetry;
      if( r_tan > 0 )
	{
	  ppath.tan_pos.resize(2);
	  ppath.tan_pos[0] = r_tan - r_geoid;
	  ppath.tan_pos[1] = ppath.pos( i_sensor, 1 );
	}

    } // do_path=1.
}




////////////////////////////////////////////////////////////////////////////
//   ppath_calculate
////////////////////////////////////////////////////////////////////////////

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
                         z_ground(row,col)>z_field(z_field.npages(),row,col) )
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
	     << "For dim = " << dim << "the length shall be " << dim*2
	     << "but it is " << cloudbox_limits.nelem() << ".";
	  throw runtime_error( os.str() );
	}
      if( !blackbody_ground && cloudbox_limits[0]!=0 )
	{
	  ostringstream os;
	  os << "The lower pressure limit for cloud box must be 0 when the\n"
             << "ground is not treated to be a blackbody, but the limit is\n"
	     << "set to " << cloudbox_limits[0] << ".";
	  throw runtime_error( os.str() );
	}
       if( cloudbox_limits[1]<=cloudbox_limits[0] || cloudbox_limits[0]<0 ||
                                           cloudbox_limits[1]>=p_grid.nelem() )
	{
	  ostringstream os;
	  os << "Incorrect value(s) for cloud box pressure limit(s) found."
	     << "\nValues are either out of range or upper limit is not "
	     << "greater than lower limit.\nWith present length of "
             << "*p_grid*, OK values are 0 - " << alpha_grid.nelem()-1
             << ".\nThe latitude index limits are set to " 
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
  e_ground      = 0.8;

  // Booleans
  Index refr_on = 0;

  // Cloud box
  Index cloudbox_on;
  ArrayOfIndex   cloudbox_limits;
  cloudbox_on = 0;
  cloudbox_limits.resize(2);
  cloudbox_limits[0] = 0;
  cloudbox_limits[0] = 4;

  // Path step length
  Numeric ppath_lmax;
  ppath_lmax = 50e3;

  // Sensor  
  Vector sensor_pos(1), sensor_los(1);
  sensor_pos[0] = r0 + 20e3;
  sensor_los[0] = 98;

  ppath_calculate( ppath, dim, p_grid, alpha_grid, beta_grid, z_field, r_geoid,
 	         z_ground, e_ground, refr_on, cloudbox_on, cloudbox_limits,
                                          ppath_lmax, sensor_pos, sensor_los );

  cout << "z = \n" << ppath.z << "\n";
  cout << "gridindex = \n" << ppath.gridindex << "\n";
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
