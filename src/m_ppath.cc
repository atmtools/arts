/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>
                            Patrick Eriksson <patrick@rss.chalmers.se>
                            
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
  \file   m_ppath.cc
  \brief  Workspace functions releated to calculation of propagation paths.

  \author Patrick Eriksson
  \date 2002-05-08 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "ppath.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////

void ppathCalc(
        // WS Output:
              Ppath&          ppath,
        // WS Input:
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Index&          blackbody_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         sensor_pos,
        const Vector&         sensor_los )
{
  //--- Check input -----------------------------------------------------------

  // Check that data sizes of grids, *z_field* and ground variables match the 
  // atmospheric dimensionality.
  //
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_ground", z_ground, atmosphere_dim, lat_grid, lon_grid );

  // Check that *z_field* has strictly increasing pages.
  //
  for( Index row=0; row<z_field.nrows(); row++ )
    {
      for( Index col=0; col<z_field.ncols(); col++ )
	{
	  // I could not get the compliler to accept a solution without dummy!!
	  Vector dummy(z_field.npages());
	  dummy = z_field(Range(joker),row,col);
	  ostringstream os;
	  os << "z_field (for latitude nr " << row << " and longitude nr " 
             << col << ")";
	  chk_if_increasing( os.str(), dummy ); 
	}
    }

  // Check that there is no gap between the ground and lowest pressure surface
  //
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
		if( atmosphere_dim > 1 )
		  os << "\nThis was found to be the case for:\n" 
		     << "latitude " << lat_grid[row];
		if( atmosphere_dim > 2 )
		  os << "\nlongitude " << lon_grid[col];
	      throw runtime_error( os.str() );
	    }
	}
    }

  // Blackbody ground and cloud box
  //
  chk_if_bool(  "blackbody_ground", blackbody_ground );
  chk_cloudbox( atmosphere_dim, p_grid, lat_grid, lon_grid, blackbody_ground, 
		cloudbox_on, cloudbox_limits );

  // Sensor position and LOS
  //
  chk_vector_length( "sensor_pos", sensor_pos, atmosphere_dim );
  chk_if_over_0( "sensor radius", sensor_pos[0] );
  if( atmosphere_dim <= 2 )
    chk_vector_length( "sensor_los", sensor_los, 1 );
  else
    {
      chk_if_in_range( "sensor latitude", sensor_pos[1], -90, 90 );
      chk_if_in_range( "sensor longitude", sensor_pos[2], -360, 360 );
      chk_vector_length( "sensor_los", sensor_los, 2 );
      chk_if_in_range( "sensor azimuth angle", sensor_los[1], -180, 180 );
    }
  chk_if_in_range( "sensor_los zenith angle", sensor_los[0], -180, 180 );
  
  //--- End: Check input ------------------------------------------------------


  //
  //

}


