/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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



/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   m_rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-11 

  \brief  Workspace functions for solving clear sky radiative transfer.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "messages.h"

extern const Numeric COSMIC_BG_TEMP;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


//! i_spaceCBR
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void i_spaceCBR(
        // WS Output:
              Matrix&         i_space,
        // WS Input:
	      const Vector&   f_grid,
              const Index&    stokes_dim )
{
  const Index n = f_grid.nelem();

  if( n == 0 )
    throw runtime_error( "The frequency grid is empty." );

  out2 << "  Setting *i_space* to hold cosmic background radiation.\n";

  i_space.resize(n,stokes_dim);
  i_space = 0;

  for( Index i=0; i<n; i++ )
    { i_space(i,0) = planck( f_grid[i], COSMIC_BG_TEMP ); }
}



//! RteCalc
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-09
*/
void RteCalc(
        // WS Output:
              Matrix&         y_rte,
	      Ppath&          ppath,
	      Ppath&          ppath_step,
	      Matrix&         i_rte,
        // WS Input:
	const Agenda&         ppath_step_agenda,
	const Agenda&         rte_agenda,
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
        const Matrix&         sensor_pos,
	const Matrix&         sensor_los,
        const Vector&         f_grid,
	const Index&          stokes_dim,
        const Index&          antenna_dim,
        const Vector&         mblock_za_grid,
	const Vector&         mblock_aa_grid )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index nmblock = sensor_pos.nrows();
  const Index nza     = mblock_za_grid.nelem();


  //--- Check input -----------------------------------------------------------

  // Basic checks of atmospheric, geoid and ground variables
  //  
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_ground", z_ground, atmosphere_dim, lat_grid, lon_grid );

  // Check that z_field has strictly increasing pages.
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

  // Frequency grid
  //
  if( nf == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );

  // Antenna
  //
  chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
  if( nza == 0 )
    throw runtime_error( "The measurement block zenith angle grid is empty." );
  chk_if_increasing( "mblock_za_grid", mblock_za_grid );
  if( antenna_dim == 1 )
    {
      if( mblock_aa_grid.nelem() != 0 )
	throw runtime_error( 
	      "For antenna_dim = 1, the azimuthal angle grid must be empty." );
    }
  else
    {
      if( mblock_aa_grid.nelem() == 0 )
	throw runtime_error( 
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }

  // Stokes
  //
  chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

  // Sensor position and LOS. 
  //
  // That the angles are inside OK ranges are checked inside ppathCalc.
  //
  if( sensor_pos.ncols() != atmosphere_dim )
    throw runtime_error( "The number of columns of sensor_pos must be equal "
                                         "to the atmospheric dimensionality" );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
  if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
    throw runtime_error( "For 3D, sensor_los shall have two columns." );
  if( sensor_los.nrows() != nmblock )
    {
      ostringstream os;
      os << "The number of rows of sensor_pos and sensor_los must be "
	 << "identical, but sensor_pos has " << nmblock << " rows,\n"
	 << "while sensor_los has " << sensor_los.nrows() << " rows.";
      throw runtime_error( os.str() );
    }

  //--- End: Check input ------------------------------------------------------


  // Number of azimuthal direction for pencil beam calculations
  Index naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }
  
  // Resize *y_rte* to have the correct size.
  y_rte.resize( nmblock * nf * nza * naa, stokes_dim );

  // Variables needed below
  Vector los( sensor_los.ncols() );   // LOS for a MPB calculation
  Index  nfdone=0;                    // Number of frequencies done


  // Loop:  measurement block / zenith angle / azimuthal angle
  //
  for( Index imblock=0; imblock<nmblock; imblock++ )
    {
      for( Index iza=0; iza<nza; iza++ )
	{
	  for( Index iaa=0; iaa<naa; iaa++ )
	    {
	      // LOS of interest
	      los     = sensor_los( imblock, Range(joker) );
	      los[0] += mblock_za_grid[iza];
	      if( antenna_dim == 2 )
		{ los[1] += mblock_aa_grid[iaa]; }

	      // Determine propagation path
	      ppathCalc( ppath, ppath_step, ppath_step_agenda, atmosphere_dim, 
		        p_grid, lat_grid, lon_grid, z_field, r_geoid, z_ground,
			    blackbody_ground, cloudbox_on, cloudbox_limits, 
			               sensor_pos(imblock,Range(joker)), los );

	      // Execute the *rte_agenda*
	      rte_agenda.execute();
	      
	      // Copy i_rte to y_rte
	      y_rte(Range(nfdone,nf),Range(joker)) = i_rte;

	      // Increase nfdone
	      nfdone += nf;
	    }
	}
    } 
}



//! RteEmissionStd
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-16
*/
void RteEmissionStd(
        // WS Output:
              Matrix&         i_rte,
	      Matrix&         i_space,
        // WS Input:
	const Agenda&         i_space_agenda,
	const Ppath&          ppath,
        const Vector&         f_grid,
        const Index&          stokes_dim )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Resize i_rte to have the correct the size
  i_rte.resize( nf, stokes_dim );

  // Initialize i_rte to the radiative background
  switch ( ppath_what_background( ppath ) )
    {
    case 1:   //--- Space ---------------------------------------------------- 
      i_space_agenda.execute();
      i_rte = i_space;
      break;
    default:
      {
	ostringstream os;
	os << "Unknown radiative background (code nr. " 
	   << ppath_what_background( ppath ) << ")";
	throw runtime_error( os.str() );
      }
    }

  // If the number of propagation path points is > 2, we are already ready,
  // the observed spectrum equals the radiative background.
  if( np >= 2 )
    {

      // Determine the atmospheric temperature at each propagation path point
      // As a temporary solution, the temperature is set to 200 K.
      Vector t_ppath(np,200);
	 
      // Determine the total absorption at each propagation path point
      // As a temporary solution, the absorption is set to 1-e6 1/m.
      Matrix abs_ppath(nf,np,1e-6);

      // Some variables used below
      Index     iv;      // Frequency index
      Numeric   trans;   // Transmission for a path step for one frequency
      Numeric   source;  // Effective source value for one step and 1 frequency

      // Loop the propagation path steps
      //
      // The number of path steps is np-1.
      // The path points are stored in such way that index 0 corresponds to
      // the point closest to the sensor.
      //
      for( Index ip=np-1; ip>0; ip-- )
	{
	  // A ground reflection at point ip?
	  if( ppath.ground  && ppath.i_ground == ip )
	    {
	      throw runtime_error("Ground reflections are not yet handled.");
	    }
	  
	  // Consider absorption and emission from point ip to ip-1 for
	  // all frequencies
	  for( iv=0; iv<nf; iv++ )
	    {
	      // Transmission along the path for the step
	      trans = exp( -ppath.l_step[ip-1] * 
			       ( abs_ppath(iv,ip) + abs_ppath(iv,ip-1) ) / 2 );

	      // The Planck function for the mean of the end temperatures
	      source = planck( f_grid[iv], (t_ppath[ip]+t_ppath[ip-1])/2 );
	      
	      // Go from point ip to point ip-1 using the constant source term
	      // approximation.
	      i_rte(iv,0) = i_rte(iv,0) * trans + source * ( 1 - trans );
	    }
	}
    }
}



//! VectorToTbByRJ
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-09
*/
void VectorToTbByRJ(
        // WS Generic Output:
              Vector&   y_out,
        // WS Generic Output Names:
        const String&   y_out_name,
        // WS Generic Input:
        const Vector&   y_in,
        const Vector&   y_in_name,
        // WS Generic Input Names:
        const String&   f_grid,
        const String&   f_grid_name )
{
  // Length of y_in
  const Index nin = y_in.nelem();

  throw runtime_error("The implementation of this function is not yet ready");

  out2 << "   " << y_out_name << " = inv_of_rj(" << y_in_name << "," 
       << f_grid_name << ")\n" ;

  // Note that y_in and y_out can be the same vector
  if ( &y_out != &y_in )
    { y_out.resize(nin); }

  
}



//! yNoPolarisation
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-09
*/
void yNoPolarisation(
        // WS Output:
              Vector&         y,
        // WS Input:
        const Matrix&         y_rte )
{
  y.resize( y_rte.nrows() );

  y = y_rte( Range(joker), 0 );
}
