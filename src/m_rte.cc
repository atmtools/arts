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
              Vector&         i_space,
        // WS Input:
	const Vector&         f_grid )
{
  const Index n = f_grid.nelem();
  if( n == 0 )
    throw runtime_error( "The frequency grid is empty." );
  i_space.resize(n);
  out2 << "  Setting i_space to hold cosmic background radiation.\n";
  for( Index i=0; i<n; i++ )
    { planck( i_space[i], f_grid[i], COSMIC_BG_TEMP ); }
}



//! RteCalc
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-29
*/
void RteCalc(
        // WS Output:
              Vector&         y,
        // WS Input:
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Matrix&         t_ground,
        const Tensor3&        e_ground,
        const Index&          blackbody_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         f_grid,
        const Vector&         i_space,
	const Matrix&         abs,
        const Index&          antenna_dim,
        const Vector&         mblock_za_grid,
        const Vector&         mblock_aa_grid,
	const Index&          stokes_dim,
        const Matrix&         sensor_pos,
        const Matrix&         sensor_los )
{
  /* Notes for the RTE calculations:

  rte_background_agenda
    ---
    Needed as to handle case with or without cloudbox/emission and to avoid
    that too many variables must be set to dummy values
    ---
    Needs   : ppath, ground variables, i_space, ...
    Creates : i_rte_background
    ---
    Versions: RteBackgroundEmission, RteBackgroundEmissionWithCloudBox

  rte_calculate_agenda
    ---
    To switch between considering emission or not. The cloudbox is not
    important here.
    ---
    Needs   : ppath, i_background, atmospheric variables, ...
    Creates : i_rte
    ----
    Versions: RteEmission

   
  The main functions for the calculations is RteCalc
  */



  // Some sizes
  //
  const Index nf      = f_grid.nelem();
  const Index nmblock = sensor_pos.nrows();
  const Index nza     = mblock_za_grid.nelem();


  //--- Check input -----------------------------------------------------------

  // There are extensives checks performed in ppathCalc and those are not
  // duplicated here.

  // Atmospheric variables
  //
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );

  // Ground variables
  //
  chk_atm_surface( "t_ground", t_ground, atmosphere_dim, lat_grid, lon_grid );
  if( e_ground.npages() != nf )
    {
      ostringstream os;
      os << "The number of pages of e_ground must be equal to the number of "
	 << "frequencies, but f_grid has the length " << nf <<",\n"
	 << "while e_ground has " << e_ground.npages() << " pages.";
      throw runtime_error( os.str() );
    }

  chk_atm_surface( "a page of e_ground", e_ground(0,Range(joker),Range(joker)),
                                          atmosphere_dim, lat_grid, lon_grid );
  
  // Frequency grid and i_space
  //
  if( nf == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  chk_vector_length( "f_grid", "i_space", f_grid, i_space );

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

  // Sensor position and LOS. That the angles are inside OK ranges are checked
  // inside ppathCalc.
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


  // More sizes
  //
  Index naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }
  
  // Matrix all the return from RTE function. 
  // This matrix must hold all Stokes elements.
  //
  Matrix i_total( nmblock * nf * nza * naa, stokes_dim );

  // Variables needed below
  //
  Ppath  ppath;
  Vector los( sensor_los.ncols() );   // LOS for a MPB direction
  Index  nfdone=0;                    // Number of frequencies done


  // Loop:  measurement block / azimuthal angle / zenith angle
  //
  for( Index imblock=0; imblock<nmblock; imblock++ )
    {
      for( Index iza=0; iza<nza; iza++ )
	{
	  for( Index iaa=0; iaa<naa; iaa++ )
	    {
	      // LOS
	      los     = sensor_los( imblock, Range(joker) );
	      los[0] += mblock_za_grid[iza];
	      if( antenna_dim == 2 )
		{ los[1] += mblock_aa_grid[iaa]; }

	      // Determine propagation path
	      ppathCalc( ppath, atmosphere_dim, p_grid, lat_grid, lon_grid, 
		    z_field, r_geoid, z_ground, blackbody_ground, cloudbox_on,
                      cloudbox_limits, sensor_pos(imblock,Range(joker)), los );

	      // Calculate spectra and WFs
	      // (this will be an agenda)
	      Matrix i_rte;
	      //rte_emission( i_rte, 
	      

	      // Copy i_rte to i_total
	      //
	      assert( i_rte.nrows() == nf );
	      assert( i_rte.ncols() == stokes_dim );
	      //
	      i_total(Range(nfdone,nf),Range(joker)) = i_rte;

	      // Increase nfdone
	      nfdone += nf;
	    }
	}
    } 


  // Apply the polarisation of the sensor
  // (so far a temporary solution)
  //
  y.resize( nmblock * nf * nza * naa );
  //
  y = i_total( Range(joker), 0 );
}
