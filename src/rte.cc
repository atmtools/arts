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
  \file   rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-29

  \brief  Functions to solve radiative transfer tasks.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "auto_md.h"
#include "check_input.h"
#include "logic.h"
#include "physics_funcs.h"
#include "rte.h"



/*===========================================================================
  === 
  ===========================================================================*/


//! rte_step_clearsky_with_emission
/*!
    Performs monochromatic radiative transfer for an atmospheric slab 
    with constant conditions, no scattering and having emission.

    The function is best explained by considering a homogenous layer. That is,
    the physical conditions inside the layer are constant. The absorption
    inside the layer is described by *ext_mat_gas* and *abs_vec_gas*,
    the blackbdoy radiation of the layer is given by *a_planck_value*
    and the propagation path length through the layer is *l_step*. 

    When calling the function, the vector *stokes_vec* shall contain the
    Stopkes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the 
    layer.

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. Scalar absorption (stokes_dim = 1). <br>
       2. The matrix ext_mat_gas is diagonal (unpolarised absorption). <br>
       3. The total general case.

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   stokes_dim         Input: As the WSV with the same name.
    \param   ext_mat_gas        Input: As the WSV with the same name.
    \param   abs_vec_gas        Input: As the WSV with the same name.
    \param   l_step             Input: The length of the RTE step.
    \param   a_planck_value     Input: Blackbody radiation.

    \author Patrick Eriksson 
    \date   2002-09-16
*/
void rte_step_clearsky_with_emission(
	      Vector&       stokes_vec,		       
        const Index&        stokes_dim,
	ConstMatrixView&    ext_mat_gas,
	ConstVectorView&    abs_vec_gas,
	const Numeric&      l_step,
	const Numeric&      a_planck_value )
{
  // Asserts
  assert( stokes_dim >= 1  &  stokes_dim <= 4 );
  assert( stokes_vec.nelem() == stokes_dim );
  assert( ext_mat_gas.nrows() == stokes_dim );
  assert( ext_mat_gas.ncols() == stokes_dim );
  assert( abs_vec_gas.nelem() == stokes_dim );
  assert( a_planck_value >= 0 );

  // Scalar case
  if( stokes_dim == 1 )
    { 
      assert( ext_mat_gas(0,0) == abs_vec_gas[0] );
      Numeric transm = exp( -l_step * abs_vec_gas[0] );
      stokes_vec[0] = stokes_vec[0] * transm + ( 1- transm ) * a_planck_value;
    }

  // Vector case
  else
    {
      // We have here two cases, diagonal or non-diagonal ext_mat_gas
      // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
      // non-zero value in position 1.
      // 
      bool diagonal = true;
      //
      for( Index i=1; diagonal  && i<stokes_dim; i++ )
	{
	  for( Index j=0; diagonal && j<i; j++ )
	    {
	      if( ext_mat_gas(i,j) != 0  ||  ext_mat_gas(j,i) )
		{ diagonal = false; }
	    }
	  assert( !diagonal  ||  ( diagonal  && abs_vec_gas[i] == 0 ) );
	}
      cout << "Diagonal = " << diagonal << "\n";

      // Unpolarised
      if( diagonal )
	{
	  // Stokes dim 1
	  assert( ext_mat_gas(0,0) == abs_vec_gas[0] );
	  Numeric transm = exp( -l_step * abs_vec_gas[0] );
	  stokes_vec[0] = stokes_vec[0] * transm + 
                                                ( 1- transm ) * a_planck_value;

	  // Stokes dims > 1
	  for( Index i=1; i<stokes_dim; i++ )
	    { stokes_vec[i] *= exp( -l_step * ext_mat_gas(i,i) ); }
	}


      // Polarised
      else
	{
	  // Here we use the general method for the cloud box.
          stokes_vecGeneral( stokes_vec, ext_mat_gas, abs_vec_gas, 
		    Vector(stokes_dim,0), l_step, a_planck_value, stokes_dim );
	}
    }
}



//! get_radiative_background
/*!
    This text is old. PE 2002-09-17

    Sets a vector (normally *i_rte*) to the radiative background for a 
    propagation path.

    The function uses *ppath* to determine the radiative background for a 
    propagation path. The possible backgrounds are listed in the header of
    the function ppath_set_background (in ppath.cc).

    The input variables beside *ppath* are the variables needed to calculate
    the magnitude of the radiative background. 

    The main purpose of the function is set *i_rte*. It is NOT needed to set 
    *i_rte* to the correct size before calling the function. The size is set
    to be [f_grid.nelem(),stokes_dim].

    The side effects of the function are: <br>
    1. The WSV *i_space* is changed if the radiative background is the 
       space. <br>
    2. The WSVs *a_pos* and *a_los* are always set to the position and LOS of 
       the last point in ppath (point np-1).

    \param   i_rte              Output: As the WSV with the same name.
    \param   i_space            Output: As the WSV with the same name.
    \param   a_pos              Output: As the WSV with the same name.
    \param   a_los              Output: As the WSV with the same name.
    \param   i_space_agenda     As the WSV with the same name.
    \param   ppath              As the WSV with the same name.
    \param   f_grid             As the WSV with the same name.
    \param   stokes_dim         As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-09-17
*/
void get_radiative_background(
              Matrix&         i_rte,
	      Ppath&          ppath_step,
	      Vector&         a_pos,
	      Vector&         a_los,
	      GridPos&        a_gp_p,
	      GridPos&        a_gp_lat,
	      GridPos&        a_gp_lon,
              Matrix&         i_space,
              Matrix&         ground_emission, 
              Matrix&         ground_los, 
	      Tensor4&        ground_refl_coeffs,
	const Ppath&          ppath,
	const Index&          mblock_index,
	const Agenda&         ppath_step_agenda,
	const Agenda&         rte_agenda,
	const Agenda&         i_space_agenda,
	const Agenda&         ground_refl_agenda,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        const Tensor3&        z_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Tensor7&        scat_i_p,
        const Tensor7&        scat_i_lat,
        const Tensor7&        scat_i_lon,
        ConstVectorView       scat_za_grid,
        ConstVectorView       scat_aa_grid,
        ConstVectorView       f_grid,
	const Index&          stokes_dim,
	const Index&          antenna_dim )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Resize i_rte to have the correct the size
  i_rte.resize( nf, stokes_dim );

  // Set a_pos, a_gp_XXX and a_los to match the last point in ppath.
  // The agendas below use different combinations of these varaibles.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  a_pos.resize( atmosphere_dim );
  a_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  a_los.resize( ppath.los.ncols() );
  a_los = ppath.los(np-1,Range(joker));
  gridpos_copy( a_gp_p, ppath.gp_p[np-1] );
  if( atmosphere_dim > 1 )
    { gridpos_copy( a_gp_lat, ppath.gp_lat[np-1] ); }
  if( atmosphere_dim > 2 )
    { gridpos_copy( a_gp_lon, ppath.gp_lon[np-1] ); }


  // Initialize i_rte to the radiative background
  switch ( ppath_what_background( ppath ) )
    {

    case 1:   //--- Space ---------------------------------------------------- 
      {
	chk_not_empty( "i_space_agenda", i_space_agenda );
	i_space_agenda.execute();

	if( i_space.nrows() != nf  ||  i_space.ncols() != stokes_dim )
	  throw runtime_error( "The size of *i_space* returned from "
                                           "*i_space_agenda* is not correct.");

	i_rte = i_space;
      }
      break;


    case 2:   //--- The ground -----------------------------------------------
      {
	chk_not_empty( "ground_refl_agenda", ground_refl_agenda );
	ground_refl_agenda.execute();

	if( ground_emission.nrows() != nf  ||  
                                        ground_emission.ncols() != stokes_dim )
	  throw runtime_error(
		  "The size of the created *ground_emission* is not correct.");

	i_rte = ground_emission;

	// Calculate the spectra hitting the ground
	//
	// Some local variables must be used not to over-write WSVs that must
	// be preserved.
	Matrix     y_rte_local;            
	Matrix     i_rte_local;            
	Ppath      ppath_local;
	Index      mblock_index_local;
	Matrix     ground_emission_local;       
	Matrix     ground_los_local;
	Tensor4    ground_refl_coeffs_local;
	// Set zenith and azimuthal angles (note that these variables are not
	// function input). Each ground los is here treated as measurement
	// block, with no za and aa grids.
	Index    nlos = ground_los.nrows();
	Matrix   sensor_pos( nlos, a_pos.nelem() );
	Matrix   sensor_los( nlos, a_los.nelem() );
	for( Index ilos=0; ilos < nlos; ilos++ )
	  {
	    sensor_pos( ilos, Range(joker) ) = a_pos;
	    sensor_los( ilos, Range(joker) ) = ground_los( ilos, Range(joker));
	  }
	Vector   mblock_za_grid(1,0);
	Vector   mblock_aa_grid(0);
	if( antenna_dim > 1 )
	  {
	    mblock_aa_grid.resize(1);
	    mblock_aa_grid = 0;
	  }
	//
	RteCalc( y_rte_local, ppath_local, ppath_step, i_rte_local, 
		 mblock_index_local, a_pos, a_los, a_gp_p, a_gp_lat, a_gp_lon,
		 i_space, ground_emission_local, ground_los_local, 
		 ground_refl_coeffs_local, ppath_step_agenda, rte_agenda, 
		 i_space_agenda, ground_refl_agenda, 
		 atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, r_geoid, 
		 z_ground, cloudbox_on, cloudbox_limits, scat_i_p, scat_i_lat,
		 scat_i_lon, scat_za_grid, scat_aa_grid, sensor_pos, 
		 sensor_los, f_grid, stokes_dim, antenna_dim, 
		 mblock_za_grid, mblock_aa_grid );

	// Decompose the calculated spectra (y_rte_local) and add the values
	// to i_rte, considering the ground reflection coeff. matrix.
	for( Index ilos=0; ilos < nlos; ilos++ )
	  {
	    Matrix i_ground(nf,stokes_dim);

	    i_ground = y_rte_local(Range(ilos*nf,nf),Range(joker));
	    
	    for( Index iv=0; iv<nf; iv++ )
	      {
		if( stokes_dim == 1 )
		  {
		    i_rte(iv,0) += ground_refl_coeffs(ilos,iv,0,0) * 
                                                                i_ground(iv,0);
		  }
		else
		  {
		    Vector stokes_vec(stokes_dim);
		    mult( stokes_vec, 
                       ground_refl_coeffs(ilos,iv,Range(joker),Range(joker)), 
                                                   i_ground(iv,Range(joker)) );
		    for( Index is=0; is < stokes_dim; is++ )
		      { i_rte(iv,is) += stokes_vec[is]; }
		  }
	      }
	  }
      }
      break;


    case 3:   //--- Cloudbox surface -----------------------------------------

      {
	CloudboxGetOutgoing( i_rte, "i_rte", scat_i_p, scat_i_lat, scat_i_lon, 
			   a_gp_p, a_gp_lat, a_gp_lon, a_los, 
                           cloudbox_on, cloudbox_limits, atmosphere_dim,
                              stokes_dim, scat_za_grid, scat_aa_grid, f_grid );
      }
      break;


    case 4:   //--- Inside of cloudbox ---------------------------------------
      {
	throw runtime_error( "RTE calculations starting inside the cloud box "
			                               "are not yet handled" );
      }
      break;


    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }
}



