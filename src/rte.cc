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

#include <stdexcept>
#include "check_input.h"
#include "logic.h"
#include "physics_funcs.h"
#include "rte.h"



/*===========================================================================
  === 
  ===========================================================================*/

//! ground_reflection_with_emission
/*!
    Updates a vector with radiances (normally *i_rte*) with the effect of a
    ground reflection.

    The assumptions for the functions shall be listed here. What are they?
    [* PE 2002-08-20 *]

    The function takes *ppath* as input to be able to determine the position
    and LOS for the ground reflection. 

    The input variables beside *ppath* are the variables needed for consistency
    checks and to determine the ground temperature and emissivity.

    The matrix *i_rte* is the incoming radiation (for all Stokes components), 
    the radiation hitting the ground at the point of the reflection.
    The size of *i_rte* must be [f_grid.nelem(),stokes_dim].

    The side effects of the function are: <br>
    1. The WSV *t_ground* and *e_ground* are set by calling the corresponding
       agendas. <br>
    2. The WSVs *a_pos* and *a_los* are set to the position and LOS, 
       respectively, for the grounmd reflection.

    \param   i_rte              Output: As the WSV with the same name.
    \param   t_ground           Output: As the WSV with the same name.
    \param   e_ground           Output: As the WSV with the same name.
    \param   a_pos              Output: As the WSV with the same name.
    \param   a_los              Output: As the WSV with the same name.
    \param   t_ground_agenda    As the WSV with the same name.
    \param   e_ground_agenda    As the WSV with the same name.
    \param   blackbody_ground   As the WSV with the same name.
    \param   ppath              As the WSV with the same name.
    \param   f_grid             As the WSV with the same name.
    \param   stokes_dim         As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-08-20
*/
void ground_reflection_with_emission(
              Matrix&         i_rte,
              Numeric&        t_ground,
	      Matrix&         e_ground,
	      Vector&         a_pos,
	      Vector&         a_los,
	const Agenda&         t_ground_agenda,
	const Agenda&         e_ground_agenda,
	const Index&          blackbody_ground,
	const Ppath&          ppath,
        const Vector&         f_grid,
	const Index&          stokes_dim )
{
  assert( !blackbody_ground );
  assert( ppath.ground );
	      
  // Set a_pos and a_los, that are used of the agendas to call
  a_pos.resize( ppath.pos.ncols() );
  a_pos = ppath.pos( ppath.i_ground, Range(joker) );
  a_los.resize( ppath.los.ncols() );
  a_los = ppath.los( ppath.i_ground, Range(joker) );
  
  // Determine t_ground and e_ground
  chk_not_empty( "t_ground_agenda", t_ground_agenda );
  chk_not_empty( "e_ground_agenda", e_ground_agenda );
  t_ground_agenda.execute();
  e_ground_agenda.execute();

  assert( e_ground.nrows() == f_grid.nelem() );
  assert( e_ground.ncols() > 0 );

  // The matrix *e_ground* has been set to have one column. This means
  // that the ground reflection does not change the polarisation.
  Numeric e;
  //
  if( e_ground.ncols() == 1 )
    {
      for( Index iv=0; iv<e_ground.nrows(); iv++ )
	{
	  e = e_ground(iv,0);
	  assert( e >= 0  &&  e <= 1 );
	  i_rte(iv,0) = i_rte(iv,0) * ( 1 - e ) +
                                            e * planck( f_grid[iv], t_ground );
	}
    }

  // Other cases are not yet handled
  else
    {
      throw runtime_error( 
                   "Only unpolarised ground reflections are handled so far." );
    }
}



//! radiation_from_blackbody_ground
/*!
    Sets the WSV *i_rte* to hold radiation from a blackbody ground.

    The radiation is assumed to be unpolarised. That is, only the first
    Stokes component will be non-zero. The radiation equals the Planck
    function for *t_ground*.

    The size of *i_rte* must be set to [f_grid.nelem(),stokes_dim] before
    calling this function.

    \param   i_rte              Output: As the WSV with the same name.
    \param   t_ground_agenda    As the WSV with the same name.
    \param   blackbody_ground   As the WSV with the same name.
    \param   f_grid             As the WSV with the same name.
    \param   stokes_dim         As the WSV with the same name.
    \param   t_ground           As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-08-20
*/
void radiation_from_blackbody_ground(
              Matrix&         i_rte,
	const Agenda&         t_ground_agenda,
	const Index&          blackbody_ground,
        const Vector&         f_grid,
	const Index&          stokes_dim,
        const Numeric&        t_ground )
{
  // Number of frequencies
  const Index nf = f_grid.nelem();
  
  // Asserts
  assert( is_size( i_rte, nf, stokes_dim ) );
  assert( blackbody_ground );

  // Determine effective emission temperature of the ground
  chk_not_empty( "t_ground_agenda", t_ground_agenda );
  t_ground_agenda.execute();

  // Set i_rte to un-polarised blackbody radiation
  //
  i_rte = 0;
  //
  for( Index iv=0; iv<nf; iv++ )
    { i_rte(iv,0) = planck( f_grid[iv], t_ground ); }
}



//! set_to_radiative_background
/*!
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
    \param   t_ground_agenda    As the WSV with the same name.
    \param   blackbody_ground   As the WSV with the same name.
    \param   ppath              As the WSV with the same name.
    \param   f_grid             As the WSV with the same name.
    \param   stokes_dim         As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-08-20
*/
void set_to_radiative_background(
              Matrix&         i_rte,
              Matrix&         i_space,
	      Vector&         a_pos,
	      Vector&         a_los,
              Numeric&        t_ground,
	const Agenda&         i_space_agenda,
	const Agenda&         t_ground_agenda,
	const Index&          blackbody_ground,
	const Ppath&          ppath,
        const Vector&         f_grid,
	const Index&          stokes_dim )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Resize i_rte to have the correct the size
  i_rte.resize( nf, stokes_dim );

  // Set a_pos and a_los to match the last point in ppath
  a_pos.resize( ppath.pos.ncols() );
  a_pos = ppath.pos(np-1,Range(joker));
  a_los.resize( ppath.los.ncols() );
  a_los = ppath.los(np-1,Range(joker));

  // Initialize i_rte to the radiative background
  switch ( ppath_what_background( ppath ) )
    {

    case 1:   //--- Space ---------------------------------------------------- 
      //
      chk_not_empty( "i_space_agenda", i_space_agenda );
      i_space_agenda.execute();
      i_rte = i_space;
      break;

    case 2:   //--- Blackbody ground -----------------------------------------
      //
      chk_not_empty( "t_ground_agenda", t_ground_agenda );
      radiation_from_blackbody_ground( i_rte, t_ground_agenda, 
                             blackbody_ground, f_grid, stokes_dim, t_ground );
      break;

    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }
}
