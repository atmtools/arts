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
#include "logic.h"
#include "physics_funcs.h"
#include "rte.h"



/*===========================================================================
  === 
  ===========================================================================*/

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
	const Index&          ip,
        const Vector&         f_grid,
	const Index&          stokes_dim )
{
  assert( !blackbody_ground );
	      
  // Set a_pos and a_los, that are used of the agendas to call
  a_pos.resize( ppath.pos.ncols() );
  a_pos = ppath.pos(ip,Range(joker));
  a_los.resize( ppath.los.ncols() );
  a_los = ppath.los(ip,Range(joker));
  
  // Determine t_ground and e_ground
  t_ground_agenda.execute();
  e_ground_agenda.execute();

  assert( e_ground.nrows() == f_grid.nelem() );
  assert( e_ground.ncols() > 0 );

  // The matrix *e_ground* has been set to have one column. This means
  // that the ground reflection does not change the polarisation.
  if( e_ground.ncols() == 1 )
    {
      for( Index iv=0; iv<e_ground.nrows(); iv++ )
	{
	  i_rte(iv,0) = i_rte(iv,0) * ( 1 - e_ground(iv,0) ) +
                              e_ground(iv,0) * planck( f_grid[iv], t_ground );
	}
    }

  // Other cases are not yet handled
  else
    {
      throw runtime_error( 
                   "Only unpolarised ground reflections are handled so far." );
    }
}



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

  // Check that *blackbody_ground* in fact is 1.
  if( !blackbody_ground )
    {
      ostringstream os;
      os << "Propagation path is calculated assuming a blackbody ground,\n"
	 << "but *blackbody_ground* is now turned off.\n";
      throw runtime_error( os.str() );
    }

  // Determine effective emission temperature of the ground
  t_ground_agenda.execute();

  // Set i_rte to un-polarised blackbody radiation
  //
  i_rte = 0;
  //
  for( Index iv=0; iv<nf; iv++ )
    { i_rte(iv,0) = planck( f_grid[iv], t_ground ); }
}



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
      i_space_agenda.execute();
      i_rte = i_space;
      break;

    case 2:   //--- Blackbody ground -----------------------------------------
      //
      radiation_from_blackbody_ground( i_rte, t_ground_agenda, 
                             blackbody_ground, f_grid, stokes_dim, t_ground );
      break;

    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }
}
