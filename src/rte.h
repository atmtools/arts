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
  \file   rte.h
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-29

  \brief  Declaration of functions in rte.cc.
*/



#ifndef rte_h
#define rte_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "agenda_class.h"
#include "arts.h"
#include "ppath.h"
#include "matpackI.h"
#include "matpackIII.h"



/*===========================================================================
  === Functions in rte.cc
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
	const Index&          stokes_dim );

void radiation_from_blackbody_ground(
              Matrix&         i_rte,
	const Agenda&         t_ground_agenda,
	const Index&          blackbody_ground,
        const Vector&         f_grid,
	const Index&          stokes_dim,
	const Numeric&        t_ground );

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
	const Index&          stokes_dim );



#endif  // rte_h
