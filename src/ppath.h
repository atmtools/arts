/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>

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



///////////////////////////////////////////////////////////////////////////////
//   File description
///////////////////////////////////////////////////////////////////////////////
/*!
  \file   ppath.h
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   Sun May  5  2002
  
  \brief  Propagation path structure and functions.
  
   This file contains the definition of the Ppath structure and the
   functions in ppath.c that are of interest elsewhere.
*/


#ifndef ppath_h
#define ppath_h


#include "arts.h"
#include "interpolation.h"
#include "matpackI.h"
#include "mystring.h"



///////////////////////////////////////////////////////////////////////////////
//   The Ppath structure
///////////////////////////////////////////////////////////////////////////////

struct Ppath {
  Index             dim;
  Index             np;
  Matrix    	    pos;
  Vector    	    z;
  Vector            l_step;
  ArrayOfGridPos    gp_p;
  ArrayOfGridPos    gp_lat;
  ArrayOfGridPos    gp_lon;
  Matrix    	    los;
  String    	    background;
  Index     	    ground;
  Index     	    i_ground;
  Vector    	    tan_pos;
  Index     	    symmetry;
  Index     	    i_symmetry;
};



///////////////////////////////////////////////////////////////////////////////
//   Functions from ppath.cc
///////////////////////////////////////////////////////////////////////////////

void ppath_init_structure( 
	      Ppath&      ppath,
	const Index&      atmosphere_dim,
	const Index&      np );

void ppath_set_background( 
	      Ppath&      ppath,
        const Index&      case_nr );

Index ppath_what_background( 
	      Ppath&      ppath );

void ppath_start_stepping(
              Ppath&          ppath,
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
	const Vector&         sensor_los );

void ppath_step_1d_geom(
	      Ppath&      ppath,
        const Index&      atmosphere_dim,
        ConstVectorView   p_grid,
        ConstVectorView   z_grid,
        const Numeric&    r_geoid,
        const Numeric&    z_ground,
        const Index&      blackbody_ground,
	const Numeric&    lmax );


#endif  // ppath_h
