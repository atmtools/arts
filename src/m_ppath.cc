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

void ppathCalc(// WS Output:
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
  // Check input
  

}


