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



/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
  \file   rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-29

  \brief  Functions to solve radiative transfer tasks.
*/



/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include "rte.h"



/*****************************************************************************
 *** High-level interpolation functions
 *****************************************************************************/



/*****************************************************************************
 *** Functions to determine background
 *****************************************************************************/

//void rte_background_emission()
//{
//}



/*****************************************************************************
 *** Core functions for RTE
 *****************************************************************************/


void rte_emission(
              Matrix&         i,
        const Ppath&          ppath, 
        const Index&          atmosphere_dim,
        ConstTensor3View      z_field,
        ConstTensor3View      t_field,
        ConstMatrixView       t_ground,
        ConstTensor3View      e_ground,
        ConstVectorView       i_space,
	ConstMatrixView       abs,
	const Index&          stokes_dim )
{
  // Sizes
  const Index   np   = z_field.npages();
  const Index   nlat = z_field.ncols();
  const Index   nlon = z_field.nrows();
  const Index   nf   = e_ground.npages();



}
