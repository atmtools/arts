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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   check_input.h

   This file contains the declaration of functions in check_input.cc.

   \author Patrick Eriksson
   \date 2002-04-15 
*/



#ifndef checkinput_h
#define checkinput_h


#include "arts.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "mystring.h"



void chk_if_bool( 
        const String&   x_name,
        const Index&    x );

void chk_if_in_range( 
	const String&   x_name,
        const Index&    x, 
        const Index&    x_low, 
        const Index&    x_high );

void chk_if_over_0( 
	const String&    x_name,
        const Numeric&   x );

void chk_if_in_range( 
	const String&    x_name,
        const Numeric&   x, 
        const Numeric&   x_low, 
        const Numeric&   x_high );

void chk_vector_length( 
	const String&      x_name,
        ConstVectorView    x,
        const Index&       l );

void chk_vector_length( 
	const String&      x1_name,
	const String&      x2_name,
        ConstVectorView    x1, 
        ConstVectorView    x2 );

void chk_if_increasing( 
	const String&      x_name,
        ConstVectorView    x );

void chk_if_decreasing( 
	const String&      x_name,
        ConstVectorView    x );

void chk_matrix_ncols( 
	const String&      x_name,
        ConstMatrixView    x,
        const Index&       l );

void chk_matrix_nrows( 
	const String&      x_name,
        ConstMatrixView    x,
        const Index&       l );

void chk_atm_grids( 
	const Index&      dim,
	ConstVectorView   p_grid,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid );

void chk_atm_field( 
	const String&     x_name,
        const Tensor3&    x, 
	const Index&      dim,
	ConstVectorView   p_grid,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid );

void chk_atm_surface( 
	const String&     x_name,
        const Matrix&     x, 
	const Index&      dim,
	ConstVectorView   alpha_grid,
	ConstVectorView   beta_grid );

void chk_cloudbox(
	const Index&          dim,
	ConstVectorView       p_grid,
	ConstVectorView       lat_grid,
	ConstVectorView       lon_grid,
        const Index&          blackbody_ground,
        const Index&          cloudbox_on, 
	const ArrayOfIndex&   cloudbox_limits );


#endif  // checkinput_h
