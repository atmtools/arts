/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   interpolation.h
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri May  3 08:54:45 2002
  
  \brief  Header file for interpolation.cc.
*/

#ifndef interpolation_h
#define interpolation_h

#include "arts.h"
#include "matpackVII.h"

//! Structure to store a grid position.
/*! 
  A grid position specifies, where an interpolation point is, relative
  to the original grid. It consists of three parts, an Index giving the
  original grid index below the interpolation point, a Numeric
  giving the fractional distance to the next original grid point, and a
  Numeric giving 1 minus this number. Of course, the last element is
  redundant. However, it is efficient to store this, since it is used
  many times over.    

  For example, idx=3 and fd=.5 means that this interpolation point is
  half-way between index 3 and 4 of the original grid.

  Grid positions for a whole new grid are stored in an Array<GridPos>
  (called ArrayOfGridPos). 
*/
struct GridPos {
   Index   idx;			/*!< Original grid index below interpolation point. */
   Numeric fd;			/*!< Fractional distance to next point (0<=fd<=1). */
   Numeric fdr; 		/*!< 1-fd. */
};

//! An Array of grid positions.
/*! 
  See \ref GridPos for details.
*/

typedef Array<GridPos> ArrayOfGridPos;

// Function documentation is in .cc file.
void gridpos( ArrayOfGridPos& gp,
              ConstVectorView old_grid,
              ConstVectorView new_grid );

#endif // interpolation_h
