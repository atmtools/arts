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
  \file   interpolation.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri May  3 08:55:51 2002
  
  \brief  Interpolation routines.
  
  There are no general single-step interpolation functions in
  ARTS. Rather, there are a set of useful utility functions that
  can be used to achieve interpolation. Roughly, you can separate these
  into functions determining grid-position Arrays, functions determining
  interpolation weight Tensors, and functions applying the
  interpolation. 

  Doing an interpolation requires a chain of function calls:

  -# gridpos (one for each interpolation dimension)
  -# interpweights
  -# interp
  
*/

#include "interpolation.h"
#include "check_input.h"
#include "logic.h"

//! Set up a grid position Array. 
/*!
 This is the simplest function to set up a grid position Array.
 As usual, gp has to have the right dimension. There could be other
 helper functions to set up grid position Arrays, but right now I
 cannot think of any. 

 The old grid has to be strictly sorted. It can be in ascending or descending
 order. But there must not be any duplicate values.

 The new grid doesn't have to be sorted, but the function will be
 faster if it is sorted or mostly sorted.

 The beauty is, that this is all it needs to do also interpolation in
 higher dimensions: You just have to call gridpos for all the
 dimensions that you want to interpolate.

 Note also, that for this step you do not need the field itself at all!

 \param gp Output: Grid position Array.
 \param old_grid The original grid.
 \param new_grid The new grid where we want to have the interpolated values. 
*/
void gridpos( ArrayOfGridPos& gp,
              ConstVectorView old_grid,
              ConstVectorView new_grid )
{
  const Index n_old = old_grid.nelem();
  const Index n_new = new_grid.nelem();

  // Assert that gp has the right size:
  assert(is_size(gp,n_new));

  // Nice trick here: If old_grid is sorted in descending order, we
  // simply use a backwards view!
  ConstVectorView og ( ( old_grid[0] < old_grid[1] ) ?
		       old_grid :
		       old_grid[Range(n_old-1,n_old,-1)]
		       );

  //  cout << "og =\n" << og << "\n";
  
  // So, og should always be sorted in strictly ascending order.
  // (No duplicate values.)
  // Assert that this is so. This may depend on user input, however,
  // inside this elementary function is not the place to check for
  // that. There must be runtime checks on higher levels to insure
  // that all grids are sorted. The assertion here is just as a last
  // safety check.
  assert(is_increasing(og));

  // Get iterator to old grid and initailize to grid start:
  ConstIterator1D po     = old_grid.begin();
  // The end of the old grid:
  ConstIterator1D po_end = old_grid.end();

  // Loop over all points in the new grid:
  for ( Index i_new=0; i_new<n_new; ++i_new )
    {
      // Get a reference to the current element of gp:
      GridPos& tgp = gp[i_new];
      // And on the current value of the new grid:
      const Numeric& tng = new_grid[i_new];

      // 
      //      while ( tng < *po )	
    }
}
