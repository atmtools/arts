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


//! Output operator for GridPos.
/*!
  This is just intended for testing and debugging.
  
  \param os Output stream.
  \param gp Grid position.

  \return The output stream.
*/
std::ostream& operator<<(std::ostream& os, const GridPos& gp)
{
  os << gp.idx << " " << gp.fd[0] << " " << gp.fd[1] << "\n";
  return os;
}


//! Set up a grid position Array. 
/*!
 This is the simplest function to set up a grid position Array.
 As usual, gp has to have the right dimension. There could be other
 helper functions to set up grid position Arrays, but right now I
 cannot think of any. 

 The old grid has to be strictly sorted. It can be in ascending or descending
 order. But there must not be any duplicate values. Furthermore, the
 old grid must contain at least two points.

 The new grid doesn't have to be sorted, but the function will be
 faster if it is sorted or mostly sorted. It is ok if the new grid
 contains only 1 point.

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
  assert( is_size(gp,n_new) );

  // Assert, that the old grid has more than one element
  assert( 1 < n_old );

  // This function hast two parts, depending on whether old_grid is
  // sorted in ascending or descending order. Maybe that's not too
  // elegant, but it's the most efficient, because in this way there
  // is no additional runtime overhead to handle both cases.

  // We use only the first two elements to decide how the grid is
  // sorted. (The rest of the grid is checked later by an assert.)
  // If both values are the same, we still assume the grid is
  // ascending. However, this will lead to an assertion fail later on,
  // because the grid has to be strictly sorted.
  bool ascending = ( old_grid[0] <= old_grid[1] );

  if (ascending)  
    {
      // So, old_grid should always be sorted in strictly ascending order.
      // (No duplicate values.)
      // Assert that this is so. This may depend on user input, however,
      // inside this elementary function is not the place to check for
      // that. There must be runtime checks on higher levels to insure
      // that all grids are sorted. The assertion here is just as a last
      // safety check.
      assert( is_increasing(old_grid) );

      // We need this to make sure that the new grid points are inside the
      // old grid:
      const Numeric og_min = old_grid[0];
      const Numeric og_max = old_grid[n_old-1];

      // We will make no firm assumptions about the new grid. But the case
      // that we have in mind is that it is either also sorted, or at
      // least partially sorted, for example like this:
      // 5 4 3 2 2.5 3 4
      // This kind of sequence should be typical if we interpolate
      // atmospheric fields along a limb propagation path.

      // Let's get some idea where the first point in the new grid is,
      // relative to the old grid. We use linear interpolation between the
      // maximum and the minimum of the old grid for this. An assertion
      // checks that frac is between 0 and 1 (otherwise we would be
      // outside the old grid).
      Numeric frac = (new_grid[0]-og_min)/(og_max-og_min);
      assert( 0<=frac && 1>=frac );

      // Initialize current_position. This statement satisfies
      // 0 <= current_position <= n_old-2
      Index   current_position = (Index) rint(frac*(n_old-2));

      // The variables lower and upper are used to remember the value of
      // the old grid at current_position and one above current_position:
      Numeric lower = old_grid[current_position];
      Numeric upper = old_grid[current_position+1];

      // Loop over all points in the new grid:
      for ( Index i_new=0; i_new<n_new; ++i_new )
	{
	  // Get a reference to the current element of gp:
	  GridPos& tgp = gp[i_new];
	  // And on the current value of the new grid:
	  const Numeric tng = new_grid[i_new];
	  assert( og_min<=tng && og_max>=tng ); // New grid point must
						// be inside old grid.

	  //       cout << "lower: " << lower << "\n";
	  //       cout << "tng:   " << tng << "\n";
	  //       cout << "upper: " << upper << "\n";

	  // Is current_position too high?
	  if ( tng < lower )
	    {
	      do
		{
		  --current_position;
		  lower = old_grid[current_position];
		}
	      while ( tng < lower );

	      upper = old_grid[current_position+1];

	      tgp.idx = current_position;
	      tgp.fd[0] = (tng-lower)/(upper-lower);
	      tgp.fd[1] = 1.0 - tgp.fd[0];
	    }
	  else
	    {
	      // Is it too low? 
	      if ( tng > upper )
		{
		  do
		    {
		      ++current_position;
		      upper = old_grid[current_position+1];
		    }
		  while ( tng > upper );

		  lower = old_grid[current_position];

		  tgp.idx = current_position;
		  tgp.fd[0] = (tng-lower)/(upper-lower);
		  tgp.fd[1] = 1.0 - tgp.fd[0];
		}
	      else
		{
		  // None of the other two conditions were true. That means:
		  // lower <= tng <= upper. The current_position is still
		  // valid.
		  //
		  // Note that it is not uniquely determined, which
		  // current position we use if a new grid point is
		  // exactly on top of an old grid point. 
		  //
		  // As it is, we safe an extra treatment for the case
		  // that the coincident point is exactly at the upper boundary
		  // of the old grid. (In this case current_position must
		  // be the second to last point, otherwise interpolation
		  // will fail later!)
		  tgp.idx = current_position;
		  tgp.fd[0] = (tng-lower)/(upper-lower);
		  tgp.fd[1] = 1.0 - tgp.fd[0];
		}
	    }      
	}
    }
  else				//   if (ascending)  
    {
      // Now we are in the "descending old grid" part. We do exatly
      // the same as in the other part, just accounting for the
      // different order of things. Comments here refer only to
      // interesting differences from the ascending case. See that
      // case for more general comments.

      // This time ensure strictly descending order:
      assert( is_decreasing(old_grid) );

      // The max is now the first point, the min the last point!
      const Numeric og_max = old_grid[0];
      const Numeric og_min = old_grid[n_old-1];

      // We have to take 1- here, because we are starting from the
      // high end.
      Numeric frac = 1 - (new_grid[0]-og_min)/(og_max-og_min);
      assert( 0<=frac && 1>=frac );

      Index   current_position = (Index) rint(frac*(n_old-2));

      // Note, that old_grid[lower] has a higher numerical value than
      // old_grid[upper]! 
      Numeric lower = old_grid[current_position];
      Numeric upper = old_grid[current_position+1];

      for ( Index i_new=0; i_new<n_new; ++i_new )
	{
	  GridPos& tgp = gp[i_new];
	  const Numeric tng = new_grid[i_new];
	  assert( og_min<=tng && og_max>=tng );

	  //       cout << "lower: " << lower << "\n";
	  //       cout << "tng:   " << tng << "\n";
	  //       cout << "upper: " << upper << "\n";

	  // Is current_position too high? (Sign of comparison changed
	  // compared to ascending case!)
	  if ( tng > lower )
	    {
	      do
		{
		  --current_position;
		  lower = old_grid[current_position];
		}
	      while ( tng > lower );

	      upper = old_grid[current_position+1];

	      tgp.idx = current_position;
	      tgp.fd[0] = (tng-lower)/(upper-lower);
	      tgp.fd[1] = 1.0 - tgp.fd[0];
	    }
	  else
	    {
	      // Is it too low? (Sign of comparison changed
	  // compared to ascending case!)
	      if ( tng < upper )
		{
		  do
		    {
		      ++current_position;
		      upper = old_grid[current_position+1];
		    }
		  while ( tng < upper );

		  lower = old_grid[current_position];

		  tgp.idx = current_position;
		  tgp.fd[0] = (tng-lower)/(upper-lower);
		  tgp.fd[1] = 1.0 - tgp.fd[0];
		}
	      else
		{
		  // None of the other two conditions were true. That means:
		  // upper <= tng <= lower. The current_position is still
		  // valid. (Note that upper and lower have switched
		  // place compared to the ascending case.)

		  tgp.idx = current_position;
		  tgp.fd[0] = (tng-lower)/(upper-lower);
		  tgp.fd[1] = 1.0 - tgp.fd[0];
		}
	    }      
	}      
    }
}


//! Compute 1D interpolation weights.
/*! 
  For this 1D case there is no distinction between "blue" and "green"
  type interpolation.

  The dimensions of itw must be consistent with cgp.

  Note that we still do not need the actual field for this step.

  
  \param itw Output: Interpolation weights.
  \param cgp The grid position Array for the column dimension.
 */
void interpweights( MatrixView itw,
               	    const ArrayOfGridPos& cgp )
{
  Index n = cgp.nelem();
  assert(is_size(itw,n,2));	// We must store 2 interpolation
				// weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tc = cgp[i];

      // Current row of interpolation weight matrix:
      VectorView ti = itw(i,Range(joker));

      // Interpolation weights are stored in this order (l=lower
      // u=upper, c=column):
      // 1. l-c
      // 2. u-c

      Index iti = 0;
	for ( Index c=1; c>=0; --c )
	  {
	    ti[iti] = tc.fd[c];
	    ++iti;
	  }
    }
}

//! Compute 2D interpolation weights for a sequence of positions.
/*! 
 Compute the weights for a "blue" type interpolation of the field,
 that means that the grid position Arrays are interpreted as defining
 a sequence of positions. ALL GRID POSITION ARRAYS MUST HAVE THE SAME LENGTH! 

 The dimensions of itw must be also consistent with this.

 Note that we still do not need the actual field for this step.

 This function can be easily distinguished from the other
 interpweights function (for "green" interpolation), because the
 output is a Matrix, whereas in the other case it is a Tensor with one
 more dimension than there are input grid position Arrays.

 \param itw Output: Interpolation weights.
 \param rgp The grid position Array for the row dimension.
 \param cgp The grid position Array for the column dimension.
 */
void interpweights( MatrixView itw,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index n = cgp.nelem();
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,4));	// We must store 4 interpolation
				// weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      // Current row of interpolation weight matrix:
      VectorView ti = itw(i,Range(joker));

      // Interpolation weights are stored in this order (l=lower
      // u=upper, r=row, c=column):
      // 1. l-r l-c
      // 2. l-r u-c
      // 3. u-r l-c
      // 4. u-r u-c

      Index iti = 0;
      // FIXME: Is there a speed gain if I use int in these loops,
      // instead of Index?
      for ( Index r=1; r>=0; --r )
	for ( Index c=1; c>=0; --c )
	  {
	    ti[iti] = tr.fd[r] * tc.fd[c];
	    ++iti;
	  }
    }
}
