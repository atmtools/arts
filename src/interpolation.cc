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

#include <math.h>
#include "array.h"
#include "check_input.h"
#include "interpolation.h"
#include "logic.h"

//! Macro for interpolation weight loops.
/*!
  We use the macro LOOPIT to make the notation for the nested for
  loops in the interpweights functions more concise, and to avoid
  typing errors.

  Should resolve to something like:

  for ( const Numeric* p=&tp.fd[1]; p>=&tp.fd[0]; --p )
*/
#define LOOPIT(x) for ( const Numeric* x=&t##x.fd[1]; x>=&t##x.fd[0]; --x )


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



//! gridpos_check_fd
/*!
   Checks that the fractional distances have a value in the range [0,1].

   This function can be used when you are sure that the fractional distances
   have been calculated correctly, but the limited numerical precision can
   give values below 0 or above 1.

   \param   gp     Output: Grid position structure.

   \author Patrick Eriksson
   \date   2002-05-21
*/
void gridpos_check_fd( GridPos&   gp )
{
  if( gp.fd[0] < 0 )
    { gp.fd[0] = 0; }
  else if( gp.fd[0] > 1 )
    { gp.fd[0] = 1; }

  if( gp.fd[1] < 0 )
    { gp.fd[1] = 0; }
  else if( gp.fd[1] > 1 )
    { gp.fd[1] = 1; }
}



//! gridpos_force_end_fd
/*!
   Forces that the fractional distances are set to 0 or 1.

   This function can be called when it is known that a position is exactly
   on a grid point. The fractional distance of the grid position is then
   0 or 1, but rounding errors can give a slightly deviating value.

   The difference between this function and gridpos_check_fd is that this
   function is only applicable for end points, while the other function can
   be called for every point.

   The input fractional distances are not allowed to deviate freom 0 and 1
   with more than 1e-6.

   \param   gp     Output: Grid position structure.

   \author Patrick Eriksson
   \date   2002-05-22
*/
void gridpos_force_end_fd( GridPos&   gp )
{
  if( gp.fd[0] < 0.5 )
    {
      assert( fabs( gp.fd[0] ) <= 1e-6 );
      assert( fabs(gp.fd[1] -1 ) <= 1e-6 );
      gp.fd[0] = 0;
      gp.fd[1] = 1;
    }
  else
    {
      assert( fabs( gp.fd[1] ) <= 1e-6 );
      assert( fabs(gp.fd[0] -1 ) <= 1e-6 );
      gp.fd[0] = 1;
      gp.fd[1] = 0;
    }    
}



//! is_gridpos_at_index_i
/*!
   Determines if a grid position is at a given grid index.

   \return         True if at index i, else false.
   \param   gp     Grid position structure.
   \param   i      The grid index of interest.

   \author Patrick Eriksson
   \date   2002-05-22
*/
bool is_gridpos_at_index_i(  
       const GridPos&   gp,
       const Index&     i )
{
  if( ( gp.fd[0] == 0  ||  gp.fd[0] == 1 )  &&  
                                            ( gp.idx + Index(gp.fd[0]) ) == i )
    { return true; }
  else
    { return false; }
}



//! gridpos2gridrange
/*!
   Determines which grid range that is of interest for a given grid position.

   The purpose of the function is to determine which two grid values that
   sorround the given point. The index of the lower grid value is returned.

   For a point exactly on a grid value it is not clear if it is the range 
   below or above that is of interest. The input argument upward is used to 
   resolve such cases, where upward != 0 means that it is the range above
   that is of interest.

   \return         The index of the lower end of the grid range.
   \param   gp     Grid position structure.
   \param   upward Direction of interest, see above.

   \author Patrick Eriksson
   \date   2002-05-20
*/
Index gridpos2gridrange(
       const GridPos&   gp,
       const bool&      upwards )
{
  assert( gp.fd[0] >= 0 );
  assert( gp.fd[0] <= 1 );

  // Not at a grid point
  if( gp.fd[0] > 0   &&  gp.fd[0] < 1 )
    {
      return gp.idx;
    }

  // Fractional distance 0
  else if( gp.fd[0] == 0 )
    {
      if( upwards )
	return gp.idx;
      else
	return gp.idx - 1;
    }

  // Fractional distance 1
  else
    {
      if( upwards )
	return gp.idx + 1;
      else
	return gp.idx;
    }
}



////////////////////////////////////////////////////////////////////////////
//			Blue interpolation
////////////////////////////////////////////////////////////////////////////

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

      // Interpolation weights are stored in this order (l=lower
      // u=upper, c=column):
      // 1. l-c
      // 2. u-c

      Index iti = 0;

      // We could use a straight-forward for loop here:
      //
      //       for ( Index c=1; c>=0; --c )
      // 	{
      // 	  ti[iti] = tc.fd[c];
      // 	  ++iti;
      // 	}
      //
      // However, it is slightly faster to use pointers (I tried it,
      // the speed gain is about 20%). So we should write something
      // like:
      //
      //       for ( const Numeric* c=&tc.fd[1]; c>=&tc.fd[0]; --c )
      // 	{
      // 	  ti[iti] = *c;
      // 	  ++iti;
      // 	}
      //
      // For higher dimensions we have to nest these loops. To avoid
      // typos and safe typing, I use the LOOPIT macro, which expands
      // to the for loop above. Note: NO SEMICOLON AFTER THE LOOPIT
      // COMMAND! 

      LOOPIT(c)
	{
	  itw(i,iti) = *c;
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

      // Interpolation weights are stored in this order (l=lower
      // u=upper, r=row, c=column):
      // 1. l-r l-c
      // 2. l-r u-c
      // 3. u-r l-c
      // 4. u-r u-c

      Index iti = 0;

      LOOPIT(r)
      LOOPIT(c)
	  {
	    itw(i,iti) = (*r) * (*c);
	    ++iti;
	  }
    }
}

//! Compute 3D interpolation weights for a sequence of positions.
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
 \param pgp The grid position Array for the page dimension.
 \param rgp The grid position Array for the row dimension.
 \param cgp The grid position Array for the column dimension.
 */
void interpweights( MatrixView itw,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index n = cgp.nelem();
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,8));	// We must store 8 interpolation
				// weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      Index iti = 0;
      LOOPIT(p)
      LOOPIT(r)
      LOOPIT(c)
	{
	  itw(i,iti) = (*p) * (*r) * (*c);
	  ++iti;
	}
    }
}

//! Compute 4D interpolation weights for a sequence of positions.
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
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 */
void interpweights( MatrixView itw,
               	    const ArrayOfGridPos& bgp,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index n = cgp.nelem();
  assert(is_size(bgp,n));	// bgp must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,16));	// We must store 16 interpolation
				// weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tb = bgp[i];
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      Index iti = 0;
      LOOPIT(b)
      LOOPIT(p)
      LOOPIT(r)
      LOOPIT(c)
	{
	  itw(i,iti) = (*b) * (*p) * (*r) * (*c);
	  ++iti;
	}
    }
}

//! Compute 5D interpolation weights for a sequence of positions.
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
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 */
void interpweights( MatrixView itw,
               	    const ArrayOfGridPos& sgp,
               	    const ArrayOfGridPos& bgp,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index n = cgp.nelem();
  assert(is_size(sgp,n));	// sgp must have same size as cgp.
  assert(is_size(bgp,n));	// bgp must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,32));	// We must store 32 interpolation
				// weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& ts = sgp[i];
      const GridPos& tb = bgp[i];
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      Index iti = 0;
      LOOPIT(s)
      LOOPIT(b)
      LOOPIT(p)
      LOOPIT(r)
      LOOPIT(c)
	{
	  itw(i,iti) = (*s) * (*b) * (*p) * (*r) * (*c);
	  ++iti;
	}
    }
}

//! Compute 6D interpolation weights for a sequence of positions.
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
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 */
void interpweights( MatrixView itw,
               	    const ArrayOfGridPos& vgp,
               	    const ArrayOfGridPos& sgp,
               	    const ArrayOfGridPos& bgp,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index n = cgp.nelem();
  assert(is_size(vgp,n));	// vgp must have same size as cgp.
  assert(is_size(sgp,n));	// sgp must have same size as cgp.
  assert(is_size(bgp,n));	// bgp must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,64));	// We must store 64 interpolation
				// weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tv = vgp[i];
      const GridPos& ts = sgp[i];
      const GridPos& tb = bgp[i];
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      Index iti = 0;
      LOOPIT(v)
      LOOPIT(s)
      LOOPIT(b)
      LOOPIT(p)
      LOOPIT(r)
      LOOPIT(c)
	{
	  itw(i,iti) = (*v) * (*s) * (*b) * (*p) * (*r) * (*c);
	  ++iti;
	}
    }
}

//! Interpolate 1D field.
/*! 
  For this 1D case there is no distinction between "blue" and "green"
  type interpolation.

  The output vector ia must have the same length as the grid position
  vector cgp. And the dimension of itw must be consistent with
  this.

  \param ia  Output: Vector containing the interpolated field values.
  \param itw Interpolation weights.
  \param a   The field to interpolate.
  \param cgp The grid position Array for the column dimension.
*/
void interp( VectorView      	   ia,
             ConstMatrixView 	   itw,
             ConstVectorView 	   a,    
             const ArrayOfGridPos& cgp)
{
  Index n = cgp.nelem();
  assert(is_size(ia,n));	//  ia must have same size as cgp.
  assert(is_size(itw,n,2));	// We need 2 interpolation
				// weights for each position.
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index c=0; c<2; ++c )
	{
	  tia += a[tc.idx+c] * itw(i,iti);
	  ++iti;
	}
    }
}

//! Interpolate 2D field to a sequence of positions.
/*! 
 This performs a "blue" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining a sequence
 of positions. ALL GRID POSITION ARRAYS MUST HAVE THE SAME LENGTH! 

 The output vector ia also must have the same length. And the
 dimension of itw must be consistent with this.

 This function can be easily distinguished from the other
 interpolation function (that creates an entire field of interpolated
 values), because of the dimension of ia and itw.

 \param ia  Output: Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param rgp The grid position Array for the row    dimension.
 \param cgp The grid position Array for the column dimension.
*/
void interp( VectorView      	   ia,
             ConstMatrixView 	   itw,
             ConstMatrixView 	   a,    
             const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index n = cgp.nelem();
  assert(is_size(ia,n));	//  ia must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,4));	// We need 4 interpolation
				// weights for each position.
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index r=0; r<2; ++r )
	for ( Index c=0; c<2; ++c )
	  {
	    tia += a(tr.idx+r,
		     tc.idx+c) * itw(i,iti);
	    ++iti;
	  }
    }
}

//! Interpolate 3D field to a sequence of positions.
/*! 
 This performs a "blue" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining a sequence
 of positions. ALL GRID POSITION ARRAYS MUST HAVE THE SAME LENGTH! 

 The output vector ia also must have the same length. And the
 dimension of itw must be consistent with this.

 This function can be easily distinguished from the other
 interpolation function (that creates an entire field of interpolated
 values), because of the dimension of ia and itw.

 \param ia  Output: Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( VectorView      	   ia,
             ConstMatrixView 	   itw,
             ConstTensor3View 	   a,    
	     const ArrayOfGridPos& pgp,
             const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index n = cgp.nelem();
  assert(is_size(ia,n));	//  ia must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,8));	// We need 8 interpolation
				// weights for each position.
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index p=0; p<2; ++p )
	for ( Index r=0; r<2; ++r )
	  for ( Index c=0; c<2; ++c )
	    {
	      tia += a(tp.idx+p,
		       tr.idx+r,
		       tc.idx+c) * itw(i,iti);
	      ++iti;
	    }
    }
}

//! Interpolate 4D field to a sequence of positions.
/*! 
 This performs a "blue" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining a sequence
 of positions. ALL GRID POSITION ARRAYS MUST HAVE THE SAME LENGTH! 

 The output vector ia also must have the same length. And the
 dimension of itw must be consistent with this.

 This function can be easily distinguished from the other
 interpolation function (that creates an entire field of interpolated
 values), because of the dimension of ia and itw.

 \param ia  Output: Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( VectorView      	   ia,
             ConstMatrixView 	   itw,
             ConstTensor4View 	   a,    
	     const ArrayOfGridPos& bgp,
	     const ArrayOfGridPos& pgp,
             const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index n = cgp.nelem();
  assert(is_size(ia,n));	//  ia must have same size as cgp.
  assert(is_size(bgp,n));	// bgp must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,16));	// We need 16 interpolation
				// weights for each position.
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tb = bgp[i];
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index b=0; b<2; ++b )
	for ( Index p=0; p<2; ++p )
	  for ( Index r=0; r<2; ++r )
	    for ( Index c=0; c<2; ++c )
	      {
		tia += a(tb.idx+b,
			 tp.idx+p,
			 tr.idx+r,
			 tc.idx+c) * itw(i,iti);
		++iti;
	      }
    }
}

//! Interpolate 5D field to a sequence of positions.
/*! 
 This performs a "blue" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining a sequence
 of positions. ALL GRID POSITION ARRAYS MUST HAVE THE SAME LENGTH! 

 The output vector ia also must have the same length. And the
 dimension of itw must be consistent with this.

 This function can be easily distinguished from the other
 interpolation function (that creates an entire field of interpolated
 values), because of the dimension of ia and itw.

 \param ia  Output: Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( VectorView      	   ia,
             ConstMatrixView 	   itw,
             ConstTensor5View 	   a,    
	     const ArrayOfGridPos& sgp,
	     const ArrayOfGridPos& bgp,
	     const ArrayOfGridPos& pgp,
             const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index n = cgp.nelem();
  assert(is_size(ia,n));	//  ia must have same size as cgp.
  assert(is_size(sgp,n));	// sgp must have same size as cgp.
  assert(is_size(bgp,n));	// bgp must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,32));	// We need 32 interpolation
				// weights for each position.
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& ts = sgp[i];
      const GridPos& tb = bgp[i];
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index s=0; s<2; ++s )
	for ( Index b=0; b<2; ++b )
	  for ( Index p=0; p<2; ++p )
	    for ( Index r=0; r<2; ++r )
	      for ( Index c=0; c<2; ++c )
		{
		  tia += a(ts.idx+s,
			   tb.idx+b,
			   tp.idx+p,
			   tr.idx+r,
			   tc.idx+c) * itw(i,iti);
		  ++iti;
		}
    }
}

//! Interpolate 6D field to a sequence of positions.
/*! 
 This performs a "blue" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining a sequence
 of positions. ALL GRID POSITION ARRAYS MUST HAVE THE SAME LENGTH! 

 The output vector ia also must have the same length. And the
 dimension of itw must be consistent with this.

 This function can be easily distinguished from the other
 interpolation function (that creates an entire field of interpolated
 values), because of the dimension of ia and itw.

 \param ia  Output: Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( VectorView      	   ia,
             ConstMatrixView 	   itw,
             ConstTensor6View 	   a,    
	     const ArrayOfGridPos& vgp,
	     const ArrayOfGridPos& sgp,
	     const ArrayOfGridPos& bgp,
	     const ArrayOfGridPos& pgp,
             const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index n = cgp.nelem();
  assert(is_size(ia,n));	//  ia must have same size as cgp.
  assert(is_size(vgp,n));	// vgp must have same size as cgp.
  assert(is_size(sgp,n));	// sgp must have same size as cgp.
  assert(is_size(bgp,n));	// bgp must have same size as cgp.
  assert(is_size(pgp,n));	// pgp must have same size as cgp.
  assert(is_size(rgp,n));	// rgp must have same size as cgp.
  assert(is_size(itw,n,64));	// We need 64 interpolation
				// weights for each position.
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPos& tv = vgp[i];
      const GridPos& ts = sgp[i];
      const GridPos& tb = bgp[i];
      const GridPos& tp = pgp[i];
      const GridPos& tr = rgp[i];
      const GridPos& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index v=0; v<2; ++v )
	for ( Index s=0; s<2; ++s )
	  for ( Index b=0; b<2; ++b )
	    for ( Index p=0; p<2; ++p )
	      for ( Index r=0; r<2; ++r )
		for ( Index c=0; c<2; ++c )
		  {
		    tia += a(tv.idx+v,
			     ts.idx+s,
			     tb.idx+b,
			     tp.idx+p,
			     tr.idx+r,
			     tc.idx+c) * itw(i,iti);
		    ++iti;
		  }
    }
}

////////////////////////////////////////////////////////////////////////////
//			Green interpolation
////////////////////////////////////////////////////////////////////////////

//! Compute 2D interpolation weights for an entire field.
/*! 
 Compute the weights for a "green" type interpolation of the field,
 that means that the grid position Arrays are interpreted as defining
 the grids for the interpolated field.

 The dimensions of itw must be consistent with this.

 Note that we still do not need the actual field for this step.

 This function can be easily distinguished from the other
 interpweights function (for "green" interpolation), because the
 output is a Tensor with one more dimension than the number of grid
 position Arrays.

 \param itw Output: Interpolation weights
 \param rgp The grid position Array for the row dimension.
 \param cgp The grid position Array for the column dimension.
 
 */
void interpweights( Tensor3View itw,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  assert(is_size(itw,nr,nc,4));	// We must store 4 interpolation
				// weights for each position.

  // We have to loop all the points in the new grid:
  for ( Index ir=0; ir<nr; ++ir )
    {
      // Current grid position:
      const GridPos& tr = rgp[ir];

      for ( Index ic=0; ic<nc; ++ic )
	{
	  // Current grid position:
	  const GridPos& tc = cgp[ic];

	  // Interpolation weights are stored in this order (l=lower
	  // u=upper, r=row, c=column):
	  // 1. l-r l-c
	  // 2. l-r u-c
	  // 3. u-r l-c
	  // 4. u-r u-c

	  Index iti = 0;

	  LOOPIT(r)
	    LOOPIT(c)
	    {
	      itw(ir,ic,iti) = (*r) * (*c);
	      ++iti;
	    }
	}
    }
}

//! Compute 3D interpolation weights for an entire field.
/*! 
 Compute the weights for a "green" type interpolation of the field,
 that means that the grid position Arrays are interpreted as defining
 the grids for the interpolated field.

 The dimensions of itw must be consistent with this.

 Note that we still do not need the actual field for this step.

 This function can be easily distinguished from the other
 interpweights function (for "green" interpolation), because the
 output is a Tensor with one more dimension than the number of grid
 position Arrays.

 \param itw Output: Interpolation weights
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 */
void interpweights( Tensor4View itw,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 8 interpolation weights for each position:
  assert(is_size(itw,np,nr,nc,8));	

  // We have to loop all the points in the new grid:
  for ( Index ip=0; ip<np; ++ip )
    {
      const GridPos& tp = pgp[ip];
      for ( Index ir=0; ir<nr; ++ir )
	{
	  const GridPos& tr = rgp[ir];
	  for ( Index ic=0; ic<nc; ++ic )
	    {
	      const GridPos& tc = cgp[ic];

	      Index iti = 0;

	      LOOPIT(p)
		LOOPIT(r)
		LOOPIT(c)
		{
		  itw(ip,ir,ic,iti) =
		    (*p) * (*r) * (*c);
		  ++iti;
		}
	    }
	}
    }
}

//! Compute 4D interpolation weights for an entire field.
/*! 
 Compute the weights for a "green" type interpolation of the field,
 that means that the grid position Arrays are interpreted as defining
 the grids for the interpolated field.

 The dimensions of itw must be consistent with this.

 Note that we still do not need the actual field for this step.

 This function can be easily distinguished from the other
 interpweights function (for "green" interpolation), because the
 output is a Tensor with one more dimension than the number of grid
 position Arrays.

 \param itw Output: Interpolation weights
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 */
void interpweights( Tensor5View itw,
               	    const ArrayOfGridPos& bgp,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 16 interpolation weights for each position:
  assert(is_size(itw,nb,np,nr,nc,16));	

  // We have to loop all the points in the new grid:
  for ( Index ib=0; ib<nb; ++ib )
    {
      const GridPos& tb = bgp[ib];
      for ( Index ip=0; ip<np; ++ip )
	{
	  const GridPos& tp = pgp[ip];
	  for ( Index ir=0; ir<nr; ++ir )
	    {
	      const GridPos& tr = rgp[ir];
	      for ( Index ic=0; ic<nc; ++ic )
		{
		  const GridPos& tc = cgp[ic];

		  Index iti = 0;

		  LOOPIT(b)
		    LOOPIT(p)
		    LOOPIT(r)
		    LOOPIT(c)
		    {
		      itw(ib,ip,ir,ic,iti) =
			(*b) * (*p) * (*r) * (*c);
		      ++iti;
		    }
		}
	    }
	}
    }
}

//! Compute 5D interpolation weights for an entire field.
/*! 
 Compute the weights for a "green" type interpolation of the field,
 that means that the grid position Arrays are interpreted as defining
 the grids for the interpolated field.

 The dimensions of itw must be consistent with this.

 Note that we still do not need the actual field for this step.

 This function can be easily distinguished from the other
 interpweights function (for "green" interpolation), because the
 output is a Tensor with one more dimension than the number of grid
 position Arrays.

 \param itw Output: Interpolation weights
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 */
void interpweights( Tensor6View itw,
               	    const ArrayOfGridPos& sgp,
               	    const ArrayOfGridPos& bgp,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 32 interpolation weights for each position:
  assert(is_size(itw,ns,nb,np,nr,nc,32));	

  // We have to loop all the points in the new grid:
  for ( Index is=0; is<ns; ++is )
    {
      const GridPos& ts = sgp[is];
      for ( Index ib=0; ib<nb; ++ib )
	{
	  const GridPos& tb = bgp[ib];
	  for ( Index ip=0; ip<np; ++ip )
	    {
	      const GridPos& tp = pgp[ip];
	      for ( Index ir=0; ir<nr; ++ir )
		{
		  const GridPos& tr = rgp[ir];
		  for ( Index ic=0; ic<nc; ++ic )
		    {
		      const GridPos& tc = cgp[ic];

		      Index iti = 0;

		      LOOPIT(s)
			LOOPIT(b)
			LOOPIT(p)
			LOOPIT(r)
			LOOPIT(c)
			{
			  itw(is,ib,ip,ir,ic,iti) =
			    (*s) * (*b) * (*p) * (*r) * (*c);
			  ++iti;
			}
		    }
		}
	    }
	}
    }
}

//! Compute 6D interpolation weights for an entire field.
/*! 
 Compute the weights for a "green" type interpolation of the field,
 that means that the grid position Arrays are interpreted as defining
 the grids for the interpolated field.

 The dimensions of itw must be consistent with this.

 Note that we still do not need the actual field for this step.

 This function can be easily distinguished from the other
 interpweights function (for "green" interpolation), because the
 output is a Tensor with one more dimension than the number of grid
 position Arrays.

 \param itw Output: Interpolation weights
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 */
void interpweights( Tensor7View itw,
               	    const ArrayOfGridPos& vgp,
               	    const ArrayOfGridPos& sgp,
               	    const ArrayOfGridPos& bgp,
               	    const ArrayOfGridPos& pgp,
               	    const ArrayOfGridPos& rgp,
               	    const ArrayOfGridPos& cgp )
{
  Index nv = vgp.nelem();
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 64 interpolation weights for each position:
  assert(is_size(itw,nv,ns,nb,np,nr,nc,64));	

  // We have to loop all the points in the new grid:
  for ( Index iv=0; iv<nv; ++iv )
    {
      const GridPos& tv = vgp[iv];
      for ( Index is=0; is<ns; ++is )
	{
	  const GridPos& ts = sgp[is];
	  for ( Index ib=0; ib<nb; ++ib )
	    {
	      const GridPos& tb = bgp[ib];
	      for ( Index ip=0; ip<np; ++ip )
		{
		  const GridPos& tp = pgp[ip];
		  for ( Index ir=0; ir<nr; ++ir )
		    {
		      const GridPos& tr = rgp[ir];
		      for ( Index ic=0; ic<nc; ++ic )
			{
			  const GridPos& tc = cgp[ic];

			  Index iti = 0;

			  LOOPIT(v)
			    LOOPIT(s)
			    LOOPIT(b)
			    LOOPIT(p)
			    LOOPIT(r)
			    LOOPIT(c)
			    {
			      itw(iv,is,ib,ip,ir,ic,iti) =
				(*v) * (*s) * (*b) * (*p) * (*r) * (*c);
			      ++iti;
			    }
			}
		    }
		}
	    }
	}
    }
}

//! Interpolate 2D field to another 2D field.
/*! 
 This performs a "green" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining the grids
 for the interpolated field.

 This function can be easily distinguished from the other
 interpolation function (that creates a sequence of interpolated
 values), because of the dimension of ia and itw.

 The size of ia and itw in all dimensions must be consistent with the grid
 position Arrays.

 \param ia  Output: Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param rgp The grid position Array for the row dimension.
 \param cgp The grid position Array for the column dimension.
 
 */
void interp( MatrixView       	   ia,
             ConstTensor3View 	   itw,
             ConstMatrixView  	   a,   
	     const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  assert(is_size(ia,nr,nc));    
  assert(is_size(itw,nr,nc,4));	// We need 4 interpolation
				// weights for each position.

  // We have to loop all the points in the new grid:
  for ( Index ir=0; ir<nr; ++ir )
    {
      // Current grid position:
      const GridPos& tr = rgp[ir];

      for ( Index ic=0; ic<nc; ++ic )
	{
	  // Current grid position:
	  const GridPos& tc = cgp[ic];

	  // Get handle to current element of output tensor and initialize
	  // it to zero:
	  Numeric& tia = ia(ir,ic);
	  tia = 0;

	  Index iti = 0;
	  for ( Index r=0; r<2; ++r )
	    for ( Index c=0; c<2; ++c )
	    {
	      tia += a(tr.idx+r,
		       tc.idx+c) * itw(ir,ic,iti);
	      ++iti;
	    }
	}
    }
}

//! Interpolate 3D field to another 3D field.
/*! 
 This performs a "green" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining the grids
 for the interpolated field.

 This function can be easily distinguished from the other
 interpolation function (that creates a sequence of interpolated
 values), because of the dimension of ia and itw.

 The size of ia and itw in all dimensions must be consistent with the grid
 position Arrays.

 \param ia  Output: Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( Tensor3View       	   ia,
             ConstTensor4View 	   itw,
             ConstTensor3View  	   a,   
	     const ArrayOfGridPos& pgp,
	     const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  assert(is_size(ia,
		 np,nr,nc));    
  assert(is_size(itw,
		 np,nr,nc,
		 8));

  // We have to loop all the points in the new grid:
  for ( Index ip=0; ip<np; ++ip )
    {
      const GridPos& tp = pgp[ip];
      for ( Index ir=0; ir<nr; ++ir )
	{
	  const GridPos& tr = rgp[ir];
	  for ( Index ic=0; ic<nc; ++ic )
	    {
	      // Current grid position:
	      const GridPos& tc = cgp[ic];

	      // Get handle to current element of output tensor and
	      // initialize it to zero:
	      Numeric& tia = ia(ip,ir,ic);
	      tia = 0;

	      Index iti = 0;
	      for ( Index p=0; p<2; ++p )
		for ( Index r=0; r<2; ++r )
		  for ( Index c=0; c<2; ++c )
		    {
		      tia += a(tp.idx+p,
			       tr.idx+r,
			       tc.idx+c) * itw(ip,ir,ic,
					       iti);
		      ++iti;
		    }
	    }
	}
    }
}

//! Interpolate 4D field to another 4D field.
/*! 
 This performs a "green" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining the grids
 for the interpolated field.

 This function can be easily distinguished from the other
 interpolation function (that creates a sequence of interpolated
 values), because of the dimension of ia and itw.

 The size of ia and itw in all dimensions must be consistent with the grid
 position Arrays.

 \param ia  Output: Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( Tensor4View       	   ia,
             ConstTensor5View 	   itw,
             ConstTensor4View  	   a,   
	     const ArrayOfGridPos& bgp,
	     const ArrayOfGridPos& pgp,
	     const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  assert(is_size(ia,
		 nb,np,nr,nc));    
  assert(is_size(itw,
		 nb,np,nr,nc,
		 16));

  // We have to loop all the points in the new grid:
  for ( Index ib=0; ib<nb; ++ib )
    {
      const GridPos& tb = bgp[ib];
      for ( Index ip=0; ip<np; ++ip )
	{
	  const GridPos& tp = pgp[ip];
	  for ( Index ir=0; ir<nr; ++ir )
	    {
	      const GridPos& tr = rgp[ir];
	      for ( Index ic=0; ic<nc; ++ic )
		{
		  // Current grid position:
		  const GridPos& tc = cgp[ic];

		  // Get handle to current element of output tensor and
		  // initialize it to zero:
		  Numeric& tia = ia(ib,ip,ir,ic);
		  tia = 0;

		  Index iti = 0;
		  for ( Index b=0; b<2; ++b )
		    for ( Index p=0; p<2; ++p )
		      for ( Index r=0; r<2; ++r )
			for ( Index c=0; c<2; ++c )
			  {
			    tia += a(tb.idx+b,
				     tp.idx+p,
				     tr.idx+r,
				     tc.idx+c) * itw(ib,ip,ir,ic,
						     iti);
			    ++iti;
			  }
		}
	    }
	}
    }
}

//! Interpolate 5D field to another 5D field.
/*! 
 This performs a "green" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining the grids
 for the interpolated field.

 This function can be easily distinguished from the other
 interpolation function (that creates a sequence of interpolated
 values), because of the dimension of ia and itw.

 The size of ia and itw in all dimensions must be consistent with the grid
 position Arrays.

 \param ia  Output: Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( Tensor5View       	   ia,
             ConstTensor6View 	   itw,
             ConstTensor5View  	   a,   
	     const ArrayOfGridPos& sgp,
	     const ArrayOfGridPos& bgp,
	     const ArrayOfGridPos& pgp,
	     const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  assert(is_size(ia,
		 ns,nb,np,nr,nc));    
  assert(is_size(itw,
		 ns,nb,np,nr,nc,
		 32));

  // We have to loop all the points in the new grid:
  for ( Index is=0; is<ns; ++is )
    {
      const GridPos& ts = sgp[is];
      for ( Index ib=0; ib<nb; ++ib )
	{
	  const GridPos& tb = bgp[ib];
	  for ( Index ip=0; ip<np; ++ip )
	    {
	      const GridPos& tp = pgp[ip];
	      for ( Index ir=0; ir<nr; ++ir )
		{
		  const GridPos& tr = rgp[ir];
		  for ( Index ic=0; ic<nc; ++ic )
		    {
		      // Current grid position:
		      const GridPos& tc = cgp[ic];

		      // Get handle to current element of output tensor and
		      // initialize it to zero:
		      Numeric& tia = ia(is,ib,ip,ir,ic);
		      tia = 0;

		      Index iti = 0;
		      for ( Index s=0; s<2; ++s )
			for ( Index b=0; b<2; ++b )
			  for ( Index p=0; p<2; ++p )
			    for ( Index r=0; r<2; ++r )
			      for ( Index c=0; c<2; ++c )
				{
				  tia += a(ts.idx+s,
					   tb.idx+b,
					   tp.idx+p,
					   tr.idx+r,
					   tc.idx+c) * itw(is,ib,ip,ir,ic,
							   iti);
				  ++iti;
				}
		    }
		}
	    }
	}
    }
}

//! Interpolate 6D field to another 6D field.
/*! 
 This performs a "green" type interpolation of the field, that means
 that the grid position Arrays are interpreted as defining the grids
 for the interpolated field.

 This function can be easily distinguished from the other
 interpolation function (that creates a sequence of interpolated
 values), because of the dimension of ia and itw.

 The size of ia and itw in all dimensions must be consistent with the grid
 position Arrays.

 \param ia  Output: Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
*/
void interp( Tensor6View       	   ia,
             ConstTensor7View 	   itw,
             ConstTensor6View  	   a,   
	     const ArrayOfGridPos& vgp,
	     const ArrayOfGridPos& sgp,
	     const ArrayOfGridPos& bgp,
	     const ArrayOfGridPos& pgp,
	     const ArrayOfGridPos& rgp,
             const ArrayOfGridPos& cgp)
{
  Index nv = vgp.nelem();
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  assert(is_size(ia,
		 nv,ns,nb,np,nr,nc));    
  assert(is_size(itw,
		 nv,ns,nb,np,nr,nc,
		 64));

  // We have to loop all the points in the new grid:
  for ( Index iv=0; iv<nv; ++iv )
    {
      const GridPos& tv = vgp[iv];
      for ( Index is=0; is<ns; ++is )
	{
	  const GridPos& ts = sgp[is];
	  for ( Index ib=0; ib<nb; ++ib )
	    {
	      const GridPos& tb = bgp[ib];
	      for ( Index ip=0; ip<np; ++ip )
		{
		  const GridPos& tp = pgp[ip];
		  for ( Index ir=0; ir<nr; ++ir )
		    {
		      const GridPos& tr = rgp[ir];
		      for ( Index ic=0; ic<nc; ++ic )
			{
			  // Current grid position:
			  const GridPos& tc = cgp[ic];

			  // Get handle to current element of output tensor and
			  // initialize it to zero:
			  Numeric& tia = ia(iv,is,ib,ip,ir,ic);
			  tia = 0;

			  Index iti = 0;
			  for ( Index v=0; v<2; ++v )
			    for ( Index s=0; s<2; ++s )
			      for ( Index b=0; b<2; ++b )
				for ( Index p=0; p<2; ++p )
				  for ( Index r=0; r<2; ++r )
				    for ( Index c=0; c<2; ++c )
				      {
					tia += a(tv.idx+v,
						 ts.idx+s,
						 tb.idx+b,
						 tp.idx+p,
						 tr.idx+r,
						 tc.idx+c) * itw(iv,is,ib,ip,ir,ic,
								 iti);
					++iti;
				      }
			}
		    }
		}
	    }
	}
    }
}
