/*!
  \file   interpolation.cc
  \author Stefan Buehler <sbuehler@ltu.se>
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
#include <cmath>
#include <iostream>
#include "array.h"
#include "check_input.h"
#include "logic.h"

// File-global constants:

//! The maximum difference from 1 that we allow for a sum check.
/*!
  The sum check makes sure that the sum of all weights is
  approximately 1.

  We cannot use a sharp comparison there, due to numerical
  noise. The value of 1e-6 is an ad-hoc value.

  This shold be ok, the main point of the test is to make sure that
  what we have really *are* interpolation weights, and not something
  else. 
*/
DEBUG_ONLY(const Numeric sum_check_epsilon = 1e-6;)

//! Allowed tolerance for fractional distance values
/*!
    Fractional distances are not allowed to be negative or exceed
    1 above the set tolerance. Values deviating less than the
    tolerance from 0 and 1 will be treated as end points, if such a
    check is performed.     

    We cannot use a sharp comparison there, due to numerical
    noise. The value of 1e-3 is an ad-hoc value.
*/
const Numeric FD_TOL = 1.5e-3;

//! Macro for interpolation weight loops.
/*!
  We use the macro LOOPIT to make the notation for the nested for
  loops in the interpweights functions more concise, and to avoid
  typing errors.

  Should resolve to something like:

  for ( const Numeric* p=&tp.fd[1]; p>=&tp.fd[0]; --p )
*/
#define LOOPIT(x) for (const Numeric* x = &t##x.fd[1]; x >= &t##x.fd[0]; --x)

//! Output operator for GridPos.
/*!
  This is just intended for testing and debugging.
  
  \param os Output stream.
  \param gp Grid position.

  \return The output stream.
*/
ostream& operator<<(ostream& os, const GridPos& gp) {
  os << gp.idx << " " << gp.fd[0] << " " << gp.fd[1] << "\n";
  return os;
}

//! Set up a grid position Array.
/*! 
 This is the function to find the position in the original grid
 and set the fd for the linear interpolation case.

 The old grid has to be strictly sorted. It can be in ascending or
 descending order. But there must not be any duplicate
 values. Furthermore, the old grid must contain at least two points.

 The new grid doesn't have to be sorted, but the function will be
 faster if it is sorted or mostly sorted. It is ok if the new grid
 contains only 1 point.

 The beauty is, that this is all it needs to do also interpolation in
 higher dimensions: You just have to call gridpos for all the
 dimensions that you want to interpolate.

 Note also, that for this step you do not need the field itself at all!

 The new grid has basically to be inside the range covered by old
 grid, but some margins are given for practical reasons. For higher
 flexibility, a selectable extrapolation is allowed. The allowed
 extrapolation is given in fractions of the distance between the
 lowermost two points for extrapolation beyond the lower grid limit and the
 distance between the uppermost two points for extrapolation beyond the upper
 grid limit, respectively. For example, if *extpolfac* is set to 0.5 an
 extrapolation of half the distance between two outermost points at the
 respective edge of the grid is allowed. In this case, if the old grid is
 given as [0, 1, ... , 3, 3.4] the new grid can be extended down to -0.5 and
 up to 3.6.

 0.5 is the default value for extpolfac. Normally, you should just use
 the function with 3 arguments, which means that the default value
 will be used.
 
 If a new grid point is exactly on top of an original grid point, then 
 the returned gp.idx is guaranteed to point exactly at this point, so 
 that gp.fd[0] = 0. The only exception from this rule is if the point 
 happens to be the last point in the original grid. In that case gp.idx 
 will point to the original grid point below. 
 
 Furthermore, you are guaranteed that for the ascending old grid case:
      old_grid[tgp.idx]<=tng || tgp.idx==0
 And for the descending old grid case: 
      old_grid[tgp.idx]>=tng || tgp.idx==0
 
 \retval gp        Grid position Array.
 \param old_grid   The original grid.
 \param new_grid   The new grid where we want to have the interpolated values. 
 \param extpolfac  Extrapolation factor. Has a default value of
                   0.5. You should normally not specify this parameter!

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void gridpos(ArrayOfGridPos& gp,
             ConstVectorView old_grid,
             ConstVectorView new_grid,
             const Numeric& extpolfac) {
  const Index n_old = old_grid.nelem();
  const Index n_new = new_grid.nelem();

  // Assert that gp has the right size:
  ARTS_ASSERT(is_size(gp, n_new));

  // Assert, that the old grid has more than one element
  ARTS_ASSERT(1 < n_old);

  // This function hast two parts, depending on whether old_grid is
  // sorted in ascending or descending order. Maybe that's not too
  // elegant, but it's the most efficient, because in this way there
  // is no additional runtime overhead to handle both cases.

  // We use only the first two elements to decide how the grid is
  // sorted. (The rest of the grid is checked later by an ARTS_ASSERT.)
  // If both values are the same, we still assume the grid is
  // ascending. However, this will lead to an assertion fail later on,
  // because the grid has to be strictly sorted.
  bool ascending = (old_grid[0] <= old_grid[1]);

  if (ascending) {
    // So, old_grid should always be sorted in strictly ascending order.
    // (No duplicate values.)
    // Assert that this is so. This may depend on user input, however,
    // inside this elementary function is not the place to check for
    // that. There must be runtime checks on higher levels to insure
    // that all grids are sorted. The assertion here is just as a last
    // safety check.
    ARTS_ASSERT(is_increasing(old_grid));

    // Limits of extrapolation.
    const Numeric og_min =
        old_grid[0] - extpolfac * (old_grid[1] - old_grid[0]);
    const Numeric og_max =
        old_grid[n_old - 1] +
        extpolfac * (old_grid[n_old - 1] - old_grid[n_old - 2]);

    // We will make no firm assumptions about the new grid. But the case
    // that we have in mind is that it is either also sorted, or at
    // least partially sorted, for example like this:
    // 5 4 3 2 2.5 3 4
    // This kind of sequence should be typical if we interpolate
    // atmospheric fields along a limb propagation path.

    // Let's get some idea where the first point in the new grid is,
    // relative to the old grid. We use linear interpolation between the
    // maximum and the minimum of the old grid for this. (Taking
    // into account the small allowed extrapolation.)
    Numeric frac = (new_grid[0] - og_min) / (og_max - og_min);

    // We are not checking if the value of frac is reasonable,
    // because there is another assertion below to catch the
    // consequences.

    // Initialize current_position:
    Index current_position = (Index)rint(frac * (Numeric)(n_old - 2));

    // The above statement should satisfy
    // 0 <= current_position <= n_old-2
    // Assert that this is indeed the case.
    //       cout << "frac             = " << frac             << "\n";
    //       cout << "current_position = " << current_position << "\n";
    ARTS_ASSERT(0 <= current_position);
    ARTS_ASSERT(current_position <= n_old - 2);

    // The variables lower and upper are used to remember the value of
    // the old grid at current_position and one above current_position:
    Numeric lower = old_grid[current_position];
    Numeric upper = old_grid[current_position + 1];

    // Loop over all points in the new grid:
    for (Index i_new = 0; i_new < n_new; ++i_new) {
      // Get a reference to the current element of gp:
      GridPos& tgp = gp[i_new];
      // And on the current value of the new grid:
      const Numeric tng = new_grid[i_new];

      // Verify that the new grid is within the limits of
      // extrapolation that we have defined by og_min and og_max:
      ARTS_ASSERT(og_min <= tng);
      ARTS_ASSERT(tng <= og_max);

      //           cout << "lower / tng / upper = "
      //                << lower << " / "
      //                << tng   << " / "
      //                << upper << "\n";

      // Is current_position too high?
      // (The current_position>0 condition is there so that the position
      // stays 0 for extrapolation.)
      if (tng < lower && current_position > 0) {
        do {
          --current_position;
          lower = old_grid[current_position];
        } while (tng < lower && current_position > 0);

        upper = old_grid[current_position + 1];

        tgp.idx = current_position;
        tgp.fd[0] = (tng - lower) / (upper - lower);
        tgp.fd[1] = 1.0 - tgp.fd[0];
      } else {
        // Is it too low?
        // (The current_position<n_old condition is there so
        // that uppers stays n_old-1 for extrapolation.)
        if (tng >= upper && current_position < n_old - 2) {
          do {
            ++current_position;
            upper = old_grid[current_position + 1];
          } while (tng >= upper && current_position < n_old - 2);

          lower = old_grid[current_position];

          tgp.idx = current_position;
          tgp.fd[0] = (tng - lower) / (upper - lower);
          tgp.fd[1] = 1.0 - tgp.fd[0];
        } else {
          // None of the other two conditions were true. That means:
          // lower <= tng < upper. The current_position is still
          // valid.

          // SAB 2010-04-28: Note that if a new grid point is exactly on
          // top of an old grid point, then you are now guaranteed to get
          // fd[0] = 0 and fd[1] = 1. (Except at the upper grid end.)

          tgp.idx = current_position;
          tgp.fd[0] = (tng - lower) / (upper - lower);
          tgp.fd[1] = 1.0 - tgp.fd[0];
        }
      }

      //          cout << tgp.idx << " " << tgp.fd[0] << " " << tgp.fd[1] << endl;

      // Safety check to ensure the above promise:
      ARTS_ASSERT(old_grid[tgp.idx] <= tng || tgp.idx == 0);
    }
  } else  //   if (ascending)
  {
    // Now we are in the "descending old grid" part. We do exactly
    // the same as in the other part, just accounting for the
    // different order of things. Comments here refer only to
    // interesting differences from the ascending case. See that
    // case for more general comments.

    // This time ensure strictly descending order:
    ARTS_ASSERT(is_decreasing(old_grid));

    // The max is now the first point, the min the last point!
    // I think the sign is right here...
    const Numeric og_max =
        old_grid[0] - extpolfac * (old_grid[1] - old_grid[0]);
    const Numeric og_min =
        old_grid[n_old - 1] +
        extpolfac * (old_grid[n_old - 1] - old_grid[n_old - 2]);

    // We have to take 1- here, because we are starting from the
    // high end.
    Numeric frac = 1 - (new_grid[0] - og_min) / (og_max - og_min);

    // We are not checking if the value of frac is reasonable,
    // because there is another assertion below to catch the
    // consequences.

    Index current_position = (Index)rint(frac * (Numeric)(n_old - 2));

    // The above statement should satisfy
    // 0 <= current_position <= n_old-2
    // Assert that this is indeed the case.
    //       cout << "frac             = " << frac             << "\n";
    //       cout << "current_position = " << current_position << "\n";
    ARTS_ASSERT(0 <= current_position);
    ARTS_ASSERT(current_position <= n_old - 2);

    // Note, that old_grid[lower] has a higher numerical value than
    // old_grid[upper]!
    Numeric lower = old_grid[current_position];
    Numeric upper = old_grid[current_position + 1];

    for (Index i_new = 0; i_new < n_new; ++i_new) {
      GridPos& tgp = gp[i_new];
      const Numeric tng = new_grid[i_new];

      // Verify that the new grid is within the limits of
      // extrapolation that we have defined by og_min and og_max:
      ARTS_ASSERT(og_min <= tng);
      ARTS_ASSERT(tng <= og_max);

      //           cout << "lower / tng / upper = "
      //                << lower << " / "
      //                << tng   << " / "
      //                << upper << "\n";

      // Is current_position too high? (Sign of comparison changed
      // compared to ascending case!)
      if (tng > lower && current_position > 0) {
        do {
          --current_position;
          lower = old_grid[current_position];
        } while (tng > lower && current_position > 0);

        upper = old_grid[current_position + 1];

        tgp.idx = current_position;
        tgp.fd[0] = (tng - lower) / (upper - lower);
        tgp.fd[1] = 1.0 - tgp.fd[0];
      } else {
        // Is it too low? (Sign of comparison changed
        // compared to ascending case!)
        if (tng <= upper && current_position < n_old - 2) {
          do {
            ++current_position;
            upper = old_grid[current_position + 1];
          } while (tng <= upper && current_position < n_old - 2);

          lower = old_grid[current_position];

          tgp.idx = current_position;
          tgp.fd[0] = (tng - lower) / (upper - lower);
          tgp.fd[1] = 1.0 - tgp.fd[0];
        } else {
          // None of the other two conditions were true. That means:
          // lower >= tng > upper. The current_position is still
          // valid. (Note that upper and lower have switched
          // place compared to the ascending case.)

          // SAB 2010-04-28: Note that if a new grid point is exactly on
          // top of an old grid point, then you are now guaranteed to get
          // fd[0] = 0 and fd[1] = 1. (Except at the upper grid end.)

          tgp.idx = current_position;
          tgp.fd[0] = (tng - lower) / (upper - lower);
          tgp.fd[1] = 1.0 - tgp.fd[0];
        }
      }

      // Safety check to ensure the above promise:
      ARTS_ASSERT(old_grid[tgp.idx] >= tng || tgp.idx == 0);
    }
  }
}

//! gridpos
/*!
   Creates a grid position structure.
  
   This is a gateway to the function for arrays of grid positions, to be
   used for e.g. "red interpolation".

   \retval  gp         The GridPos structure. 
   \param   old_grid   The original grid.
   \param   new_grid   The position where we want to have the interpolated 
                       value.
   \param   extpolfac  Extrapolation factor. Default value is 0.5,
                       which means that extrapolation of half of the
                       last grid distance is allowed.
                       You don't have to specify this.

   \author Patrick Eriksson
   \date   2002-12-31

   \author Stefan Buehler
   \date   2008-03-03
*/
void gridpos(GridPos& gp,
             ConstVectorView old_grid,
             const Numeric& new_grid,
             const Numeric& extpolfac) {
  ArrayOfGridPos agp(1);
  Vector v(1, new_grid);
  gridpos(agp, old_grid, v, extpolfac);
  gridpos_copy(gp, agp[0]);
}

//! gridpos_1to1
/*!
   Creates a position array for a 1-to-1 mapping.

   This function treats interpolaton back to the same grid. The function
   matches: gridpos( gp, grid, grid ), but is faster.

   The function gives gp correct size.

   \retval  gp         Array of GridPos structure. 
   \param   old_grid   The grid.

   \author Patrick Eriksson
   \date   2012-06-22
*/
void gridpos_1to1(ArrayOfGridPos& gp, ConstVectorView grid) {
  const Index n = grid.nelem();
  gp.resize(n);

  for (Index i = 0; i < n - 1; i++) {
    gp[i].idx = i;
    gp[i].fd[0] = 0;
    gp[i].fd[1] = 1;
  }
  gp[n - 1].idx = n - 2;
  gp[n - 1].fd[0] = 1;
  gp[n - 1].fd[1] = 0;
}

//! gridpos_copy
/*!
   Copies the content of a GridPos structure.

   \retval  gp_new   The GridPos structure to be filled. 
   \param   gp_old   The GridPos structure to be copied. 

   \author Patrick Eriksson
   \date   2002-09-17
*/
void gridpos_copy(GridPos& gp_new, const GridPos& gp_old) {
  gp_new.idx = gp_old.idx;
  gp_new.fd[0] = gp_old.fd[0];
  gp_new.fd[1] = gp_old.fd[1];
}

//! fractional_gp
/*!
   Returns the complete fractional grid position.

   The fractional grid position is the sum of grid index and fraction
   distance, and is not an integer value.

   \return         Fractional grid position.
   \param   gp     Input: Grid position structure.

   \author Patrick Eriksson
   \date   2004-09-28
*/
Numeric fractional_gp(const GridPos& gp) {
  return (Numeric(gp.idx) + gp.fd[0]);
}

//! gridpos_check_fd
/*!
   Checks that the fractional distances have a value in the range [0,1].

   This function can be used when you are sure that the fractional distances
   have been calculated correctly, but the limited numerical precision can
   give values below 0 or above 1.

   \retval   gp     Grid position structure.

   \author Patrick Eriksson
   \date   2002-05-21
*/
void gridpos_check_fd(GridPos& gp) {
  // Catch values that "must" be wrong
  ARTS_ASSERT(gp.fd[0] > -FD_TOL);
  ARTS_ASSERT(gp.fd[0] < 1.0 + FD_TOL);
  ARTS_ASSERT(gp.fd[1] > -FD_TOL);
  ARTS_ASSERT(gp.fd[1] < 1.0 + FD_TOL);

  if (gp.fd[0] < 0.0) {
    gp.fd[0] = 0.0;
    gp.fd[1] = 1.0;
  } else if (gp.fd[0] > 1.0) {
    gp.fd[0] = 1.0;
    gp.fd[1] = 0.0;
  }

  if (gp.fd[1] < 0.0) {
    gp.fd[1] = 0.0;
    gp.fd[0] = 1.0;
  } else if (gp.fd[1] > 1.0) {
    gp.fd[1] = 1.0;
    gp.fd[0] = 0.0;
  }
}

//! gridpos_force_end_fd
/*!
   Forces that the fractional distances is set to OK with respect to
   interpolation.

   This function can be called when it is known that a position is exactly
   on a grid point. The fractional distance of the grid position is then
   0 or 1, but rounding errors can give a slightly deviating value.

   The fractional distance is set to 0, beside for end point of the complete
   grid range where is set to 1.

   The difference between this function and gridpos_check_fd is that this
   function is only applicable for end points, while the other function can
   be called for every point.

   \retval   gp     Grid position structure.
   \param     n     Number of points of the grid.

   \author Patrick Eriksson
   \date   2002-05-22
*/
void gridpos_force_end_fd(GridPos& gp, const Index& n) {
  ARTS_ASSERT(gp.idx >= 0);

  // If fd=1, shift to grid index above
  if (gp.fd[0] > 0.5) {
    gp.idx += 1;
  }
  gp.fd[0] = 0;
  gp.fd[1] = 1;

  ARTS_ASSERT(gp.idx < n);

  // End of complete grid range must be handled separately
  if (gp.idx == n - 1) {
    gp.idx -= 1;
    gp.fd[0] = 1;
    gp.fd[1] = 0;
  }
}

//! gridpos_upperend_check
/*!
   A function to handle the upper end after shifts of grid positions

   The function ensures that a shift of grid positions not results in a too
   high gp.idx. That happens for positions exactly at the upper cloudbox
   boundaries.

   \retval   gp     Grid position structure.
   \param    ie     Index of upper limit, shifted

   \author Patrick Eriksson
   \date   2010-11-18
*/
void gridpos_upperend_check(GridPos& gp, const Index& ie) {
  if (gp.idx == ie) {
    ARTS_ASSERT(gp.fd[0] < 0.005);  // To capture obviously bad cases
    gp.idx -= 1;
    gp.fd[0] = 1.0;
    gp.fd[1] = 0.0;
  }
}

//! gridpos_upperend_check
/*!
   A function to handle the upper end after shifts of grid positions

   The function ensures that a shift of grid positions not results in a too
   high gp.idx. That happens for positions exactly at the upper cloudbox
   boundaries.

   \retval   gp     Array of grid position structure.
   \param    ie     Index of upper limit, shifted

   \author Patrick Eriksson
   \date   2011-02-28
*/
void gridpos_upperend_check(ArrayOfGridPos& gp, const Index& ie) {
  for (Index i = 0; i < gp.nelem(); i++) {
    if (gp[i].idx == ie) {
      ARTS_ASSERT(gp[i].fd[0] < 0.005);  // To capture obviously bad cases
      gp[i].idx -= 1;
      gp[i].fd[0] = 1.0;
      gp[i].fd[1] = 0.0;
    }
  }
}

//! Grid position matching a grid of length 1
/*!
  The function sets *gp* to the expected values if there would exist a grid of
  length 1.

  Grids of length 1 are normally not allowed, but the jacobian part makes use
  of such grids to flag a completely flat field in the diemsnion of concern. 

  \param[in/out] gp   A grid positions.

  \author Patrick Eriksson 
  \date   2015-09-10
*/
void gp4length1grid(ArrayOfGridPos& gp) {
  for (Index i = 0; i < gp.nelem(); i++) {
    gp[i].idx = 0;
    gp[i].fd[0] = 0;
    gp[i].fd[1] = 1;
  }
}

//! is_gridpos_at_index_i
/*!
   Determines if a grid position is at a given grid index.

   \return           True if at index i, else false.
   \param   gp       Grid position structure.
   \param   i        The grid index of interest.
   \param   strict   If true, fractional distances are demanded to be exactly 0
                     or 1. Otherwise FD_TOL is considered. 

   \author Patrick Eriksson
   \date   2002-05-22
*/
bool is_gridpos_at_index_i(const GridPos& gp,
                           const Index& i,
                           const bool& strict) {
  if (strict) {
    // Assume that gridpos_force_end_fd has been used. The expression 0==0
    // should be safer than 1==1.
    if (gp.idx == i && gp.fd[0] == 0) {
      return true;
    } else if (gp.idx == i - 1 && gp.fd[1] == 0) {
      return true;
    }
  } else {
    if (gp.idx == i && gp.fd[0] < FD_TOL) {
      return true;
    } else if (gp.idx == i - 1 && gp.fd[1] < FD_TOL) {
      return true;
    }
  }
  return false;
}

//! gridpos2gridrange
/*!
   Determines which grid range that is of interest for a given grid position.

   The purpose of the function is to determine which two grid values that
   surround the given point. The index of the lower grid value is returned.

   For a point exactly on a grid value it is not clear if it is the range 
   below or above that is of interest. The input argument upward is used to 
   resolve such cases, where upward == 1 means that it is the range above
   that is of interest.

   \param   gp      Grid position structure.
   \param   upwards Direction of interest, see above.
   \return          The index of the lower end of the grid range.

   \author Patrick Eriksson
   \date   2002-05-20
*/
Index gridpos2gridrange(const GridPos& gp, const bool& upwards) {
  ARTS_ASSERT(gp.fd[0] >= 0);
  ARTS_ASSERT(gp.fd[1] >= 0);

  // Not at a grid point
  if (gp.fd[0] > 0 && gp.fd[1] > 0) {
    return gp.idx;
  }

  // Fractional distance 0
  else if (gp.fd[0] == 0) {
    if (upwards)
      return gp.idx;
    else
      return gp.idx - 1;
  }

  // Fractional distance 1
  else {
    if (upwards)
      return gp.idx + 1;
    else
      return gp.idx;
  }
}

////////////////////////////////////////////////////////////////////////////
//                      Red Interpolation
////////////////////////////////////////////////////////////////////////////

//! Red 1D interpolation weights.
/*! 
  "Red" interpolation returns just a scalar, so the weights are stored
  in a Vector.   

  The length of itw must be consistent with the dimension of the
  field to be interpolated (2^n).

  \retval itw Interpolation weights.
  \param  tc  The grid position for the column dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
void interpweights(VectorView itw, const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 2));  // We must store 2 interpolation
                            // weights.

  // Interpolation weights are stored in this order (l=lower
  // u=upper, c=column):
  // 1. l-c
  // 2. u-c

  Index iti = 0;

  // We could use a straight-forward for loop here:
  //
  //       for ( Index c=1; c>=0; --c )
  //    {
  //      ti[iti] = tc.fd[c];
  //      ++iti;
  //    }
  //
  // However, it is slightly faster to use pointers (I tried it,
  // the speed gain is about 20%). So we should write something
  // like:
  //
  //       for ( const Numeric* c=&tc.fd[1]; c>=&tc.fd[0]; --c )
  //    {
  //      ti[iti] = *c;
  //      ++iti;
  //    }
  //
  // For higher dimensions we have to nest these loops. To avoid
  // typos and safe typing, I use the LOOPIT macro, which expands
  // to the for loop above. Note: NO SEMICOLON AFTER THE LOOPIT
  // COMMAND!

  LOOPIT(c) {
    itw[iti] = *c;
    ++iti;
  }
}

//! Red 2D interpolation weights.
/*! 
  "Red" interpolation returns just a scalar, so the weights are stored
  in a Vector.   

  The length of itw must be consistent with the dimension of the
  field to be interpolated (2^n).

  \retval itw Interpolation weights.
  \param tr   The grid position for the row dimension.
  \param tc   The grid position for the column dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
void interpweights(VectorView itw, const GridPos& tr, const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 4));  // We must store 4 interpolation
                            // weights.
  Index iti = 0;

  LOOPIT(r)
  LOOPIT(c) {
    itw[iti] = (*r) * (*c);
    ++iti;
  }
}

//! Red 3D interpolation weights.
/*! 
  "Red" interpolation returns just a scalar, so the weights are stored
  in a Vector.   

  The length of itw must be consistent with the dimension of the
  field to be interpolated (2^n).

  \retval itw Interpolation weights.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
void interpweights(VectorView itw,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 8));  // We must store 8 interpolation
                            // weights.
  Index iti = 0;

  LOOPIT(p)
  LOOPIT(r)
  LOOPIT(c) {
    itw[iti] = (*p) * (*r) * (*c);
    ++iti;
  }
}

//! Red 4D interpolation weights.
/*! 
  "Red" interpolation returns just a scalar, so the weights are stored
  in a Vector.   

  The length of itw must be consistent with the dimension of the
  field to be interpolated (2^n).

  \retval itw Interpolation weights.
  \param tb   The grid position for the book    dimension.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
void interpweights(VectorView itw,
                   const GridPos& tb,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 16));  // We must store 16 interpolation
                             // weights.
  Index iti = 0;

  LOOPIT(b)
  LOOPIT(p)
  LOOPIT(r)
  LOOPIT(c) {
    itw[iti] = (*b) * (*p) * (*r) * (*c);
    ++iti;
  }
}

//! Red 5D interpolation weights.
/*! 
  "Red" interpolation returns just a scalar, so the weights are stored
  in a Vector.   

  The length of itw must be consistent with the dimension of the
  field to be interpolated (2^n).

  \retval itw Interpolation weights.
  \param ts   The grid position for the shelf   dimension.
  \param tb   The grid position for the book    dimension.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
void interpweights(VectorView itw,
                   const GridPos& ts,
                   const GridPos& tb,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 32));  // We must store 32 interpolation
                             // weights.
  Index iti = 0;

  LOOPIT(s)
  LOOPIT(b)
  LOOPIT(p)
  LOOPIT(r)
  LOOPIT(c) {
    itw[iti] = (*s) * (*b) * (*p) * (*r) * (*c);
    ++iti;
  }
}

//! Red 6D interpolation weights.
/*! 
  "Red" interpolation returns just a scalar, so the weights are stored
  in a Vector.   

  The length of itw must be consistent with the dimension of the
  field to be interpolated (2^n).

  \retval itw Interpolation weights.
  \param  tv  The grid position for the vitrine dimension.
  \param  ts  The grid position for the shelf   dimension.
  \param  tb  The grid position for the book    dimension.
  \param  tp  The grid position for the page    dimension.
  \param  tr  The grid position for the row     dimension.
  \param  tc  The grid position for the column  dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
void interpweights(VectorView itw,
                   const GridPos& tv,
                   const GridPos& ts,
                   const GridPos& tb,
                   const GridPos& tp,
                   const GridPos& tr,
                   const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 64));  // We must store 64 interpolation
                             // weights.
  Index iti = 0;

  LOOPIT(v)
  LOOPIT(s)
  LOOPIT(b)
  LOOPIT(p)
  LOOPIT(r)
  LOOPIT(c) {
    itw[iti] = (*v) * (*s) * (*b) * (*p) * (*r) * (*c);
    ++iti;
  }
}

//! Red 1D Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param tc   The grid position for the column dimension.

  \return Interpolated value.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
Numeric interp(ConstVectorView itw, ConstVectorView a, const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 2));  // We need 2 interpolation
                            // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  ARTS_ASSERT(is_same_within_epsilon(sum(itw), 1, sum_check_epsilon));

  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for (Index c = 0; c < 2; ++c) {
    tia += a[tc.idx + c] * itw[iti];
    ++iti;
  }

  return tia;
}

//! Red 2D Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \return Interpolated value.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
Numeric interp(ConstVectorView itw,
               ConstMatrixView a,
               const GridPos& tr,
               const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 4));  // We need 4 interpolation
                            // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  ARTS_ASSERT(is_same_within_epsilon(sum(itw), 1, sum_check_epsilon));

  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for (Index r = 0; r < 2; ++r)
    for (Index c = 0; c < 2; ++c) {
      tia += a(tr.idx + r, tc.idx + c) * itw[iti];
      ++iti;
    }

  return tia;
}

//! Red 3D Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \return Interpolated value.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
Numeric interp(ConstVectorView itw,
               ConstTensor3View a,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 8));  // We need 8 interpolation
                            // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  ARTS_ASSERT(is_same_within_epsilon(sum(itw), 1, sum_check_epsilon));

  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for (Index p = 0; p < 2; ++p)
    for (Index r = 0; r < 2; ++r)
      for (Index c = 0; c < 2; ++c) {
        tia += a(tp.idx + p, tr.idx + r, tc.idx + c) * itw[iti];
        ++iti;
      }

  return tia;
}

//! Red 4D Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param tb   The grid position for the book    dimension.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \return Interpolated value.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
Numeric interp(ConstVectorView itw,
               ConstTensor4View a,
               const GridPos& tb,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 16));  // We need 16 interpolation
                             // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  ARTS_ASSERT(is_same_within_epsilon(sum(itw), 1, sum_check_epsilon));

  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for (Index b = 0; b < 2; ++b)
    for (Index p = 0; p < 2; ++p)
      for (Index r = 0; r < 2; ++r)
        for (Index c = 0; c < 2; ++c) {
          tia += a(tb.idx + b, tp.idx + p, tr.idx + r, tc.idx + c) *
                 itw[iti];
          ++iti;
        }

  return tia;
}

//! Red 5D Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param ts   The grid position for the shelf   dimension.
  \param tb   The grid position for the book    dimension.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \return Interpolated value.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
Numeric interp(ConstVectorView itw,
               ConstTensor5View a,
               const GridPos& ts,
               const GridPos& tb,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 32));  // We need 32 interpolation
                             // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  ARTS_ASSERT(is_same_within_epsilon(sum(itw), 1, sum_check_epsilon));

  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for (Index s = 0; s < 2; ++s)
    for (Index b = 0; b < 2; ++b)
      for (Index p = 0; p < 2; ++p)
        for (Index r = 0; r < 2; ++r)
          for (Index c = 0; c < 2; ++c) {
            tia += a(ts.idx + s,
                         tb.idx + b,
                         tp.idx + p,
                         tr.idx + r,
                         tc.idx + c) *
                   itw[iti];
            ++iti;
          }

  return tia;
}

//! Red 6D Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (2^n).

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param tv   The grid position for the vitrine dimension.
  \param ts   The grid position for the shelf   dimension.
  \param tb   The grid position for the book    dimension.
  \param tp   The grid position for the page    dimension.
  \param tr   The grid position for the row     dimension.
  \param tc   The grid position for the column  dimension.

  \return Interpolated value.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri Jun 28 10:53:32 2002
*/
Numeric interp(ConstVectorView itw,
               ConstTensor6View a,
               const GridPos& tv,
               const GridPos& ts,
               const GridPos& tb,
               const GridPos& tp,
               const GridPos& tr,
               const GridPos& tc) {
  ARTS_ASSERT(is_size(itw, 64));  // We need 64 interpolation
                             // weights.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  ARTS_ASSERT(is_same_within_epsilon(sum(itw), 1, sum_check_epsilon));

  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for (Index v = 0; v < 2; ++v)
    for (Index s = 0; s < 2; ++s)
      for (Index b = 0; b < 2; ++b)
        for (Index p = 0; p < 2; ++p)
          for (Index r = 0; r < 2; ++r)
            for (Index c = 0; c < 2; ++c) {
              tia += a(tv.idx + v,
                           ts.idx + s,
                           tb.idx + b,
                           tp.idx + p,
                           tr.idx + r,
                           tc.idx + c) *
                     itw[iti];
              ++iti;
            }

  return tia;
}

////////////////////////////////////////////////////////////////////////////
//                      Blue interpolation
////////////////////////////////////////////////////////////////////////////

//! Compute 1D interpolation weights.
/*! 
  For this 1D case there is no distinction between "blue" and "green"
  type interpolation.

  The dimensions of itw must be consistent with cgp.

  Note that we still do not need the actual field for this step.

  
  \retval itw Interpolation weights.
  \param  cgp The grid position Array for the column dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri May  3 08:55:51 2002
*/
void interpweights(MatrixView itw, const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(itw, n, 2));  // We must store 2 interpolation
                               // weights for each position.

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    //        {
    //          ti[iti] = tc.fd[c];
    //          ++iti;
    //        }
    //
    // However, it is slightly faster to use pointers (I tried it,
    // the speed gain is about 20%). So we should write something
    // like:
    //
    //       for ( const Numeric* c=&tc.fd[1]; c>=&tc.fd[0]; --c )
    //        {
    //          ti[iti] = *c;
    //          ++iti;
    //        }
    //
    // For higher dimensions we have to nest these loops. To avoid
    // typos and safe typing, I use the LOOPIT macro, which expands
    // to the for loop above. Note: NO SEMICOLON AFTER THE LOOPIT
    // COMMAND!

    LOOPIT(c) {
      itw(i, iti) = *c;
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

 \retval itw Interpolation weights.
 \param  rgp The grid position Array for the row dimension.
 \param  cgp The grid position Array for the column dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(MatrixView itw,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(rgp, n));     // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 4));  // We must store 4 interpolation
                               // weights for each position.

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    LOOPIT(c) {
      itw(i, iti) = (*r) * (*c);
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

 \retval itw Interpolation weights.
 \param  pgp The grid position Array for the page dimension.
 \param  rgp The grid position Array for the row dimension.
 \param  cgp The grid position Array for the column dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(MatrixView itw,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(pgp, n));     // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));     // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 8));  // We must store 8 interpolation
                               // weights for each position.

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
    // Current grid positions:
    const GridPos& tp = pgp[i];
    const GridPos& tr = rgp[i];
    const GridPos& tc = cgp[i];

    Index iti = 0;
    LOOPIT(p)
    LOOPIT(r)
    LOOPIT(c) {
      itw(i, iti) = (*p) * (*r) * (*c);
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

 \retval itw Interpolation weights.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(MatrixView itw,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(bgp, n));      // bgp must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));      // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));      // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 16));  // We must store 16 interpolation
                                // weights for each position.

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
    // Current grid positions:
    const GridPos& tb = bgp[i];
    const GridPos& tp = pgp[i];
    const GridPos& tr = rgp[i];
    const GridPos& tc = cgp[i];

    Index iti = 0;
    LOOPIT(b)
    LOOPIT(p)
    LOOPIT(r)
    LOOPIT(c) {
      itw(i, iti) = (*b) * (*p) * (*r) * (*c);
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

 \retval itw Interpolation weights.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(MatrixView itw,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(sgp, n));      // sgp must have same size as cgp.
  ARTS_ASSERT(is_size(bgp, n));      // bgp must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));      // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));      // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 32));  // We must store 32 interpolation
                                // weights for each position.

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    LOOPIT(c) {
      itw(i, iti) = (*s) * (*b) * (*p) * (*r) * (*c);
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

 \retval itw Interpolation weights.
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(MatrixView itw,
                   const ArrayOfGridPos& vgp,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(vgp, n));      // vgp must have same size as cgp.
  ARTS_ASSERT(is_size(sgp, n));      // sgp must have same size as cgp.
  ARTS_ASSERT(is_size(bgp, n));      // bgp must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));      // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));      // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 64));  // We must store 64 interpolation
                                // weights for each position.

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    LOOPIT(c) {
      itw(i, iti) = (*v) * (*s) * (*b) * (*p) * (*r) * (*c);
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

  \retval ia  Vector containing the interpolated field values.
  \param itw Interpolation weights.
  \param a   The field to interpolate.
  \param cgp The grid position Array for the column dimension.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri May  3 08:55:51 2002
*/
void interp(VectorView ia,
            ConstMatrixView itw,
            ConstVectorView a,
            const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(ia, n));      //  ia must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 2));  // We need 2 interpolation
                               // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(
      is_same_within_epsilon(sum(itw(0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
    // Current grid positions:
    const GridPos& tc = cgp[i];

    // Get handle to current element of output vector and initialize
    // it to zero:
    Numeric& tia = ia[i];
    tia = 0;

    Index iti = 0;
    for (Index c = 0; c < 2; ++c) {
      ARTS_ASSERT(tc.idx + c < a.nelem());  // Temporary !?
      tia += a(tc.idx + c) * itw(i, iti);
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

 \retval ia  Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param rgp The grid position Array for the row    dimension.
 \param cgp The grid position Array for the column dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(VectorView ia,
            ConstMatrixView itw,
            ConstMatrixView a,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(ia, n));      //  ia must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));     // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 4));  // We need 4 interpolation
                               // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(
      is_same_within_epsilon(sum(itw(0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
    // Current grid positions:
    const GridPos& tr = rgp[i];
    const GridPos& tc = cgp[i];

    // Get handle to current element of output vector and initialize
    // it to zero:
    Numeric& tia = ia[i];
    tia = 0;

    Index iti = 0;
    for (Index r = 0; r < 2; ++r)
      for (Index c = 0; c < 2; ++c) {
        tia += a(tr.idx + r, tc.idx + c) * itw(i, iti);
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

 \retval ia  Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(VectorView ia,
            ConstMatrixView itw,
            ConstTensor3View a,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(ia, n));      //  ia must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));     // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));     // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 8));  // We need 8 interpolation
                               // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(
      is_same_within_epsilon(sum(itw(0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
    // Current grid positions:
    const GridPos& tp = pgp[i];
    const GridPos& tr = rgp[i];
    const GridPos& tc = cgp[i];

    // Get handle to current element of output vector and initialize
    // it to zero:
    Numeric& tia = ia[i];
    tia = 0;

    Index iti = 0;
    for (Index p = 0; p < 2; ++p)
      for (Index r = 0; r < 2; ++r)
        for (Index c = 0; c < 2; ++c) {
          tia += a(tp.idx + p, tr.idx + r, tc.idx + c) * itw(i, iti);
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

 \retval ia  Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(VectorView ia,
            ConstMatrixView itw,
            ConstTensor4View a,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(ia, n));       //  ia must have same size as cgp.
  ARTS_ASSERT(is_size(bgp, n));      // bgp must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));      // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));      // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 16));  // We need 16 interpolation
                                // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(
      is_same_within_epsilon(sum(itw(0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    for (Index b = 0; b < 2; ++b)
      for (Index p = 0; p < 2; ++p)
        for (Index r = 0; r < 2; ++r)
          for (Index c = 0; c < 2; ++c) {
            tia += a(tb.idx + b, tp.idx + p, tr.idx + r, tc.idx + c) *
                   itw(i, iti);
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

 \retval ia  Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(VectorView ia,
            ConstMatrixView itw,
            ConstTensor5View a,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(ia, n));       //  ia must have same size as cgp.
  ARTS_ASSERT(is_size(sgp, n));      // sgp must have same size as cgp.
  ARTS_ASSERT(is_size(bgp, n));      // bgp must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));      // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));      // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 32));  // We need 32 interpolation
                                // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(
      is_same_within_epsilon(sum(itw(0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    for (Index s = 0; s < 2; ++s)
      for (Index b = 0; b < 2; ++b)
        for (Index p = 0; p < 2; ++p)
          for (Index r = 0; r < 2; ++r)
            for (Index c = 0; c < 2; ++c) {
              tia += a(ts.idx + s,
                           tb.idx + b,
                           tp.idx + p,
                           tr.idx + r,
                           tc.idx + c) *
                     itw(i, iti);
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

 \retval ia  Vector containing the interpolated field values.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(VectorView ia,
            ConstMatrixView itw,
            ConstTensor6View a,
            const ArrayOfGridPos& vgp,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index n = cgp.nelem();
  ARTS_ASSERT(is_size(ia, n));       //  ia must have same size as cgp.
  ARTS_ASSERT(is_size(vgp, n));      // vgp must have same size as cgp.
  ARTS_ASSERT(is_size(sgp, n));      // sgp must have same size as cgp.
  ARTS_ASSERT(is_size(bgp, n));      // bgp must have same size as cgp.
  ARTS_ASSERT(is_size(pgp, n));      // pgp must have same size as cgp.
  ARTS_ASSERT(is_size(rgp, n));      // rgp must have same size as cgp.
  ARTS_ASSERT(is_size(itw, n, 64));  // We need 64 interpolation
                                // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(
      is_same_within_epsilon(sum(itw(0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the sequence:
  for (Index i = 0; i < n; ++i) {
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
    for (Index v = 0; v < 2; ++v)
      for (Index s = 0; s < 2; ++s)
        for (Index b = 0; b < 2; ++b)
          for (Index p = 0; p < 2; ++p)
            for (Index r = 0; r < 2; ++r)
              for (Index c = 0; c < 2; ++c) {
                tia += a(tv.idx + v,
                             ts.idx + s,
                             tb.idx + b,
                             tp.idx + p,
                             tr.idx + r,
                             tc.idx + c) *
                       itw(i, iti);
                ++iti;
              }
  }
}

////////////////////////////////////////////////////////////////////////////
//                      Green interpolation
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

 \retval itw Interpolation weights
 \param rgp The grid position Array for the row dimension.
 \param cgp The grid position Array for the column dimension.
 
 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(Tensor3View itw,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  ARTS_ASSERT(is_size(itw, nr, nc, 4));  // We must store 4 interpolation
                                    // weights for each position.

  // We have to loop all the points in the new grid:
  for (Index ir = 0; ir < nr; ++ir) {
    // Current grid position:
    const GridPos& tr = rgp[ir];

    for (Index ic = 0; ic < nc; ++ic) {
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
      LOOPIT(c) {
        itw(ir, ic, iti) = (*r) * (*c);
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

 \retval itw Interpolation weights
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(Tensor4View itw,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 8 interpolation weights for each position:
  ARTS_ASSERT(is_size(itw, np, nr, nc, 8));

  // We have to loop all the points in the new grid:
  for (Index ip = 0; ip < np; ++ip) {
    const GridPos& tp = pgp[ip];
    for (Index ir = 0; ir < nr; ++ir) {
      const GridPos& tr = rgp[ir];
      for (Index ic = 0; ic < nc; ++ic) {
        const GridPos& tc = cgp[ic];

        Index iti = 0;

        LOOPIT(p)
        LOOPIT(r)
        LOOPIT(c) {
          itw(ip, ir, ic, iti) = (*p) * (*r) * (*c);
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

 \retval itw Interpolation weights
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(Tensor5View itw,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 16 interpolation weights for each position:
  ARTS_ASSERT(is_size(itw, nb, np, nr, nc, 16));

  // We have to loop all the points in the new grid:
  for (Index ib = 0; ib < nb; ++ib) {
    const GridPos& tb = bgp[ib];
    for (Index ip = 0; ip < np; ++ip) {
      const GridPos& tp = pgp[ip];
      for (Index ir = 0; ir < nr; ++ir) {
        const GridPos& tr = rgp[ir];
        for (Index ic = 0; ic < nc; ++ic) {
          const GridPos& tc = cgp[ic];

          Index iti = 0;

          LOOPIT(b)
          LOOPIT(p)
          LOOPIT(r)
          LOOPIT(c) {
            itw(ib, ip, ir, ic, iti) = (*b) * (*p) * (*r) * (*c);
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

 \retval itw Interpolation weights
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(Tensor6View itw,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 32 interpolation weights for each position:
  ARTS_ASSERT(is_size(itw, ns, nb, np, nr, nc, 32));

  // We have to loop all the points in the new grid:
  for (Index is = 0; is < ns; ++is) {
    const GridPos& ts = sgp[is];
    for (Index ib = 0; ib < nb; ++ib) {
      const GridPos& tb = bgp[ib];
      for (Index ip = 0; ip < np; ++ip) {
        const GridPos& tp = pgp[ip];
        for (Index ir = 0; ir < nr; ++ir) {
          const GridPos& tr = rgp[ir];
          for (Index ic = 0; ic < nc; ++ic) {
            const GridPos& tc = cgp[ic];

            Index iti = 0;

            LOOPIT(s)
            LOOPIT(b)
            LOOPIT(p)
            LOOPIT(r)
            LOOPIT(c) {
              itw(is, ib, ip, ir, ic, iti) =
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

 \retval itw Interpolation weights
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.
 
 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interpweights(Tensor7View itw,
                   const ArrayOfGridPos& vgp,
                   const ArrayOfGridPos& sgp,
                   const ArrayOfGridPos& bgp,
                   const ArrayOfGridPos& pgp,
                   const ArrayOfGridPos& rgp,
                   const ArrayOfGridPos& cgp) {
  Index nv = vgp.nelem();
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  // We must store 64 interpolation weights for each position:
  ARTS_ASSERT(is_size(itw, nv, ns, nb, np, nr, nc, 64));

  // We have to loop all the points in the new grid:
  for (Index iv = 0; iv < nv; ++iv) {
    const GridPos& tv = vgp[iv];
    for (Index is = 0; is < ns; ++is) {
      const GridPos& ts = sgp[is];
      for (Index ib = 0; ib < nb; ++ib) {
        const GridPos& tb = bgp[ib];
        for (Index ip = 0; ip < np; ++ip) {
          const GridPos& tp = pgp[ip];
          for (Index ir = 0; ir < nr; ++ir) {
            const GridPos& tr = rgp[ir];
            for (Index ic = 0; ic < nc; ++ic) {
              const GridPos& tc = cgp[ic];

              Index iti = 0;

              LOOPIT(v)
              LOOPIT(s)
              LOOPIT(b)
              LOOPIT(p)
              LOOPIT(r)
              LOOPIT(c) {
                itw(iv, is, ib, ip, ir, ic, iti) =
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

 \retval ia  Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param rgp The grid position Array for the row dimension.
 \param cgp The grid position Array for the column dimension.
 
 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(MatrixView ia,
            ConstTensor3View itw,
            ConstMatrixView a,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  ARTS_ASSERT(is_size(ia, nr, nc));
  ARTS_ASSERT(is_size(itw, nr, nc, 4));  // We need 4 interpolation
                                    // weights for each position.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(is_same_within_epsilon(
      sum(itw(0, 0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the new grid:
  for (Index ir = 0; ir < nr; ++ir) {
    // Current grid position:
    const GridPos& tr = rgp[ir];

    for (Index ic = 0; ic < nc; ++ic) {
      // Current grid position:
      const GridPos& tc = cgp[ic];

      // Get handle to current element of output tensor and initialize
      // it to zero:
      Numeric& tia = ia(ir, ic);
      tia = 0;

      Index iti = 0;
      for (Index r = 0; r < 2; ++r)
        for (Index c = 0; c < 2; ++c) {
          ARTS_ASSERT(tr.idx + r < a.nrows());  // Temporary !?
          ARTS_ASSERT(tc.idx + c < a.ncols());  // Temporary !?
          tia += a(tr.idx + r, tc.idx + c) * itw(ir, ic, iti);
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

 \retval ia  Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(Tensor3View ia,
            ConstTensor4View itw,
            ConstTensor3View a,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  ARTS_ASSERT(is_size(ia, np, nr, nc));
  ARTS_ASSERT(is_size(itw, np, nr, nc, 8));

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(is_same_within_epsilon(
      sum(itw(0, 0, 0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the new grid:
  for (Index ip = 0; ip < np; ++ip) {
    const GridPos& tp = pgp[ip];
    for (Index ir = 0; ir < nr; ++ir) {
      const GridPos& tr = rgp[ir];
      for (Index ic = 0; ic < nc; ++ic) {
        // Current grid position:
        const GridPos& tc = cgp[ic];

        // Get handle to current element of output tensor and
        // initialize it to zero:
        Numeric& tia = ia(ip, ir, ic);
        tia = 0;

        Index iti = 0;
        for (Index p = 0; p < 2; ++p)
          for (Index r = 0; r < 2; ++r)
            for (Index c = 0; c < 2; ++c) {
              ARTS_ASSERT(tp.idx + p < a.npages());  // Temporary !?
              ARTS_ASSERT(tr.idx + r < a.nrows());   // Temporary !?
              ARTS_ASSERT(tc.idx + c < a.ncols());   // Temporary !?
              tia += a(tp.idx + p, tr.idx + r, tc.idx + c) *
                     itw(ip, ir, ic, iti);
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

 \retval ia  Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(Tensor4View ia,
            ConstTensor5View itw,
            ConstTensor4View a,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  ARTS_ASSERT(is_size(ia, nb, np, nr, nc));
  ARTS_ASSERT(is_size(itw, nb, np, nr, nc, 16));

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(is_same_within_epsilon(
      sum(itw(0, 0, 0, 0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the new grid:
  for (Index ib = 0; ib < nb; ++ib) {
    const GridPos& tb = bgp[ib];
    for (Index ip = 0; ip < np; ++ip) {
      const GridPos& tp = pgp[ip];
      for (Index ir = 0; ir < nr; ++ir) {
        const GridPos& tr = rgp[ir];
        for (Index ic = 0; ic < nc; ++ic) {
          // Current grid position:
          const GridPos& tc = cgp[ic];

          // Get handle to current element of output tensor and
          // initialize it to zero:
          Numeric& tia = ia(ib, ip, ir, ic);
          tia = 0;

          Index iti = 0;
          for (Index b = 0; b < 2; ++b)
            for (Index p = 0; p < 2; ++p)
              for (Index r = 0; r < 2; ++r)
                for (Index c = 0; c < 2; ++c) {
                  tia += a(tb.idx + b, tp.idx + p, tr.idx + r, tc.idx + c) *
                         itw(ib, ip, ir, ic, iti);
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

 \retval ia  Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(Tensor5View ia,
            ConstTensor6View itw,
            ConstTensor5View a,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  ARTS_ASSERT(is_size(ia, ns, nb, np, nr, nc));
  ARTS_ASSERT(is_size(itw, ns, nb, np, nr, nc, 32));

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(is_same_within_epsilon(
      sum(itw(0, 0, 0, 0, 0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the new grid:
  for (Index is = 0; is < ns; ++is) {
    const GridPos& ts = sgp[is];
    for (Index ib = 0; ib < nb; ++ib) {
      const GridPos& tb = bgp[ib];
      for (Index ip = 0; ip < np; ++ip) {
        const GridPos& tp = pgp[ip];
        for (Index ir = 0; ir < nr; ++ir) {
          const GridPos& tr = rgp[ir];
          for (Index ic = 0; ic < nc; ++ic) {
            // Current grid position:
            const GridPos& tc = cgp[ic];

            // Get handle to current element of output tensor and
            // initialize it to zero:
            Numeric& tia = ia(is, ib, ip, ir, ic);
            tia = 0;

            Index iti = 0;
            for (Index s = 0; s < 2; ++s)
              for (Index b = 0; b < 2; ++b)
                for (Index p = 0; p < 2; ++p)
                  for (Index r = 0; r < 2; ++r)
                    for (Index c = 0; c < 2; ++c) {
                      tia += a(ts.idx + s,
                                   tb.idx + b,
                                   tp.idx + p,
                                   tr.idx + r,
                                   tc.idx + c) *
                             itw(is, ib, ip, ir, ic, iti);
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

 \retval ia  Interpolated field.
 \param itw Interpolation weights.
 \param a   The field to interpolate.
 \param vgp The grid position Array for the vitrine dimension.
 \param sgp The grid position Array for the shelf   dimension.
 \param bgp The grid position Array for the book    dimension.
 \param pgp The grid position Array for the page    dimension.
 \param rgp The grid position Array for the row     dimension.
 \param cgp The grid position Array for the column  dimension.

 \author Stefan Buehler <sbuehler@ltu.se>
 \date   Fri May  3 08:55:51 2002
*/
void interp(Tensor6View ia,
            ConstTensor7View itw,
            ConstTensor6View a,
            const ArrayOfGridPos& vgp,
            const ArrayOfGridPos& sgp,
            const ArrayOfGridPos& bgp,
            const ArrayOfGridPos& pgp,
            const ArrayOfGridPos& rgp,
            const ArrayOfGridPos& cgp) {
  Index nv = vgp.nelem();
  Index ns = sgp.nelem();
  Index nb = bgp.nelem();
  Index np = pgp.nelem();
  Index nr = rgp.nelem();
  Index nc = cgp.nelem();
  ARTS_ASSERT(is_size(ia, nv, ns, nb, np, nr, nc));
  ARTS_ASSERT(is_size(itw, nv, ns, nb, np, nr, nc, 64));

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  ARTS_ASSERT(is_same_within_epsilon(
      sum(itw(0, 0, 0, 0, 0, 0, Range(joker))), 1, sum_check_epsilon));

  // We have to loop all the points in the new grid:
  for (Index iv = 0; iv < nv; ++iv) {
    const GridPos& tv = vgp[iv];
    for (Index is = 0; is < ns; ++is) {
      const GridPos& ts = sgp[is];
      for (Index ib = 0; ib < nb; ++ib) {
        const GridPos& tb = bgp[ib];
        for (Index ip = 0; ip < np; ++ip) {
          const GridPos& tp = pgp[ip];
          for (Index ir = 0; ir < nr; ++ir) {
            const GridPos& tr = rgp[ir];
            for (Index ic = 0; ic < nc; ++ic) {
              // Current grid position:
              const GridPos& tc = cgp[ic];

              // Get handle to current element of output tensor and
              // initialize it to zero:
              Numeric& tia = ia(iv, is, ib, ip, ir, ic);
              tia = 0;

              Index iti = 0;
              for (Index v = 0; v < 2; ++v)
                for (Index s = 0; s < 2; ++s)
                  for (Index b = 0; b < 2; ++b)
                    for (Index p = 0; p < 2; ++p)
                      for (Index r = 0; r < 2; ++r)
                        for (Index c = 0; c < 2; ++c) {
                          tia += a(tv.idx + v,
                                       ts.idx + s,
                                       tb.idx + b,
                                       tp.idx + p,
                                       tr.idx + r,
                                       tc.idx + c) *
                                 itw(iv, is, ib, ip, ir, ic, iti);
                          ++iti;
                        }
            }
          }
        }
      }
    }
  }
}


//! Polynomial interpolation.
/*! 
  This function performs a polinomial interpolation. Given arrays xa(n) and
  ya(n), and a given value x, this function returns a value y and an error 
  estimate dy.
  This function is (almost) copied from: 
  Numerical Recipies in C, pages 108-110.
   
  \param y_int interpolated value
  \param dy_int error estimate
  \param xa original grid (n elements)
  \param ya corresponding values (n elements)
  \param n order of polynom
  \param x requested grid point

  \return interpolated value

  \author Claudia Emde
  \date 2004-03-12
*/
void polint(Numeric& y_int,
            Numeric& dy_int,
            ConstVectorView xa,
            ConstVectorView ya,
            const Index& n,
            const Numeric& x) {
  Index ns = 1;
  Numeric den, dif, dift, ho, hp, w;

  dif = abs(x - xa[0]);

  Vector c(n);
  Vector d(n);

  // Find the index of the closest table entry
  for (Index i = 0; i < n; i++) {
    if ((dift = abs(x - xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    // Initialize c and d
    c[i] = ya[i];
    d[i] = ya[i];
  }
  // Initial approximation to y
  y_int = ya[ns--];

  for (Index m = 1; m < n; m++) {
    for (Index i = 0; i < n - m; i++) {
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = c[i + 1] - d[i];
      den = ho - hp;
      // This error occurs when two input xa's are identical.
      ARTS_ASSERT(den != 0.);
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    y_int += (dy_int = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
  }
}

//! Polynomial interpolation.
/*! 
  This function performs a polynomial interpolation. Given two vectors x, y 
  the value of y corresponding to x_i is returned. The function uses the
  common gridpositions (which are also used for linear interpolation).
   
  \param x  original grid
  \param y  values corresponding to x-values
  \param x_i a value x, for which the corresponding y value is requested
  \param gp grid position of x_i in relation to x

  \return interpolated value

  \author Claudia Emde
  \date 2004-03-17
*/
Numeric interp_poly(ConstVectorView x,
                    ConstVectorView y,
                    const Numeric& x_i,
                    const GridPos& gp) {
  Index N_x = x.nelem();

  ARTS_ASSERT(N_x == y.nelem());
  ARTS_ASSERT(N_x > 2);

  Vector xa(4), ya(4);
  Numeric y_int;
  Numeric dy_int;
  y_int = 0.;

  // 1 - polynomial interpolation (3 points) with grid position search
  // 2 - polynomial interpolation (3 points) without grid position search
  // 3 - polynomial interpolation (4 points)

  Index interp_method = 1;

  if (interp_method == 1) {
    // Pick out three points for interpolation
    if ((gp.fd[0] <= 0.5 && gp.idx > 0) || gp.idx == N_x - 2) {
      xa[0] = x[gp.idx - 1];
      xa[1] = x[gp.idx];
      xa[2] = x[gp.idx + 1];

      ya[0] = y[gp.idx - 1];
      ya[1] = y[gp.idx];
      ya[2] = y[gp.idx + 1];
    }

    else if ((gp.fd[0] > 0.5 && gp.idx < N_x - 2) || gp.idx == 0) {
      xa[0] = x[gp.idx];
      xa[1] = x[gp.idx + 1];
      xa[2] = x[gp.idx + 2];

      ya[0] = y[gp.idx];
      ya[1] = y[gp.idx + 1];
      ya[2] = y[gp.idx + 2];
    }

    else if (gp.idx == N_x - 1) {
      xa[0] = x[N_x - 2];
      xa[1] = x[N_x - 1];
      xa[2] = x[N_x];

      ya[0] = y[N_x - 2];
      ya[1] = y[N_x - 1];
      ya[2] = y[N_x];
    } else {
      ARTS_ASSERT(false);
      arts_exit();
    }

    polint(y_int, dy_int, xa, ya, 3, x_i);

  }

  else if (interp_method == 2) {
    if (gp.idx == 0) {
      xa[0] = x[gp.idx];
      xa[1] = x[gp.idx + 1];
      xa[2] = x[gp.idx + 2];

      ya[0] = y[gp.idx];
      ya[1] = y[gp.idx + 1];
      ya[2] = y[gp.idx + 2];
    } else if (gp.idx == N_x - 1) {
      xa[0] = x[gp.idx - 2];
      xa[1] = x[gp.idx - 1];
      xa[2] = x[gp.idx];

      ya[0] = y[gp.idx - 2];
      ya[1] = y[gp.idx - 1];
      ya[2] = y[gp.idx];
    } else {
      xa[0] = x[gp.idx - 1];
      xa[1] = x[gp.idx];
      xa[2] = x[gp.idx + 1];

      ya[0] = y[gp.idx - 1];
      ya[1] = y[gp.idx];
      ya[2] = y[gp.idx + 1];
    }

    // Polynominal interpolation, n = 3
    polint(y_int, dy_int, xa, ya, 3, x_i);
  }

  else if (interp_method == 3) {
    // Take 4 points
    if (gp.idx == 0) {
      xa[0] = -x[gp.idx + 1];
      xa[1] = x[gp.idx + 0];
      xa[2] = x[gp.idx + 1];
      xa[3] = x[gp.idx + 2];

      ya[0] = y[gp.idx + 1];
      ya[1] = y[gp.idx + 0];
      ya[2] = y[gp.idx + 1];
      ya[3] = y[gp.idx + 2];
    } else if (gp.idx == N_x - 1) {
      xa[0] = x[gp.idx - 1];
      xa[1] = x[gp.idx - 0];
      xa[2] = 2 * x[gp.idx] - x[gp.idx - 1];
      xa[3] = 2 * x[gp.idx] - x[gp.idx - 2];

      ya[0] = y[gp.idx - 1];
      ya[1] = y[gp.idx - 0];
      ya[2] = y[gp.idx - 1];
      ya[3] = y[gp.idx - 2];
    } else if (gp.idx == N_x - 2) {
      xa[0] = x[gp.idx - 2];
      xa[1] = x[gp.idx - 1];
      xa[2] = x[gp.idx];
      xa[3] = x[gp.idx + 1];

      ya[0] = y[gp.idx - 2];
      ya[1] = y[gp.idx - 1];
      ya[2] = y[gp.idx];
      ya[3] = y[gp.idx + 1];
    } else {
      xa[0] = x[gp.idx - 1];
      xa[1] = x[gp.idx];
      xa[2] = x[gp.idx + 1];
      xa[3] = x[gp.idx + 2];

      ya[0] = y[gp.idx - 1];
      ya[1] = y[gp.idx];
      ya[2] = y[gp.idx + 1];
      ya[3] = y[gp.idx + 2];
    }
    // Polinominal interpolation, n = 4
    polint(y_int, dy_int, xa, ya, 4, x_i);
  }

  return y_int;
}
