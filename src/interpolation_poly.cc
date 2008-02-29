/* Copyright (C) 2008 Stefan Buehler <sbuehler(at)ltu.se>

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
  \file   interpolation_poly.cc
  \author Stefan Buehler <sbuehler(at)ltu.se>
  \date   2008-02-04
  
  \brief  Interpolation routines for cubic and higher order interpolation. 
  
  The data structures and functions provided here follow the same
  philosophy as those for linear interpolation in
  interpolation{.h,.cc}. You will need a sequence of three steps to
  perform an interpolation: 

  -# gridpos_poly (one for each interpolation dimension)
  -# interpweights_poly
  -# interp_poly
  
  Not only is the philosophy the same, these higher order functions
  also make direct use of the linear functions in some important
  cases. 
*/

#include <iostream>
#include "interpolation_poly.h"
#include "interpolation.h"
#include "logic.h"

// These two macros here are from Numerical Receipes. They give the
// maxium and minimum of two Index arguments.
static Index imaxarg1, imaxarg2;
#define IMAX(a,b) (imaxarg1=(a), imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
                   (imaxarg1) : (imaxarg2))

static Index iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a), iminarg2=(b),(iminarg1) < (iminarg2) ?\
                   (iminarg1) : (iminarg2))

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
const Numeric sum_check_epsilon = 1e-6;

//! Set up grid positions for higher order interpolation.
/*!
  This function performs the same task as gridpos, but for arbitrary
  orders of interpolation. (Linear is also possible, as a special
  case.) 
  
  The formula for calculating the weights w is taken from Numerical
  Recipes, 2nd edition, section 3.1, eq. 3.1.1.

  \param gp Output: An array of grid positions.
  \param old_grid Original grid.
  \param new_grid New grid.
  \param order Interpolation order.
               1 = linear, 2 = quadratic, etc..
               The number of points used in the
               interpolation will be order+1.
  \param extpolfac Extrapolation fraction. Should normally not be
                   specified, then the default of 0.5 is used. 
*/
void gridpos_poly(ArrayOfGridPosPoly& gp,
                  ConstVectorView old_grid,
                  ConstVectorView new_grid,
                  const Index order,
                  const Numeric  extpolfac)
{
  // Number of points used in the interpolation (order + 1):
  Index m=order+1;

  const Index n_old = old_grid.nelem();
  const Index n_new = new_grid.nelem();

  // Consistently with gridpos, the array size of gp has to be set
  // outside. Here, we only assert that it is correct:
  assert( is_size(gp,n_new) );
  
  // First call the traditional gridpos to find the grid positions:
  ArrayOfGridPos gp_trad(n_new);
  gridpos_extpol( gp_trad, old_grid, new_grid, extpolfac );

  for (Index s=0; s<n_new; ++s)
    {
                   
      // Here we calculate the index of the first of the range of
      // points used for interpolation. For linear interpolation this
      // is identical to j. The idea for this expression is from
      // Numerical Receipes (Chapter 3, section "after the hunt"), but
      // there is is for 1-based arrays.
      Index k = IMIN(IMAX(gp_trad[s].idx-(m-1)/2, 0),
                     n_old-m);

      cout << "m: "<< m << ", gp[s].idx: " << gp[s].idx << ", k: " << k << endl;

      gp[s].idx = k;

      // Make gp[s].w the right size:
      gp[s].w.resize(m);
      
      // Calculate w for each interpolation point. In the linear case
      // these are just the fractional distances to each interpolation
      // point. The w here correspond exactly to the terms in fron of
      // the yi in Numerical Recipes, 2nd edition, section 3.1,
      // eq. 3.1.1.
      for (Index i=0; i<m; ++i)
        {
          //  Numerical Recipes, 2nd edition, section 3.1, eq. 3.1.1.

          // Numerator:
          Numeric num = 1;
          for (Index j=0; j<m; ++j)
            if (j!=i)
              num *= new_grid[s] - old_grid[k+j];
      
          // Denominator:
          Numeric denom = 1;
          for (Index j=0; j<m; ++j)
            if (j!=i)
              denom *= old_grid[k+i] - old_grid[k+j];

          gp[s].w[i] = num / denom;
        }

      // Debugging: Test if sum of all w is 1, as it should be:
      Numeric testsum = 0;
      for (Index i=0; i<m; ++i) testsum += gp[s].w[i];
      cout << "Testsum = " << testsum << endl;        
      
    }
}


////////////////////////////////////////////////////////////////////////////
//                      Red Interpolation
////////////////////////////////////////////////////////////////////////////


//! Red 1D polynomial interpolation weights.
/*!
  This is like the corresponding *interpweights* function, but works
  for interpolation with arbitrary polynomial order. Interpolation
  order 1 should give the same result as our traditional linear
  interpolation routines (although the code is different). I did this
  on purpose, so that consistency can be checked. Order 0 (nearest
  neighbor) is not implemented, but could easily be, if anybody ever
  needs it.

  In contrast to the traditional *interpweights*, the output vector itw
  is resized automatically, so that it is easy to switch interpolation
  orders without changing a lot of code.

  The interpolation order is determined from the size of the w vector
  in the grid positions tc.

  \retval itw Interpolation weights.
  \param  tc Grid position (of the interpolation point in old_grid).
*/
void interpweights_poly( Vector& itw,
                         const GridPosPoly& tc )
{
  // We need the number of interpolation points, which is the
  // interpolation order plus one.
  Index m=tc.w.nelem();
  
  // In the linear case we need 2 weights, in the quadratic case 3, etc..
  itw.resize(m);

  // This loop is over all points used in the interpolation:
  for (Index i=0; i<m; ++i)
    {
      itw[i] = tc.w[i];
    }

}

//! Red 1D Polynomial Interpolate.
/*! 
  "Red" interpolation returns just a scalar.

  The dimension of itw must be consistent with the dimension of the
  interpolation (m^n), where m is the number of points in the
  interpolation scheme (2 for linear), and n is the dimension of the field.

  \param itw  Interpolation weights.
  \param a    The field to interpolate.
  \param tc   The grid position for the column dimension.

  \return Interpolated value.
*/
Numeric interp_poly( ConstVectorView    itw,
                     ConstVectorView    a,    
                     const GridPosPoly& tc )
{
  // Number of points in interpolation scheme:
  Index m = itw.nelem();

  // Dimensions of itw and tc must be consistent:
  assert( is_size(tc.w,m) );

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one.
  assert( is_same_within_epsilon( itw.sum(),
                                  1,
                                  sum_check_epsilon ) );
  
  // To store interpolated value:
  Numeric tia = 0;

  Index iti = 0;
  for ( Index c=0; c<m; ++c )
    {
      tia += a[tc.idx+c] * itw[iti];
      ++iti;
    }

  return tia;
}


////////////////////////////////////////////////////////////////////////////
//                      Blue interpolation
////////////////////////////////////////////////////////////////////////////

//! Compute 1D polynomial interpolation weights.
/*! 
  For this 1D case there is no distinction between "blue" and "green"
  type interpolation.

  This is like the corresponding *interpweights* function, but works
  for interpolation with arbitrary polynomial order. Interpolation
  order 1 should give the same result as our traditional linear
  interpolation routines (although the code is different). I did this
  on purpose, so that consistency can be checked. Order 0 (nearest
  neighbor) is not implemented, but could easily be, if anybody ever
  needs it.

  Note that we still do not need the actual field for this step.

  In contrast to the linear *interpweights* function, itw is sized
  automatically, to allow easy switching between interpolation orders.
  
  The interpolation order is determined from the size of the w vector
  in the grid positions tc.

  \retval itw Interpolation weights.
  \param cgp  The grid position Array for the column dimension.
*/
void interpweights_poly( Matrix& itw,
                         const ArrayOfGridPosPoly& cgp )
{
  Index n = cgp.nelem();

  // We need the number of interpolation points, which is the
  // interpolation order plus one.
  Index m=cgp[0].w.nelem();

  itw.resize(n,m);      // We must store m interpolation
                        // weights for each position.

  // We have to loop all the points in the sequence:
  for ( Index s=0; s<n; ++s )
    {
      // Current grid positions:
      const GridPosPoly& tc = cgp[s];
      
      // Check that the number of interpolation points to use is the
      // same for all points:
      assert(is_size(tc.w,m));
      
      // This loop is over all points used in the interpolation:
      for (Index i=0; i<m; ++i)
        itw(s,i) = tc.w[i];
    }
}

//! Polynomial interpolation of 1D field.
/*! 
  For this 1D case there is no distinction between "blue" and "green"
  type interpolation.

  This is like the corresponding *interp* function, but works
  for interpolation with arbitrary polynomial order. Interpolation
  order 1 should give the same result as our traditional linear
  interpolation routines (although the code is different). 

  The output vector ia must have the same length as the grid position
  vector cgp. And the dimension of itw must be consistent with
  this.

  \retval ia Vector containing the interpolated field values.
  \param itw Interpolation weights.
  \param a   The field to interpolate.
  \param cgp The grid position Array for the column dimension.
*/
void interp_poly( VectorView            ia,
                  ConstMatrixView       itw,
                  ConstVectorView       a,    
                  const ArrayOfGridPosPoly& cgp)
{
  // Number of point in the sequence:
  Index n = itw.nrows();
  // Number of interpolation points (interpolation order plus one):
  Index m=itw.ncols();

  assert(is_size(cgp,n));       // cgp must have one element for each point.
  assert(is_size(ia,n));        // ia must have same size as cgp.

  // Check that interpolation weights are valid. The sum of all
  // weights (last dimension) must always be approximately one. We
  // only check the first element.
  assert( is_same_within_epsilon( itw(0,Range(joker)).sum(),
                                  1,
                                  sum_check_epsilon ) );
  
  // We have to loop all the points in the sequence:
  for ( Index i=0; i<n; ++i )
    {
      // Current grid positions:
      const GridPosPoly& tc = cgp[i];

      // Get handle to current element of output vector and initialize
      // it to zero:
      Numeric& tia = ia[i];
      tia = 0;

      Index iti = 0;
      for ( Index c=0; c<m; ++c )
        {
          tia += a[tc.idx+c] * itw(i,iti);
          ++iti;
        }
    }
}
