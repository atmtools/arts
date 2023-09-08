/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
  \file   check_input.cc
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date 2002-04-15 

  General functions to check the size and logic of input to functions.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "check_input.h"
#include <cfloat>
#include <cmath>
#include <stdexcept>
#include "array.h"
#include <workspace.h>
#include "gridded_fields.h"
#include "logic.h"

using GriddedFieldGrids::GFIELD2_LAT_GRID;
using GriddedFieldGrids::GFIELD2_LON_GRID;
using GriddedFieldGrids::GFIELD3_P_GRID;
using GriddedFieldGrids::GFIELD3_LAT_GRID;
using GriddedFieldGrids::GFIELD3_LON_GRID;

/*===========================================================================
  === Functions for Index
  ===========================================================================*/

//! chk_if_bool
/*! 
    Checks that a variable of type Index has the value 0 or 1.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Index.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_bool(const String& x_name, const Index& x) {
  ARTS_USER_ERROR_IF (!is_bool(x),
      "The variable *", x_name, "* must be a boolean (0 or 1).\n"
      "The present value of *", x_name, "* is ", x, ".")
}

//! chk_if_in_range
/*! 
    Checks that a variable of type Index has a value inside the specified
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Index.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_in_range(const String& x_name,
                     const Index& x,
                     const Index& x_low,
                     const Index& x_high) {
  ARTS_USER_ERROR_IF ((x < x_low) || (x > x_high),
    "The variable *", x_name, "* must fulfill:\n"
      , "   ", x_low, " <= ", x_name, " <= ", x_high, "\n"
      , "The present value of *", x_name, "* is ", x, ".")
}

//! chk_if_increasing
/*! 
    Checks if an ArrayOfIndex is strictly increasing. Cloned from
    Patricks similar function for Vector.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type ArrayOfIndex.

    \author Stefan Buehler
    \date   2007-05-18
*/
void chk_if_increasing(const String& x_name, const ArrayOfIndex& x) {
  ARTS_USER_ERROR_IF (!is_increasing(x),
      "The ArrayOfIndex *", x_name, "* must have strictly\n"
      "increasing values, but this is not the case.\n"
      "x = ", x, "\n")
}

/*===========================================================================
  === Functions for Numeric
  ===========================================================================*/

//! chk_if_positive
/*! 
    Checks that a variable of type Numeric is >0
    The function gives an error message if this is not the case.
    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_positive(const String& x_name, const Numeric& x) {
  ARTS_USER_ERROR_IF (x <= 0,
      "The variable *", x_name, "* must be > 0.\n"
      "The present value of *", x_name, "* is ", x, ".")
}

//! chk_not_negative
/*! 
    Checks that a variable of type Numeric is 0 or positive.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_not_negative(const String& x_name, const Numeric& x) {
  ARTS_USER_ERROR_IF (x < 0,
      "The variable *", x_name, "* must be >= 0.\n"
      "The present value of *", x_name, "* is ", x, ".")
}

//! chk_if_in_range
/*! 
    Checks that a variable of type Numeric has a value inside the specified
    range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_in_range(const String& x_name,
                     const Numeric& x,
                     const Numeric& x_low,
                     const Numeric& x_high) {
  ARTS_USER_ERROR_IF ((x < x_low) || (x > x_high),
      "The variable *", x_name, "* must fulfill:\n"
      "   ", x_low, " <= ", x_name, " <= ", x_high, "\n"
      "The present value of *", x_name, "* is ", x, ".")
}

//! chk_if_in_range_exclude_low
/*! 
    Checks that a variable of type Numeric has a value inside the specified
    range. The low value is excluded from the valid range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Oliver Lemke
    \date   2016-05-10
*/
void chk_if_in_range_exclude_low(const String& x_name,
                                 const Numeric& x,
                                 const Numeric& x_low,
                                 const Numeric& x_high) {
  ARTS_USER_ERROR_IF ((x <= x_low) || (x > x_high),
      "The variable *", x_name, "* must fulfill:\n"
      "   ", x_low, " < ", x_name, " <= ", x_high, "\n"
      "The present value of *", x_name, "* is ", x, ".")
}

//! chk_if_in_range_exclude_high
/*! 
    Checks that a variable of type Numeric has a value inside the specified
    range. The high value is excluded from the valid range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Oliver Lemke
    \date   2016-05-10
*/
void chk_if_in_range_exclude_high(const String& x_name,
                                  const Numeric& x,
                                  const Numeric& x_low,
                                  const Numeric& x_high) {
  ARTS_USER_ERROR_IF ((x < x_low) || (x >= x_high),
      "The variable *", x_name, "* must fulfill:\n"
      "   ", x_low, " <= ", x_name, " < ", x_high, "\n"
      "The present value of *", x_name, "* is ", x, ".")
}

//! chk_if_in_range_exclude
/*! 
    Checks that a variable of type Numeric has a value inside the specified
    range. The low and high values are excluded from the valid range.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Numeric.
    \param    x_low    Lowest allowed value for x.
    \param    x_high   Highest allowed value for x.

    \author Oliver Lemke
    \date   2016-05-10
*/
void chk_if_in_range_exclude(const String& x_name,
                             const Numeric& x,
                             const Numeric& x_low,
                             const Numeric& x_high) {
  ARTS_USER_ERROR_IF ((x <= x_low) || (x >= x_high),
      "The variable *", x_name, "* must fulfill:\n"
      "   ", x_low, " < ", x_name, " < ", x_high, "\n"
      "The present value of *", x_name, "* is ", x, ".")
}

/*===========================================================================
  === Functions for Vector
  ===========================================================================*/

//! chk_vector_length
/*! 
    Checks that a vector has the specified length.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_vector_length(const String& x_name,
                       ConstVectorView x,
                       const Index& l) {
  ARTS_USER_ERROR_IF (x.nelem() != l,
      "The vector *", x_name, "* must have the length ", l, ".\n"
      "The present length of *", x_name, "* is ", x.nelem(), ".")
}

//! chk_vector_length
/*! 
    Checks if two vectors have the same length.

    The function gives an error message if this is not the case.

    \param    x1_name   The name of the first vector
    \param    x2_name   The name of the second vector
    \param    x1        The first vector.
    \param    x2        The second vector.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_vector_length(const String& x1_name,
                       const String& x2_name,
                       ConstVectorView x1,
                       ConstVectorView x2) {
  ARTS_USER_ERROR_IF (x1.nelem() != x2.nelem(),
      "The vectors *", x1_name, "* and *", x2_name,
      "* must have the same length.\n"
      "The length of *", x1_name, "* is ", x1.nelem(), ".\n"
      "The length of *", x2_name, "* is ", x2.nelem(), ".")
}

//! chk_if_increasing
/*! 
    Checks if a vector is strictly increasing.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_increasing(const String& x_name, ConstVectorView x) {
  ARTS_USER_ERROR_IF (!is_increasing(x),
      "The vector *", x_name, "* must have strictly\n"
      "increasing values, but this is not the case.\n"
      "x = ", x, "\n")
}

//! chk_if_decreasing
/*! 
    Checks if a vector is strictly decreasing.

    Duplicated values are not allowed.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A variable of type Vector.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_if_decreasing(const String& x_name, ConstVectorView x) {
  ARTS_USER_ERROR_IF (!is_decreasing(x),
      "The vector *", x_name, "* must have strictly\ndecreasing "
      "values, but this is not the case.\n")
}

//! chk_if_equal
/*!
 * Checks if two vectors are equal within a margin.
 *
 * \param   x1_name The name of the first variable (used in error message).
 * \param   x2_name The name of the second variable (used in error message).
 * \param   v1      First vector
 * \param   v2      Second vector
 * \param   margin  uncertainty margin. Default: 1e-6
 *
 * \author  Gerrit Holl
 * \date    2011-05-04
*/

void chk_if_equal(const String& x1_name,
                  const String& x2_name,
                  ConstVectorView v1,
                  ConstVectorView v2,
                  Numeric margin) {
  chk_vector_length(x1_name, x2_name, v1, v2);

  for (Index i = 0; i < v1.nelem(); i++) {
    ARTS_USER_ERROR_IF (abs(v1[i] - v2[i]) > margin,
        "Vectors ", x1_name, " and ", x2_name, " differ.\n",
        x1_name, "[", i, "]"
        , " = ", v1[i], "\n"
        , x2_name, "[", i, "]"
        , " = ", v2[i], "\n"
        , "Difference should not exceed ", margin, "\n")
  }
}

/*===========================================================================
  === Functions for Matrix
  ===========================================================================*/

//! chk_matrix_ncols
/*! 
    Checks that a matrix has the specified number of columns.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A matrix.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2002-05-16
*/
void chk_matrix_ncols(const String& x_name, ConstMatrixView x, const Index& l) {
  ARTS_USER_ERROR_IF (x.ncols() != l,
    "The matrix *", x_name, "* must have ", l, " columns,\n"
      , "but the number of columns is ", x.ncols(), ".")
}

//! chk_matrix_nrows
/*! 
    Checks that a matrix has the specified number of rows.

    The function gives an error message if this is not the case.

    \param    x_name   The name of the variable.
    \param    x        A matrix.
    \param    l        The expected length of x.

    \author Patrick Eriksson 
    \date   2002-05-16
*/
void chk_matrix_nrows(const String& x_name, ConstMatrixView x, const Index& l) {
  ARTS_USER_ERROR_IF (x.nrows() != l,
      "The matrix *", x_name, "* must have ", l, " rows,\n"
      , "but the number of rows is ", x.nrows(), ".")
}

/*===========================================================================
  === Functions for Tensors
  ===========================================================================*/

//! Runtime check for size of Vector.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Vector.
  \param    c        Required number of elements

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name, ConstVectorView x, const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, c),
    "The object *", x_name, "* does not have the right size.\n"
    "Dimension should be:" " ", c,  ",\nbut it is:          " " ",
    x.nelem(), ".")
}

//! Runtime check for size of Matrix.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Matrix.
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name,
              ConstMatrixView x,
              const Index& r,
              const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, r, c),
      "The object *", x_name, "* does not have the right size.\n"
      , "Dimensions should be:"
      , " ", r, " ", c, ",\nbut they are:         "
      , " ", x.nrows(), " ", x.ncols(), ".")
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Tensor3.
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name,
              ConstTensor3View x,
              const Index& p,
              const Index& r,
              const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, p, r, c),
      "The object *", x_name, "* does not have the right size.\n"
      , "Dimensions should be:"
      , " ", p, " ", r, " ", c, ",\nbut they are:         "
      , " ", x.npages(), " ", x.nrows(), " ", x.ncols(), ".")
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Tensor4.
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name,
              ConstTensor4View x,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, b, p, r, c),
      "The object *", x_name, "* does not have the right size.\n"
      , "Dimensions should be:"
      , " ", b, " ", p, " ", r, " ", c
      , ",\nbut they are:         "
      , " ", x.nbooks(), " ", x.npages(), " ", x.nrows(), " "
      , x.ncols(), ".")
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Tensor5.
  \param    s        Required number of shelves
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name,
              ConstTensor5View x,
              const Index& s,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, s, b, p, r, c),
      "The object *", x_name, "* does not have the right size.\n"
      , "Dimensions should be:"
      , " ", s, " ", b, " ", p, " ", r, " ", c
      , ",\nbut they are:         "
      , " ", x.nshelves(), " ", x.nbooks(), " ", x.npages(), " "
      , x.nrows(), " ", x.ncols(), ".")
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Tensor6.
  \param    v        Required number of vitrines
  \param    s        Required number of shelves
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name,
              ConstTensor6View x,
              const Index& v,
              const Index& s,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, v, s, b, p, r, c),
      "The object *", x_name, "* does not have the right size.\n"
      , "Dimensions should be:"
      , " ", v, " ", s, " ", b, " ", p, " ", r, " ", c
      , ",\nbut they are:         "
      , " ", x.nvitrines(), " ", x.nshelves(), " ", x.nbooks()
      , " ", x.npages(), " ", x.nrows(), " ", x.ncols(), ".")
}

//! Runtime check for size of Tensor.
/*! 
  This is the runtime version of is_size. An appropriate error message
  is generated if the size is not correct.

  \param    x_name   The name of the agenda.
  \param    x        A variable of type Tensor7.
  \param    l        Required number of libraries
  \param    v        Required number of vitrines
  \param    s        Required number of shelves
  \param    b        Required number of books
  \param    p        Required number of pages
  \param    r        Required number of rows
  \param    c        Required number of columns

  \author Stefan Buehler
  \date   2002-11-29
*/
void chk_size(const String& x_name,
              ConstTensor7View x,
              const Index& l,
              const Index& v,
              const Index& s,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c) {
  ARTS_USER_ERROR_IF (!is_size(x, l, v, s, b, p, r, c),
      "The object *", x_name, "* does not have the right size.\n"
      , "Dimensions should be:"
      , " ", l, " ", v, " ", s, " ", b, " ", p, " ", r
      , " ", c, ",\nbut they are:         "
      , " ", x.nlibraries(), " ", x.nvitrines(), " ", x.nshelves()
      , " ", x.nbooks(), " ", x.npages(), " ", x.nrows(), " "
      , x.ncols(), ".")
}

/*===========================================================================
  === Functions for interpolation grids
  ===========================================================================*/

//! Check interpolation grids
/*!
 This function checks if old and new grid for an interpolation are
 ok. If not, it throws a detailed runtime error message. This is
 intended for workspace method input variable checking. If the original grid does
 not have to cover the whole new grid. The returned ing_min and ing_max give
 the positions in the new grid of the first values that are outside the old grid.
 This is only allowed if the boundary value in the input data is 0.

 \param[out] ing_min             Index in the new grid with first value covered
                                 by the old grid.
 \param[out] ing_max             Index in the new grid with last value covered
                                 by the old grid.
 \param[in]  which_interpolation A string describing the interpolation for
                                 which the grids are intended.
 \param[in]  old_grid            The original grid.
 \param[in]  new_grid            The new grid.
 \param[in]  data                The data for the interpolation.
 \param[in]  order               Interpolation order. (Default value is 1.)
 \param[in]  extpolfac           The extrapolation fraction. See gridpos function
                                 for details. Has a default value, which is
                                 consistent with gridpos.

 \author Oliver Lemke
 \date   2012-07-11
 */
void chk_interpolation_grids_loose(Index& ing_min,
                                   Index& ing_max,
                                   const String& which_interpolation,
                                   ConstVectorView old_grid,
                                   ConstVectorView new_grid,
                                   ConstVectorView data,
                                   const Index order) {
  chk_interpolation_grids_loose_no_data_check(
      ing_min, ing_max, which_interpolation, old_grid, new_grid, order);

  chk_interpolation_grids_loose_check_data(
      ing_min, ing_max, which_interpolation, old_grid, new_grid, data);
}

//! Check interpolation grids
/*!
 This function checks if old and new grid for an interpolation are
 ok. If not, it throws a detailed runtime error message. This is
 intended for workspace method input variable checking. The original grid does
 not have to cover the whole new grid. The returned ing_min and ing_max give
 the positions in the new grid of the first values that are outside the old grid.
 
 \param[out] ing_min             Index in the new grid with first value covered
                                 by the old grid.
 \param[out] ing_max             Index in the new grid with last value covered
                                 by the old grid.
 \param[in]  which_interpolation A string describing the interpolation for
                                 which the grids are intended. 
 \param[in]  old_grid            The original grid.
 \param[in]  new_grid            The new grid.
 \param[in]  order               Interpolation order. (Default value is 1.)
 \param[in]  extpolfac           The extrapolation fraction. See gridpos function
                                 for details. Has a default value, which is
                                 consistent with gridpos.  
 
 \author Oliver Lemke (based on chk_interpolation_grids by Stefan)
 \date   2012-03-28
 */
void chk_interpolation_grids_loose_no_data_check(
    Index& ing_min,
    Index& ing_max,
    const String& which_interpolation,
    ConstVectorView old_grid,
    ConstVectorView new_grid,
    const Index order) {
  const Index n_old = old_grid.nelem();

  ARTS_USER_ERROR_IF (!new_grid.nelem(),
                      "The new grid is not allowed to be empty.");

  // Old grid must have at least order+1 elements:
  ARTS_USER_ERROR_IF (n_old < order + 1, "There is a problem with the grids for the following interpolation:\n",
                      which_interpolation, "\n",
    "The original grid must have at least ", order + 1, " elements.")

  // Decide whether we have an ascending or descending grid:
  const bool ascending = (old_grid[0] <= old_grid[1]);

  // Minimum and maximum allowed value from old grid. (Will include
  // extrapolation tolerance.)
  Numeric og_min, og_max;

  ing_min = 0;
  ing_max = new_grid.nelem() - 1;
  if (ascending) {
    // Old grid must be strictly increasing (no duplicate values.)
    ARTS_USER_ERROR_IF (!is_increasing(old_grid), "There is a problem with the grids for the following interpolation:\n",
                        which_interpolation, "\n",
      "The original grid must be strictly sorted\n"
        , "(no duplicate values). Yours is:\n"
        , old_grid, ".")

    // Limits of extrapolation.
    og_min = old_grid[0];
    og_max = old_grid[n_old - 1];
  } else {
    // Old grid must be strictly decreasing (no duplicate values.)
    ARTS_USER_ERROR_IF (!is_decreasing(old_grid), "There is a problem with the grids for the following interpolation:\n",
                        which_interpolation, "\n",
      "The original grid must be strictly sorted\n"
        , "(no duplicate values). Yours is:\n"
        , old_grid, ".")

    // The max is now the first point, the min the last point!
    og_max = old_grid[0];
    og_min = old_grid[n_old - 1];
  }

  // Min and max of new grid:
  const Numeric ng_min = min(new_grid);
  const Numeric ng_max = max(new_grid);

  // If new grid is not inside old grid, determine the indexes of the range
  // that is.

  const Index iog_min = 0;
  const Index iog_max = old_grid.nelem() - 1;

  ing_min = 0;
  ing_max = new_grid.nelem() - 1;

  if (ascending) {
    if (ng_max > og_max) {
      while (ing_max > 0 && new_grid[ing_max] > old_grid[iog_max]) ing_max--;
    }

    if (ng_min < og_min) {
      while (ing_min < new_grid.nelem() - 1 &&
             new_grid[ing_min] < old_grid[iog_min])
        ing_min++;
    }
  } else {
    if (ng_min < og_min) {
      while (ing_max > 0 && new_grid[ing_max] < old_grid[iog_max]) ing_max--;
    }

    if (ng_max > og_max) {
      while (ing_min < new_grid.nelem() - 1 &&
             new_grid[ing_min] > old_grid[iog_min])
        ing_min++;
    }
  }
}

//! Check log pressure interpolation grids
/*!
 This function checks if old and new grid for an interpolation are
 ok. If not, it throws a detailed runtime error message. This is
 intended for workspace method input variable checking. The original grid does
 not have to cover the whole new grid. The returned ing_min and ing_max give
 the positions in the new grid of the first values that are outside the old grid.

 \param[out] ing_min             Index in the new grid with first value covered
                                 by the old grid.
 \param[out] ing_max             Index in the new grid with last value covered
                                 by the old grid.
 \param[in]  which_interpolation A string describing the interpolation for
                                 which the grids are intended.
 \param[in]  old_grid            The original grid.
 \param[in]  new_grid            The new grid.
 \param[in]  order               Interpolation order. (Default value is 1.)
 \param[in]  extpolfac           The extrapolation fraction. See gridpos function
                                 for details. Has a default value, which is
                                 consistent with gridpos.

 \author Oliver Lemke (based on chk_interpolation_grids by Stefan)
 \date   2012-03-28
 */
void chk_interpolation_pgrids_loose_no_data_check(
    Index& ing_min,
    Index& ing_max,
    const String& which_interpolation,
    ConstVectorView old_pgrid,
    ConstVectorView new_pgrid,
    const Index order) {
  // Local variable to store log of the pressure grids
  Vector logold(old_pgrid.nelem());
  Vector lognew(new_pgrid.nelem());

  transform(logold, log, old_pgrid);
  transform(lognew, log, new_pgrid);

  chk_interpolation_grids_loose_no_data_check(
      ing_min, ing_max, which_interpolation, logold, lognew, order);
}

//! Check interpolation grids
/*!
 This function checks if old and new grid for an interpolation are
 ok. If not, it throws a detailed runtime error message. This is
 intended for workspace method input variable checking. If the original grid does
 not have to cover the whole new grid. The returned ing_min and ing_max give
 the positions in the new grid of the first values that are outside the old grid.
 This is only allowed if the boundary value in the input data is 0.

 \param[out] ing_min             Index in the new grid with first value covered
                                 by the old grid.
 \param[out] ing_max             Index in the new grid with last value covered
                                 by the old grid.
 \param[in]  which_interpolation A string describing the interpolation for
                                 which the grids are intended.
 \param[in]  old_grid            The original grid.
 \param[in]  new_grid            The new grid.
 \param[in]  data                The data for the interpolation.

 \author Oliver Lemke
 \date   2012-03-28
 */
void chk_interpolation_grids_loose_check_data(Index& ing_min,
                                              Index& ing_max,
                                              const String& which_interpolation,
                                              ConstVectorView old_grid,
                                              ConstVectorView new_grid,
                                              ConstVectorView data) {
  ARTS_USER_ERROR_IF (!new_grid.nelem(),
                      "The new grid is not allowed to be empty.");

  // Decide whether we have an ascending or descending grid:
  const bool ascending = (old_grid[0] <= old_grid[1]);

  // If new grid is not inside old grid, determine the indexes of the range
  // that is.

  const Index iog_min = ascending ? old_grid.nelem() - 1 : 0;
  const Index iog_max = ascending ? 0 : old_grid.nelem() - 1;

  ARTS_USER_ERROR_IF (ing_min > 0 && data[iog_min] != 0, "There is a problem with the grids for the following interpolation:\n"
  , which_interpolation, "\n",
    "\nThe new grid is not fully inside the original grid.\n"
      , "This is allowed if the corresponding boundary value of raw data is 0.\n"
      , "New grid point: ", new_grid[ing_min], "\n"
      , "Old grid point: ", old_grid[iog_min], "\n"
      , "Boundary value: ", data[iog_min])

  ARTS_USER_ERROR_IF (ing_max < new_grid.nelem() - 1 && data[iog_max] != 0, "There is a problem with the grids for the following interpolation:\n"
  , which_interpolation, "\n",
    "\nThe the new grid is not fully inside the original grid.\n"
      , "This is allowed if the corresponding boundary value of raw data is 0.\n"
      , "New grid point: ", new_grid[ing_max], "\n"
      , "Old grid point: ", old_grid[iog_max], "\n"
      , "Boundary value: ", data[iog_max])
}

//! Check interpolation grids
/*!
  This function checks if old and new grid for an interpolation are
  ok. If not, it throws a detailed runtime error message. This is
  intended for workspace method input variable checking. 
  
  \param[in] which_interpolation A string describing the interpolation for
                                 which the grids are intended. 
  \param[in] old_grid            The original grid.
  \param[in] new_grid            The new grid.
  \param[in] order               Interpolation order. (Default value is 1.)
  \param[in] extpolfac           The extrapolation fraction. See gridpos function
                                 for details. Has a default value, which is
                                 consistent with gridpos.  
  
  \author Stefan Buehler
  \date   2008-11-24 
*/
void chk_interpolation_grids(const String& which_interpolation,
                             ConstVectorView old_grid,
                             ConstVectorView new_grid,
                             const Index order,
                             const Numeric& extpolfac,
                             const bool islog) {
  const Index n_old = old_grid.nelem();

  ARTS_USER_ERROR_IF (!new_grid.nelem(),
                      "The new grid is not allowed to be empty.");

  ARTS_USER_ERROR_IF (order < 0,
      "There is a problem with the grids for the following "
      "interpolation:\n",
      which_interpolation, "\n"
      "Interpolation order must be 0 or larger (but your's is ", order,
      ").")

  // Old grid must have at least order+1 elements:
  ARTS_USER_ERROR_IF (n_old < order + 1,
    "There is a problem with the grids for the following "
      , "interpolation:\n"
      , which_interpolation, "\n"
      , "For interpolation order ", order
      , ", the original grid must have at least\n"
      , order + 1, " elements (but your's has only ", n_old, ").")

  // Decide whether we have an ascending or descending grid:
  const bool ascending = ((n_old > 1) ? (old_grid[0] <= old_grid[1]) : true);

  // Minimum and maximum allowed value from old grid. (Will include
  // extrapolation tolerance.)
  Numeric og_min = old_grid[0], og_max = old_grid[0];

  if (ascending) {
    // Old grid must be strictly increasing (no duplicate values.)
    ARTS_USER_ERROR_IF (!is_increasing(old_grid),
      "There is a problem with the grids for the "
        , "following interpolation:\n"
        , which_interpolation, "\n"
        , "The original grid must be strictly sorted\n"
        , "(no duplicate values). Yours is:\n"
        , old_grid, ".")

    // Limits of extrapolation.
    if (n_old > 1) {
      og_min = old_grid[0] - extpolfac * (old_grid[1] - old_grid[0]);
      og_max = old_grid[n_old - 1] +
               extpolfac * (old_grid[n_old - 1] - old_grid[n_old - 2]);
    }
  } else {
    // Old grid must be strictly decreasing (no duplicate values.)
    ARTS_USER_ERROR_IF (!is_decreasing(old_grid),
      "There is a problem with the grids for the "
        , "following interpolation:\n"
        , which_interpolation, "\n"
        , "The original grid must be strictly sorted\n"
        , "(no duplicate values). Yours is:\n"
        , old_grid, ".")

    // The max is now the first point, the min the last point!
    // I think the sign is right here...
    if (n_old > 1) {
      og_max = old_grid[0] - extpolfac * (old_grid[1] - old_grid[0]);
      og_min = old_grid[n_old - 1] +
               extpolfac * (old_grid[n_old - 1] - old_grid[n_old - 2]);
    }
  }

  // Min and max of new grid:
  const Numeric ng_min = min(new_grid);
  const Numeric ng_max = max(new_grid);

  // New grid must be inside old grid (plus extpolfac).
  // (Values right on the edge (ng_min==og_min) are still allowed.)

  if (n_old > 1) {
    ARTS_USER_ERROR_IF (ng_min < og_min,
      "There is a problem with the grids for the "
        , "following interpolation:\n"
        , which_interpolation, "\n"
        , "The minimum of the new grid must be inside "
        , "the original grid.\n(We allow a bit of extrapolation, "
        , "but not so much).\n"
        , "Minimum of original grid:           ", min(old_grid),
        islog ? var_string(" (", exp(min(old_grid)), ")") : var_string(),
        "\nMinimum allowed value for new grid: ", og_min,
        islog ? var_string(" (", exp(og_min), ")") : var_string(),
        "\nActual minimum of new grid:         ", ng_min,
        islog ? var_string(" (", exp(ng_min), ")") : var_string())

    ARTS_USER_ERROR_IF (ng_max > og_max,
      "There is a problem with the grids for the "
        , "following interpolation:\n"
        , which_interpolation, "\n"
        , "The maximum of the new grid must be inside\n"
        , "the original grid. (We allow a bit of extrapolation,\n"
        , "but not so much).\n"
        , "Maximum of original grid:           ", max(old_grid),
        islog ? var_string(" (", exp(max(old_grid)), ")") : var_string(),
        "\nMaximum allowed value for new grid: ", og_max,
        islog ? var_string(" (", exp(og_max), ")") : var_string(),
        "\nActual maximum of new grid:         ", ng_max,
        islog ? var_string(" (", exp(ng_max), ")") : var_string())
  }

  // If we get here, than everything should be fine.
}

//! Check interpolation grids
/*!
  This function checks if old and new grid for an interpolation are
  ok. If not, it throws a detailed runtime error message. This is
  intended for workspace method input variable checking. 
  
  This is for the special case that the new grid is just a single
  Numeric, instead of a Vector. ("Red" interpolation.)
  It just calles the other more general chk_interpolation_grids
  function for which both grid arguments are vectors.

  \param[in] which_interpolation A string describing the interpolation for
                                 which the grids are intended. 
  \param[in] old_grid            The original grid.
  \param[in] new_grid            The new grid.
  \param[in] order               Interpolation order. (Default value is 1.)
  \param[in] extpolfac           The extrapolation fraction. See gridpos function
                                 for details. Has a default value, which is
                                 consistent with gridpos.  
  
  \author Stefan Buehler
  \date   2008-11-24 
*/
void chk_interpolation_grids(const String& which_interpolation,
                             ConstVectorView old_grid,
                             const Numeric& new_grid,
                             const Index order,
                             const Numeric& extpolfac) {
  const Vector v(1, new_grid);
  chk_interpolation_grids(which_interpolation, old_grid, v, order, extpolfac);
}

//! Check log pressure interpolation grids
/*!
 This function checks if old and new grid for an interpolation are
 ok. If not, it throws a detailed runtime error message. This is
 intended for workspace method input variable checking.

 \param[in] which_interpolation A string describing the interpolation for
 which the grids are intended.
 \param[in] old_pgrid            The original grid.
 \param[in] new_pgrid            The new grid.
 \param[in] order               Interpolation order. (Default value is 1.)
 \param[in] extpolfac           The extrapolation fraction. See gridpos function
 for details. Has a default value, which is
 consistent with gridpos.

 \author Oliver Lemke
 \date   2012-07-11
 */
void chk_interpolation_pgrids(const String& which_interpolation,
                              ConstVectorView old_pgrid,
                              ConstVectorView new_pgrid,
                              const Index order,
                              const Numeric& extpolfac) {
  // Local variable to store log of the pressure grids
  Vector logold(old_pgrid.nelem());
  Vector lognew(new_pgrid.nelem());

  transform(logold, log, old_pgrid);
  transform(lognew, log, new_pgrid);

  chk_interpolation_grids(
      which_interpolation, logold, lognew, order, extpolfac, true);
}

/*===========================================================================
  === Functions related to atmospheric and surface grids and fields.
  ===========================================================================*/

//! chk_atm_grids
/*! 
    Checks if the atmospheric grids and the specified atmospheric 
    dimensionality match, and if the grids are ordered correctly.

    The function gives an error message if this is not the case.

    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_atm_grids(const Index& dim,
                   ConstVectorView p_grid,
                   ConstVectorView lat_grid,
                   ConstVectorView lon_grid) {
  // p_grid
  ARTS_USER_ERROR_IF (p_grid.nelem() < 2,
                      "The length of *p_grid* must be >= 2.");
  chk_if_decreasing("p_grid", p_grid);

  // lat_grid
  if (dim == 1) {
    ARTS_USER_ERROR_IF (lat_grid.nelem() > 0,
                        "For dim=1, the length of *lat_grid* must be 0.");
  } else {
    ARTS_USER_ERROR_IF (lat_grid.nelem() < 2,
                        "For dim>1, the length of *lat_grid* must be >= 2.");
    chk_if_increasing("lat_grid", lat_grid);
  }

  // lon_grid
  if (dim < 3) {
    ARTS_USER_ERROR_IF (lon_grid.nelem() > 0,
                        "For dim<3, the length of *lon_grid* must be 0.");
  } else {
    ARTS_USER_ERROR_IF (lon_grid.nelem() < 2,
                        "For dim=3, the length of *lon_grid* must be >= 2.");
    chk_if_increasing("lon_grid", lon_grid);
  }

  // Check that latitude and longitude grids are inside OK ranges for 3D
  if (dim == 3) {
    ARTS_USER_ERROR_IF (lat_grid[0] < -90,
          "The latitude grid cannot extend below -90 degrees for 3D");
    ARTS_USER_ERROR_IF (lat_grid[lat_grid.nelem() - 1] > 90,
          "The latitude grid cannot extend above 90 degrees for 3D");
    ARTS_USER_ERROR_IF (lon_grid[0] < -360,
          "No longitude (in lon_grid) can be below -360 degrees.");
    ARTS_USER_ERROR_IF (lon_grid[lon_grid.nelem() - 1] > 360,
          "No longitude (in lon_grid) can be above 360 degrees.");
    ARTS_USER_ERROR_IF (lon_grid[lon_grid.nelem() - 1] - lon_grid[0] > 360,
          "The longitude grid is not allowed to cover more than 360 degrees.");
  }
}

//! chk_atm_field (simple fields)
/*! 
    Checks if an atmospheric field matches the dimensionality and the grids.

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric field.
    \param    x            A variable holding an atmospheric field.
    \param    dim          The atmospheric dimensionality.
    \param    p_grid       The pressure grid.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.
    \param    chk_lat90    Flag whether pole consistency check to be done (only
                           relevant for dim==3.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_atm_field(const String& x_name,
                   ConstTensor3View x,
                   const Index& dim,
                   ConstVectorView p_grid,
                   ConstVectorView lat_grid,
                   ConstVectorView lon_grid,
                   const bool& chk_lat90) {
  // It is assumed that grids OK-ed through chk_atm_grids
  Index npages = p_grid.nelem(), nrows = 1, ncols = 1;
  if (dim > 1) nrows = lat_grid.nelem();
  if (dim > 2) ncols = lon_grid.nelem();
  ARTS_USER_ERROR_IF (x.ncols() != ncols || x.nrows() != nrows || x.npages() != npages,
    "The atmospheric field *", x_name, "* has wrong size.\n"
      , "Expected size is ", npages, " x ", nrows, " x ", ncols
      , ", while actual size is ", x.npages(), " x ", x.nrows(), " x "
      , x.ncols(), ".")

  // NaNs are not allowed
  for (Index ip = 0; ip < npages; ip++) {
    for (Index ir = 0; ir < nrows; ir++) {
      for (Index ic = 0; ic < ncols; ic++) {
        ARTS_USER_ERROR_IF (std::isnan(x(ip, ir, ic)),
          "The variable *", x_name, "* contains one or "
            , "several NaNs. This is not allowed!")
      }
    }
  }

  // Special 3D checks:
  if (dim == 3) {
    // If all lons are covered, check if cyclic
    if (is_lon_cyclic(lon_grid)) {
      const Index ic = ncols - 1;
      for (Index ip = 0; ip < npages; ip++) {
        for (Index ir = 0; ir < nrows; ir++) {
          ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                  x(ip, ir, ic), x(ip, ir, 0), 4 * DBL_EPSILON),
            "The variable *", x_name, "* covers 360 "
              , "degrees in the longitude direction, but the field "
              , "seems to deviate between first and last longitude "
              , "point. The field must be \"cyclic\".\n"
              , "Difference: ", std::setprecision(16)
              , x(ip, ir, ic) - x(ip, ir, 0), "\n"
              , "Epsilon   : "
              , 4 * DBL_EPSILON * std::max(x(ip, ir, ic), x(ip, ir, 0)))
        }
      }
    }

    chk_if_bool("chk_lat90", chk_lat90);
    if (chk_lat90) {
      // No variation at the South pole!
      if (lat_grid[0] == -90) {
        for (Index ip = 0; ip < npages; ip++) {
          for (Index ic = 1; ic < ncols; ic++) {
            ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                    x(ip, 0, ic), x(ip, 0, ic - 1), 2 * DBL_EPSILON),
              "The variable *", x_name, "* covers the South\n"
                , "pole. The data corresponding to the pole can not\n"
                , "vary with longitude, but this appears to be the\n"
                , "case.")
              /*                            , "case: at ", ip, ".th  pressure it has val\n"
                            , x(ip,0,ic-1), ", but ", x(ip,0,ic)
                            , " (i.e., a diff of ", fabs(x(ip,0,ic)-x(ip,0,ic-1))
                            , ") at ", ic-1, "th and ", ic, "th longitudes!\n";
*/
          }
        }
      }
      // No variation at the North pole!
      if (lat_grid[nrows - 1] == 90) {
        const Index ir = nrows - 1;
        for (Index ip = 0; ip < npages; ip++) {
          for (Index ic = 1; ic < ncols; ic++) {
            ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                    x(ip, ir, ic), x(ip, ir, ic - 1), 2 * DBL_EPSILON),
              "The variable *", x_name, "* covers the North\n"
                , "pole. The data corresponding to the pole can not\n"
                , "vary with longitude, but this appears to be the "
                , "case.")
              /*                            , "case: at ", ip, ".th  pressure it has val\n"
                            , x(ip,ir,ic-1), ", but ", x(ip,ir,ic)
                            , " (i.e., a diff of ", fabs(x(ip,ir,ic)-x(ip,ir,ic-1))
                            , ") at ", ic-1, "th and ", ic, "th longitudes!\n";
*/
          }
        }
      }
    }
  }
}

//! chk_atm_field (fields with one more dimension)
/*! 
    Checks if an atmospheric field matches the dimensionality and the
    grids. This is the version for fields like vmr_field, which are a
    Tensor4, not a Tensor3. (First dimension is the gas species.)

    The function gives an error message if this is not the case.

    \param    x_name       The name of the atmospheric field.
    \param    x            A variable holding an atmospheric field.
    \param    dim          The atmospheric dimensionality.
    \param    nspecies     Number of species.
    \param    p_grid       The pressure grid.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Stefan Buehler, cloned from Patrick Eriksson 
    \date   2002-12-20
*/
void chk_atm_field(const String& x_name,
                   ConstTensor4View x,
                   const Index& dim,
                   const Index& nspecies,
                   ConstVectorView p_grid,
                   ConstVectorView lat_grid,
                   ConstVectorView lon_grid,
                   const bool& check_nan) {
  const Index nbooks = nspecies;
  //
  if (nbooks == 0) {
    ARTS_USER_ERROR_IF (x.nbooks(),
      "The atmospheric field *", x_name, "* should be empty.\n")
    return;
  }

  Index npages = p_grid.nelem(), nrows = 1, ncols = 1;
  if (dim > 1) nrows = lat_grid.nelem();
  if (dim > 2) ncols = lon_grid.nelem();

  ARTS_USER_ERROR_IF (x.ncols() != ncols || x.nrows() != nrows || x.npages() != npages ||
      x.nbooks() != nbooks,
    "The atmospheric field *", x_name, "* has wrong size.\n"
      , "Expected size is ", nbooks, " x ", npages, " x ", nrows
      , " x ", ncols, ",\n"
      , "while actual size is ", x.nbooks(), " x ", x.npages(), " x "
      , x.nrows(), " x ", x.ncols(), ".")

  if (check_nan)
  // NaNs are not allowed
  {
    for (Index ib = 0; ib < nbooks; ib++) {
      for (Index ip = 0; ip < npages; ip++) {
        for (Index ir = 0; ir < nrows; ir++) {
          for (Index ic = 0; ic < ncols; ic++) {
            ARTS_USER_ERROR_IF (std::isnan(x(ib, ip, ir, ic)),
              "The variable *", x_name, "* contains one or "
                , "several NaNs. This is not allowed!")
          }
        }
      }
    }
  }

  // Special 3D checks:
  if (dim == 3) {
    // If all lons are covered, check if cyclic
    if ((lon_grid[ncols - 1] - lon_grid[0]) == 360) {
      const Index ic = ncols - 1;
      for (Index is = 0; is < nspecies; is++) {
        for (Index ip = 0; ip < npages; ip++) {
          for (Index ir = 0; ir < nrows; ir++) {
            ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                    x(is, ip, ir, ic), x(is, ip, ir, 0), 2 * DBL_EPSILON),
              "The variable *", x_name, "* covers 360 "
                , "degrees in the longitude direction, but at least "
                , "one field seems to deviate between first and last "
                , "longitude point. The field must be \"cyclic\". "
                , "This was found for field with index ", is
                , " (0-based).")
          }
        }
      }
    }
    // No variation at the South pole!
    if (lat_grid[0] == -90) {
      for (Index is = 0; is < nspecies; is++) {
        for (Index ip = 0; ip < npages; ip++) {
          for (Index ic = 1; ic < ncols; ic++) {
            ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                    x(is, ip, 0, ic), x(is, ip, 0, ic - 1), 2 * DBL_EPSILON),
              "The variable *", x_name, "* covers the South "
                , "pole. The data corresponding to the pole can not "
                , "vary with longitude, but this appears to be the "
                , "case. This was found for field with index ", is
                , " (0-based).")
          }
        }
      }
    }
    // No variation at the North pole!
    if (lat_grid[nrows - 1] == 90) {
      const Index ir = nrows - 1;
      for (Index is = 0; is < nspecies; is++) {
        for (Index ip = 0; ip < npages; ip++) {
          for (Index ic = 1; ic < ncols; ic++) {
            ARTS_USER_ERROR_IF (!is_same_within_epsilon(x(is, ip, ir, ic),
                                        x(is, ip, ir, ic - 1),
                                        2 * DBL_EPSILON),
              "The variable *", x_name, "* covers the North "
                , "pole. The data corresponding to the pole can not "
                , "vary with longitude, but this appears to be the "
                , "case. This was found for field with index ", is
                , " (0-based).")
          }
        }
      }
    }
  }
}

//! chk_atm_vecfield_lat90
/*! 
    Checks if a two-compnent vector atmospheric field is consistant at the poles.

    Similar to the field-constant at poles check of chk_atm_field, but checking
    the total vector instead of each component to be constant (since this
    involves some numerics, we allow a deviation threshold instead of perfect
    match).
    Note that each component of the vector is stored in a separate atmospheric
    field. Intended for variables that are supposed to be two horizontal components of
    a 3D vector field (e.g., winds, magnetic field).
    It is assumed that individual fields have passed chk_atm_field.
    The function gives an error message if a mismatch is encountered.

    \param    x1_name      The name of the atmospheric field.
    \param    x1           A variable holding an atmospheric field.
    \param    x2_name      The name of the atmospheric field.
    \param    x2           A variable holding an atmospheric field.
    \param    dim          The atmospheric dimensionality.
    \param    lat_grid     The latitude grid.
    \param    threshold    The percentage threshold the total vector lengths
                           along the pole are allowed to deviate.

    \author Jana Mendrok
    \date   2012-06-29
*/
void chk_atm_vecfield_lat90(const String& x1_name,
                            ConstTensor3View x1,
                            const String& x2_name,
                            ConstTensor3View x2,
                            const Index& dim,
                            ConstVectorView lat_grid,
                            const Numeric& threshold) {
  // It is assumed that grids OK-ed through chk_atm_grids and fields
  // individually OK-ed.

  // We only need to check 3D cases. Else there is no variation in longitude
  // anyways.
  if (dim == 3) {
    Index npages = x1.npages();
    Index nrows = x1.nrows();
    Index ncols = x1.ncols();

    // For safety check that both fields have identical dimensions
    ARTS_USER_ERROR_IF (x2.ncols() != ncols || x2.nrows() != nrows || x2.npages() != npages,
      "The atmospheric fields *", x1_name, "* and *", x2_name
        , "* do not match in size.\n"
        , "*", x1_name, "*'s size is ", npages, " x ", nrows
        , " x ", ncols, ", while *", x1_name, "*'s size is "
        , x2.npages(), " x ", x2.nrows(), " x ", x2.ncols(), ".")

    // redefine ratio deviation threshold of vector lengths to ratio of
    // squared vector lengths, cause don't want to calc squareroot everytime.
    // (val1**2/val2**2 - 1) / (val1/val2 - 1) = val1/val2 + 1
    // and with val1~val2                      = 2
    // i.e., (val1**2/val2**2 - 1) ~ 2 * (val1/val2 - 1)
    //
    // with val1/1 = sqrt(vec1/2), hence val1/2**2 = vec1/2
    //       (vec1/vec2 - 1) ~ 2 * (sqrt(vec1)/sqrt(vec2) - 1)
    //       (sqrt(vec1)/sqrt(vec2) - 1) ~ (vec1/vec2 - 1) / 2
    //
    // we want to check: sqrt(vec1)/sqrt(vec2) - 1. < threshold
    // i.e., with the above,
    //       (vec1/vec2 - 1) / 2 < threshold
    //       (vec1/vec2 - 1) < threshold*2
    Numeric th = threshold * 2.;

    // No variation at the South pole!
    if (lat_grid[0] == -90) {
      Numeric vec1, vec2;
      for (Index ip = 0; ip < npages; ip++) {
        for (Index ic = 1; ic < ncols; ic++) {
          vec1 = x1(ip, 0, ic) * x1(ip, 0, ic) + x2(ip, 0, ic) * x2(ip, 0, ic);
          vec2 = x1(ip, 0, ic - 1) * x1(ip, 0, ic - 1) +
                 x2(ip, 0, ic - 1) * x2(ip, 0, ic - 1);
          ARTS_USER_ERROR_IF (fabs(vec1 / vec2 - 1.) > th,
            "The variables *", x1_name, "* and *", x2_name
              , "* are assumed\n"
              , "to be two horizontal components of a vector field.\n"
              , "At the pole, the data (here: the total length of\n"
              , "the horizontal vector) can NOT vary with longitude,\n"
              , "but this appears to be the case on the South pole.\n"
              , "The threshold is ", threshold, ", but the actual\n"
              , "deviation at pressure level ", ip, " and longitude\n"
              , "points ", ic - 1, " and ", ic, " is "
              , sqrt(vec1) / sqrt(vec2) - 1.)
        }
      }
    }
    // No variation at the North pole!
    if (lat_grid[nrows - 1] == 90) {
      Numeric vec1, vec2;
      const Index ir = nrows - 1;
      for (Index ip = 0; ip < npages; ip++) {
        for (Index ic = 1; ic < ncols; ic++) {
          vec1 =
              x1(ip, ir, ic) * x1(ip, ir, ic) + x2(ip, ir, ic) * x2(ip, ir, ic);
          vec2 = x1(ip, ir, ic - 1) * x1(ip, ir, ic - 1) +
                 x2(ip, ir, ic - 1) * x2(ip, ir, ic - 1);
          ARTS_USER_ERROR_IF (fabs(vec1 / vec2 - 1.) > th,
            "The variables *", x1_name, "* and *", x2_name
              , "* are assumed\n"
              , "to be two horizontal components of a vector field.\n"
              , "At the pole, the data (here: the total length of\n"
              , "the horizontal vector) can NOT vary with longitude,\n"
              , "but this appears to be the case on the North pole.\n"
              , "The threshold is ", threshold, ", but the actual\n"
              , "deviation at pressure level ", ip, " and longitude\n"
              , "points ", ic - 1, " and ", ic, " is "
              , sqrt(vec1) / sqrt(vec2) - 1.)
        }
      }
    }
  }
}

//! chk_atm_surface
/*! 
    Checks if a surface-type variable matches the dimensionality and the grids.

    The function gives an error message if this is not the case.

    \param    x_name       The name of the surface-type variable.
    \param    x            The variable holding the surface data.
    \param    dim          The atmospheric dimensionality.
    \param    lat_grid     The latitude grid.
    \param    lon_grid     The longitude grid.

    \author Patrick Eriksson 
    \date   2002-04-15
*/
void chk_atm_surface(const String& x_name,
                     const Matrix& x,
                     const Index& dim,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid) {
  Index ncols = 1, nrows = 1;
  if (dim > 1) nrows = lat_grid.nelem();
  if (dim > 2) ncols = lon_grid.nelem();
  ARTS_USER_ERROR_IF (x.ncols() != ncols || x.nrows() != nrows,
    "The surface variable *", x_name, "* has wrong size.\n"
      , "Expected size is ", nrows, " x ", ncols, ","
      , " while actual size is ", x.nrows(), " x ", x.ncols(), ".")

  // Special 3D checks:
  if (dim == 3) {
    // If all lons are covered, check if cyclic
    if ((lon_grid[ncols - 1] - lon_grid[0]) == 360) {
      const Index ic = ncols - 1;
      for (Index ir = 0; ir < nrows; ir++) {
        ARTS_USER_ERROR_IF (!is_same_within_epsilon(x(ir, ic), x(ir, 0), 2 * DBL_EPSILON),
          "The variable *", x_name, "* covers 360 "
            , "degrees in the longitude direction, but the field "
            , "seems to deviate between first and last longitude "
            , "point. The field must be \"cyclic\".")
      }
    }

    // No variation at the South pole!
    if (lat_grid[0] == -90) {
      for (Index ic = 1; ic < ncols; ic++) {
        ARTS_USER_ERROR_IF (!is_same_within_epsilon(x(0, ic), x(0, ic - 1), 2 * DBL_EPSILON),
          "The variable *", x_name, "* covers the South "
            , "pole. The data corresponding to the pole can not "
            , "vary with longitude, but this appears to be the "
            , "case.")
      }
    }
    // No variation at the North pole!
    if (lat_grid[nrows - 1] == 90) {
      const Index ir = nrows - 1;
      for (Index ic = 1; ic < ncols; ic++) {
        ARTS_USER_ERROR_IF (!is_same_within_epsilon(
                x(ir, ic), x(ir, ic - 1), 2 * DBL_EPSILON),
          "The variable *", x_name, "* covers the North "
            , "pole. The data corresponding to the pole can not "
            , "vary with longitude, but this appears to be the "
            , "case.")
      }
    }
  }
}

/*===========================================================================
  === Functions related to specif WSV
  ===========================================================================*/

void chk_rte_los(const String& name,
                 ConstVectorView los)
{
  ARTS_USER_ERROR_IF(los.nelem() != 2,
                     "The vector *", name, "* must have length 2.");
  ARTS_USER_ERROR_IF(los[0] < 0 || los[0] > 180,
      "The zenith angle in *", name, "* must be in the range [0,180].");
  ARTS_USER_ERROR_IF(los[1] < -180 || los[1] > 180,
      "The azimuth angle in *", name, "* must be in the range [-180,180].");
}


void chk_rte_pos(const String& name,
                 ConstVectorView pos)
{
  ARTS_USER_ERROR_IF(pos.nelem() != 3,
                      "The vector *", name, "* must have length 3.")
  ARTS_USER_ERROR_IF(pos[1] < -90 || pos[1] > 90,
      "The latitude in *", name, "* must be in the range [-90,90].")
  ARTS_USER_ERROR_IF(pos[2] < -180 || pos[2] >= 180,
      "The longitude in *", name, "* must be in the range [-180,180[].")
}


void chk_refellipsoid(const Vector2 refellipsoid)
{
  ARTS_USER_ERROR_IF(refellipsoid[0] <= 0 || refellipsoid[1] <= 0,
                     "All elements of *refellipsoid* must be > 0.");
  ARTS_USER_ERROR_IF(abs(refellipsoid[1]/refellipsoid[0]-1) > 0.5,
      "The ratio of the two radii in *refellipsoid* is outisde of [0.5,1.5].\n"
      "Do you really want to have such a flat reference ellipsoid?");
}


void chk_sensor_pos(const String& name,
                    ConstMatrixView sensor_pos)
{
  ARTS_USER_ERROR_IF(sensor_pos.ncols() != 3,
                     "*", name, "* must have three columns.");
  ARTS_USER_ERROR_IF(sensor_pos.nrows() == 0,
                     "*", name, "*must have at least one row.");
  for (Index i=0; i<sensor_pos.nrows(); i++) {
    ARTS_USER_ERROR_IF(sensor_pos(i,1) < -90 || sensor_pos(i,1) > 90,
                       "Unvalid latitude in *", name, "*.\n",
                       "Latitudes must be inside the range [-90,90],\n",
                       "but ", name, "(", i, ",1) is ", sensor_pos(i,1));
    ARTS_USER_ERROR_IF(sensor_pos(i,2) < -180 || sensor_pos(i,1) >= 180,
                       "Unvalid longitude in *", name, "*.\n",
                       "Longitudes must be inside the range [-180,180[,\n",
                       "but ", name, "(", i, ",2) is ", sensor_pos(i,2));
  }
}


void chk_sensor_los(const String& name,
                    ConstMatrixView sensor_los)
{
  ARTS_USER_ERROR_IF(sensor_los.ncols() != 2,
                     "*", name, "* must have two columns.");
  ARTS_USER_ERROR_IF(sensor_los.nrows() == 0,
                     "*", name, "* must have at least one row.");
  for (Index i=0; i<sensor_los.nrows(); i++) {
    ARTS_USER_ERROR_IF(sensor_los(i,0) < 0 || sensor_los(i,0) > 180,
                       "Unvalid zenith angle in *", name, "*.\n",
                       "Zenith angles must be inside the range [0,180],\n"
                       "but ", name, "(", i, ",0) is ", sensor_los(i,0));
    ARTS_USER_ERROR_IF(sensor_los(i,1) < -180 || sensor_los(i,1) > 180,
                       "Unvalid azimuth angle in *", name, "*.\n",
                       "Azimuth angles must be inside the range [-180,180],\n"
                       "but ", name, "(", i, ",1) is ", sensor_los(i,1));
  }
}


void chk_sensor_poslos(const String& name1,
                       ConstMatrixView sensor_pos,
                       const String& name2,
                       ConstMatrixView sensor_los) {
  chk_sensor_pos(name1, sensor_pos);
  chk_sensor_los(name2, sensor_los);
  ARTS_USER_ERROR_IF(sensor_los.nrows() != sensor_pos.nrows(),
                     "*", name1, "* and *", name2,
                     "* must have the same number of rows.");
}


void chk_surface_elevation(const GriddedField2& surface_elevation)
{
  ARTS_USER_ERROR_IF(surface_elevation.get_grid_name(0) != "Latitude",
                     "Name of first grid must be \"Latitude\".");
  ARTS_USER_ERROR_IF(surface_elevation.get_grid_name(1) != "Longitude",
                     "Name of second grid must be \"Longitude\".");
  const Vector& lat_grid = surface_elevation.get_numeric_grid(0);
  ARTS_USER_ERROR_IF(surface_elevation.data.nrows() != lat_grid.nelem(),
                     "Inconsistent latitude size in *surface_elevation*\n"
                     "Length of latitude grid: ", lat_grid.nelem(), "\n"
                     "Latitude size of data: ", surface_elevation.data.nrows());
  const Vector& lon_grid = surface_elevation.get_numeric_grid(1);
  ARTS_USER_ERROR_IF(surface_elevation.data.ncols() != lon_grid.nelem(),
                     "Inconsistent longitude size in *surface_elevation*\n"
                     "Length of longitude grid: ", lon_grid.nelem(), "\n"
                     "Longitude size of data: ", surface_elevation.data.ncols());
  ARTS_USER_ERROR_IF(surface_elevation.data.empty(),
                     "The data in *surface_elevation* are empty. Not allowed!");
}












// The ones below should be deleted ZZZ

//! chk_rte_pos
/*! 
    Performs all needed checks of rte_pos and rte_pos2.

    The function gives an error message if this is not the case.

    \param    rte_pos          As WSV rte_pos or rte_pos2.
    \param    is_rte_pos2      True if rte_pos actually is rte_pos2.

    \author Patrick Eriksson 
    \date   2012-03-26
*/
void chk_rte_pos(ConstVectorView rte_pos,
                 const bool& is_rte_pos2) {
  String vname = "*rte_pos*";
  if (is_rte_pos2) {
    vname = "*rte_pos2*";
  }
  ARTS_USER_ERROR_IF (rte_pos.nelem() != 3,
    "For 3D, ", vname, " must have length 3.")
  ARTS_USER_ERROR_IF (rte_pos[1] < -90 || rte_pos[1] > 90,
    "The (3D) latitude in ", vname, " must be in the "
      , "range [-90,90].")
  ARTS_USER_ERROR_IF (rte_pos[2] < -360 || rte_pos[2] > 360,
    "The longitude in ", vname, " must be in the "
      , "range [-360,360].")
}

//! chk_rte_los
/*! 
    Performs all needed checks of rte_los

    The function gives an error message if this is not the case.

    \param    rte_los          As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2012-03-26
*/
void chk_rte_los(ConstVectorView rte_los) {
  ARTS_USER_ERROR_IF (rte_los.nelem() != 2,
        "For 3D, los-vectors must have length 2.");
  ARTS_USER_ERROR_IF (rte_los[0] < 0 || rte_los[0] > 180,
        "For 3D, the zenith angle of a los-vector must "
        "be in the range [0,180].");
  ARTS_USER_ERROR_IF (rte_los[1] < -180 || rte_los[1] > 180,
        "For 3D, the azimuth angle of a los-vector must "
        "be in the range [-180,180].");
}


/*===========================================================================
 === Functions related to GriddedFields.
 ===========================================================================*/

//! Check name of grid in GriddedField.
/**
 Does a case-insensitive check to verify that the name of the grid at the
 given index in the GriddedField has the expected name.

 \param gf         GriddedField to check.
                   Can be a GriddedField of any dimension.
 \param gridindex  Index of grid.
 \param gridname   Expected name of grid.
 */
void chk_griddedfield_gridname(const GriddedField& gf,
                               const Index gridindex,
                               const String& gridname) {
  ARTS_USER_ERROR_IF (gf.get_dim() - 1 < gridindex,
    "Grid index ", gridindex, " exceeds dimension of GriddedField",
    (gf.get_name().size()) ? var_string(" \"", gf.get_name(), "\"") : var_string())

  String gfgridnameupper = gf.get_grid_name(gridindex);
  toupper(gfgridnameupper);

  String gridnameupper = gridname;
  toupper(gridnameupper);

  ARTS_USER_ERROR_IF (gfgridnameupper != gridnameupper,
    "Name of grid ", gridindex, " in GriddedField",
    (gf.get_name().size()) ? var_string(" \"", gf.get_name(), "\"") : var_string(),
    " is \"", gf.get_grid_name(gridindex), "\".\n"
      , "The expected name is \"", gridname, "\".")
}

/*===========================================================================
 === Functions checking sensor
 ===========================================================================*/

/** Check met_mm_backend.

 Verifies that the backend description matrix has the correct size and format.

 \param[in] mmb met_mm_backend

 \throws std::runtime_error

 \author Oliver Lemke
 */
void chk_met_mm_backend(const Matrix& mmb) {
  ARTS_USER_ERROR_IF (!mmb.nrows(),
    "No channels defined in *met_mm_backend*.");

  ARTS_USER_ERROR_IF (mmb.ncols() != 4,
    "*met_mm_backend* must have 4 columns.");

  for (Index ch = 0; ch < mmb.nrows(); ch++) {
    Numeric lo = mmb(ch, 0);
    Numeric offset1 = mmb(ch, 1);
    Numeric offset2 = mmb(ch, 2);
    Numeric bandwidth = mmb(ch, 3);

    // Negative LO
    ARTS_USER_ERROR_IF (lo < 0.,
      "Error in channel ", ch + 1, " at row ", ch
        , " in *met_mm_backend*.\n"
        , "Center frequency is negative: ", mmb(ch, 0), " Hz")

    // Negative offsets
    ARTS_USER_ERROR_IF (offset1 < 0. || offset2 < 0.,
      "Error in channel ", ch + 1, " at row ", ch
        , " in *met_mm_backend*.\n"
        , "Offset is negative:\n"
        , "offset1: ", offset1, " Hz\n"
        , "offset2: ", offset2, " Hz\n")

    // First offset is smaller than second offset
    ARTS_USER_ERROR_IF (offset1 != 0. && offset1 <= offset2,
      "Error in channel ", ch + 1, " at row ", ch
        , " in *met_mm_backend*.\n"
        , "First passband offset is smaller than/equal to the second offset:\n"
        , "offset1: ", offset1, " Hz\n"
        , "offset2: ", offset2, " Hz\n")

    // Bandwidth too wide, overlap with LO
    ARTS_USER_ERROR_IF (offset1 > 0 && offset1 - offset2 <= bandwidth / 2.,
      "Error in channel ", ch + 1, " at row ", ch
        , " in *met_mm_backend*.\n"
        , "Band touches or overlaps with the center frequency:\n"
        , "offset1                        : ", offset1, " Hz\n"
        , "offset2                        : ", offset2, " Hz\n"
        , "bandwidth                      : ", bandwidth, " Hz\n"
        , "offset1 - offset2 - bandwidth/2: "
        , offset1 - offset2 - bandwidth / 2., " Hz\n")

    // Bandwidth too wide, passbands overlap
    ARTS_USER_ERROR_IF (offset2 > 0 && offset2 <= bandwidth / 2.,
      "Error in channel ", ch + 1, " at row ", ch
        , " in *met_mm_backend*.\n"
        , "Bands overlap or touch, offset2 > bandwidth/2:\n"
        , "offset2    : ", offset2, " Hz\n"
        , "bandwidth/2: ", bandwidth / 2., " Hz\n")

    // Channel too wide, goes negative
    ARTS_USER_ERROR_IF (lo - offset1 - offset2 - bandwidth / 2. <= 0,
      "Error in channel ", ch + 1, " at row ", ch
        , " in *met_mm_backend*.\n"
        , "Band too wide, reaches/exceeds 0 Hz:\n"
        , "LO                                  : ", lo, " Hz\n"
        , "offset1                             : ", offset1, " Hz\n"
        , "offset2                             : ", offset2, " Hz\n"
        , "bandwidth                           : ", bandwidth, " Hz\n"
        , "LO - offset1 - offset2 - bandwidth/2: "
        , lo - offset1 - offset2 - bandwidth / 2., " Hz\n")
  }
}
