/*===========================================================================
  ===  File description 
  ===========================================================================*/

/*!
   \file   check_input.h
   \author Patrick Eriksson <patrick.eriksson@chalmers.se>
   \date   2002-04-15 

   This file contains the declaration of functions in check_input.cc.
*/

#ifndef checkinput_h
#define checkinput_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "agenda_class.h"
#include "exceptions.h"
#include "gridded_fields.h"
#include "matpack_data.h"
#include "mystring.h"
//#include <cfloat>

/*===========================================================================
  === Functions in check_input.cc
  ===========================================================================*/

void chk_if_bool(const String& x_name, const Index& x);

void chk_if_in_range(const String& x_name,
                     const Index& x,
                     const Index& x_low,
                     const Index& x_high);

void chk_if_positive(const String& x_name, const Numeric& x);

void chk_if_increasing(const String& x_name, const ArrayOfIndex& x);

void chk_not_negative(const String& x_name, const Numeric& x);

void chk_if_in_range(const String& x_name,
                     const Numeric& x,
                     const Numeric& x_low,
                     const Numeric& x_high);

void chk_if_in_range_exclude_low(const String& x_name,
                                 const Numeric& x,
                                 const Numeric& x_low,
                                 const Numeric& x_high);

void chk_if_in_range_exclude_high(const String& x_name,
                                  const Numeric& x,
                                  const Numeric& x_low,
                                  const Numeric& x_high);

void chk_if_in_range_exclude(const String& x_name,
                             const Numeric& x,
                             const Numeric& x_low,
                             const Numeric& x_high);

void chk_vector_length(const String& x_name, ConstVectorView x, const Index& l);

void chk_vector_length(const String& x1_name,
                       const String& x2_name,
                       ConstVectorView x1,
                       ConstVectorView x2);

void chk_if_increasing(const String& x_name, ConstVectorView x);

void chk_if_decreasing(const String& x_name, ConstVectorView x);

void chk_if_equal(const String& x1_name,
                  const String& x2_name,
                  ConstVectorView v1,
                  ConstVectorView v2,
                  Numeric margin = 1e-6);

void chk_matrix_ncols(const String& x_name, ConstMatrixView x, const Index& l);

void chk_matrix_nrows(const String& x_name, ConstMatrixView x, const Index& l);

/** Subclasses of runtime_error.
 
 This is so that I can distinguish what went wrong in chk_contains.
 
 \author Stefan Buehler
 \date   2013-04-23 */
class runtime_error_not_found : public runtime_error {
 public:
  runtime_error_not_found(const string& s) : runtime_error(s) {}
};

/** Subclasses of runtime_error.
 
 This is so that I can distinguish what went wrong in chk_contains.
 
 \author Stefan Buehler
 \date   2013-04-23 */
class runtime_error_not_unique : public runtime_error {
 public:
  runtime_error_not_unique(const string& s) : runtime_error(s) {}
};

/*===========================================================================
  === Template Functions for Arrays
  ===========================================================================*/

//! Check if an array contains a value.
/*!
  This makes sure that the array *x* contains the element with
  value *what* exactly once.

  As a bonus, it returns the index of *what* in *x*.

  This template function can be used for arrays of anything, provided
  that the "==" operator is defined.

  \return The index of the thing we looked for.
  \param x_name Name of the array to check
  \param x The array to check
  \param what The value to look for.
 
  \throw runtime_error_not_found, runtime_error_not_unique

  \author Stefan Buehler
  \date   2002-11-28
*/
template <class T>
Index chk_contains(const String& x_name, const Array<T>& x, const T& what) {
  // To generate error messages:
  ostringstream os;

  // To store the positions:
  ArrayOfIndex pos;

  // Find all positions of what in x and store in pos:
  find_all(pos, x, what);

  switch (pos.nelem()) {
    case 0:
      // Not found.
      os << "The array *" << x_name << "* must contain the element " << what
         << ",\n"
         << "but it does not.";
      throw runtime_error_not_found(os.str());
      break;

    case 1:
      // Found once, this is what we want!
      return pos[0];

    default:
      // Found more than once.
      os << "The array *" << x_name << "* must contain the element " << what
         << "\n"
         << "exactly once, but it does contain it " << pos.nelem() << " times.";
      throw runtime_error_not_unique(os.str());
      break;
  }

  return -1;
}

//! Check the size of an array.
/*! 
    Checks the size of an Array. Cloned from
    Patricks similar function for Vector.

    The function throws a runtime_error if the size is not correct.  

    This is a template function that works for any array type.

    \param    x_name   The name of the variable.
    \param    x        A variable of type ArrayOfIndex.
    \param    c        The size to match

    \author Stefan Buehler
    \date   2007-05-18
*/
template <class T>
void chk_size(const String& x_name, const Array<T>& x, const Index& c) {
  ARTS_USER_ERROR_IF (x.nelem() != c,
      "The array *", x_name, "*\n"
      "does not have the right size.\n"
      "The size should be: ", c,"\n"
      "but it is:          ", x.nelem())
}

/*===========================================================================
  === Functions for Tensors
  ===========================================================================*/

void chk_size(const String& x_name, ConstVectorView x, const Index& c);

void chk_size(const String& x_name,
              ConstMatrixView x,
              const Index& r,
              const Index& c);

void chk_size(const String& x_name,
              ConstTensor3View x,
              const Index& p,
              const Index& r,
              const Index& c);

void chk_size(const String& x_name,
              ConstTensor4View x,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c);

void chk_size(const String& x_name,
              ConstTensor5View x,
              const Index& s,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c);

void chk_size(const String& x_name,
              ConstTensor6View x,
              const Index& v,
              const Index& s,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c);

void chk_size(const String& x_name,
              ConstTensor7View x,
              const Index& l,
              const Index& v,
              const Index& s,
              const Index& b,
              const Index& p,
              const Index& r,
              const Index& c);

void chk_not_empty(const String& x_name, const Agenda& x);

void chk_interpolation_grids_loose(Index& ing_min,
                                   Index& ing_max,
                                   const String& which_interpolation,
                                   ConstVectorView old_grid,
                                   ConstVectorView new_grid,
                                   ConstVectorView data,
                                   const Index order = 1);

void chk_interpolation_grids_loose_no_data_check(
    Index& ing_min,
    Index& ing_max,
    const String& which_interpolation,
    ConstVectorView old_grid,
    ConstVectorView new_grid,
    const Index order = 1);

void chk_interpolation_pgrids_loose_no_data_check(
    Index& ing_min,
    Index& ing_max,
    const String& which_interpolation,
    ConstVectorView old_pgrid,
    ConstVectorView new_pgrid,
    const Index order = 1);

void chk_interpolation_grids_loose_check_data(Index& ing_min,
                                              Index& ing_max,
                                              const String& which_interpolation,
                                              ConstVectorView old_grid,
                                              ConstVectorView new_grid,
                                              ConstVectorView data);

void chk_interpolation_grids(const String& which_interpolation,
                             ConstVectorView old_grid,
                             ConstVectorView new_grid,
                             const Index order = 1,
                             const Numeric& extpolfac = 0.5,
                             const bool islog = false);

void chk_interpolation_grids(const String& which_interpolation,
                             ConstVectorView old_grid,
                             const Numeric& new_grid,
                             const Index order = 1,
                             const Numeric& extpolfac = 0.5);

void chk_interpolation_pgrids(const String& which_interpolation,
                              ConstVectorView old_pgrid,
                              ConstVectorView new_pgrid,
                              const Index order = 1,
                              const Numeric& extpolfac = 0.5);

void chk_atm_grids(const Index& dim,
                   ConstVectorView p_grid,
                   ConstVectorView lat_grid,
                   ConstVectorView lon_grid);

void chk_atm_field(const String& x_name,
                   ConstTensor3View x,
                   const Index& dim,
                   ConstVectorView p_grid,
                   ConstVectorView lat_grid,
                   ConstVectorView lon_grid,
                   const bool& chk_lat90 = 1);

void chk_atm_field(const String& x_name,
                   ConstTensor4View x,
                   const Index& dim,
                   const Index& nspecies,
                   ConstVectorView p_grid,
                   ConstVectorView lat_grid,
                   ConstVectorView lon_grid,
                   const bool& check_nan = 1);

void chk_atm_vecfield_lat90(const String& x1_name,
                            ConstTensor3View x1,
                            const String& x2_name,
                            ConstTensor3View x2,
                            const Index& dim,
                            ConstVectorView lat_grid,
                            const Numeric& threshold = 1e-3);
//        const Numeric&    threshold = 2*DBL_EPSILON );

void chk_latlon_true(const Index& atmosphere_dim,
                     ConstVectorView lat_grid,
                     ConstVectorView lat_true,
                     ConstVectorView lon_true);

void chk_atm_surface(const String& x_name,
                     const Matrix& x,
                     const Index& dim,
                     ConstVectorView lat_grid,
                     ConstVectorView lon_grid);

void chk_rte_pos(const Index& atmosphere_dim,
                 ConstVectorView rte_pos,
                 const bool& is_rte_pos2 = false);

void chk_rte_los(const Index& atmosphere_dim, ConstVectorView rte_los);

void chk_refellipsoid(ConstVectorView refellipsoid);
void chk_refellipsoidZZZ(ConstVectorView refellipsoid);
                           
void chk_griddedfield_gridname(const GriddedField& gf,
                               const Index gridindex,
                               const String& gridname);

void chk_met_mm_backend(const Matrix& bdsp);

#endif  // checkinput_h
