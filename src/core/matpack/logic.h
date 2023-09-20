/*!
  \file   logic.h
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Fri May  3 19:10:04 2002
  
  \brief  Header file for logic.cc
*/

#ifndef logic_h
#define logic_h

#include <array.h>
#include "matpack_data.h"

bool is_bool(const Index& x);

bool is_multiple(const Index& x, const Index& y);

bool is_size(ConstVectorView x, const Index& l);

bool is_size(ConstMatrixView x, const Index& r, const Index& c);

bool is_size(ConstTensor3View x,
             const Index& p,
             const Index& r,
             const Index& c);

bool is_size(ConstTensor4View x,
             const Index& b,
             const Index& p,
             const Index& r,
             const Index& c);

bool is_size(ConstTensor5View x,
             const Index& s,
             const Index& b,
             const Index& p,
             const Index& r,
             const Index& c);

bool is_size(ConstTensor6View x,
             const Index& v,
             const Index& s,
             const Index& b,
             const Index& p,
             const Index& r,
             const Index& c);

bool is_size(ConstTensor7View x,
             const Index& l,
             const Index& v,
             const Index& s,
             const Index& b,
             const Index& p,
             const Index& r,
             const Index& c);

bool is_sorted(ConstVectorView x);

bool is_increasing(ConstVectorView x);

/*! Checks if the vector is increasing with a regular interval

  The check is abs((x[i] - x[i-1]) - (x[1] - x[0])) < epsilon

  @param[in] x A vector
  @param[in] epsilon A comparison
*/
bool is_regularly_increasing_within_epsilon(ConstVectorView x, const Numeric epsilon=1e-8);

bool is_increasing(const ArrayOfIndex& x);

bool is_decreasing(ConstVectorView x);

bool is_unique(const ArrayOfIndex& x);

bool is_singular(ConstMatrixView A);

bool is_diagonal(ConstMatrixView A);

bool is_same_within_epsilon(const Numeric& a,
                            const Numeric& b,
                            const Numeric& epsilon);

bool is_lon_cyclic(ConstVectorView grid, const Numeric& epsilon = 0.001);

////////////////////////////////////////////////////////////////////////////
//   Template functions (have to be here in the .h file).
////////////////////////////////////////////////////////////////////////////

//! Verifies that the size of x is n.
/*! 
  This function is supposed to be used together with ARTS_ASSERT like this:
  ARTS_ASSERT(is_size(x,n)). It works for any array type.

  \param  x The Array to check.
  \param  n The desired length.
  \return True if the size of x is n.
*/
template <class T>
bool is_size(const Array<T>& x, const Size& n) {
  return (n == x.size());
}

#endif  // logic_h
