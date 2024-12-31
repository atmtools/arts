#pragma once

#include <configtypes.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

/*! Check, if two numbers agree within a given epsilon.

  This logical function verifies if two numbers are the same for the
  desired number of digits. The comparison statement comes from
  Oliver: ( abs(a-b) <= epsilon * max(a,b) )

  Modified to make sure that negative numbers are also treated correctly.

  The variable epsilon gives the number of digits used for the
  comparison. (epsilon = 0.0001 for a comparison up to the 5th digit)

  @param a A number.
  @param b Another number.
  @param epsilon The epsilon of the required agreement.

  @return True if the two numbers are the same.
 */
template <typename T>
constexpr bool is_same_within_epsilon(T a,
                                      T b,
                                      T e = std::numeric_limits<T>::epsilon()) {
  return std::abs(a - b) <= e * std::max(std::abs(a), std::abs(b));
}

template <typename T>
constexpr bool in_range(T x, T x_low, T x_high) {
  return (x >= x_low) and (x <= x_high);
}

/*! Checks if an ArrayOfIndex is unique, i.e., has no duplicate values
  
  This only returns true if the array does not contain any duplicate
  values.  
  
  \return      True if unique, otherwise false.
  \param   x   An ArrayOfIndex.
  
  \author Stefan Buehler
  \date   2008-08-24

*/
template <typename T>
constexpr bool is_unique(const std::vector<T>& x) {
  // We simply compare the second element to the first,
  // the third to the first and second, and so on.

  for (Size i = 1; i < x.size(); ++i)
    for (Size s = 0; s < i; ++s)
      if (x[i] == x[s]) return false;

  return true;
}
