#pragma once

#include <algorithm>
#include <array>
#include <climits>
#include <compare>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "configtypes.h"

/** An array of Index. */
template <typename base>
using Array = std::vector<base>;

/** An array of Index. */
using ArrayOfIndex = Array<Index>;

using ArrayOfArrayOfIndex = Array<ArrayOfIndex>;

/** An array of Numeric. */
using ArrayOfNumeric = Array<Numeric>;

namespace std {
std::ostream& operator<<(std::ostream& os, const ArrayOfIndex& x);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfIndex& x);
}  // namespace std

std::ostream& operator<<(std::ostream& os, const ArrayOfNumeric& x);

/** Max function. */
template <class base>
constexpr base max(const Array<base>& x) {
  return *std::ranges::max_element(x);
}

/** Min function. */
template <class base>
constexpr base min(const Array<base>& x) {
  return *std::ranges::min_element(x);
}

//! Find first occurance.
/*!
  This returns the index of the first occurance of w in 
  array x.  

  A return value of -1 indicates that no matching element was found.

  \return The index of the thing we looked for.
  \param  x   The array to search.
  \param w The value to look for.

  \author Stefan Buehler
  \date   2002-11-28
*/
template <class base>
Index find_first(const Array<base>& x, const base& w) {
  for (Size i = 0; i < x.size(); ++i)
    if (w == x[i]) return i;

  return -1;
}

//! Find all occurances.
/*!
  This calculates an array of indices of all occurances of w in
  array x.

  An empty output array means that no occurance was found.

  \retval pos Array with positions of w in the array.
  \param  x   The array to search.
  \param  w   The value to look for.

  \author Stefan Buehler
  \date   2002-11-28
*/
template <class base>
void find_all(ArrayOfIndex& pos, const Array<base>& x, const base& w) {
  pos.resize(0);
  for (Size i = 0; i < x.size(); ++i)
    if (w == x[i]) pos.push_back(i);
}

/** Helper comparison class to sort an array or vector based on an ArrayOfNumeric.

 Usage:

 Vector v1
 ArrayOfNumeric v2;
 ...
 std::sort(v1.begin(), v1.end(), CmpArrayOfNumeric(v2));

 Source: http://stackoverflow.com/questions/8147911/locking-two-vectors-and-sorting-them
 */
class CmpArrayOfNumeric {
 public:
  CmpArrayOfNumeric(const ArrayOfNumeric& vec) : values(vec) {}
  bool operator()(const int& a, const int& b) const {
    return values[a] < values[b];
  }

  const ArrayOfNumeric& values;
};

//! Determine total number of elements in an ArrayOfArray
template <class base>
Index TotalNumberOfElements(const Array<Array<base> >& aa) {
  Index N_aa = 0;
  for (Size i = 0; i < aa.size(); i++) {
    N_aa += aa[i].size();
  }

  return N_aa;
}

//! Determine the index of an element in a flattened version of the array
template <class base>
Index FlattenedIndex(const Array<Array<base> >& aa,
                     Size outer,
                     Size inner = 0) {
  ARTS_ASSERT(outer < aa.size());
  ARTS_ASSERT(inner < aa[outer].size());

  Index N = 0;
  for (Size i = 0; i < outer; i++) {
    N += aa[i].size();
  }

  return N + inner;
}

// It is not a good idea to put all the predefined array types in one
// place. If I do this than a file cannot use one without defining all
// the others.

//! Make a std::array of a list of variables (must be 1-long at least)
template <typename T, typename... Ts>
constexpr std::array<T, 1 + sizeof...(Ts)> stdarrayify(const T& first,
                                                       const Ts&... the_rest) {
  return {first, T(the_rest)...};
}

template <typename T>
std::string stringify(const Array<T>& list,
                      const char* const sep = " ",
                      const char* const beg = "") {
  std::ostringstream os;
  for (auto& x : list) os << beg << x << sep;
  return os.str();
}
