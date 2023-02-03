/* Copyright (C) 2001-2012 Stefan Buehler <sbuehler@ltu.se>

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

/**
    \file   array.h

   This file contains the definition of Array.

   \author Stefan Buehler
   \date   2001-09-05
*/

#ifndef array_h
#define array_h

#include <array>
#include <climits>
#include <iomanip>
#include <iostream>
#include <vector>

#include "matpack_concepts.h"

/** This can be used to make arrays out of anything.

    A simple \#define does not do for this, since I have to implement
    member functions like nelem, in order to be consistent with
    Vector.

    Because constructors are not inherited, I have to re-define all
    constructors.
*/
template <class base>
class Array : public std::vector<base> {
 public:
  // Constructors:
  Array() : std::vector<base>() {}

  explicit Array(Index n) : std::vector<base>(n) {}

  Array(Index n, const base& fillvalue) : std::vector<base>(n) {
    std::fill(this->begin(), this->end(), fillvalue);
  }

  Array(const Array& A) : std::vector<base>(A) {}

  Array(Array&& A) noexcept : std::vector<base>(std::move(A)) {}

  Array(std::initializer_list<base> init) : std::vector<base>(init) {}

  template <class base2, size_t N>
  explicit Array(const std::array<base2, N>& input)
      : std::vector<base>(input.begin(), input.end()) {
    static_assert(std::is_convertible<base2, base>::value,
                  "Must be convertible");
  }

  Array(std::vector<base> x) : std::vector<base>(std::move(x)) {}

  // Assignment operators:
  Array& operator=(base x) {
    std::fill(this->begin(), this->end(), x);
    return *this;
  }

  Array& operator=(const Array& A) {
    this->resize(A.size());
    std::copy(A.begin(), A.end(), this->begin());
    return *this;
  }

  Array& operator=(Array&& A) noexcept {
    std::vector<base>::operator=(std::move(A));
    return *this;
  }

  // Number of elements:
  [[nodiscard]] Index nelem() const ARTS_NOEXCEPT {
    size_t s = this->size();
    ARTS_ASSERT(s < LONG_MAX);
    return static_cast<Index>(s);
  }

  // Index operators:
  const base& operator[](const Index n) const {
    ARTS_ASSERT(0 <= n);
    ARTS_ASSERT(n < nelem());
    return std::vector<base>::operator[](n);
  }

  base& operator[](const Index n) {
    ARTS_ASSERT(0 <= n);
    ARTS_ASSERT(n < nelem());
    return std::vector<base>::operator[](n);
  }

  // Helper functions
  void push_back_n(const base& elem, const Index n) {
    for (Index i = 0; i < n; i++) std::vector<base>::push_back(elem);
  }

  virtual ~Array() = default;

  friend std::ostream& operator<<(std::ostream& os, const Array& v) {
    typename Array::const_iterator i = v.begin();
    const typename Array::const_iterator end = v.end();

    if (i != end) {
      os << std::setw(3) << *i;
      ++i;
    }

    for (; i != end; ++i) {
      os << " " << std::setw(3) << *i;
    }

    return os;
  }
};

/** An array of Index. */
using ArrayOfIndex = Array<Index>;

using ArrayOfArrayOfIndex = Array<ArrayOfIndex>;

/** An array of Numeric. */
using ArrayOfNumeric = Array<Numeric>;

/** Max function. */
template <class base>
inline base max(const Array<base>& x) {
  // Initial value for max:
  base max = x[0];

  typename Array<base>::const_iterator xi = x.begin();
  const typename Array<base>::const_iterator xe = x.end();

  for (; xi != xe; ++xi) {
    if (*xi > max) max = *xi;
  }

  return max;
}

/** Min function. */
template <class base>
inline base min(const Array<base>& x) {
  // Initial value for min:
  base min = x[0];

  typename Array<base>::const_iterator xi = x.begin();
  const typename Array<base>::const_iterator xe = x.end();

  for (; xi != xe; ++xi) {
    if (*xi < min) min = *xi;
  }

  return min;
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
  for (Index i = 0; i < x.nelem(); ++i)
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
  for (Index i = 0; i < x.nelem(); ++i)
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
  for (Index i = 0; i < aa.nelem(); i++) {
    N_aa += aa[i].nelem();
  }

  return N_aa;
}

//! Determine the index of an element in a flattened version of the array
template <class base>
Index FlattenedIndex(const Array<Array<base> >& aa,
                     Index outer,
                     Index inner = 0) {
  ARTS_ASSERT(outer < aa.nelem());
  ARTS_ASSERT(inner < aa[outer].nelem());

  Index N = 0;
  for (Index i = 0; i < outer; i++) {
    N += aa[i].nelem();
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

#endif  // array_h
