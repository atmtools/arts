////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
  \file   sorting.h

  Contains sorting routines.

  \author Oliver Lemke
  \date 2003-08-20
  */

#ifndef sorting_h
#define sorting_h

#include <algorithm>
#include <functional>

#include "array.h"
#include "matpack_concepts.h"


/** get_sorted_indexes
 *
 * The output array contains the sorted indexes of the input data.
 *
 * The following member functions must be defined for the input data type:
 * T.begin, T.end, T.operator[], T.operator<
 *
 * Furthermore, it must provide the type T::const_iterator.
 *
 * \param sorted  Output array with sorted indexes
 * \param data    Data to sort
 *
 * \author Oliver Lemke <olemke@core-dump.info>
 * \date   2003-08-20
 */
template <typename T>
void get_sorted_indexes(ArrayOfIndex& sorted, const T& data) {
  sorted.resize(data.size());
  std::iota(sorted.begin(), sorted.end(), 0);
  sort(sorted.begin(), sorted.end(), [&data](const Index a, const Index b) {
    return data[a] < data[b];
  });
}

#endif /* sorting_h */
