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
#include <numeric>

#include "array.h"

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

template <typename T>
concept bubble_sorting_function = requires(T a) {
  { a(Size{}, Size{}) } -> std::same_as<bool>;
};

template <typename T>
concept bubble_sortable = requires(T& a) {
  std::swap(a[Size{}], a[Size{}]);
  { a.size() } -> std::convertible_to<Size>;
};

template <bubble_sorting_function SortFn,
          bubble_sortable T,
          bubble_sortable... Ts>
void bubble_sort_by(const SortFn& sort_fn,
                    T& data1,
                    Ts&... data_n) ARTS_NOEXCEPT {
  const Size n = static_cast<Size>(data1.size());
  ARTS_ASSERT(n == static_cast<Size>(data_n.size()) && ...,
              "All data must have the same size.");
  for (Size i = 0; i < n; i++) {
    for (Size j = i + 1; j < n; j++) {
      if (sort_fn(i, j)) {
        (std::swap(data_n[i], data_n[j]), ...);
      }
    }
  }
}

#endif /* sorting_h */
