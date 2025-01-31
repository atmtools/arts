#pragma once

#include <configtypes.h>
#include <debug.h>

#include <algorithm>
#include <vector>

/** An array. */
template <typename base>
using Array = std::vector<base>;

/** An array of Index. */
using ArrayOfIndex = Array<Index>;

using ArrayOfArrayOfIndex = Array<ArrayOfIndex>;

/** An array of Numeric. */
using ArrayOfNumeric = Array<Numeric>;

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

//! Determine total number of elements in an ArrayOfArray
template <class base>
Index TotalNumberOfElements(const Array<Array<base>>& aa) {
  Index N_aa = 0;
  for (Size i = 0; i < aa.size(); i++) {
    N_aa += aa[i].size();
  }

  return N_aa;
}

namespace arr {
enum class CheckStatus : char { NotUnique, Unique, NotFound };

/*! Return the possible position and status of a value in an array.

  The first return argument is -1 if the value is not found, otherwise
  it points to the first occurance of the value. The second return
  argument indicates if the value is unique in the array or not.

  @param x The array to check
  @param what The value to look for.

  @return The index of the thing we looked for and the status of the value.

  @author Richard Larsson
  @date   2024-12-23
*/
template <class T>
constexpr std::pair<Index, CheckStatus> contains(const Array<T>& x,
                                                 const T& what) {
  auto pos1 = std::ranges::find(x, what);

  if (pos1 == x.end()) return {-1, CheckStatus::NotFound};

  auto pos2 = std::ranges::find(std::next(pos1), x.end(), what);

  const auto indpos = std::distance(x.begin(), pos1);
  if (pos2 == x.end()) return {indpos, CheckStatus::Unique};
  return {indpos, CheckStatus::NotUnique};
}

template <typename T>
concept sizeable = requires(T& x) { x.size(); };

/** Check that all have the same size. */
template <sizeable T, sizeable... Rest>
constexpr bool same_size(const T& x, const Rest&... rest)
  requires(sizeof...(rest) > 0)
{
  return ((x.size() == rest.size()) && ...);
}

/** Check that all have the same size. */
template <sizeable T, sizeable... Rest>
constexpr bool each_same_size(const std::vector<T>& x,
                              const std::vector<Rest>&... rest)
  requires(sizeof...(rest) > 0)
{
  if (not same_size(x, rest...)) return false;

  for (Size i = 0; i < x.size(); i++) {
    if (not same_size(x[i], rest[i]...)) return false;
  }
  return true;
}

template <typename T>
concept elemwise_sizeable = sizeable<T> and requires(T& x) { x[0].size(); };

/** Check that all and their elements have the same size. */
template <elemwise_sizeable T, elemwise_sizeable... Rest>
constexpr bool elemwise_same_size(const T& x, const Rest&... rest) {
  if constexpr (sizeof...(rest) != 0) {
    if (not same_size(x, rest...)) return false;
  }

  for (Size i = 0; i < x.size(); i++) {
    if (not same_size(x[i], x[0], rest[i]...)) return false;
  }

  return true;
}

template <typename T, typename U>
concept resizeable = requires(T& x, U sz) { x.resize(sz); };

//! Resize all the vectors/arrays using .resize(sz)
template <typename SZ, resizeable<SZ> T, resizeable<SZ>... Rest>
void resize(const SZ& sz, T& x, Rest&... rest) {
  x.resize(sz);

  if constexpr (sizeof...(rest) > 0) {
    resize(sz, rest...);
  }
}

template <typename T, typename U>
concept reserveable = requires(T& x, U sz) { x.reserve(sz); };

//! Reserve all the vectors/arrays using .reserve(sz)
template <typename SZ, reserveable<SZ> T, reserveable<SZ>... Rest>
void reserve(const SZ& sz, T& x, Rest&... rest) {
  x.reserve(sz);

  if constexpr (sizeof...(rest) > 0) {
    reserve(sz, rest...);
  }
}

template <typename T, typename U>
concept elemwise_resizeable = requires(T& x, U sz) { (*x.begin()).resize(sz); };

//! Resize all the sub-vectors/sub-arrays using .resize(sz)
template <typename SZ,
          elemwise_resizeable<SZ> T,
          elemwise_resizeable<SZ>... Rest>
void elemwise_resize(const SZ& sz, T& x, Rest&... rest) {
  for (auto& v : x) v.resize(sz);

  if constexpr (sizeof...(rest) > 0) {
    elemwise_resize(sz, rest...);
  }
}
}  // namespace arr
