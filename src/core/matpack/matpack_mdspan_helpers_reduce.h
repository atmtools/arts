#pragma once

#include <compare.h>
#include <nonstd.h>

#include <algorithm>
#include <cmath>

#include "matpack_mdspan_common_types.h"

namespace matpack {
/** Sum all elements in the range using it's reduce operation
 * 
 * @param self The range to sum
 * @return The sum of all elements in the range
 */
template <any_md Self>
constexpr auto sum(const Self& self) {
  return std::reduce(self.elem_begin(), self.elem_end());
}

/** Sum the absolute of all elements in the range
 * 
 * @param self The range to sum
 * @return The absolute sum of all elements in the range
 */
template <any_md Self>
constexpr auto abssum(const Self& self) {
  const auto sum_abs = [](auto& a) -> value_type<Self> {
    return nonstd::abs(a);
  };

  return std::transform_reduce(self.elem_begin(),
                               self.elem_end(),
                               value_type<Self>{},
                               std::plus<>{},
                               sum_abs);
}

/** The mean of all elements in the range using it's reduce operation after scaling the value by the size of the range
 * 
 * This is effectively the same as sum(self / len(self)).
 *
 * Note that this is numerically different than fastmean(self), which is the faster version that
 * assumes that there are no numerical issues doing the scaling only once.
 * 
 * @param self The range to sum
 * @return The sum of all elements in the range
 */
template <any_md Self>
constexpr auto mean(const Self& self) {
  const auto scale = [n = static_cast<Numeric>(self.size())](auto& a) {
    return a / n;
  };
  return std::transform_reduce(self.elem_begin(),
                               self.elem_end(),
                               value_type<Self>{},
                               std::plus<>{},
                               scale);
}

/** The mean of all elements in the range using it's reduce operation after scaling the value by the size of the range
 * 
 * This is effectively the same as sum(self) / len(self).
 *
 * Note that this is numerically different than mean(self), which performs the scaling for each element.
 * The difference is that this version is faster, but may have numerical issues if the range is very large.
 * 
 * @param self The range to sum
 * @return The sum of all elements in the range
 */
template <any_md Self>
constexpr auto fastmean(const Self& self) {
  const auto len = static_cast<Numeric>(self.size());
  return sum(self) / len;
}

/** Get the minimum value in the range
 * 
 * @param self The range to find the minimum value in
 * @return The smallest value in the range or the maximum value of the type if the range is empty
 */
template <any_md Self>
constexpr auto min(const Self& self) {
  return stdr::min(elemwise_range(self));
}

/** Get the maximum value in the range
 * 
 * @param self The range to find the maximum value in
 * @return The largest value in the range or the minimum value of the type if the range is empty
 */
template <any_md Self>
constexpr auto max(const Self& self) {
  return stdr::max(elemwise_range(self));
}

/** Get the maximum value in the range
 * 
 * @param self The range to find the maximum value in
 * @return The largest value in the range or the minimum value of the type if the range is empty
 */
template <any_md Self>
constexpr auto minmax(const Self& self) {
  auto e      = elemwise_range(self);
  auto [a, b] = stdr::minmax_element(e);
  assert(a != stdr::end(e));
  assert(b != stdr::end(e));
  return std::pair{*a, *b};
}

/** Get the number of non-NaN elements in the range
 * 
 * @param self The range to count the number of non-NaN elements in
 * @return The number of elements that are not NaN
 */
template <any_md Self>
constexpr Size nansize(const Self& self) {
  Size out = 0;

  for (auto& v : elemwise_range(self)) {
    if (not nonstd::isnan(v)) out++;
  }

  return out;
}

/** Get the sum of all elements in the range that are not NaN
 * 
 * @param self The range to sum
 * @return The sum of all elements in the range that are not NaN
 */
template <any_md Self>
constexpr auto nansum(const Self& self) {
  const auto sum_nonnan = [](auto& a) -> value_type<Self> {
    if (nonstd::isnan(a)) return {};
    return a;
  };

  return std::transform_reduce(self.elem_begin(),
                               self.elem_end(),
                               value_type<Self>{},
                               std::plus<>{},
                               sum_nonnan);
}

/** Find the mean of all elements in the range that are not NaN
 * 
 * This is effectively the same as nansum(self / nansize(self)).
 * 
 * See explanation in mean(self) for the difference between this and nanfastmean(self).
 * 
 * @param self The range to find the mean of
 * @return The mean of all elements in the range that are not NaN
 */
template <any_md Self>
constexpr auto nanmean(const Self& self) {
  const auto scale_nonnan =
      [n = static_cast<Numeric>(nansize(self))](auto& a) -> value_type<Self> {
    if (nonstd::isnan(a)) return {};
    return a / n;
  };

  return std::transform_reduce(self.elem_begin(),
                               self.elem_end(),
                               value_type<Self>{},
                               std::plus<>{},
                               scale_nonnan);
}

/** Find the mean of all elements in the range that are not NaN
 * 
 * This is effectively the same as nansum(self) / nansize(self).
 * 
 * See explanation in fastmean(self) for the difference between this and nanmean(self).
 * 
 * @param self The range to find the mean of
 * @return The mean of all elements in the range that are not NaN
 */
template <any_md Self>
constexpr auto nanfastmean(const Self& self) {
  auto n = static_cast<Numeric>(nansize(self));
  return nansum(self) / n;
}

/** Find the minimum value in the range that is not NaN
 * 
 * @param self The range to find the minimum value in
 * @return The smallest value in the range filtering out nan
 */
template <any_md Self>
constexpr auto nanmin(const Self& self) {
  const auto not_isnan = [](auto& a) { return not nonstd::isnan(a); };
  return stdr::min(elemwise_range(self) | stdv::filter(not_isnan));
}

/** Find the maximum value in the range that is not NaN
 * 
 * @param self The range to find the maximum value in
 * @return The largest value in the range filtering out nan
 */
template <any_md Self>
constexpr auto nanmax(const Self& self) {
  const auto not_isnan = [](auto& a) { return not nonstd::isnan(a); };
  return stdr::max(elemwise_range(self) | stdv::filter(not_isnan));
}

/** Get the sum of the elementwise product of two ranges
 * 
 * This is effectively the dot-product of two vectors, thus the name.
 * 
 * @param self One of the ranges
 * @param other The other range
 * @param init The initial value, defaults to 0 for most common types.  Beware the type.
 * @return The dot-product of the two ranges
 */
template <any_md Self, any_md Other, class T = value_type<Self>>
constexpr T dot(const Self& self, const Other& other, T init = {}) {
  return std::transform_reduce(
      self.elem_begin(), self.elem_end(), other.elem_begin(), init);
}

/** Get the sum of the elementwise product of two ranges
 * 
 * This is effectively the square root of the dot-product of two vectors, thus the name.
 * 
 * @param self One of the ranges
 * @param other The other range
 * @return The hypothenuse of the range
 */
template <any_md Self>
constexpr auto hypot(const Self& self) {
  return std::sqrt(dot(self, self));
}

template <any_md Self, any_md Other>
constexpr bool operator==(const Self& self, const Other& other) {
  return stdr::equal(elemwise_range(self), elemwise_range(other));
}

template <any_md Self, any_md Other>
constexpr bool operator!=(const Self& self, const Other& other) {
  return not stdr::equal(elemwise_range(self), elemwise_range(other));
}

template <any_md Self>
constexpr bool operator==(const Self& self, const value_type<Self>& other) {
  for (auto& v : elemwise_range(self)) {
    if (v != other) return false;
  }
  return true;
}

template <any_md Self>
constexpr bool operator!=(const Self& self, const value_type<Self>& other) {
  for (auto& v : elemwise_range(self)) {
    if (v == other) return false;
  }
  return true;
}

template <any_md Self>
constexpr bool any_negative(const Self& self) {
  return stdr::any_of(elemwise_range(self), Cmp::lt(value_type<Self>{}));
}

template <any_md Self>
constexpr bool any_nan(const Self& self) {
  const auto isnan = [](auto& a) { return nonstd::isnan(a); };
  return stdr::any_of(elemwise_range(self), isnan);
}
}  // namespace matpack
