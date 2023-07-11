#pragma once

#include "lbl_concepts.h"

namespace fwd::lbl {
namespace internal {
/** Sums up the contribution of input lines [first, last)
 *
 * Internal function, not to be used directly, if you want to sum up the
 * contribution of a list of lines, use main sumup() functions
 * 
 * @param first Iterator to the first element
 * @param last Iterator to the last element
 * @param f A frequency
 * @return constexpr Complex Some sort of absorption coefficient
 */
constexpr Complex sumup(const std::forward_iterator auto& first,
                        const std::forward_iterator auto& last,
                        Numeric f)
  requires(singleable<decltype(*first)> and singleable<decltype(*last)>)
{
  return std::transform_reduce(
      first, last, Complex{0, 0}, std::plus<>{}, [f](const auto& l) {
        return l.at(f);
      });
}
}  // namespace internal

/** Sums up the contribution of input lines
 * 
 * @param lines List of absorption lines
 * @param f A frequency
 * @return constexpr Complex Some sort of absorption coefficient
 */
constexpr Complex sumup(const list_singleable auto& lines, Numeric f) {
  return internal::sumup(lines.begin(), lines.end(), f);
}

/** Sums up the contribution of input lines in range [f - fc, f + fc]
 * 
 * @param lines List of absorption lines
 * @param f A frequency
 * @param fc A cutoff frequency
 * @return constexpr Complex Some sort of absorption coefficient
 */
constexpr Complex sumup(const list_singleable auto& lines,
                        Numeric f,
                        Numeric fc) {
  const auto first = std::lower_bound(
      lines.begin(), lines.end(), f - fc, [](const auto& l, const Numeric fm) {
        return l.F0 < fm;
      });
  const auto last = std::upper_bound(
      first, lines.end(), f + fc, [](const Numeric fm, const auto& l) {
        return fm < l.F0;
      });
  return internal::sumup(first, last, f);
}
}  // namespace fwd::lbl
