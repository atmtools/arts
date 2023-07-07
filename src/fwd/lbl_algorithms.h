#pragma once

#include "lbl_concepts.h"

namespace lbl {
struct no_cutoff_range {
  static constexpr Numeric limit = std::numeric_limits<Numeric>::infinity();
};

template <cutoffable cutoff = no_cutoff_range>
[[nodiscard]] constexpr Complex sumup(const list_singleable auto& lines,
                                      Numeric f) {
  if constexpr (has_cutoff(cutoff{})) {
    const auto first = std::lower_bound(
        lines.begin(),
        lines.end(),
        f - 750e9,
        [](const auto& l, const Numeric fm) { return l.F0 < fm; });
    const auto last = std::upper_bound(
        first, lines.end(), f + 750e9, [](const Numeric fm, const auto& l) {
          return fm < l.F0;
        });

    return std::transform_reduce(
        first, last, Complex{0, 0}, std::plus<>{}, [f](const auto& l) {
          return l.at(f);
        });
  } else {
    return std::transform_reduce(lines.begin(),
                                 lines.end(),
                                 Complex{0, 0},
                                 std::plus<>(),
                                 [f](const auto& l) { return l.at(f); });
  }
}
}  // namespace lbl
