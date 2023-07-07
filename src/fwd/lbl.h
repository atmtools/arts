#pragma once

#include <functional>

#include "fwd/lbl_concepts.h"
#include "lbl_mtckd_voigt.h"
#include "matpack_concepts.h"

namespace lbl {
class full {
  std::vector<std::variant<mtckd::band_lm, mtckd::band>> models{};

public:
  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ComplexVector& out, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;

  template <bandable bandable_t>
  full& add(bandable_t&& model) {
    models.emplace_back(std::forward<bandable_t>(model));
    return *this;
  }
};
}  // namespace lbl
