#pragma once

#include <functional>

#include "fwd_lbl_mtckd_voigt.h"
#include "matpack_concepts.h"

namespace fwd::lbl {
using band_models = std::variant<mtckd::band_lm, mtckd::band>;

struct full {
  std::vector<band_models> models{};

  full() = default;

  full(const AtmPoint& atm_point,
       const SpeciesIsotopologueRatios& isotopologue_ratios,
       const ArrayOfArrayOfAbsorptionLines& specbands);

  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ExhaustiveComplexVectorView out, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
  void at_par(ExhaustiveComplexVectorView out, const Vector& fs) const;
  [[nodiscard]] ComplexVector at_par(const Vector& fs) const;

  template <typename bandable_t>
  full& add(bandable_t&& model) {
    models.emplace_back(std::forward<bandable_t>(model));
    return *this;
  }
};
}  // namespace fwd::lbl
