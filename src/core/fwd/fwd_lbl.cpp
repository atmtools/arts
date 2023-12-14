#include "fwd_lbl.h"

#include <algorithm>
#include <numeric>
#include <utility>
#include <variant>

template <Index... ints>
std::vector<fwd::lbl::band_models> all_models(
    const AtmPoint& atm_point,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ArrayOfArrayOfAbsorptionLines& specbands,
    std::integer_sequence<Index, ints...>) {
  std::vector<fwd::lbl::band_models> out{
      std::variant_alternative_t<ints, fwd::lbl::band_models>(
          atm_point, isotopologue_ratios, specbands)...};
  out.erase(
      std::remove_if(out.begin(),
                     out.end(),
                     [](auto& m) {
                       return 0 == std::visit(
                                       [](auto& mod) { return mod.size(); }, m);
                     }),
      out.end());
  return out;
}

fwd::lbl::full::full(const AtmPoint& atm_point,
                     const SpeciesIsotopologueRatios& isotopologue_ratios,
                     const ArrayOfArrayOfAbsorptionLines& specbands)
    : models(all_models(atm_point,
                        isotopologue_ratios,
                        specbands,
                        std::make_integer_sequence<
                            Index,
                            std::variant_size_v<lbl::band_models>>{})) {
  ARTS_USER_ERROR_IF(
      size() not_eq static_cast<std::size_t>(Absorption::size(specbands)),
      "Size mismatch between specbands and models");
}

Complex fwd::lbl::full::at(Numeric f) const {
  return std::transform_reduce(models.begin(),
                               models.end(),
                               Complex{},
                               std::plus<>{},
                               [f](const auto& m) {
                                 return std::visit(
                                     [f](auto& mod) { return mod.at(f); }, m);
                               });
}

std::size_t fwd::lbl::full::size() const {
  return std::transform_reduce(
      models.begin(),
      models.end(),
      std::size_t{},
      std::plus<>{},
      [](const auto& m) {
        return std::visit([](const auto& mod) { return mod.size(); }, m);
      });
}

ComplexVector fwd::lbl::full::at(const Vector& f) const {
  ComplexVector out(f.size());
  at(out, f);
  return out;
}

void fwd::lbl::full::at(ExhaustiveComplexVectorView out,
                        const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), out.begin(), [this](const auto& f) {
    return at(f);
  });
}

ComplexVector fwd::lbl::full::at_par(const Vector& f) const {
  ComplexVector out(f.size());
  at_par(out, f);
  return out;
}

void fwd::lbl::full::at_par(ExhaustiveComplexVectorView out,
                            const Vector& fs) const {
  const Index n = fs.size();
#pragma omp parallel for
  for (Index i = 0; i < n; ++i) {
    out[i] = at(fs[i]);
  }
}
