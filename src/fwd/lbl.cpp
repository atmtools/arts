#include "lbl.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>
#include <variant>

template <Index... ints>
std::vector<lbl::band_models> all_models(
    Numeric t,
    Numeric p,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ArrayOfArrayOfSpeciesTag& allspecs,
    const Vector& allvmrs,
    const ArrayOfArrayOfAbsorptionLines& specbands,
    std::integer_sequence<Index, ints...>) {
  std::vector<lbl::band_models> out{
      std::variant_alternative_t<ints, lbl::band_models>(
          t, p, isotopologue_ratios, allspecs, allvmrs, specbands)...};
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

lbl::full::full(Numeric t,
                Numeric p,
                const SpeciesIsotopologueRatios& isotopologue_ratios,
                const ArrayOfArrayOfSpeciesTag& allspecs,
                const Vector& allvmrs,
                const ArrayOfArrayOfAbsorptionLines& specbands)
    : models(all_models(t,
                        p,
                        isotopologue_ratios,
                        allspecs,
                        allvmrs,
                        specbands,
                        std::make_integer_sequence<
                            Index,
                            std::variant_size_v<lbl::band_models>>{})) {
  ARTS_USER_ERROR_IF(size() not_eq static_cast<std::size_t>(nelem(specbands)),
                     "Size mismatch between specbands and models");
}

Complex lbl::full::at(Numeric f) const {
  return std::transform_reduce(models.begin(),
                               models.end(),
                               Complex{},
                               std::plus<>{},
                               [f](const auto& m) {
                                 return std::visit(
                                     [f](auto& mod) { return mod.at(f); }, m);
                               });
}

std::size_t lbl::full::size() const {
  return std::transform_reduce(
      models.begin(),
      models.end(),
      std::size_t{},
      std::plus<>{},
      [](const auto& m) {
        return std::visit([](const auto& mod) { return mod.size(); }, m);
      });
}

ComplexVector lbl::full::at(const Vector& f) const {
  ComplexVector out(f.size());
  at(out, f);
  return out;
}

void lbl::full::at(ExhaustiveComplexVectorView out, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), out.begin(), [this](const auto& f) {
    return at(f);
  });
}

ComplexVector lbl::full::at_par(const Vector& f) const {
  ComplexVector out(f.size());
  at_par(out, f);
  return out;
}

void lbl::full::at_par(ExhaustiveComplexVectorView out, const Vector& fs) const {
  const Index n = fs.size();
  #pragma omp parallel for
  for (Index i=0; i<n; ++i) {
    out[i] = at(fs[i]);
  }
}

static_assert(lbl::bandable<lbl::full>,
              "lbl::mtckd::single is not representative of a band of lines");
