#pragma once

#include <Faddeeva/Faddeeva.hh>
#include <algorithm>
#include <cmath>
#include <ratio>
#include <vector>

#include "absorptionlines.h"
#include "lbl_concepts.h"
#include "lineshapemodel.h"
#include "matpack_concepts.h"
#include "species_tags.h"

namespace fwd::lbl::mtckd {
static constexpr Numeric cutoff_freq = 750e9;

struct single {
  Numeric scl{};
  Numeric invGD{};
  Numeric F0{};
  Numeric z_imag{};
  Complex cutoff{};

  single(Numeric T,
         Numeric P,
         const SpeciesIsotopologueRatios& isotopologue_ratios,
         const ArrayOfArrayOfSpeciesTag& allspecs,
         const Vector& allvmrs,
         const AbsorptionLines& band,
         Index line);

  // This times numdens * f * (1 - exp(hf / kT)) is the absorption coeff
  template <bool no_cutoff = false>
  [[nodiscard]] Complex at(Numeric f) const {
    if constexpr (no_cutoff) {
      return scl * Faddeeva::w(Complex{invGD * (f - F0), z_imag});
    } else {
      return scl * Faddeeva::w(Complex{invGD * (f - F0), z_imag}) - cutoff;
    }
  }

  [[nodiscard]] static bool is_valid(const AbsorptionLines& band);
};

struct single_lm {
  Complex scl{};
  Numeric invGD{};
  Numeric F0{};
  Numeric z_imag{};
  Complex cutupp{};
  Complex cutlow{};

  single_lm(Numeric T,
            Numeric P,
            const SpeciesIsotopologueRatios& isotopologue_ratios,
            const ArrayOfArrayOfSpeciesTag& allspecs,
            const Vector& allvmrs,
            const AbsorptionLines& band,
            Index line);

  [[nodiscard]] constexpr Complex cutoff(Numeric f) const {
    f = std::clamp<Numeric>(0.5 + (f - F0) / (2 * cutoff_freq), 0.0, 1.0);
    return {std::lerp(cutlow.real(), cutupp.real(), f),
            std::lerp(cutlow.imag(), cutupp.imag(), f)};
  }

  // This times numdens * f * (1 - exp(hf / kT)) is the absorption coeff
  template <bool no_cutoff = false>
  [[nodiscard]] Complex at(Numeric f) const {
    if constexpr (no_cutoff) {
      return scl * Faddeeva::w(Complex{invGD * (f - F0), z_imag});
    } else {
      return scl * Faddeeva::w(Complex{invGD * (f - F0), z_imag}) - cutoff(f);
    }
  }

  [[nodiscard]] static bool is_valid(const AbsorptionLines& band);
};

struct band {
  std::vector<single> lines;
  Numeric T;
  Numeric P;

  band() = default;

  band(Numeric T,
       Numeric P,
       const SpeciesIsotopologueRatios& isotopologue_ratios,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const Vector& allvmrs,
       const ArrayOfArrayOfAbsorptionLines& specbands);

  [[nodiscard]] static std::size_t validity_count(
      const ArrayOfArrayOfAbsorptionLines& band);
  [[nodiscard]] std::size_t size() const { return lines.size(); }

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ExhaustiveComplexVectorView out, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};

struct band_lm {
  std::vector<std::vector<single_lm>> bands;
  Numeric T;
  Numeric P;

  band_lm() = default;

  band_lm(Numeric T,
          Numeric P,
          const SpeciesIsotopologueRatios& isotopologue_ratios,
          const ArrayOfArrayOfSpeciesTag& allspecs,
          const Vector& allvmrs,
          const ArrayOfArrayOfAbsorptionLines& specbands);

  [[nodiscard]] static std::size_t validity_count(
      const ArrayOfArrayOfAbsorptionLines& band);
  [[nodiscard]] std::size_t size() const;

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ExhaustiveComplexVectorView out, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};
}  // namespace fwd::lbl::mtckd
