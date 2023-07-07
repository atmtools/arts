#include "lbl_mtckd_voigt.h"

#include <algorithm>
#include <iterator>
#include <numeric>

#include "absorptionlines.h"
#include "arts_constants.h"
#include "arts_conversions.h"
#include "lbl_algorithms.h"
#include "linescaling.h"
#include "lineshapemodel.h"
#include "physics_funcs.h"

namespace lbl::mtckd {
single::single(Numeric T,
               Numeric P,
               const SpeciesIsotopologueRatios& isotopologue_ratios,
               const ArrayOfArrayOfSpeciesTag& allspecs,
               const Vector& allvmrs,
               const AbsorptionLines& band,
               Index line) {
  using Constant::inv_sqrt_pi;
  using Conversion::hz2joule;
  using Conversion::kelvin2joule;

  const Numeric QT = single_partition_function(T, band.Isotopologue());
  const Numeric QT0 = single_partition_function(band.T0, band.Isotopologue());

  const auto& ln = band.lines[line];

  const auto X = band.ShapeParameters(
      line, T, P, band.BroadeningSpeciesVMR(allvmrs, allspecs));

  F0 = ln.F0 + X.D0;
  invGD = 1.0 / (band.DopplerConstant(T) * F0);
  z_imag = invGD * X.G0;

  const Numeric selfvmr = band.SelfVMR(allvmrs, allspecs);
  const Numeric isotrat =
      isotopologue_ratios[band.quantumidentity.Isotopologue()];

  scl = -inv_sqrt_pi * invGD * QT0 * isotrat * selfvmr * ln.I0 *
        boltzman_ratio(T, band.T0, ln.E0) /
        (F0 * std::expm1(-hz2joule(F0) / kelvin2joule(band.T0)) * QT);

  cutoff = at<true>(F0 + cutoff_range::limit);
}

bool single::is_valid(const AbsorptionLines& band) {
  return band.lineshapetype == LineShapeType::VP and
         band.cutoff == AbsorptionCutoffType::ByLine and
         band.cutofffreq == cutoff_range::limit and
         band.population == AbsorptionPopulationType::LTE and
         band.normalization == AbsorptionNormalizationType::SFS and
         not band.AnyLinemixing();
}

single_lm::single_lm(Numeric T,
                     Numeric P,
                     const SpeciesIsotopologueRatios& isotopologue_ratios,
                     const ArrayOfArrayOfSpeciesTag& allspecs,
                     const Vector& allvmrs,
                     const AbsorptionLines& band,
                     Index line) {
  using Constant::inv_sqrt_pi;
  using Conversion::hz2joule;
  using Conversion::kelvin2joule;

  const Numeric QT = single_partition_function(T, band.Isotopologue());
  const Numeric QT0 = single_partition_function(band.T0, band.Isotopologue());

  const auto& ln = band.lines[line];

  const auto X = band.ShapeParameters(
      line, T, P, band.BroadeningSpeciesVMR(allvmrs, allspecs));

  F0 = ln.F0 + X.D0 + X.DV;
  invGD = 1.0 / (band.DopplerConstant(T) * F0);
  z_imag = invGD * X.G0;

  const Numeric selfvmr = band.SelfVMR(allvmrs, allspecs);
  const Numeric isotrat =
      isotopologue_ratios[band.quantumidentity.Isotopologue()];

  scl = (-inv_sqrt_pi * invGD * QT0 * isotrat * selfvmr * ln.I0 *
         boltzman_ratio(T, band.T0, ln.E0) /
         (F0 * std::expm1(-hz2joule(F0) / kelvin2joule(band.T0)) * QT)) *
        Complex{1 + X.G, -X.Y};

  cutupp = at<true>(F0 + cutoff_range::limit);
  cutlow = at<true>(F0 - cutoff_range::limit);
}

bool single_lm::is_valid(const AbsorptionLines& band) {
  return band.lineshapetype == LineShapeType::VP and
         band.cutoff == AbsorptionCutoffType::ByLine and
         band.cutofffreq == cutoff_range::limit and
         band.population == AbsorptionPopulationType::LTE and
         band.normalization == AbsorptionNormalizationType::SFS and
         band.AnyLinemixing();
}

std::size_t band::validity_count(
    const ArrayOfArrayOfAbsorptionLines& specbands) {
  std::size_t n = 0;
  for (auto& bands : specbands) {
    for (auto& band : bands) {
      if (single::is_valid(band)) n += band.lines.size();
    }
  }
  return n;
}

band::band(Numeric t,
           Numeric p,
           const SpeciesIsotopologueRatios& isotopologue_ratios,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const Vector& allvmrs,
           const ArrayOfArrayOfAbsorptionLines& specbands)
    : T(t), P(p) {
  lines.reserve(validity_count(specbands));

  for (auto& bs : specbands) {
    for (auto& b : bs) {
      if (single::is_valid(b)) {
        for (std::size_t line = 0; line < b.lines.size(); ++line) {
          lines.emplace_back(
              t, P, isotopologue_ratios, allspecs, allvmrs, b, line);
        }
      }
    }
  }

  std::sort(lines.begin(), lines.end(), [](const auto& l1, const auto& l2) {
    return l1.F0 < l2.F0;
  });
}

Complex band::at(Numeric f) const {
  using Conversion::hz2joule;
  using Conversion::kelvin2joule;

  const Numeric fscl = -f * std::expm1(-hz2joule(f) / kelvin2joule(T));
  const Numeric nscl = number_density(P, T);

  const Complex sum = sumup<cutoff_range>(lines, f);
  return sum.real() < 0 ? Complex{0, 0} : sum * fscl * nscl;
}

void band::at(ComplexVector& out, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), out.begin(), [this](const auto& f) {
    return at(f);
  });
}

ComplexVector band::at(const Vector& f) const {
  ComplexVector out(f.size());
  at(out, f);
  return out;
}

std::size_t band_lm::validity_count(
    const ArrayOfArrayOfAbsorptionLines& specbands) {
  std::size_t n = 0;
  for (auto& bands : specbands) {
    for (auto& band : bands) {
      if (single_lm::is_valid(band)) {
        n += 1;
        break;
      }
    }
  }
  return n;
}

std::size_t band_lm::size() const {
  return std::transform_reduce(bands.begin(),
                               bands.end(),
                               std::size_t{},
                               std::plus<>{},
                               [](const auto& b) { return b.size(); });
}

band_lm::band_lm(Numeric t,
                 Numeric p,
                 const SpeciesIsotopologueRatios& isotopologue_ratios,
                 const ArrayOfArrayOfSpeciesTag& allspecs,
                 const Vector& allvmrs,
                 const ArrayOfArrayOfAbsorptionLines& specbands)
    : T(t), P(p) {
  bands.reserve(validity_count(specbands));
  for (auto& bs : specbands) {
    for (auto& b : bs) {
      if (single_lm::is_valid(b)) {
        bands.emplace_back();
        bands.back().reserve(b.lines.size());
        for (std::size_t line = 0; line < b.lines.size(); ++line) {
          bands.back().emplace_back(
              t, P, isotopologue_ratios, allspecs, allvmrs, b, line);
        }
      }
    }
  }

  for (auto& lines : bands) {
    std::sort(lines.begin(), lines.end(), [](const auto& l1, const auto& l2) {
      return l1.F0 < l2.F0;
    });
  }
}

Complex band_lm::at(Numeric f) const {
  using Conversion::hz2joule;
  using Conversion::kelvin2joule;

  const Numeric fscl = -f * std::expm1(-hz2joule(f) / kelvin2joule(T));
  const Numeric nscl = number_density(P, T);

  const auto allsum = [f](auto& band) {
    const Complex sum = sumup<cutoff_range>(band, f);
    return sum.real() < 0 ? Complex{0, 0} : sum;
  };

  const auto sum = std::transform_reduce(
      bands.begin(), bands.end(), Complex{0, 0}, std::plus<>{}, allsum);

  return fscl * nscl * sum;
}

void band_lm::at(ComplexVector& out, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), out.begin(), [this](const auto& f) {
    return at(f);
  });
}

ComplexVector band_lm::at(const Vector& f) const {
  ComplexVector out(f.size());
  at(out, f);
  return out;
}

std::size_t full::size() const { return bands.size() + bands_lm.size(); }

Complex full::at(Numeric f) const { return bands.at(f) + bands_lm.at(f); }

ComplexVector full::at(const Vector& f) const {
  ComplexVector out(f.size());
  at(out, f);
  return out;
}

void full::at(ComplexVector& out, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), out.begin(), [this](const auto& f) {
    return at(f);
  });
}
}  // namespace lbl::mtckd

static_assert(lbl::singleable<lbl::mtckd::single>,
              "lbl::mtckd::single is not representative of a single line");

static_assert(lbl::bandable<lbl::mtckd::band>,
              "lbl::mtckd::single is not representative of a band of lines");

static_assert(lbl::singleable<lbl::mtckd::single_lm>,
              "lbl::mtckd::single is not representative of a single line");

static_assert(lbl::bandable<lbl::mtckd::band_lm>,
              "lbl::mtckd::single is not representative of a band of lines");

static_assert(lbl::bandable<lbl::mtckd::full>,
              "lbl::mtckd::single is not representative of a band of lines");
