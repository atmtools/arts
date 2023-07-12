#include "fwd_radiance.h"

#include <predefined/predef_data.h>

#include <cstddef>
#include <exception>
#include <memory>
#include <stdexcept>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "cia.h"
#include "debug.h"
#include "matpack_concepts.h"
#include "physics_funcs.h"

namespace fwd::profile {
radiance::radiance(const Vector& z,
                   const Vector& p,
                   const Vector& t,
                   const std::vector<Vector>& allvmrs,
                   const ArrayOfArrayOfSpeciesTag& allspecs,
                   const PredefinedModelData& predef_data,
                   const ArrayOfCIARecord& cia_data,
                   const SpeciesIsotopologueRatios& isotopologue_ratios,
                   const ArrayOfArrayOfAbsorptionLines& lbl_data,
                   Numeric cia_extrap,
                   Index cia_robust,
                   Verbosity cia_verb)
    : refl(0), altitude(z.begin(), z.end()), temperature(t.begin(), t.end()) {
  const std::vector<std::shared_ptr<CIARecord>> cia_records = [&]() {
    std::vector<std::shared_ptr<CIARecord>> ciar;
    ciar.reserve(cia_data.size());
    for (const auto& cia : cia_data) {
      ciar.push_back(std::make_shared<CIARecord>(cia));
    }
    return ciar;
  }();
  const std::shared_ptr<PredefinedModelData> predef_model_data =
      std::make_shared<PredefinedModelData>(predef_data);

  const std::size_t n = altitude.size();
  const std::size_t m = allspecs.size();
  ARTS_USER_ERROR_IF(n not_eq static_cast<std::size_t>(p.size()),
                     "Size mismatch: altitude and pressure (",
                     n,
                     " vs ",
                     p.size(),
                     ')')
  ARTS_USER_ERROR_IF(n not_eq temperature.size(),
                     "Size mismatch: altitude and temperature (",
                     n,
                     " vs ",
                     t.size(),
                     ')')
  ARTS_USER_ERROR_IF(n not_eq allvmrs.size(),
                     "Size mismatch: altitude and vmr (",
                     n,
                     " vs ",
                     allvmrs.size(),
                     ')')
  ARTS_USER_ERROR_IF(
      std::any_of(allvmrs.begin(),
                  allvmrs.end(),
                  [m](auto& vmr) {
                    return static_cast<std::size_t>(vmr.size()) not_eq m;
                  }),
      "Size mismatch: species and vmr")
  ARTS_USER_ERROR_IF(m not_eq lbl_data.size(),
                     "Size mismatch: species and absorption lines data")
  ARTS_USER_ERROR_IF(not std::is_sorted(altitude.begin(), altitude.end()),
                     "Altitude must be sorted in ascending order")
  ARTS_USER_ERROR_IF(n == 0, "Must have some atmospheric size: [\n")

  models.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    models.emplace_back(p[i],
                        t[i],
                        allvmrs[i],
                        allspecs,
                        predef_model_data,
                        cia_records,
                        isotopologue_ratios,
                        lbl_data,
                        cia_extrap,
                        cia_robust,
                        cia_verb);
  }
}

ExhaustiveVectorView radiance::planar(ExhaustiveVectorView rad,
                                      Numeric f,
                                      Numeric za) const try {
  const std::size_t n = altitude.size();
  ARTS_USER_ERROR_IF(n not_eq static_cast<std::size_t>(rad.nelem()),
                     "Bad size absorption input vector view'\n")
  ARTS_USER_ERROR_IF(
      za > 89 and za < 91,
      "You cannot look sideways in a plane-parallel atmosphere.\n"
      "The zenith angle must be above or below 1 degree of the limb.\n")

  const Numeric z_scl = 1.0 / std::abs(Conversion::cosd(za));

  const bool looking_down = za > 90;
  if (looking_down) {
    const Numeric Bbg = refl == 1.0 ? 0.0 : planck(f, temperature[0]);
    const Numeric Ibg = refl == 0.0 ? 0.0 : planar(rad, f, 180 - za)[0];

    rad[0] = std::lerp(Bbg, Ibg, refl);

    Numeric abs_past = models.front().at(f).real();

    for (std::size_t i = 1; i < n; ++i) {
      const Numeric abs_this = models[i].at(f).real();
      const Numeric B = planck(f, temperature[i - 1]);
      const Numeric T = std::exp(-z_scl * std::midpoint(abs_this, abs_past) *
                                 (altitude[i] - altitude[i - 1]));

      rad[i] = std::lerp(B, rad[i - 1], T);
      abs_past = abs_this;
    }
  } else {
    rad[n - 1] = planck(f, Constant::cosmic_microwave_background_temperature);

    Numeric abs_past = models.front().at(f).real();
    for (std::size_t i = n - 2; i < n; --i) {  // wraps around
      const Numeric abs_this = models[i].at(f).real();
      const Numeric B = planck(f, temperature[i + 1]);
      const Numeric T = std::exp(-z_scl * std::midpoint(abs_this, abs_past) *
                                 (altitude[i + 1] - altitude[i]));

      rad[i] = std::lerp(B, rad[i + 1], T);
      abs_past = abs_this;
    }
  }

  return rad;
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Input: f=", f, "; za=", za, '\n', e.what()));
}

Vector radiance::planar(Numeric f, Numeric za) const {
  Vector rad(altitude.size());
  planar(rad, f, za);
  return rad;
}

void radiance::planar_par(ExhaustiveMatrixView rad,
                          const Vector& fs,
                          Numeric za) const {
  const Index n = fs.size();

#pragma omp parallel for
  for (Index i = 0; i < n; i++) {
    planar(rad[i], fs[i], za);
  }
}

Matrix radiance::planar_par(const Vector& fs, Numeric za) const {
  Matrix rad(fs.size(), altitude.size());
  planar_par(rad, fs, za);
  return rad;
}

std::ostream& operator<<(std::ostream& os, const radiance&) {
  return os << "No useful output operator yet";
}
}  // namespace fwd::profile
