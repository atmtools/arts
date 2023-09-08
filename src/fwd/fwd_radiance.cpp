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
spectral_radiance::spectral_radiance(
    const Vector& z,
    const std::vector<AtmPoint>& atm_points_,
    const ArrayOfArrayOfSpeciesTag& allspecs,
    const PredefinedModelData& predef_data,
    const ArrayOfCIARecord& cia_data,
    const ArrayOfXsecRecord& xsec_data,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ArrayOfArrayOfAbsorptionLines& lbl_data,
    Numeric cia_extrap,
    Index cia_robust)
    : refl(0), altitude(z.begin(), z.end()), atm_points(atm_points_.begin(), atm_points_.end()) {
  const std::vector<std::shared_ptr<CIARecord>> cia_records = [&]() {
    std::vector<std::shared_ptr<CIARecord>> ciar;
    ciar.reserve(cia_data.size());
    for (const auto& cia : cia_data) {
      ciar.push_back(std::make_shared<CIARecord>(cia));
    }
    return ciar;
  }();
  const std::vector<std::shared_ptr<XsecRecord>> xsec_records = [&]() {
    std::vector<std::shared_ptr<XsecRecord>> xsecr;
    xsecr.reserve(cia_data.size());
    for (const auto& xsec : xsec_data) {
      xsecr.push_back(std::make_shared<XsecRecord>(xsec));
    }
    return xsecr;
  }();
  const std::shared_ptr<PredefinedModelData> predef_model_data =
      std::make_shared<PredefinedModelData>(predef_data);

  const std::size_t n = altitude.size();
  const std::size_t m = allspecs.size();
  
  ARTS_USER_ERROR_IF(n == 0, "Must have some atmospheric size")

  ARTS_USER_ERROR_IF(n not_eq atm_points.size(),
                     "Size mismatch: altitude and vmr (",
                     n,
                     " vs ",
                     atm_points.size(),
                     ')')
                     
  ARTS_USER_ERROR_IF(m not_eq lbl_data.size(),
                     "Size mismatch: species and absorption lines data")

  ARTS_USER_ERROR_IF(not std::is_sorted(altitude.begin(), altitude.end()),
                     "Altitude must be sorted in ascending order")

  String error_msg;
  models.resize(n);
#pragma omp parallel for
  for (std::size_t i = 0; i < n; ++i) {
    if (error_msg.size()) continue;
    try {
      models[i] = full_absorption(atm_points[i],
                                  allspecs,
                                  predef_model_data,
                                  cia_records,
                                  xsec_records,
                                  isotopologue_ratios,
                                  lbl_data,
                                  cia_extrap,
                                  cia_robust);
    } catch (std::exception& e) {
#pragma omp critical
      if (error_msg.size() == 0) error_msg = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error_msg.size(), error_msg, '\n')
}

ExhaustiveVectorView spectral_radiance::planar(ExhaustiveVectorView rad,
                                               Numeric f,
                                               Numeric za) const try {
  using Conversion::cosd;
  using std::abs;
  using std::exp;
  using std::lerp;
  using std::midpoint;
  using std::size_t;
  constexpr Numeric Tcmb = Constant::cosmic_microwave_background_temperature;

  const std::size_t n = altitude.size();
  ARTS_USER_ERROR_IF(n not_eq static_cast<size_t>(rad.nelem()),
                     "Bad size absorption input vector view'\n")
  ARTS_USER_ERROR_IF(
      za == 90.0,
      "You cannot look sideways in a plane-parallel atmosphere.\n"
      "The zenith angle must be above or below 1 degree of the limb.\n")

  const Numeric z_scl = 1.0 / abs(cosd(za));

  const bool looking_down = za > 90;
  if (looking_down) {
    const Numeric Bbg = refl == 1.0 ? 0.0 : planck(f, atm_points.front().temperature);
    const Numeric Ibg = refl == 0.0 ? 0.0 : planar(rad, f, 180 - za)[0];
    rad[0] = lerp(Bbg, Ibg, refl);

    Numeric abs_past = models.front().at(f).real();
    Numeric B_past = planck(f, atm_points.front().temperature);
    for (size_t i = 1; i < n; ++i) {
      const Numeric abs_this = models[i].at(f).real();
      const Numeric B_this = planck(f, atm_points[i].temperature);
      const Numeric B = midpoint(B_past, B_this);
      const Numeric T = exp(-z_scl * midpoint(abs_this, abs_past) *
                            (altitude[i] - altitude[i - 1]));

      rad[i] = lerp(B, rad[i - 1], T);
      abs_past = abs_this;
      B_past = B_this;
    }
  } else {
    rad[n - 1] = planck(f, Tcmb);

    Numeric abs_past = models.back().at(f).real();
    Numeric B_past = planck(f, atm_points.back().temperature);
    for (size_t i = n - 2; i < n; --i) {  // wraps around
      const Numeric abs_this = models[i].at(f).real();
      const Numeric B_this = planck(f, atm_points[i].temperature);
      const Numeric B = midpoint(B_past, B_this);
      const Numeric T = exp(-z_scl * midpoint(abs_this, abs_past) *
                            (altitude[i + 1] - altitude[i]));

      rad[i] = lerp(B, rad[i + 1], T);
      abs_past = abs_this;
      B_past = B_this;
    }
  }

  return rad;
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("Input: f=", f, "; za=", za, '\n', e.what()));
}

Vector spectral_radiance::planar(Numeric f, Numeric za) const {
  Vector rad(altitude.size());
  planar(rad, f, za);
  return rad;
}

void spectral_radiance::planar_par(ExhaustiveMatrixView rad,
                                   const Vector& fs,
                                   Numeric za) const {
  const Index n = fs.size();

#pragma omp parallel for
  for (Index i = 0; i < n; i++) {
    planar(rad[i], fs[i], za);
  }
}

Matrix spectral_radiance::planar_par(const Vector& fs, Numeric za) const {
  Matrix rad(fs.size(), altitude.size());
  planar_par(rad, fs, za);
  return rad;
}

std::ostream& operator<<(std::ostream& os, const spectral_radiance&) {
  return os << "No useful output operator yet";
}
}  // namespace fwd::profile
