#include "fwd_irradiance.h"

#include <algorithm>

#include "arts_constants.h"
#include "arts_conversions.h"
#include "debug.h"
#include "matpack_view.h"

namespace fwd::profile {
irradiance::irradiance(const Vector& z,
                       const Vector& p,
                       const Vector& t,
                       const std::vector<Vector>& allvmrs,
                       const ArrayOfArrayOfSpeciesTag& allspecs,
                       const PredefinedModelData& predef_data,
                       const ArrayOfCIARecord& cia_data,
                       const ArrayOfXsecRecord& xsec_data,
                       const SpeciesIsotopologueRatios& isotopologue_ratios,
                       const ArrayOfArrayOfAbsorptionLines& lbl_data,
                       Numeric cia_extrap,
                       Index cia_robust)
    : rad(z,
          p,
          t,
          allvmrs,
          allspecs,
          predef_data,
          cia_data,
          xsec_data,
          isotopologue_ratios,
          lbl_data,
          cia_extrap,
          cia_robust) {}

void irradiance::planar(ExhaustiveVectorView irr,
                        Numeric f,
                        const Index streams) const {
  const std::size_t n = rad.altitude.size();
  ARTS_USER_ERROR_IF(streams < 1 or streams % 2 not_eq 0,
                     "streams must be a positive even number, got: ",
                     streams)
  ARTS_USER_ERROR_IF(irr.size() not_eq static_cast<Index>(n),
                     "Size mismatch: irradiance and altitude (",
                     irr.size(),
                     " vs ",
                     n,
                     ')')

  const auto cos_za = [streams](Index i) {
    return (i == streams - 1) ? -1.0
                              : 1.0 - 2.0 * static_cast<Numeric>(i) /
                                          (static_cast<Numeric>(streams) - 1.0);
  };

  const Numeric dza_daa = Constant::pi * (1.0 - cos_za(1));

  Vector rads(n);
  for (Index i = 0; i < streams; i++) {
    const Numeric za = Conversion::acosd(cos_za(i));
    rad.planar(rads, f, za);

    if (i == 0) {
      for (std::size_t j = 0; j < n; j++) {
        irr[j] = 0.5 * rads[j] * dza_daa;
      }
    } else if (i == streams - 1) {
      for (std::size_t j = 0; j < n; j++) {
        irr[j] += 0.5 * rads[j] * dza_daa;
      }
    } else {
      for (std::size_t j = 0; j < n; j++) {
        irr[j] += rads[j] * dza_daa;
      }
    }
  }
}

Vector irradiance::planar(Numeric f, const Index streams) const {
  Vector irr(rad.altitude.size());
  planar(irr, f, streams);
  return irr;
}

Matrix irradiance::planar_par(const Vector& f, const Index streams) const {
  const Index n = f.size();
  Matrix irr(n, rad.altitude.size());
#pragma omp parallel for
  for (Index i = 0; i < n; i++) planar(irr[i], f[i], streams);
  return irr;
}

std::ostream& operator<<(std::ostream& os, const irradiance&) {
  return os << "No useful output operator yet";
}
}  // namespace fwd::profile
