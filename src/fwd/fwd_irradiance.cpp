#include "fwd_irradiance.h"
#include <algorithm>
#include "arts_constants.h"
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
                       Index cia_robust,
                       Verbosity verb)
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
          cia_robust,
          verb) {}

void irradiance::planar(ExhaustiveVectorView irr,
                        Numeric f,
                        const Vector& za) const {
  const std::size_t n = rad.altitude.size();
  const Index m = za.size();
  ARTS_USER_ERROR_IF(
      not static_cast<bool>(m) or za[0] not_eq 0.0 or za[m - 1] not_eq 180.0,
      "We need zenith angles, [0, ..., 180], got: [",
      za,
      "]\n")
  ARTS_USER_ERROR_IF(irr.size() not_eq static_cast<Index>(n),
                     "Size mismatch: irradiance and altitude (",
                     irr.size(),
                     " vs ",
                     n,
                     ')')
  ARTS_USER_ERROR_IF(not std::is_sorted(za.begin(), za.end()),
                     "Zenith angles must be sorted in ascending order")

  // Do first step manually to reset input irradiance
  Numeric h_past = 1.0, h_this = Conversion::cosd(za[1]);
  Vector rads_past=rad.planar(f, za[0]), rads_this=rad.planar(f, za[1]);
  Numeric scl = 0.5 * (h_past - h_this);
  for (std::size_t j = 0; j < n; j++) {
    irr[j] = scl * std::midpoint(rads_this[j], rads_past[j]);
  }
  h_past = h_this;
  rads_past.swap(rads_this);

  for (Index i = 2; i < m; i++) {
    h_this = Conversion::cosd(za[i]);
    rad.planar(rads_this, f, za[i]);

    scl = 0.5 * (h_past - h_this);
    for (std::size_t j = 0; j < n; j++) {
      irr[j] += scl * std::midpoint(rads_this[j], rads_past[j]);
    }

    h_past = h_this;
    rads_past.swap(rads_this);
  }
}

Vector irradiance::planar(Numeric f, const Vector& za) const {
  Vector irr(rad.altitude.size());
  planar(irr, f, za);
  return irr;
}

std::ostream& operator<<(std::ostream& os, const irradiance&) {
  return os << "No useful output operator yet";
}
}  // namespace fwd::profile
