#pragma once

#include <utility>

#include "fwd_radiance.h"

namespace fwd::profile {
struct irradiance {
  spectral_radiance rad;

  irradiance() = default;

  irradiance(const Vector& z,
             const ArrayOfAtmPoint& atm_points,
             const ArrayOfArrayOfSpeciesTag& allspecs,
             const PredefinedModelData& predef_data,
             const ArrayOfCIARecord& cia_data,
             const ArrayOfXsecRecord& xsec_data,
             const SpeciesIsotopologueRatios& isotopologue_ratios,
             const ArrayOfArrayOfAbsorptionLines& lbl_data,
             Numeric cia_extrap = {},
             Index cia_robust = {});
  
  irradiance(spectral_radiance  fwd_rad) : rad(std::move(fwd_rad)) {}

  //! FIXME: Will change in arts-3
  [[nodiscard]] Vector planar(Numeric f, const Index streams) const;

  //! FIXME: Will change in arts-3
  void planar(ExhaustiveVectorView irr, Numeric f, const Index streams) const;
  [[nodiscard]] Matrix planar_par(const Vector& f, const Index streams) const;

  friend std::ostream& operator<<(std::ostream&, const irradiance&);
};
}  // namespace fwd::profile
