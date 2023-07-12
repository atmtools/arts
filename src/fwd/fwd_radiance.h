#pragma once

#include <ostream>

#include "fwd_abs.h"

namespace fwd::profile {
struct radiance {
  Numeric refl;
  std::vector<Numeric> altitude;
  std::vector<Numeric> temperature;  // FIXME: Should be AtmPoint...
  std::vector<full_absorption> models;

  radiance() = default;

  radiance(const Vector& z,
           const Vector& p,
           const Vector& t,
           const std::vector<Vector>& allvmrs,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const PredefinedModelData& predef_data,
           const ArrayOfCIARecord& cia_data,
           const SpeciesIsotopologueRatios& isotopologue_ratios,
           const ArrayOfArrayOfAbsorptionLines& lbl_data,
           Numeric cia_extrap = {},
           Index cia_robust = {},
           Verbosity cia_verb = {});

  //! FIXME: Will change in arts-3
  [[nodiscard]] Vector planar(Numeric f, Numeric za) const;

  //! FIXME: Will change in arts-3
  ExhaustiveVectorView planar(ExhaustiveVectorView rad,
                              Numeric f,
                              Numeric za) const;

  //! FIXME: Will change in arts-3
  [[nodiscard]] Matrix planar_par(const Vector& f, Numeric za) const;
  void planar_par(ExhaustiveMatrixView rad, const Vector& f, Numeric za) const;

  friend std::ostream& operator<<(std::ostream&, const radiance&);
};
}  // namespace fwd::profile
