#include "fwd_abs.h"
#include <ostream>

namespace fwd::profile {
struct full {
  Numeric refl;
  std::vector<Numeric> altitude;
  std::vector<Numeric> temperature;  // FIXME: Should be AtmPoint...
  std::vector<full_absorption> models;

  full() = default;

  full(const Vector& z,
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
  [[nodiscard]] Vector plane_par(Numeric f, Numeric za) const;

  //! FIXME: Will change in arts-3
  ExhaustiveVectorView plane_par(ExhaustiveVectorView, Numeric f, Numeric za) const;

  friend std::ostream& operator<<(std::ostream&, const full&);
};
}  // namespace fwd::profile
