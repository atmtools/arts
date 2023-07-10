#include "abs.h"

namespace profile {
struct plane_parallel_planck_surface_no_polarization {
  std::vector<Numeric> altitude;
  std::vector<Numeric> temperature;
  std::vector<fwd::full> models;

  plane_parallel_planck_surface_no_polarization() = default;

  plane_parallel_planck_surface_no_polarization(
      const Vector& z,
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
  
  [[nodiscard]] Vector profile_at(Numeric f, Numeric za) const;
};
}  // namespace profile
