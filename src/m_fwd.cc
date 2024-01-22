#include <fwd.h>
#include <workspace.h>

#include <algorithm>

#include "atm.h"
#include "atm_path.h"
#include "debug.h"
#include "matpack_data.h"

void spectral_radiance_profile_operatorPlaneParallel(
    SpectralRadianceProfileOperator& spectral_radiance_profile_operator,
    const AtmField& atm_field,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const PredefinedModelData& predefined_model_data,
    const ArrayOfCIARecord& abs_cia_data,
    const ArrayOfXsecRecord& xsec_fit_data,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Numeric& cia_extrap,
    const Index& cia_robust,
    const Vector& z_grid,
    const Vector& lat_grid,
    const Vector& lon_grid) {
  const Index n = z_grid.size();

  ARTS_USER_ERROR_IF(n == 0, "Must have z_grid.size() > 0")
  ARTS_USER_ERROR_IF(not std::ranges::is_sorted(z_grid),
                     "z_grid must be sorted in ascending order\nz_grid: ",
                     z_grid)

  const ArrayOfAtmPoint ppvar_atm =
      extract1D(atm_field, z_grid, lat_grid, lon_grid);

  spectral_radiance_profile_operator =
      SpectralRadianceProfileOperator(z_grid,
                                      ppvar_atm,
                                      abs_species,
                                      predefined_model_data,
                                      abs_cia_data,
                                      xsec_fit_data,
                                      isotopologue_ratios,
                                      abs_lines_per_species,
                                      cia_extrap,
                                      cia_robust);
}

void spectral_radiance_fieldPlaneParallelSpectralRadianceOperator(
    Tensor7& spectral_radiance_field,
    const Vector& f_grid,
    const Vector& za_grid,
    const SpectralRadianceProfileOperator& spectral_radiance_profile_operator) {
  const Index n = f_grid.size();
  const Index m = za_grid.size();

  spectral_radiance_field.resize(
      n, spectral_radiance_profile_operator.altitude.size(), 1, 1, m, 1, 1);

  String error_msg;
#pragma omp parallel for collapse(2)
  for (Index i = 0; i < n; ++i) {
    for (Index j = 0; j < m; ++j) {
      if (error_msg.size()) continue;
      try {
        spectral_radiance_field(i, joker, 0, 0, j, 0, 0) =
            spectral_radiance_profile_operator.planar(f_grid(i), za_grid(j));
      } catch (std::exception& e) {
#pragma omp critical
        if (error_msg.size() == 0) error_msg = e.what();
      }
    }
  }

  ARTS_USER_ERROR_IF(error_msg.size(), error_msg, '\n')
}
