#include <fwd.h>

#include <algorithm>

#include "agenda_set.h"
#include "auto_md.h"
#include "messages.h"
#include "rte.h"

void spectral_radiance_profile_operatorPlaneParallel(
    Workspace& ws,
    SpectralRadianceProfileOperator& spectral_radiance_profile_operator,
    const Tensor3& z_field,
    const Numeric& ppath_lmax,
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Tensor3& t_field,
    const EnergyLevelMap& nlte_field,
    const Tensor4& vmr_field,
    const Tensor3& wind_u_field,
    const Tensor3& wind_v_field,
    const Tensor3& wind_w_field,
    const Tensor3& mag_u_field,
    const Tensor3& mag_v_field,
    const Tensor3& mag_w_field,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const PredefinedModelData& predefined_model_data,
    const ArrayOfCIARecord& abs_cia_data,
    const ArrayOfXsecRecord& xsec_fit_data,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const Numeric& cia_extrap,
    const Index& cia_robust,
    const Verbosity& verbosity) {
  const Agenda ppath_agenda =
      AgendaManip::get_ppath_agenda(ws, "PlaneParallel");

  Vector rte_pos{max(z_field) + 10};
  Vector rte_los{180};

  Ppath ppath;
  ppath_agendaExecute(ws,
                      ppath,
                      ppath_lmax,
                      ppath_lmax,
                      rte_pos,
                      rte_los,
                      {},
                      0,
                      0,
                      {},
                      ppath_agenda);

  Vector ppvar_p, ppvar_t;
  EnergyLevelMap ppvar_nlte;
  Matrix ppvar_vmr, ppvar_wind, ppvar_mag;
  get_ppath_atmvars(ppvar_p,
                    ppvar_t,
                    ppvar_nlte,
                    ppvar_vmr,
                    ppvar_wind,
                    ppvar_mag,
                    ppath,
                    atmosphere_dim,
                    p_grid,
                    t_field,
                    nlte_field,
                    vmr_field,
                    wind_u_field,
                    wind_v_field,
                    wind_w_field,
                    mag_u_field,
                    mag_v_field,
                    mag_w_field);

  const Vector z = reverse(Vector{ppath.pos(joker, 0)});
  const Vector p = reverse(ppvar_p);
  const Vector t = reverse(ppvar_t);

  ppvar_vmr = transpose(ppvar_vmr);
  std::vector<Vector> allvmrs{ppvar_vmr.begin(), ppvar_vmr.end()};
  std::reverse(allvmrs.begin(), allvmrs.end());

  spectral_radiance_profile_operator =
      SpectralRadianceProfileOperator(z,
                                      p,
                                      t,
                                      allvmrs,
                                      abs_species,
                                      predefined_model_data,
                                      abs_cia_data,
                                      xsec_fit_data,
                                      isotopologue_ratios,
                                      abs_lines_per_species,
                                      cia_extrap,
                                      cia_robust,
                                      verbosity);
}

void spectral_radiance_fieldPlaneParallelSpectralRadianceOperator(
    Tensor7& spectral_radiance_field,
    const SpectralRadianceProfileOperator& spectral_radiance_profile_operator,
    const Vector& f_grid,
    const Vector& za_grid,
    const Verbosity&) {
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
}
