#include <arts_omp.h>
#include <disort.h>
#include <lagrange_interp.h>
#include <vdisort_settings.h>
#include <workspace.h>

/* The VDISORT counterparts of the DISORT solver methods (m_disort.cc), one for one, with Stokes vectors. */

void vdisort_spectral_rad_fieldCalc(VDisortRadiance&       vdisort_spectral_rad_field,
                                    ZenGriddedField1&      quadrature,
                                    const VDisortSettings& vdisort_settings,
                                    const AziGrid&         phis) try {
  ARTS_TIME_REPORT

  const Index nv    = vdisort_settings.frequency_count();
  const Index np    = vdisort_settings.layer_count();
  const Index nquad = vdisort_settings.quadrature_dimension;
  const Index nphi  = phis.size();

  vdisort::main_data dis = vdisort_settings.init();

  //! Supplementary outputs
  quadrature = vdisort_quadrature_weights(dis);

  //! Main output
  vdisort_spectral_rad_field.resize(
      vdisort_settings.freq_grid, vdisort_settings.alt_grid, phis, ZenGrid{quadrature.grid<0>()});

  // The propagation azimuths as VDISORT's
  Vector phi(nphi);
  for (Index p = 0; p < nphi; p++) phi[p] = vdisort_azimuth(phis[p]);

  Tensor4 u(np, nphi, nquad, 4);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis, u)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      vdisort_settings.set(dis, iv);
      dis.gridded_u(u, phi);
      for (Index l = 0; l < np; l++)
        for (Index p = 0; p < nphi; p++)
          for (Index i = 0; i < nquad; i++)
            vdisort_spectral_rad_field.data[iv, l, p, i] = Stokvec{u[l, p, i, 0], u[l, p, i, 1], u[l, p, i, 2], u[l, p, i, 3]};
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  vdisort_spectral_rad_field.sort(dis.mu());

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in vdisort-spectral:\n{}", error);
}
ARTS_METHOD_ERROR_CATCH

void vdisort_spectral_flux_fieldCalc(SpectralFluxField&     spectral_flux_field,
                                     const VDisortSettings& vdisort_settings) try {
  ARTS_TIME_REPORT

  const Index nv = vdisort_settings.frequency_count();

  spectral_flux_field.resize(vdisort_settings.freq_grid, vdisort_settings.alt_grid);

  vdisort::main_data dis = vdisort_settings.init();

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      vdisort_settings.set(dis, iv);

      dis.gridded_flux(spectral_flux_field.up[iv],
                       spectral_flux_field.down_diffuse[iv],
                       spectral_flux_field.down_direct[iv],
                       spectral_flux_field.dfdt[iv]);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in vdisort:\n{}", error);
}
ARTS_METHOD_ERROR_CATCH

void spectral_radFromVDisort(StokvecVector&              spectral_rad,
                             const VDisortRadiance&      vdisort_spectral_rad_field,
                             const PropagationPathPoint& ray_point) {
  using namespace lagrange_interp;
  using aa_cyc_t = cycler<0.0, 360.0>;

  ARTS_TIME_REPORT

  const auto& f_grid   = vdisort_spectral_rad_field.freq_grid;
  const auto& alt_grid = vdisort_spectral_rad_field.alt_grid;
  const auto& azi_grid = vdisort_spectral_rad_field.azi_grid;
  const auto& zen_grid = vdisort_spectral_rad_field.zen_grid;
  const auto& data     = vdisort_spectral_rad_field.data;

  const Size nf = f_grid.size();
  spectral_rad.resize(nf);

  if (nf == 0) return;

  ARTS_USER_ERROR_IF(alt_grid.size() < 2, "VDISORT altitude grid must have at least two points")

  const Numeric z     = ray_point.altitude();
  const bool    above = alt_grid.size() < 2 or z >= alt_grid[1];
  const bool    below = z <= alt_grid.back();

  // One Stokes component at a time, as spectral_radFromDisort interpolates I
  Tensor4 component(data.shape());
  for (Index s = 0; s < 4; s++) {
    for (Size i = 0; i < data.size(); i++) component.data_handle()[i] = data.data_handle()[i][s];

    if (above or below) {
      std::visit(
          [&, idx = above ? 0 : alt_grid.size() - 2](auto&& aa, auto&& za) {
            for (Size v = 0; v < nf; v++) spectral_rad[v][s] = interp(component[v, idx], aa, za);
          },
          variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]),
          variant_lag(zen_grid, ray_point.los[0]));
    } else {
      std::visit(
          [&](auto&& alt, auto&& aa, auto&& za) {
            for (Size v = 0; v < nf; v++) spectral_rad[v][s] = interp(component[v], alt, aa, za);
          },
          variant_lag(std::span{alt_grid}.subspan(1), z),
          variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]),
          variant_lag(zen_grid, ray_point.los[0]));
    }
  }
}

void vdisort_spectral_rad_fieldApplyUnit(VDisortRadiance&                         vdisort_spectral_rad_field,
                                         const PropagationPathPoint&              ray_point,
                                         const SpectralRadianceTransformOperator& spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecVector spectral_rad(vdisort_spectral_rad_field.freq_grid.size());
  StokvecMatrix spectral_rad_jac(0, vdisort_spectral_rad_field.freq_grid.size());

  const auto [nv, np, nz, nq] = vdisort_spectral_rad_field.data.shape();

  std::string error{};

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3) firstprivate(spectral_rad, spectral_rad_jac)
  for (Index i = 0; i < np; i++) {
    for (Index j = 0; j < nz; j++) {
      for (Index k = 0; k < nq; k++) {
        try {
          for (Index v = 0; v < nv; v++) spectral_rad[v] = vdisort_spectral_rad_field.data[v, i, j, k];

          spectral_rad_transform_operator(
              spectral_rad, spectral_rad_jac, vdisort_spectral_rad_field.freq_grid, ray_point);

          for (Index v = 0; v < nv; v++) vdisort_spectral_rad_field.data[v, i, j, k] = spectral_rad[v];
        } catch (std::exception& e) {
#pragma omp critical
          if (error.empty()) error = e.what();
        }
      }
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), error)
}
ARTS_METHOD_ERROR_CATCH
