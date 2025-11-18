#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <format>

////////////////////////////////////////////////////////////////////////
// Core Disort
////////////////////////////////////////////////////////////////////////

void disort_spectral_rad_fieldCalc(DisortRadiance& disort_spectral_rad_field,
                                   ZenGriddedField1& disort_quadrature,
                                   const DisortSettings& disort_settings,
                                   const AziGrid& phis) {
  ARTS_TIME_REPORT

  const Index nv    = disort_settings.frequency_count();
  const Index np    = disort_settings.layer_count();
  const Index nquad = disort_settings.quadrature_dimension;

  disort::main_data dis = disort_settings.init();
  Tensor3 ims(np, phis.size(), nquad / 2);
  Tensor3 tms(np, phis.size(), nquad);

  //! Supplementary outputs
  disort_quadrature = dis.gridded_weights();

  //! Main output
  disort_spectral_rad_field.resize(disort_settings.freq_grid,
                                   disort_settings.alt_grid,
                                   phis,
                                   ZenGrid{disort_quadrature.grid<0>()});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) \
    firstprivate(dis, tms, ims)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      disort_settings.set(dis, iv);
      dis.gridded_u_corr(disort_spectral_rad_field.data[iv], tms, ims, phis);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  //! FIXME: It would be nice to remove this if the internal angles can be solved
  disort_spectral_rad_field.sort(dis.mu());

  ARTS_USER_ERROR_IF(
      error.size(), "Error occurred in disort-spectral:\n{}", error);
}

void disort_spectral_flux_fieldCalc(DisortFlux& disort_spectral_flux_field,
                                    const DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  const Index nv = disort_settings.frequency_count();

  disort_spectral_flux_field.resize(disort_settings.freq_grid,
                                    disort_settings.alt_grid);

  disort::main_data dis = disort_settings.init();

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      disort_settings.set(dis, iv);

      dis.gridded_flux(disort_spectral_flux_field.up[iv],
                       disort_spectral_flux_field.down_diffuse[iv],
                       disort_spectral_flux_field.down_direct[iv]);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in disort:\n{}", error);
}

////////////////////////////////////////////////////////////////////////
// Integrate functions
////////////////////////////////////////////////////////////////////////

void spectral_radFromDisort(StokvecVector& spectral_rad,
                            const DisortRadiance& disort_spectral_rad_field,
                            const PropagationPathPoint& ray_point) {
  using namespace lagrange_interp;

  ARTS_TIME_REPORT

  const auto& f_grid   = disort_spectral_rad_field.freq_grid;
  const auto& alt_grid = disort_spectral_rad_field.alt_grid;
  const auto& azi_grid = disort_spectral_rad_field.azi_grid;
  const auto& zen_grid = disort_spectral_rad_field.zen_grid;
  const auto& data     = disort_spectral_rad_field.data;

  const Size nf = f_grid.size();
  spectral_rad.resize(nf);

  if (nf == 0) return;

  ARTS_USER_ERROR_IF(alt_grid.size() < 2,
                     "DISORT altitude grid must have at least two points")

  using aa_cyc_t    = cycler<0.0, 360.0>;
  const auto aa_lag = variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]);
  const auto za_lag = variant_lag<identity>(zen_grid, ray_point.los[0]);

  const Numeric z  = ray_point.altitude();
  const bool above = alt_grid.size() < 2 or z >= alt_grid[1];
  const bool below = z <= alt_grid.back();

  if (above or below) {
    std::visit(
        [&, idx = above ? 0 : alt_grid.size() - 2](auto& aa, auto& za) {
          std::transform(
              data.begin(), data.end(), spectral_rad.begin(), [&](auto&& d) {
                return interp(d[idx], aa, za);
              });
        },
        aa_lag,
        za_lag);
  } else {
    const auto alt_lag = variant_lag<>(std::span{alt_grid}.subspan(1), z);

    std::visit(
        [&](auto& alt, auto& aa, auto& za) {
          std::transform(
              data.begin(), data.end(), spectral_rad.begin(), [&](auto&& d) {
                return interp(d, alt, aa, za);
              });
        },
        alt_lag,
        aa_lag,
        za_lag);
  }
}

void spectral_radIntegrateDisort(
    StokvecVector& /*spectral_rad*/,
    const DisortRadiance& /*disort_spectral_rad_field*/,
    const ZenGriddedField1& /*disort_quadrature*/) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR("Not implemented")
}

void SpectralFluxDisort(Matrix& /*spectral_flux_field_up*/,
                        Matrix& /*spectral_flux_field_down*/,
                        const DisortFlux& /*disort_spectral_flux_field*/) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR("Not implemented")
}

void disort_spectral_rad_fieldApplyUnit(
    DisortRadiance& disort_spectral_rad_field,
    const PropagationPathPoint& ray_point,
    const SpectralRadianceTransformOperator&
        spectral_rad_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecVector spectral_rad(disort_spectral_rad_field.freq_grid.size());
  StokvecMatrix spectral_rad_jac(0, disort_spectral_rad_field.freq_grid.size());

  const Index nv = disort_spectral_rad_field.data.nbooks();
  const Index np = disort_spectral_rad_field.data.npages();
  const Index nz = disort_spectral_rad_field.data.nrows();
  const Index nq = disort_spectral_rad_field.data.ncols();

  std::string error{};

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3) \
    firstprivate(spectral_rad, spectral_rad_jac)
  for (Index i = 0; i < np; i++) {
    for (Index j = 0; j < nz; j++) {
      for (Index k = 0; k < nq; k++) {
        try {
          for (Index v = 0; v < nv; v++) {
            spectral_rad[v][0] = disort_spectral_rad_field.data[v, i, j, k];
          }

          spectral_rad_transform_operator(spectral_rad,
                                          spectral_rad_jac,
                                          disort_spectral_rad_field.freq_grid,
                                          ray_point);

          for (Index v = 0; v < nv; v++) {
            disort_spectral_rad_field.data[v, i, j, k] = spectral_rad[v][0];
          }
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
