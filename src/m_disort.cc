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
                                   const AziGrid& phis) try {
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
ARTS_METHOD_ERROR_CATCH

void disort_spectral_flux_fieldCalc(DisortFlux& disort_spectral_flux_field,
                                    const DisortSettings& disort_settings) try {
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
ARTS_METHOD_ERROR_CATCH

////////////////////////////////////////////////////////////////////////
// Coupled Disort
////////////////////////////////////////////////////////////////////////

void disort_spectral_rad_fieldCoupledCalc(
    DisortRadiance& disort_spectral_rad_field,
    ZenGriddedField1& disort_quadrature,
    const DisortSettings& atm_disort_settings,
    const DisortSettings& subsurf_disort_settings,
    const AziGrid& phis,
    const Numeric& tolerance,
    const Index& max_iterations,
    const Numeric& relaxation) try {
  ARTS_TIME_REPORT

  const Index nv    = atm_disort_settings.frequency_count();
  const Index np    = atm_disort_settings.layer_count();
  const Index nquad = atm_disort_settings.quadrature_dimension;

  disort::main_data dis_atm     = atm_disort_settings.init();
  disort::main_data dis_subsurf = subsurf_disort_settings.init();
  Tensor3 ims(np, phis.size(), nquad / 2);
  Tensor3 tms(np, phis.size(), nquad);

  //! Supplementary outputs
  disort_quadrature = dis_atm.gridded_weights();

  //! Main outputs
  DisortRadiance disort_atm_spectral_rad_field,
      disort_subsurf_spectral_rad_field;
  disort_atm_spectral_rad_field.resize(atm_disort_settings.freq_grid,
                                       atm_disort_settings.alt_grid,
                                       phis,
                                       disort_quadrature.grid<0>());
  disort_subsurf_spectral_rad_field.resize(subsurf_disort_settings.freq_grid,
                                           subsurf_disort_settings.alt_grid,
                                           phis,
                                           disort_quadrature.grid<0>());

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) \
    firstprivate(dis_atm, dis_subsurf, tms, ims)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      atm_disort_settings.set(dis_atm, iv);
      subsurf_disort_settings.set(dis_subsurf, iv);

      const auto _ = disort::couple(
          dis_atm, dis_subsurf, tolerance, max_iterations, relaxation);

      dis_atm.gridded_u_corr(
          disort_atm_spectral_rad_field.data[iv], tms, ims, phis);
      dis_subsurf.gridded_u_corr(
          disort_subsurf_spectral_rad_field.data[iv], tms, ims, phis);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(
      error.size(), "Error occurred in disort-spectral:\n{}", error);

  //! FIXME: It would be nice to remove this if the internal angles can be solved
  disort_atm_spectral_rad_field.sort(dis_atm.mu());
  disort_subsurf_spectral_rad_field.sort(dis_subsurf.mu());

  disort_spectral_rad_field =
      disort_atm_spectral_rad_field.combine(disort_subsurf_spectral_rad_field);
}
ARTS_METHOD_ERROR_CATCH

void disort_spectral_flux_fieldCoupledCalc(
    DisortFlux& disort_spectral_flux_field,
    const DisortSettings& atm_disort_settings,
    const DisortSettings& subsurf_disort_settings,
    const Numeric& tolerance,
    const Index& max_iterations,
    const Numeric& relaxation) try {
  ARTS_TIME_REPORT

  const Index nv = atm_disort_settings.frequency_count();

  DisortFlux disort_atm_spectral_flux_field, disort_subsurf_spectral_flux_field;
  disort_atm_spectral_flux_field.resize(atm_disort_settings.freq_grid,
                                        atm_disort_settings.alt_grid);

  disort_subsurf_spectral_flux_field.resize(subsurf_disort_settings.freq_grid,
                                            subsurf_disort_settings.alt_grid);

  disort::main_data dis_atm, dis_subsurf;
  try {
    dis_atm = atm_disort_settings.init();
  } catch (std::exception& e) {
    throw std::runtime_error("Error in atmosphere disort settings:\n" +
                             std::string(e.what()));
  }
  try {
    dis_subsurf = subsurf_disort_settings.init();
  } catch (std::exception& e) {
    throw std::runtime_error("Error in subsurface disort settings:\n" +
                             std::string(e.what()));
  }

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) \
    firstprivate(dis_atm, dis_subsurf)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      atm_disort_settings.set(dis_atm, iv);
      subsurf_disort_settings.set(dis_subsurf, iv);

      const auto _ = disort::couple(
          dis_atm, dis_subsurf, tolerance, max_iterations, relaxation);

      dis_atm.gridded_flux(disort_atm_spectral_flux_field.up[iv],
                           disort_atm_spectral_flux_field.down_diffuse[iv],
                           disort_atm_spectral_flux_field.down_direct[iv]);
      dis_subsurf.gridded_flux(
          disort_subsurf_spectral_flux_field.up[iv],
          disort_subsurf_spectral_flux_field.down_diffuse[iv],
          disort_subsurf_spectral_flux_field.down_direct[iv]);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in disort:\n{}", error);

  disort_spectral_flux_field = disort_atm_spectral_flux_field.combine(
      disort_subsurf_spectral_flux_field);
}
ARTS_METHOD_ERROR_CATCH

////////////////////////////////////////////////////////////////////////
// Integrate functions
////////////////////////////////////////////////////////////////////////

void spectral_radFromDisort(StokvecVector& spectral_rad,
                            const DisortRadiance& disort_spectral_rad_field,
                            const PropagationPathPoint& ray_point) {
  using namespace lagrange_interp;
  using aa_cyc_t = cycler<0.0, 360.0>;

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

  const Numeric z  = ray_point.altitude();
  const bool above = alt_grid.size() < 2 or z >= alt_grid[1];
  const bool below = z <= alt_grid.back();

  if (above or below) {
    std::visit(
        [&, idx = above ? 0 : alt_grid.size() - 2](auto&& aa, auto&& za) {
          std::transform(
              data.begin(), data.end(), spectral_rad.begin(), [&](auto&& d) {
                return interp(d[idx], aa, za);
              });
        },
        variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]),
        variant_lag(zen_grid, ray_point.los[0]));
  } else {
    std::visit(
        [&](auto&& alt, auto&& aa, auto&& za) {
          std::transform(
              data.begin(), data.end(), spectral_rad.begin(), [&](auto&& d) {
                return interp(d, alt, aa, za);
              });
        },
        variant_lag(std::span{alt_grid}.subspan(1), z),
        variant_lag<aa_cyc_t>(azi_grid, ray_point.los[1]),
        variant_lag(zen_grid, ray_point.los[0]));
  }
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
