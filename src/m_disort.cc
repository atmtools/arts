#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <format>

#include "lagrange_interp.h"
#include "path_point.h"

////////////////////////////////////////////////////////////////////////
// Core Disort
////////////////////////////////////////////////////////////////////////

void disort_spectral_radiance_fieldCalc(
    DisortRadiance& disort_spectral_radiance_field,
    ZenithGriddedField1& disort_quadrature,
    const DisortSettings& disort_settings,
    const AzimuthGrid& phis) {
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
  disort_spectral_radiance_field.resize(
      disort_settings.frequency_grid,
      disort_settings.altitude_grid,
      phis,
      ZenithGrid{disort_quadrature.grid<0>()});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) \
    firstprivate(dis, tms, ims)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      disort_settings.set(dis, iv);
      dis.gridded_u_corr(
          disort_spectral_radiance_field.data[iv], tms, ims, phis);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  //! FIXME: It would be nice to remove this if the internal angles can be solved
  disort_spectral_radiance_field.sort(dis.mu());

  ARTS_USER_ERROR_IF(
      error.size(), "Error occurred in disort-spectral:\n{}", error);
}

void disort_spectral_flux_fieldCalc(DisortFlux& disort_spectral_flux_field,
                                    const DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  const Index nv = disort_settings.frequency_count();

  disort_spectral_flux_field.resize(disort_settings.frequency_grid,
                                    disort_settings.altitude_grid);

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

void spectral_radianceFromReverseDisort(
    StokvecVector& spectral_radiance,
    const DisortRadiance& disort_spectral_radiance_field,
    const PropagationPathPoint& ray_path_point) {
  ARTS_TIME_REPORT

  const auto& f_grid   = disort_spectral_radiance_field.frequency_grid;
  const auto& alt_grid = disort_spectral_radiance_field.altitude_grid;
  const auto& aa_grid  = disort_spectral_radiance_field.azimuth_grid;
  const auto& za_grid  = disort_spectral_radiance_field.zenith_grid;
  const auto& data     = disort_spectral_radiance_field.data;

  const Size nf = f_grid.size();
  spectral_radiance.resize(nf);

  if (nf == 0) return;

  ARTS_USER_ERROR_IF(alt_grid.size() < 2,
                     "DISORT altitude grid must have at least two points")

  //! FIXME: Is the azimuth angle supposed to be mirrored or not?
  const Vector2 mirror_los = path::mirror(ray_path_point.los);

  using aa_cyc_t = lagrange_interp::cycler<0.0, 360.0>;
  const auto aa_lag =
      lagrange_interp::variant_lag<aa_cyc_t>(aa_grid, mirror_los[1]);
  const auto za_lag = lagrange_interp::variant_lag<lagrange_interp::identity>(
      za_grid, mirror_los[0]);

  const Numeric z  = ray_path_point.altitude();
  const bool above = alt_grid.size() < 2 or z >= alt_grid[1];
  const bool below = z <= alt_grid.back();

  if (above or below) {
    std::visit(
        [&, idx = above ? 0 : alt_grid.size() - 2](auto& aa, auto& za) {
          for (Size iv = 0; iv < nf; iv++) {
            spectral_radiance[iv] = interp(data[iv, idx], aa, za);
          }
        },
        aa_lag,
        za_lag);
  } else {
    const auto alt_lag =
        lagrange_interp::variant_lag<>(std::span{alt_grid}.subspan(1), z);

    std::visit(
        [&](auto& alt, auto& aa, auto& za) {
          for (Size iv = 0; iv < nf; iv++) {
            spectral_radiance[iv] = interp(data[iv], alt, aa, za);
          }
        },
        alt_lag,
        aa_lag,
        za_lag);
  }
}

void spectral_radianceIntegrateDisort(
    StokvecVector& /*spectral_radiance*/,
    const DisortRadiance& /*disort_spectral_radiance_field*/,
    const ZenithGriddedField1& /*disort_quadrature*/) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR("Not implemented")
}

void SpectralFluxDisort(Matrix& /*spectral_flux_field_up*/,
                        Matrix& /*spectral_flux_field_down*/,
                        const DisortFlux& /*disort_spectral_flux_field*/) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR("Not implemented")
}

void disort_spectral_radiance_fieldApplyUnit(
    DisortRadiance& disort_spectral_radiance_field,
    const PropagationPathPoint& ray_path_point,
    const SpectralRadianceTransformOperator&
        spectral_radiance_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecVector spectral_radiance(
      disort_spectral_radiance_field.frequency_grid.size());
  StokvecMatrix spectral_radiance_jacobian(
      0, disort_spectral_radiance_field.frequency_grid.size());

  const Index nv = disort_spectral_radiance_field.data.nbooks();
  const Index np = disort_spectral_radiance_field.data.npages();
  const Index nz = disort_spectral_radiance_field.data.nrows();
  const Index nq = disort_spectral_radiance_field.data.ncols();

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3) \
    firstprivate(spectral_radiance, spectral_radiance_jacobian)
  for (Index i = 0; i < np; i++) {
    for (Index j = 0; j < nz; j++) {
      for (Index k = 0; k < nq; k++) {
        for (Index v = 0; v < nv; v++) {
          spectral_radiance[v][0] =
              disort_spectral_radiance_field.data[v, i, j, k];
        }
        spectral_radiance_transform_operator(
            spectral_radiance,
            spectral_radiance_jacobian,
            disort_spectral_radiance_field.frequency_grid,
            ray_path_point);
        for (Index v = 0; v < nv; v++) {
          disort_spectral_radiance_field.data[v, i, j, k] =
              spectral_radiance[v][0];
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
