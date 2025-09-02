#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <format>

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
    Tensor4& disort_spectral_radiance_field,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point,
    const SpectralRadianceTransformOperator&
        spectral_radiance_transform_operator) try {
  ARTS_TIME_REPORT

  StokvecVector spectral_radiance(frequency_grid.size());
  StokvecMatrix spectral_radiance_jacobian(0, frequency_grid.size());

  const Index nv = disort_spectral_radiance_field.nbooks();
  const Index np = disort_spectral_radiance_field.npages();
  const Index nz = disort_spectral_radiance_field.nrows();
  const Index nq = disort_spectral_radiance_field.ncols();

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3) \
    firstprivate(spectral_radiance, spectral_radiance_jacobian)
  for (Index i = 0; i < np; i++) {
    for (Index j = 0; j < nz; j++) {
      for (Index k = 0; k < nq; k++) {
        for (Index v = 0; v < nv; v++) {
          spectral_radiance[v][0] = disort_spectral_radiance_field[v, i, j, k];
        }
        spectral_radiance_transform_operator(spectral_radiance,
                                             spectral_radiance_jacobian,
                                             frequency_grid,
                                             ray_path_point);
        for (Index v = 0; v < nv; v++) {
          disort_spectral_radiance_field[v, i, j, k] = spectral_radiance[v][0];
        }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH
