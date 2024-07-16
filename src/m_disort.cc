#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include "atm.h"
#include "enums.h"
#include "physics_funcs.h"
#include "rtepack.h"
#include "sorted_grid.h"

ArrayOfAscendingGrid ray_path_tau_arrFromPath(
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfPropmatVector& ray_path_propagation_matrix) {
  const Size N = ray_path.size();
  ARTS_ASSERT(N > 1)
  ARTS_ASSERT(N == ray_path_propagation_matrix.size());

  const Vector r = [n = N - 1, &ray_path]() {
    Vector out(n);
    for (Size i = 0; i < n; i++)
      out[i] = ray_path[i + 1].altitude() - ray_path[i].altitude();

    return out;
  }();

  const Index nv = ray_path_propagation_matrix.front().size();

  ArrayOfAscendingGrid ray_path_tau_arr(nv);

  for (Index iv = 0; iv < nv; iv++) {
    Vector tau(N - 1);
    for (Size i = 0; i < N - 1; i++) {
      tau[i] = r[i] * std::midpoint(ray_path_propagation_matrix[i + 1][iv].A(),
                                    ray_path_propagation_matrix[i + 0][iv].A());
    }

    for (Size i = 1; i < N - 1; i++) {
      tau[i] += tau[i - 1];
    }

    ray_path_tau_arr[iv] = std::move(tau);
  }

  return ray_path_tau_arr;
}

void ray_pathGeometricUplooking(ArrayOfPropagationPathPoint& ray_path,
                                const AtmField& atmospheric_field,
                                const SurfaceField& surface_field,
                                const Numeric& latitude,
                                const Numeric& longitude,
                                const Numeric& max_step) {
  ray_pathGeometric(
      ray_path,
      atmospheric_field,
      surface_field,
      {surface_field.single_value(SurfaceKey::h, latitude, longitude),
       latitude,
       longitude},
      {0, 0},
      max_step,
      1.0,
      true,
      false,
      true,
      true,
      false);
}

void disort_intensitiesClearskyDisort(
    Tensor3& disort_intensities,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfAscendingGrid& ray_path_frequency_grid,
    const Index& NQuad,
    const Index& NLeg_,
    const Index& NFourier_) {
  const Size N = ray_path.size();
  ARTS_USER_ERROR_IF(N < 2, "Cannot determine sizes from 0-, or 1-sized input")

  // Input has same array sizes
  ARTS_USER_ERROR_IF(N != ray_path_propagation_matrix.size(),
                     "Incorrect number of propagation matrices.");
  ARTS_USER_ERROR_IF(N != ray_path_frequency_grid.size(),
                     "Incorrect number of frequency grids.");

  const Index nv = ray_path_frequency_grid.front().size();

  // Input has same frequency sizes
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          ray_path_frequency_grid, Cmp::ne(nv), &AscendingGrid::size),
      "Incorrect size for propagation matrices and frequency grids.");
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          ray_path_propagation_matrix, Cmp::ne(nv), &PropmatVector::size),
      "Incorrect size for propagation matrices and frequency grids.");

  // No polarization allowed
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_propagation_matrix,
                          [](const PropmatVector& pms) {
                            return std::ranges::any_of(
                                pms, Cmp::eq(true), &Propmat::is_polarized);
                          }),
      "No implementation for polarized propagation matrices.");

  disort_intensities.resize(nv, N - 1, NQuad);
  disort_intensities = 0.0;

  // Altitude is increasing
  ARTS_USER_ERROR_IF(
      std::ranges::is_sorted(
          ray_path,
          [](const PropagationPathPoint& a, const PropagationPathPoint& b) {
            return a.altitude() > b.altitude();
          }),
      "Ray path points must be sorted by increasing altitude.");

  const Index NLeg     = NLeg_ < 0 ? NQuad : NLeg_;
  const Index NFourier = NFourier_ < 0 ? NLeg : NFourier_;

  ARTS_USER_ERROR_IF(NLeg < 1, "Must be at least 1")

  disort::main_data dis(N - 1, NQuad, NLeg, NFourier, 2, NLeg, 1);

  dis.solar_zenith() = 1.0;
  dis.beam_azimuth() = 0.0;
  dis.brdf_modes()[0] =
      disort::BDRF{[](ExhaustiveMatrixView x,
                      const ExhaustiveConstVectorView&,
                      const ExhaustiveConstVectorView&) { x = 0.0; }};
  dis.omega()                         = 0.0;
  dis.f()                             = 0.0;
  dis.all_legendre_coeffs()           = 0.0;
  dis.all_legendre_coeffs()(joker, 0) = 1.0;
  dis.positive_boundary()             = 0.0;
  dis.negative_boundary()             = 0.0;

  const ArrayOfAscendingGrid ray_path_tau_arr =
      ray_path_tau_arrFromPath(ray_path, ray_path_propagation_matrix);

#pragma omp parallel for firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    dis.tau() = ray_path_tau_arr[iv];
    for (Size i = 0; i < N - 1; i++) {
      const Numeric& f0 = ray_path_frequency_grid[i + 0][iv];
      const Numeric& f1 = ray_path_frequency_grid[i + 1][iv];

      const Numeric& t0 = ray_path_atmospheric_point[i + 0].temperature;
      const Numeric& t1 = ray_path_atmospheric_point[i + 1].temperature;

      const Numeric y0 = planck(f0, t0);
      const Numeric y1 = planck(f1, t1);

      dis.source_poly()(i, 0) = std::midpoint(y0, y1);
      dis.source_poly()(i, 1) = 0;
    }

    dis.update_all(0.0);

    dis.gridded_u(disort_intensities[iv].reshape_as(N - 1, 1, NQuad), {0.0});
  }
}
