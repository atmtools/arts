#include <algorithm>
#include <chrono>
#include <numbers>
#include <ranges>

#include <iostream>

#include "arts_conversions.h"
#include "interpolation.h"
#include "matpack/matpack.h"
#include "scattering/mie.h"
#include "scattering/phase_matrix.h"
#include "test_utils.h"
#include "matpack/matpack_mdspan_helpers_eigen.h"

using namespace scattering;

using PhaseMatrixTROGridded =
    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded>;

using PhaseMatrixTROSpectral =
    PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded>;

using PhaseMatrixAROGridded =
    PhaseMatrixData<Numeric, Format::ARO, Representation::Gridded>;

using PhaseMatrixAROSpectral =
    PhaseMatrixData<Numeric, Format::ARO, Representation::Spectral>;

/** Create a TRO phase matrix for testing
 *
 * Creates a phase matrix with Legendre polynomials with the degree
 * similar to the frequency index.
 */
PhaseMatrixTROGridded make_phase_matrix(
    std::shared_ptr<const Vector> t_grid,
    std::shared_ptr<const Vector> f_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid) {
  PhaseMatrixTROGridded phase_matrix(t_grid, f_grid, za_scat_grid);
  Vector za_grid  = Vector{grid_vector(*za_scat_grid)};
  za_grid        *= Conversion::deg2rad(1.0);
  Vector aa_grid(1);
  aa_grid = 0.0;

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
        phase_matrix[i_t, Range(i_f, 1), joker, i_s] =
            evaluate_spherical_harmonic(i_f, 0, aa_grid, za_grid);
      }
    }
  }
  return phase_matrix;
}

/** Create a TRO phase matrix for a liquid sphere.
 *
 * Creates phase matrix data for a liquid sphere with a radius of 100um.
 */
PhaseMatrixTROGridded make_phase_matrix_liquid_sphere(
    std::shared_ptr<const Vector> t_grid,
    std::shared_ptr<const Vector> f_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid) {
  PhaseMatrixTROGridded phase_matrix(t_grid, f_grid, za_scat_grid);

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      auto scat_data = scattering::MieSphere<Numeric>::Liquid(
          (*f_grid)[i_f], (*t_grid)[i_t], 1e-3, grid_vector(*za_scat_grid));
      auto scat_matrix = scat_data.get_scattering_matrix_compact();
      for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
        phase_matrix[i_t, i_f, joker, i_s] = scat_matrix[joker, i_s];
      }
    }
  }
  return phase_matrix;
}

/** Create a ARO phase matrix for testing
 *
 * Creates a phase matrix with Legendre polynomials with the degree
 * similar to the frequency index and order similar to the temperature
 * index.
 */
PhaseMatrixAROGridded make_phase_matrix(
    std::shared_ptr<const Vector> t_grid,
    std::shared_ptr<const Vector> f_grid,
    std::shared_ptr<const Vector> za_inc_grid,
    std::shared_ptr<const Vector> delta_aa_grid,
    std::shared_ptr<const ZenithAngleGrid> za_scat_grid) {
  PhaseMatrixAROGridded phase_matrix(
      t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
  Vector aa_grid  = *delta_aa_grid;
  aa_grid        *= Conversion::deg2rad(1.0);
  Vector za_grid  = Vector{grid_vector(*za_scat_grid)};
  za_grid        *= Conversion::deg2rad(1.0);

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      for (Size i_za_inc = 0; i_za_inc < za_inc_grid->size(); ++i_za_inc) {
        for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
          Index l = i_t;
          Index m = std::min(i_t, i_f);
          phase_matrix[i_t, i_f, i_za_inc, joker, joker, i_s] =
              evaluate_spherical_harmonic(l, m, aa_grid, za_grid);
        }
      }
    }
  }
  return phase_matrix;
}

bool test_phase_matrix_tro() {
  auto sht          = sht::provider.get_instance(1, 32);
  auto t_grid       = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid       = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid = std::make_shared<ZenithAngleGrid>(*sht->get_za_grid_ptr());
  auto phase_matrix_gridded =
      make_phase_matrix(t_grid, f_grid, za_scat_grid);

  //
  // Test conversion between spectral and gridded format.
  //
  auto phase_matrix_spectral = phase_matrix_gridded.to_spectral(sht);
  ComplexVector coeffs_ref(phase_matrix_spectral.extent(2));
  for (Index i_t = 0; i_t < phase_matrix_spectral.extent(0); ++i_t) {
    for (Index i_f = 0; i_f < phase_matrix_spectral.extent(1); ++i_f) {
      coeffs_ref      = std::complex<Numeric>(0.0, 0.0);
      coeffs_ref[i_f] = std::complex<Numeric>(1.0, 0.0);
      for (Index i_s = 0; i_s < phase_matrix_spectral.extent(3); ++i_s) {
        Numeric err = max_error<ComplexVector>(
            coeffs_ref,
            static_cast<ComplexVector>(
                phase_matrix_spectral[i_t, i_f, joker, i_s]));
        if (err > 1e-6) return false;
      }
    }
  }
  auto phase_matrix_gridded_2 = phase_matrix_spectral.to_gridded();
  Numeric err = max_error(phase_matrix_gridded, phase_matrix_gridded_2);
  if (err > 1e-6) {
    return false;
  }

  //
  // Test conversion to lab frame.
  //

  std::shared_ptr<Vector> delta_aa_grid =
      std::make_shared<Vector>(Vector({0.0, 180}));
  std::shared_ptr<Vector> za_inc_grid =
      std::make_shared<Vector>(Vector({90.0}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_new =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({90.0})));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_liquid =
      std::make_shared<ZenithAngleGrid>(
          IrregularZenithAngleGrid(Vector({0.0, 10, 20, 160, 180.0})));

  auto phase_matrix_liquid =
      make_phase_matrix_liquid_sphere(t_grid, f_grid, za_scat_grid_liquid);
  auto phase_matrix_lab = phase_matrix_liquid.to_lab_frame(
      za_inc_grid, delta_aa_grid, za_scat_grid_new);

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      // Backward scattering direction. Only two independent elements.
      // Off-diagonal elements must be close to 0.
      auto pm_f           = phase_matrix_lab[i_t, i_f, 0, 0, 0, joker];
      Numeric coeff_max_f = pm_f[0];
      auto pm_b           = phase_matrix_lab[i_t, i_f, 0, 1, 0, joker];
      Numeric coeff_max_b = pm_b[0];

      // For the forward direction we have:
      // Z22 == Z33
      Numeric delta = std::abs(pm_f[1 * 4 + 1] - pm_f[2 * 4 + 2]) / coeff_max_f;
      if (delta > 1e-3) return false;

      // For the backward direction we have:
      // Z11 == -Z33
      delta = std::abs(pm_b[1 * 4 + 1] + pm_b[2 * 4 + 2]) / coeff_max_b;
      if (delta > 1e-3) return false;
      // Z44 == Z11 - 2 * Z22
      delta =
          std::abs(pm_b[0 * 4 + 0] - 2.0 * pm_b[1 * 4 + 1] - pm_b[3 * 4 + 3]) /
          coeff_max_b;
      if (delta > 1e-3) return false;

      // And all off-diagonal elements should be zero.
      for (Index i_s1 = 0; i_s1 < 4; ++i_s1) {
        for (Index i_s2 = 0; i_s2 < 4; ++i_s2) {
          if (i_s1 != i_s2) {
            Numeric c = pm_b[i_s1 * 4 + i_s2];
            if (std::abs(c) / coeff_max_b > 1e-3) {
              return false;

              c = pm_f[i_s1 * 4 + i_s2];
              if (std::abs(c) / coeff_max_f > 1e-3) {
                return false;
              }
            }
          }
        }
      }
    }
  }        


  // Test extraction of backscatter matrix.
  auto backscatter_matrix = phase_matrix_liquid.extract_backscatter_matrix();
  err                     = max_error<Tensor3>(
      backscatter_matrix,
      static_cast<Tensor3>(phase_matrix_liquid[joker, joker, 4, joker]));
  if (err > 1e-15) return false;

  auto forwardscatter_matrix =
      phase_matrix_liquid.extract_forwardscatter_matrix();
  err = max_error<Tensor3>(
      forwardscatter_matrix,
      static_cast<Tensor3>(phase_matrix_liquid[joker, joker, 0, joker]));
  if (err > 1e-15) return false;

  auto phase_matrix_liquid_sht =
      make_phase_matrix_liquid_sphere(t_grid, f_grid, za_scat_grid);
  auto phase_matrix_liquid_spectral = phase_matrix_liquid_sht.to_spectral(sht);
  auto backscatter_matrix_2 =
      phase_matrix_liquid_spectral.extract_backscatter_matrix();
  err            = max_error<Tensor3>(backscatter_matrix, backscatter_matrix_2);
  Tensor3 delta  = backscatter_matrix;
  delta         -= backscatter_matrix_2;
  if (err > 1e-15) return false;

  auto forwardscatter_matrix_2 =
      phase_matrix_liquid_spectral.extract_forwardscatter_matrix();
  err = max_error<Tensor3>(forwardscatter_matrix, forwardscatter_matrix_2);
  if (err > 1e-15) return false;

  // Test reduction of stokes elements.
  auto phase_matrix_liquid_1 = phase_matrix_liquid.extract_stokes_coeffs();
  err                        = max_error<ConstTensor4View>(
      phase_matrix_liquid_1,
      phase_matrix_liquid);
  // Extraction of stokes parameters should be exact.
  if (err > 0.0) return false;

  auto phase_matrix_liquid_spectral_1 =
      phase_matrix_liquid_spectral.extract_stokes_coeffs();
  err = max_error<matpack::strided_view_t<const std::complex<Numeric>, 4>>(
      phase_matrix_liquid_spectral_1,
      phase_matrix_liquid_spectral);
  // Extraction of stokes parameters should be exact.
  if (err > 0.0) return false;

  return true;
}

bool test_phase_matrix_copy_const_tro() {
  auto sht          = sht::provider.get_instance(1, 32);
  auto t_grid       = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid       = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid = sht->get_za_grid_ptr();
  auto phase_matrix_gridded =
      make_phase_matrix(t_grid, f_grid, za_scat_grid);
  auto phase_matrix_spectral = phase_matrix_gridded.to_spectral(sht);

  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded>
      phase_matrix_gridded_2(phase_matrix_gridded);
  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded>
      phase_matrix_gridded_3(phase_matrix_spectral);
  Numeric err = max_error(phase_matrix_gridded_2, phase_matrix_gridded_3);
  if (err > 1e-6) {
    return false;
  }

  PhaseMatrixData<Numeric, Format::TRO, Representation::Spectral>
      phase_matrix_spectral_2 = phase_matrix_gridded_2;
  err = max_error(phase_matrix_spectral.extract_stokes_coeffs(),
                  phase_matrix_spectral_2);
  if (err > 1e-6) {
    return false;
  }

  return true;
}

/** Test regridding of TRO phase matrices.
 *
 * This method ensures the regridding of phase matrix in TRO format in both
 * gridded and spectral representation yield the expected results.
 *
 * @return true if all tests passed, false otherwise.
 */
bool test_phase_matrix_regrid_tro() {
  auto sht          = sht::provider.get_instance(1, 32);
  auto t_grid       = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid       = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_scat_grid = sht->get_za_grid_ptr();
  auto phase_matrix_gridded =
      make_phase_matrix(t_grid, f_grid, za_scat_grid);
  auto phase_matrix_spectral = phase_matrix_gridded.to_spectral(sht);

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto t_grid_new       = std::make_shared<Vector>(Vector({210}));
  auto f_grid_new       = std::make_shared<Vector>(Vector({1e9}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_new = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0})));

  ScatteringDataGrids grids{t_grid_new, f_grid_new, za_scat_grid};
  auto weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid, grids);

  auto phase_matrix_gridded_interp =
      phase_matrix_gridded.regrid(grids, weights);
  Numeric err = max_error(
      static_cast<MatrixView>(phase_matrix_gridded[0, 0, joker, joker]),
      static_cast<MatrixView>(phase_matrix_gridded_interp[0, 0, joker, joker]));
  if (err > 0.0) {
    return false;
  }

  //
  // Do the same for data in spectral representation. Here, however, all
  // scattering zenith angles are extracted because there's no way to perform
  // angle interpolation in spectral space.
  //

  auto phase_matrix_spectral_interp =
      phase_matrix_spectral.regrid(grids, weights);
  phase_matrix_gridded_interp = phase_matrix_spectral_interp.to_gridded();

  err = max_error(
      static_cast<MatrixView>(phase_matrix_gridded[0, 0, joker, joker]),
      static_cast<MatrixView>(phase_matrix_gridded_interp[0, 0, joker, joker]));
  if (err > 1e-10) {
    return false;
  }

  //
  // Test interpolation at mid-point between first and second elements
  // along temperature, frequency and scattering zenith angle.
  //

  (*t_grid_new)[0]       = 230.0;
  (*f_grid_new)[0]       = 5.5e9;
  grid_vector(*za_scat_grid_new)[0] = 0.5 * (grid_vector(*za_scat_grid)[0] + (grid_vector(*za_scat_grid))[1]);
  grids   = ScatteringDataGrids{t_grid_new, f_grid_new, za_scat_grid_new};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid, grids);
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);

  Vector phase_matrix_ref(6);
  phase_matrix_ref = 0.0;
  for (Index i_t = 0; i_t < 2; ++i_t) {
    for (Index i_f = 0; i_f < 2; ++i_f) {
      for (Index i_za_scat = 0; i_za_scat < 2; ++i_za_scat) {
        phase_matrix_ref += static_cast<matpack::data_t<double, 1>>(
            0.125 * phase_matrix_gridded[i_t, i_f, i_za_scat, joker]);
      }
    }
  }
  err = max_error(
      static_cast<VectorView>(phase_matrix_ref),
      static_cast<VectorView>(phase_matrix_gridded_interp[0, 0, 0, joker]));
  if (err > 1e-15) {
    return false;
  }

  //
  // Do the same in spectral space.
  //

  auto phase_matrix_interp =
      phase_matrix_spectral.regrid(grids, weights).to_gridded();
  Matrix phase_matrix_spectral_ref(grid_size(*za_scat_grid), 6);
  phase_matrix_spectral_ref = 0.0;
  for (Index i_t = 0; i_t < 2; ++i_t) {
    for (Index i_f = 0; i_f < 2; ++i_f) {
      phase_matrix_spectral_ref +=
          static_cast<matpack::data_t<double, 2>>(
              0.25 * phase_matrix_gridded[i_t, i_f, joker, joker]);
    }
  }
  err = max_error(
      static_cast<MatrixView>(phase_matrix_spectral_ref),
      static_cast<MatrixView>(phase_matrix_interp[0, 0, joker, joker]));
  if (err > 1e-10) {
    return false;
  }

  //
  // Test interpolation for arbitrary values along axes.
  //

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_inc =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector(std::views::iota(0, grid_size(*za_scat_grid)))));

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 4> &>(
      phase_matrix_gridded));

  (*t_grid_new)[0]       = 1.2345;
  (*f_grid_new)[0]       = 1.2345;
  grid_vector(*za_scat_grid_new)[0] = 1.2345;
  grids   = ScatteringDataGrids{t_grid_new, f_grid_new, za_scat_grid_new};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid_inc, grids);
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<0>(
      reinterpret_cast<matpack::data_t<std::complex<Numeric>, 4> &>(
          phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along f-axis.

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 4> &>(
      phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<1>(
      reinterpret_cast<matpack::data_t<std::complex<Numeric>, 4> &>(
          phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along za-scat-axis.

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<Numeric, 4> &>(
      phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

bool test_backscatter_matrix_regrid_tro() {
  auto sht          = sht::provider.get_instance(1, 32);
  auto t_grid       = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid       = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_scat_grid = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(sht->get_zenith_angle_grid().angles));
  auto phase_matrix = make_phase_matrix(t_grid, f_grid, za_scat_grid);
  auto backscatter_matrix = phase_matrix.extract_backscatter_matrix();

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto t_grid_new       = std::make_shared<Vector>(Vector({210}));
  auto f_grid_new       = std::make_shared<Vector>(Vector({1e9}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_new = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0})));

  ScatteringDataGrids grids{t_grid_new, f_grid_new, za_scat_grid};
  auto weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid, grids);

  auto backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  Numeric err                    = max_error(
      static_cast<VectorView>(backscatter_matrix[0, 0, joker]),
      static_cast<VectorView>(backscatter_matrix_interp[0, 0, joker]));
  if (err > 0.0) {
    return false;
  }

  //
  // Test interpolation for arbitrary values along axes.
  //

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 3> &>(
      backscatter_matrix));

  (*t_grid_new)[0]       = 1.2345;
  (*f_grid_new)[0]       = 1.2345;
  grid_vector(*za_scat_grid_new)[0] = 1.2345;
  grids   = ScatteringDataGrids{t_grid_new, f_grid_new, za_scat_grid_new};
  weights = calc_regrid_weights(
      t_grid, f_grid, nullptr, nullptr, nullptr, nullptr, grids);
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err = std::abs(backscatter_matrix_interp[0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along f-axis.
  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 3> &>(
      backscatter_matrix));
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err = std::abs(backscatter_matrix_interp[0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }
  return true;
}

bool test_phase_matrix_aro() {
  Index l_max        = 128;
  Index m_max        = 128;
  auto sht           = sht::provider.get_instance(l_max, m_max);
  auto t_grid        = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid        = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_inc_grid   = std::make_shared<Vector>(Vector({20.0}));
  auto za_scat_grid  = sht->get_za_grid_ptr();
  auto delta_aa_grid = sht->get_aa_grid_ptr();

  auto phase_matrix_gridded = make_phase_matrix(
      t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);

  //
  // Test conversion between spectral and gridded format.
  //
  auto phase_matrix_spectral = phase_matrix_gridded.to_spectral(sht);
  ComplexVector coeffs_ref(phase_matrix_spectral.extent(3));
  for (Index i_t = 0; i_t < phase_matrix_spectral.extent(0); ++i_t) {
    for (Index i_f = 0; i_f < phase_matrix_spectral.extent(1); ++i_f) {
      Index l                                = i_t;
      Index m                                = std::min(i_t, i_f);
      coeffs_ref                             = std::complex<Numeric>(0.0, 0.0);
      coeffs_ref[sht->get_coeff_index(l, m)] = std::complex<Numeric>(1.0, 0.0);
      for (Index i_za_inc = 0; i_za_inc < phase_matrix_spectral.extent(2);
           ++i_za_inc) {
        for (Index i_s = 0; i_s < phase_matrix_spectral.extent(4); ++i_s) {
          Numeric err = max_error<ComplexVector>(
              coeffs_ref,
              static_cast<ComplexVector>(
                  phase_matrix_spectral[i_t, i_f, i_za_inc, joker, i_s]));
          if (err > 1e-6) return false;
        }
      }
    }
  }
  auto phase_matrix_gridded_2 = phase_matrix_spectral.to_gridded();
  Numeric err = max_error(phase_matrix_gridded, phase_matrix_gridded_2);
  if (err > 1e-6) {
    return false;
  }

  auto backscatter_matrix = phase_matrix_gridded.extract_backscatter_matrix();
  auto backscatter_matrix_2 =
      phase_matrix_spectral.extract_backscatter_matrix();
  err = max_error<Tensor4>(backscatter_matrix, backscatter_matrix_2);
  if (err > 1e-3) return false;

  auto forwardscatter_matrix =
      phase_matrix_gridded.extract_forwardscatter_matrix();
  auto forwardscatter_matrix_2 =
      phase_matrix_spectral.extract_forwardscatter_matrix();
  err = max_error<Tensor4>(forwardscatter_matrix, forwardscatter_matrix_2);
  if (err > 1e-3) return false;
  auto phase_matrix_gridded_1 = phase_matrix_gridded.extract_stokes_coeffs();
  err                         = max_error<Tensor6View>(
      phase_matrix_gridded_1,
      phase_matrix_gridded);
  if (err > 0) return false;

  return true;
}

/** Test regridding of ARO phase matrices.
 *
 * This method ensures the regridding of phase matrix in ARO format in both
 * gridded and spectral representation yield the expected results.
 *
 * @return true if all tests passed, false otherwise.
 */
bool test_phase_matrix_regrid_aro() {
  auto sht           = sht::provider.get_instance(1, 32);
  auto t_grid        = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid        = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_scat_grid  = sht->get_za_grid_ptr();
  auto za_inc_grid   = std::make_shared<Vector>(Vector({0.0, 20.0, 40.0}));
  auto delta_aa_grid = std::make_shared<Vector>(std::views::iota(0, 180));
  auto phase_matrix_gridded = make_phase_matrix(
      t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
  auto phase_matrix_spectral = phase_matrix_gridded.to_spectral();

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto t_grid_new       = std::make_shared<Vector>(Vector({210}));
  auto f_grid_new       = std::make_shared<Vector>(Vector({1e9}));
  auto za_inc_grid_new  = std::make_shared<Vector>(Vector{0.0});
  auto aa_scat_grid_new = std::make_shared<Vector>(Vector({0.0}));
  auto za_scat_grid_new = std::make_shared<ZenithAngleGrid>(
      IrregularZenithAngleGrid(Vector{grid_vector(*za_scat_grid)[0]}));

  ScatteringDataGrids grids{t_grid_new,
                            f_grid_new,
                            za_inc_grid_new,
                            aa_scat_grid_new,
                            za_scat_grid_new};
  auto weights = calc_regrid_weights(t_grid,
                                     f_grid,
                                     nullptr,
                                     std::make_shared<Vector>(grid_vector(*za_scat_grid)),
                                     delta_aa_grid,
                                     za_scat_grid,
                                     grids);

  auto phase_matrix_gridded_interp =
      phase_matrix_gridded.regrid(grids, weights);
  Numeric err = max_error(
      static_cast<VectorView>(phase_matrix_gridded[0, 0, 0, 0, 0, joker]),
      static_cast<VectorView>(
          phase_matrix_gridded_interp[0, 0, 0, 0, 0, joker]));
  if (err > 1e-10) {
    return false;
  }

  //
  // Do the same for data in spectral representation. Here, however, all
  // scattering angles are extracted because there's no way to perform
  // angle interpolation in spectral space.
  //

  auto phase_matrix_spectral_interp =
      phase_matrix_spectral.regrid(grids, weights);
  phase_matrix_gridded_interp = phase_matrix_spectral_interp.to_gridded();

  err = max_error(static_cast<Tensor3View>(
                      phase_matrix_gridded[0, 0, 0, joker, joker, joker]),
                  static_cast<Tensor3View>(phase_matrix_gridded_interp[
                      0, 0, 0, joker, joker, joker]));
  if (err > 1e-10) {
    return false;
  }

  //
  // Test interpolation for arbitrary values along axes.
  //
  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  std::shared_ptr<const Vector> za_inc_grid_inc = std::make_shared<Vector>(
      std::from_range_t{}, std::views::iota(Size{0}, za_inc_grid->size()));
  std::shared_ptr<const Vector> aa_scat_grid_inc = std::make_shared<Vector>(
      std::from_range_t{}, std::views::iota(Size{0}, delta_aa_grid->size()));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_inc =
      std::make_shared<ZenithAngleGrid>(
          IrregularZenithAngleGrid(Vector{std::views::iota(0, grid_size(*za_scat_grid))})
                                        );

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 6> &>(
      phase_matrix_gridded));

  (*t_grid_new)[0]       = 1.2345;
  (*f_grid_new)[0]       = 1.2345;
  (*za_inc_grid_new)[0]  = 1.2345;
  (*aa_scat_grid_new)[0] = 1.2345;
  grid_vector(*za_scat_grid_new)[0] = 1.2345;

  grids                       = ScatteringDataGrids{t_grid_new,
                              f_grid_new,
                              za_inc_grid_new,
                              aa_scat_grid_new,
                              za_scat_grid_new};
  weights                     = calc_regrid_weights(t_grid,
                                f_grid,
                                nullptr,
                                za_inc_grid_inc,
                                aa_scat_grid_inc,
                                za_scat_grid_inc,
                                grids);
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);

  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<0>(
      reinterpret_cast<matpack::data_t<std::complex<Numeric>, 5> &>(
          phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along f-axis.

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 6> &>(
      phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<1>(
      reinterpret_cast<matpack::data_t<std::complex<Numeric>, 5> &>(
          phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along za_inc axis.

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<Numeric, 6> &>(
      phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  fill_along_axis<2>(
      reinterpret_cast<matpack::data_t<std::complex<Numeric>, 5> &>(
          phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along aa_scat axis.
  fill_along_axis<3>(reinterpret_cast<matpack::data_t<Numeric, 6> &>(
      phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along aa_scat axis.
  fill_along_axis<4>(reinterpret_cast<matpack::data_t<Numeric, 6> &>(
      phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

bool test_backscatter_matrix_regrid_aro() {
  auto sht           = sht::provider.get_instance(1, 32);
  auto t_grid        = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid        = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_scat_grid  = sht->get_za_grid_ptr();
  auto za_inc_grid   = std::make_shared<Vector>(Vector({0.0, 20.0, 40.0}));
  auto delta_aa_grid = std::make_shared<Vector>(std::views::iota(0, 180));
  auto phase_matrix_gridded = make_phase_matrix(
      t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
  auto backscatter_matrix = phase_matrix_gridded.extract_backscatter_matrix();

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto t_grid_new       = std::make_shared<Vector>(Vector({210}));
  auto f_grid_new       = std::make_shared<Vector>(Vector({1e9}));
  auto za_inc_grid_new  = std::make_shared<Vector>(Vector{0.0});
  auto aa_scat_grid_new = std::make_shared<Vector>(Vector({0.0}));
  auto za_scat_grid_new = std::make_shared<ZenithAngleGrid>(
      IrregularZenithAngleGrid(Vector{grid_vector(*za_scat_grid)[0]}));

  ScatteringDataGrids grids{t_grid_new,
                            f_grid_new,
                            za_inc_grid_new,
                            aa_scat_grid_new,
                            za_scat_grid_new};
  auto weights = calc_regrid_weights(t_grid,
                                     f_grid,
                                     nullptr,
                                     std::make_shared<Vector>(grid_vector(*za_scat_grid)),
                                     delta_aa_grid,
                                     za_scat_grid,
                                     grids);

  auto backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  Numeric err                    = max_error(
      static_cast<VectorView>(backscatter_matrix[0, 0, 0, joker]),
      static_cast<VectorView>(backscatter_matrix_interp[0, 0, 0, joker]));
  if (err > 1e-10) {
    return false;
  }

  //
  // Test interpolation for arbitrary values along axes.
  //
  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  std::shared_ptr<const Vector> za_inc_grid_inc = std::make_shared<Vector>(
      std::from_range_t{}, std::views::iota(Size{0}, za_inc_grid->size()));
  std::shared_ptr<const Vector> aa_scat_grid_inc = std::make_shared<Vector>(
      std::from_range_t{}, std::views::iota(Size{0}, delta_aa_grid->size()));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_inc =
      std::make_shared<ZenithAngleGrid>(
          IrregularZenithAngleGrid(Vector{std::views::iota(0, grid_size(*za_scat_grid))}));

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 4> &>(
      backscatter_matrix));

  (*t_grid_new)[0]       = 1.2345;
  (*f_grid_new)[0]       = 1.2345;
  (*za_inc_grid_new)[0]  = 1.2345;
  (*aa_scat_grid_new)[0] = 1.2345;
  grid_vector(*za_scat_grid_new)[0] = 1.2345;

  grids                     = ScatteringDataGrids{t_grid_new,
                              f_grid_new,
                              za_inc_grid_new,
                              aa_scat_grid_new,
                              za_scat_grid_new};
  weights                   = calc_regrid_weights(t_grid,
                                f_grid,
                                nullptr,
                                za_inc_grid_inc,
                                aa_scat_grid_inc,
                                za_scat_grid_inc,
                                grids);
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err = std::abs(backscatter_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along f-axis.

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 4> &>(
      backscatter_matrix));
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err = std::abs(backscatter_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  // Test interpolation along za_inc axis.

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<Numeric, 4> &>(
      backscatter_matrix));
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err = std::abs(backscatter_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) {
    return false;
  }

  return true;
}

int main() {
#ifndef ARTS_NO_SHTNS
  bool passed = false;
  std::cout << "Testing phase matrix (TRO): ";
  passed = test_phase_matrix_tro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing phase matrix copy constructor (TRO): ";
  passed = test_phase_matrix_copy_const_tro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing phase matrix regridding (TRO): ";
  passed = test_phase_matrix_regrid_tro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing backscatter matrix regridding (TRO): ";
  passed = test_backscatter_matrix_regrid_tro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing phase matrix (ARO): ";
  passed = test_phase_matrix_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing phase matrix regridding (ARO): ";
  passed = test_phase_matrix_regrid_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing backscatter matrix regridding (ARO): ";
  passed = test_backscatter_matrix_regrid_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
#endif

  return 0;
}
