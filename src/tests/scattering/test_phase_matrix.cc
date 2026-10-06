#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <numbers>
#include <ranges>

#include "arts_conversions.h"
#include "interpolation.h"
#include "matpack/matpack.h"
#include "matpack/matpack_mdspan_helpers_eigen.h"
#include "scattering/mie.h"
#include "scattering/phase_matrix.h"
#include "test_utils.h"

using namespace scattering;

using PhaseMatrixTROGridded = PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded>;

using PhaseMatrixTROSpectral = PhaseMatrixData<Numeric, Format::TRO, Representation::Spectral>;

using PhaseMatrixAROGridded = PhaseMatrixData<Numeric, Format::ARO, Representation::Gridded>;

using PhaseMatrixAROSpectral = PhaseMatrixData<Numeric, Format::ARO, Representation::Spectral>;

using PhaseMatrixAROFourier = PhaseMatrixData<Numeric, Format::ARO, Representation::Fourier>;

/** Create a TRO phase matrix for testing
 *
 * Creates a phase matrix with Legendre polynomials with the degree
 * similar to the frequency index.
 */
PhaseMatrixTROGridded make_phase_matrix(std::shared_ptr<const Vector>          t_grid,
                                        std::shared_ptr<const Vector>          f_grid,
                                        std::shared_ptr<const ZenithAngleGrid> za_scat_grid) {
  PhaseMatrixTROGridded phase_matrix(t_grid, f_grid, za_scat_grid);
  Vector                za_grid  = Vector{grid_vector(*za_scat_grid)};
  za_grid                       *= Conversion::deg2rad(1.0);
  Vector aa_grid(1);
  aa_grid = 0.0;

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
        phase_matrix[i_t, Range(i_f, 1), joker, i_s] = evaluate_spherical_harmonic(i_f, 0, aa_grid, za_grid);
      }
    }
  }
  return phase_matrix;
}

/** Create a TRO phase matrix for a liquid sphere.
 *
 * Creates phase matrix data for a liquid sphere with a radius of 100um.
 */
PhaseMatrixTROGridded make_phase_matrix_liquid_sphere(std::shared_ptr<const Vector>          t_grid,
                                                      std::shared_ptr<const Vector>          f_grid,
                                                      std::shared_ptr<const ZenithAngleGrid> za_scat_grid) {
  PhaseMatrixTROGridded phase_matrix(t_grid, f_grid, za_scat_grid);

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      auto scat_data =
          scattering::MieSphere<Numeric>::Liquid((*f_grid)[i_f], (*t_grid)[i_t], 1e-3, grid_vector(*za_scat_grid));
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
PhaseMatrixAROGridded make_phase_matrix(std::shared_ptr<const Vector>          t_grid,
                                        std::shared_ptr<const Vector>          f_grid,
                                        std::shared_ptr<const Vector>          za_inc_grid,
                                        std::shared_ptr<const Vector>          delta_aa_grid,
                                        std::shared_ptr<const ZenithAngleGrid> za_scat_grid) {
  PhaseMatrixAROGridded phase_matrix(t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
  Vector                aa_grid  = *delta_aa_grid;
  aa_grid                       *= Conversion::deg2rad(1.0);
  Vector za_grid                 = Vector{grid_vector(*za_scat_grid)};
  za_grid                       *= Conversion::deg2rad(1.0);

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      for (Size i_za_inc = 0; i_za_inc < za_inc_grid->size(); ++i_za_inc) {
        for (Index i_s = 0; i_s < phase_matrix.n_stokes_coeffs; ++i_s) {
          Index l                                             = i_t;
          Index m                                             = std::min(i_t, i_f);
          phase_matrix[i_t, i_f, i_za_inc, joker, joker, i_s] = evaluate_spherical_harmonic(l, m, aa_grid, za_grid);
        }
      }
    }
  }
  return phase_matrix;
}

/** The TRO series with coefficient 1 at degree i_f (for every temperature and element), of degree 15
 *
 * It is the series of make_phase_matrix, exactly.
 */
PhaseMatrixTROSpectral make_unit_series(std::shared_ptr<const Vector> t_grid, std::shared_ptr<const Vector> f_grid) {
  PhaseMatrixTROSpectral series(t_grid, f_grid, 15);
  for (Size i_t = 0; i_t < t_grid->size(); ++i_t)
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f)
      for (Index i_s = 0; i_s < series.n_stokes_coeffs; ++i_s) series[i_t, i_f, i_f, i_s] = 1.0;
  return series;
}

bool test_phase_matrix_tro() {
  auto                             t_grid               = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto                             f_grid               = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid         = std::make_shared<ZenithAngleGrid>(GaussLegendreGrid(32));
  auto                             phase_matrix_gridded = make_phase_matrix(t_grid, f_grid, za_scat_grid);

  //
  // A series evaluates exactly on any grid, and the projection of gridded
  // data converges to the series they sample (second order in the spacing).
  //
  auto    phase_matrix_spectral  = make_unit_series(t_grid, f_grid);
  auto    phase_matrix_gridded_2 = phase_matrix_spectral.to_gridded(za_scat_grid);
  Numeric err                    = max_error(phase_matrix_gridded, phase_matrix_gridded_2);
  if (err > 1e-12) { return false; }

  auto fine = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(nlinspace(0.0, 180.0, 3601)));
  auto projected = make_phase_matrix(t_grid, f_grid, fine).to_spectral(15);
  err            = max_error<matpack::strided_view_t<const std::complex<Numeric>, 4>>(projected, phase_matrix_spectral);
  if (err > 1e-5) { return false; }

  //
  // Test conversion to lab frame.
  //

  std::shared_ptr<Vector>          delta_aa_grid = std::make_shared<Vector>(Vector({0.0, 180}));
  std::shared_ptr<Vector>          za_inc_grid   = std::make_shared<Vector>(Vector({90.0}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_new =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({90.0})));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_liquid =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0.0, 10, 20, 160, 180.0})));

  auto phase_matrix_liquid = make_phase_matrix_liquid_sphere(t_grid, f_grid, za_scat_grid_liquid);
  auto phase_matrix_lab    = phase_matrix_liquid.to_lab_frame(za_inc_grid, delta_aa_grid, za_scat_grid_new);

  for (Size i_t = 0; i_t < t_grid->size(); ++i_t) {
    for (Size i_f = 0; i_f < f_grid->size(); ++i_f) {
      // Backward scattering direction. Only two independent elements.
      // Off-diagonal elements must be close to 0.
      auto    pm_f        = phase_matrix_lab[i_t, i_f, 0, 0, 0, joker];
      Numeric coeff_max_f = pm_f[0];
      auto    pm_b        = phase_matrix_lab[i_t, i_f, 0, 1, 0, joker];
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
      delta = std::abs(pm_b[0 * 4 + 0] - 2.0 * pm_b[1 * 4 + 1] - pm_b[3 * 4 + 3]) / coeff_max_b;
      if (delta > 1e-3) return false;

      // And all off-diagonal elements should be zero.
      for (Index i_s1 = 0; i_s1 < 4; ++i_s1) {
        for (Index i_s2 = 0; i_s2 < 4; ++i_s2) {
          if (i_s1 != i_s2) {
            Numeric c = pm_b[i_s1 * 4 + i_s2];
            if (std::abs(c) / coeff_max_b > 1e-3) {
              return false;

              c = pm_f[i_s1 * 4 + i_s2];
              if (std::abs(c) / coeff_max_f > 1e-3) { return false; }
            }
          }
        }
      }
    }
  }

  // Test extraction of backscatter matrix.
  auto backscatter_matrix = phase_matrix_liquid.extract_backscatter_matrix();
  err = max_error<Tensor3>(backscatter_matrix, static_cast<Tensor3>(phase_matrix_liquid[joker, joker, 4, joker]));
  if (err > 1e-15) return false;

  auto forwardscatter_matrix = phase_matrix_liquid.extract_forwardscatter_matrix();
  err = max_error<Tensor3>(forwardscatter_matrix, static_cast<Tensor3>(phase_matrix_liquid[joker, joker, 0, joker]));
  if (err > 1e-15) return false;

  // The back- and forward-scatter matrices of a series are the series at 180 and 0 deg
  const auto ends = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{0.0, 180.0}));
  const auto at_ends = phase_matrix_spectral.to_gridded(ends);
  err = max_error<Tensor3>(phase_matrix_spectral.extract_backscatter_matrix(),
                           static_cast<Tensor3>(at_ends[joker, joker, 1, joker]));
  if (err > 1e-14) return false;
  err = max_error<Tensor3>(phase_matrix_spectral.extract_forwardscatter_matrix(),
                           static_cast<Tensor3>(at_ends[joker, joker, 0, joker]));
  if (err > 1e-14) return false;
  auto phase_matrix_liquid_spectral = phase_matrix_spectral;

  // Test reduction of stokes elements.
  auto phase_matrix_liquid_1 = phase_matrix_liquid.extract_stokes_coeffs();
  err                        = max_error<ConstTensor4View>(phase_matrix_liquid_1, phase_matrix_liquid);
  // Extraction of stokes parameters should be exact.
  if (err > 0.0) return false;

  auto phase_matrix_liquid_spectral_1 = phase_matrix_liquid_spectral.extract_stokes_coeffs();
  err = max_error<matpack::strided_view_t<const std::complex<Numeric>, 4>>(phase_matrix_liquid_spectral_1,
                                                                           phase_matrix_liquid_spectral);
  // Extraction of stokes parameters should be exact.
  if (err > 0.0) return false;

  return true;
}

bool test_phase_matrix_copy_const_tro() {
  auto                                   t_grid               = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto                                   f_grid               = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid         = std::make_shared<ZenithAngleGrid>(GaussLegendreGrid(32));
  auto                                   phase_matrix_gridded = make_phase_matrix(t_grid, f_grid, za_scat_grid);
  auto                                   phase_matrix_spectral = make_unit_series(t_grid, f_grid);

  // A gridded copy of a series of degree 15 is the series at the 32 Gauss-Legendre nodes
  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> phase_matrix_gridded_2(phase_matrix_gridded);
  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> phase_matrix_gridded_3(phase_matrix_spectral);
  Numeric err = max_error(phase_matrix_gridded_2, phase_matrix_gridded_3);
  if (err > 1e-12) { return false; }

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
  auto t_grid                = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid                = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid = std::make_shared<ZenithAngleGrid>(GaussLegendreGrid(32));
  auto phase_matrix_gridded  = make_phase_matrix(t_grid, f_grid, za_scat_grid);
  auto phase_matrix_spectral = make_unit_series(t_grid, f_grid);

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto                             t_grid_new = std::make_shared<Vector>(Vector({210}));
  auto                             f_grid_new = std::make_shared<Vector>(Vector({1e9}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_new =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0})));

  ScatteringDataGrids grids{t_grid_new, f_grid_new, za_scat_grid};
  auto                weights = calc_regrid_weights(t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid, grids);

  auto    phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  Numeric err = max_error(static_cast<MatrixView>(phase_matrix_gridded[0, 0, joker, joker]),
                          static_cast<MatrixView>(phase_matrix_gridded_interp[0, 0, joker, joker]));
  if (err > 0.0) { return false; }

  //
  // Do the same for data in spectral representation. Here, however, all
  // scattering zenith angles are extracted because there's no way to perform
  // angle interpolation in spectral space.
  //

  auto phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  phase_matrix_gridded_interp       = phase_matrix_spectral_interp.to_gridded();

  err = max_error(static_cast<MatrixView>(phase_matrix_gridded[0, 0, joker, joker]),
                  static_cast<MatrixView>(phase_matrix_gridded_interp[0, 0, joker, joker]));
  if (err > 1e-10) { return false; }

  //
  // Test interpolation at mid-point between first and second elements
  // along temperature, frequency and scattering zenith angle.
  //

  (*t_grid_new)[0] = 230.0;
  (*f_grid_new)[0] = 5.5e9;
  *za_scat_grid_new =
      IrregularZenithAngleGrid(Vector{0.5 * (grid_vector(*za_scat_grid)[0] + (grid_vector(*za_scat_grid))[1])});
  grids                       = ScatteringDataGrids{t_grid_new, f_grid_new, za_scat_grid_new};
  weights                     = calc_regrid_weights(t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid, grids);
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);

  Vector phase_matrix_ref(6);
  phase_matrix_ref = 0.0;
  for (Index i_t = 0; i_t < 2; ++i_t) {
    for (Index i_f = 0; i_f < 2; ++i_f) {
      for (Index i_za_scat = 0; i_za_scat < 2; ++i_za_scat) {
        phase_matrix_ref +=
            static_cast<matpack::data_t<double, 1>>(0.125 * phase_matrix_gridded[i_t, i_f, i_za_scat, joker]);
      }
    }
  }
  err = max_error(static_cast<VectorView>(phase_matrix_ref),
                  static_cast<VectorView>(phase_matrix_gridded_interp[0, 0, 0, joker]));
  if (err > 1e-15) { return false; }

  //
  // Do the same in spectral space.
  //

  auto   phase_matrix_interp = phase_matrix_spectral.regrid(grids, weights).to_gridded();
  Matrix phase_matrix_spectral_ref(grid_size(*za_scat_grid), 6);
  phase_matrix_spectral_ref = 0.0;
  for (Index i_t = 0; i_t < 2; ++i_t) {
    for (Index i_f = 0; i_f < 2; ++i_f) {
      phase_matrix_spectral_ref +=
          static_cast<matpack::data_t<double, 2>>(0.25 * phase_matrix_gridded[i_t, i_f, joker, joker]);
    }
  }
  err = max_error(static_cast<MatrixView>(phase_matrix_spectral_ref),
                  static_cast<MatrixView>(phase_matrix_interp[0, 0, joker, joker]));
  if (err > 1e-10) { return false; }

  //
  // Test interpolation for arbitrary values along axes.
  //

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_inc =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector(stdv::iota(0, grid_size(*za_scat_grid)))));

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 4>&>(phase_matrix_gridded));

  (*t_grid_new)[0]            = 1.2345;
  (*f_grid_new)[0]            = 1.2345;
  *za_scat_grid_new           = IrregularZenithAngleGrid(Vector{1.2345});
  grids                       = ScatteringDataGrids{t_grid_new, f_grid_new, za_scat_grid_new};
  weights                     = calc_regrid_weights(t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid_inc, grids);
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<std::complex<Numeric>, 4>&>(phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err                          = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along f-axis.

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 4>&>(phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<std::complex<Numeric>, 4>&>(phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err                          = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along za-scat-axis.

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<Numeric, 4>&>(phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  return true;
}

bool test_backscatter_matrix_regrid_tro() {
  auto t_grid = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_scat_grid =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(GaussLegendreGrid(32).angles.vec()));
  auto phase_matrix       = make_phase_matrix(t_grid, f_grid, za_scat_grid);
  auto backscatter_matrix = phase_matrix.extract_backscatter_matrix();

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto                             t_grid_new = std::make_shared<Vector>(Vector({210}));
  auto                             f_grid_new = std::make_shared<Vector>(Vector({1e9}));
  std::shared_ptr<ZenithAngleGrid> za_scat_grid_new =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0})));

  ScatteringDataGrids grids{t_grid_new, f_grid_new, za_scat_grid};
  auto                weights = calc_regrid_weights(t_grid, f_grid, nullptr, nullptr, nullptr, za_scat_grid, grids);

  auto    backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  Numeric err                       = max_error(static_cast<VectorView>(backscatter_matrix[0, 0, joker]),
                                                static_cast<VectorView>(backscatter_matrix_interp[0, 0, joker]));
  if (err > 0.0) { return false; }

  //
  // Test interpolation for arbitrary values along axes.
  //

  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 3>&>(backscatter_matrix));

  (*t_grid_new)[0]          = 1.2345;
  (*f_grid_new)[0]          = 1.2345;
  *za_scat_grid_new         = IrregularZenithAngleGrid(Vector{1.2345});
  grids                     = ScatteringDataGrids{t_grid_new, f_grid_new, za_scat_grid_new};
  weights                   = calc_regrid_weights(t_grid, f_grid, nullptr, nullptr, nullptr, nullptr, grids);
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err                       = std::abs(backscatter_matrix_interp[0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along f-axis.
  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 3>&>(backscatter_matrix));
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err                       = std::abs(backscatter_matrix_interp[0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }
  return true;
}

bool test_phase_matrix_aro() {
  Index l_max         = 128;
  Index m_max         = 128;
  auto  sht           = sht::provider.get_instance(l_max, m_max);
  auto  t_grid        = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto  f_grid        = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto  za_inc_grid   = std::make_shared<Vector>(Vector({20.0}));
  auto  za_scat_grid  = sht->get_za_grid_ptr();
  auto  delta_aa_grid = sht->get_aa_grid_ptr();

  auto phase_matrix_gridded = make_phase_matrix(t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);

  //
  // Test conversion between spectral and gridded format.
  //
  auto          phase_matrix_spectral = phase_matrix_gridded.to_spectral(sht);
  ComplexVector coeffs_ref(phase_matrix_spectral.extent(3));
  for (Index i_t = 0; i_t < phase_matrix_spectral.extent(0); ++i_t) {
    for (Index i_f = 0; i_f < phase_matrix_spectral.extent(1); ++i_f) {
      Index l                                = i_t;
      Index m                                = std::min(i_t, i_f);
      coeffs_ref                             = std::complex<Numeric>(0.0, 0.0);
      coeffs_ref[sht->get_coeff_index(l, m)] = std::complex<Numeric>(1.0, 0.0);
      for (Index i_za_inc = 0; i_za_inc < phase_matrix_spectral.extent(2); ++i_za_inc) {
        for (Index i_s = 0; i_s < phase_matrix_spectral.extent(4); ++i_s) {
          Numeric err = max_error<ComplexVector>(
              coeffs_ref, static_cast<ComplexVector>(phase_matrix_spectral[i_t, i_f, i_za_inc, joker, i_s]));
          if (err > 1e-6) return false;
        }
      }
    }
  }
  auto    phase_matrix_gridded_2 = phase_matrix_spectral.to_gridded();
  Numeric err                    = max_error(phase_matrix_gridded, phase_matrix_gridded_2);
  if (err > 1e-6) { return false; }

  auto backscatter_matrix   = phase_matrix_gridded.extract_backscatter_matrix();
  auto backscatter_matrix_2 = phase_matrix_spectral.extract_backscatter_matrix();
  err                       = max_error<Tensor4>(backscatter_matrix, backscatter_matrix_2);
  if (err > 1e-3) return false;

  auto forwardscatter_matrix   = phase_matrix_gridded.extract_forwardscatter_matrix();
  auto forwardscatter_matrix_2 = phase_matrix_spectral.extract_forwardscatter_matrix();
  err                          = max_error<Tensor4>(forwardscatter_matrix, forwardscatter_matrix_2);
  if (err > 1e-3) return false;
  auto phase_matrix_gridded_1 = phase_matrix_gridded.extract_stokes_coeffs();
  err                         = max_error<Tensor6View>(phase_matrix_gridded_1, phase_matrix_gridded);
  if (err > 0) return false;

  return true;
}

/** ARO azimuthal Fourier modes.
 *
 * Gridded data linear in the azimuth difference are their own interpolant,
 * so their modes are exact: f = Delta [rad] on [-pi, pi] has C_m = 0 and
 * S_m = 2 (-1)^(m + 1) / m.  A series evaluates exactly on any grid, and the
 * modes of finely gridded samples converge to it.
 */
bool test_phase_matrix_aro_fourier() {
  auto t_grid       = std::make_shared<Vector>(Vector({210.0, 250.0}));
  auto f_grid       = std::make_shared<Vector>(Vector({1e9}));
  auto za_inc_grid  = std::make_shared<Vector>(Vector({20.0, 140.0}));
  auto za_scat_grid = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({10.0, 90.0})));
  const Index M     = 6;

  auto          two = std::make_shared<Vector>(Vector({-180.0, 180.0}));
  PhaseMatrixAROGridded linear(t_grid, f_grid, za_inc_grid, two, za_scat_grid);
  for (Index e = 0; e < 16; ++e) {
    linear[joker, joker, joker, 0, joker, e] = -std::numbers::pi;
    linear[joker, joker, joker, 1, joker, e] = std::numbers::pi;
  }
  auto    modes = linear.to_fourier(M);
  Numeric err   = 0.0;
  for (Index m = 0; m <= M; ++m) {
    const Numeric S = m == 0 ? 0.0 : 2.0 * (m % 2 == 1 ? 1.0 : -1.0) / static_cast<Numeric>(m);
    err             = std::max(err, std::abs(modes[1, 0, 1, 0, m, 0, 5]));
    err             = std::max(err, std::abs(modes[1, 0, 1, 0, m, 1, 5] - S));
  }
  if (err > 1e-13) return false;
  // The scattering zenith angles 10 and 90 deg do not cover all directions: no phase integral
  if (not std::isnan(modes.get_phase_integral()[0, 0, 0])) return false;

  // The phase integral of Z11 = 1 and Z11 = za [rad], linear in za between nodes spanning [0, 180] deg:
  // 2 pi int sin = 4 pi and 2 pi int x sin(x) = 2 pi^2, exactly
  auto full = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0.0, 60.0, 180.0})));
  PhaseMatrixAROGridded flat(t_grid, f_grid, za_inc_grid, two, full);
  for (Index i_s = 0; i_s < 3; ++i_s) {
    flat[0, joker, joker, joker, i_s, 0] = 1.0;
    flat[1, joker, joker, joker, i_s, 0] = Conversion::deg2rad(grid_vector(*full)[i_s]);
  }
  const auto flat_modes = flat.to_fourier(M);
  if (std::abs(flat_modes.get_phase_integral()[0, 0, 1] - 4.0 * std::numbers::pi) > 1e-13) return false;
  if (std::abs(flat_modes.get_phase_integral()[1, 0, 0] - 2.0 * std::numbers::pi * std::numbers::pi) > 1e-13)
    return false;

  // Data that do not span one period cannot give the modes
  try {
    PhaseMatrixAROGridded half(t_grid, f_grid, za_inc_grid, std::make_shared<Vector>(Vector({0.0, 180.0})), za_scat_grid);
    (void)half.to_fourier(M);
    return false;
  } catch (const std::exception&) {
  }

  // A series, evaluated finely and back
  PhaseMatrixAROFourier series(t_grid, f_grid, za_inc_grid, za_scat_grid, M);
  for (Index m = 0; m <= M; ++m) {
    series[joker, joker, joker, joker, m, 0, joker] = 1.0 / (1.0 + static_cast<Numeric>(m));
    if (m > 0) series[joker, joker, joker, joker, m, 1, joker] = 0.5 / static_cast<Numeric>(m);
  }
  auto fine = std::make_shared<Vector>(nlinspace(-180.0, 180.0, 3601));
  auto back = series.to_gridded(fine).to_fourier(M);
  err       = max_error<matpack::strided_view_t<const Numeric, 7>>(back, series);
  if (err > 1e-4) return false;

  // Truncation, and its limit
  if (series.to_fourier(2).get_max_mode() != 2) return false;
  try {
    (void)series.to_fourier(M + 1);
    return false;
  } catch (const std::exception&) {
  }
  return true;
}

/** The Fourier modes of a TRO laboratory-frame phase matrix against tro_lab_frame.
 *
 * The Rayleigh scattering matrix is a polynomial in cos(Theta) that vanishes
 * where a physical one must, so its laboratory-frame phase matrix is a
 * trigonometric polynomial in the azimuth difference of degree 2: the modes
 * from the converged trapezoidal rule must reproduce tro_lab_frame at any
 * azimuth to rounding, and the modes above 2 must vanish.  The directions
 * include the zenith, the nadir and the forward direction (za_inc = za_scat).
 */
bool test_tro_lab_frame_fourier_modes() {
  auto       t_grid       = std::make_shared<Vector>(Vector({250.0}));
  auto       f_grid       = std::make_shared<Vector>(Vector({1e9, 2e9}));
  auto       za_inc_grid  = std::make_shared<Vector>(Vector({0.0, 35.0, 120.0, 180.0}));
  auto       za_scat_grid = std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector({0.0, 35.0, 77.0, 180.0})));
  const auto rayleigh     = [](Numeric theta, matpack::data_t<Numeric, 3>& F) {
    const Numeric c = std::cos(theta);
    for (Index i_f = 0; i_f < F.extent(1); ++i_f) {
      const Numeric s = 1.0 + static_cast<Numeric>(i_f);
      F[0, i_f, 0]    = s * 0.75 * (1.0 + c * c);
      F[0, i_f, 1]    = s * -0.75 * (1.0 - c * c);
      F[0, i_f, 2]    = s * 0.75 * (1.0 + c * c);
      F[0, i_f, 3]    = s * 1.5 * c;
      F[0, i_f, 4]    = 0.0;
      F[0, i_f, 5]    = s * 1.5 * c;
    }
  };
  const Index M     = 5;
  Matrix      integral(1, 2);  // 2 pi int F11 dcos(Theta) = 4 pi s
  integral[0, 0] = 4.0 * std::numbers::pi;
  integral[0, 1] = 8.0 * std::numbers::pi;
  auto modes = tro_lab_frame_fourier_modes<Numeric>(t_grid, f_grid, za_inc_grid, za_scat_grid, M, integral, rayleigh);
  for (Index i = 0; i < 4; ++i)
    if (modes.get_phase_integral()[0, 0, i] != integral[0, 0] or modes.get_phase_integral()[0, 1, i] != integral[0, 1])
      return false;
  auto        delta = std::make_shared<Vector>(Vector({0.0, 13.0, 90.0, 179.0, 180.0, 181.0, 271.0, 359.0}));
  auto        lab   = tro_lab_frame<Numeric>(t_grid, f_grid, za_inc_grid, delta, za_scat_grid, rayleigh);
  Numeric     err   = max_error<matpack::strided_view_t<const Numeric, 6>>(modes.to_gridded(delta), lab);
  if (err > 1e-12) return false;
  // Azimuth differences of any period: -90 and -1 deg are 270 and 359 deg
  auto negative = std::make_shared<Vector>(Vector({-90.0, -1.0}));
  auto positive = std::make_shared<Vector>(Vector({270.0, 359.0}));
  err           = max_error<matpack::strided_view_t<const Numeric, 6>>(
      tro_lab_frame<Numeric>(t_grid, f_grid, za_inc_grid, negative, za_scat_grid, rayleigh),
      tro_lab_frame<Numeric>(t_grid, f_grid, za_inc_grid, positive, za_scat_grid, rayleigh));
  if (err > 1e-12) return false;

  Numeric high = 0.0;
  for (Index m = 3; m <= M; ++m)
    for (Numeric v : modes[joker, joker, joker, joker, m, joker, joker] | by_elem) high = std::max(high, std::abs(v));
  return high < 1e-13;
}

/** Test regridding of ARO phase matrices.
 *
 * This method ensures the regridding of phase matrix in ARO format in both
 * gridded and spectral representation yield the expected results.
 *
 * @return true if all tests passed, false otherwise.
 */
bool test_phase_matrix_regrid_aro() {
  auto sht                   = sht::provider.get_instance(1, 32);
  auto t_grid                = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid                = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  auto za_scat_grid          = sht->get_za_grid_ptr();
  auto za_inc_grid           = std::make_shared<Vector>(Vector({0.0, 20.0, 40.0}));
  auto delta_aa_grid         = std::make_shared<Vector>(stdv::iota(0, 180));
  auto phase_matrix_gridded  = make_phase_matrix(t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
  auto phase_matrix_spectral = phase_matrix_gridded.to_spectral();

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto t_grid_new       = std::make_shared<Vector>(Vector({210}));
  auto f_grid_new       = std::make_shared<Vector>(Vector({1e9}));
  auto za_inc_grid_new  = std::make_shared<Vector>(Vector{0.0});
  auto aa_scat_grid_new = std::make_shared<Vector>(Vector({0.0}));
  auto za_scat_grid_new =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{grid_vector(*za_scat_grid)[0]}));

  ScatteringDataGrids grids{t_grid_new, f_grid_new, za_inc_grid_new, aa_scat_grid_new, za_scat_grid_new};
  auto                weights = calc_regrid_weights(t_grid,
                                                    f_grid,
                                                    nullptr,
                                                    std::make_shared<Vector>(grid_vector(*za_scat_grid)),
                                                    delta_aa_grid,
                                                    za_scat_grid,
                                                    grids);

  auto    phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  Numeric err = max_error(static_cast<VectorView>(phase_matrix_gridded[0, 0, 0, 0, 0, joker]),
                          static_cast<VectorView>(phase_matrix_gridded_interp[0, 0, 0, 0, 0, joker]));
  if (err > 1e-10) { return false; }

  //
  // Do the same for data in spectral representation. Here, however, all
  // scattering angles are extracted because there's no way to perform
  // angle interpolation in spectral space.
  //

  auto phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  phase_matrix_gridded_interp       = phase_matrix_spectral_interp.to_gridded();

  err = max_error(static_cast<Tensor3View>(phase_matrix_gridded[0, 0, 0, joker, joker, joker]),
                  static_cast<Tensor3View>(phase_matrix_gridded_interp[0, 0, 0, joker, joker, joker]));
  if (err > 1e-10) { return false; }

  //
  // Test interpolation for arbitrary values along axes.
  //
  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  std::shared_ptr<const Vector> za_inc_grid_inc =
      std::make_shared<Vector>(std::from_range_t{}, stdv::iota(Size{0}, za_inc_grid->size()));
  std::shared_ptr<const Vector> aa_scat_grid_inc =
      std::make_shared<Vector>(std::from_range_t{}, stdv::iota(Size{0}, delta_aa_grid->size()));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_inc =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{stdv::iota(0, grid_size(*za_scat_grid))}));

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 6>&>(phase_matrix_gridded));

  (*t_grid_new)[0]       = 1.2345;
  (*f_grid_new)[0]       = 1.2345;
  (*za_inc_grid_new)[0]  = 1.2345;
  (*aa_scat_grid_new)[0] = 1.2345;
  *za_scat_grid_new      = IrregularZenithAngleGrid(Vector{1.2345});

  grids   = ScatteringDataGrids{t_grid_new, f_grid_new, za_inc_grid_new, aa_scat_grid_new, za_scat_grid_new};
  weights = calc_regrid_weights(t_grid, f_grid, nullptr, za_inc_grid_inc, aa_scat_grid_inc, za_scat_grid_inc, grids);
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);

  if (err > 1e-10) { return false; }

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<std::complex<Numeric>, 5>&>(phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err                          = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along f-axis.

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 6>&>(phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<std::complex<Numeric>, 5>&>(phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err                          = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along za_inc axis.

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<Numeric, 6>&>(phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<std::complex<Numeric>, 5>&>(phase_matrix_spectral));
  phase_matrix_spectral_interp = phase_matrix_spectral.regrid(grids, weights);
  err                          = std::abs(phase_matrix_spectral_interp[0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along aa_scat axis.
  fill_along_axis<3>(reinterpret_cast<matpack::data_t<Numeric, 6>&>(phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along aa_scat axis.
  fill_along_axis<4>(reinterpret_cast<matpack::data_t<Numeric, 6>&>(phase_matrix_gridded));
  phase_matrix_gridded_interp = phase_matrix_gridded.regrid(grids, weights);
  err                         = std::abs(phase_matrix_gridded_interp[0, 0, 0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  return true;
}

bool test_backscatter_matrix_regrid_aro() {
  auto t_grid               = std::make_shared<Vector>(Vector({210.0, 250.0, 270.0}));
  auto f_grid               = std::make_shared<Vector>(Vector({1e9, 10e9, 100e9}));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid = std::make_shared<ZenithAngleGrid>(GaussLegendreGrid(32));
  auto za_inc_grid          = std::make_shared<Vector>(Vector({0.0, 20.0, 40.0}));
  auto delta_aa_grid        = std::make_shared<Vector>(stdv::iota(0, 180));
  auto phase_matrix_gridded = make_phase_matrix(t_grid, f_grid, za_inc_grid, delta_aa_grid, za_scat_grid);
  auto backscatter_matrix   = phase_matrix_gridded.extract_backscatter_matrix();

  //
  // First test: Extract element at lowest temp, freq and za_scat angle.
  //

  auto t_grid_new       = std::make_shared<Vector>(Vector({210}));
  auto f_grid_new       = std::make_shared<Vector>(Vector({1e9}));
  auto za_inc_grid_new  = std::make_shared<Vector>(Vector{0.0});
  auto aa_scat_grid_new = std::make_shared<Vector>(Vector({0.0}));
  auto za_scat_grid_new =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{grid_vector(*za_scat_grid)[0]}));

  ScatteringDataGrids grids{t_grid_new, f_grid_new, za_inc_grid_new, aa_scat_grid_new, za_scat_grid_new};
  auto                weights = calc_regrid_weights(t_grid,
                                                    f_grid,
                                                    nullptr,
                                                    std::make_shared<Vector>(grid_vector(*za_scat_grid)),
                                                    delta_aa_grid,
                                                    za_scat_grid,
                                                    grids);

  auto    backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  Numeric err                       = max_error(static_cast<VectorView>(backscatter_matrix[0, 0, 0, joker]),
                                                static_cast<VectorView>(backscatter_matrix_interp[0, 0, 0, joker]));
  if (err > 1e-10) { return false; }

  //
  // Test interpolation for arbitrary values along axes.
  //
  fill_along_axis<0>(*t_grid);
  fill_along_axis<0>(*f_grid);
  std::shared_ptr<const Vector> za_inc_grid_inc =
      std::make_shared<Vector>(std::from_range_t{}, stdv::iota(Size{0}, za_inc_grid->size()));
  std::shared_ptr<const Vector> aa_scat_grid_inc =
      std::make_shared<Vector>(std::from_range_t{}, stdv::iota(Size{0}, delta_aa_grid->size()));
  std::shared_ptr<const ZenithAngleGrid> za_scat_grid_inc =
      std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{stdv::iota(0, grid_size(*za_scat_grid))}));

  // Test interpolation along temperature axis.

  fill_along_axis<0>(reinterpret_cast<matpack::data_t<Numeric, 4>&>(backscatter_matrix));

  (*t_grid_new)[0]       = 1.2345;
  (*f_grid_new)[0]       = 1.2345;
  (*za_inc_grid_new)[0]  = 1.2345;
  (*aa_scat_grid_new)[0] = 1.2345;
  *za_scat_grid_new      = IrregularZenithAngleGrid(Vector{1.2345});

  grids   = ScatteringDataGrids{t_grid_new, f_grid_new, za_inc_grid_new, aa_scat_grid_new, za_scat_grid_new};
  weights = calc_regrid_weights(t_grid, f_grid, nullptr, za_inc_grid_inc, aa_scat_grid_inc, za_scat_grid_inc, grids);
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err                       = std::abs(backscatter_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along f-axis.

  fill_along_axis<1>(reinterpret_cast<matpack::data_t<Numeric, 4>&>(backscatter_matrix));
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err                       = std::abs(backscatter_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  // Test interpolation along za_inc axis.

  fill_along_axis<2>(reinterpret_cast<matpack::data_t<Numeric, 4>&>(backscatter_matrix));
  backscatter_matrix_interp = backscatter_matrix.regrid(grids, weights);
  err                       = std::abs(backscatter_matrix_interp[0, 0, 0, 0] - 1.2345);
  if (err > 1e-10) { return false; }

  return true;
}

/** Forward and backward scattering in the laboratory frame.
 *
 * Exactly forward (Theta = 0) or backward (Theta = pi) the scattering plane
 * is undefined and rotation_coefficients() has dedicated branches.  The
 * reference is the general branch itself: approaching the exact direction
 * from several directions (along the azimuth from both sides, along the
 * zenith in the principal plane, and obliquely), with a physical scattering
 * matrix (F12 = F34 = 0, F22 = F33 forward, F22 = -F33 backward; constant in
 * Theta, so that only the rotations matter), Z must converge to F, and the
 * exact direction must give that limit.  In particular forward-scattered Q
 * stays Q: Z22 = F22.  The general branch is evaluated at Theta = 0.05 and
 * 0.005 rad, and its distance from F must shrink at least like Theta (the
 * meridional planes of the two rays turn by O(Theta)).
 */
bool test_forward_backward_limit() {
  const auto expand = [](const Vector& f) {
    Matrix F(4, 4, 0.0);
    F[0, 0] = f[0];
    F[0, 1] = F[1, 0] = f[1];
    F[1, 1]           = f[2];
    F[2, 2]           = f[3];
    F[2, 3]           = f[4];
    F[3, 2]           = -f[4];
    F[3, 3]           = f[5];
    return F;
  };
  const auto lab = [](const Vector& f, Numeric za_inc, Numeric delta_aa, Numeric za_scat) {
    Vector     z(16);
    const auto rc = detail::rotation_coefficients<Numeric>(0.0, za_inc, delta_aa, za_scat);
    detail::expand_and_transform<Numeric>(z, rtepack::compact_planar_muelmat{f}, rc, delta_aa > 180.0);
    Matrix Z(4, 4);
    for (Index i = 0; i < 4; i++)
      for (Index j = 0; j < 4; j++) Z[i, j] = z[4 * i + j];
    return Z;
  };
  const auto distance = [](const Matrix& a, const Matrix& b) {
    Numeric d = 0.0;
    for (Index i = 0; i < 4; i++)
      for (Index j = 0; j < 4; j++) d = std::max(d, std::abs(a[i, j] - b[i, j]));
    return d;
  };

  const Vector forward{1.0, 0.0, 0.8, 0.8, 0.0, 0.6}, backward{1.0, 0.0, 0.8, -0.8, 0.0, -0.6};
  for (Numeric za : {30.0, 60.0, 90.0, 135.0}) {
    // Exact directions
    if (distance(lab(forward, za, 0.0, za), expand(forward)) > 1e-12) return false;
    if (distance(lab(backward, za, 180.0, 180.0 - za), expand(backward)) > 1e-12) return false;

    // Approaches, as (delta_aa, za_scat) for a scattering angle of about t rad
    for (Numeric t : {0.05}) {
      const Numeric                               s   = std::sin(Conversion::deg2rad(za));
      const Numeric                               deg = Conversion::rad2deg(t);
      const std::array<std::array<Numeric, 2>, 4> to_forward{
          {{deg / s, za}, {360.0 - deg / s, za}, {0.0, za + deg}, {0.6 * deg / s, za - 0.8 * deg}}};
      const std::array<std::array<Numeric, 2>, 4> to_backward{{{180.0 + deg / s, 180.0 - za},
                                                               {180.0 - deg / s, 180.0 - za},
                                                               {180.0, 180.0 - za + deg},
                                                               {180.0 + 0.6 * deg / s, 180.0 - za - 0.8 * deg}}};
      for (const auto& [f, approach] : {std::pair{&forward, &to_forward}, std::pair{&backward, &to_backward}}) {
        for (const auto& [daa, zs] : *approach) {
          const auto    shrink = [&](Numeric x, Numeric centre) { return centre + 0.1 * (x - centre); };
          const Numeric centre = f == &forward ? (daa > 180.0 ? 360.0 : 0.0) : 180.0;
          const Numeric far    = distance(lab(*f, za, daa, zs), expand(*f));
          const Numeric near =
              distance(lab(*f, za, shrink(daa, centre), shrink(zs, f == &forward ? za : 180.0 - za)), expand(*f));
          if (far > 0.5 or near > 0.15 * far + 1e-12) return false;
        }
      }
    }
  }
  return true;
}

/** Rays at the poles (za = 0 or 180 deg) in the laboratory frame.
 *
 * At a pole the meridional basis of a ray is that of its azimuth, the
 * limit along its meridian, and rotation_coefficients() has dedicated
 * branches.  The reference is the limit along the meridian: moving one
 * pole ray at a time off its pole (za = 0 + e or 180 - e at the same
 * azimuth), Z must converge to the pole value, at least like e, for every
 * delta_aa (both sides of 180 deg).  One pole ray takes any F (F12, F34 !=
 * 0); two pole rays scatter exactly forward or backward and take a
 * physical F, for which the two meridional bases differ by delta_aa.  e is
 * 1 and 0.1 deg, above the 1e-6 rad at which a ray snaps to its pole and
 * the 1e-6 rad scattering angle below which a pair is exactly forward or
 * backward.
 */
bool test_pole_limit() {
  const auto lab = [](const Vector& f, Numeric za_inc, Numeric delta_aa, Numeric za_scat) {
    Vector     z(16);
    const auto rc = detail::rotation_coefficients<Numeric>(0.0, za_inc, delta_aa, za_scat);
    detail::expand_and_transform<Numeric>(z, rtepack::compact_planar_muelmat{f}, rc, delta_aa > 180.0);
    return z;
  };
  const auto distance = [](const Vector& a, const Vector& b) {
    Numeric d = 0.0;
    for (Index i = 0; i < 16; i++) d = std::max(d, std::abs(a[i] - b[i]));
    return d;
  };
  //! A pole za moved by e along its meridian
  const auto off     = [](Numeric za, Numeric e) { return za == 0.0 ? e : 180.0 - e; };
  const auto is_pole = [](Numeric za) { return za == 0.0 or za == 180.0; };

  const Vector generic{1.0, -0.3, 0.7, 0.5, 0.2, 0.4};
  const Vector forward{1.0, 0.0, 0.8, 0.8, 0.0, 0.6}, backward{1.0, 0.0, 0.8, -0.8, 0.0, -0.6};
  const std::array<std::tuple<Numeric, Numeric, const Vector*>, 8> cases{{{0.0, 50.0, &generic},
                                                                          {180.0, 50.0, &generic},
                                                                          {50.0, 0.0, &generic},
                                                                          {50.0, 180.0, &generic},
                                                                          {0.0, 0.0, &forward},
                                                                          {180.0, 180.0, &forward},
                                                                          {0.0, 180.0, &backward},
                                                                          {180.0, 0.0, &backward}}};
  for (const auto& [za_inc, za_scat, f] : cases) {
    for (Numeric daa : {10.0, 60.0, 120.0, 170.0, 190.0, 250.0, 300.0, 350.0}) {
      const Vector pole = lab(*f, za_inc, daa, za_scat);
      for (bool move_inc : {true, false}) {
        if (not is_pole(move_inc ? za_inc : za_scat)) continue;
        const auto moved = [&](Numeric e) {
          return move_inc ? lab(*f, off(za_inc, e), daa, za_scat) : lab(*f, za_inc, daa, off(za_scat, e));
        };
        const Numeric far = distance(moved(1.0), pole), near = distance(moved(0.1), pole);
        if (far > 0.1 or near > 0.15 * far + 1e-12) return false;
      }
    }
  }
  return true;
}

int main() {
  std::cout << "Testing rays at the poles in the laboratory frame: ";
  if (test_pole_limit()) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing forward and backward scattering in the laboratory frame: ";
  if (test_forward_backward_limit()) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

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

#ifndef ARTS_NO_SHTNS
  std::cout << "Testing phase matrix (ARO): ";
  passed = test_phase_matrix_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
#endif

  std::cout << "Testing azimuthal Fourier modes (ARO): ";
  passed = test_phase_matrix_aro_fourier();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  std::cout << "Testing laboratory-frame Fourier modes (TRO): ";
  passed = test_tro_lab_frame_fourier_modes();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

#ifndef ARTS_NO_SHTNS
  std::cout << "Testing phase matrix regridding (ARO): ";
  passed = test_phase_matrix_regrid_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }
#endif

  std::cout << "Testing backscatter matrix regridding (ARO): ";
  passed = test_backscatter_matrix_regrid_aro();
  if (passed) {
    std::cout << "PASSED." << '\n';
  } else {
    std::cout << "FAILED." << '\n';
    return 1;
  }

  return 0;
}
