/*===========================================================================
  === File description
  ===========================================================================*/

/*!
  \file   m_oem.cc
  \author Patrick Eriksson <patrick.eriksson@chalmers.se>
  \date   2015-09-08

  \brief  Workspace functions related to making OEM inversions.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <array.h>
#include <atm.h>
#include <config.h>
#include <debug.h>
#include <jacobian.h>
#include <oem_settings.h>
#include <workspace.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef OEM_SUPPORT
#include "oem.h"
#endif

namespace {
enum class OEMAlgorithm : char { Linear, GaussNewton, LevenbergMarquardt };

struct OEMMethod {
  OEMAlgorithm algorithm;
  bool         conjugate_gradient = false;
  bool         measurement_space  = false;

  bool linear() const { return algorithm == OEMAlgorithm::Linear; }
  bool damped() const { return algorithm == OEMAlgorithm::LevenbergMarquardt; }
};

OEMMethod parse_oem_method(const String& method) {
  if (method == "li") return {.algorithm = OEMAlgorithm::Linear};
  if (method == "li_cg") return {.algorithm = OEMAlgorithm::Linear, .conjugate_gradient = true};
  if (method == "li_cg_m")
    return {.algorithm = OEMAlgorithm::Linear, .conjugate_gradient = true, .measurement_space = true};
  if (method == "gn") return {.algorithm = OEMAlgorithm::GaussNewton};
  if (method == "gn_cg") return {.algorithm = OEMAlgorithm::GaussNewton, .conjugate_gradient = true};
  if (method == "gn_cg_m")
    return {.algorithm = OEMAlgorithm::GaussNewton, .conjugate_gradient = true, .measurement_space = true};
  if (method == "lm" || method == "ml") return {.algorithm = OEMAlgorithm::LevenbergMarquardt};
  if (method == "lm_cg" || method == "ml_cg")
    return {.algorithm = OEMAlgorithm::LevenbergMarquardt, .conjugate_gradient = true};
  if (method == "li_m") return {.algorithm = OEMAlgorithm::Linear, .measurement_space = true};
  if (method == "gn_m") return {.algorithm = OEMAlgorithm::GaussNewton, .measurement_space = true};
  ARTS_USER_ERROR(
      "Unknown OEM method '{}'. Supported methods: li, li_m, li_cg, li_cg_m, gn, gn_m, gn_cg, gn_cg_m, lm, lm_cg; aliases: ml, ml_cg.",
      method)
}

// Both solvers consume the same validated, named damping settings.
template <typename Optimizer> void configure_lm(Optimizer&                        optimizer,
                                                const LevenbergMarquardtSettings& settings,
                                                Numeric                           tolerance,
                                                unsigned int                      iterations) {
  optimizer.set_tolerance(tolerance);
  optimizer.set_maximum_iterations(iterations);
  optimizer.set_lambda(settings.initial_damping);
  optimizer.set_lambda_decrease(settings.decrease_factor);
  optimizer.set_lambda_increase(settings.increase_factor);
  optimizer.set_lambda_maximum(settings.maximum_damping);
  optimizer.set_lambda_threshold(settings.damping_threshold);
  optimizer.set_lambda_constraint(settings.convergence_damping_limit);
}

// Validation must not run the user's forward model or invert a covariance.
void check_oem_inputs(const Vector&           x,
                      const Vector&           yf,
                      const Matrix&           jacobian,
                      const Vector&           xa,
                      const CovarianceMatrix& covmat_sx,
                      const Vector&           y,
                      const CovarianceMatrix& covmat_se,
                      const OEMMethod&        method,
                      const Vector&           normalization,
                      Index                   max_iter,
                      Numeric                 stop_dx,
                      Numeric                 max_start_cost,
                      Index                   clear_matrices,
                      Index                   display_progress) {
  const Size n = xa.size();
  const Size m = y.size();
  ARTS_USER_ERROR_IF(stdr::any_of(y, [](Numeric value) { return !std::isfinite(value); }),
                     "measurement_vec values must be finite.")

  ARTS_USER_ERROR_IF(
      (x.size() != n) && (x.size() != 0),
      "The length of *model_state_vec* must be either the same as *model_state_vec_apriori* or 0. x.size(): {}, xa.size(): {}",
      x.size(),
      xa.size());
  ARTS_USER_ERROR_IF(covmat_sx.ncols() != covmat_sx.nrows(),
                     "*model_state_covmat* must be a square matrix. covmat_sx shape: {}x{}",
                     covmat_sx.nrows(),
                     covmat_sx.ncols());
  ARTS_USER_ERROR_IF(
      static_cast<Size>(covmat_sx.ncols()) != n,
      "Inconsistency in size between *model_state_vec* and *model_state_covmat*. x.size(): {}, covmat_sx size: {}x{}",
      x.size(),
      covmat_sx.nrows(),
      covmat_sx.ncols());
  ARTS_USER_ERROR_IF(
      (yf.size() != m) && (yf.size() != 0),
      "The length of *measurement_vec_fit* must be either the same as *measurement_vec* or 0. yf.size(): {}, y.size(): {}",
      yf.size(),
      y.size());
  ARTS_USER_ERROR_IF(covmat_se.ncols() != covmat_se.nrows(),
                     "*measurement_vec_error_covmat* must be a square matrix. covmat_se shape: {}x{}",
                     covmat_se.nrows(),
                     covmat_se.ncols());
  ARTS_USER_ERROR_IF(
      static_cast<Size>(covmat_se.ncols()) != m,
      "Inconsistency in size between *measurement_vec* and *measurement_vec_error_covmat*. y.size(): {}, covmat_se size: {}x{}",
      y.size(),
      covmat_se.nrows(),
      covmat_se.ncols());
  ARTS_USER_ERROR_IF(
      (static_cast<Size>(jacobian.nrows()) != m) && (!jacobian.empty()),
      "The number of rows of the jacobian must be either the number of elements in *measurement_vec* or 0. jacobian.nrows(): {}, y.size(): {}",
      jacobian.nrows(),
      y.size());
  ARTS_USER_ERROR_IF(
      (static_cast<Size>(jacobian.ncols()) != n) && (!jacobian.empty()),
      "The number of cols of the jacobian must be either the number of elements in *model_state_vec_apriori* or 0. jacobian.ncols(): {}, xa.size(): {}",
      jacobian.ncols(),
      xa.size());

  ARTS_USER_ERROR_IF(n == 0 || m == 0, "OEM requires nonempty state and measurement vectors.")
  ARTS_USER_ERROR_IF(n > std::numeric_limits<unsigned int>::max() || m > std::numeric_limits<unsigned int>::max(),
                     "OEM state and measurement dimensions exceed the solver index range.")
  ARTS_USER_ERROR_IF(!normalization.empty() && normalization.size() != n,
                     "model_state_covmat_normalization must be empty or have {} elements.",
                     n)
  ARTS_USER_ERROR_IF(stdr::any_of(normalization, [](Numeric v) { return !std::isfinite(v) || v <= 0; }),
                     "model_state_covmat_normalization values must be finite and > 0.")
  ARTS_USER_ERROR_IF(
      method.measurement_space && !normalization.empty(),
      "model_state_covmat_normalization is not supported for measurement-space methods; use li_cg or gn_cg with normalization.")
  ARTS_USER_ERROR_IF(max_iter <= 0 || max_iter >= std::numeric_limits<unsigned int>::max(),
                     "max_iter must be positive and less than {}.",
                     std::numeric_limits<unsigned int>::max())
  ARTS_USER_ERROR_IF(!std::isfinite(stop_dx) || stop_dx <= 0, "stop_dx must be finite and > 0.")
  ARTS_USER_ERROR_IF(std::isnan(max_start_cost), "max_start_cost must not be NaN.")
  ARTS_USER_ERROR_IF(clear_matrices < 0 || clear_matrices > 1, "clear_matrices must be 0 or 1.")
  ARTS_USER_ERROR_IF(display_progress < 0 || display_progress > 1, "display_progress must be 0 or 1.")
}

// Shared iteration dispatch for full and reduced forward-model adapters.
template <typename Forward> void oem_compute(Forward&                          aw,
                                             oem::Vector&                      x_oem,
                                             OptimalEstimationDiagnostics&     oem_diagnostics,
                                             const Vector&                     model_state_vec_apriori,
                                             const Vector&                     measurement_vec,
                                             const CovarianceMatrix&           model_state_covmat,
                                             const CovarianceMatrix&           measurement_vec_error_covmat,
                                             const OEMMethod&                  selected,
                                             const Vector&                     model_state_covmat_normalization,
                                             const Vector&                     measurement_vec_normalization,
                                             Index                             max_iter,
                                             Numeric                           stop_dx,
                                             const LevenbergMarquardtSettings& lm_ga_settings,
                                             Index                             display_progress,
                                             const BlockMatrix*                projected_damping = nullptr) {
  const Index n = model_state_vec_apriori.size(), m = measurement_vec.size();
  auto&       lm_ga_history = oem_diagnostics.lm_ga_history;
  auto&       errors        = oem_diagnostics.errors;
  bool        apply_norm    = false;
  oem::Matrix T{};
  if (model_state_covmat_normalization.size() == static_cast<Size>(n)) {
    T.resize(n, n);
    static_cast<::Matrix&>(T) = 0.0;
    diagonal(T)               = model_state_covmat_normalization;
    apply_norm                = true;
  }

  oem::CovarianceMatrix Se(measurement_vec_error_covmat), Sa(model_state_covmat);
  oem::Vector           xa_oem(model_state_vec_apriori), y_oem(measurement_vec);
  const auto            iterations = static_cast<unsigned int>(selected.linear() ? 1 : max_iter);
  const auto            verbosity  = static_cast<unsigned int>(display_progress);

  // Read diagnostics from the formulation that actually ran, including when
  // the forward model throws. Costs always use the same measurement scaling.
  auto run = [&]<typename Retrieval, typename Optimizer>(Retrieval& retrieval, Optimizer& optimizer) {
    auto diagnostics = [&] {
      oem_diagnostics.final_cost       = retrieval.cost / static_cast<Numeric>(m);
      oem_diagnostics.measurement_cost = retrieval.cost_y / static_cast<Numeric>(m);
      oem_diagnostics.iterations       = static_cast<Index>(retrieval.iterations);
    };
    retrieval.iterations = 0;
    try {
      const auto status = retrieval.template compute<Optimizer&, oem::ArtsLog>(
          x_oem, y_oem, optimizer, verbosity, lm_ga_history, selected.linear());
      oem_diagnostics.status =
          status == 0 ? OptimalEstimationStatus::Converged : OptimalEstimationStatus::IterationLimit;
    } catch (...) {
      diagnostics();
      throw;
    }
    diagnostics();
  };

  auto solve = [&]<typename Solver>(Solver& solver) {
    if (selected.damped()) {
      // D = diag(Sa^-1), not diag(Sa)^-1 when the prior is correlated.
      CovarianceMatrix damping;
      if (projected_damping) {
        damping.add_correlation_inverse(Block(Range(0, n), Range(0, n), {0, 0}, *projected_damping));
      } else {
        damping.add_correlation_inverse(
            Block(Range(0, n),
                  Range(0, n),
                  std::make_pair(0, 0),
                  std::make_shared<Sparse>(Sparse::diagonal(model_state_covmat.inverse_diagonal()))));
      }
      oem::CovarianceMatrix precision = inv(oem::CovarianceMatrix(damping));
      invlib::LevenbergMarquardt<Numeric, oem::CovarianceMatrix, Solver> optimizer(precision, solver);
      configure_lm(optimizer, lm_ga_settings, stop_dx, iterations);
      oem::OEM_STANDARD<Forward> retrieval(aw, xa_oem, Sa, Se);
      run(retrieval, optimizer);
      if (optimizer.get_stop_reason() == invlib::LMStopReason::DampingLimit)
        oem_diagnostics.status = OptimalEstimationStatus::DampingLimit;
    } else {
      invlib::GaussNewton<Numeric, Solver> optimizer(stop_dx, iterations, solver);
      // Both measurement solvers accept the lazy system.
      if constexpr (std::is_same_v<Solver, oem::CG> or std::is_same_v<Solver, oem::DirectMeasurementSolver>) {
        if (selected.measurement_space) {
          oem::OEM_MFORM<Forward> retrieval(aw, xa_oem, Sa, Se);
          run(retrieval, optimizer);
          return;
        }
      }
      oem::OEM_STANDARD<Forward> retrieval(aw, xa_oem, Sa, Se);
      run(retrieval, optimizer);
    }
  };

  if (selected.conjugate_gradient) {
    oem::CG solver(T, apply_norm, 1e-10, 0);
    solver.measurement_scales = measurement_vec_normalization;
    bool warned               = false;
    solver.set_iteration_limit_warning([&] {
      if (not warned) {
        errors.emplace_back("Warning: CG iteration limit reached; OEM continued with the last linear-solver iterate.");
        warned = true;
      }
    });
    solve(solver);
  } else if (selected.measurement_space) {
    oem::DirectMeasurementSolver solver{measurement_vec_normalization};
    solve(solver);
  } else {
    oem::Std solver(T, apply_norm);
    solve(solver);
  }
}
}  // namespace

void model_state_vec_aprioriFromState(Vector& xa, const Vector& x) {
  ARTS_TIME_REPORT

  xa = x;
}

void measurement_vec_fitFromMeasurement(Vector& yf, const Vector& y) {
  ARTS_TIME_REPORT

  yf = y;
}

void measurement_vec_error_covmatNormalization(Vector& normalization, const CovarianceMatrix& covariance) {
  covariance.validate(covariance.nrows());
  Vector scales = covariance.diagonal();
  for (auto& value : scales) {
    ARTS_USER_ERROR_IF(not std::isfinite(value) or value <= 0,
                       "Measurement noise variances must be finite and strictly positive.")
    value = std::sqrt(value);
  }
  normalization = std::move(scales);
}

void OEM(const Workspace&                  ws,
         Vector&                           model_state_vec,
         Vector&                           measurement_vec_fit,
         Matrix&                           measurement_jac,
         AtmField&                         atm_field,
         AbsorptionBands&                  abs_bands,
         ArrayOfSensorObsel&               measurement_sensor,
         SurfaceField&                     surf_field,
         SubsurfaceField&                  subsurf_field,
         Matrix&                           measurement_gain_mat,
         OptimalEstimationDiagnostics&     oem_diagnostics,
         const JacobianTargets&            jac_targets,
         const Vector&                     model_state_vec_apriori,
         const CovarianceMatrix&           model_state_covmat_input,
         const Vector&                     measurement_vec,
         const CovarianceMatrix&           measurement_vec_error_covmat_input,
         const Agenda&                     inversion_iterate_agenda,
         const String&                     method,
         const Numeric&                    max_start_cost,
         const Vector&                     model_state_covmat_normalization,
         const Vector&                     measurement_vec_normalization,
         const Index&                      max_iter,
         const Numeric&                    stop_dx,
         const LevenbergMarquardtSettings& lm_ga_settings,
         const Index&                      clear_matrices,
         const Index&                      display_progress) {
  ARTS_TIME_REPORT

  const OEMMethod selected = parse_oem_method(method);
  // Freeze covariance values and finish cache preparation before iteration.
  const auto  state_snapshot = model_state_covmat_input.prepared(not selected.measurement_space or clear_matrices == 0);
  const auto  measurement_snapshot         = measurement_vec_error_covmat_input.prepared();
  const auto& model_state_covmat           = *state_snapshot;
  const auto& measurement_vec_error_covmat = *measurement_snapshot;

  ARTS_USER_ERROR_IF(not measurement_vec_normalization.empty() and not selected.measurement_space,
                     "measurement_vec_normalization is only supported for measurement-space methods.")
  ARTS_USER_ERROR_IF(
      not measurement_vec_normalization.empty() and measurement_vec_normalization.size() != measurement_vec.size(),
      "measurement_vec_normalization must be empty or have {} elements.",
      measurement_vec.size())
  for (const auto scale : measurement_vec_normalization)
    ARTS_USER_ERROR_IF(not std::isfinite(scale) or scale <= 0,
                       "measurement_vec_normalization values must be finite and > 0.")

  check_oem_inputs(model_state_vec,
                   measurement_vec_fit,
                   measurement_jac,
                   model_state_vec_apriori,
                   model_state_covmat,
                   measurement_vec,
                   measurement_vec_error_covmat,
                   selected,
                   model_state_covmat_normalization,
                   max_iter,
                   stop_dx,
                   max_start_cost,
                   clear_matrices,
                   display_progress);
  if (selected.damped()) lm_ga_settings.validate();

  const Index n = model_state_covmat.nrows();
  const Index m = measurement_vec.size();
  // Covariance snapshots are prepared before iteration. Explicit precision
  // consumers still request inverse storage when needed.

  measurement_gain_mat.resize(0, 0);
  oem_diagnostics     = {};
  auto& lm_ga_history = oem_diagnostics.lm_ga_history;
  auto& errors        = oem_diagnostics.errors;
  lm_ga_history.resize(selected.damped() ? max_iter + 1 : 0);
  lm_ga_history = NAN;

  // A cached simulation/Jacobian pair must correspond to the supplied start.
  if (model_state_vec.empty()) {
    model_state_vec = model_state_vec_apriori;
    measurement_vec_fit.resize(0);
    measurement_jac.resize(0, 0);
  }
  if (measurement_vec_fit.empty() || measurement_jac.empty()) {
    inversion_iterate_agendaExecute(ws,
                                    atm_field,
                                    abs_bands,
                                    measurement_sensor,
                                    surf_field,
                                    subsurf_field,
                                    measurement_vec_fit,
                                    measurement_jac,
                                    jac_targets,
                                    jac_targets,
                                    model_state_vec,
                                    inversion_iterate_agenda);
  }
  ARTS_USER_ERROR_IF(measurement_jac.nrows() != m || measurement_jac.ncols() != n,
                     "inversion_iterate_agenda must return an {} by {} measurement_jac; got {} by {}.",
                     m,
                     n,
                     measurement_jac.nrows(),
                     measurement_jac.ncols())

  ARTS_USER_ERROR_IF(
      measurement_vec_fit.size() not_eq measurement_vec.size(),
      "Mismatch between simulated y and input y.\n"
      "Input y is size {}"
      " but simulated y is {}"
      "\n"
      "Use your frequency grid vector and your sensor response matrix to match simulations with measurements.\n",
      measurement_vec.size(),
      measurement_vec_fit.size())

  // TODO: Get this from invlib log.
  // Start value of cost function
  Numeric cost_start = NAN;
  if (selected.damped() || display_progress || max_start_cost > 0) {
    Vector dy   = measurement_vec;
    dy         -= measurement_vec_fit;
    Vector sdy  = measurement_vec;
    mult_inv(sdy.view_as(sdy.size(), 1), measurement_vec_error_covmat, dy.view_as(dy.size(), 1));
    Vector dx   = model_state_vec;
    dx         -= model_state_vec_apriori;
    Vector sdx  = model_state_vec;
    mult_inv(sdx.view_as(sdx.size(), 1), model_state_covmat, dx.view_as(dx.size(), 1));
    cost_start  = dot(dx, sdx) + dot(dy, sdy);
    cost_start /= static_cast<Numeric>(m);
  }
  oem_diagnostics.initial_cost = cost_start;

  // Handle cases with too large start cost
  if (max_start_cost > 0 && cost_start > max_start_cost) {
    // No inversion; retain the starting state and its simulated measurement.
    oem_diagnostics.status = OptimalEstimationStatus::StartCostLimit;
    if (clear_matrices) measurement_jac.resize(0, 0);
    //
    if (display_progress) {
      std::cout << "\n   No OEM inversion, too high start cost:\n"
                << "        Set limit : " << max_start_cost << '\n'
                << "      Found value : " << cost_start << '\n'
                << '\n';
    }
  }
  // Otherwise do inversion
  else {
    oem::Vector        x_oem(model_state_vec);
    oem::AgendaWrapper aw(&ws,
                          static_cast<unsigned int>(m),
                          static_cast<unsigned int>(n),
                          measurement_jac,
                          measurement_vec_fit,
                          model_state_vec,
                          &atm_field,
                          &abs_bands,
                          &measurement_sensor,
                          &surf_field,
                          &subsurf_field,
                          &jac_targets,
                          &inversion_iterate_agenda);
    try {
      oem_compute(aw,
                  x_oem,
                  oem_diagnostics,
                  model_state_vec_apriori,
                  measurement_vec,
                  model_state_covmat,
                  measurement_vec_error_covmat,
                  selected,
                  model_state_covmat_normalization,
                  measurement_vec_normalization,
                  max_iter,
                  stop_dx,
                  lm_ga_settings,
                  display_progress);
      // Ensure that the returned gain/Jacobian describe the retrieved state.
      // An already current Jacobian needs neither an agenda call nor a fit copy.
      if (!selected.linear() && !clear_matrices) aw.ensure_jacobian(x_oem);
    } catch (const std::exception& e) {
      oem_diagnostics.status        = OptimalEstimationStatus::Error;
      static_cast<::Vector&>(x_oem) = NAN;
      for (const auto& message : oem::handle_nested_exception(e)) {
        std::stringstream stream{message};
        for (std::string line; std::getline(stream, line);) errors.push_back(line);
      }
    }

    model_state_vec = x_oem;

    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices) {
      measurement_jac.resize(0, 0);
      measurement_gain_mat.resize(0, 0);
    } else if ((oem_diagnostics.status == OptimalEstimationStatus::Converged or
                oem_diagnostics.status == OptimalEstimationStatus::IterationLimit or
                oem_diagnostics.status == OptimalEstimationStatus::DampingLimit)) {
      measurement_gain_mat.resize(n, m);
      Matrix tmp1(n, m), tmp2(n, n), tmp3(n, n);
      mult_inv(tmp1, transpose(measurement_jac), measurement_vec_error_covmat);
      mult(tmp2, tmp1, measurement_jac);
      add_inv(tmp2, model_state_covmat);
      inv(tmp3, tmp2);
      mult(measurement_gain_mat, tmp3, tmp1);
    }
  }
}

namespace {
// Apply diagonal precision without changing basis storage.
BlockMatrix scale_basis_rows(const BlockMatrix& basis, const Vector& scale) {
  return std::visit(
      [&]<typename T>(const std::shared_ptr<T>& matrix) -> BlockMatrix {
        auto out = std::make_shared<T>(*matrix);
        if constexpr (std::same_as<T, Sparse>) {
          for (auto [row, col, value] : *out | by_elem) value *= scale[row];
        } else {
          for (Index row = 0; row < out->nrows(); ++row)
            for (Index col = 0; col < out->ncols(); ++col) (*out)[row, col] *= scale[row];
        }
        return out;
      },
      basis.data);
}

BlockMatrix basis_precision(const BlockMatrix& basis, const BlockMatrix& weighted) {
  const Index r = basis.ncols();
  if (weighted.is_sparse()) {
    assert(basis.is_sparse());
    auto out    = std::make_shared<Sparse>(r, r);
    out->matrix = basis.sparse().matrix.transpose() * weighted.sparse().matrix;
    return out;
  }
  auto out = std::make_shared<Matrix>(r, r);
  basis.multiply_right(transpose(*out), transpose(weighted.dense()));
  return out;
}

CovarianceMatrix basis_covariance(BlockMatrix values) {
  const Index      size = values.nrows();
  CovarianceMatrix out;
  out.add_correlation(Block(Range(0, size), Range(0, size), {0, 0}, std::move(values)));
  return out;
}

// This is preparation, never an operation in an OEM iteration. Keep sparse
// projections and sparse/diagonal noise sparse all the way to the factor cache.
CovarianceMatrix measurement_covariance_projection(const BlockMatrix& C, const CovarianceMatrix& noise) {
  const Index m = C.ncols(), q = C.nrows();
  BlockMatrix values;
  if (C.is_sparse()) {
    Sparse weighted;
    if (const auto diagonal = noise.diagonal_if_diagonal()) {
      weighted = C.sparse();
      for (auto [row, col, value] : weighted | by_elem) value *= (*diagonal)[col];
    } else if (stdr::all_of(noise.get_blocks(), [](const auto& block) { return block.is_sparse(); })) {
      std::vector<Eigen::Triplet<Numeric>> entries;
      for (const auto& block : noise.get_blocks()) {
        const auto r0 = block.get_row_range().offset, c0 = block.get_column_range().offset;
        const auto [i, j] = block.get_indices();
        for (const auto [row, col, value] : block.get_sparse() | by_elem) {
          entries.emplace_back(r0 + row, c0 + col, value);
          if (i != j) entries.emplace_back(c0 + col, r0 + row, value);
        }
      }
      Sparse matrix(m, m);
      matrix.matrix.setFromTriplets(entries.begin(), entries.end());
      weighted.resize(q, m);
      mult(weighted, C.sparse(), matrix);
    }
    if (weighted.nrows() == q) {
      auto sparse    = std::make_shared<Sparse>(q, q);
      sparse->matrix = weighted.matrix * C.sparse().matrix.transpose();
      values         = std::move(sparse);
    }
  }
  if (not values.not_null()) {
    Matrix weighted(m, q);
    if (C.is_dense())
      mult(weighted, noise, transpose(C.dense()));
    else {
      Matrix ct(m, q, 0.);
      for (const auto [row, col, value] : C.sparse() | by_elem) ct[col, row] = value;
      mult(weighted, noise, ct);
    }
    auto dense = std::make_shared<Matrix>(q, q);
    C.multiply_left(*dense, weighted);
    values = std::move(dense);
  }
  CovarianceMatrix result;
  result.add_correlation(Block(Range(0, q), Range(0, q), {0, 0}, std::move(values)));
  return result;
}
}  // namespace

void measurement_basis_matCalc(BlockMatrix&            measurement_basis_mat,
                               const Matrix&           measurement_jac,
                               const CovarianceMatrix& measurement_vec_error_covmat) {
  ARTS_TIME_REPORT
  const auto& J = measurement_jac;
  const Index m = J.nrows(), n = J.ncols();
  ARTS_USER_ERROR_IF(m <= 0 or n <= 0 or measurement_vec_error_covmat.nrows() != m,
                     "A nonempty measurement_jac and matching measurement covariance are required.")
  ARTS_USER_ERROR_IF(stdr::any_of(J | by_elem, [](Numeric v) { return not std::isfinite(v); }),
                     "measurement_jac must be finite.")
  const auto                            noise = measurement_vec_error_covmat.prepared();
  std::map<std::vector<Numeric>, Index> directions;
  ArrayOfIndex                          group(m);
  Vector                                amplitude(m);
  // Canonicalize the complete row, retaining the sign in its amplitude.
  // Exact keys deliberately do not merge merely similar sensitivities.
  for (Index i = 0; i < m; ++i) {
    Index pivot = 0;
    for (Index k = 1; k < n; ++k)
      if (std::abs(J[i, k]) > std::abs(J[i, pivot])) pivot = k;
    amplitude[i] = J[i, pivot] == 0 ? 1 : J[i, pivot];
    std::vector<Numeric> row(n);
    for (Index k = 0; k < n; ++k) {
      row[k] = J[i, k] / amplitude[i];
      ARTS_USER_ERROR_IF(
          (row[k] == 0 and J[i, k] != 0), "Jacobian row {} spans too many orders of magnitude to group reliably.", i)
    }
    const auto [it, inserted] = directions.try_emplace(std::move(row), directions.size());
    group[i]                  = it->second;
  }
  const Index q = directions.size();
  // No grouping: an identity projection avoids any noise solve or dense basis.
  if (q == m) {
    auto identity = std::make_shared<Sparse>(m, m);
    id_mat(*identity);
    measurement_basis_mat = std::move(identity);
    return;
  }
  Vector scale(q, 0.);
  for (Index i = 0; i < m; ++i) scale[group[i]] = std::max(scale[group[i]], std::abs(amplitude[i]));
  for (Index i = 0; i < m; ++i) {
    amplitude[i] /= scale[group[i]];
    ARTS_USER_ERROR_IF(amplitude[i] == 0, "Channel scaling underflow during measurement grouping.")
  }
  if (const auto diagonal = noise->diagonal_if_diagonal()) {
    Vector norm(q, 0.);
    for (Index i = 0; i < m; ++i) norm[group[i]] = std::hypot(norm[group[i]], amplitude[i] / std::sqrt((*diagonal)[i]));
    std::vector<Eigen::Triplet<Numeric>> entries;
    entries.reserve(m);
    for (Index i = 0; i < m; ++i) {
      const Numeric sigma = std::sqrt((*diagonal)[i]);
      const Numeric value = (amplitude[i] / sigma / norm[group[i]]) / sigma;
      ARTS_USER_ERROR_IF(not std::isfinite(value), "Measurement basis exceeds numerical range.")
      entries.emplace_back(group[i], i, value);
    }
    auto sparse = std::make_shared<Sparse>(q, m);
    sparse->matrix.setFromTriplets(entries.begin(), entries.end());
    measurement_basis_mat = std::move(sparse);
  } else {
    // J = T R. C = T^T Se^-1 retains all state-dependent likelihood terms,
    // including information carried by noise correlations outside each group.
    Matrix T(m, q, 0.), weighted(m, q);
    for (Index i = 0; i < m; ++i) T[i, group[i]] = amplitude[i];
    mult_inv(weighted, *noise, T);
    ARTS_USER_ERROR_IF(stdr::any_of(weighted | by_elem, [](Numeric v) { return not std::isfinite(v); }),
                       "Measurement basis exceeds numerical range.")
    measurement_basis_mat = std::make_shared<Matrix>(transpose(weighted));
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReducedOEMBasisCalc(BlockMatrix&            model_state_basis_mat,
                         BlockMatrix&            measurement_basis_mat,
                         Vector&                 oem_basis_singular_values,
                         const Matrix&           measurement_jac,
                         const CovarianceMatrix& model_state_covmat,
                         const CovarianceMatrix& measurement_vec_error_covmat) {
  ARTS_TIME_REPORT
  const Index m = measurement_jac.nrows(), n = measurement_jac.ncols();
  ARTS_USER_ERROR_IF(m <= 0 or n <= 0, "ReducedOEMBasisCalc requires a nonempty measurement_jac.")
  ARTS_USER_ERROR_IF(model_state_covmat.nrows() != n or measurement_vec_error_covmat.nrows() != m,
                     "Covariance sizes must match the {} measurement_jac columns and {} rows.",
                     n,
                     m)
  ARTS_USER_ERROR_IF(stdr::any_of(measurement_jac | by_elem, [](auto v) { return not std::isfinite(v); }),
                     "measurement_jac must be finite.")

  const CovarianceSquareRoot prior(model_state_covmat);
  const CovarianceSquareRoot noise(measurement_vec_error_covmat);
  Matrix                     scaled(m, n), whitened(m, n);
  prior.multiply_left(transpose(scaled), transpose(measurement_jac), true);
  noise.solve_left(whitened, scaled);
  ARTS_USER_ERROR_IF(stdr::any_of(whitened | by_elem, [](auto v) { return not std::isfinite(v); }),
                     "Whitened Jacobian is not finite; check covariance and Jacobian scales.")

  Matrix u, v;
  Vector singular_values;
  // Preserve both null spaces; choosing which modes to discard is a separate step.
  svd(u, singular_values, v, whitened);
  ARTS_USER_ERROR_IF(stdr::any_of(singular_values, [](auto x) { return not std::isfinite(x); }),
                     "Information spectrum exceeds numerical range.")

  Matrix B(n, n), C(m, m);
  prior.multiply_left(B, v);
  noise.solve_left(transpose(C), u, true);
  ARTS_USER_ERROR_IF(stdr::any_of(B | by_elem, [](auto x) { return not std::isfinite(x); }) or
                         stdr::any_of(C | by_elem, [](auto x) { return not std::isfinite(x); }),
                     "Basis matrices exceed numerical range.")

  // Publish the matched bases and spectrum only after all transformations succeed.
  model_state_basis_mat     = std::make_shared<Matrix>(std::move(B));
  measurement_basis_mat     = std::make_shared<Matrix>(std::move(C));
  oem_basis_singular_values = std::move(singular_values);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void ReducedOEMBasisReduce(BlockMatrix&   model_state_basis_mat,
                           BlockMatrix&   measurement_basis_mat,
                           Numeric&       oem_basis_lost_dofs,
                           Numeric&       oem_basis_lost_information_bits,
                           const Vector&  oem_basis_singular_values,
                           const Index&   rank,
                           const Numeric& max_lost_dofs,
                           const Numeric& max_lost_information_bits) {
  ARTS_TIME_REPORT
  const auto& B               = model_state_basis_mat;
  const auto& C               = measurement_basis_mat;
  const auto& singular_values = oem_basis_singular_values;
  const Index n = B.nrows(), m = C.ncols(), p = std::min(m, n);
  ARTS_USER_ERROR_IF(
      n <= 0 or m <= 0 or B.ncols() <= 0 or B.ncols() > n or C.nrows() <= 0 or C.nrows() > m,
      "Basis matrices must be nonempty with at most the original number of modes; use ReducedOEMBasisCalc first.")
  ARTS_USER_ERROR_IF(singular_values.size() != static_cast<std::size_t>(p) or
                         stdr::any_of(singular_values, [](auto x) { return not std::isfinite(x) or x < 0; }) or
                         not std::is_sorted(singular_values.begin(), singular_values.end(), std::greater<>{}),
                     "oem_basis_singular_values must contain {} finite, nonnegative values in descending order "
                     "from the same ReducedOEMBasisCalc call as the full bases.",
                     p)
  ARTS_USER_ERROR_IF(rank != -1 and (rank < 1 or rank > n),
                     "rank must be -1 for automatic selection or between 1 and {} full state variables.",
                     n)

  const auto valid_loss = [](Numeric value) { return value == -1 or (std::isfinite(value) and value >= 0); };
  ARTS_USER_ERROR_IF(not valid_loss(max_lost_dofs) or not valid_loss(max_lost_information_bits),
                     "Information-loss limits must be finite and nonnegative, or -1 to leave a limit unset.")
  ARTS_USER_ERROR_IF(rank != -1 and (max_lost_dofs != -1 or max_lost_information_bits != -1),
                     "Supply either an explicit rank or information-loss limits, not both.")

  const auto contribution = [](Numeric s) {
    const Numeric fraction = s / std::hypot(1, s);
    const Numeric bits =
        (s <= 1 ? .5 * std::log1p(s * s) : std::log(s) + .5 * std::log1p((1 / s) * (1 / s))) / std::log(2.);
    return std::pair{fraction * fraction, bits};
  };
  Index   retained  = rank;
  Numeric lost_dofs = 0, lost_bits = 0;
  if (retained == -1) {
    const Numeric dofs_limit = max_lost_dofs == -1 ? std::numeric_limits<Numeric>::infinity() : max_lost_dofs;
    const Numeric bits_limit = max_lost_information_bits == -1
                                   ? (max_lost_dofs == -1 ? 0 : std::numeric_limits<Numeric>::infinity())
                                   : max_lost_information_bits;
    retained                 = p;
    // Sum from the weakest end: subtracting from a large total loses weak tails.
    // Keep one coefficient even when no mode is informative, as required by ReducedOEM.
    while (retained > 1) {
      const auto [dofs, bits] = contribution(singular_values[retained - 1]);
      const Numeric next_dofs = lost_dofs + dofs;
      const Numeric next_bits = lost_bits + bits;
      if (next_dofs > dofs_limit or next_bits > bits_limit) break;
      lost_dofs = next_dofs;
      lost_bits = next_bits;
      --retained;
    }
  } else {
    for (Index i = p; i > retained; --i) {
      const auto [dofs, bits]  = contribution(singular_values[i - 1]);
      lost_dofs               += dofs;
      lost_bits               += bits;
    }
  }

  const Index q = std::min(retained, m);
  ARTS_USER_ERROR_IF(
      retained > B.ncols() or q > C.nrows(),
      "Requested modes have already been removed. Restore saved bases or rerun ReducedOEMBasisCalc before increasing rank or tightening loss limits.")
  const auto leading = [](const BlockMatrix& basis, Index rows, Index cols) {
    return std::visit(
        [&]<typename T>(const std::shared_ptr<T>& matrix) -> BlockMatrix {
          if constexpr (std::same_as<T, Matrix>)
            return std::make_shared<Matrix>((*matrix)[Range(0, rows), Range(0, cols)]);
          else {
            auto out    = std::make_shared<Sparse>(rows, cols);
            out->matrix = matrix->matrix.topLeftCorner(rows, cols);
            return out;
          }
        },
        basis.data);
  };
  BlockMatrix reduced_B = leading(B, n, retained);
  BlockMatrix reduced_C = leading(C, q, m);
  // Materialize both slices before replacing either input; retain the full spectrum for total losses.
  model_state_basis_mat           = std::move(reduced_B);
  measurement_basis_mat           = std::move(reduced_C);
  oem_basis_lost_dofs             = lost_dofs;
  oem_basis_lost_information_bits = lost_bits;
}

void ReducedOEM(const Workspace&                  ws,
                Vector&                           model_state_vec,
                Vector&                           measurement_vec_fit,
                Matrix&                           measurement_jac,
                AtmField&                         atm_field,
                AbsorptionBands&                  abs_bands,
                ArrayOfSensorObsel&               measurement_sensor,
                SurfaceField&                     surf_field,
                SubsurfaceField&                  subsurf_field,
                Matrix&                           measurement_gain_mat,
                OptimalEstimationDiagnostics&     oem_diagnostics,
                const JacobianTargets&            jac_targets,
                const Vector&                     model_state_vec_apriori,
                const CovarianceMatrix&           model_state_covmat_input,
                const Vector&                     measurement_vec,
                const CovarianceMatrix&           measurement_vec_error_covmat_input,
                const Agenda&                     inversion_iterate_agenda,
                const BlockMatrix&                model_state_basis_mat,
                const BlockMatrix&                measurement_basis_mat,
                const String&                     method,
                const Numeric&                    max_start_cost,
                const Vector&                     model_state_covmat_normalization,
                const Vector&                     measurement_vec_normalization,
                const Index&                      max_iter,
                const Numeric&                    stop_dx,
                const LevenbergMarquardtSettings& lm_ga_settings,
                const Index&                      clear_matrices,
                const Index&                      display_progress) {
  ARTS_TIME_REPORT
  const auto  selected = parse_oem_method(method);
  const auto& B        = model_state_basis_mat;
  const auto& C        = measurement_basis_mat;
  const Index n = model_state_vec_apriori.size(), m = measurement_vec.size();
  const Index r = B.ncols(), q = C.nrows();
  ARTS_USER_ERROR_IF(
      B.nrows() != n or r <= 0 or r > n, "model_state_basis_mat must have {} rows and between 1 and {} columns.", n, n)
  ARTS_USER_ERROR_IF(
      C.ncols() != m or q <= 0 or q > m, "measurement_basis_mat must have {} columns and between 1 and {} rows.", m, m)
  ARTS_USER_ERROR_IF(not B.is_finite() or not C.is_finite(), "ReducedOEM reduction matrices must be finite.")
  ARTS_USER_ERROR_IF(stdr::any_of(model_state_vec_apriori, [](auto v) { return not std::isfinite(v); }),
                     "ReducedOEM prior state must be finite.")

  const bool state_identity = B.is_identity();
  const auto prior          = model_state_covmat_input.prepared(
      state_identity ? (not selected.measurement_space or clear_matrices == 0) : selected.damped());
  const auto noise = measurement_vec_error_covmat_input.prepared();
  check_oem_inputs(model_state_vec,
                   measurement_vec_fit,
                   measurement_jac,
                   model_state_vec_apriori,
                   *prior,
                   measurement_vec,
                   *noise,
                   selected,
                   {},
                   max_iter,
                   stop_dx,
                   max_start_cost,
                   clear_matrices,
                   display_progress);
  if (selected.damped()) lm_ga_settings.validate();

  std::shared_ptr<const CovarianceMatrix> reduced_prior, reduced_noise;
  Vector                                  za(r, 0), start(r, 0), reduced_y(q);
  BlockMatrix                             damping;
  {
    // Factor and validate once; release preparation scratch before iteration.
    BlockMatrix weighted;
    if (state_identity) {
      reduced_prior = prior;
    } else {
      const auto diagonal = B.is_sparse() ? prior->diagonal_if_diagonal() : std::nullopt;
      if (diagonal) {
        Vector inverse = *diagonal;
        for (auto& value : inverse) value = 1 / value;
        weighted = scale_basis_rows(B, inverse);
      } else {
        auto values = std::make_shared<Matrix>(n, r);
        std::visit(
            [&]<typename T>(const std::shared_ptr<T>& matrix) {
              if constexpr (std::same_as<T, Matrix>)
                mult_inv(*values, *prior, *matrix);
              else
                mult_inv(*values, *prior, Matrix(*matrix));
            },
            B.data);
        weighted = std::move(values);
      }
      const auto  reduced_precision = basis_covariance(basis_precision(B, weighted)).prepared();
      BlockMatrix covariance;
      if (auto diagonal = reduced_precision->diagonal_if_diagonal()) {
        for (auto& value : *diagonal) value = 1 / value;
        covariance = std::make_shared<Sparse>(Sparse::diagonal(*diagonal));
      } else {
        Matrix identity(r, r);
        id_mat(identity);
        auto values = std::make_shared<Matrix>(r, r);
        mult_inv(*values, *reduced_precision, identity);
        covariance = std::move(values);
      }
      reduced_prior =
          basis_covariance(std::move(covariance)).prepared(not selected.measurement_space or clear_matrices == 0);
    }

    reduced_noise = measurement_covariance_projection(C, *noise).prepared();

    C.multiply_left(reduced_y, measurement_vec);
    if (not model_state_vec.empty()) {
      Vector delta  = model_state_vec;
      delta        -= model_state_vec_apriori;
      Vector rhs(r), represented(n);
      if (state_identity) {
        start       = delta;
        represented = delta;
      } else {
        weighted.multiply_right(rhs.view_as(1, r), delta.view_as(1, n));
        mult(start, *reduced_prior, rhs);
        B.multiply_left(represented, start);
      }
      for (Index i = 0; i < n; ++i)
        ARTS_USER_ERROR_IF(not std::isfinite(delta[i]) or not std::isfinite(represented[i]) or
                               std::abs(represented[i] - delta[i]) > 1e-8 * (1 + std::abs(delta[i])),
                           "ReducedOEM starting state must lie in the supplied affine subspace.")
    }
    check_oem_inputs(start,
                     {},
                     {},
                     za,
                     *reduced_prior,
                     reduced_y,
                     *reduced_noise,
                     selected,
                     model_state_covmat_normalization,
                     max_iter,
                     stop_dx,
                     max_start_cost,
                     clear_matrices,
                     display_progress);
    ARTS_USER_ERROR_IF(
        not measurement_vec_normalization.empty() and
            (not selected.measurement_space or measurement_vec_normalization.size() != static_cast<Size>(q)),
        "ReducedOEM measurement_vec_normalization requires a measurement-space method and {} elements.",
        q)
    ARTS_USER_ERROR_IF(
        stdr::any_of(measurement_vec_normalization, [](auto v) { return not std::isfinite(v) or v <= 0; }),
        "measurement_vec_normalization values must be finite and > 0.")

    if (selected.damped() and not state_identity)
      damping = basis_precision(B, scale_basis_rows(B, prior->inverse_diagonal()));
  }

  measurement_gain_mat.resize(0, 0);
  oem_diagnostics = {};
  oem_diagnostics.lm_ga_history.resize(selected.damped() ? max_iter + 1 : 0);
  oem_diagnostics.lm_ga_history = NAN;
  if (model_state_vec.empty()) {
    model_state_vec = model_state_vec_apriori;
    measurement_vec_fit.resize(0);
    measurement_jac.resize(0, 0);
  }
  oem::AgendaWrapper        full(&ws,
                                 m,
                                 n,
                                 measurement_jac,
                                 measurement_vec_fit,
                                 model_state_vec,
                                 &atm_field,
                                 &abs_bands,
                                 &measurement_sensor,
                                 &surf_field,
                                 &subsurf_field,
                                 &jac_targets,
                                 &inversion_iterate_agenda);
  oem::ReducedAgendaWrapper reduced(full, model_state_vec_apriori, B, C, measurement_jac, state_identity);
  oem::Vector               z(start);
  model_state_vec = reduced.expand(z);
  reduced.ensure_jacobian(z);

  // Keep diagnostics in the original measurement space, including discarded residuals.
  const auto full_cost = [&] {
    Vector dx = model_state_vec, dy = measurement_vec;
    dx -= model_state_vec_apriori;
    dy -= measurement_vec_fit;
    Vector sx(n), sy(m);
    mult_inv(sx.view_as(n, 1), *prior, dx.view_as(n, 1));
    mult_inv(sy.view_as(m, 1), *noise, dy.view_as(m, 1));
    const Numeric measurement_cost = dot(dy, sy) / static_cast<Numeric>(m);
    return std::pair{dot(dx, sx) / static_cast<Numeric>(m) + measurement_cost, measurement_cost};
  };
  oem_diagnostics.initial_cost = full_cost().first;
  if (max_start_cost > 0 and oem_diagnostics.initial_cost > max_start_cost) {
    oem_diagnostics.status = OptimalEstimationStatus::StartCostLimit;
    if (clear_matrices) measurement_jac.resize(0, 0);
    return;
  }

  try {
    oem_compute(reduced,
                z,
                oem_diagnostics,
                za,
                reduced_y,
                *reduced_prior,
                *reduced_noise,
                selected,
                model_state_covmat_normalization,
                measurement_vec_normalization,
                max_iter,
                stop_dx,
                lm_ga_settings,
                display_progress,
                damping.not_null() ? &damping : nullptr);
    model_state_vec = reduced.expand(z);
    // LI and rejected LM trials also need the physical outputs at the returned state.
    if (clear_matrices) {
      full.ensure_measurement(reduced.expand(z));
    } else {
      reduced.ensure_jacobian(z);
      const auto& jac = reduced.get_jacobian();
      Matrix      rhs(r, q), hessian(r, r), posterior(r, r), gain(r, q);
      mult_inv(rhs, transpose(jac), *reduced_noise);
      mult(hessian, rhs, jac);
      add_inv(hessian, *reduced_prior);
      inv(posterior, hessian);
      mult(gain, posterior, rhs);
      measurement_gain_mat.resize(n, m);
      if (state_identity)
        C.multiply_right(measurement_gain_mat, gain);
      else {
        Matrix expanded_gain(n, q);
        B.multiply_left(expanded_gain, gain);
        C.multiply_right(measurement_gain_mat, expanded_gain);
      }
    }
    const auto [cost, measurement_cost] = full_cost();
    oem_diagnostics.final_cost          = cost;
    oem_diagnostics.measurement_cost    = measurement_cost;
  } catch (const std::exception& error) {
    oem_diagnostics.status           = OptimalEstimationStatus::Error;
    oem_diagnostics.final_cost       = NAN;
    oem_diagnostics.measurement_cost = NAN;
    model_state_vec                  = NAN;
    measurement_jac.resize(0, 0);
    measurement_gain_mat.resize(0, 0);
    for (const auto& message : oem::handle_nested_exception(error)) {
      std::stringstream stream{message};
      for (std::string line; std::getline(stream, line);) oem_diagnostics.errors.push_back(line);
    }
  }
  if (clear_matrices) measurement_jac.resize(0, 0);
}

void measurement_vec_error_covmat_observation_systemCalc(Matrix&       measurement_vec_error_covmat_observation_system,
                                                         const Matrix& measurement_gain_mat,
                                                         const CovarianceMatrix& measurement_vec_error_covmat) {
  ARTS_TIME_REPORT

  Index  n(measurement_gain_mat.nrows()), m(measurement_gain_mat.ncols());
  Matrix tmp1(m, n);

  ARTS_USER_ERROR_IF(
      (m == 0) || (n == 0),
      "The gain matrix *measurement_gain_mat* is required to compute the observation error covariance matrix.");

  measurement_vec_error_covmat_observation_system.resize(n, n);
  mult(tmp1, measurement_vec_error_covmat, transpose(measurement_gain_mat));
  mult(measurement_vec_error_covmat_observation_system, measurement_gain_mat, tmp1);
}

void model_state_covmat_smoothing_errorCalc(Matrix&                 model_state_covmat_smoothing_error,
                                            const Matrix&           measurement_averaging_kernel,
                                            const CovarianceMatrix& model_state_covmat) {
  ARTS_TIME_REPORT

  Index  n(measurement_averaging_kernel.ncols());
  Matrix tmp1(n, n), tmp2(n, n);

  ARTS_USER_ERROR_IF(
      n == 0,
      "The averaging kernel matrix *measurement_gain_mat* is required to compute the smoothing error covariance matrix.");

  model_state_covmat_smoothing_error.resize(n, n);

  // Sign doesn't matter since we're dealing with a quadratic form.
  id_mat(tmp1);
  tmp1 -= measurement_averaging_kernel;

  mult(tmp2, model_state_covmat, transpose(tmp1));
  mult(model_state_covmat_smoothing_error, tmp1, tmp2);
}

void measurement_averaging_kernelCalc(Matrix&       measurement_averaging_kernel,
                                      const Matrix& measurement_gain_mat,
                                      const Matrix& measurement_jac) {
  ARTS_TIME_REPORT

  Index n(measurement_jac.ncols());

  ARTS_USER_ERROR_IF(measurement_jac.empty(), "The Jacobian matrix is empty.");

  measurement_averaging_kernel.resize(n, n);
  mult(measurement_averaging_kernel, measurement_gain_mat, measurement_jac);
}
