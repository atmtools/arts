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
#include <limits>
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
  if (method == "li") return {.algorithm=OEMAlgorithm::Linear};
  if (method == "li_cg") return {.algorithm=OEMAlgorithm::Linear, .conjugate_gradient=true};
  if (method == "li_cg_m") return {.algorithm=OEMAlgorithm::Linear, .conjugate_gradient=true, .measurement_space=true};
  if (method == "gn") return {.algorithm=OEMAlgorithm::GaussNewton};
  if (method == "gn_cg") return {.algorithm=OEMAlgorithm::GaussNewton, .conjugate_gradient=true};
  if (method == "gn_cg_m") return {.algorithm=OEMAlgorithm::GaussNewton, .conjugate_gradient=true, .measurement_space=true};
  if (method == "lm" || method == "ml") return {.algorithm=OEMAlgorithm::LevenbergMarquardt};
  if (method == "lm_cg" || method == "ml_cg") return {.algorithm=OEMAlgorithm::LevenbergMarquardt, .conjugate_gradient=true};
  ARTS_USER_ERROR_IF(method == "li_m" || method == "gn_m",
                     "OEM method '{}' is not supported. Use 'li_cg_m' or 'gn_cg_m' for measurement-space solves.",
                     method)
  ARTS_USER_ERROR(
      "Unknown OEM method '{}'. Supported methods: li, li_cg, li_cg_m, gn, gn_cg, gn_cg_m, lm, lm_cg; aliases: ml, ml_cg.",
      method)
}

// Both solvers consume the same validated, named damping settings.
template <typename Optimizer>
void configure_lm(Optimizer& optimizer, const OEMLMSettings& settings, Numeric tolerance, unsigned int iterations) {
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
}  // namespace

void model_state_vec_aprioriFromState(Vector& xa, const Vector& x) {
  ARTS_TIME_REPORT

  xa = x;
}

void measurement_vec_fitFromMeasurement(Vector& yf, const Vector& y) {
  ARTS_TIME_REPORT

  yf = y;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void OEM(const Workspace&        ws,
         Vector&                 model_state_vec,
         Vector&                 measurement_vec_fit,
         Matrix&                 measurement_jac,
         AtmField&               atm_field,
         AbsorptionBands&        abs_bands,
         ArrayOfSensorObsel&     measurement_sensor,
         SurfaceField&           surf_field,
         SubsurfaceField&        subsurf_field,
         Matrix&                 measurement_gain_mat,
         Vector&                 oem_diagnostics,
         Vector&                 lm_ga_history,
         ArrayOfString&          errors,
         const JacobianTargets&  jac_targets,
         const Vector&           model_state_vec_apriori,
         const CovarianceMatrix& model_state_covmat,
         const Vector&           measurement_vec,
         const CovarianceMatrix& measurement_vec_error_covmat,
         const Agenda&           inversion_iterate_agenda,
         const String&           method,
         const Numeric&          max_start_cost,
         const Vector&           model_state_covmat_normalization,
         const Index&            max_iter,
         const Numeric&          stop_dx,
         const Vector&           lm_ga_settings,
         const Index&            clear_matrices,
         const Index&            display_progress) {
  ARTS_TIME_REPORT

  const OEMMethod selected = parse_oem_method(method);
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
  const std::optional<OEMLMSettings> lm_settings =
      selected.damped() ? std::optional{OEMLMSettings::from_vector(lm_ga_settings)} : std::nullopt;

  const Index n = model_state_covmat.nrows();
  const Index m = measurement_vec.size();
  model_state_covmat.compute_inverse();
  measurement_vec_error_covmat.compute_inverse();

  errors.clear();
  measurement_gain_mat.resize(0, 0);
  oem_diagnostics.resize(5);
  oem_diagnostics = NAN;
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
                                    model_state_vec,
                                    1,
                                    0,
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
  oem_diagnostics[1] = cost_start;

  // Handle cases with too large start cost
  if (max_start_cost > 0 && cost_start > max_start_cost) {
    // No inversion; retain the starting state and its simulated measurement.
    oem_diagnostics[0] = 99;
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
    bool        apply_norm = false;
    oem::Matrix T{};
    if (model_state_covmat_normalization.size() == static_cast<Size>(n)) {
      T.resize(n, n);
      static_cast<::Matrix&>(T) = 0.0;
      diagonal(T)               = model_state_covmat_normalization;
      apply_norm                = true;
    }

    oem::CovarianceMatrix Se(measurement_vec_error_covmat), Sa(model_state_covmat);
    oem::Vector           xa_oem(model_state_vec_apriori), y_oem(measurement_vec), x_oem(model_state_vec);
    oem::AgendaWrapper    aw(&ws,
                             static_cast<unsigned int>(m),
                             static_cast<unsigned int>(n),
                             measurement_jac,
                             measurement_vec_fit,
                             &atm_field,
                             &abs_bands,
                             &measurement_sensor,
                             &surf_field,
                             &subsurf_field,
                             &jac_targets,
                             &inversion_iterate_agenda);
    const auto            iterations = static_cast<unsigned int>(selected.linear() ? 1 : max_iter);
    const auto            verbosity  = static_cast<unsigned int>(display_progress);

    // Read diagnostics from the formulation that actually ran, including when
    // the forward model throws. Costs always use the same measurement scaling.
    auto run = [&]<typename Retrieval, typename Optimizer>(Retrieval& retrieval, Optimizer& optimizer) {
      auto diagnostics = [&] {
        oem_diagnostics[2] = retrieval.cost / static_cast<Numeric>(m);
        oem_diagnostics[3] = retrieval.cost_y / static_cast<Numeric>(m);
        oem_diagnostics[4] = static_cast<Numeric>(retrieval.iterations);
      };
      retrieval.iterations = 0;
      try {
        oem_diagnostics[0] = retrieval.template compute<Optimizer&, oem::ArtsLog>(
            x_oem, y_oem, optimizer, verbosity, lm_ga_history, selected.linear());
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
        damping.add_correlation_inverse(
            Block(Range(0, n),
                  Range(0, n),
                  std::make_pair(0, 0),
                  std::make_shared<Sparse>(Sparse::diagonal(model_state_covmat.inverse_diagonal()))));
        oem::CovarianceMatrix precision = inv(oem::CovarianceMatrix(damping));
        invlib::LevenbergMarquardt<Numeric, oem::CovarianceMatrix, Solver> optimizer(precision, solver);
        configure_lm(optimizer, *lm_settings, stop_dx, iterations);
        oem::OEM_STANDARD<oem::AgendaWrapper> retrieval(aw, xa_oem, Sa, Se);
        run(retrieval, optimizer);
        if (optimizer.get_lambda() > optimizer.get_lambda_maximum()) oem_diagnostics[0] = 2;
      } else {
        invlib::GaussNewton<Numeric, Solver> optimizer(stop_dx, iterations, solver);
        // Only CG supports the lazy matrix expression in the m formulation.
        if constexpr (std::is_same_v<Solver, oem::CG>) {
          if (selected.measurement_space) {
            oem::OEM_MFORM<oem::AgendaWrapper> retrieval(aw, xa_oem, Sa, Se);
            run(retrieval, optimizer);
            return;
          }
        }
        oem::OEM_STANDARD<oem::AgendaWrapper> retrieval(aw, xa_oem, Sa, Se);
        run(retrieval, optimizer);
      }
    };

    try {
      if (selected.conjugate_gradient) {
        oem::CG solver(T, apply_norm, 1e-10, 0);
        solve(solver);
      } else {
        oem::Std solver(T, apply_norm);
        solve(solver);
      }
      // invlib can stop with the Jacobian from the preceding state. Refresh it
      // so the returned gain/Jacobian describe the actual retrieved state.
      if (!selected.linear() && !clear_matrices) {
        oem::Vector fitted;
        aw.Jacobian(x_oem, fitted);
      }
    } catch (const std::exception& e) {
      oem_diagnostics[0]            = 9;
      static_cast<::Vector&>(x_oem) = NAN;
      for (const auto& message : oem::handle_nested_exception(e)) {
        std::stringstream stream{message};
        for (std::string line; std::getline(stream, line);) errors.push_back(line);
      }
    }

    model_state_vec     = x_oem;
    measurement_vec_fit = aw.get_measurement_vec();

    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices) {
      measurement_jac.resize(0, 0);
      measurement_gain_mat.resize(0, 0);
    } else if (oem_diagnostics[0] <= 2) {
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
