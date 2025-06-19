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

#include <workspace.h>

#include <cmath>
#include <string>

#include "array.h"
#include "atm.h"
#include "config.h"
#include "debug.h"
#include "jacobian.h"

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef OEM_SUPPORT
#include "oem.h"
#endif

void model_state_vector_aprioriFromState(Vector& xa, const Vector& x) {
  ARTS_TIME_REPORT

  xa = x;
}

void measurement_vector_fittedFromMeasurement(Vector& yf, const Vector& y) {
  ARTS_TIME_REPORT

  yf = y;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void OEM(const Workspace& ws,
         Vector& model_state_vector,
         Vector& measurement_vector_fitted,
         Matrix& measurement_jacobian,
         AtmField& atmospheric_field,
         AbsorptionBands& absorption_bands,
         ArrayOfSensorObsel& measurement_sensor,
         SurfaceField& surface_field,
         SubsurfaceField& subsurface_field,
         Matrix& measurement_gain_matrix,
         Vector& oem_diagnostics,
         Vector& lm_ga_history,
         ArrayOfString& errors,
         const JacobianTargets& jacobian_targets,
         const Vector& model_state_vector_apriori,
         const CovarianceMatrix& model_state_covariance_matrix,
         const Vector& measurement_vector,
         const CovarianceMatrix& measurement_vector_error_covariance_matrix,
         const Agenda& inversion_iterate_agenda,
         const String& method,
         const Numeric& max_start_cost,
         const Vector& model_state_covariance_matrix_normalization,
         const Index& max_iter,
         const Numeric& stop_dx,
         const Vector& lm_ga_settings,
         const Index& clear_matrices,
         const Index& display_progress) {
  ARTS_TIME_REPORT

  // Main sizes
  const Index n = model_state_covariance_matrix.nrows();
  const Index m = measurement_vector.size();

  // Checks
  model_state_covariance_matrix.compute_inverse();
  measurement_vector_error_covariance_matrix.compute_inverse();

  OEM_checks(ws,
             model_state_vector,
             measurement_vector_fitted,
             measurement_jacobian,
             atmospheric_field,
             absorption_bands,
             measurement_sensor,
             surface_field,
             subsurface_field,
             jacobian_targets,
             inversion_iterate_agenda,
             model_state_vector_apriori,
             model_state_covariance_matrix,
             measurement_vector,
             measurement_vector_error_covariance_matrix,
             method,
             model_state_covariance_matrix_normalization,
             max_iter,
             stop_dx,
             lm_ga_settings,
             clear_matrices,
             display_progress);

  // Size diagnostic output and init with NaNs
  oem_diagnostics.resize(5);
  oem_diagnostics = NAN;
  //
  if (method == "ml" || method == "lm" || method == "ml_cg" ||
      method == "lm_cg") {
    lm_ga_history.resize(max_iter + 1);
    lm_ga_history = NAN;
  } else {
    lm_ga_history.resize(0);
  }

  // Check for start vector and precomputed yf, jacobian
  if (model_state_vector.size() != static_cast<Size>(n)) {
    model_state_vector = model_state_vector_apriori;
    measurement_vector_fitted.resize(0);
    measurement_jacobian.resize(0, 0);
  }

  // If no precomputed value given, we compute yf and jacobian to
  // compute initial cost (and use in the first OEM iteration).
  if (measurement_vector_fitted.size() == 0) {
    inversion_iterate_agendaExecute(ws,
                                    atmospheric_field,
                                    absorption_bands,
                                    measurement_sensor,
                                    surface_field,
                                    subsurface_field,
                                    measurement_vector_fitted,
                                    measurement_jacobian,
                                    jacobian_targets,
                                    model_state_vector_apriori,
                                    1,
                                    0,
                                    inversion_iterate_agenda);
  }

  ARTS_USER_ERROR_IF(
      measurement_vector_fitted.size() not_eq measurement_vector.size(),
      "Mismatch between simulated y and input y.\n"
      "Input y is size {}"
      " but simulated y is {}"
      "\n"
      "Use your frequency grid vector and your sensor response matrix to match simulations with measurements.\n",
      measurement_vector.size(),
      measurement_vector_fitted.size())

  // TODO: Get this from invlib log.
  // Start value of cost function
  Numeric cost_start = NAN;
  if (method == "ml" || method == "lm" || display_progress ||
      max_start_cost > 0) {
    Vector dy   = measurement_vector;
    dy         -= measurement_vector_fitted;
    Vector sdy  = measurement_vector;
    mult_inv(sdy.view_as(sdy.size(), 1),
             measurement_vector_error_covariance_matrix,
             dy.view_as(dy.size(), 1));
    Vector dx   = model_state_vector;
    dx         -= model_state_vector_apriori;
    Vector sdx  = model_state_vector;
    mult_inv(sdx.view_as(sdx.size(), 1),
             model_state_covariance_matrix,
             dx.view_as(dx.size(), 1));
    cost_start  = dot(dx, sdx) + dot(dy, sdy);
    cost_start /= static_cast<Numeric>(m);
  }
  oem_diagnostics[1] = cost_start;

  // Handle cases with too large start cost
  if (max_start_cost > 0 && cost_start > max_start_cost) {
    // Flag no inversion in oem_diagnostics, and let x to be undefined
    oem_diagnostics[0] = 99;
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
    bool apply_norm = false;
    oem::Matrix T{};
    if (model_state_covariance_matrix_normalization.size() ==
        static_cast<Size>(n)) {
      T.resize(n, n);
      T           *= 0.0;
      diagonal(T)  = model_state_covariance_matrix_normalization;
      for (Index i = 0; i < n; i++) {
        T(i, i) = model_state_covariance_matrix_normalization[i];
      }
      apply_norm = true;
    }

    oem::CovarianceMatrix Se(measurement_vector_error_covariance_matrix),
        Sa(model_state_covariance_matrix);
    oem::Vector xa_oem(model_state_vector_apriori), y_oem(measurement_vector),
        x_oem(model_state_vector);
    oem::AgendaWrapper aw(&ws,
                          (unsigned int)m,
                          (unsigned int)n,
                          measurement_jacobian,
                          measurement_vector_fitted,
                          &atmospheric_field,
                          &absorption_bands,
                          &measurement_sensor,
                          &surface_field,
                          &subsurface_field,
                          &jacobian_targets,
                          &inversion_iterate_agenda);
    oem::OEM_STANDARD<oem::AgendaWrapper> oem(aw, xa_oem, Sa, Se);
    oem::OEM_MFORM<oem::AgendaWrapper> oem_m(aw, xa_oem, Sa, Se);
    int oem_verbosity = static_cast<int>(display_progress);

    int return_code = 0;

    try {
      if (method == "li") {
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, 1, s);  // Linear case, only one step.
        return_code = oem.compute<oem::GN, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "li_m") {
        ARTS_USER_ERROR("{} is not supported", method)
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, 1, s);  // Linear case, only one step.
        //        return_code = oem_m.compute<oem::GN, oem::ArtsLog>(
        //            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "li_cg") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, 1, cg);  // Linear case, only one step.
        return_code = oem.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "li_cg_m") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, 1, cg);  // Linear case, only one step.
        return_code = oem_m.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history, true);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn") {
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, (unsigned int)max_iter, s);
        return_code = oem.compute<oem::GN, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn_m") {
        ARTS_USER_ERROR("{} is not supported", method)
        oem::Std s(T, apply_norm);
        oem::GN gn(stop_dx, (unsigned int)max_iter, s);
        //        return_code = oem_m.compute<oem::GN, oem::ArtsLog>(
        //            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn_cg") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
        return_code = oem.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if (method == "gn_cg_m") {
        oem::CG cg(T, apply_norm, 1e-10, 0);
        oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
        return_code = oem_m.compute<oem::GN_CG, oem::ArtsLog>(
            x_oem, y_oem, gn, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
      } else if ((method == "lm") || (method == "ml")) {
        oem::Std s(T, apply_norm);
        Sparse diagonal =
            Sparse::diagonal(model_state_covariance_matrix.inverse_diagonal());
        CovarianceMatrix SaDiag{};
        SaDiag.add_correlation_inverse(
            Block(Range(0, n),
                  Range(0, n),
                  std::make_pair(0, 0),
                  std::make_shared<Sparse>(diagonal)));
        oem::CovarianceMatrix SaInvLM = inv(oem::CovarianceMatrix(SaDiag));
        oem::LM lm(SaInvLM, s);

        lm.set_tolerance(stop_dx);
        lm.set_maximum_iterations((unsigned int)max_iter);
        lm.set_lambda(lm_ga_settings[0]);
        lm.set_lambda_decrease(lm_ga_settings[1]);
        lm.set_lambda_increase(lm_ga_settings[2]);
        lm.set_lambda_maximum(lm_ga_settings[3]);
        lm.set_lambda_threshold(lm_ga_settings[4]);
        lm.set_lambda_constraint(lm_ga_settings[5]);

        return_code = oem.compute<oem::LM&, oem::ArtsLog>(
            x_oem, y_oem, lm, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
        if (lm.get_lambda() > lm.get_lambda_maximum()) {
          oem_diagnostics[0] = 2;
        }
      } else if ((method == "lm_cg") || (method == "ml_cg")) {
        oem::CG cg(T, apply_norm, 1e-10, 0);

        Sparse diagonal =
            Sparse::diagonal(model_state_covariance_matrix.inverse_diagonal());
        CovarianceMatrix SaDiag{};
        SaDiag.add_correlation_inverse(
            Block(Range(0, n),
                  Range(0, n),
                  std::make_pair(0, 0),
                  std::make_shared<Sparse>(diagonal)));
        oem::LM_CG lm(SaDiag, cg);

        lm.set_maximum_iterations((unsigned int)max_iter);
        lm.set_lambda(lm_ga_settings[0]);
        lm.set_lambda_decrease(lm_ga_settings[1]);
        lm.set_lambda_increase(lm_ga_settings[2]);
        lm.set_lambda_threshold(lm_ga_settings[3]);
        lm.set_lambda_maximum(lm_ga_settings[4]);

        return_code = oem.compute<oem::LM_CG&, oem::ArtsLog>(
            x_oem, y_oem, lm, oem_verbosity, lm_ga_history);
        oem_diagnostics[0] = static_cast<Index>(return_code);
        if (lm.get_lambda() > lm.get_lambda_maximum()) {
          oem_diagnostics[0] = 2;
        }
      }

      oem_diagnostics[2] = oem.cost / static_cast<Numeric>(m);
      oem_diagnostics[3] = oem.cost_y / static_cast<Numeric>(m);
      oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);
    } catch (const std::exception& e) {
      oem_diagnostics[0]           = 9;
      oem_diagnostics[2]           = oem.cost;
      oem_diagnostics[3]           = oem.cost_y;
      oem_diagnostics[4]           = static_cast<Numeric>(oem.iterations);
      x_oem                       *= NAN;
      std::vector<std::string> sv  = oem::handle_nested_exception(e);
      for (auto& s : sv) {
        std::stringstream ss{s};
        std::string t{};
        while (std::getline(ss, t)) {
          errors.push_back(t.c_str());
        }
      }
    } catch (...) {
      throw;
    }

    model_state_vector        = x_oem;
    measurement_vector_fitted = aw.get_measurement_vector();

    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices) {
      measurement_jacobian.resize(0, 0);
      measurement_gain_matrix.resize(0, 0);
    } else if (oem_diagnostics[0] <= 2) {
      measurement_gain_matrix.resize(n, m);
      Matrix tmp1(n, m), tmp2(n, n), tmp3(n, n);
      mult_inv(tmp1,
               transpose(measurement_jacobian),
               measurement_vector_error_covariance_matrix);
      mult(tmp2, tmp1, measurement_jacobian);
      add_inv(tmp2, model_state_covariance_matrix);
      inv(tmp3, tmp2);
      mult(measurement_gain_matrix, tmp3, tmp1);
    }
  }
}

void measurement_vector_error_covariance_matrix_observation_systemCalc(
    Matrix& measurement_vector_error_covariance_matrix_observation_system,
    const Matrix& measurement_gain_matrix,
    const CovarianceMatrix& measurement_vector_error_covariance_matrix) {
  ARTS_TIME_REPORT

  Index n(measurement_gain_matrix.nrows()), m(measurement_gain_matrix.ncols());
  Matrix tmp1(m, n);

  ARTS_USER_ERROR_IF(
      (m == 0) || (n == 0),
      "The gain matrix *measurement_gain_matrix* is required to compute the observation error covariance matrix.");
  ARTS_USER_ERROR_IF(
      (measurement_vector_error_covariance_matrix.nrows() != m) ||
          (measurement_vector_error_covariance_matrix.ncols() != m),
      "The covariance matrix measurement_vector_error_covariance_matrix has invalid dimensions.");

  measurement_vector_error_covariance_matrix_observation_system.resize(n, n);
  mult(tmp1,
       measurement_vector_error_covariance_matrix,
       transpose(measurement_gain_matrix));
  mult(measurement_vector_error_covariance_matrix_observation_system,
       measurement_gain_matrix,
       tmp1);
}

void model_state_covariance_matrix_smoothing_errorCalc(
    Matrix& model_state_covariance_matrix_smoothing_error,
    const Matrix& measurement_averaging_kernel,
    const CovarianceMatrix& model_state_covariance_matrix) {
  ARTS_TIME_REPORT

  Index n(measurement_averaging_kernel.ncols());
  Matrix tmp1(n, n), tmp2(n, n);

  ARTS_USER_ERROR_IF(
      n == 0,
      "The averaging kernel matrix *measurement_gain_matrix* is required to compute the smoothing error covariance matrix.");
  ARTS_USER_ERROR_IF(
      (model_state_covariance_matrix.nrows() != n) ||
          (model_state_covariance_matrix.ncols() != n),
      "The covariance matrix *model_state_covariance_matrix* invalid dimensions.");

  model_state_covariance_matrix_smoothing_error.resize(n, n);

  // Sign doesn't matter since we're dealing with a quadratic form.
  id_mat(tmp1);
  tmp1 -= measurement_averaging_kernel;

  mult(tmp2, model_state_covariance_matrix, transpose(tmp1));
  mult(model_state_covariance_matrix_smoothing_error, tmp1, tmp2);
}

void measurement_averaging_kernelCalc(Matrix& measurement_averaging_kernel,
                                      const Matrix& measurement_gain_matrix,
                                      const Matrix& measurement_jacobian) {
  ARTS_TIME_REPORT

  Index m(measurement_jacobian.nrows()), n(measurement_jacobian.ncols());

  ARTS_USER_ERROR_IF(measurement_jacobian.empty(),
                     "The Jacobian matrix is empty.");

  ARTS_USER_ERROR_IF((measurement_gain_matrix.shape() != std::array{n, m}),
                     R"(Matrices have inconsistent sizes.

measurement_gain_matrix: {:B,},
measurement_jacobian:    {:B,}
)",
                     measurement_gain_matrix.shape(),
                     measurement_jacobian.shape());

  measurement_averaging_kernel.resize(n, n);
  mult(measurement_averaging_kernel,
       measurement_gain_matrix,
       measurement_jacobian);
}
