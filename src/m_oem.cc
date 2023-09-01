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

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include "array.h"
#include "arts_omp.h"
#include "atm.h"
#include <workspace.h>
#include "debug.h"
#include "jacobian.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "rte.h"
#include "special_interp.h"
#include "surf.h"
#include "surface.h"
#include "check_input.h"

#pragma GCC diagnostic ignored "-Wconversion"

#ifdef OEM_SUPPORT
#include "oem.h"
#endif

/* Workspace method: Doxygen documentation will be auto-generated */
void vmr_fieldClip(Tensor4& vmr_field,
                   const ArrayOfArrayOfSpeciesTag& abs_species,
                   const String& species,
                   const Numeric& limit_low,
                   const Numeric& limit_high) {
  Index iq = -1;
  if (species == "ALL") {
  }

  else {
    for (Index i = 0; i < abs_species.nelem(); i++) {
      if (abs_species[i].Species() == SpeciesTag(species).Spec()) {
        iq = i;
        break;
      }
    }
    ARTS_USER_ERROR_IF (iq < 0,
      "Could not find ", species, " in abs_species.\n")
  }

  Tensor4Clip(vmr_field, iq, limit_low, limit_high);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void xClip(Vector& x,
           const ArrayOfRetrievalQuantity& jacobian_quantities,
           const Index& ijq,
           const Numeric& limit_low,
           const Numeric& limit_high) {
  // Sizes
  const Index nq = jacobian_quantities.nelem();

  ARTS_USER_ERROR_IF (ijq < -1, "Argument *ijq* must be >= -1.");
  ARTS_USER_ERROR_IF (ijq >= nq,
      "Argument *ijq* is too high.\n"
      "You have selected index: ", ijq, "\n"
      "but the number of quantities is only: ", nq, "\n"
      "(Note that zero-based indexing is used)\n")

  // Jacobian indices
  ArrayOfArrayOfIndex ji;
  {
    bool any_affine;
    jac_ranges_indices(ji, any_affine, jacobian_quantities);
  }

  Index ifirst = 0, ilast = x.nelem() - 1;
  if (ijq > -1) {
    ifirst = ji[ijq][0];
    ilast = ji[ijq][1];
  }

  if (!std::isinf(limit_low)) {
    for (Index i = ifirst; i <= ilast; i++) {
      if (x[i] < limit_low) x[i] = limit_low;
    }
  }
  if (!std::isinf(limit_high)) {
    for (Index i = ifirst; i <= ilast; i++) {
      if (x[i] > limit_high) x[i] = limit_high;
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void x2artsSensor(const Workspace& ws,
                  Matrix& sensor_los,
                  Vector& f_backend,
                  Vector& y_baseline,
                  Sparse& sensor_response,
                  Vector& sensor_response_f,
                  ArrayOfIndex& sensor_response_pol,
                  Matrix& sensor_response_dlos,
                  Vector& sensor_response_f_grid,
                  ArrayOfIndex& sensor_response_pol_grid,
                  Matrix& sensor_response_dlos_grid,
                  Matrix& mblock_dlos,
                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                  const Vector& x,
                  const Agenda& sensor_response_agenda,
                  const Index& sensor_checked,
                  const ArrayOfTime& sensor_time) {
  // Basics
  //
  ARTS_USER_ERROR_IF (sensor_checked != 1,
        "The sensor response must be flagged to have "
        "passed a consistency check (sensor_checked=1).");

  // Revert transformation
  Vector x_t(x);
  transform_x_back(x_t, jacobian_quantities);

  // Main sizes
  const Index nq = jacobian_quantities.nelem();

  // Jacobian indices
  ArrayOfArrayOfIndex ji;
  {
    bool any_affine;
    jac_ranges_indices(ji, any_affine, jacobian_quantities, true);
  }

  // Check input
  ARTS_USER_ERROR_IF (x_t.nelem() != ji[nq - 1][1] + 1,
        "Length of *x* does not match length implied by "
        "*jacobian_quantities*.");

  // Flag indicating that y_baseline is not set
  bool yb_set = false;

  // Shall sensor responses be calculed?
  bool do_sensor = false;

  // Loop retrieval quantities
  for (Index q = 0; q < nq; q++) {
    // Index range of this retrieval quantity
    const Index np = ji[q][1] - ji[q][0] + 1;

    // Pointing off-set
    // ----------------------------------------------------------------------------
    if (jacobian_quantities[q].Target().isPointing()) {
      // Handle pointing "jitter" seperately
      if (jacobian_quantities[q].Grids()[0][0] == -1) {
        ARTS_USER_ERROR_IF (sensor_los.nrows() != np,
              "Mismatch between pointing jacobian and *sensor_los* found.");
        // Simply add retrieved off-set(s) to za column of *sensor_los*
        for (Index i = 0; i < np; i++) {
          sensor_los(i, 0) += x_t[ji[q][0] + i];
        }
      }
      // Polynomial representation
      else {
        ARTS_USER_ERROR_IF (sensor_los.nrows() != sensor_time.nelem(),
              "Sizes of *sensor_los* and *sensor_time* do not match.");
        Vector w;
        for (Index c = 0; c < np; c++) {
          polynomial_basis_func(w, time_vector(sensor_time), c);
          for (Index i = 0; i < w.nelem(); i++) {
            sensor_los(i, 0) += w[i] * x_t[ji[q][0] + c];
          }
        }
      }
    }

    // Frequency shift or stretch
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q].Target().isFrequency()) {
      if (jacobian_quantities[q] == Jacobian::Sensor::FrequencyShift) {
        ARTS_ASSERT(np == 1);
        if (x_t[ji[q][0]] != 0) {
          do_sensor = true;
          f_backend += x_t[ji[q][0]];
        }
      } else if (jacobian_quantities[q] == Jacobian::Sensor::FrequencyStretch) {
        ARTS_ASSERT(np == 1);
        if (x_t[ji[q][0]] != 0) {
          do_sensor = true;
          Vector w;
          polynomial_basis_func(w, f_backend, 1);
          for (Index i = 0; i < w.nelem(); i++) {
            f_backend[i] += w[i] * x_t[ji[q][0]];
          }
        }
      } else {
        ARTS_ASSERT(0);
      }
    }

    // Baseline fit: polynomial or sinusoidal
    // ----------------------------------------------------------------------------
    else if (jacobian_quantities[q] == Jacobian::Sensor::Polyfit ||
             jacobian_quantities[q] == Jacobian::Sensor::Sinefit) {
      if (!yb_set) {
        yb_set = true;
        Index y_size = sensor_los.nrows() * sensor_response_f_grid.nelem() *
                       sensor_response_pol_grid.nelem() *
                       sensor_response_dlos_grid.nrows();
        y_baseline.resize(y_size);
        y_baseline = 0;
      }

      for (Index mb = 0; mb < sensor_los.nrows(); ++mb) {
        calcBaselineFit(y_baseline,
                        x_t,
                        mb,
                        sensor_response,
                        sensor_response_pol_grid,
                        sensor_response_f_grid,
                        sensor_response_dlos_grid,
                        jacobian_quantities[q],
                        q,
                        ji);
      }
    }
  }

  // *y_baseline* not yet set?
  if (!yb_set) {
    y_baseline.resize(1);
    y_baseline[0] = 0;
  }

  // Recalculate sensor_response?
  if (do_sensor) {
    sensor_response_agendaExecute(ws,
                                  sensor_response,
                                  sensor_response_f,
                                  sensor_response_f_grid,
                                  sensor_response_pol,
                                  sensor_response_pol_grid,
                                  sensor_response_dlos,
                                  sensor_response_dlos_grid,
                                  mblock_dlos,
                                  f_backend,
                                  sensor_response_agenda);
  }
}


#ifdef OEM_SUPPORT
/* Workspace method: Doxygen documentation will be auto-generated */
void OEM(const Workspace& ws,
         Vector& x,
         Vector& yf,
         Matrix& jacobian,
         Matrix& dxdy,
         Vector& oem_diagnostics,
         Vector& lm_ga_history,
         ArrayOfString& errors,
         const Vector& xa,
         const CovarianceMatrix& covmat_sx,
         const Vector& y,
         const CovarianceMatrix& covmat_se,
         const ArrayOfRetrievalQuantity& jacobian_quantities,
         const Agenda& inversion_iterate_agenda,
         const String& method,
         const Numeric& max_start_cost,
         const Vector& x_norm,
         const Index& max_iter,
         const Numeric& stop_dx,
         const Vector& lm_ga_settings,
         const Index& clear_matrices,
         const Index& display_progress) {
  // Main sizes
  const Index n = covmat_sx.nrows();
  const Index m = y.nelem();

  // Checks
  covmat_sx.compute_inverse();
  covmat_se.compute_inverse();

  OEM_checks(ws,
             x,
             yf,
             jacobian,
             inversion_iterate_agenda,
             xa,
             covmat_sx,
             y,
             covmat_se,
             jacobian_quantities,
             method,
             x_norm,
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
  if (x.nelem() != n) {
    x = xa;
    yf.resize(0);
    jacobian.resize(0, 0);
  }

  // If no precomputed value given, we compute yf and jacobian to
  // compute initial cost (and use in the first OEM iteration).
  if (yf.nelem() == 0) {
    inversion_iterate_agendaExecute(
        ws, yf, jacobian, xa, 1, 0, inversion_iterate_agenda);
  }

  ARTS_USER_ERROR_IF (yf.nelem() not_eq y.nelem(),
      "Mismatch between simulated y and input y.\n"
      "Input y is size ", y.nelem(), " but simulated y is ",
      yf.nelem(), "\n"
      "Use your frequency grid vector and your sensor response matrix to match simulations with measurements.\n")

  // TODO: Get this from invlib log.
  // Start value of cost function
  Numeric cost_start = NAN;
  if (method == "ml" || method == "lm" || display_progress ||
      max_start_cost > 0) {
    Vector dy = y;
    dy -= yf;
    Vector sdy = y;
    mult_inv(ExhaustiveMatrixView{sdy}, covmat_se, ExhaustiveMatrixView{dy});
    Vector dx = x;
    dx -= xa;
    Vector sdx = x;
    mult_inv(ExhaustiveMatrixView{sdx}, covmat_sx, ExhaustiveMatrixView{dx});
    cost_start = dx * sdx + dy * sdy;
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
           << "        Set limit : " << max_start_cost << std::endl
           << "      Found value : " << cost_start << std::endl
           << std::endl;
    }
  }
  // Otherwise do inversion
  else {
    bool apply_norm = false;
    oem::Matrix T{};
    if (x_norm.nelem() == n) {
      T.resize(n, n);
      T *= 0.0;
      T.diagonal() = x_norm;
      for (Index i = 0; i < n; i++) {
        T(i, i) = x_norm[i];
      }
      apply_norm = true;
    }

    oem::CovarianceMatrix Se(covmat_se), Sa(covmat_sx);
    oem::Vector xa_oem(xa), y_oem(y), x_oem(x);
    oem::AgendaWrapper aw(const_cast<Workspace*>(&ws), // FIXME: WHAT IS HAPPENING HERE????
                          (unsigned int)m,
                          (unsigned int)n,
                          jacobian,
                          yf,
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
        ARTS_USER_ERROR(method, " is not supported")
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
        ARTS_USER_ERROR(method, " is not supported")
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
        Sparse diagonal = Sparse::diagonal(covmat_sx.inverse_diagonal());
        CovarianceMatrix SaDiag{};
        SaDiag.add_correlation_inverse(Block(Range(0, n),
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

        Sparse diagonal = Sparse::diagonal(covmat_sx.inverse_diagonal());
        CovarianceMatrix SaDiag{};
        SaDiag.add_correlation_inverse(Block(Range(0, n),
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
      oem_diagnostics[0] = 9;
      oem_diagnostics[2] = oem.cost;
      oem_diagnostics[3] = oem.cost_y;
      oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);
      x_oem *= NAN;
      std::vector<std::string> sv = oem::handle_nested_exception(e);
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

    x = x_oem;
    yf = aw.get_measurement_vector();

    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices) {
      jacobian.resize(0, 0);
      dxdy.resize(0, 0);
    } else if (oem_diagnostics[0] <= 2) {
      dxdy.resize(n, m);
      Matrix tmp1(n, m), tmp2(n, n), tmp3(n, n);
      mult_inv(tmp1, transpose(jacobian), covmat_se);
      mult(tmp2, tmp1, jacobian);
      add_inv(tmp2, covmat_sx);
      inv(tmp3, tmp2);
      mult(dxdy, tmp3, tmp1);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void covmat_soCalc(Matrix& covmat_so,
                   const Matrix& dxdy,
                   const CovarianceMatrix &covmat_se) {
  Index n(dxdy.nrows()), m(dxdy.ncols());
  Matrix tmp1(m, n);

  ARTS_USER_ERROR_IF ((m == 0) || (n == 0),
        "The gain matrix *dxdy* is required to compute the observation error covariance matrix.");
  ARTS_USER_ERROR_IF ((covmat_se.nrows() != m) || (covmat_se.ncols() != m),
        "The covariance matrix covmat_se has invalid dimensions.");

  covmat_so.resize(n, n);
  mult(tmp1, covmat_se, transpose(dxdy));
  mult(covmat_so, dxdy, tmp1);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void covmat_ssCalc(Matrix& covmat_ss,
                   const Matrix& avk,
                   const CovarianceMatrix &covmat_sx) {
  Index n(avk.ncols());
  Matrix tmp1(n, n), tmp2(n, n);

  ARTS_USER_ERROR_IF (n == 0,
        "The averaging kernel matrix *dxdy* is required to compute the smoothing error covariance matrix.");
  ARTS_USER_ERROR_IF ((covmat_sx.nrows() != n) || (covmat_sx.ncols() != n),
        "The covariance matrix *covmat_sx* invalid dimensions.");

  covmat_ss.resize(n, n);

  // Sign doesn't matter since we're dealing with a quadratic form.
  id_mat(tmp1);
  tmp1 -= avk;

  mult(tmp2, covmat_sx, transpose(tmp1));
  mult(covmat_ss, tmp1, tmp2);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixFromCovarianceMatrix(Matrix& S,
                                const CovarianceMatrix &Sc) {
  S = Matrix(Sc);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void avkCalc(Matrix& avk,
             const Matrix& dxdy,
             const Matrix &jacobian) {
  Index m(jacobian.nrows()), n(jacobian.ncols());
  ARTS_USER_ERROR_IF ((m == 0) || (n == 0),
                      "The Jacobian matrix is empty.");
  ARTS_USER_ERROR_IF ((dxdy.nrows() != n) || (dxdy.ncols() != m),
      "Matrices have inconsistent sizes.\n"
      "  Size of gain matrix: ", dxdy.nrows(), " x ", dxdy.ncols(),
      "\n"
      "     Size of Jacobian: ", jacobian.nrows(), " x ",
      jacobian.ncols(), "\n")

  avk.resize(n, n);
  mult(avk, dxdy, jacobian);
}

#else

void covmat_soCalc(Matrix& /* covmat_so */,
                   const Matrix& /* dxdy */,
                   const CovarianceMatrix& /* covmat_ /*v*/) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}

void covmat_ssCalc(Matrix& /*covmat_ss*/,
                   const Matrix& /*avk*/,
                   const CovarianceMatrix& /*covmat_ /*v*/) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}

void avkCalc(Matrix& /* avk */,
             const Matrix& /* dxdy */,
             const Matrix& /* jacobia /*v*/) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}

void OEM(Workspace&,
         Vector&,
         Vector&,
         Matrix&,
         Matrix&,
         Vector&,
         Vector&,
         ArrayOfString&,
         const Vector&,
         const CovarianceMatrix&,
         const Vector&,
         const CovarianceMatrix&,
         const Index&,
         const ArrayOfRetrievalQuantity&,
         const ArrayOfArrayOfIndex&,
         const Agenda&,
         const String&,
         const Numeric&,
         const Vector&,
         const Index&,
         const Numeric&,
         const Vector&,
         const Index&,
         const Index&) {
  ARTS_USER_ERROR (
      "WSM is not available because ARTS was compiled without "
      "OEM support.");
}
#endif

#if defined(OEM_SUPPORT) && 0

#include "agenda_wrapper_mpi.h"
#include "oem_mpi.h"

//
// Performs manipulations of workspace variables necessary for distributed
// retrievals with MPI:
//
//   - Splits up sensor positions evenly over processes
//   - Splits up inverse covariance matrices.
//
void MPI_Initialize(Matrix& sensor_los,
                    Matrix& sensor_pos,
                    Vector& sensor_time) {
  int initialized;

  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(nullptr, nullptr);
  }

  int rank, nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int nmblock = (int)sensor_pos.nrows();
  int mblock_range = nmblock / nprocs;
  int mblock_start = mblock_range * rank;
  int remainder = nmblock % std::max(mblock_range, nprocs);

  //
  // Split up sensor positions.
  //

  if (rank < remainder) {
    mblock_range += 1;
    mblock_start += rank;
  } else {
    mblock_start += remainder;
  }

  if (nmblock > 0) {
    Range range = Range(mblock_start, mblock_range);

    Matrix tmp_m = sensor_los(range, joker);
    sensor_los = tmp_m;

    tmp_m = sensor_pos(range, joker);
    sensor_pos = tmp_m;

    Vector tmp_v = sensor_time[range];
    sensor_time = tmp_v;
  } else {
    sensor_los.resize(0, 0);
    sensor_pos.resize(0, 0);
    sensor_time.resize(0);
  }
}

void OEM_MPI(const Workspace& ws,
             Vector& x,
             Vector& yf,
             Matrix& jacobian,
             Matrix& dxdy,
             Vector& oem_diagnostics,
             Vector& lm_ga_history,
             Matrix& sensor_los,
             Matrix& sensor_pos,
             Vector& sensor_time,
             CovarianceMatrix& covmat_sx,
             CovarianceMatrix& covmat_se,
             const Vector& xa,
             const Vector& y,
             const ArrayOfRetrievalQuantity& jacobian_quantities,
             const Agenda& inversion_iterate_agenda,
             const String& method,
             const Numeric& max_start_cost,
             const Vector& x_norm,
             const Index& max_iter,
             const Numeric& stop_dx,
             const Vector& lm_ga_settings,
             const Index& clear_matrices,
             const Index /*v*/) {
  // Main sizes
  const Index n = covmat_sx.nrows();
  const Index m = y.nelem();

  // Check WSVs
  OEM_checks(ws,
             x,
             yf,
             jacobian,
             inversion_iterate_agenda,
             xa,
             covmat_sx,
             covmat_se,
             jacobian_quantities,
             method,
             x_norm,
             max_iter,
             stop_dx,
             lm_ga_settings,
             clear_matrices,
             display_progress);

  // Calculate spectrum and Jacobian for a priori state
  // Jacobian is also input to the agenda, and to flag this is this first
  // call, this WSV must be set to be empty.
  jacobian.resize(0, 0);

  // Initialize MPI environment.
  MPI_Initialize(sensor_los, sensor_pos, sensor_time);

  // Setup distributed matrices.
  MPICovarianceMatrix SeInvMPI(covmat_se);
  MPICovarianceMatrix SaInvMPI(covmat_sx);

  // Create temporary MPI vector from local results and use conversion to
  // standard vector to broadcast results to all processes.
  oem::Vector tmp;
  inversion_iterate_agendaExecute(
      ws, tmp, jacobian, xa, 1, inversion_iterate_agenda);
  yf = MPIVector(tmp);

  // Size diagnostic output and init with NaNs
  oem_diagnostics.resize(5);
  oem_diagnostics = NAN;
  //
  if (method == "ml" || method == "lm") {
    lm_ga_history.resize(max_iter);
    lm_ga_history = NAN;
  } else {
    lm_ga_history.resize(0);
  }

  // Start value of cost function. Covariance matrices are already distributed
  // over processes, so we need to use invlib matrix algebra.
  Numeric cost_start = NAN;
  if (method == "ml" || method == "lm" || display_progress ||
      max_start_cost > 0) {
    oem::Vector dy = y;
    dy -= yf;
    cost_start = dot(dy, SeInvMPI * dy);
  }
  oem_diagnostics[1] = cost_start;

  // Handle cases with too large start cost
  if (max_start_cost > 0 && cost_start > max_start_cost) {
    // Flag no inversion in oem_diagnostics, and let x to be undefined
    oem_diagnostics[0] = 99;
    //
    if (display_progress) {
      cout << "\n   No OEM inversion, too high start cost:\n"
           << "        Set limit : " << max_start_cost << endl
           << "      Found value : " << cost_start << endl
           << endl;
    }
  }

  // Otherwise do inversion
  else {
    // Size remaining output arguments
    x.resize(n);
    dxdy.resize(n, m);

    OEMVector xa_oem(xa), y_oem(y), x_oem;
    oem::AgendaWrapperMPI aw(&ws, &inversion_iterate_agenda, m, n);

    OEM_PS_PS_MPI<AgendaWrapperMPI> oem(aw, xa_oem, SaInvMPI, SeInvMPI);

    // Call selected method
    int return_value = 99;

    if (method == "li") {
      oem::CG cg(1e-12, 0);
      oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
      return_value = oem.compute<oem::GN_CG, invlib::MPILog>(
          x_oem, y_oem, gn, 2 * (int)display_progress);
    } else if (method == "gn") {
	    oem::CG cg(1e-12, 0);
      oem::GN_CG gn(stop_dx, (unsigned int)max_iter, cg);
      return_value = oem.compute<oem::GN_CG, invlib::MPILog>(
          x_oem, y_oem, gn, 2 * (int)display_progress);
    } else if ((method == "lm") || (method == "ml")) {
	    oem::CG cg(1e-12, 0);
      LM_CG_MPI lm(SaInvMPI, cg);

      lm.set_tolerance(stop_dx);
      lm.set_maximum_iterations((unsigned int)max_iter);
      lm.set_lambda(lm_ga_settings[0]);
      lm.set_lambda_decrease(lm_ga_settings[1]);
      lm.set_lambda_increase(lm_ga_settings[2]);
      lm.set_lambda_threshold(lm_ga_settings[3]);
      lm.set_lambda_maximum(lm_ga_settings[4]);

      return_value = oem.compute<oem::LM_CG_MPI, invlib::MPILog>(
          x_oem, y_oem, lm, 2 * (int)display_progress);
    }

    oem_diagnostics[0] = return_value;
    oem_diagnostics[2] = oem.cost;
    oem_diagnostics[3] = oem.cost_y;
    oem_diagnostics[4] = static_cast<Numeric>(oem.iterations);

    x = x_oem;
    // Shall empty jacobian and dxdy be returned?
    if (clear_matrices && (oem_diagnostics[0])) {
      jacobian.resize(0, 0);
      dxdy.resize(0, 0);
    }
  }
  MPI_Finalize();
}

#else

void OEM_MPI(Workspace&,
             Vector&,
             Vector&,
             Matrix&,
             Matrix&,
             Vector&,
             Vector&,
             Matrix&,
             Matrix&,
             Vector&,
             CovarianceMatrix&,
             CovarianceMatrix&,
             const Vector&,
             const Vector&,
             const Index&,
             const ArrayOfRetrievalQuantity&,
             const Agenda&,
             const String&,
             const Numeric&,
             const Vector&,
             const Index&,
             const Numeric&,
             const Vector&,
             const Index&,
             const Index&) {
  ARTS_USER_ERROR (
      "You have to compile ARTS with OEM support "
      " and enable MPI to use OEM_MPI.");
}

#endif  // OEM_SUPPORT && ENABLE_MPI
