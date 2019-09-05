/**
  @file   oem.h
  @author Simon Pfreundschuh <simonpf@chalmers.se>
  @date   Fri March 25 15:53:54 2016

  @brief Defines the ARTS interface to the invlib library.

  Since invlib is a template library, the interface is defined mostly
  through type definitions that instantiate the generic invlib classes
  with the corresponding ARTS types.
*/
#ifndef _ARTS_OEM_H_
#define _ARTS_OEM_H_

#include <type_traits>

#include "invlib/algebra.h"
#include "invlib/algebra/precision_matrix.h"
#include "invlib/algebra/solvers.h"
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/map.h"
#include "invlib/optimization.h"
#include "invlib/profiling/timer.h"

////////////////////////////////////////////////////////////////////////////////
//  Type Aliases
////////////////////////////////////////////////////////////////////////////////

namespace oem {

/** invlib wrapper type for ARTS vectors.*/
using Vector = invlib::Vector<ArtsVector>;
/** invlib wrapper type for ARTS matrices.*/
using Matrix = invlib::Matrix<ArtsMatrix>;
/** invlib wrapper type for ARTS matrices to be passed by reference.*/
using MatrixReference = invlib::Matrix<ArtsMatrixReference<::Matrix>>;
/** invlib wrapper type for ARTS the ARTS covariance class.*/
using CovarianceMatrix = invlib::Matrix<ArtsCovarianceMatrixWrapper>;
using Identity = invlib::MatrixIdentity<Matrix>;

////////////////////////////////////////////////////////////////////////////////
// OEM Formulations
////////////////////////////////////////////////////////////////////////////////

/** Different formulations of the OEM optimization according to Rodgers 2000.*/
using invlib::Formulation;

/** OEM standard form.
 *
 * Class template implementing the standard form of the OEM
 * optimization defined according to Chapter 4 in Rodgers (2000).
 *
 * In this formulation, each iteration requires the solution
 * of a system of linear equations of size n-times-n.
 */
template <typename ForwardModel>
using OEM_STANDARD = invlib::MAP<ForwardModel,
                                 Matrix,
                                 CovarianceMatrix,
                                 CovarianceMatrix,
                                 Vector,
                                 Formulation::STANDARD,
                                 invlib::Rodgers531>;

/** OEM n form.
 *
 * Class template implementing the n-form of the OEM
 * optimization defined according to Chapter 4 in Rodgers (2000).
 *
 * In this formulation, each iteration requires the solution
 * of a system of linear equations of size n-times-n.
 */
template <typename ForwardModel>
using OEM_NFORM = invlib::MAP<ForwardModel,
                              Matrix,
                              CovarianceMatrix,
                              CovarianceMatrix,
                              Vector,
                              Formulation::NFORM>;

/** OEM m form.
 *
 * Class template implementing the m-form of the OEM
 * optimization defined according to Chapter 4 in Rodgers (2000).
 *
 * In this formulation, each iteration requires the solution
 * of a system of linear equations of size m-times-m.
 */
template <typename ForwardModel>
using OEM_MFORM = invlib::MAP<ForwardModel,
                              Matrix,
                              CovarianceMatrix,
                              CovarianceMatrix,
                              Vector,
                              Formulation::MFORM>;

////////////////////////////////////////////////////////////////////////////////
// Solvers
////////////////////////////////////////////////////////////////////////////////

/** Normalizing solver.
 * 
 * Solver class that wraps around a given solver and transforms the linear
 * system from left and right with the given transformation matrix. This is used
 * to implement the normalization from qpack.
 * 
 * @tparam TransformationMatrixType The type of the transformation matrix.
 * @tparam SolverType The underlying solver type used to solve the linear
 * system.
 */
template <typename TransformationMatrixType,
          typename SolverType = invlib::Standard>
class NormalizingSolver : SolverType {
 public:
  template <typename... Params>
  NormalizingSolver(const TransformationMatrixType &trans,
                    bool apply,
                    Params... params)
      : SolverType(params...), apply_(apply), trans_(trans) {}

  /** Solve linear system.
   *
   * Solves the transformed linear system using the
   * solve(...) method of the underlying solver type.
   *
   * @param[in] A Matrix defining the  linear system.
   * @param[in] v RHS vector of the linear system.
   *
   * @return The solution vector of the linear system.
   */
  template <typename MatrixType, typename VectorType>
  auto solve(const MatrixType &A, const VectorType &v) ->
      typename VectorType::ResultType {
    typename VectorType::ResultType w;
    if (apply_) {
      typename VectorType::ResultType vv = trans_ * v;
      auto &&ww = SolverType::solve(trans_ * A * trans_, vv);
      w = trans_ * ww;
    } else {
      w = SolverType::solve(A, v);
      VectorType u = v - A * w;
    }
    return w;
  }

 private:
  /** Whether or not to apply the transformation.*/
  const bool apply_ = false;
  /** The transformation matrix.*/
  const TransformationMatrixType &trans_;
};

/** The invlib standard solver
 *
 * This solver uses the built-in ARTS QR solver to solve a
 * given linear system.
 */
using Std = NormalizingSolver<Matrix, invlib::Standard>;

/** The invlib CG solver.
 *
 *  The invlib Conjugate Grdient (CG) solver. The solver only
 *  performs matrix-vector multiplication and is therefore better
 *  suited for large linear systems.
 */
using CG = NormalizingSolver<Matrix, invlib::ConjugateGradient<>>;

/** OEM Gauss-Newton optimization using normed ARTS QR solver.*/
using GN = invlib::GaussNewton<Numeric, Std>;
/** Gauss-Newton (GN) optimization using normed CG solver.*/
using GN_CG = invlib::GaussNewton<Numeric, CG>;
/** Levenberg-Marquardt (LM) optimization using normed ARTS QR solver.*/
using LM = invlib::LevenbergMarquardt<Numeric, CovarianceMatrix, Std>;
/** Levenberg-Marquardt (LM) optimization using normed CG solver.*/
using LM_CG = invlib::LevenbergMarquardt<Numeric, CovarianceMatrix, CG>;

////////////////////////////////////////////////////////////////////////////////
//  Custom Log Class
////////////////////////////////////////////////////////////////////////////////

/** Log customization for different optimization methods.
 *
 * This type trait is used to customize the log for different optimizers.
 * Its purpose is to allow the log to look different for the Gauss-Newton
 * method and the Levenberg-Marquardt method.
 */
template <typename T>
struct OptimizerLog;

/** Log customization for LM method
 *
 * This essentially adds a line for the gamma parameter to
 * the output.
 */
template <typename RealType, typename DampingMatrix, typename Solver>
struct OptimizerLog<
    invlib::LevenbergMarquardt<RealType, DampingMatrix, Solver>> {
  /** Method name */
  static constexpr auto name = "Levenberg-Marquardt";

  /** Name to append to header line. */
  static std::string header() {
    std::string out = "Gamma Factor";
    return out;
  }

  /** Returns the string to append to the log of a single step. */
  static std::string log(
      const invlib::LevenbergMarquardt<RealType, DampingMatrix, Solver> &g,
      Vector &gamma_history_,
      size_t i) {
    std::string lambda = std::to_string(g.get_lambda());
    std::string out(15 - std::min<size_t>(lambda.size(), 15), ' ');
    out += lambda;
    gamma_history_[i] = g.get_lambda();
    return out;
  }
};

/** Log customization for GN method. */
template <typename RealType, typename Solver>
struct OptimizerLog<invlib::GaussNewton<RealType, Solver>> {
  /** Method name */
  static constexpr auto name = "Gauss-Newton";

  /** Name to append to header line. */
  static std::string header() { return ""; }

  static std::string log(const invlib::GaussNewton<RealType, Solver> &,
                         Vector &,
                         size_t) {
    return "";
  }
};

/** OEM log output
 *
 * This class takes care of formatting the OEM iteration information
 * and displaying in the command line.
 *
 * @tparam The invlib log type defining which type of logging to perform.
 */
template <invlib::LogType type>
class ArtsLog {
 public:
  /** Create log.
   * 
   * @param verbosity Verbosity level 0 for silent, 2 for verbose
   * @param gamma_history Reference to vector in which to store gamma
   * values of LM iteration
   * @param linear Flag indicating whether forward model is linear.
   */
  ArtsLog(unsigned int v, ::Vector &g, bool l = false)
      : verbosity_(v), gamma_history_(g), linear_(l), finalized_(false) {}

  /** Finalizes log output if necessary.*/
  ~ArtsLog() {
    if ((verbosity_ >= 1) && (!finalized_)) {
      std::cout << invlib::separator() << std::endl << std::endl;
      std::cout << "Error during OEM computation." << std::endl;
      std::cout << std::endl;
      std::cout << invlib::center("----") << std::endl;
      std::cout << std::endl;
    }
  }

  /** Initialize log output.
   *
   * This function is called from within invlib to initialize the log
   * output. Prints general information on the OEM settings.
   */
  template <typename... Params>
  void init(Params &... params) {
    if (verbosity_ >= 1) {
      std::tuple<Params &...> tuple(params...);

      auto &y = std::get<4>(tuple);
      scaling_factor_ = 1.0 / static_cast<Numeric>(y.nelem());
      std::cout << std::endl;
      std::cout << invlib::center("MAP Computation") << std::endl;

      // Print formulation.
      int formulation = static_cast<int>(std::get<6>(tuple));
      switch (formulation) {
        case 0:
          std::cout << "Formulation: Standard" << std::endl;
          break;
        case 1:
          std::cout << "Formulation: N-Form" << std::endl;
          break;

        case 2:
          std::cout << "Formulation: M-Form" << std::endl;
          break;
      }

      // Print optimization method.
      using OptimizationType = typename std::decay<
          typename std::tuple_element<5, decltype(tuple)>::type>::type;
      std::cout << "Method:      "
                << invlib::OptimizerLog<OptimizationType>::name;
      std::cout << std::endl;

      std::cout << std::endl;
      std::cout << std::setw(5) << "Step" << std::setw(15) << "Total Cost";
      std::cout << std::setw(15) << "x-Cost" << std::setw(15) << "y-Cost";
      std::cout << std::setw(15) << "Conv. Crit.";
      std::cout << std::setw(15) << OptimizerLog<OptimizationType>::header();
      std::cout << std::endl << invlib::separator() << std::endl;
    }
  }

  /** Print step to command line.
   *
   * This function is called from invlib to log a new step to the command
   * line.
   */
  template <typename... Params>
  void step(const Params &... params) {
    if (verbosity_ >= 1) {
      std::tuple<const Params &...> tuple(params...);
      using OptimizationType = typename std::decay<
          typename std::tuple_element<5, decltype(tuple)>::type>::type;

      auto step_number = std::get<0>(tuple);
      std::cout << std::setw(5) << step_number;
      if (step_number == 0) {
        start_cost_ = std::get<1>(tuple);
      }
      std::cout << std::setw(15) << scaling_factor_ * std::get<1>(tuple);
      std::cout << std::setw(15) << scaling_factor_ * std::get<2>(tuple);
      std::cout << std::setw(15) << scaling_factor_ * std::get<3>(tuple);

      if (std::isnan(std::get<4>(tuple))) {
        std::cout << std::setw(15) << " ";
      } else {
        std::cout << std::setw(15) << std::get<4>(tuple);
      }
      std::cout << OptimizerLog<OptimizationType>::log(
          std::get<5>(tuple), gamma_history_, std::get<0>(tuple));
      std::cout << std::endl;
    }
  }

  /** Finalize log output.
   *
   * This function is called from within invlib to finalize the log
   * output.
   */
  template <typename... Params>
  void finalize(const Params &... params) {
    if (verbosity_ >= 1) {
      std::cout << invlib::separator() << std::endl;

      std::tuple<const Params &...> tuple(params...);
      std::cout << std::endl;

      std::cout << "Total number of steps:            ";
      std::cout << std::get<1>(tuple) << std::endl;
      std::cout << "Final scaled cost function value: ";
      std::cout << std::get<2>(tuple) * scaling_factor_ << std::endl;

      bool converged = std::get<0>(tuple);
      if (converged) {
        std::cout << "OEM computation converged." << std::endl;
      } else if (linear_) {
        std::cout << "Linear OEM computation finished." << std::endl;
      } else {
        std::cout << "OEM computation DID NOT converge!" << std::endl;
      }
    }

    finalized_ = true;
  }

  /** Print timing information to command line.*/
  template <typename... Params>
  void time(const Params &... params) {
    if (verbosity_ >= 1) {
      std::tuple<const Params &...> tuple(params...);
      std::cout << std::endl;
      std::cout << "Elapsed Time for Retrieval:                       ";
      std::cout << std::get<0>(tuple) << std::endl;
      std::cout << "Time in inversion_iterate Agenda (No Jacobian):   ";
      std::cout << std::get<1>(tuple) << std::endl;
      std::cout << "Time in inversion_iterate Agenda (With Jacobian): ";
      std::cout << std::get<2>(tuple) << std::endl;

      std::cout << std::endl;
      std::cout << invlib::center("----") << std::endl;
      std::cout << std::endl;
    }
  }

 private:
  /** Verbosity level of logger */
  int verbosity_;
  /** Reference to ARTS vector holding the LM gamma history*/
  Vector gamma_history_;
  /** Scaling factor for the cost.*/
  Numeric scaling_factor_ = 0.0;
  /** Start cost (not computed by invlib)*/
  Numeric start_cost_ = 0.0;
  /** Flag indicating whether forward model is linear.*/
  bool linear_ = false;
  /** Flag indicating whether output has been finalized.*/
  bool finalized_ = false;
};

////////////////////////////////////////////////////////////////////////////////
//  Exception Handling
////////////////////////////////////////////////////////////////////////////////

/** Handle exception encountered within invlib.
 *
 * During OEM iteration invlib executes the ARTS inversion_iterate_agenda
 * multiple times during which error can occur. This function converts
 * nested exceptions to a vector of strings suitable for printing.
 *
 * @tparam E The exception type which to handle 
 * @param[in] e The specific exception type to handle
 * @param The nesting level, should be 0 when called.
 */
template <typename E>
std::vector<std::string> handle_nested_exception(const E &e, int level = 0) {
  const std::exception *re;
  std::vector<std::string> errors{};

  re = dynamic_cast<const std::exception *>(&e);
  if (re) {
    std::string s{};

    // If invlib level, extend error description.
    if (level == 0) {
      s = "Run-time error in oem computation: ";
    }

    s += re->what();
    errors.push_back(s);
  }

  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &ne) {
    std::vector<std::string> sv(handle_nested_exception(ne, level + 1));
    errors.insert(errors.end(), sv.begin(), sv.end());
  } catch (...) {
  }
  return errors;
}

////////////////////////////////////////////////////////////////////////////////
// Forward model interface
////////////////////////////////////////////////////////////////////////////////

/** Interface to ARTS inversion_iterate_agenda
 *  
 *  This wrapper class implements the invlib-to-ARTS interface to the
 *  inversion_iterate_agendaExecute function, which implements the forward
 *  model used in the invlib iteration.
 */
class AgendaWrapper {
 public:
  /** Dimension of the measurement space.*/
  const unsigned int m = 0;
  /** Dimension of the state space.*/
  const unsigned int n = 0;

  /** Create inversion_iterate_agendaExecute wrapper.
   *
   * Initializes the wrapper object for the inversion_iterate_agendaExecute
   * method. The object forwards the evaluate() and evaluate_jacobian() calls
   * made by the iterative OEM methods to inversion_iterate_agendaExecute using
   * the arguments provided to the constructor.
   * 
   * \param[in] ws Pointer to the current ARTS workspace.
   * \param[in] measurment_space_dimension Dimension of the measurement space
   * \param[in] arts_jacobian Reference to the jacobian WSV of the workspace.
   * \param[in] arts_y Reference to the arts y WSV.
   * \param[in] inversion_iterate_agenda Pointer to the x argument of the agenda
   * execution function.
   */
  AgendaWrapper(Workspace *ws,
                unsigned int measurement_space_dimension,
                unsigned int state_space_dimension,
                ::Matrix &arts_jacobian,
                ::Vector &arts_y,
                const Agenda *inversion_iterate_agenda)
      : m(measurement_space_dimension),
        n(state_space_dimension),
        inversion_iterate_agenda_(inversion_iterate_agenda),
        iteration_counter_(0),
        jacobian_(arts_jacobian),
        reuse_jacobian_((arts_jacobian.nrows() != 0) &&
                        (arts_jacobian.ncols() != 0) && (arts_y.nelem() != 0)),
        ws_(ws),
        yi_(arts_y) {}

  /** Return most recently simulated measurement vector.
   *
   * @return The simulated observation vector.
   */
  ArtsVector get_measurement_vector() { return yi_; }

  AgendaWrapper(const AgendaWrapper &) = delete;
  AgendaWrapper(AgendaWrapper &&) = delete;
  AgendaWrapper &operator=(const AgendaWrapper &) = delete;
  AgendaWrapper &operator=(AgendaWrapper &&) = delete;

  /** Evaluate forward model and compute Jacobian.
   *
   * Forwards the call to evaluate_jacobian() and evaluate() that is made by
   * Gauss-Newton and Levenberg-Marquardt OEM methods using the variables pointed
   * to by the pointers provided to the constructor as arguments.

   * \param[out] y The measurement vector y = K(x) for the current state vector x
   * as computed by the forward model.
   * \param[out] J The Jacobian Ki=d/dx(K(x)) of the forward model.
   * \param[in] x The current state vector x.
   */
  MatrixReference Jacobian(const Vector &xi, Vector &yi) {
    if (!reuse_jacobian_) {
      inversion_iterate_agendaExecute(
          *ws_, yi_, jacobian_, xi, 1, 0, *inversion_iterate_agenda_);
      yi = yi_;
      iteration_counter_ += 1;
    } else {
      reuse_jacobian_ = false;
      yi = yi_;
    }
    return jacobian_;
  }

  /** Evaluate the ARTS forward model.
   *
   * Call the ARTS forward model defined by inversion_iterate_agenda
   * and return resulting observation vector.
   * 
   * @param[in] xi The current state vector of the OEM iteration.
   * @return The observation vector y contained in the yf WSV after
   *   executing the inversion_iterate_agenda.
   */
  Vector evaluate(const Vector &xi) {
    if (!reuse_jacobian_) {
      Matrix dummy;
      inversion_iterate_agendaExecute(*ws_,
                                      yi_,
                                      dummy,
                                      xi,
                                      0,
                                      iteration_counter_,
                                      *inversion_iterate_agenda_);
    } else {
      reuse_jacobian_ = false;
    }
    return yi_;
  }

 private:
  /** Pointer to the inversion_iterate_agenda of the workspace. */
  const Agenda *inversion_iterate_agenda_;
  unsigned int iteration_counter_;
  /** Reference to the jacobian WSV.*/
  MatrixReference jacobian_;
  /** Flag whether to reuse Jacobian from previous calculation. */
  bool reuse_jacobian_;
  /** Pointer to current ARTS workspace */
  Workspace *ws_;
  /** Cached simulation result. */
  Vector yi_;
};
}  // namespace oem


/** Clip Tensor4
 *
 * @param[in] The tensor to which to apply the clipping.
 * @param[in] The book index to which to apply the clipping.
 * @param[in] limit_low Lower limit below which to clip values.
 * @param[in] limit_high Upper limit below which to clip values.
 */
void Tensor4Clip(Tensor4& x,
                 const Index& iq,
                 const Numeric& limit_low,
                 const Numeric& limit_high) {
  // Sizes
  const Index nq = x.nbooks();

  if (iq < -1) throw runtime_error("Argument *iq* must be >= -1.");
  if (iq >= nq) {
    ostringstream os;
    os << "Argument *iq* is too high.\n"
       << "You have selected index: " << iq << "\n"
       << "but the number of quantities is only: " << nq << "\n"
       << "(Note that zero-based indexing is used)\n";
    throw runtime_error(os.str());
  }

  Index ifirst = 0, ilast = nq - 1;
  if (iq > -1) {
    ifirst = iq;
    ilast = iq;
  }

  if (!std::isinf(limit_low)) {
    for (Index i = ifirst; i <= ilast; i++) {
      for (Index p = 0; p < x.npages(); p++) {
        for (Index r = 0; r < x.nrows(); r++) {
          for (Index c = 0; c < x.ncols(); c++) {
            if (x(i, p, r, c) < limit_low) x(i, p, r, c) = limit_low;
          }
        }
      }
    }
  }

  if (!std::isinf(limit_high)) {
    for (Index i = ifirst; i <= ilast; i++) {
      for (Index p = 0; p < x.npages(); p++) {
        for (Index r = 0; r < x.nrows(); r++) {
          for (Index c = 0; c < x.ncols(); c++) {
            if (x(i, p, r, c) > limit_high) x(i, p, r, c) = limit_high;
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// OEM error checking
////////////////////////////////////////////////////////////////////////////////

/** Error checking for OEM method.
 *
 * @param[in, out] x Checked to have size consistent with xa or zero.
 * Set to xa if empty.
 * @param[in, out] yf Checked to have size consistent with y or zero.
 * Computed by execution inversion_iterate_agenda if necessary.
 * @param[in, out] jacobian Checked to be consistent with xa and covmat_se
 * or empty. Computed by executin inversion_iterate_agenda if necessary.
 * @param[in] inversion_iterate_agenda The inversion_iterate_agenda to execute
 * to compute yf.
 * @param[in] xa The a priori vector
 * @param[in] covmat_sx The state-space covariance matrix. Checked to be 
 * square and consistent with xa.
 * @param[in] y The observation vector to fit.
 * @param[in] covmat_se The observation error covariance matrix. Checked to
 * by square and consistent with y.
 * @param[in] jacobian_do The jacobian_do flag of the workspace. Checked to
 * be 1
 * @param[in] jacobian_quantities: The Jacobian quantities array checked to
 * be consistent with jacobian_indices
 * @param[in] method The method string. Checked to be a valid OEM method
 * string.
 * @param[in] x_norm Vector to use to normalize linear systems occurring
 * in the OEM minimization. Checked to be same size as x or empty.
 * @param[in] max_iter Maximum number of OEM iteration. Checked to be positive.
 * @param[in] stop_dx The convergence criterion for the OEM iteration. Checked
 * to be positive.
 * @param[in] lm_ga_settings Vector containint setting for the Levenberg-Marquardt
 * method. Checked to contain 6 elements that are all greater or equal zero.
 * @param clear_matrices Flag whether or not to clear matrices after OEM run.
 * Checked to be 1 or 0.
 * @param display_progress Whether or not to display iteration progress. Checked
 * to be 1 or 0.
 */
void OEM_checks(Workspace& ws,
                Vector& x,
                Vector& yf,
                Matrix& jacobian,
                const Agenda& inversion_iterate_agenda,
                const Vector& xa,
                const CovarianceMatrix& covmat_sx,
                const Vector& y,
                const CovarianceMatrix& covmat_se,
                const Index& jacobian_do,
                const ArrayOfRetrievalQuantity& jacobian_quantities,
                const String& method,
                const Vector& x_norm,
                const Index& max_iter,
                const Numeric& stop_dx,
                const Vector& lm_ga_settings,
                const Index& clear_matrices,
                const Index& display_progress) {
  const Index nq = jacobian_quantities.nelem();
  const Index n = xa.nelem();
  const Index m = y.nelem();

  if ((x.nelem() != n) && (x.nelem() != 0))
    throw runtime_error(
        "The length of *x* must be either the same as *xa* or 0.");
  if (covmat_sx.ncols() != covmat_sx.nrows())
    throw runtime_error("*covmat_sx* must be a square matrix.");
  if (covmat_sx.ncols() != n)
    throw runtime_error("Inconsistency in size between *x* and *covmat_sx*.");
  if ((yf.nelem() != m) && (yf.nelem() != 0))
    throw runtime_error(
        "The length of *yf* must be either the same as *y* or 0.");
  if (covmat_se.ncols() != covmat_se.nrows())
    throw runtime_error("*covmat_se* must be a square matrix.");
  if (covmat_se.ncols() != m)
    throw runtime_error("Inconsistency in size between *y* and *covmat_se*.");
  if (!jacobian_do)
    throw runtime_error(
        "Jacobian calculations must be turned on (but jacobian_do=0).");
  if ((jacobian.nrows() != m) && (!jacobian.empty()))
    throw runtime_error(
        "The number of rows of the jacobian must be either the number of elements in *y* or 0.");
  if ((jacobian.ncols() != n) && (!jacobian.empty()))
    throw runtime_error(
        "The number of cols of the jacobian must be either the number of elements in *xa* or 0.");

  ArrayOfArrayOfIndex jacobian_indices;
  bool any_affine;
  jac_ranges_indices(jacobian_indices, any_affine, jacobian_quantities);
  if (jacobian_indices.nelem() != nq)
    throw runtime_error(
        "Different number of elements in *jacobian_quantities* "
        "and *jacobian_indices*.");
  if (nq && jacobian_indices[nq - 1][1] + 1 != n)
    throw runtime_error(
        "Size of *covmat_sx* do not agree with Jacobian "
        "information (*jacobian_indices*).");

  // Check GINs
  if (!(method == "li" || method == "gn" || method == "li_m" ||
        method == "gn_m" || method == "ml" || method == "lm" ||
        method == "li_cg" || method == "gn_cg" || method == "li_cg_m" ||
        method == "gn_cg_m" || method == "lm_cg" || method == "ml_cg")) {
    throw runtime_error(
        "Valid options for *method* are \"nl\", \"gn\" and "
        "\"ml\" or \"lm\".");
  }

  if (!(x_norm.nelem() == 0 || x_norm.nelem() == n)) {
    throw runtime_error(
        "The vector *x_norm* must have length 0 or match "
        "*covmat_sx*.");
  }

  if (x_norm.nelem() > 0 && min(x_norm) <= 0) {
    throw runtime_error("All values in *x_norm* must be > 0.");
  }

  if (max_iter <= 0) {
    throw runtime_error("The argument *max_iter* must be > 0.");
  }

  if (stop_dx <= 0) {
    throw runtime_error("The argument *stop_dx* must be > 0.");
  }

  if ((method == "ml") || (method == "lm") || (method == "lm_cg") ||
      (method == "ml_cg")) {
    if (lm_ga_settings.nelem() != 6) {
      throw runtime_error(
          "When using \"ml\", *lm_ga_setings* must be a "
          "vector of length 6.");
    }
    if (min(lm_ga_settings) < 0) {
      throw runtime_error(
          "The vector *lm_ga_setings* can not contain any "
          "negative value.");
    }
  }

  if (clear_matrices < 0 || clear_matrices > 1)
    throw runtime_error("Valid options for *clear_matrices* are 0 and 1.");
  if (display_progress < 0 || display_progress > 1)
    throw runtime_error("Valid options for *display_progress* are 0 and 1.");

  // If necessary compute yf and jacobian.
  if (x.nelem() == 0) {
    x = xa;
    inversion_iterate_agendaExecute(
        ws, yf, jacobian, xa, 1, 0, inversion_iterate_agenda);
  }
  if ((yf.nelem() == 0) || (jacobian.empty())) {
    inversion_iterate_agendaExecute(
        ws, yf, jacobian, x, 1, 0, inversion_iterate_agenda);
  }
}

#endif  // _ARTS_OEM_H_
