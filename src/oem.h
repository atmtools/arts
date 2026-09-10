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

#include <algorithm>
#include <type_traits>

#include "invlib/algebra.h"
#include "invlib/algebra/precision_matrix.h"
#include "invlib/algebra/solvers.h"
#include "invlib/interfaces/arts_wrapper.h"
#include "invlib/map.h"
#include "invlib/optimization.h"
#include "jacobian.h"

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
using Identity         = invlib::MatrixIdentity<Matrix>;

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
template <typename ForwardModel> using OEM_STANDARD = invlib::
    MAP<ForwardModel, Matrix, CovarianceMatrix, CovarianceMatrix, Vector, Formulation::STANDARD, invlib::Rodgers531>;

/** OEM n form.
 *
 * Class template implementing the n-form of the OEM
 * optimization defined according to Chapter 4 in Rodgers (2000).
 *
 * In this formulation, each iteration requires the solution
 * of a system of linear equations of size n-times-n.
 */
template <typename ForwardModel> using OEM_NFORM =
    invlib::MAP<ForwardModel, Matrix, CovarianceMatrix, CovarianceMatrix, Vector, Formulation::NFORM>;

/** OEM m form.
 *
 * Class template implementing the m-form of the OEM
 * optimization defined according to Chapter 4 in Rodgers (2000).
 *
 * In this formulation, each iteration requires the solution
 * of a system of linear equations of size m-times-m. Use the state-step
 * criterion (Rodgers 5.30): the m-form residual is not a state-space gradient
 * and cannot be used with Rodgers 5.31, especially when m != n.
 */
template <typename ForwardModel> using OEM_MFORM = invlib::
    MAP<ForwardModel, Matrix, CovarianceMatrix, CovarianceMatrix, Vector, Formulation::MFORM, invlib::Rodgers530>;

////////////////////////////////////////////////////////////////////////////////
// Solvers
////////////////////////////////////////////////////////////////////////////////

// Matrix-free symmetric scaling by inverse measurement standard deviations.
template <typename MatrixType> struct NoiseScaledSystem {
  const MatrixType& matrix;
  const ::Vector& scales;

  template <typename V> typename V::ResultType operator*(const V& value) const {
    typename V::ResultType input = value;
    for (Index i = 0; i < input.rows(); ++i) input(i) /= scales[i];
    typename V::ResultType output = matrix * input;
    for (Index i = 0; i < output.rows(); ++i) output(i) /= scales[i];
    return output;
  }
};

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
template <typename TransformationMatrixType, typename SolverType = invlib::Standard> class NormalizingSolver
    : SolverType {
 public:
  template <typename... Params> NormalizingSolver(const TransformationMatrixType &trans, bool apply, Params... params)
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
  template <typename MatrixType, typename VectorType> auto solve(const MatrixType &A, const VectorType &v) ->
      typename VectorType::ResultType {
    // invlib's relative CG residual is undefined for an exactly zero RHS.
    // The solution is zero for the positive-definite OEM systems.
    bool zero_rhs = true;
    for (Index i = 0; i < v.rows(); ++i) zero_rhs = zero_rhs && v(i) == 0;
    if (zero_rhs) return v;

    if constexpr (std::is_same_v<SolverType, invlib::ConjugateGradient<>>) {
      if (not measurement_scales.empty()) {
        typename VectorType::ResultType rhs = v;
        for (Index i = 0; i < rhs.rows(); ++i) rhs(i) /= measurement_scales[i];
        auto result = SolverType::solve(NoiseScaledSystem<MatrixType>{A, measurement_scales}, rhs);
        for (Index i = 0; i < result.rows(); ++i) result(i) /= measurement_scales[i];
        return result;
      }
    }
    typename VectorType::ResultType w;
    if (apply_) {
      typename VectorType::ResultType vv = trans_ * v;
      auto                          &&ww = SolverType::solve(trans_ * A * trans_, vv);
      w                                  = trans_ * ww;
    } else {
      w = SolverType::solve(A, v);
    }
    return w;
  }

  void set_iteration_limit_warning(std::function<void()> warning)
    requires std::is_same_v<SolverType, invlib::ConjugateGradient<>> {
    SolverType::iteration_limit_warning = std::move(warning);
  }

  ::Vector measurement_scales;

 private:
  /** Whether or not to apply the transformation.*/
  const bool apply_ = false;
  /** The transformation matrix.*/
  const TransformationMatrixType &trans_;
};

// Request matrix-matrix assembly in MFORM, then factor the measurement-sized
// system. CG does not opt into this policy and remains matrix-free.
struct DirectMeasurementSolver {
  static constexpr bool dense_measurement_system = true;
  ::Vector measurement_scales;

  template <typename M, typename V> typename V::ResultType solve(const M& system, const V& rhs) {
    Matrix dense = system;
    const Index size = rhs.rows();
    typename V::ResultType scaled_rhs = rhs;
    if (not measurement_scales.empty()) {
      for (Index i = 0; i < size; ++i) {
        scaled_rhs(i) /= measurement_scales[i];
        for (Index j = 0; j < size; ++j)
          dense(i, j) = (dense(i, j) / measurement_scales[i]) / measurement_scales[j];
      }
    }
    typename V::ResultType result = invlib::Standard{}.solve(dense, scaled_rhs);
    if (not measurement_scales.empty())
      for (Index i = 0; i < size; ++i) result(i) /= measurement_scales[i];
    return result;
  }
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
template <typename T> struct OptimizerLog;

/** Log customization for LM method
 *
 * This essentially adds a line for the gamma parameter to
 * the output.
 */
template <typename RealType, typename DampingMatrix, typename Solver>
struct OptimizerLog<invlib::LevenbergMarquardt<RealType, DampingMatrix, Solver>> {
  /** Method name */
  static constexpr auto name = "Levenberg-Marquardt";

  /** Name to append to header line. */
  static std::string header() {
    std::string out = "Gamma Factor";
    return out;
  }

  /** Returns the string to append to the log of a single step. */
  static std::string log(const invlib::LevenbergMarquardt<RealType, DampingMatrix, Solver> &g,
                         ::Vector                                                          &gamma_history_,
                         size_t                                                             i) {
    std::string lambda = std::to_string(g.get_lambda());
    std::string out(15 - std::min<size_t>(lambda.size(), 15), ' ');
    out               += lambda;
    gamma_history_[i]  = g.get_lambda();
    return out;
  }
};

/** Log customization for GN method. */
template <typename RealType, typename Solver> struct OptimizerLog<invlib::GaussNewton<RealType, Solver>> {
  /** Method name */
  static constexpr auto name = "Gauss-Newton";

  /** Name to append to header line. */
  static std::string header() { return ""; }

  static std::string log(const invlib::GaussNewton<RealType, Solver> &, ::Vector &, size_t) { return ""; }
};

/** OEM log output
 *
 * This class takes care of formatting the OEM iteration information
 * and displaying in the command line.
 *
 * @tparam The invlib log type defining which type of logging to perform.
 */
template <invlib::LogType type> class ArtsLog {
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
      std::cout << invlib::separator() << '\n' << '\n';
      std::cout << "Error during OEM computation." << '\n';
      std::cout << '\n';
      std::cout << invlib::center("----") << '\n';
      std::cout << '\n';
    }
  }

  /** Initialize log output.
   *
   * This function is called from within invlib to initialize the log
   * output. Prints general information on the OEM settings.
   */
  template <typename... Params> void init(Params &...params) {
    if (verbosity_ >= 1) {
      std::tuple<Params &...> tuple(params...);

      auto &y         = std::get<4>(tuple);
      scaling_factor_ = 1.0 / static_cast<Numeric>(y.size());
      std::cout << '\n';
      std::cout << invlib::center("MAP Computation") << '\n';

      // Print formulation.
      int formulation = static_cast<int>(std::get<6>(tuple));
      switch (formulation) {
        case 0: std::cout << "Formulation: Standard" << '\n'; break;
        case 1: std::cout << "Formulation: N-Form" << '\n'; break;

        case 2: std::cout << "Formulation: M-Form" << '\n'; break;
      }

      // Print optimization method.
      using OptimizationType = typename std::decay<typename std::tuple_element<5, decltype(tuple)>::type>::type;
      std::cout << "Method:      " << invlib::OptimizerLog<OptimizationType>::name;
      std::cout << '\n';

      std::cout << '\n';
      std::cout << std::setw(5) << "Step" << std::setw(15) << "Total Cost";
      std::cout << std::setw(15) << "x-Cost" << std::setw(15) << "y-Cost";
      std::cout << std::setw(15) << "Conv. Crit.";
      std::cout << std::setw(15) << OptimizerLog<OptimizationType>::header();
      std::cout << '\n' << invlib::separator() << '\n';
    }
  }

  /** Print step to command line.
   *
   * This function is called from invlib to log a new step to the command
   * line.
   */
  template <typename... Params> void step(const Params &...params) {
    std::tuple<const Params &...> tuple(params...);
    using OptimizationType = std::remove_cvref_t<decltype(std::get<5>(tuple))>;
    // History is an output of OEM, independent of console verbosity.
    const auto optimizer_log =
        OptimizerLog<OptimizationType>::log(std::get<5>(tuple), gamma_history_, std::get<0>(tuple));
    if (verbosity_ >= 1) {
      auto step_number = std::get<0>(tuple);
      std::cout << std::setw(5) << step_number;
      if (step_number == 0) { start_cost_ = std::get<1>(tuple); }
      std::cout << std::setw(15) << scaling_factor_ * std::get<1>(tuple);
      std::cout << std::setw(15) << scaling_factor_ * std::get<2>(tuple);
      std::cout << std::setw(15) << scaling_factor_ * std::get<3>(tuple);

      if (std::isnan(std::get<4>(tuple))) {
        std::cout << std::setw(15) << " ";
      } else {
        std::cout << std::setw(15) << std::get<4>(tuple);
      }
      std::cout << optimizer_log;
      std::cout << '\n';
    }
  }

  /** Finalize log output.
   *
   * This function is called from within invlib to finalize the log
   * output.
   */
  template <typename... Params> void finalize(const Params &...params) {
    if (verbosity_ >= 1) {
      std::cout << invlib::separator() << '\n';

      std::tuple<const Params &...> tuple(params...);
      std::cout << '\n';

      std::cout << "Total number of steps:            ";
      std::cout << std::get<1>(tuple) << '\n';
      std::cout << "Final scaled cost function value: ";
      std::cout << std::get<2>(tuple) * scaling_factor_ << '\n';

      bool converged = std::get<0>(tuple);
      if (converged) {
        std::cout << "OEM computation converged." << '\n';
      } else if (linear_) {
        std::cout << "Linear OEM computation finished." << '\n';
      } else {
        std::cout << "OEM computation DID NOT converge!" << '\n';
      }
    }

    finalized_ = true;
  }

  /** Print timing information to command line.*/
  template <typename... Params> void time(const Params &...params) {
    if (verbosity_ >= 1) {
      std::tuple<const Params &...> tuple(params...);
      std::cout << '\n';
      std::cout << "Elapsed Time for Retrieval:                       ";
      std::cout << std::get<0>(tuple) << '\n';
      std::cout << "Time in inversion_iterate Agenda (No Jacobian):   ";
      std::cout << std::get<1>(tuple) << '\n';
      std::cout << "Time in inversion_iterate Agenda (With Jacobian): ";
      std::cout << std::get<2>(tuple) << '\n';

      std::cout << '\n';
      std::cout << invlib::center("----") << '\n';
      std::cout << '\n';
    }
  }

 private:
  /** Verbosity level of logger */
  int verbosity_;
  /** Reference to ARTS vector holding the LM gamma history*/
  ::Vector &gamma_history_;
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
template <typename E> std::vector<std::string> handle_nested_exception(const E &e, int level = 0) {
  const std::exception    *re;
  std::vector<std::string> errors{};

  re = dynamic_cast<const std::exception *>(&e);
  if (re) {
    std::string s{};

    // If invlib level, extend error description.
    if (level == 0) { s = "Run-time error in oem computation: "; }

    s += re->what();
    errors.push_back(s);
  }

  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &ne) {
    std::vector<std::string> sv(handle_nested_exception(ne, level + 1));
    errors.insert(errors.end(), sv.begin(), sv.end());
  } catch (...) {}
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
  AgendaWrapper(const Workspace *const ws,
                unsigned int           measurement_space_dimension,
                unsigned int           state_space_dimension,
                ::Matrix              &arts_jacobian,
                ::Vector              &arts_y,
                const ::Vector        &initial_state,
                AtmField              *atm_field,
                AbsorptionBands       *abs_bands,
                ArrayOfSensorObsel    *measurement_sensor,
                SurfaceField          *surf_field,
                SubsurfaceField       *subsurf_field,
                const JacobianTargets *jac_targets,
                const Agenda          *inversion_iterate_agenda)
      : m(measurement_space_dimension),
        n(state_space_dimension),
        inversion_iterate_agenda_(inversion_iterate_agenda),
        jacs(jac_targets),
        atm(atm_field),
        absdata(abs_bands),
        sensor(measurement_sensor),
        surf(surf_field),
        subsurf(subsurf_field),
        jacobian_(arts_jacobian),
        ws_(ws),
        yi_(arts_y),
        measurement_state_(initial_state),
        jacobian_state_(initial_state),
        measurement_valid_(arts_y.size() == m),
        jacobian_valid_(measurement_valid_ && arts_jacobian.nrows() == m && arts_jacobian.ncols() == n) {}

  /** Return most recently simulated measurement vector.
   *
   * @return The simulated observation vector.
   */
  const ::Vector &get_measurement_vec() const { return yi_; }

  AgendaWrapper(const AgendaWrapper &)            = delete;
  AgendaWrapper(AgendaWrapper &&)                 = delete;
  AgendaWrapper &operator=(const AgendaWrapper &) = delete;
  AgendaWrapper &operator=(AgendaWrapper &&)      = delete;

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
    ensure_jacobian(xi);
    // Assign directly to the ARTS storage, avoiding an intermediate invlib vector.
    static_cast<::Vector &>(yi) = yi_;
    return jacobian_;
  }

  /** Ensure that the Jacobian and physical model describe xi.
   * A value-only trial preserves the previous Jacobian, but may change the
   * atmosphere or other inouts. Restore those before returning a cached matrix.
   */
  void ensure_jacobian(const Vector &xi) {
    if (!jacobian_valid_ || !same_state(jacobian_state_, xi)) {
      execute(xi, true);
    } else if (!measurement_valid_ || !same_state(measurement_state_, xi)) {
      execute(xi, false);
    }
  }

  /** Evaluate the forward model, reusing only the most recent successful state.
   * Results from before a failed or rejected trial cannot restore physical
   * inouts; those require another agenda execution at the accepted state.
   */
  Vector evaluate(const Vector &xi) {
    if (!measurement_valid_ || !same_state(measurement_state_, xi)) execute(xi, false);
    return yi_;
  }

 private:
  static bool same_state(const ::Vector &cached, const Vector &state) {
    return cached.size() == state.size() && std::equal(cached.elem_begin(), cached.elem_end(), state.elem_begin());
  }

  void execute(const Vector &xi, bool with_jacobian) {
    // The agenda can partially modify its outputs before throwing. Invalidate
    // first, and publish a new state tag only after all output checks succeed.
    measurement_valid_ = false;
    if (with_jacobian) jacobian_valid_ = false;
    ::Matrix               dummy;
    auto                  &jacobian = with_jacobian ? static_cast<::Matrix &>(jacobian_) : dummy;
    const JacobianTargets  no_targets{};
    const JacobianTargets &derivative_targets = with_jacobian ? *jacs : no_targets;
    inversion_iterate_agendaExecute(*ws_,
                                    *atm,
                                    *absdata,
                                    *sensor,
                                    *surf,
                                    *subsurf,
                                    yi_,
                                    jacobian,
                                    *jacs,
                                    derivative_targets,
                                    xi,
                                    *inversion_iterate_agenda_);
    ARTS_USER_ERROR_IF(yi_.size() != m, "inversion_iterate_agenda must return {} measurements; got {}.", m, yi_.size())
    if (with_jacobian) {
      ARTS_USER_ERROR_IF(jacobian.nrows() != m || jacobian.ncols() != n,
                         "inversion_iterate_agenda must return an {} by {} Jacobian.",
                         m,
                         n)
      jacobian_state_ = static_cast<const ::Vector &>(xi);
      jacobian_valid_ = true;
    }
    measurement_state_ = static_cast<const ::Vector &>(xi);
    measurement_valid_ = true;
  }

  const Agenda          *inversion_iterate_agenda_;
  const JacobianTargets *jacs;
  AtmField              *atm;
  AbsorptionBands       *absdata;
  ArrayOfSensorObsel    *sensor;
  SurfaceField          *surf;
  SubsurfaceField       *subsurf;
  MatrixReference        jacobian_;
  const Workspace *const ws_;
  // Borrow the workspace output rather than copying it into and out of the adapter.
  ::Vector &yi_;
  // O(n) state tags; the potentially much larger Jacobian stays in its original storage.
  ::Vector measurement_state_, jacobian_state_;
  bool     measurement_valid_, jacobian_valid_;
};
}  // namespace oem

#endif  // _ARTS_OEM_H_
