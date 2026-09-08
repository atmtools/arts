#include <callback.h>
#include <oem_settings.h>
#include <workspace.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string_view>

#include "invlib/algebra.h"
#include "invlib/algebra/solvers.h"
#include "invlib/archetypes/matrix_archetype.h"
#include "invlib/optimization/levenberg_marquardt.h"
#include "invlib/optimization/minimize.h"

namespace {

constexpr std::array<std::string_view, 10> methods{
    "li", "li_cg", "li_cg_m", "gn", "gn_cg", "gn_cg_m", "lm", "ml", "lm_cg", "ml_cg"};

void require(bool condition, std::string_view message) {
  if (not condition) throw std::runtime_error(std::string(message));
}

void close(Numeric actual, Numeric expected, Numeric tolerance, std::string_view name) {
  if (not std::isfinite(actual) or std::abs(actual - expected) > tolerance) {
    throw std::runtime_error(std::format("{}: got {}, expected {} (tolerance {})", name, actual, expected, tolerance));
  }
}

Matrix matrix(Index rows, Index cols, std::initializer_list<Numeric> values) {
  Matrix result(rows, cols);
  auto   value = values.begin();
  for (Index i = 0; i < rows; ++i) {
    for (Index j = 0; j < cols; ++j) result[i, j] = *value++;
  }
  return result;
}

CovarianceMatrix covariance(Matrix value) {
  const Index      n = value.nrows();
  CovarianceMatrix result;
  result.add_correlation(Block(Range(0, n), Range(0, n), {0, 0}, std::make_shared<Matrix>(std::move(value))));
  return result;
}

bool is_lm(std::string_view method) { return method.starts_with("lm") or method.starts_with("ml"); }

// This fixture uses the public workspace interface, including agenda validation,
// covariance inversion, optimizer dispatch, diagnostics, and gain calculation.
// No atmospheric data files or Python/Matlab installation are required.
struct Retrieval {
  Workspace          ws{WorkspaceInitialization::Empty};
  Vector             x, yf;
  Matrix             jac, gain;
  AtmField           atm;
  AbsorptionBands    bands;
  ArrayOfSensorObsel sensor;
  SurfaceField       surf;
  SubsurfaceField    subsurf;
  Vector             diagnostics, history;
  ArrayOfString      errors;
  JacobianTargets    targets;
  Vector             xa{0.5, -0.25};
  Vector             y{2.0, -1.0, 1.5};
  CovarianceMatrix   sa = covariance(matrix(2, 2, {4, 1, 1, 2}));
  CovarianceMatrix   se = covariance(matrix(3, 3, {1, 0.2, 0, 0.2, 2, 0.3, 0, 0.3, 0.5}));
  Agenda             agenda{"inversion_iterate_agenda"};
  Vector             normalization;
  Vector             settings{10, 3, 2, 1e8, 0.1, 0};
  Index              max_iter = 40;
  // Accuracy is also checked against independent state and cost oracles.
  Numeric stop_dx          = 1e-9;
  Numeric max_start_cost   = std::numeric_limits<Numeric>::infinity();
  Index   clear_matrices   = 0;
  Index   display_progress = 0;
  Index   calls            = 0;
  Vector  first_state;

  std::function<void(const Vector&, Vector&, Matrix&, bool)> forward;

  Retrieval() {
    forward = [](const Vector& state, Vector& fit, Matrix& jacobian, bool do_jac) {
      fit = Vector{state[0] + 2 * state[1] + 0.25, 2 * state[0] - state[1] - 0.5, state[0] + state[1] + 1};
      if (do_jac)
        jacobian = matrix(3, 2, {1, 2, 2, -1, 1, 1});
      else
        jacobian.resize(0, 0);
    };
    set_target_size(2);

    CallbackOperator callback;
    callback.inputs     = {"model_state_vec", "do_jac"};
    callback.outputs    = {"measurement_vec_fit", "measurement_jac"};
    callback.callback.f = [this](Workspace& local) {
      const auto& state = local.get<Vector>("model_state_vec");
      if (calls++ == 0) first_state = state;
      forward(state,
              local.get<Vector>("measurement_vec_fit"),
              local.get<Matrix>("measurement_jac"),
              local.get<Index>("do_jac") != 0);
    };
    agenda.add(Method("analytical_forward_model", Wsv{callback}));
    agenda.finalize(true);
  }

  void set_target_size(Size n) {
    targets.atm       = {Jacobian::AtmTarget{.type = AtmKey::t, .target_pos = 0, .x_start = 0, .x_size = n}};
    targets.finalized = true;
  }

  void run(std::string_view method) {
    OEM(ws,
        x,
        yf,
        jac,
        atm,
        bands,
        sensor,
        surf,
        subsurf,
        gain,
        diagnostics,
        history,
        errors,
        targets,
        xa,
        sa,
        y,
        se,
        agenda,
        String{method},
        max_start_cost,
        normalization,
        max_iter,
        stop_dx,
        settings,
        clear_matrices,
        display_progress);
  }
};

void check_history(const Retrieval& r, std::string_view method) {
  if (not is_lm(method)) {
    require(r.history.empty(), "Non-LM method returned gamma history");
    return;
  }
  require(r.history.size() == static_cast<Size>(r.max_iter + 1), "Gamma history has wrong size");
  close(r.history[0], r.settings[0], 0, "Initial gamma");
  const auto iterations = static_cast<Index>(r.diagnostics[4]);
  for (Index i = 0; i <= iterations; ++i) {
    require(std::isfinite(r.history[i]) and r.history[i] >= 0, "Missing silent-mode gamma history");
  }
  for (Index i = iterations + 1; i <= r.max_iter; ++i) {
    require(std::isnan(r.history[i]), "Unused gamma history must remain NaN");
  }
}

void check_affine_solution(const Retrieval& r, std::string_view method) {
  require(r.errors.empty(), "Successful retrieval returned errors");
  require(r.diagnostics.size() == 5, "Diagnostics must contain five values");
  require(r.x.size() == 2 and r.yf.size() == 3, "Wrong retrieval output dimensions");

  // Exact rational oracle for K = [[1,2],[2,-1],[1,1]], offset =
  // [1/4,-1/2,1], xa = [1/2,-1/4], y = [2,-1,3/2], and the
  // correlated covariances above. Obtained by rational Gaussian elimination:
  // G = (Sa^-1 + K^T Se^-1 K)^-1 K^T Se^-1,
  // x = xa + G (y - K xa - offset). No OEM or ARTS solver supplies the oracle.
  close(r.x[0], 3112.0 / 18575, 2e-7, "MAP state[0]");
  close(r.x[1], 21727.0 / 37150, 2e-7, "MAP state[1]");
  close(r.yf[0], 117931.0 / 74300, 4e-7, "Fitted measurement[0]");
  close(r.yf[1], -13927.0 / 18575, 4e-7, "Fitted measurement[1]");
  close(r.yf[2], 65101.0 / 37150, 4e-7, "Fitted measurement[2]");
  close(r.diagnostics[2], 245309.0 / 891600, 1e-10, "Normalized MAP cost");
  close(r.diagnostics[3], 581955097.0 / 5520490000, 2e-7, "Normalized measurement cost");
  if (method.starts_with("li")) {
    require(r.diagnostics[0] == 0 or r.diagnostics[0] == 1, "Linear retrieval failed");
    close(r.diagnostics[4], 1, 0, "Linear iteration count");
  } else {
    close(r.diagnostics[0], 0, 0, "Convergence status");
    require(r.diagnostics[4] > 0 and r.diagnostics[4] <= static_cast<Numeric>(r.max_iter), "Invalid iteration count");
  }

  if (r.clear_matrices) {
    require(r.jac.empty() and r.gain.empty(), "clear_matrices did not clear Jacobian and gain");
  } else {
    require(r.jac.nrows() == 3 and r.jac.ncols() == 2, "Wrong Jacobian dimensions");
    require(r.gain.nrows() == 2 and r.gain.ncols() == 3, "Wrong gain dimensions");
    const Matrix expected_gain =
        matrix(2, 3, {1276.0 / 18575, 1092.0 / 3715, 4586.0 / 18575, 4323.0 / 18575, -784.0 / 3715, 4328.0 / 18575});
    for (Index i = 0; i < 2; ++i) {
      for (Index j = 0; j < 3; ++j) close(r.gain[i, j], expected_gain[i, j], 1e-11, "Gain");
    }
  }
  check_history(r, method);
}

void test_affine(std::string_view method) {
  for (const bool cached : {false, true}) {
    for (const bool normalized : {false, true}) {
      if (normalized and method.ends_with("_m")) continue;
      Retrieval r;
      if (normalized) r.normalization = Vector{0.25, 4};
      if (cached) {
        // A user-provided starting state is allowed to differ from the prior.
        r.x = Vector{1.5, -0.75};
        r.forward(r.x, r.yf, r.jac, true);
      }
      r.run(method);
      check_affine_solution(r, method);
      if (not cached) {
        close(r.diagnostics[1], 2863.0 / 1424, 1e-12, "Normalized initial cost");
        close(r.first_state[0], r.xa[0], 0, "Initial agenda state[0]");
        close(r.first_state[1], r.xa[1], 0, "Initial agenda state[1]");
      }
    }
  }

  Retrieval fresh_start;
  fresh_start.x = Vector{1.5, -0.75};
  fresh_start.run(method);
  check_affine_solution(fresh_start, method);
  close(fresh_start.first_state[0], 1.5, 0, "Uncached initial state[0]");
  close(fresh_start.first_state[1], -0.75, 0, "Uncached initial state[1]");

  Retrieval cleared;
  cleared.clear_matrices = 1;
  cleared.run(method);
  check_affine_solution(cleared, method);
}

void test_exact_start(std::string_view method) {
  Retrieval r;
  // Both residual and prior departure are exactly zero. In particular, CG
  // must return a zero step without attempting a 0/0 residual normalization.
  r.forward(r.xa, r.y, r.jac, true);
  r.jac.resize(0, 0);
  r.run(method);
  require(r.errors.empty(), "Exact initial solution returned errors");
  close(r.x[0], r.xa[0], 0, "Exact initial state[0]");
  close(r.x[1], r.xa[1], 0, "Exact initial state[1]");
  close(r.diagnostics[0], 0, 0, "Exact initial solution status");
  close(r.diagnostics[2], 0, 0, "Exact initial solution cost");
}

void test_disabled_start_cost(std::string_view method) {
  Retrieval r;
  r.max_start_cost = -1;
  r.run(method);
  check_affine_solution(r, method);
  if (is_lm(method))
    close(r.diagnostics[1], 2863.0 / 1424, 1e-12, "LM initial cost with disabled limit");
  else
    require(std::isnan(r.diagnostics[1]), "Unrequested initial cost should be NaN");
}

void test_runtime_failure(std::string_view method) {
  Retrieval  r;
  const auto working_forward = r.forward;
  r.forward                  = [&](const Vector& state, Vector& fit, Matrix& jacobian, bool do_jac) {
    if (r.calls > 1) throw std::runtime_error("deliberate regression forward-model failure");
    working_forward(state, fit, jacobian, do_jac);
  };
  r.settings[3] = 20;
  r.gain        = matrix(1, 1, {999});
  r.errors      = {"stale error"};
  r.run(method);
  close(r.diagnostics[0], 9, 0, "Forward-model failure status");
  close(r.diagnostics[2], 2863.0 / 1424, 1e-12, "Normalized failure cost");
  close(r.diagnostics[3], 2863.0 / 1424, 1e-12, "Normalized failure measurement cost");
  close(r.diagnostics[4], 0, 0, "Failure iteration count");
  require(stdr::all_of(r.x, [](Numeric value) { return std::isnan(value); }), "Failed state must be NaN");
  require(r.gain.empty(), "Failed retrieval returned stale gain");
  require(stdr::none_of(r.errors, [](const String& value) { return value == "stale error"; }),
          "Failed retrieval retained previous errors");
  require(stdr::any_of(r.errors,
                       [](const String& value) {
                         return std::string_view(value).contains("deliberate regression forward-model failure");
                       }),
          "Forward-model exception was not reported");

  r.forward = working_forward;
  r.x.resize(0);
  r.run(method);
  check_affine_solution(r, method);
}

void quadratic_model(Retrieval& r) {
  r.xa = Vector{1};
  r.y  = Vector{4};
  r.sa = covariance(matrix(1, 1, {4}));
  r.se = covariance(matrix(1, 1, {0.25}));
  r.set_target_size(1);
  r.forward = [](const Vector& state, Vector& fit, Matrix& jacobian, bool do_jac) {
    fit = Vector{state[0] * state[0]};
    if (do_jac)
      jacobian = matrix(1, 1, {2 * state[0]});
    else
      jacobian.resize(0, 0);
  };
}

void test_nonlinear(std::string_view method) {
  if (method.starts_with("li")) return;
  Retrieval r;
  quadratic_model(r);
  r.run(method);
  require(r.errors.empty(), "Nonlinear retrieval returned errors");
  close(r.diagnostics[0], 0, 0, "Nonlinear convergence status");

  // Positive MAP solution for F(x)=x^2 solves 32*x^3 - 127*x - 1=0.
  // Bisection on [1,2] supplies an independent, deterministic scalar oracle.
  Numeric low = 1, high = 2;
  for (Index i = 0; i < 60; ++i) {
    const Numeric mid = (low + high) / 2;
    if (32 * mid * mid * mid - 127 * mid - 1 < 0)
      low = mid;
    else
      high = mid;
  }
  const Numeric expected = (low + high) / 2;
  close(r.x[0], expected, 2e-7, "Nonlinear MAP state");
  close(r.yf[0], r.x[0] * r.x[0], 1e-12, "Nonlinear fitted measurement");
  close(r.jac[0, 0], 2 * r.x[0], 1e-12, "Jacobian at returned state");
  close(r.gain[0, 0], 8 * r.x[0] / (0.25 + 16 * r.x[0] * r.x[0]), 1e-12, "Nonlinear gain");
  const Numeric cost_y = 4 * std::pow(4 - r.yf[0], 2);
  close(r.diagnostics[2], 0.25 * std::pow(r.x[0] - 1, 2) + cost_y, 1e-12, "Nonlinear MAP cost");
  close(r.diagnostics[3], cost_y, 1e-12, "Nonlinear measurement cost");
  check_history(r, method);
}

void test_underdetermined(std::string_view method) {
  Retrieval r;
  r.y       = Vector{2};
  r.sa      = covariance(matrix(2, 2, {4, 0, 0, 9}));
  r.se      = covariance(matrix(1, 1, {1}));
  r.forward = [](const Vector& state, Vector& fit, Matrix& jacobian, bool do_jac) {
    fit = Vector{state[0] + 2 * state[1]};
    if (do_jac)
      jacobian = matrix(1, 2, {1, 2});
    else
      jacobian.resize(0, 0);
  };
  r.run(method);
  require(r.errors.empty(), "Underdetermined retrieval returned errors");
  close(r.x[0], 57.0 / 82, 2e-7, "Underdetermined state[0]");
  close(r.x[1], 103.0 / 164, 2e-7, "Underdetermined state[1]");
  close(r.yf[0], 80.0 / 41, 4e-7, "Underdetermined fit");
  close(r.gain[0, 0], 4.0 / 41, 1e-12, "Underdetermined gain[0]");
  close(r.gain[1, 0], 18.0 / 41, 1e-12, "Underdetermined gain[1]");
  close(r.diagnostics[2], 4.0 / 41, 1e-10, "Underdetermined cost");
  close(r.diagnostics[3], 4.0 / 1681, 2e-7, "Underdetermined measurement cost");
}

void test_lm_settings() {
  for (const auto method : {"lm", "ml", "lm_cg", "ml_cg"}) {
    Retrieval r;
    r.max_iter = 1;
    r.stop_dx  = 1e3;
    r.settings = OEMLMSettings{.initial_damping           = 12,
                               .decrease_factor           = 3,
                               .increase_factor           = 2,
                               .maximum_damping           = 1e6,
                               .damping_threshold         = 0.01,
                               .convergence_damping_limit = 10}
                     .as_vector();
    r.run(method);
    require(r.errors.empty(), "Damped step returned errors");
    // One step with damping 12*diag(Sa^-1). The non-diagonal prior
    // distinguishes diag(Sa^-1) from inverse(diag(Sa)).
    close(r.x[0], 31732.0 / 69551, 1e-11, "Damped state[0]");
    close(r.x[1], 24127.0 / 139102, 1e-11, "Damped state[1]");
    close(r.history[0], 12, 0, "Configured initial gamma");
    close(r.history[1], 4, 0, "Configured gamma decrease");
    // The supplied tolerance and convergence gamma limit allow convergence
    // now. Ignoring either setting requires further iterations.
    close(r.diagnostics[0], 0, 0, "Configured convergence settings");
    close(r.diagnostics[4], 1, 0, "Damped step iteration count");

    Retrieval rejected_trials;
    quadratic_model(rejected_trials);
    rejected_trials.x        = Vector{0.1};
    rejected_trials.max_iter = 1;
    rejected_trials.settings = Vector{0, 3, 2, 100, 0.1, 0};
    rejected_trials.run(method);
    require(rejected_trials.errors.empty(), "LM failed to recover from rejected trials");
    // Starting at x=0.1 gives g=-3.417, H=0.41, D=0.25. The cost
    // rejects gamma=0,0.1,0.2,0.4,0.8,1.6,3.2. At gamma=6.4 the
    // accepted step is 3.417/(0.41+6.4*0.25)=1.7.
    close(rejected_trials.x[0], 1.8, 1e-12, "State after rejected LM trials");
    close(rejected_trials.history[1], 6.4, 1e-12, "Configured gamma threshold and increase");
    close(rejected_trials.diagnostics[0], 1, 0, "LM iteration limit");

    Retrieval gamma_limit;
    quadratic_model(gamma_limit);
    gamma_limit.x        = Vector{0.1};
    gamma_limit.settings = Vector{0, 3, 2, 1, 0.1, 0};
    gamma_limit.run(method);
    require(gamma_limit.errors.empty(), "Gamma exhaustion is not an agenda exception");
    close(gamma_limit.diagnostics[0], 2, 0, "Configured gamma maximum");
    close(gamma_limit.x[0], 0.1, 1e-12, "Rejected LM step must not change state");
    close(gamma_limit.yf[0], 0.01, 1e-12, "Fit after gamma exhaustion");
  }
}

void test_lm_outcomes() {
  for (const auto method : {"lm", "ml", "lm_cg", "ml_cg"}) {
    Retrieval over_damped;
    over_damped.settings = OEMLMSettings{.initial_damping = 1e20, .maximum_damping = 1e20}.as_vector();
    over_damped.run(method);
    // Adding one to this maximum rounds back to the maximum. A numeric
    // sentinel used to turn this rejected, unchanged step into convergence.
    require(over_damped.errors.empty(), "Damping exhaustion is not an agenda exception");
    close(over_damped.diagnostics[0], 2, 0, "Huge damping must report exhaustion");
    close(over_damped.diagnostics[4], 1, 0, "Huge damping iteration count");
    close(over_damped.x[0], over_damped.xa[0], 0, "Huge damping state[0]");
    close(over_damped.x[1], over_damped.xa[1], 0, "Huge damping state[1]");
    close(over_damped.diagnostics[2], over_damped.diagnostics[1], 0, "Huge damping unchanged cost");
    close(over_damped.history[1], 1e20, 0, "Huge damping history retains actual damping");
    // Initial Jacobian, initial cost, one trial, restoration of the accepted
    // state's fit, and the public interface's final Jacobian refresh.
    require(over_damped.calls <= 5,
            std::format("Huge damping needlessly retried an exhausted step ({} agenda calls)", over_damped.calls));

    for (const Numeric tolerance : {1e-12, 1e-16, 1e-20}) {
      Retrieval tight;
      tight.stop_dx = tolerance;
      tight.run(method);
      check_affine_solution(tight, method);
      // Once the predicted reduction is beneath the cost's floating-point
      // resolution, comparing two rounded cost values must not exhaust damping.
      require(tight.calls < 30, "Stationary affine retrieval exhausted damping trials");
    }

    for (const Numeric damping : {10., 1e20}) {
      Retrieval stationary;
      stationary.xa = Vector{0};
      stationary.x  = Vector{1};
      stationary.y  = Vector{2};
      stationary.sa = covariance(matrix(1, 1, {1}));
      stationary.se = covariance(matrix(1, 1, {1}));
      stationary.set_target_size(1);
      stationary.settings = OEMLMSettings{.initial_damping = damping, .maximum_damping = damping}.as_vector();
      stationary.stop_dx  = 1e-20;
      stationary.forward  = [](const Vector& state, Vector& fit, Matrix& jacobian, bool do_jac) {
        fit = state;
        if (do_jac)
          jacobian = matrix(1, 1, {1});
        else
          jacobian.resize(0, 0);
      };
      stationary.run(method);
      // J=x^2+(2-x)^2 has an exactly representable stationary point at x=1,
      // although neither the measurement residual nor the prior departure is zero.
      require(stationary.errors.empty(), "Stationary nonzero-cost retrieval returned errors");
      close(stationary.diagnostics[0], 0, 0, "Stationary nonzero-cost convergence status");
      close(stationary.x[0], 1, 0, "Stationary nonzero-cost state");
      close(stationary.yf[0], 1, 0, "Stationary nonzero-cost fit");
      close(stationary.diagnostics[2], 2, 0, "Stationary nonzero cost");
      close(stationary.gain[0, 0], 0.5, 0, "Stationary nonzero-cost gain");
      require(stationary.calls <= 4, "Stationary point required repeated damping trials");
    }
  }
}

void test_named_settings() {
  const Vector legacy{12, 3, 2, 1e6, 0.01, 10};
  auto         named = OEMLMSettings::from_vector(legacy);
  close(named.initial_damping, 12, 0, "Named initial damping");
  close(named.decrease_factor, 3, 0, "Named decrease divisor");
  close(named.increase_factor, 2, 0, "Named increase multiplier");
  close(named.maximum_damping, 1e6, 0, "Named maximum damping");
  close(named.damping_threshold, 0.01, 0, "Named damping threshold");
  close(named.convergence_damping_limit, 10, 0, "Named convergence damping limit");
  const auto roundtrip = named.as_vector();
  for (Index i = 0; i < 6; ++i) close(roundtrip[i], legacy[i], 0, "Named settings round trip");

  const auto   defaults = OEMLMSettings{}.as_vector();
  const Vector expected_defaults{10, 2, 2, 100, 1, 0};
  for (Index i = 0; i < 6; ++i) close(defaults[i], expected_defaults[i], 0, "Visible named defaults");

  // Validation also occurs after edits, at the conversion boundary.
  named.decrease_factor = 0.5;
  bool rejected         = false;
  try {
    static_cast<void>(named.as_vector());
  } catch (const std::exception& error) { rejected = std::string_view(error.what()).contains("decrease_factor"); }
  require(rejected, "Invalid edited settings must identify the named control");
}

template <typename Configure> void rejects_before_agenda(std::string_view method, Configure configure) {
  Retrieval r;
  configure(r);
  bool rejected = false;
  try {
    r.run(method);
  } catch (const std::exception&) { rejected = true; }
  require(rejected, std::format("Invalid {} configuration was accepted", method));
  require(r.calls == 0, "Invalid configuration executed the agenda");
}

void test_validation() {
  for (const auto method : {"li_m", "gn_m", "not-a-method"}) {
    rejects_before_agenda(method, [](Retrieval&) {});
  }
  for (const auto method : {"li_cg_m", "gn_cg_m"}) {
    rejects_before_agenda(method, [](Retrieval& r) { r.normalization = Vector{1, 1}; });
  }
  for (const Numeric bad :
       {0.0, -1.0, std::numeric_limits<Numeric>::infinity(), std::numeric_limits<Numeric>::quiet_NaN()}) {
    rejects_before_agenda("gn", [bad](Retrieval& r) { r.stop_dx = bad; });
    rejects_before_agenda("gn", [bad](Retrieval& r) { r.normalization = Vector{1, bad}; });
  }
  for (const auto method : {"lm", "ml", "lm_cg", "ml_cg"}) {
    rejects_before_agenda(method, [](Retrieval& r) { r.settings.resize(5); });
    for (Index i = 0; i < 6; ++i) {
      rejects_before_agenda(method, [i](Retrieval& r) { r.settings[i] = -1; });
      rejects_before_agenda(method, [i](Retrieval& r) { r.settings[i] = std::numeric_limits<Numeric>::quiet_NaN(); });
    }
    rejects_before_agenda(method, [](Retrieval& r) { r.settings[1] = 1; });
    rejects_before_agenda(method, [](Retrieval& r) { r.settings[2] = 1; });
    rejects_before_agenda(method, [](Retrieval& r) { r.settings[4] = 0; });
    rejects_before_agenda(method, [](Retrieval& r) { r.settings[0] = r.settings[3] + 1; });
    rejects_before_agenda(method, [](Retrieval& r) { r.settings[4] = r.settings[3] + 1; });
  }
}

void test_skipped_and_reused_outputs() {
  for (Index clear : {0, 1}) {
    Retrieval r;
    r.gain           = matrix(1, 1, {999});
    r.errors         = {"stale error"};
    r.max_start_cost = 0.1;
    r.clear_matrices = clear;
    r.run("lm");
    close(r.diagnostics[0], 99, 0, "Start-cost rejection");
    close(r.diagnostics[1], 2863.0 / 1424, 1e-12, "Rejected initial cost");
    require(r.gain.empty(), "Skipped retrieval returned stale gain");
    require(r.errors.empty(), "Skipped retrieval returned stale errors");
    if (clear) require(r.jac.empty(), "Skipped retrieval ignored clear_matrices");
    for (Index i = 2; i < 5; ++i) require(std::isnan(r.diagnostics[i]), "Skipped diagnostics must remain NaN");

    r.max_start_cost = -1;
    r.run("lm");
    check_affine_solution(r, "lm");
  }
}

template <typename Operation> void rejects_with(Operation operation, std::string_view message) {
  try {
    operation();
  } catch (const std::exception& error) {
    require(std::string_view(error.what()).contains(message),
            std::format("Expected error containing '{}', got '{}'", message, error.what()));
    return;
  }
  throw std::runtime_error(std::format("Expected failure containing '{}'", message));
}

using SolverVector = invlib::Vector<invlib::VectorArchetype<Numeric>>;
using SolverMatrix = invlib::Matrix<invlib::MatrixArchetype<Numeric>>;

SolverVector solver_vector(std::initializer_list<Numeric> values) {
  SolverVector result;
  result.resize(static_cast<unsigned int>(values.size()));
  unsigned int i = 0;
  for (Numeric value : values) result(i++) = value;
  return result;
}

SolverMatrix solver_matrix(unsigned int rows, unsigned int cols, std::initializer_list<Numeric> values) {
  SolverMatrix result;
  result.resize(rows, cols);
  auto value = values.begin();
  for (unsigned int i = 0; i < rows; ++i) {
    for (unsigned int j = 0; j < cols; ++j) result(i, j) = *value++;
  }
  return result;
}

struct IdentityPreconditioner {
  IdentityPreconditioner() = default;
  template <typename MatrixType> explicit IdentityPreconditioner(const MatrixType&) {}
  SolverVector operator()(const SolverVector& value) const { return value; }
};

template <typename Factory> void check_cg_termination(Factory make_solver) {
  const SolverMatrix diagonal = solver_matrix(2, 2, {1, 0, 0, 2});
  const SolverMatrix identity = solver_matrix(2, 2, {1, 0, 0, 1});
  const SolverVector rhs      = solver_vector({1, 1});
  const SolverVector zero     = solver_vector({0, 0});

  // Two distinct eigenvalues require two CG steps for this RHS. Hitting
  // the budget must report failure, not return the first inaccurate iterate.
  for (const Numeric bad :
       {0., -1., std::numeric_limits<Numeric>::quiet_NaN(), std::numeric_limits<Numeric>::infinity()}) {
    rejects_with([&] { static_cast<void>(make_solver(bad, 1)); }, "tolerance");
  }
  for (int bad : {0, -1}) {
    rejects_with([&] { static_cast<void>(make_solver(1e-12, bad)); }, "max_iterations");
  }
  auto limited = make_solver(1e-12, 1);
  rejects_with([&] { static_cast<void>(limited.solve(diagonal, rhs)); }, "iteration limit");
  auto two_steps = make_solver(1e-12, 2);
  for (Index run = 0; run < 2; ++run) {
    const auto solution = two_steps.solve(diagonal, rhs);
    close(solution(0), 1, 1e-14, "CG solution at iteration budget[0]");
    close(solution(1), 0.5, 1e-14, "CG solution at iteration budget[1]");
  }
  // Budgets reset after failures and successes, and a zero RHS is already
  // solved. Test the native solver directly, bypassing OEM's zero shortcut.
  for (Index run = 0; run < 2; ++run) {
    const auto solution = limited.solve(identity, rhs);
    close(solution(0), 1, 0, "Repeated CG solution[0]");
    close(solution(1), 1, 0, "Repeated CG solution[1]");
    const auto zero_solution = limited.solve(identity, zero);
    close(zero_solution(0), 0, 0, "Native CG zero RHS[0]");
    close(zero_solution(1), 0, 0, "Native CG zero RHS[1]");
  }

  for (const Numeric bad : {std::numeric_limits<Numeric>::quiet_NaN(), std::numeric_limits<Numeric>::infinity()}) {
    const SolverVector nonfinite_rhs = solver_vector({bad, 1});
    rejects_with([&] { static_cast<void>(limited.solve(identity, nonfinite_rhs)); }, "finite");
    const SolverMatrix nonfinite_matrix = solver_matrix(2, 2, {bad, 0, 0, 1});
    rejects_with([&] { static_cast<void>(limited.solve(nonfinite_matrix, rhs)); }, "finite");
  }
  for (const Numeric curvature : {0., -1.}) {
    const SolverMatrix breakdown = solver_matrix(2, 2, {curvature, 0, 0, curvature});
    rejects_with([&] { static_cast<void>(limited.solve(breakdown, rhs)); }, "curvature");
  }
}

struct NeverConvergedCGSettings {
  explicit NeverConvergedCGSettings(double) {}
  SolverVector start_vector(const SolverVector& rhs) const { return 0.0 * rhs; }
  bool         converged(const SolverVector&, const SolverVector&) const { return false; }
};

void test_cg_termination() {
  check_cg_termination([](Numeric tolerance, int budget) { return invlib::ConjugateGradient<>(tolerance, 0, budget); });
  const IdentityPreconditioner identity;
  check_cg_termination([&](Numeric tolerance, int budget) {
    return invlib::PreconditionedConjugateGradient<IdentityPreconditioner, true>(identity, tolerance, 0, budget);
  });
  check_cg_termination([](Numeric tolerance, int budget) {
    return invlib::PreconditionedConjugateGradient<IdentityPreconditioner, false>(tolerance, 0, budget);
  });

  // The safety budget belongs to the solve loop, even when a custom settings
  // functor replaces the default convergence predicate.
  invlib::ConjugateGradient<NeverConvergedCGSettings> custom(1e-12, 0, 1);
  const SolverMatrix                                  diagonal = solver_matrix(2, 2, {1, 0, 0, 2});
  const SolverVector                                  rhs      = solver_vector({1, 1});
  rejects_with([&] { static_cast<void>(custom.solve(diagonal, rhs)); }, "iteration limit");

  // A fixed-step policy must also stop when the exact solution is reached,
  // before the next conjugate-direction update attempts a 0/0 division.
  invlib::ConjugateGradient<invlib::CGStepLimit<3>> fixed_steps(1e-12);
  const SolverMatrix                                identity_matrix = solver_matrix(2, 2, {1, 0, 0, 1});
  const auto                                        exact           = fixed_steps.solve(identity_matrix, rhs);
  close(exact(0), 1, 0, "Fixed-step CG exact solution[0]");
  close(exact(1), 1, 0, "Fixed-step CG exact solution[1]");

  for (const auto method : {"li_cg", "li_cg_m", "gn_cg", "gn_cg_m", "lm_cg", "ml_cg"}) {
    for (const Numeric bad : {std::numeric_limits<Numeric>::quiet_NaN(), std::numeric_limits<Numeric>::infinity()}) {
      Retrieval r;
      r.y[0]     = bad;
      r.max_iter = 1;
      r.run(method);
      close(r.diagnostics[0], 9, 0, std::format("{} with measurement {} must fail", method, bad));
      require(not r.errors.empty(), "Nonfinite CG input lost the failure reason");
      require(r.gain.empty(), "Failed CG retrieval returned gain");
    }
  }
}

void test_lm_trial_limit() {
  const SolverMatrix curvature = solver_matrix(1, 1, {1});
  const SolverVector initial   = solver_vector({1});
  using Optimizer              = invlib::LevenbergMarquardt<Numeric, SolverMatrix>;
  struct RejectTrials {
    Index   calls = 0;
    Numeric cost_function(const SolverVector&, bool = false) { return calls++ == 0 ? 0 : 1; }
  } rejected;

  Optimizer limited(curvature);
  rejects_with([&] { limited.set_maximum_trials(0); }, "trial");
  limited.set_maximum_trials(3);
  require(limited.get_maximum_trials() == 3, "Configured LM trial budget was lost");
  limited.set_lambda(1);
  limited.set_lambda_increase(1.01);
  rejects_with([&] { static_cast<void>(limited.step(initial, initial, curvature, rejected)); },
               "Levenberg-Marquardt trial limit");
  require(rejected.calls == 4, "LM exceeded its three-trial budget or stopped before spending it");
  require(limited.get_stop_reason() == invlib::LMStopReason::TrialLimit,
          "LM trial exhaustion lost its explicit stop reason");
  require(limited.stop_iteration() and not limited.converged(), "LM trial exhaustion reported convergence");

  struct QuadraticCost {
    Numeric cost_function(const SolverVector& x, bool = false) { return 0.5 * x(0) * x(0); }
  } quadratic;
  Optimizer one_trial(curvature);
  one_trial.set_maximum_trials(1);
  one_trial.set_lambda(1);
  one_trial.set_lambda_threshold(0.01);
  SolverVector state = initial;
  // Both outer steps succeed on their last allowed trial. Sharing the trial
  // counter between steps would wrongly reject the second call.
  for (Index step = 0; step < 2; ++step) {
    const Numeric expected  = state(0) * one_trial.get_lambda() / (1 + one_trial.get_lambda());
    const auto    dx        = one_trial.step(state, state, curvature, quadratic);
    state                  += dx;
    close(state(0), expected, 1e-14, "LM per-step trial budget");
  }

  for (const auto method : {"lm", "ml", "lm_cg", "ml_cg"}) {
    for (const bool stalled : {false, true}) {
      Retrieval r;
      quadratic_model(r);
      r.x        = Vector{0.1};
      r.max_iter = 1;
      r.settings = OEMLMSettings{.initial_damping   = 0,
                                 .increase_factor   = stalled ? std::nextafter(1., 2.) : 1.0001,
                                 .damping_threshold = stalled ? std::nextafter(0., 1.) : 0.1}
                       .as_vector();
      r.run(method);
      close(r.diagnostics[0], 9, 0, "LM retry failure must not report convergence");
      const std::string_view reason = stalled ? "damping did not increase" : "Levenberg-Marquardt trial limit";
      require(stdr::any_of(r.errors, [&](const String& error) { return std::string_view(error).contains(reason); }),
              "LM retry failure lost its reason");
      require(r.gain.empty(), "Failed LM retrieval returned gain");
      // Initial Jacobian and initial LM cost precede the trial evaluations.
      require(r.calls <= 102, "OEM exceeded the default LM trial budget");
    }
  }
}

void test_lm_stop_reasons() {
  const SolverMatrix curvature = solver_matrix(1, 1, {1});
  const SolverVector initial   = solver_vector({1});
  using Optimizer              = invlib::LevenbergMarquardt<Numeric, SolverMatrix>;
  struct QuadraticCost {
    Numeric cost_function(const SolverVector& x, bool = false) { return 0.5 * x(0) * x(0); }
  } quadratic;

  Optimizer ordinary(curvature);
  ordinary.set_lambda(1);
  const auto accepted = ordinary.step(initial, initial, curvature, quadratic);
  close(accepted(0), -0.5, 0, "Ordinary LM step");
  require(ordinary.get_stop_reason() == invlib::LMStopReason::None,
          "An accepted LM step prematurely acquired a stop reason");
  require(not ordinary.stop_iteration() and not ordinary.converged(), "An accepted LM step prematurely stopped");

  Optimizer exhausted(curvature);
  exhausted.set_lambda(1e20);
  exhausted.set_lambda_maximum(1e20);
  const auto rejected = exhausted.step(initial, initial, curvature, quadratic);
  close(rejected(0), 0, 0, "Exhausted LM returns no rejected displacement");
  close(exhausted.get_lambda(), 1e20, 0, "Exhausted LM retains actual damping");
  require(exhausted.get_stop_reason() == invlib::LMStopReason::DampingLimit,
          "LM maximum damping lost its explicit stop reason");
  require(exhausted.stop_iteration() and not exhausted.converged(), "LM maximum damping reported convergence");

  struct LargeBaselineCost {
    Numeric cost_function(const SolverVector& x, bool = false) { return 1e20 + 0.5 * x(0) * x(0); }
  } large_baseline;
  Optimizer unresolved(curvature);
  unresolved.set_lambda(1);
  unresolved.set_lambda_maximum(1);
  const auto unresolved_step = unresolved.step(initial, initial, curvature, large_baseline);
  // Both damped and undamped cost changes round to zero. A large constant
  // objective offset must not disguise the remaining undamped state error.
  close(unresolved_step(0), 0, 0, "Unresolved nonstationary step must not be applied");
  require(unresolved.get_stop_reason() == invlib::LMStopReason::DampingLimit,
          "A large constant cost baseline incorrectly certified stationarity");
  require(not unresolved.converged(), "Unresolved nonstationary objective reported convergence");

  struct InaccurateModelCost {
    Numeric cost_function(const SolverVector& x, bool = false) { return 0.125 * x(0) * x(0); }
  } inaccurate_model;
  Optimizer weak_reduction(curvature);
  weak_reduction.set_lambda(1);
  weak_reduction.set_lambda_maximum(1);
  const auto weak_step = weak_reduction.step(initial, initial, curvature, inaccurate_model);
  // Actual reduction is positive but only one quarter of the supplied model's
  // prediction. The rejected step used to leak out because its ratio was >0.
  close(weak_step(0), 0, 0, "Positive but inadequate reduction must not be applied");
  require(weak_reduction.get_stop_reason() == invlib::LMStopReason::DampingLimit,
          "Positive but inadequate reduction lost its failure reason");
  require(not weak_reduction.converged(), "Positive but inadequate reduction reported convergence");

  struct NonfiniteTrials {
    Index   calls = 0;
    Numeric cost_function(const SolverVector&, bool = false) {
      return calls++ == 0 ? 1 : std::numeric_limits<Numeric>::quiet_NaN();
    }
  } nonfinite_trials;
  Optimizer nonfinite(curvature);
  nonfinite.set_lambda(1);
  nonfinite.set_lambda_maximum(1);
  const auto nonfinite_step = nonfinite.step(initial, initial, curvature, nonfinite_trials);
  close(nonfinite_step(0), 0, 0, "A NaN-cost trial must not be applied");
  require(nonfinite.get_stop_reason() == invlib::LMStopReason::DampingLimit,
          "A NaN-cost trial escaped without a failure reason");
  require(nonfinite.stop_iteration() and not nonfinite.converged(), "A NaN-cost trial reported convergence");
  require(nonfinite_trials.calls == 2, "A NaN-cost trial exceeded the maximum damping budget");

  struct OffsetQuadraticCost {
    Numeric cost_function(const SolverVector& x, bool = false) { return 1 + 0.5 * x(0) * x(0); }
  } offset_quadratic;
  for (const Numeric state : {0., 1e-9}) {
    Optimizer stationary(curvature);
    stationary.set_lambda(0);
    const auto x  = solver_vector({state});
    const auto dx = stationary.step(x, x, curvature, offset_quadratic);
    require(stationary.get_stop_reason() == invlib::LMStopReason::Stationary,
            "Stationary LM point lost its explicit stop reason");
    require(stationary.stop_iteration() and stationary.converged(), "Stationary LM point did not converge");
    close(x(0) + dx(0), 0, 1e-8, "LM stationary state");
  }

  struct RejectTrials {
    Index   calls = 0;
    Numeric cost_function(const SolverVector&, bool = false) { return calls++ == 0 ? 0 : 1; }
  } reject_all;
  Optimizer stalled(curvature);
  stalled.set_lambda(std::nextafter(0., 1.));
  stalled.set_lambda_threshold(std::nextafter(0., 1.));
  stalled.set_lambda_increase(std::nextafter(1., 2.));
  rejects_with([&] { static_cast<void>(stalled.step(initial, initial, curvature, reject_all)); },
               "damping did not increase");
  require(stalled.get_stop_reason() == invlib::LMStopReason::DampingStalled,
          "LM unchanged damping lost its explicit stop reason");
  require(stalled.stop_iteration() and not stalled.converged(), "LM unchanged damping reported convergence");
}

void test_generic_minimize_outcomes() {
  const SolverMatrix curvature = solver_matrix(1, 1, {1});
  const SolverVector initial   = solver_vector({1});
  using Optimizer              = invlib::LevenbergMarquardt<Numeric, SolverMatrix>;
  struct QuadraticCost {
    Numeric      baseline        = 0;
    Index        criterion_calls = 0;
    Numeric      cost_function(const SolverVector& x, bool = false) { return baseline + 0.5 * x(0) * x(0); }
    SolverVector gradient(const SolverVector& x) { return x; }
    SolverMatrix Hessian(const SolverVector&) { return solver_matrix(1, 1, {1}); }
    Numeric      criterion(const SolverVector&, const SolverVector& dx) {
      ++criterion_calls;
      return dx.norm();
    }
  } quadratic;
  SolverVector result;
  Optimizer    ordinary(curvature);
  ordinary.set_lambda(1);
  require(invlib::minimize(quadratic, ordinary, initial, result, 10, 1e-12) == 0,
          "Generic LM minimization failed to converge");
  close(result(0), 0, 0, "Generic LM solution");
  require(ordinary.get_stop_reason() == invlib::LMStopReason::Stationary,
          "Generic minimization failed to respect stationary termination");

  Optimizer exhausted(curvature);
  exhausted.set_lambda(1e20);
  exhausted.set_lambda_maximum(1e20);
  quadratic.criterion_calls = 0;
  require(invlib::minimize(quadratic, exhausted, initial, result, 10, 1e-12) == 1,
          "Generic minimization reported damping exhaustion as convergence");
  close(result(0), 1, 0, "Generic failed LM state");
  require(quadratic.criterion_calls == 0, "Generic minimization tested a rejected zero step for convergence");

  for (const unsigned int iterations : {0u, 1u}) {
    Optimizer limited(curvature);
    limited.set_lambda(1);
    require(invlib::minimize(quadratic, limited, initial, result, iterations, 1e-12) == 1,
            "Generic minimization concealed iteration exhaustion");
    close(result(0), iterations == 0 ? 1 : 0.5, 0, "Generic iteration-limited state");
  }

  const SolverVector near_solution = solver_vector({1e-9});
  quadratic.baseline               = 1;
  for (const Numeric tolerance : {1e-8, 1e-12}) {
    Optimizer stationary(curvature);
    stationary.set_lambda(0);
    const auto status = invlib::minimize(quadratic, stationary, near_solution, result, 10, tolerance);
    require(stationary.converged(), "Generic stationary fixture did not reach an LM stop reason");
    require(status == (tolerance > 1e-9 ? 0 : 1),
            "Generic minimization ignored its own tolerance after optimizer termination");
    close(result(0), 0, 0, "Generic minimization applies a verified stationary step");
  }

  // Existing custom minimizers need not implement the optional stop hooks.
  struct CustomMinimizer {
    SolverVector step(const SolverVector&, const SolverVector& gradient, const SolverMatrix&, QuadraticCost&) {
      return -0.5 * gradient;
    }
  } custom;
  require(invlib::minimize(quadratic, custom, initial, result, 1, 1.) == 0,
          "Generic minimization broke a custom minimizer without stop hooks");
  close(result(0), 0.5, 0, "Custom minimizer state");
  require(invlib::minimize(quadratic, custom, initial, result, 1, 1e-12) == 1,
          "Generic custom minimization concealed iteration exhaustion");
}

}  // namespace

int main(int argc, char** argv) try {
  require(argc == 2, "Usage: test_oem_methods METHOD|settings|validation|termination|outcomes");
  const std::string_view selected{argv[1]};
  if (selected == "settings") {
    test_lm_settings();
    test_named_settings();
    test_skipped_and_reused_outputs();
  } else if (selected == "validation") {
    test_validation();
  } else if (selected == "termination") {
    test_cg_termination();
    test_lm_trial_limit();
  } else if (selected == "outcomes") {
    test_lm_outcomes();
    test_lm_stop_reasons();
    test_generic_minimize_outcomes();
  } else {
    require(stdr::find(methods, selected) != methods.end(), "Unknown test method");
    test_affine(selected);
    test_exact_start(selected);
    test_disabled_start_cost(selected);
    test_runtime_failure(selected);
    test_nonlinear(selected);
    test_underdetermined(selected);
  }
  return EXIT_SUCCESS;
} catch (const std::exception& error) {
  std::println(stderr, "OEM regression ({}): {}", argc > 1 ? argv[1] : "", error.what());
  return EXIT_FAILURE;
}
