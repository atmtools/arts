#pragma once

#include <covariance_matrix.h>
#include <enums.h>
#include <matpack.h>
#include <retrieval_target.h>
#include <xml_io_stream_aggregate.h>

#include <string>

/** Named controls for OEM Levenberg--Marquardt damping.
 *
 * These defaults are an explicit starting configuration, not a convergence
 * guarantee. Call validate() after editing fields. This type does not depend on the OEM solver.
 */
struct LevenbergMarquardtSettings {
  /** Damping for the first trial; zero starts with a Gauss--Newton step. */
  Numeric initial_damping = 10;
  /** Divisor for a damping reduction, greater than one. */
  Numeric decrease_factor = 2;
  /** Multiplier for a damping increase, greater than one. */
  Numeric increase_factor = 2;
  /** Damping limit at which another rejected trial stops the retrieval. */
  Numeric maximum_damping = 100;
  /** Restart value after a low-damping rejection; reductions below it use zero. */
  Numeric damping_threshold = 1;
  /** Upper updated damping at which the ordinary stop_dx test is enabled. */
  Numeric convergence_damping_limit = 0;

  /** Maximum linear-solve trials per outer iteration, including stationarity checks. */
  Index maximum_trials = 100;

  /** Reject invalid individual values and inconsistent combinations. */
  void validate() const;

  /** Show every field, including defaults, without requiring valid settings. */
  [[nodiscard]] std::string repr() const;

  /** Validate and explain the configured controls and their interactions. */
  [[nodiscard]] std::string describe() const;
};

template <> struct std::formatter<LevenbergMarquardtSettings> {
  format_tags                      tags;
  [[nodiscard]] constexpr auto&    inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto&    inner_fmt() const { return *this; }
  constexpr auto                   parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class FmtContext> auto format(const LevenbergMarquardtSettings& value, FmtContext& ctx) const {
    return tags.format(ctx, value.repr());
  }
};

template <> struct xml_io_stream_name<LevenbergMarquardtSettings> {
  static constexpr std::string_view name = "LevenbergMarquardtSettings";
};

template <> struct xml_io_stream_aggregate<LevenbergMarquardtSettings> {
  static constexpr bool value = true;
};

/** Calculation controls shared by full and reduced optimal estimation.
 * Validate once at the calculation boundary, after editing the configuration.
 */
struct OptimalEstimationSettings {
  OptimalEstimationMethod    method           = OptimalEstimationMethod::gn;
  Index                      max_iter         = 10;
  Numeric                    stop_dx          = 0.01;
  Numeric                    max_start_cost   = std::numeric_limits<Numeric>::infinity();
  Numeric                    cg_tolerance     = 1e-10;
  Index                      cg_max_iter      = 0;
  LevenbergMarquardtSettings lm               = {};
  Index                      display_progress = 0;
  bool                       clear_matrices   = false;

  void                      validate() const;
  [[nodiscard]] std::string repr() const;
};

template <> struct std::formatter<OptimalEstimationSettings> {
  format_tags                   tags;
  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }
  constexpr auto                parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class Context> auto format(const OptimalEstimationSettings& v, Context& ctx) const {
    return tags.format(ctx, v.repr());
  }
};
template <> struct xml_io_stream_name<OptimalEstimationSettings> {
  static constexpr std::string_view name = "OptimalEstimationSettings";
};
template <> struct xml_io_stream_aggregate<OptimalEstimationSettings> {
  static constexpr bool value = true;
};

/** OEM results. Costs are per measurement and NaN when unavailable.
 * Iterations counts completed outer iterations, including zero when not run.
 */
struct OptimalEstimationDiagnostics {
  OptimalEstimationStatus status           = OptimalEstimationStatus::NotRun;
  Numeric                 initial_cost     = std::numeric_limits<Numeric>::quiet_NaN();
  Numeric                 final_cost       = std::numeric_limits<Numeric>::quiet_NaN();
  Numeric                 measurement_cost = std::numeric_limits<Numeric>::quiet_NaN();
  Index                   iterations       = 0;
  Vector                  lm_ga_history;
  ArrayOfString           errors;
};

template <> struct std::formatter<OptimalEstimationDiagnostics> {
  format_tags                      tags;
  [[nodiscard]] constexpr auto&    inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto&    inner_fmt() const { return *this; }
  constexpr auto                   parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class FmtContext> auto format(const OptimalEstimationDiagnostics& v, FmtContext& ctx) const {
    return tags.format(ctx,
                       "OptimalEstimationDiagnostics(status="sv,
                       toString(v.status),
                       ", initial_cost="sv,
                       v.initial_cost,
                       ", final_cost="sv,
                       v.final_cost,
                       ", measurement_cost="sv,
                       v.measurement_cost,
                       ", iterations="sv,
                       v.iterations,
                       ", lm_ga_history="sv,
                       v.lm_ga_history,
                       ", errors="sv,
                       v.errors,
                       ")"sv);
  }
};

template <> struct xml_io_stream_name<OptimalEstimationDiagnostics> {
  static constexpr std::string_view name = "OptimalEstimationDiagnostics";
};
template <> struct xml_io_stream_aggregate<OptimalEstimationDiagnostics> {
  static constexpr bool value = true;
};

/** Owning numerical problem and results for oemCalc/oemCalcReduced.
 * Physical forward-model data and target metadata remain in the workspace.
 * Per-call solver scratch has automatic lifetime; covariance caches belong
 * to the covariance objects and can be released by clear_auxiliary().
 */
struct OptimalEstimationData {
  Vector           measurement_vec;
  Vector           model_state_vec_apriori;
  CovarianceMatrix model_state_covmat;
  CovarianceMatrix measurement_vec_error_covmat;

  Vector      model_state_vec;
  Vector      measurement_vec_fit;
  Matrix      measurement_jac;
  BlockMatrix model_state_basis_mat{Matrix{}};
  BlockMatrix measurement_basis_mat{Matrix{}};
  Vector      model_state_covmat_normalization;
  Vector      measurement_vec_normalization;

  Matrix                       measurement_gain_mat;
  Matrix                       measurement_averaging_kernel;
  Matrix                       observation_error_covmat;
  Matrix                       smoothing_error_covmat;
  OptimalEstimationDiagnostics diagnostics;
  Vector                       basis_singular_values;
  Numeric                      basis_lost_dofs             = std::numeric_limits<Numeric>::quiet_NaN();
  Numeric                      basis_lost_information_bits = std::numeric_limits<Numeric>::quiet_NaN();

  /** Pending per-target covariance blocks used during target-based setup. */
  JacobianTargetsDiagonalCovarianceMatrixMap covmat_diagonal_blocks;

  [[nodiscard]] bool checked() const noexcept { return checked_; }
  void               uncheck() noexcept { checked_ = false; }
  void               check(const JacobianTargets* targets = nullptr);
  /** Validate once, reusing the checked status until invalidated.
   *
   * A checked object deliberately skips revalidation, including the comparison
   * against jac_targets. Since jac_targets is a separate workspace variable,
   * editing it cannot clear this object's checked status; a resulting mismatch
   * is reported by oemCalc when the agenda returns a Jacobian of the wrong shape.
   */
  void ensure_checked(const JacobianTargets& targets) {
    if (not checked_) check(&targets);
  }
  void require_unchecked(std::string_view operation) const;

  /** Release recomputable products while preserving the problem, state,
   * selected bases, normalization and diagnostics. */
  void clear_auxiliary();
  /** Reset the entire retrieval, including inputs and diagnostics. */
  void clear();

 private:
  bool checked_ = false;
};

template <> struct std::formatter<OptimalEstimationData> {
  format_tags                   tags;
  constexpr auto                parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class Context> auto format(const OptimalEstimationData& v, Context& ctx) const {
    return tags.format(ctx,
                       "OptimalEstimationData(measurements="sv,
                       v.measurement_vec.size(),
                       ", states="sv,
                       v.model_state_vec_apriori.size(),
                       ", diagnostics="sv,
                       v.diagnostics,
                       ")"sv);
  }
};
template <> struct xml_io_stream_name<OptimalEstimationData> {
  static constexpr std::string_view name = "OptimalEstimationData";
};
template <> struct xml_io_stream<OptimalEstimationData> {
  static constexpr std::string_view type_name = "OptimalEstimationData";
  static void write(std::ostream&, const OptimalEstimationData&, bofstream* = nullptr, std::string_view = "");
  static void read(std::istream&, OptimalEstimationData&, bifstream* = nullptr);
};
