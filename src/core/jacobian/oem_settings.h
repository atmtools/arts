#pragma once

#include <enums.h>
#include <matpack.h>
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
