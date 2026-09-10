#pragma once

#include <matpack.h>
#include <xml_io_stream_aggregate.h>

#include <string>

/** Named controls for OEM Levenberg--Marquardt damping.
 *
 * These defaults are an explicit starting configuration, not a convergence
 * guarantee. The existing OEM vector argument retains its empty default.
 * Call validate() after editing fields, or as_vector() to validate and convert
 * to the workspace interface. This type does not depend on the OEM solver.
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

  /** Validate and return the legacy six-element lm_ga_settings vector. */
  [[nodiscard]] Vector as_vector() const;

  /** Parse and validate the legacy six-element lm_ga_settings vector. */
  [[nodiscard]] static LevenbergMarquardtSettings from_vector(const Vector& values);

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
