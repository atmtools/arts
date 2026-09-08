#pragma once

#include <matpack.h>

#include <string>

/** Named controls for OEM Levenberg--Marquardt damping.
 *
 * These defaults are an explicit starting configuration, not a convergence
 * guarantee. The existing OEM vector argument retains its empty default.
 * Call validate() after editing fields, or as_vector() to validate and convert
 * to the workspace interface. This type does not depend on the OEM solver.
 */
struct OEMLMSettings {
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
  [[nodiscard]] static OEMLMSettings from_vector(const Vector& values);

  /** Show every field, including defaults, without requiring valid settings. */
  [[nodiscard]] std::string repr() const;

  /** Validate and explain the configured controls and their interactions. */
  [[nodiscard]] std::string describe() const;
};
