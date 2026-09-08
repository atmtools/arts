#include "oem_settings.h"

#include <debug.h>

#include <array>
#include <cmath>
#include <format>
#include <string_view>
#include <utility>

void OEMLMSettings::validate() const {
  const std::array<std::pair<std::string_view, Numeric>, 6> fields{{
      {"initial_damping", initial_damping},
      {"decrease_factor", decrease_factor},
      {"increase_factor", increase_factor},
      {"maximum_damping", maximum_damping},
      {"damping_threshold", damping_threshold},
      {"convergence_damping_limit", convergence_damping_limit},
  }};
  for (const auto& [name, value] : fields) {
    ARTS_USER_ERROR_IF(
        !std::isfinite(value) || value < 0, "OEMLMSettings.{} must be finite and nonnegative; got {}.", name, value)
  }
  ARTS_USER_ERROR_IF(
      decrease_factor <= 1, "OEMLMSettings.decrease_factor divides damping and must be > 1; got {}.", decrease_factor)
  ARTS_USER_ERROR_IF(increase_factor <= 1,
                     "OEMLMSettings.increase_factor multiplies damping and must be > 1; got {}.",
                     increase_factor)
  ARTS_USER_ERROR_IF(damping_threshold <= 0,
                     "OEMLMSettings.damping_threshold must be > 0 so rejected undamped steps can restart; got {}.",
                     damping_threshold)
  ARTS_USER_ERROR_IF(damping_threshold > maximum_damping,
                     "OEMLMSettings.damping_threshold ({}) must not exceed maximum_damping ({}).",
                     damping_threshold,
                     maximum_damping)
  ARTS_USER_ERROR_IF(initial_damping > maximum_damping,
                     "OEMLMSettings.initial_damping ({}) must not exceed maximum_damping ({}).",
                     initial_damping,
                     maximum_damping)
}

Vector OEMLMSettings::as_vector() const {
  validate();
  return {
      initial_damping, decrease_factor, increase_factor, maximum_damping, damping_threshold, convergence_damping_limit};
}

OEMLMSettings OEMLMSettings::from_vector(const Vector& values) {
  ARTS_USER_ERROR_IF(values.size() != 6,
                     "lm_ga_settings must contain 6 values in this order: initial_damping, decrease_factor, "
                     "increase_factor, maximum_damping, damping_threshold, convergence_damping_limit; got {}.",
                     values.size())
  OEMLMSettings settings{.initial_damping           = values[0],
                         .decrease_factor           = values[1],
                         .increase_factor           = values[2],
                         .maximum_damping           = values[3],
                         .damping_threshold         = values[4],
                         .convergence_damping_limit = values[5]};
  settings.validate();
  return settings;
}

std::string OEMLMSettings::repr() const {
  return std::format(
      "OEMLMSettings(initial_damping={}, decrease_factor={}, increase_factor={}, "
      "maximum_damping={}, damping_threshold={}, convergence_damping_limit={})",
      initial_damping,
      decrease_factor,
      increase_factor,
      maximum_damping,
      damping_threshold,
      convergence_damping_limit);
}

std::string OEMLMSettings::describe() const {
  validate();
  return std::format(
      "{}\n\n"
      "The first trial uses damping {}. Zero gives a Gauss-Newton step; increasing damping constrains the step.\n"
      "When damping is reduced, it is divided by {}. A reduction below {} switches it to zero.\n"
      "A rejected trial below {} restarts at {}; otherwise damping is multiplied by {}, up to {}.\n"
      "Another rejected trial at that maximum stops the retrieval.\n"
      "Accepted trials do not always reduce damping, particularly after a rejection.\n"
      "The ordinary stop_dx test is enabled when updated damping is <= {}. This limit is not a state-step tolerance.\n"
      "{}\n"
      "These controls determine trial steps. The prior and measurement covariances determine statistical weights.\n",
      repr(),
      initial_damping,
      decrease_factor,
      damping_threshold,
      damping_threshold,
      damping_threshold,
      increase_factor,
      maximum_damping,
      convergence_damping_limit,
      convergence_damping_limit == 0
          ? "A zero convergence_damping_limit waits for damping to reach zero before enabling the ordinary stopping test."
          : "A positive convergence_damping_limit allows stopping while damping still constrains the step; inspect the costs and residuals.");
}
