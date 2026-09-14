#include "oem_settings.h"

#include <debug.h>

#include <array>
#include <cmath>
#include <format>
#include <string_view>
#include <utility>

void LevenbergMarquardtSettings::validate() const {
  const std::array<std::pair<std::string_view, Numeric>, 6> fields{{
      {"initial_damping", initial_damping},
      {"decrease_factor", decrease_factor},
      {"increase_factor", increase_factor},
      {"maximum_damping", maximum_damping},
      {"damping_threshold", damping_threshold},
      {"convergence_damping_limit", convergence_damping_limit},
  }};
  for (const auto& [name, value] : fields) {
    ARTS_USER_ERROR_IF(!std::isfinite(value) || value < 0,
                       "LevenbergMarquardtSettings.{} must be finite and nonnegative; got {}.",
                       name,
                       value)
  }
  ARTS_USER_ERROR_IF(decrease_factor <= 1,
                     "LevenbergMarquardtSettings.decrease_factor divides damping and must be > 1; got {}.",
                     decrease_factor)
  ARTS_USER_ERROR_IF(increase_factor <= 1,
                     "LevenbergMarquardtSettings.increase_factor multiplies damping and must be > 1; got {}.",
                     increase_factor)
  ARTS_USER_ERROR_IF(
      damping_threshold <= 0,
      "LevenbergMarquardtSettings.damping_threshold must be > 0 so rejected undamped steps can restart; got {}.",
      damping_threshold)
  ARTS_USER_ERROR_IF(damping_threshold > maximum_damping,
                     "LevenbergMarquardtSettings.damping_threshold ({}) must not exceed maximum_damping ({}).",
                     damping_threshold,
                     maximum_damping)
  ARTS_USER_ERROR_IF(initial_damping > maximum_damping,
                     "LevenbergMarquardtSettings.initial_damping ({}) must not exceed maximum_damping ({}).",
                     initial_damping,
                     maximum_damping)
}

std::string LevenbergMarquardtSettings::repr() const {
  return std::format(
      "LevenbergMarquardtSettings(initial_damping={}, decrease_factor={}, increase_factor={}, "
      "maximum_damping={}, damping_threshold={}, convergence_damping_limit={})",
      initial_damping,
      decrease_factor,
      increase_factor,
      maximum_damping,
      damping_threshold,
      convergence_damping_limit);
}

std::string LevenbergMarquardtSettings::describe() const {
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

void OptimalEstimationData::clear() { *this = OptimalEstimationData{}; }

void OptimalEstimationData::clear_auxiliary() {
  covmat_diagonal_blocks       = JacobianTargetsDiagonalCovarianceMatrixMap{};
  measurement_vec_fit          = Vector{};
  measurement_jac              = Matrix{};
  measurement_gain_mat         = Matrix{};
  measurement_averaging_kernel = Matrix{};
  observation_error_covmat     = Matrix{};
  smoothing_error_covmat       = Matrix{};
  basis_singular_values        = Vector{};
  basis_lost_dofs              = NAN;
  basis_lost_information_bits  = NAN;
  model_state_covmat.clear_cache();
  measurement_vec_error_covmat.clear_cache();
}

void OptimalEstimationData::require_unchecked(std::string_view operation) const {
  ARTS_USER_ERROR_IF(checked_,
                     "Cannot {} on checked OptimalEstimationData. Call oem.uncheck() before replacing members, "
                     "then oem.check() after editing. Reading members and editing their existing values are allowed.",
                     operation)
}

void OptimalEstimationData::check(const JacobianTargets* targets) {
  uncheck();
  std::string errors;
  const auto  error = [&](std::string message) { errors += "\n - " + message; };
  const Size  n = model_state_vec_apriori.size(), m = measurement_vec.size();
  if (n == 0) error("model_state_vec_apriori must be nonempty.");
  if (m == 0) error("measurement_vec must be nonempty.");
  const auto vector = [&](std::string_view name, const Vector& value, Size size, bool optional, bool positive = false) {
    if ((!optional or !value.empty()) and value.size() != size)
      error(std::format("{} has {} elements; expected {}{}.", name, value.size(), size, optional ? " or empty" : ""));
    for (Size i = 0; i < value.size(); ++i)
      if (!std::isfinite(value[i]) or (positive and value[i] <= 0)) {
        error(
            std::format("{}[{}] = {} must be finite{}.", name, i, value[i], positive ? " and strictly positive" : ""));
        break;
      }
  };
  vector("model_state_vec_apriori", model_state_vec_apriori, n, false);
  vector("measurement_vec", measurement_vec, m, false);
  vector("model_state_vec", model_state_vec, n, true);
  vector("measurement_vec_fit", measurement_vec_fit, m, true);
  // This preflight is shared by full and reduced OEM. The selected calculation
  // checks which normalization dimension applies to its solve.
  const auto normalization = [&](std::string_view name, const Vector& value, Size full_size, Size reduced_size) {
    if (value.empty()) return;
    if (value.size() != full_size and value.size() != reduced_size)
      error(std::format("{} has {} elements; expected {} (full space){} or empty.",
                        name,
                        value.size(),
                        full_size,
                        reduced_size ? std::format(" or {} (reduced space)", reduced_size) : ""));
    vector(name, value, value.size(), false, true);
  };
  normalization("model_state_covmat_normalization",
                model_state_covmat_normalization,
                n,
                model_state_basis_mat.not_null() ? model_state_basis_mat.ncols() : 0);
  normalization("measurement_vec_normalization",
                measurement_vec_normalization,
                m,
                measurement_basis_mat.not_null() ? measurement_basis_mat.nrows() : 0);
  if (!measurement_jac.empty()) {
    if (measurement_jac.nrows() != static_cast<Index>(m) or measurement_jac.ncols() != static_cast<Index>(n))
      error(std::format("measurement_jac is {} by {}; expected {} by {} or empty.",
                        measurement_jac.nrows(),
                        measurement_jac.ncols(),
                        m,
                        n));
    for (auto value : stdr::subrange(measurement_jac.elem_begin(), measurement_jac.elem_end()))
      if (!std::isfinite(value)) {
        error("measurement_jac contains a non-finite value.");
        break;
      }
  }
  const auto covariance = [&](std::string_view name, const CovarianceMatrix& value, Size size) {
    try {
      value.validate(static_cast<Index>(size));
    } catch (const std::exception& e) { error(std::format("{}: {}", name, e.what())); }
  };
  covariance("model_state_covmat", model_state_covmat, n);
  covariance("measurement_vec_error_covmat", measurement_vec_error_covmat, m);
  if (model_state_basis_mat.not_null() and (model_state_basis_mat.nrows() != 0 or model_state_basis_mat.ncols() != 0)) {
    if (model_state_basis_mat.nrows() != static_cast<Index>(n) or model_state_basis_mat.ncols() <= 0 or
        model_state_basis_mat.ncols() > static_cast<Index>(n))
      error(std::format("model_state_basis_mat must have {} rows and between 1 and {} columns.", n, n));
    if (!model_state_basis_mat.is_finite()) error("model_state_basis_mat contains non-finite values.");
  }
  if (measurement_basis_mat.not_null() and (measurement_basis_mat.nrows() != 0 or measurement_basis_mat.ncols() != 0)) {
    if (measurement_basis_mat.ncols() != static_cast<Index>(m) or measurement_basis_mat.nrows() <= 0 or
        measurement_basis_mat.nrows() > static_cast<Index>(m))
      error(std::format("measurement_basis_mat must have {} columns and between 1 and {} rows.", m, m));
    if (!measurement_basis_mat.is_finite()) error("measurement_basis_mat contains non-finite values.");
  }
  if (targets) {
    if (!targets->finalized)
      error("jac_targets must be finalized.");
    else if (targets->x_size() != n)
      error(std::format("jac_targets maps {} state elements but the prior contains {}.", targets->x_size(), n));
  }
  ARTS_USER_ERROR_IF(
      !errors.empty(),
      "OptimalEstimationData check failed:{}\nCorrect these inputs and call check() again. The object remains unchecked.",
      errors)
  checked_ = true;
}

namespace {
auto oem_xml_members(auto& v) {
  return std::tie(v.measurement_vec,
                  v.model_state_vec_apriori,
                  v.model_state_covmat,
                  v.measurement_vec_error_covmat,
                  v.model_state_vec,
                  v.measurement_vec_fit,
                  v.measurement_jac,
                  v.model_state_basis_mat,
                  v.measurement_basis_mat,
                  v.model_state_covmat_normalization,
                  v.measurement_vec_normalization,
                  v.measurement_gain_mat,
                  v.measurement_averaging_kernel,
                  v.observation_error_covmat,
                  v.smoothing_error_covmat,
                  v.diagnostics,
                  v.basis_singular_values,
                  v.basis_lost_dofs,
                  v.basis_lost_information_bits,
                  v.covmat_diagonal_blocks);
}
}  // namespace
void xml_io_stream<OptimalEstimationData>::write(std::ostream&                os,
                                                 const OptimalEstimationData& value,
                                                 bofstream*                   binary,
                                                 std::string_view             name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  std::apply([&]<typename... T>(const T&... members) { (xml_io_stream<T>::write(os, members, binary), ...); },
             oem_xml_members(value));
  tag.write_to_end_stream(os);
}
void xml_io_stream<OptimalEstimationData>::read(std::istream& is, OptimalEstimationData& value, bifstream* binary) {
  value.require_unchecked("read XML");
  OptimalEstimationData next;
  XMLTag                tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  std::apply([&]<typename... T>(T&... members) { (xml_io_stream<T>::read(is, members, binary), ...); },
             oem_xml_members(next));
  tag.read_from_stream(is);
  tag.check_end_name(type_name);
  value = std::move(next);
}
