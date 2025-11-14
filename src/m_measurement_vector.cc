#include <jac_polyfit.h>
#include <workspace.h>

#include <memory>

void measurement_vector_errorFromModelState(
    Vector& measurement_vector_error,
    Matrix& measurement_jacobian_error,
    const ArrayOfSensorObsel& measurement_sensor,
    const JacobianTargets& jac_targets,
    const Vector& model_state_vector) try {
  ARTS_TIME_REPORT

  measurement_vector_error.resize(measurement_sensor.size());
  measurement_vector_error = 0.0;

  measurement_jacobian_error.resize(measurement_sensor.size(),
                                    jac_targets.x_size());
  measurement_jacobian_error = 0.0;

  for (auto& elem : jac_targets.error) {
    elem.update_model(measurement_vector_error, model_state_vector);
    elem.update_jac(measurement_jacobian_error,
                    model_state_vector,
                    measurement_vector_error);
  }
}
ARTS_METHOD_ERROR_CATCH

void measurement_vectorConditionalAddError(
    Vector& measurement_vector,
    Matrix& measurement_jacobian,
    const Vector& measurement_vector_error,
    const Matrix& measurement_jacobian_error,
    const Index& do_jac) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      measurement_vector.shape() != measurement_vector_error.shape(),
      R"(Mismatched shapes:

measurement_vector.shape()       : {:B,}
measurement_vector_error.shape() : {:B,}
)",
      measurement_vector.shape(),
      measurement_vector_error.shape())

  measurement_vector += measurement_vector_error;

  if (do_jac != 0) {
    ARTS_USER_ERROR_IF(
        measurement_jacobian.shape() != measurement_jacobian_error.shape(),
        R"(Mismatched shapes:

measurement_jacobian.shape()       : {:B,}
measurement_jacobian_error.shape() : {:B,}
)",
        measurement_jacobian.shape(),
        measurement_jacobian_error.shape())

    measurement_jacobian += measurement_jacobian_error;
  }
}
ARTS_METHOD_ERROR_CATCH

void jac_targetsAddErrorPolyFit(JacobianTargets& jac_targets,
                                const ArrayOfSensorObsel& measurement_sensor,
                                const Vector& t,
                                const Index& sensor_elem,
                                const Index& polyorder) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      polyorder < 0, "Polyorder must be non-negative: {}", polyorder)

  std::vector<std::pair<Size, const void*>> sensor_grid_ptrs;
  for (Size i = 0; i < measurement_sensor.size(); i++) {
    const void* ptr =
        reinterpret_cast<const void*>(measurement_sensor[i].f_grid_ptr().get());
    if (std::ranges::none_of(sensor_grid_ptrs | std::views::values,
                             Cmp::eq(ptr))) {
      sensor_grid_ptrs.emplace_back(i, ptr);
    }
  }

  ARTS_USER_ERROR_IF(
      sensor_grid_ptrs.size() <= static_cast<Size>(sensor_elem),
      "There are only {0} independent sensor frequency grids.  "
      "Cannot select index: {1}, please choose an index less than {0}.",
      sensor_grid_ptrs.size(),
      sensor_elem)

  const Size y_start = sensor_grid_ptrs[sensor_elem].first;
  const Size y_size  = [&]() -> Size {
    for (Size i = y_start + 1; i < measurement_sensor.size(); i++) {
      if (not measurement_sensor[y_start].same_freqs(measurement_sensor[i])) {
        return i;
      }
    }

    return measurement_sensor.size();
  }() - y_start;

  ARTS_USER_ERROR_IF(y_size <= static_cast<Size>(polyorder),
                     R"(Not enough data points for the given polyorder.

Expected y-size: {}
Poly-order:      {}
)",
                     y_size,
                     polyorder)

  ARTS_USER_ERROR_IF(y_size != static_cast<Size>(t.size()),
                     R"(Not enough grid points for the given data size.
Expected y-size: {}
Grid (size: {}): {:B,}
)",
                     y_size,
                     t.size(),
                     t)

  make_polyfit(
      jac_targets.emplace_back(ErrorKey{.y_start = y_start, .y_size = y_size}),
      static_cast<Size>(polyorder),
      t);
}
