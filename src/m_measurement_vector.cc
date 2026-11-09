#include <jac_polyfit.h>
#include <workspace.h>

#include <memory>

void measurement_vec_errorFromModelState(Vector&                measurement_vec_error,
                                         Matrix&                measurement_jac_error,
                                         const ArrayOfSensor&   measurement_sensor,
                                         const JacobianTargets& model_state_targets,
                                         const JacobianTargets& jac_targets,
                                         const Vector&          model_state_vec) try {
  ARTS_TIME_REPORT

  const Size ny = measurement_size(measurement_sensor);

  measurement_vec_error.resize(ny);
  measurement_vec_error = 0.0;

  measurement_jac_error.resize(ny, jac_targets.x_size());
  measurement_jac_error = 0.0;

  for (const auto& elem : model_state_targets.error) { elem.update_model(measurement_vec_error, model_state_vec); }
  for (const auto& elem : jac_targets.error) {
    elem.update_jac(measurement_jac_error, model_state_vec, measurement_vec_error);
  }
}
ARTS_METHOD_ERROR_CATCH

void measurement_vecAddError(Vector&       measurement_vec,
                             Matrix&       measurement_jac,
                             const Vector& measurement_vec_error,
                             const Matrix& measurement_jac_error) try {
  ARTS_TIME_REPORT

  if (measurement_vec.shape() != measurement_vec_error.shape())
    arts::user_error(R"(Mismatched shapes:

measurement_vec.shape()       : {:B,}
measurement_vec_error.shape() : {:B,}
)",
                     measurement_vec.shape(),
                     measurement_vec_error.shape());

  measurement_vec += measurement_vec_error;

  if (not measurement_jac.empty()) {
    if (measurement_jac.shape() != measurement_jac_error.shape())
      arts::user_error(R"(Mismatched shapes:

measurement_jac.shape()       : {:B,}
measurement_jac_error.shape() : {:B,}
)",
                       measurement_jac.shape(),
                       measurement_jac_error.shape());

    measurement_jac += measurement_jac_error;
  }
}
ARTS_METHOD_ERROR_CATCH

void jac_targetsAddErrorPolyFit(JacobianTargets&     jac_targets,
                                const ArrayOfSensor& measurement_sensor,
                                const Vector&        t,
                                const Index&         sensor_elem,
                                const Index&         polyorder) {
  ARTS_TIME_REPORT

  if (polyorder < 0) arts::user_error("Polyorder must be non-negative: {}", polyorder);

  if (sensor_elem < 0 or measurement_sensor.size() <= static_cast<Size>(sensor_elem))
    arts::user_error("There are only {0} sensors.  Cannot select index: {1}, please choose an index less than {0}.",
                     measurement_sensor.size(),
                     sensor_elem);

  const Sensor& sensor  = measurement_sensor[sensor_elem];
  const Size    y_start = measurement_start(measurement_sensor, static_cast<Size>(sensor_elem));
  const Size    y_size  = sensor.elements.size();

  if (y_size <= static_cast<Size>(polyorder))
    arts::user_error(R"(Sensor "{}": not enough data points for the given polyorder.

Expected y-size: {}
Poly-order:      {}
)",
                     sensor.name(),
                     y_size,
                     polyorder);

  if (y_size != static_cast<Size>(t.size()))
    arts::user_error(R"(Sensor "{}": not enough grid points for the given data size.
Expected y-size: {}
Grid (size: {}): {:B,}
)",
                     sensor.name(),
                     y_size,
                     t.size(),
                     t);

  make_polyfit(
      jac_targets.emplace_back(ErrorKey{.y_start = y_start, .y_size = y_size}), static_cast<Size>(polyorder), t);
}
