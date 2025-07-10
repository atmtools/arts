#include <workspace.h>

#include <limits>
#include <memory>

#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_helpers_grid_t.h"

void measurement_vector_errorFromModelState(
    Vector& measurement_vector_error,
    Matrix& measurement_jacobian_error,
    const ArrayOfSensorObsel& measurement_sensor,
    const JacobianTargets& jacobian_targets,
    const Vector& model_state_vector) {
  ARTS_TIME_REPORT

  measurement_vector_error.resize(measurement_sensor.size());
  measurement_vector_error = 0.0;

  measurement_jacobian_error.resize(measurement_sensor.size(),
                                    jacobian_targets.x_size());
  measurement_jacobian_error = 0.0;

  for (auto& elem : jacobian_targets.error()) {
    elem.update_model(measurement_vector_error, model_state_vector);
    elem.update_jac(measurement_jacobian_error, model_state_vector);
  }
}

void measurement_vectorConditionalAddError(
    Vector& measurement_vector,
    Matrix& measurement_jacobian,
    const Vector& measurement_vector_error,
    const Matrix& measurement_jacobian_error,
    const Index& do_jacobian) {
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

  if (do_jacobian != 0) {
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

void model_state_vectorUpdateError(Vector& model_state_vector,
                                   const JacobianTargets& jacobian_targets,
                                   const Vector& measurement_vector,
                                   const Vector& measurement_vector_fitted) {
  ARTS_TIME_REPORT

  for (auto& elem : jacobian_targets.error()) {
    const Range r(elem.type.y_size, elem.type.y_size);
    Vector meas{measurement_vector[r]};
    meas -= measurement_vector_fitted[r];
    elem.update_state(model_state_vector, meas);
  }
}

struct polyfit {
  std::shared_ptr<Vector> st;

  Vector operator()(Vector y, const ConstVectorView p) const {
    ARTS_USER_ERROR_IF(not st, "No t-vector provided for polyinv.")
    const Vector& t = *st;

    ARTS_USER_ERROR_IF(y.size() != t.size(), "Mismatched y and t sizes.")

    for (Size j = 0; j < t.size(); j++) {
      const Numeric tj = t[j];
      Numeric xn       = 1.0;
      for (Size i = 0; i < p.size(); i++) {
        y[j] += p[i] * xn;
        xn   *= tj;
      }
    }

    return y;
  }

  Matrix operator()(Matrix dy, const ConstVectorView p) const {
    ARTS_USER_ERROR_IF(not st, "No t-vector provided for polyinv.")
    const Vector& t = *st;

    ARTS_USER_ERROR_IF(static_cast<Index>(t.size()) != dy.nrows(),
                       "Mismatched y and dy sizes.")
    ARTS_USER_ERROR_IF(dy.ncols() != static_cast<Index>(p.size()),
                       "Mismatched dy and p sizes.")

    for (Size j = 0; j < t.size(); j++) {
      const Numeric tj = t[j];
      Numeric xn       = 1.0;
      for (Size i = 0; i < p.size(); i++) {
        dy[j, i]  = xn;
        xn       *= tj;
      }
    }

    return dy;
  }
};

struct polyinv {
  std::shared_ptr<Vector> st;

  Vector operator()(Vector p, const ConstVectorView y) const {
    ARTS_USER_ERROR_IF(not st, "No t-vector provided for polyfit.")
    const Vector& t = *st;

    ARTS_USER_ERROR_IF(y.size() != t.size(), "Mismatched y and t sizes.")

    Jacobian::polyfit(p, t, y);

    return p;
  }
};

void jacobian_targetsAddErrorPolyFit(
    JacobianTargets& jacobian_targets,
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

  auto st = std::make_shared<Vector>(t);
  const polyfit pf{st};
  const polyinv pi{st};

  auto& x = jacobian_targets.emplace_back(
      ErrorKey{.y_start = y_start, .y_size = y_size}, {});
  x.x_size             = polyorder + 1;
  x.inverse_state      = pi;
  x.transform_state    = pf;
  x.transform_jacobian = pf;
}
