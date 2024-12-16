#include <jacobian.h>
#include <matpack.h>

#include <limits>

void measurement_vector_errorInitStandard(Vector& measurement_vector_error,
                                          Matrix& measuremen_jacobian_error,
                                          const Vector& measurement_vector,
                                          const Matrix& measuremen_jacobian) {
  measurement_vector_error.resize(measurement_vector.shape());
  measuremen_jacobian_error.resize(measuremen_jacobian.shape());
  measurement_vector_error  = 0.0;
  measuremen_jacobian_error = 0.0;
}

void measurement_vector_errorAddErrorState(
    Vector& measurement_vector_error,
    Matrix& measuremen_jacobian_error,
    const JacobianTargets& jacobian_targets,
    const Vector& model_state_vector) {
  for (auto& elem : jacobian_targets.error()) {
    elem.update_y(measurement_vector_error,
                  measuremen_jacobian_error,
                  model_state_vector);
  }
}

void measurement_vectorAddError(Vector& measurement_vector,
                                Matrix& measuremen_jacobian,
                                const Vector& measurement_vector_error,
                                const Matrix& measuremen_jacobian_error) {
  ARTS_USER_ERROR_IF(
      measurement_vector.shape() != measurement_vector_error.shape(),
      R"(Mismatched shapes:

measurement_vector.shape()       : {:B,}
measurement_vector_error.shape() : {:B,}
)",
      measurement_vector.shape(),
      measurement_vector_error.shape())
  ARTS_USER_ERROR_IF(
      measuremen_jacobian.shape() != measuremen_jacobian_error.shape(),
      R"(Mismatched shapes:

measuremen_jacobian.shape()       : {:B,}
measuremen_jacobian_error.shape() : {:B,}
)",
      measuremen_jacobian.shape(),
      measuremen_jacobian_error.shape())

  measurement_vector  += measurement_vector_error;
  measuremen_jacobian += measuremen_jacobian_error;
}

void model_state_vectorUpdateError(Vector& model_state_vector,
                                   const JacobianTargets& jacobian_targets,
                                   const Vector& measurement_vector,
                                   const Vector& measurement_vector_fitted) {
  for (auto& elem : jacobian_targets.error()) {
    elem.update_x(
        model_state_vector, measurement_vector, measurement_vector_fitted);
  }
}

struct polyfit {
  Vector t;

  void operator()(ExhaustiveVectorView y,
                  MatrixView dy,
                  const ExhaustiveConstVectorView p) const {
    ARTS_USER_ERROR_IF(y.size() != t.size(), "Mismatched y and t sizes.")
    ARTS_USER_ERROR_IF(y.size() != dy.nrows(), "Mismatched y and dy sizes.")
    ARTS_USER_ERROR_IF(dy.ncols() != p.size(), "Mismatched dy and p sizes.")

    for (Index j = 0; j < t.size(); j++) {
      for (Index i = 0; i < p.size(); i++) {
        const Numeric xn  = std::pow(t[j], i);
        y[j]             += p[i] * xn;
        dy[j, i]          = xn;
      }
    }
  }

  void operator()(ExhaustiveVectorView p,
                  const ExhaustiveConstVectorView y) const {
    ARTS_USER_ERROR_IF(y.size() != t.size(), "Mismatched y and t sizes.")
    Jacobian::polyfit(p, t, y);
  }
};

void jacobian_targetsAddErrorPolyFit(
    JacobianTargets& jacobian_targets,
    const ArrayOfSensorObsel& measurement_sensor,
    const Vector& t,
    const Index& sensor_elem,
    const Index& polyorder) {
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

  const polyfit p{t};

  jacobian_targets.emplace_back(ErrorKey{.y_start = y_start, .y_size = y_size},
                                static_cast<Size>(polyorder + 1),
                                p);
}
