#include <arts_omp.h>
#include <obsel.h>
#include <rtepack.h>
#include <workspace.h>

#include <boost/math/distributions/normal.hpp>
#include <numeric>
#include <stdexcept>

#include "debug.h"

void measurement_sensorInit(ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  measurement_sensor = {};
}

void measurement_sensorAddSimple(ArrayOfSensorObsel& measurement_sensor,
                                 const AscendingGrid& frequency_grid,
                                 const Vector3& pos,
                                 const Vector2& los,
                                 const Stokvec& pol) try {
  ARTS_TIME_REPORT

  const Index n = frequency_grid.size();
  const Size sz = measurement_sensor.size();

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  for (Index i = 0; i < n; i++) {
    StokvecMatrix w(1, n, 0);
    w[0, i]                    = pol;
    measurement_sensor[i + sz] = {f, p, std::move(w)};
  }
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorAddVectorGaussian(ArrayOfSensorObsel& measurement_sensor,
                                         const AscendingGrid& frequency_grid,
                                         const Vector& stds,
                                         const Vector3& pos,
                                         const Vector2& los,
                                         const Stokvec& pol) try {
  ARTS_TIME_REPORT

  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  const Size n       = frequency_grid.size();
  const Size sz      = measurement_sensor.size();
  const Size nonzero = pol.nonzero_components();

  ARTS_USER_ERROR_IF(n < 2, "Must have a frequency grid")
  ARTS_USER_ERROR_IF(n != stds.size(),
                     "Must have a standard deviation for each frequency point")
  ARTS_USER_ERROR_IF(nonzero == 0, "pol is 0")

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n; i++) {
    try {
      StokvecMatrix w(1, n, pol);

      const gauss dist(frequency_grid[i], stds[i]);
      for (Size j = 0; j < n; j++) {
        w[0, j] *= pdf(dist, frequency_grid[j]);
      }

      measurement_sensor[i + sz] = {f, p, std::move(w)};
      measurement_sensor[i + sz].normalize(pol);
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "{}", error)
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorAddSimpleGaussian(ArrayOfSensorObsel& measurement_sensor,
                                         const AscendingGrid& frequency_grid,
                                         const Numeric& std,
                                         const Vector3& pos,
                                         const Vector2& los,
                                         const Stokvec& pol) {
  ARTS_TIME_REPORT

  const Vector stds(frequency_grid.size(), std);
  measurement_sensorAddVectorGaussian(
      measurement_sensor, frequency_grid, stds, pos, los, pol);
}
