#include <arts_omp.h>
#include <matpack.h>
#include <obsel.h>
#include <rtepack.h>

#include <numeric>
#include <stdexcept>

#include "debug.h"
#include "matpack_data.h"

void measurement_sensorInit(ArrayOfSensorObsel& measurement_sensor) {
  measurement_sensor = {};
}

void measurement_sensorAddSimple(ArrayOfSensorObsel& measurement_sensor,
                                 const AscendingGrid& frequency_grid,
                                 const Vector3& pos,
                                 const Vector2& los,
                                 const Stokvec& pol) try {
  const Index n = frequency_grid.size();
  const Size sz = measurement_sensor.size();

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  for (Index i = 0; i < n; i++) {
    StokvecMatrix w(1, n, 0);
    w(0, i)                    = pol;
    measurement_sensor[i + sz] = {f, p, std::move(w)};
  }
}
ARTS_METHOD_ERROR_CATCH

Numeric gauss(Numeric f0, Numeric f, Numeric fwhm) {
  return std::exp(-4 * Constant::ln_2 * Math::pow2((f - f0) / fwhm));
}

void measurement_sensorAddVectorGaussian(ArrayOfSensorObsel& measurement_sensor,
                                         const AscendingGrid& frequency_grid,
                                         const Vector& fwhm,
                                         const Vector3& pos,
                                         const Vector2& los,
                                         const Stokvec& pol) try {
  const Index n = frequency_grid.size();
  const Size sz = measurement_sensor.size();

  ARTS_USER_ERROR_IF(n < 2, "Must have a frequency grid")
  ARTS_USER_ERROR_IF(n != fwhm.size(), "Must have a frequency grid")

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index i = 0; i < n; i++) {
    try {
      StokvecMatrix w(1, n, pol);

      for (Index j = 0; j < n; j++) {
        w(0, j) *= gauss(frequency_grid[i], frequency_grid[j], fwhm[i]);
      }

      measurement_sensor[i + sz] = {f, p, std::move(w)};
      measurement_sensor[i + sz].normalize(pol, sum(pol));
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
                                         const Numeric& fwhm,
                                         const Vector3& pos,
                                         const Vector2& los,
                                         const Stokvec& pol) {
  measurement_sensorAddVectorGaussian(measurement_sensor,
                                      frequency_grid,
                                      Vector(frequency_grid.size(), fwhm),
                                      pos,
                                      los,
                                      pol);
}
