#include <arts_omp.h>
#include <matpack.h>
#include <obsel.h>
#include <rtepack.h>

#include <numeric>
#include <stdexcept>

#include "debug.h"
#include "matpack_data.h"

void measurement_sensorSimple(ArrayOfSensorObsel& measurement_sensor,
                              const AscendingGrid& frequency_grid,
                              const Vector3& pos,
                              const Vector2& los,
                              const Stokvec& pol) try {
  measurement_sensor.resize(frequency_grid.size());

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  for (Size i = 0; i < measurement_sensor.size(); i++) {
    StokvecMatrix w(1, frequency_grid.size(), 0);

    w(0, i)               = pol;
    measurement_sensor[i] = SensorObsel(f, p, std::move(w));
  }
}
ARTS_METHOD_ERROR_CATCH

Numeric gauss(Numeric f0, Numeric f, Numeric fwhm) {
  return std::exp(-4 * Constant::ln_2 * Math::pow2((f - f0) / fwhm));
}

void measurement_sensorSimpleGaussian(ArrayOfSensorObsel& measurement_sensor,
                                      const AscendingGrid& frequency_grid,
                                      const Vector& fwhm,
                                      const Vector3& pos,
                                      const Vector2& los,
                                      const Stokvec& pol) try {
  const Index n = frequency_grid.size();
  ARTS_USER_ERROR_IF(n < 2, "Must have a frequency grid")
  ARTS_USER_ERROR_IF(n != fwhm.size(), "Must have a frequency grid")

  measurement_sensor.resize(n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < measurement_sensor.size(); i++) {
    try {
      StokvecMatrix w(1, frequency_grid.size(), pol);

      for (Index j = 0; j < n; j++) {
        w(0, j) *= gauss(frequency_grid[i], frequency_grid[j], fwhm[i]);
      }

      measurement_sensor[i] = SensorObsel(f, p, std::move(w));
      measurement_sensor[i].normalize(pol, sum(pol));
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "{}", error)
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorSimpleGaussian(ArrayOfSensorObsel& measurement_sensor,
                                      const AscendingGrid& frequency_grid,
                                      const Numeric& fwhm,
                                      const Vector3& pos,
                                      const Vector2& los,
                                      const Stokvec& pol) {
  measurement_sensorSimpleGaussian(measurement_sensor,
                                   frequency_grid,
                                   Vector(frequency_grid.size(), fwhm),
                                   pos,
                                   los,
                                   pol);
}
