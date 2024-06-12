#include <matpack.h>
#include <obsel.h>
#include <rtepack.h>

#include <numeric>

#include "debug.h"

//! Ensures sum of w is 1.0
void normalize(Vector& w) {
  std::transform(w.begin(),
                 w.end(),
                 w.begin(),
                 [w_sum = std::reduce(w.begin(), w.end(), 0.0)](
                     const auto& wx) { return wx / w_sum; });
}

void gaussian(Vector& w,
              const AscendingGrid& frequency_grid,
              const Numeric& f0,
              const Numeric& fwhm) {
  std::transform(frequency_grid.begin(),
                 frequency_grid.end(),
                 w.begin(),
                 [f0, fwhm](const auto& f) {
                   return std::exp(-0.5 * Math::pow2((f - f0) / fwhm));
                 });
}

void measurement_vector_sensorSimple(
    ArrayOfSensorObsel& measurement_vector_sensor,
    const AscendingGrid& frequency_grid,
    const Vector3& pos,
    const Vector2& los,
    const Stokvec& pol) try {
  measurement_vector_sensor.resize(frequency_grid.size());

  std::transform(
      frequency_grid.begin(),
      frequency_grid.end(),
      measurement_vector_sensor.begin(),
      [pos, los, pol](const auto& f) {
        return SensorObsel{
            .f_grid_w      = Vector{1},
            .f_grid        = {f},
            .poslos_grid_w = MuelmatVector{1},
            .poslos_grid   = SensorPosLosVector{{.pos = pos, .los = los}},
            .polarization  = pol};
      });
}
ARTS_METHOD_ERROR_CATCH

void measurement_vector_sensorGaussianFrequencyGrid(
    ArrayOfSensorObsel& measurement_vector_sensor,
    const AscendingGrid& frequency_grid,
    const ArrayOfVector2& f0_fwmh,
    const Numeric& weight_cutoff,
    const Vector3& pos,
    const Vector2& los,
    const Stokvec& pol) try {
  const Index n = frequency_grid.size();
  ARTS_USER_ERROR_IF(n < 2, "Must have a frequency grid")

  measurement_vector_sensor.resize(0);
  measurement_vector_sensor.reserve(f0_fwmh.size());

  for (const auto [f0, fwhm] : f0_fwmh) {
    auto& obsel = measurement_vector_sensor.emplace_back(SensorObsel{
        .f_grid_w      = Vector(n, 1.0 / static_cast<Numeric>(n)),
        .f_grid        = frequency_grid,
        .poslos_grid_w = MuelmatVector{1},
        .poslos_grid   = SensorPosLosVector{{.pos = pos, .los = los}},
        .polarization  = pol,
    });
    ARTS_USER_ERROR_IF(f0 <= 0.0, "Must have a positive frequency")
    ARTS_USER_ERROR_IF(fwhm < 0.0, "Must have a non-negative FWHM")

    if (fwhm == 0.0) {
      auto ptr       = std::ranges::find(frequency_grid, f0);
      obsel.f_grid   = AscendingGrid{f0};
      obsel.f_grid_w = {1.0};
      ARTS_USER_ERROR_IF(
          ptr == frequency_grid.end(),
          "Frequency not found in the grid, must be exact to give dirac-channels")
      continue;
    }

    gaussian(obsel.f_grid_w, obsel.f_grid, f0, fwhm);

    if (weight_cutoff > 0.0) {
      //! Remove weights above and below the cutoff
      while (obsel.f_grid_w.size() > 0 and
             obsel.f_grid_w.front() < weight_cutoff) {
        obsel.f_grid_w.erase(obsel.f_grid_w.begin());
        obsel.f_grid.erase(obsel.f_grid.begin());
      }
      while (obsel.f_grid_w.size() > 0 and
             obsel.f_grid_w.back() < weight_cutoff) {
        obsel.f_grid_w.erase(obsel.f_grid_w.end() - 1);
        obsel.f_grid.erase(obsel.f_grid.end() - 1);
      }
    }

    normalize(obsel.f_grid_w);
  }
}
ARTS_METHOD_ERROR_CATCH

void measurement_vector_sensorGaussian(
    ArrayOfSensorObsel& measurement_vector_sensor,
    const ArrayOfVector3& f0_fwmh_df,
    const Numeric& weight_cutoff,
    const Vector3& pos,
    const Vector2& los,
    const Stokvec& pol) try {
  ARTS_USER_ERROR_IF(weight_cutoff >= 1.0 or weight_cutoff <= 0.0,
                     "Weight cutoff must be in (0,1]")

  const Numeric wc = std::sqrt(-2.0 * std::log(weight_cutoff));

  measurement_vector_sensor.resize(0);
  measurement_vector_sensor.reserve(f0_fwmh_df.size());

  for (auto& guass_params : f0_fwmh_df) {
    const Numeric f0 = guass_params[0], fwhm = guass_params[1],
                  df = guass_params[2];
    ARTS_USER_ERROR_IF(f0 <= 0.0, "Must have a positive frequency")
    ARTS_USER_ERROR_IF(fwhm <= 0.0, "Must have a positive FWHM")
    ARTS_USER_ERROR_IF(df <= 0.0, "Must have a positive frequency step")

    const auto n = 1 + static_cast<Index>(wc * fwhm / df);
    Vector f_grid(2 * n + 1);
    for (Index i = -n; i <= n; i++) {
      f_grid[i + n] = f0 + static_cast<Numeric>(i) * df;
    }

    auto& obsel = measurement_vector_sensor.emplace_back(SensorObsel{
        .f_grid_w      = Vector(2 * n + 1),
        .f_grid        = std::move(f_grid),
        .poslos_grid_w = MuelmatVector{1},
        .poslos_grid   = SensorPosLosVector{{.pos = pos, .los = los}},
        .polarization  = pol,
    });

    gaussian(obsel.f_grid_w, obsel.f_grid, f0, fwhm);
    normalize(obsel.f_grid_w);
  }
}
ARTS_METHOD_ERROR_CATCH
