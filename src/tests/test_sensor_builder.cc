#include <arts_conversions.h>
#include <sensor_builder.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>

#include "arts_constants.h"

namespace {
void assert_close(Numeric actual,
                  Numeric expected,
                  Numeric tol,
                  std::string_view name) {
  ARTS_USER_ERROR_IF(std::abs(actual - expected) > tol,
                     "{} mismatch: actual {} expected {}",
                     name,
                     actual,
                     expected)
}

void assert_stokvec(Stokvec actual,
                    Stokvec expected,
                    Numeric tol,
                    std::string_view name) {
  for (Size i = 0; i < 4; ++i) {
    ARTS_USER_ERROR_IF(std::abs(actual[i] - expected[i]) > tol,
                       "{} mismatch at {}: actual {} expected {}",
                       name,
                       i,
                       actual,
                       expected)
  }
}

Numeric gaussian_airy_expected_gain(Numeric zenith_deg,
                                    Numeric frequency,
                                    Numeric aperture_diameter) {
  constexpr Numeric gaussian_airy_hwhm_factor =
      3.8317059702075123156 / Constant::pi;
  const Numeric wavelength = Constant::speed_of_light / frequency;
  const Numeric hwhm_deg = Conversion::rad2deg(gaussian_airy_hwhm_factor *
                                               wavelength / aperture_diameter);
  const Numeric ratio    = zenith_deg / hwhm_deg;

  return std::exp(-Constant::ln_2 * ratio * ratio);
}

void test_sensor_builder_returns_meta_per_geometry() {
  sensor::SensorBuilder builder(
      {sensor::BoxChannel{AscendingGrid{100.0, 101.0}},
       sensor::DiracChannel{200.0}},
      sensor::PencilBeamAntenna{});

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};

  const auto [obsels, meta] = builder(pos, los);

  ARTS_USER_ERROR_IF(
      obsels.size() != 4, "Expected 4 obsels, got {}", obsels.size())
  ARTS_USER_ERROR_IF(
      meta.size() != 2, "Expected 2 meta entries, got {}", meta.size())
  ARTS_USER_ERROR_IF(meta[0].count() != 2 or meta[1].count() != 2,
                     "Each meta block must describe one 2-channel geometry")

  const auto& gf0 = std::get<SortedGriddedField1>(meta[0].data);
  const auto& gf1 = std::get<SortedGriddedField1>(meta[1].data);

  ARTS_USER_ERROR_IF(
      gf0.gridname<0>() != "channel" or gf1.gridname<0>() != "channel",
      "SensorBuilder meta grid must be the channel axis")
  assert_close(gf0.grid<0>()[0], 0.0, 0.0, "meta[0] channel index 0");
  assert_close(gf0.grid<0>()[1], 1.0, 0.0, "meta[0] channel index 1");

  ARTS_USER_ERROR_IF(obsels[0].poslos_grid()[0].pos != pos[0],
                     "First geometry position mismatch")
  ARTS_USER_ERROR_IF(obsels[2].poslos_grid()[0].pos != pos[1],
                     "Second geometry position mismatch")
  ARTS_USER_ERROR_IF(not obsels[0].same_poslos(obsels[1]),
                     "Obsels from the same geometry must share poslos")
  ARTS_USER_ERROR_IF(not obsels[2].same_poslos(obsels[3]),
                     "Obsels from the same geometry must share poslos")
  ARTS_USER_ERROR_IF(not obsels[0].same_freqs(obsels[2]),
                     "Obsels for the same channel must share frequencies")
  ARTS_USER_ERROR_IF(not obsels[1].same_freqs(obsels[3]),
                     "Obsels for the same channel must share frequencies")
}

void test_sensor_builder_rejects_mismatched_geometry_counts() {
  sensor::SensorBuilder builder({sensor::DiracChannel{}},
                                sensor::PencilBeamAntenna{});

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};

  bool threw = false;
  try {
    static_cast<void>(builder(pos, los));
  } catch (const std::runtime_error&) {
    threw = true;
  }

  ARTS_USER_ERROR_IF(
      not threw,
      "SensorBuilder must reject mismatching position and LOS counts")
}

void test_sensor_builder_uses_gaussian_airy_frequency_dependence() {
  const Stokvec peak_weight{2.0, 0.0, 0.0, 0.0};
  sensor::SensorBuilder builder(
      {sensor::BoxChannel{AscendingGrid{100.0e9, 200.0e9}}},
      sensor::GaussianAiryAntenna{
          ZenGrid{{0.0, 0.2}}, AziGrid{{0.0}}, 1.0, peak_weight});

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 1> los{{{45.0, 30.0}}};

  const auto [obsels, meta] = builder(pos, los);

  ARTS_USER_ERROR_IF(
      obsels.size() != 1, "Expected 1 obsel, got {}", obsels.size())
  ARTS_USER_ERROR_IF(
      meta.size() != 1, "Expected 1 meta entry, got {}", meta.size())

  const Numeric low_gain  = gaussian_airy_expected_gain(0.2, 100.0e9, 1.0);
  const Numeric high_gain = gaussian_airy_expected_gain(0.2, 200.0e9, 1.0);

  assert_stokvec(obsels[0].weight_matrix()[1, 0],
                 0.5 * low_gain * peak_weight,
                 1e-12,
                 "builder gaussian airy off-axis low frequency");
  assert_stokvec(obsels[0].weight_matrix()[1, 1],
                 0.5 * high_gain * peak_weight,
                 1e-12,
                 "builder gaussian airy off-axis high frequency");
}
}  // namespace

int main() {
  test_sensor_builder_returns_meta_per_geometry();
  test_sensor_builder_rejects_mismatched_geometry_counts();
  test_sensor_builder_uses_gaussian_airy_frequency_dependence();
  return 0;
}