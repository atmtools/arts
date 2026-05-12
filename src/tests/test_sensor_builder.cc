#include <sensor_builder.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>

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

void test_sensor_builder_returns_meta_per_geometry() {
  sensor::SensorBuilder builder{
      .channels = {sensor::BoxChannel{AscendingGrid{100.0, 101.0}},
                   sensor::DiracChannel{200.0}},
      .antenna  = sensor::PencilBeamAntenna{},
  };

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
  sensor::SensorBuilder builder{
      .channels = {sensor::DiracChannel{}},
      .antenna  = sensor::PencilBeamAntenna{},
  };

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
}  // namespace

int main() {
  test_sensor_builder_returns_meta_per_geometry();
  test_sensor_builder_rejects_mismatched_geometry_counts();
  return 0;
}