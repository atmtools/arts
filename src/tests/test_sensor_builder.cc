#include <arts_constants.h>
#include <arts_conversions.h>
#include <obsel.h>
#include <planet_data.h>
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

void test_sensor_builder_returns_meta_per_geometry() try {
  sensor::Builder builder({sensor::BoxChannel{AscendingGrid{100.0, 101.0}},
                           sensor::DiracChannel{200.0}},
                          std::make_shared<const sensor::PencilBeamAntenna>());

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const Vector2 ell{Body::Earth::a, Body::Earth::b};

  const auto [obsels, meta] = builder(pos, los, ell);

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
      "Builder meta grid must be the channel axis")
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
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_sensor_builder_returns_meta_per_geometry failed:\n{}", e.what()));
}

void test_sensor_builder_rejects_mismatched_geometry_counts() try {
  sensor::Builder builder({sensor::DiracChannel{}},
                          std::make_shared<const sensor::PencilBeamAntenna>());

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const Vector2 ell{Body::Earth::a, Body::Earth::b};

  bool threw = false;
  try {
    static_cast<void>(builder(pos, los, ell));
  } catch (const std::runtime_error&) {
    threw = true;
  }

  ARTS_USER_ERROR_IF(not threw,
                     "Builder must reject mismatching position and LOS counts")
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_sensor_builder_rejects_mismatched_geometry_counts failed:\n{}", e.what()));
}

void test_unflatten_updates_shared_poslos_grids() try {
  auto freq_grid = std::make_shared<const AscendingGrid>(AscendingGrid{100.0});

  SensorPosLosVector poslos_grid(1);
  poslos_grid[0] = {.pos = {600e3, 10.0, 20.0}, .los = {30.0, 40.0}};
  auto shared_poslos =
      std::make_shared<const SensorPosLosVector>(std::move(poslos_grid));

  sensor::SparseStokvecMatrix weights(1, 1);
  weights[0, 0] = {1.0, 0.0, 0.0, 0.0};

  ArrayOfSensorObsel obsels(2);
  obsels[0] = {freq_grid, shared_poslos, weights};
  obsels[1] = {freq_grid, shared_poslos, weights};

  const auto original_poslos = obsels[0].poslos_grid_ptr();
  const Vector zeniths{15.0};

  unflatten(obsels, zeniths, obsels[0], SensorKeyType::zen);

  ARTS_USER_ERROR_IF(obsels[0].poslos_grid_ptr() == original_poslos,
                     "unflatten must replace the shared poslos grid")
  ARTS_USER_ERROR_IF(
      not obsels[0].same_poslos(obsels[1]),
      "obsels sharing a geometry must still share it after unflatten")
  assert_close(obsels[0].poslos_grid()[0].los[0],
               15.0,
               0.0,
               "updated zenith for first obsel");
  assert_close(obsels[1].poslos_grid()[0].los[0],
               15.0,
               0.0,
               "updated zenith for second obsel");
  assert_close(obsels[0].poslos_grid()[0].los[1],
               40.0,
               0.0,
               "azimuth must remain unchanged");

  const auto original_freq = obsels[0].f_grid_ptr();
  const Vector freqs{101.0};

  unflatten(obsels, freqs, obsels[0], SensorKeyType::freq);

  ARTS_USER_ERROR_IF(obsels[0].f_grid_ptr() == original_freq,
                     "unflatten must replace the shared frequency grid")
  ARTS_USER_ERROR_IF(
      not obsels[0].same_freqs(obsels[1]),
      "obsels sharing frequencies must still share them after unflatten")
  assert_close(
      obsels[0].f_grid()[0], 101.0, 0.0, "updated frequency for first obsel");
  assert_close(
      obsels[1].f_grid()[0], 101.0, 0.0, "updated frequency for second obsel");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "test_unflatten_updates_shared_poslos_grids failed:\n{}", e.what()));
}
}  // namespace

int main() try {
  test_sensor_builder_returns_meta_per_geometry();
  test_sensor_builder_rejects_mismatched_geometry_counts();
  test_unflatten_updates_shared_poslos_grids();
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "test_sensor_builder failed:\n{}", e.what());
  return EXIT_FAILURE;
}
