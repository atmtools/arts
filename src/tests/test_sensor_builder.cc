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
void assert_close(Numeric actual, Numeric expected, Numeric tol, std::string_view name) {
  if (std::abs(actual - expected) > tol) arts::user_error("{} mismatch: actual {} expected {}", name, actual, expected);
}

void test_sensor_builder_returns_meta_per_geometry() try {
  sensor::Builder builder({sensor::BoxChannel{AscendingGrid{100.0, 101.0}}, sensor::DiracChannel{200.0}},
                          std::make_shared<const sensor::PencilBeamAntenna>());

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const std::array<Time, 2>    time{{Time{"2001-01-01 00:00:00"}, Time{"2001-01-01 00:01:00"}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const auto [obsels, meta] = builder(pos, los, time, ell);

  if (obsels.size() != 4) arts::user_error("Expected 4 obsels, got {}", obsels.size());
  if (meta.size() != 2) arts::user_error("Expected 2 meta entries, got {}", meta.size());
  if (meta[0].count() != 2 or meta[1].count() != 2)
    arts::user_error("Each meta block must describe one 2-channel geometry");
  if (not meta[0].time or *meta[0].time != time[0]) arts::user_error("First metadata time mismatch");
  if (not meta[1].time or *meta[1].time != time[1]) arts::user_error("Second metadata time mismatch");

  const auto& gf0 = std::get<SortedGriddedField1>(meta[0].data);
  const auto& gf1 = std::get<SortedGriddedField1>(meta[1].data);

  if (gf0.gridname<0>() != "channel" or gf1.gridname<0>() != "channel")
    arts::user_error("Builder meta grid must be the channel axis");
  assert_close(gf0.grid<0>()[0], 0.0, 0.0, "meta[0] channel index 0");
  assert_close(gf0.grid<0>()[1], 1.0, 0.0, "meta[0] channel index 1");

  if (obsels[0].poslos_grid()[0].pos != pos[0]) arts::user_error("First geometry position mismatch");
  if (obsels[2].poslos_grid()[0].pos != pos[1]) arts::user_error("Second geometry position mismatch");
  if (not obsels[0].same_poslos(obsels[1])) arts::user_error("Obsels from the same geometry must share poslos");
  if (not obsels[2].same_poslos(obsels[3])) arts::user_error("Obsels from the same geometry must share poslos");
  if (not obsels[0].same_freqs(obsels[2])) arts::user_error("Obsels for the same channel must share frequencies");
  if (not obsels[1].same_freqs(obsels[3])) arts::user_error("Obsels for the same channel must share frequencies");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_sensor_builder_returns_meta_per_geometry failed:\n{}", e.what()));
}

void test_sensor_builder_rejects_mismatched_geometry_counts() try {
  sensor::Builder builder({sensor::DiracChannel{}}, std::make_shared<const sensor::PencilBeamAntenna>());

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  bool threw = false;
  try {
    static_cast<void>(builder(pos, los, ell));
  } catch (const std::runtime_error&) { threw = true; }

  if (not threw) arts::user_error("Builder must reject mismatching position and LOS counts");

  const std::array<Vector2, 1> matching_los{{{20.0, 30.0}}};
  const std::array<Time, 2>    mismatched_time{{Time{}, Time{}}};
  threw = false;
  try {
    static_cast<void>(builder(pos, matching_los, mismatched_time, ell));
  } catch (const std::runtime_error&) { threw = true; }
  if (not threw) arts::user_error("Builder must reject mismatching geometry and time counts");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_sensor_builder_rejects_mismatched_geometry_counts failed:\n{}", e.what()));
}

void test_unflatten_updates_shared_poslos_grids() try {
  auto freq_grid = std::make_shared<const AscendingGrid>(AscendingGrid{100.0});

  SensorPosLosVector poslos_grid(1);
  poslos_grid[0]     = {.pos = {600e3, 10.0, 20.0}, .los = {30.0, 40.0}};
  auto shared_poslos = std::make_shared<const SensorPosLosVector>(std::move(poslos_grid));

  sensor::SparseStokvecMatrix weights(1, 1);
  weights[0, 0] = {1.0, 0.0, 0.0, 0.0};

  ArrayOfSensorObsel obsels(2);
  obsels[0] = {freq_grid, shared_poslos, weights};
  obsels[1] = {freq_grid, shared_poslos, weights};

  const auto   original_poslos = obsels[0].poslos_grid_ptr();
  const Vector zeniths{15.0};

  unflatten(obsels, zeniths, obsels[0], SensorKeyType::zen);

  if (obsels[0].poslos_grid_ptr() == original_poslos) arts::user_error("unflatten must replace the shared poslos grid");
  if (not obsels[0].same_poslos(obsels[1]))
    arts::user_error("obsels sharing a geometry must still share it after unflatten");
  assert_close(obsels[0].poslos_grid()[0].los[0], 15.0, 0.0, "updated zenith for first obsel");
  assert_close(obsels[1].poslos_grid()[0].los[0], 15.0, 0.0, "updated zenith for second obsel");
  assert_close(obsels[0].poslos_grid()[0].los[1], 40.0, 0.0, "azimuth must remain unchanged");

  const auto   original_freq = obsels[0].f_grid_ptr();
  const Vector freqs{101.0};

  unflatten(obsels, freqs, obsels[0], SensorKeyType::freq);

  if (obsels[0].f_grid_ptr() == original_freq) arts::user_error("unflatten must replace the shared frequency grid");
  if (not obsels[0].same_freqs(obsels[1]))
    arts::user_error("obsels sharing frequencies must still share them after unflatten");
  assert_close(obsels[0].f_grid()[0], 101.0, 0.0, "updated frequency for first obsel");
  assert_close(obsels[1].f_grid()[0], 101.0, 0.0, "updated frequency for second obsel");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_unflatten_updates_shared_poslos_grids failed:\n{}", e.what()));
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
