#include <arts_constants.h>
#include <arts_conversions.h>
#include <frequency_bandpass_filters.h>
#include <frequency_range_selection.h>
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
  sensor::Builder builder{.channels = {sensor::BoxChannel{AscendingGrid{100.0, 101.0}}, sensor::DiracChannel{200.0}}, .antenna = sensor::PencilBeamAntenna{}};

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const auto [obsels, meta] = builder(pos, los, ell);

  if (obsels.size() != 4) arts::user_error("Expected 4 obsels, got {}", obsels.size());
  if (meta.size() != 2) arts::user_error("Expected 2 meta entries, got {}", meta.size());
  if (meta[0].count() != 2 or meta[1].count() != 2)
    arts::user_error("Each meta block must describe one 2-channel geometry");

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
  sensor::Builder builder{.channels = {sensor::DiracChannel{}}, .antenna = sensor::PencilBeamAntenna{}};

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  bool threw = false;
  try {
    static_cast<void>(builder(pos, los, ell));
  } catch (const std::runtime_error&) { threw = true; }

  if (not threw) arts::user_error("Builder must reject mismatching position and LOS counts");
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

void test_sensor_builder_polarization_per_channel() try {
  sensor::Channel vertical   = sensor::DiracChannel{100.0};
  sensor::Channel horizontal = sensor::DiracChannel{200.0};
  vertical.polarization      = {.state = {1.0, 1.0, 0.0, 0.0}};
  horizontal.polarization    = {.state = {1.0, -1.0, 0.0, 0.0}};

  const sensor::Builder builder{.channels = {vertical, horizontal}, .antenna = sensor::PencilBeamAntenna{}};

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 1> los{{{150.0, 40.0}}};
  const auto [obsels, meta] = builder(pos, los, Vector2{Body::Earth::a, Body::Earth::b});

  const Stokvec v = obsels[0].weight_matrix()[0, 0];
  const Stokvec h = obsels[1].weight_matrix()[0, 0];
  for (Size i = 0; i < 4; i++) {
    assert_close(v[i], vertical.polarization.state[i], 0.0, "vertical channel weight");
    assert_close(h[i], horizontal.polarization.state[i], 0.0, "horizontal channel weight");
  }

  // The heterodyne mapping keeps the polarization of the channel it maps
  const sensor::HeterodyneFrequencyRange backend(50e9, Vector2{0.0, sensor::FrequencyRange::inf});
  const auto mapped = sensor::make_bandpass_channels(backend, std::span<const sensor::Channel>{&horizontal, 1});
  for (Size i = 0; i < 4; i++) {
    assert_close(mapped.front().polarization.state[i], horizontal.polarization.state[i], 0.0, "mapped polarization");
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_sensor_builder_polarization_per_channel failed:\n{}", e.what()));
}

//! A cross-track scan: the boresight is in the plane across the flight direction, and the polarization turns with
//! the signed scan angle; a bare line of sight cannot be built for such a channel because its side is unknown
void test_sensor_builder_cross_track() try {
  using Conversion::cosd, Conversion::sind;

  sensor::Channel channel = sensor::DiracChannel{100e9};
  channel.polarization    = {.state = {1.0, 1.0, 0.0, 0.0}, .scan_rotation = 1.0};

  const sensor::Builder builder{.channels = {channel}, .antenna = sensor::PencilBeamAntenna{}};

  const std::array<Vector3, 1> pos{{{800e3, 10.0, 20.0}}};
  const std::array<Numeric, 1> heading{{30.0}};
  const Vector                 scan{-20.0, 0.0, 20.0};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const auto [obsels, meta] = builder.cross_track(pos, heading, scan, ell);
  if (obsels.size() != 3 or meta.size() != 3) arts::user_error("Expected 3 obsels and meta blocks");

  const std::array<Vector2, 3> expected_los{{{160.0, -60.0}, {180.0, 120.0}, {160.0, 120.0}}};
  const std::array<Stokvec, 3> expected_weight{
      {{1.0, cosd(40.0), -sind(40.0), 0.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, cosd(40.0), sind(40.0), 0.0}}};
  for (Size i = 0; i < 3; i++) {
    const auto& los = obsels[i].poslos_grid()[0].los;
    assert_close(los[0], expected_los[i][0], 1e-12, "cross-track zenith");
    assert_close(los[1], expected_los[i][1], 1e-12, "cross-track azimuth");
    const Stokvec w = obsels[i].weight_matrix()[0, 0];
    for (Size k = 0; k < 4; k++) assert_close(w[k], expected_weight[i][k], 1e-14, "cross-track polarization");
  }

  bool threw = false;
  try {
    const std::array<Vector2, 1> los{{{160.0, 120.0}}};
    static_cast<void>(builder(pos, los, ell));
  } catch (const std::runtime_error&) { threw = true; }
  if (not threw) arts::user_error("A scan-dependent polarization needs a scan geometry");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_sensor_builder_cross_track failed:\n{}", e.what()));
}
}  // namespace

//! A camera gives one obsel per pixel and channel, by pixel and then channel, and an image of metadata per geometry
void test_sensor_builder_with_camera() try {
  const sensor::CameraAntenna camera{.grid = {.rows = AscendingGrid{-1.0, 1.0}, .cols = AscendingGrid{-1.0, 0.0, 1.0}}};
  const sensor::Builder       builder{.channels = {sensor::BoxChannel{AscendingGrid{100.0, 101.0}}, sensor::DiracChannel{200.0}},
                                      .antenna  = camera};

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{120.0, 30.0}, {140.0, 50.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const auto [obsels, meta] = builder(pos, los, ell);

  constexpr Size npix = 6, nchan = 2;
  if (obsels.size() != 2 * npix * nchan) arts::user_error("Expected {} obsels, got {}", 2 * npix * nchan, obsels.size());
  if (meta.size() != 2) arts::user_error("Expected one meta block per geometry, got {}", meta.size());

  for (Size g = 0; g < 2; g++) {
    const auto& image = std::get<sensor::CameraGriddedField>(meta[g].data);
    if (meta[g].count() != npix * nchan) arts::user_error("Meta block {} must count every pixel and channel", g);
    if (image.grid<0>().size() != 2 or image.grid<1>().size() != 3 or image.grid<2>().size() != nchan)
      arts::user_error("Meta block {} must be rows by columns by channels", g);
    for (Size i = 0; i < 2; i++) assert_close(image.grid<0>()[i], camera.grid.rows[i], 0.0, "meta row offset");
    for (Size j = 0; j < 3; j++) assert_close(image.grid<1>()[j], camera.grid.cols[j], 0.0, "meta column offset");
    for (Size c = 0; c < nchan; c++) assert_close(image.grid<2>()[c], static_cast<Numeric>(c), 0.0, "meta channel");

    const auto pixels = camera(sensor::DiracChannel{200.0}, Stokvec{1.0, 0.0, 0.0, 0.0}, pos[g], los[g], ell);
    for (Size p = 0; p < npix; p++) {
      const auto& box   = obsels[g * npix * nchan + p * nchan + 0];
      const auto& dirac = obsels[g * npix * nchan + p * nchan + 1];

      if (box.f_grid().size() != 2 or dirac.f_grid().size() != 1)
        arts::user_error("Obsels of geometry {} pixel {} must be the box and then the Dirac channel", g, p);
      if (box.poslos_grid_ptr() != dirac.poslos_grid_ptr())
        arts::user_error("The channels of one pixel must share its lines of sight");
      if (box.f_grid_ptr() != obsels[p * nchan].f_grid_ptr())
        arts::user_error("A channel must have the same frequency grid at every pixel and geometry");

      const auto& expected = pixels[p].poslos_grid()[0].los;
      assert_close(box.poslos_grid()[0].los[0], expected[0], 1e-12, "pixel zenith");
      assert_close(box.poslos_grid()[0].los[1], expected[1], 1e-12, "pixel azimuth");
    }
  }

  bool threw = false;
  try {
    const std::array<Numeric, 1> heading{{0.0}};
    static_cast<void>(builder.cross_track(std::span{pos}.first(1), heading, Vector{0.0}, ell));
  } catch (const std::runtime_error&) { threw = true; }
  if (not threw) arts::user_error("A cross-track scan of a camera must be refused");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_sensor_builder_with_camera failed:\n{}", e.what()));
}

int main() try {
  test_sensor_builder_returns_meta_per_geometry();
  test_sensor_builder_rejects_mismatched_geometry_counts();
  test_unflatten_updates_shared_poslos_grids();
  test_sensor_builder_polarization_per_channel();
  test_sensor_builder_cross_track();
  test_sensor_builder_with_camera();
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "test_sensor_builder failed:\n{}", e.what());
  return EXIT_FAILURE;
}
