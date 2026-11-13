#include <arts_constants.h>
#include <arts_conversions.h>
#include <frequency_bandpass_filters.h>
#include <frequency_range_selection.h>
#include <planet_data.h>
#include <sensor.h>
#include <sensor_builder.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string_view>

namespace {
void assert_close(Numeric actual, Numeric expected, Numeric tol, std::string_view name) {
  if (std::abs(actual - expected) > tol) arts::user_error("{} mismatch: actual {} expected {}", name, actual, expected);
}

//! The rows and columns of the weights of an element
std::pair<std::vector<Size>, std::vector<Size>> rows_and_cols(const sensor::SparseStokvecMatrix& w) {
  std::vector<Size> rows, cols;
  for (const auto& e : w) {
    rows.push_back(e.irow);
    cols.push_back(e.icol);
  }
  return {rows, cols};
}

//! One sensor for all placements: the frequencies of all channels, every placement's line of sight, and the elements
//! by placement and then channel, as the placement by channel meta
void test_sensor_builder_meta_has_placement_axis() try {
  sensor::Builder builder{.channels = {sensor::BoxChannel{AscendingGrid{100.0, 101.0}}, sensor::DiracChannel{200.0}},
                          .antenna  = sensor::PencilBeamAntenna{}};

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{20.0, 30.0}, {40.0, 50.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const Sensor sensor = builder(pos, los, ell);

  if (sensor.elements.size() != 4) arts::user_error("Expected 4 elements, got {}", sensor.elements.size());
  if (sensor.freq_grid.size() != 3) arts::user_error("Expected the 3 frequencies of both channels");
  if (sensor.geometry.size() != 2) arts::user_error("Expected one line of sight per placement");
  if (sensor::count(sensor.meta) != 4) arts::user_error("The meta must describe every element");

  const auto& gf = std::get<SortedGriddedField2>(sensor.meta);
  if (gf.gridname<0>() != "placement" or gf.gridname<1>() != "channel")
    arts::user_error("Builder meta must be placement by channel");
  for (Size i = 0; i < 2; i++) assert_close(gf.grid<0>()[i], static_cast<Numeric>(i), 0.0, "meta placement");
  for (Size i = 0; i < 2; i++) assert_close(gf.grid<1>()[i], static_cast<Numeric>(i), 0.0, "meta channel");

  if (sensor.geometry[0].pos != pos[0]) arts::user_error("First placement position mismatch");
  if (sensor.geometry[1].pos != pos[1]) arts::user_error("Second placement position mismatch");

  const std::array<std::vector<Size>, 4> expected_rows{{{0, 0}, {0}, {1, 1}, {1}}};
  const std::array<std::vector<Size>, 4> expected_cols{{{0, 1}, {2}, {0, 1}, {2}}};
  for (Size i = 0; i < 4; i++) {
    const auto [rows, cols] = rows_and_cols(sensor.elements[i]);
    if (rows != expected_rows[i] or cols != expected_cols[i])
      arts::user_error("Element {} must weight its placement's line of sight and its channel's frequencies", i);
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_sensor_builder_meta_has_placement_axis failed:\n{}", e.what()));
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

void test_unflatten_updates_the_sensor() try {
  sensor::SparseStokvecMatrix weights(1, 1);
  weights[0, 0] = {1.0, 0.0, 0.0, 0.0};

  SortedGriddedField1 meta;
  meta.data_name     = "test";
  meta.gridname<0>() = "channel";
  meta.grid<0>()     = AscendingGrid{0.0, 1.0};
  meta.data.resize(2);
  meta.data = 0.0;

  Sensor sensor{.freq_grid = AscendingGrid{100.0},
                .geometry  = SensorPosLosVector{{.pos = {600e3, 10.0, 20.0}, .los = {30.0, 40.0}}},
                .elements  = {weights, weights},
                .meta      = std::move(meta)};
  sensor.check();

  sensor.unflatten(Vector{15.0}, SensorKeyType::zen);
  assert_close(sensor.geometry[0].los[0], 15.0, 0.0, "updated zenith");
  assert_close(sensor.geometry[0].los[1], 40.0, 0.0, "azimuth must remain unchanged");
  assert_close(sensor.flat(SensorKeyType::zen)[0], 15.0, 0.0, "flat zenith");

  sensor.unflatten(Vector{101.0}, SensorKeyType::freq);
  assert_close(sensor.freq_grid[0], 101.0, 0.0, "updated frequency");
  assert_close(sensor.flat(SensorKeyType::freq)[0], 101.0, 0.0, "flat frequency");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_unflatten_updates_the_sensor failed:\n{}", e.what()));
}

//! The measurement elements of the sensors, in order, are the elements of the measurement vector
void test_measurement_elements() try {
  const sensor::Builder        builder{.channels = {sensor::DiracChannel{100.0}, sensor::DiracChannel{200.0}},
                                       .antenna  = sensor::PencilBeamAntenna{}};
  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 1> los{{{150.0, 40.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const ArrayOfSensor sensors{builder(pos, los, ell), builder(pos, los, ell), builder(pos, los, ell)};

  if (measurement_size(sensors) != 6) arts::user_error("Expected 6 measurement elements");
  for (Size i = 0; i <= 3; i++) {
    if (measurement_start(sensors, i) != 2 * i) arts::user_error("Sensor {} must start at {}", i, 2 * i);
  }
  for (Size iy = 0; iy < 6; iy++) {
    if (measurement_element(sensors, iy) != std::pair<Size, Size>{iy / 2, iy % 2})
      arts::user_error("Measurement element {} must be element {} of sensor {}", iy, iy % 2, iy / 2);
  }

  bool threw = false;
  try {
    static_cast<void>(measurement_element(sensors, 6));
  } catch (const std::runtime_error&) { threw = true; }
  if (not threw) arts::user_error("A measurement element beyond the sensors must be refused");
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_measurement_elements failed:\n{}", e.what()));
}

//! The index lists, for each geometry, the elements with weights there and where those weights are
void test_geometry_index() try {
  constexpr Size ngeom = 4, nfreq = 4;

  //! Two placements of two geometries each, and two channels of two frequencies each; then an element
  //! without weights, and one with weights at the first and last geometries only
  const std::array<std::vector<std::pair<Size, Size>>, 6> weights{{
      {{0, 0}, {0, 1}, {1, 0}, {1, 1}},
      {{0, 2}, {0, 3}, {1, 2}, {1, 3}},
      {{2, 0}, {2, 1}, {3, 0}, {3, 1}},
      {{2, 2}, {2, 3}, {3, 2}, {3, 3}},
      {},
      {{0, 1}, {3, 2}},
  }};

  Sensor sensor{.freq_grid = AscendingGrid{1.0, 2.0, 3.0, 4.0}, .geometry = SensorPosLosVector(ngeom)};
  for (const auto& w : weights) {
    sensor::SparseStokvecMatrix& elem = sensor.elements.emplace_back(ngeom, nfreq);
    for (const auto& [g, f] : w) elem[g, f] = Stokvec{1.0, 0.0, 0.0, 0.0};
  }

  struct Expected {
    Size element, begin, end;
  };
  const std::array<std::vector<Expected>, ngeom> expected{{
      {{0, 0, 2}, {1, 0, 2}, {5, 0, 1}},
      {{0, 2, 4}, {1, 2, 4}},
      {{2, 0, 2}, {3, 0, 2}},
      {{2, 2, 4}, {3, 2, 4}, {5, 1, 2}},
  }};
  const std::array<Size, 6>                      expected_counts{2, 2, 2, 2, 0, 2};

  const sensor::GeometryIndex index(sensor);

  for (Size ip = 0; ip < ngeom; ip++) {
    const auto got = index[ip];
    if (got.size() != expected[ip].size())
      arts::user_error("Geometry {} has {} elements, expected {}", ip, got.size(), expected[ip].size());

    for (Size k = 0; k < got.size(); k++) {
      const auto& [element, begin, end] = expected[ip][k];
      const auto& all                   = sensor.elements[element].vector();
      if (got[k].element != element or got[k].weights.data() != all.data() + begin or
          got[k].weights.size() != end - begin)
        arts::user_error(
            "Entry {} of geometry {} must be the weights {} to {} of element {}", k, ip, begin, end, element);
    }
  }

  for (Size ie = 0; ie < expected_counts.size(); ie++) {
    if (index.geometry_count(ie) != expected_counts[ie])
      arts::user_error(
          "Element {} has weights at {} geometries, expected {}", ie, index.geometry_count(ie), expected_counts[ie]);
  }

  //! A weight outside the geometries or the frequencies of the sensor is refused
  for (const auto& [g, f] : std::array<std::pair<Size, Size>, 2>{{{ngeom, 0}, {0, nfreq}}}) {
    Sensor bad            = sensor;
    bad.elements[4][g, f] = Stokvec{1.0, 0.0, 0.0, 0.0};

    bool threw = false;
    try {
      static_cast<void>(sensor::GeometryIndex(bad));
    } catch (const std::runtime_error&) { threw = true; }
    if (not threw) arts::user_error("A weight at geometry {} and frequency {} must be refused", g, f);
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("test_geometry_index failed:\n{}", e.what()));
}

void test_sensor_builder_polarization_per_channel() try {
  sensor::Channel vertical   = sensor::DiracChannel{100.0};
  sensor::Channel horizontal = sensor::DiracChannel{200.0};
  vertical.polarization      = {.state = {1.0, 1.0, 0.0, 0.0}};
  horizontal.polarization    = {.state = {1.0, -1.0, 0.0, 0.0}};

  const sensor::Builder builder{.channels = {vertical, horizontal}, .antenna = sensor::PencilBeamAntenna{}};

  const std::array<Vector3, 1> pos{{{600e3, 10.0, 20.0}}};
  const std::array<Vector2, 1> los{{{150.0, 40.0}}};
  const Sensor                 sensor = builder(pos, los, Vector2{Body::Earth::a, Body::Earth::b});

  const Stokvec v = sensor.elements[0][0, 0];
  const Stokvec h = sensor.elements[1][0, 1];
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

  const Sensor sensor = builder.cross_track(pos, heading, scan, ell);
  if (sensor.elements.size() != 3 or sensor.geometry.size() != 3)
    arts::user_error("Expected 3 elements and 3 lines of sight");

  const std::array<Vector2, 3> expected_los{{{160.0, -60.0}, {180.0, 120.0}, {160.0, 120.0}}};
  const std::array<Stokvec, 3> expected_weight{
      {{1.0, cosd(40.0), -sind(40.0), 0.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, cosd(40.0), sind(40.0), 0.0}}};
  for (Size i = 0; i < 3; i++) {
    const auto& los = sensor.geometry[i].los;
    assert_close(los[0], expected_los[i][0], 1e-12, "cross-track zenith");
    assert_close(los[1], expected_los[i][1], 1e-12, "cross-track azimuth");
    const Stokvec w = sensor.elements[i][i, 0];
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

//! A camera gives one element per placement, pixel and channel, in that order, as its placement by rows by columns
//! by channels meta
void test_sensor_builder_with_camera() try {
  const sensor::CameraAntenna camera{.grid = {.rows = AscendingGrid{-1.0, 1.0}, .cols = AscendingGrid{-1.0, 0.0, 1.0}}};
  const sensor::Builder       builder{
      .channels = {sensor::BoxChannel{AscendingGrid{100.0, 101.0}}, sensor::DiracChannel{200.0}}, .antenna = camera};

  const std::array<Vector3, 2> pos{{{600e3, 10.0, 20.0}, {601e3, 11.0, 21.0}}};
  const std::array<Vector2, 2> los{{{120.0, 30.0}, {140.0, 50.0}}};
  const Vector2                ell{Body::Earth::a, Body::Earth::b};

  const Sensor sensor = builder(pos, los, ell);

  constexpr Size npix = 6, nchan = 2;
  if (sensor.elements.size() != 2 * npix * nchan)
    arts::user_error("Expected {} elements, got {}", 2 * npix * nchan, sensor.elements.size());
  if (sensor.geometry.size() != 2 * npix) arts::user_error("Expected one line of sight per placement and pixel");

  const auto& image = std::get<SortedGriddedField4>(sensor.meta);
  if (image.grid<0>().size() != 2 or image.grid<1>().size() != 2 or image.grid<2>().size() != 3 or
      image.grid<3>().size() != nchan)
    arts::user_error("The meta must be placements by rows by columns by channels");
  for (Size i = 0; i < 2; i++) assert_close(image.grid<1>()[i], camera.grid.rows[i], 0.0, "meta row offset");
  for (Size j = 0; j < 3; j++) assert_close(image.grid<2>()[j], camera.grid.cols[j], 0.0, "meta column offset");
  for (Size c = 0; c < nchan; c++) assert_close(image.grid<3>()[c], static_cast<Numeric>(c), 0.0, "meta channel");

  for (Size g = 0; g < 2; g++) {
    const Sensor pixels =
        sensor::antenna_sensor(camera, sensor::DiracChannel{200.0}, Stokvec{1.0, 0.0, 0.0, 0.0}, pos[g], los[g], ell);
    for (Size p = 0; p < npix; p++) {
      const auto [box_rows, box_cols]     = rows_and_cols(sensor.elements[g * npix * nchan + p * nchan + 0]);
      const auto [dirac_rows, dirac_cols] = rows_and_cols(sensor.elements[g * npix * nchan + p * nchan + 1]);

      const Size row = g * npix + p;
      if (box_rows != std::vector<Size>{row, row} or dirac_rows != std::vector<Size>{row})
        arts::user_error("Placement {} pixel {} must weight its own line of sight", g, p);
      if (box_cols != std::vector<Size>{0, 1} or dirac_cols != std::vector<Size>{2})
        arts::user_error("Placement {} pixel {} must be the box and then the Dirac channel", g, p);

      const auto& expected = pixels.geometry[p].los;
      assert_close(sensor.geometry[row].los[0], expected[0], 1e-12, "pixel zenith");
      assert_close(sensor.geometry[row].los[1], expected[1], 1e-12, "pixel azimuth");
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
  test_sensor_builder_meta_has_placement_axis();
  test_sensor_builder_rejects_mismatched_geometry_counts();
  test_unflatten_updates_the_sensor();
  test_measurement_elements();
  test_geometry_index();
  test_sensor_builder_polarization_per_channel();
  test_sensor_builder_cross_track();
  test_sensor_builder_with_camera();
  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "test_sensor_builder failed:\n{}", e.what());
  return EXIT_FAILURE;
}
