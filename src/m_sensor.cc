#include <arts_conversions.h>
#include <arts_omp.h>
#include <debug.h>
#include <rtepack.h>
#include <sensor.h>
#include <workspace.h>

#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <format>
#include <stdexcept>

namespace {
/** Appends the sensor that make makes for its name
 *
 *  The name is the given name, or "sensor-N" with N the position of the sensor
 *  in the list if the given name is empty.  Errors while making the sensor
 *  name it.
 */
template <typename Make> void append_sensor(ArrayOfSensor& measurement_sensors, const String& name, Make&& make) {
  const String sensor_name = name.empty() ? std::format("sensor-{}", measurement_sensors.size() + 1) : name;

  Sensor sensor;
  try {
    sensor = make(sensor_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::format("Sensor \"{}\":\n{}", sensor_name, e.what()));
  }

  sensor.check();
  measurement_sensors.push_back(std::move(sensor));
}

//! The structure of the elements of a spectrometer, one per frequency, named as the sensor
SortedGriddedField1 frequency_meta(const String& name, const AscendingGrid& freq_grid) {
  SortedGriddedField1 gf;
  gf.data_name     = name;
  gf.gridname<0>() = "frequency";
  gf.grid<0>()     = freq_grid;
  gf.data.resize(freq_grid.size());
  gf.data = 0.0;
  return gf;
}
}  // namespace

void measurement_sensorsInit(ArrayOfSensor& measurement_sensors) {
  ARTS_TIME_REPORT

  measurement_sensors = {};
}

void measurement_sensorsAddSimple(ArrayOfSensor&       measurement_sensors,
                                  const AscendingGrid& freq_grid,
                                  const Vector3&       pos,
                                  const Vector2&       los,
                                  const Stokvec&       pol,
                                  const String&        name) try {
  ARTS_TIME_REPORT

  append_sensor(measurement_sensors, name, [&](const String& sensor_name) {
    const Size n = freq_grid.size();

    Sensor sensor{.freq_grid = freq_grid,
                  .geometry  = SensorPosLosVector{{.pos = pos, .los = los}},
                  .elements  = std::vector<sensor::SparseStokvecMatrix>(n, sensor::SparseStokvecMatrix(1, n)),
                  .meta      = frequency_meta(sensor_name, freq_grid)};

    for (Size i = 0; i < n; i++) sensor.elements[i][0, i] = pol;

    return sensor;
  });
}
ARTS_METHOD_ERROR_CATCH

namespace {
Sensor gaussian_sensor(const String&        name,
                       const AscendingGrid& freq_grid,
                       const Vector&        stds,
                       const Vector3&       pos,
                       const Vector2&       los,
                       const Stokvec&       pol) {
  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  const Size n       = freq_grid.size();
  const Size nonzero = pol.nonzero_components();

  if (n < 2) arts::user_error("Must have a frequency grid");
  if (n != stds.size()) arts::user_error("Must have a standard deviation for each frequency point");
  if (stdr::any_of(stds, Cmp::le<0.0>())) arts::user_error("Standard deviation must be positive.\nstds := {:B,}", stds);
  if (nonzero == 0) arts::user_error("The polarization pol must have at least one non-zero Stokes component");

  Sensor sensor{.freq_grid = freq_grid,
                .geometry  = SensorPosLosVector{{.pos = pos, .los = los}},
                .elements  = std::vector<sensor::SparseStokvecMatrix>(n),
                .meta      = frequency_meta(name, freq_grid)};

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n; i++) {
    try {
      StokvecMatrix w(1, n, pol);

      const gauss dist(freq_grid[i], stds[i]);
      for (Size j = 0; j < n; j++) { w[0, j] *= pdf(dist, freq_grid[j]); }

      sensor.elements[i] = sensor::SparseStokvecMatrix(w);
      sensor::normalize(sensor.elements[i], pol);
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (error.size()) arts::user_error(error);

  return sensor;
}
}  // namespace

void measurement_sensorsAddVectorGaussian(ArrayOfSensor&       measurement_sensors,
                                          const AscendingGrid& freq_grid,
                                          const Vector&        stds,
                                          const Vector3&       pos,
                                          const Vector2&       los,
                                          const Stokvec&       pol,
                                          const String&        name) try {
  ARTS_TIME_REPORT

  append_sensor(measurement_sensors, name, [&](const String& sensor_name) {
    return gaussian_sensor(sensor_name, freq_grid, stds, pos, los, pol);
  });
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorsAddSimpleGaussian(ArrayOfSensor&       measurement_sensors,
                                          const AscendingGrid& freq_grid,
                                          const Numeric&       std,
                                          const Vector3&       pos,
                                          const Vector2&       los,
                                          const Stokvec&       pol,
                                          const String&        name) {
  ARTS_TIME_REPORT

  const Vector stds(freq_grid.size(), std);
  measurement_sensorsAddVectorGaussian(measurement_sensors, freq_grid, stds, pos, los, pol, name);
}

namespace {
Sensor zenith_response_sensor(const String&        name,
                              const AscendingGrid& freq_grid,
                              const Vector3&       pos,
                              const Vector2&       los,
                              const Stokvec&       pol,
                              const AscendingGrid& dzen_grid,
                              const Vector&        zen_weights) {
  const Size n_freq = freq_grid.size();
  const Size n_za   = dzen_grid.size();

  if (n_za != zen_weights.size()) arts::user_error("dzen_grid and zen_weights must have same size");
  const Size nonzero = pol.nonzero_components();
  if (nonzero == 0) arts::user_error("The polarization pol must have at least one non-zero Stokes component");

  Sensor sensor{.freq_grid = freq_grid,
                .geometry  = SensorPosLosVector(n_za),
                .elements  = std::vector<sensor::SparseStokvecMatrix>(n_freq),
                .meta      = frequency_meta(name, freq_grid)};

  for (Size i = 0; i < n_za; ++i) { sensor.geometry[i] = SensorPosLos(pos, los + Vector2{dzen_grid[i], 0.0}); }

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n_freq; i++) {
    try {
      sensor::SparseStokvecMatrix w(n_za, n_freq);
      for (Size j = 0; j < n_za; j++) w[j, i] = zen_weights[j] * pol;
      sensor::normalize(w, pol);
      sensor.elements[i] = std::move(w);
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }
  if (error.size()) arts::user_error(error);

  return sensor;
}
}  // namespace

void measurement_sensorsAddGaussianZenith(ArrayOfSensor&       measurement_sensors,
                                          const AscendingGrid& freq_grid,
                                          const Vector3&       pos,
                                          const Vector2&       los,
                                          const Stokvec&       pol,
                                          const AscendingGrid& dzen_grid,
                                          const Numeric&       std_zen,
                                          const String&        name) try {
  ARTS_TIME_REPORT

  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  append_sensor(measurement_sensors, name, [&](const String& sensor_name) {
    if (std_zen <= 0) arts::user_error("Standard deviation must be positive.");
    if (dzen_grid.size() == 0) arts::user_error("dzen_grid cannot be empty.");

    Vector      zen_weights(dzen_grid.size());
    const gauss dist(0.0, std_zen);
    for (Size i = 0; i < dzen_grid.size(); ++i) { zen_weights[i] = pdf(dist, dzen_grid[i]); }

    return zenith_response_sensor(sensor_name, freq_grid, pos, los, pol, dzen_grid, zen_weights);
  });
}
ARTS_METHOD_ERROR_CATCH

namespace {
template <typename T, typename... Grids>
Sensor raw_sensor_sensor(const String&                               name,
                         const AscendingGrid&                        freq_grid,
                         const Vector3&                              pos,
                         const Vector2&                              los,
                         const matpack::gridded_data_t<T, Grids...>& raw_sensor,
                         const bool                                  normalize)
    requires(sizeof...(Grids) <= 6 and (std::same_as<Grids, AscendingGrid> and ...) and
             (std::same_as<Numeric, T> or std::same_as<Stokvec, T>)) {
  constexpr Size idf   = 0;
  constexpr Size idza  = 1;
  constexpr Size idaa  = 2;
  constexpr Size idalt = 3;
  constexpr Size idlat = 4;
  constexpr Size idlon = 5;

  struct InternalGrid {
    const AscendingGrid* grid{nullptr};
    Size                 i{};

    [[nodiscard]] bool    exists() const { return grid != nullptr; }
    [[nodiscard]] Size    size() const { return exists() ? grid->size() : 1; }
    [[nodiscard]] Numeric get(Size i) const { return exists() ? (*grid)[i] : 0.0; }

    InternalGrid() = default;
    InternalGrid(const char* name, const matpack::gridded_data_t<T, Grids...>& raw_sensor) {
      auto ptr = stdr::find(raw_sensor.grid_names, name);
      i        = stdr::distance(stdr::begin(raw_sensor.grid_names), ptr);
      if constexpr (sizeof...(Grids) >= 1)
        if (i == 0) grid = &raw_sensor.template grid<0>();
      if constexpr (sizeof...(Grids) >= 2)
        if (i == 1) grid = &raw_sensor.template grid<1>();
      if constexpr (sizeof...(Grids) >= 3)
        if (i == 2) grid = &raw_sensor.template grid<2>();
      if constexpr (sizeof...(Grids) >= 4)
        if (i == 3) grid = &raw_sensor.template grid<3>();
      if constexpr (sizeof...(Grids) >= 5)
        if (i == 4) grid = &raw_sensor.template grid<4>();
      if constexpr (sizeof...(Grids) >= 6)
        if (i == 5) grid = &raw_sensor.template grid<5>();
    }
  };

  std::array<InternalGrid, 6> grids;
  grids[idf]   = InternalGrid("dfreq", raw_sensor);
  grids[idza]  = InternalGrid("dzen", raw_sensor);
  grids[idaa]  = InternalGrid("dazi", raw_sensor);
  grids[idalt] = InternalGrid("dalt", raw_sensor);
  grids[idlat] = InternalGrid("dlat", raw_sensor);
  grids[idlon] = InternalGrid("dlon", raw_sensor);

  const Size NF   = grids[idf].size();
  const Size NZA  = grids[idza].size();
  const Size NAA  = grids[idaa].size();
  const Size NALT = grids[idalt].size();
  const Size NLAT = grids[idlat].size();
  const Size NLON = grids[idlon].size();

  if (not raw_sensor.ok() or raw_sensor.data.size() != NF * NZA * NAA * NALT * NLAT * NLON)
    arts::user_error(R"(The raw sensor weights data have an unexpected size or are not OK.

The internal method wants to view the data as a 6-dimensional
array of shape:  [{0:B}, {1:}, {2:}, {3:}, {4:}, {5:}]
where the dimensions are:
  Frequency perturbance     ("dfreq"): {0:}
  Zenith angle perturbance  ("dzen"):  {1:}
  Azimuth angle perturbance ("dazi"):  {2:}
  Altitude perturbance      ("dalt"):  {3:}
  Latitude perturbance      ("dlat"):  {4:}
  Longitude perturbance     ("dlon"):  {5:}

But the array has the following shape: {6:B,}

Note that the grid names must match the quoted grid names above, and that
the order of the grids is enforced.  The raw sensor data sent in has the
following grid names: {7:B,}
)",
                     NF,
                     NZA,
                     NAA,
                     NALT,
                     NLAT,
                     NLON,
                     raw_sensor.data.shape(),
                     raw_sensor.grid_names);

  Index lowest = -1;
  for (Size i = 0; i < grids.size(); i++) {
    if (not grids[i].exists()) continue;

    if (lowest != -1 and static_cast<Index>(grids[i].i) <= lowest)
      arts::user_error(R"(The raw sensor input grids must be sorted as follows:

Frequency grid perturbance (named "dfreq"),
Zenith angle perturbance   (named "dzen"),
Azimuth angle perturbance  (named "dazi"),
Altitude perturbance       (named "dalt"),
Latitude perturbance       (named "dlat"),
Longitude perturbance      (named "dlon")

Note that only the grid names above are allowed, and the order is
enforced.

Your sorting is not correct, and instead reads:
{:B,}
)",
                       raw_sensor.grid_names);
    lowest = grids[i].i;
  }

  auto view = raw_sensor.data.view_as(NF, NZA, NAA, NALT, NLAT, NLON);

  // The perturbed frequencies of every frequency of the grid
  std::vector<Numeric> all_freqs;
  all_freqs.reserve(freq_grid.size() * NF);
  for (auto f : freq_grid) {
    for (Size ifreq = 0; ifreq < NF; ifreq++) all_freqs.push_back(f + grids[idf].get(ifreq));
  }
  stdr::sort(all_freqs);
  all_freqs.erase(std::unique(all_freqs.begin(), all_freqs.end()), all_freqs.end());

  Sensor sensor{.freq_grid = AscendingGrid{std::move(all_freqs)},
                .geometry  = SensorPosLosVector(NZA * NAA * NALT * NLAT * NLON),
                .elements  = {},
                .meta      = frequency_meta(name, freq_grid)};

  auto wposlos = sensor.geometry.view_as(NZA, NAA, NALT, NLAT, NLON);
  for (Size iza = 0; iza < NZA; iza++) {
    for (Size iaa = 0; iaa < NAA; iaa++) {
      for (Size ialt = 0; ialt < NALT; ialt++) {
        for (Size ilat = 0; ilat < NLAT; ilat++) {
          for (Size ilon = 0; ilon < NLON; ilon++) {
            const Vector3 dpos{grids[idalt].get(ialt), grids[idlat].get(ilat), grids[idlon].get(ilon)};
            const Vector2 dlos{grids[idza].get(iza), grids[idaa].get(iaa)};

            wposlos[iza, iaa, ialt, ilat, ilon] = SensorPosLos(pos + dpos, los + dlos);
          }
        }
      }
    }
  }

  StokvecMatrix w(sensor.geometry.size(), sensor.freq_grid.size(), Stokvec{0, 0, 0, 0});
  auto          wview = w.view_as(NZA, NAA, NALT, NLAT, NLON, sensor.freq_grid.size());

  sensor.elements.reserve(freq_grid.size());
  for (auto f : freq_grid) {
    w = Stokvec{0, 0, 0, 0};

    for (Size ifreq = 0; ifreq < NF; ifreq++) {
      const auto iv = static_cast<Size>(
          stdr::distance(sensor.freq_grid.begin(), stdr::lower_bound(sensor.freq_grid, f + grids[idf].get(ifreq))));

      for (Size iza = 0; iza < NZA; iza++) {
        for (Size iaa = 0; iaa < NAA; iaa++) {
          for (Size ialt = 0; ialt < NALT; ialt++) {
            for (Size ilat = 0; ilat < NLAT; ilat++) {
              for (Size ilon = 0; ilon < NLON; ilon++) {
                wview[iza, iaa, ialt, ilat, ilon, iv] = Stokvec{view[ifreq, iza, iaa, ialt, ilat, ilon]};
              }
            }
          }
        }
      }
    }

    sensor.elements.emplace_back(w);
    if (normalize) sensor::normalize(sensor.elements.back());
  }

  return sensor;
}
}  // namespace

void measurement_sensorsAddRawSensor(ArrayOfSensor&                                  measurement_sensors,
                                     const AscendingGrid&                            freq_grid,
                                     const Vector3&                                  pos,
                                     const Vector2&                                  los,
                                     const Generic<const SortedGriddedField1,
                                                   const SortedGriddedField2,
                                                   const SortedGriddedField3,
                                                   const SortedGriddedField4,
                                                   const SortedGriddedField5,
                                                   const SortedGriddedField6,
                                                   const StokvecSortedGriddedField1,
                                                   const StokvecSortedGriddedField2,
                                                   const StokvecSortedGriddedField3,
                                                   const StokvecSortedGriddedField4,
                                                   const StokvecSortedGriddedField5,
                                                   const StokvecSortedGriddedField6> raw_sensor,
                                     const Index&                                    normalize,
                                     const String&                                   name) try {
  ARTS_TIME_REPORT

  append_sensor(measurement_sensors, name, [&](const String& sensor_name) {
    return std::visit(
        [&](const auto& selected) {
          return raw_sensor_sensor(sensor_name, freq_grid, pos, los, *selected, static_cast<bool>(normalize));
        },
        raw_sensor);
  });
}
ARTS_METHOD_ERROR_CATCH

namespace {
Sensor camera_sensor(const String&        name,
                     const AscendingGrid& freq_grid,
                     const Vector3&       pos,
                     const Vector2&       los,
                     const Stokvec&       pol,
                     const Index&         n_h,
                     const Index&         n_w,
                     const Numeric&       ccd_h,
                     const Numeric&       ccd_w,
                     const Numeric&       focal_length,
                     const Numeric&       aperture_diameter,
                     const Numeric&       focus_distance) {
  const Size n_freq   = freq_grid.size();
  const Size n_pixels = static_cast<Size>(n_h) * static_cast<Size>(n_w);

  if (n_freq == 0) arts::user_error("freq_grid must not be empty");
  if (n_h < 1) arts::user_error("n_h must be at least 1, got {}", n_h);
  if (n_w < 1) arts::user_error("n_w must be at least 1, got {}", n_w);
  if (ccd_h <= 0) arts::user_error("ccd_h must be positive, got {}", ccd_h);
  if (ccd_w <= 0) arts::user_error("ccd_w must be positive, got {}", ccd_w);
  if (focal_length <= 0) arts::user_error("focal_length must be positive, got {}", focal_length);
  if (aperture_diameter <= 0) arts::user_error("aperture_diameter must be positive, got {}", aperture_diameter);
  if (focus_distance <= 0) arts::user_error("focus_distance must be positive, got {}", focus_distance);
  if (pol.nonzero_components() == 0)
    arts::user_error("The polarization pol must have at least one non-zero Stokes component");

  // Pixel pitch (distance between pixel centers)
  const Numeric pitch_h = ccd_h / static_cast<Numeric>(n_h);
  const Numeric pitch_w = ccd_w / static_cast<Numeric>(n_w);

  // Thin-lens adjusted image distance for the given focus distance.
  // Thin lens equation: 1/f = 1/d_o + 1/d_i
  //   => d_i = f * d_o / (d_o - f)
  if (focus_distance <= focal_length)
    arts::user_error("focus_distance ({}) must be greater than focal_length ({})", focus_distance, focal_length);

  const Numeric image_distance = focal_length * focus_distance / (focus_distance - focal_length);

  // Camera coordinate system via 3D rotation.
  //
  // We build an orthonormal basis (e_up, e_right, g_hat) where g_hat is the
  // optical axis.  The CCD plane is perpendicular to g_hat.  The lens
  // inverts the image, so a CCD pixel at (dx, dy) maps to scene direction:
  //
  //   v = image_distance * g_hat + dy * e_zen - dx * e_azi
  //
  // where e_zen points toward increasing zenith and e_azi toward increasing
  // azimuth.  The signs account for the lens inversion: CCD "up" (dy>0)
  // captures the scene below the optical axis (toward increasing zenith
  // for a horizontal camera).
  //
  // This approach replaces the additive dzen/dazi flat-plane approximation
  // and works correctly at all zenith angles including 0° and 180°.

  const Numeric zen0 = Conversion::deg2rad(los[0]);
  const Numeric azi0 = Conversion::deg2rad(los[1]);

  const Numeric szen = std::sin(zen0);
  const Numeric czen = std::cos(zen0);
  const Numeric sazi = std::sin(azi0);
  const Numeric cazi = std::cos(azi0);

  // Optical axis unit vector
  const Numeric gx = szen * cazi;
  const Numeric gy = szen * sazi;
  const Numeric gz = czen;

  // e_zen: unit vector toward increasing zenith (on the sphere)
  const Numeric ezx = czen * cazi;
  const Numeric ezy = czen * sazi;
  const Numeric ezz = -szen;

  // e_azi: unit vector toward increasing azimuth (on the sphere)
  const Numeric eax = -sazi;
  const Numeric eay = cazi;
  // eaz = 0

  SensorPosLosVector all_poslos(n_pixels);
  for (Index ih = 0; ih < n_h; ++ih) {
    for (Index iw = 0; iw < n_w; ++iw) {
      const Numeric dy = (static_cast<Numeric>(ih) - 0.5 * static_cast<Numeric>(n_h - 1)) * pitch_h;
      const Numeric dx = (static_cast<Numeric>(iw) - 0.5 * static_cast<Numeric>(n_w - 1)) * pitch_w;

      const Numeric vx = image_distance * gx + dy * ezx - dx * eax;
      const Numeric vy = image_distance * gy + dy * ezy - dx * eay;
      const Numeric vz = image_distance * gz + dy * ezz;

      const Numeric vnorm = std::sqrt(vx * vx + vy * vy + vz * vz);
      const Numeric nz    = vz / vnorm;

      const Numeric zen = Conversion::rad2deg(std::acos(std::clamp(nz, -1.0, 1.0)));
      const Numeric azi = Conversion::rad2deg(std::atan2(vy, vx));

      const Size pix_idx  = static_cast<Size>(ih) * static_cast<Size>(n_w) + static_cast<Size>(iw);
      all_poslos[pix_idx] = {.pos = pos, .los = {zen, azi}};
    }
  }

  // Each pixel (ih,iw) at frequency ifreq is an element with a sparse
  // weight matrix of shape (n_pixels, n_freq) with a single nonzero
  // entry at (pix_idx, ifreq).
  std::vector<sensor::SparseStokvecMatrix> elements(n_pixels * n_freq);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Index ih = 0; ih < n_h; ++ih) {
    for (Index iw = 0; iw < n_w; ++iw) {
      try {
        const Size pix_idx = static_cast<Size>(ih) * static_cast<Size>(n_w) + static_cast<Size>(iw);

        for (Size ifreq = 0; ifreq < n_freq; ++ifreq) {
          sensor::SparseStokvecMatrix w(n_pixels, n_freq);
          w[pix_idx, ifreq]                  = pol;
          elements[pix_idx * n_freq + ifreq] = std::move(w);
        }
      } catch (std::runtime_error& e) {
#pragma omp critical
        if (error.empty()) error = e.what();
      }
    }
  }

  if (error.size()) arts::user_error("{}", error);

  // Build camera metadata: gridded_data_t<Numeric, AscendingGrid, AscendingGrid, Vector>
  // grid<0> = row angular offsets from boresight (degrees)
  // grid<1> = col angular offsets from boresight (degrees)
  // grid<2> = frequency / color
  sensor::CameraGriddedField gf;
  gf.data_name     = name;
  gf.gridname<0>() = "row_offset_deg";
  gf.gridname<1>() = "col_offset_deg";
  gf.gridname<2>() = "frequency";

  // Row angular offsets
  Vector row_vec(n_h);
  for (Index ih = 0; ih < n_h; ++ih) {
    const Numeric dy = (static_cast<Numeric>(ih) - 0.5 * static_cast<Numeric>(n_h - 1)) * pitch_h;
    row_vec[ih]      = Conversion::rad2deg(std::atan2(dy, image_distance));
  }
  gf.grid<0>() = AscendingGrid{std::move(row_vec)};

  // Col angular offsets
  Vector col_vec(n_w);
  for (Index iw = 0; iw < n_w; ++iw) {
    const Numeric dx = (static_cast<Numeric>(iw) - 0.5 * static_cast<Numeric>(n_w - 1)) * pitch_w;
    col_vec[iw]      = Conversion::rad2deg(std::atan2(dx, image_distance));
  }
  gf.grid<1>() = AscendingGrid{std::move(col_vec)};

  // Frequency / color grid
  gf.grid<2>() = Vector(freq_grid);

  gf.data.resize(n_h, n_w, n_freq);
  gf.data = 0.0;

  return Sensor{.freq_grid = freq_grid,
                .geometry  = std::move(all_poslos),
                .elements  = std::move(elements),
                .meta      = std::move(gf)};
}
}  // namespace

void measurement_sensorsAddCamera(ArrayOfSensor&       measurement_sensors,
                                  const AscendingGrid& freq_grid,
                                  const Vector3&       pos,
                                  const Vector2&       los,
                                  const Stokvec&       pol,
                                  const Index&         n_h,
                                  const Index&         n_w,
                                  const Numeric&       ccd_h,
                                  const Numeric&       ccd_w,
                                  const Numeric&       focal_length,
                                  const Numeric&       aperture_diameter,
                                  const Numeric&       focus_distance,
                                  const String&        name) try {
  ARTS_TIME_REPORT

  append_sensor(measurement_sensors, name, [&](const String& sensor_name) {
    return camera_sensor(
        sensor_name, freq_grid, pos, los, pol, n_h, n_w, ccd_h, ccd_w, focal_length, aperture_diameter, focus_distance);
  });
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorsAddSensorBuilder(ArrayOfSensor&       measurement_sensors,
                                         const SurfaceField&  surf_field,
                                         const Vector3&       pos,
                                         const Vector2&       los,
                                         const SensorBuilder& sensor_builder,
                                         const String&        name) try {
  ARTS_TIME_REPORT

  append_sensor(measurement_sensors, name, [&](const String& sensor_name) {
    return sensor_builder(std::span{&pos, 1}, std::span{&los, 1}, surf_field.ellipsoid, sensor_name);
  });
}
ARTS_METHOD_ERROR_CATCH
