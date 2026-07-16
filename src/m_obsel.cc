#include <arts_conversions.h>
#include <arts_omp.h>
#include <debug.h>
#include <obsel.h>
#include <rtepack.h>
#include <sensor_meta_info.h>
#include <workspace.h>

#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include <memory>
#include <stdexcept>

void measurement_sensorInit(ArrayOfSensorObsel& measurement_sensor, ArrayOfSensorMetaInfo& measurement_sensor_meta) {
  ARTS_TIME_REPORT

  measurement_sensor      = {};
  measurement_sensor_meta = {};
}

void measurement_sensorAddSimple(ArrayOfSensorObsel&    measurement_sensor,
                                 ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                 const AscendingGrid&   freq_grid,
                                 const Vector3&         pos,
                                 const Vector2&         los,
                                 const Stokvec&         pol) try {
  ARTS_TIME_REPORT

  const Index n  = freq_grid.size();
  const Size  sz = measurement_sensor.size();

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(freq_grid);
  auto p = std::make_shared<const SensorPosLosVector>(SensorPosLosVector{{.pos = pos, .los = los}});

  StokvecMatrix w(1, n, 0);
  for (Index i = 0; i < n; i++) {
    w[0, i]                    = pol;
    measurement_sensor[i + sz] = {f, p, w};
    w[0, i]                    = 0.0;
  }

  // Metadata: 1D sorted gridded field with frequency grid
  SortedGriddedField1 gf;
  gf.data_name     = "ffts";
  gf.gridname<0>() = "frequency";
  gf.grid<0>()     = freq_grid;
  gf.data.resize(n);
  gf.data = 0.0;
  measurement_sensor_meta.push_back(SensorMetaInfo{std::move(gf)});
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorAddVectorGaussian(ArrayOfSensorObsel&    measurement_sensor,
                                         ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                         const AscendingGrid&   freq_grid,
                                         const Vector&          stds,
                                         const Vector3&         pos,
                                         const Vector2&         los,
                                         const Stokvec&         pol) try {
  ARTS_TIME_REPORT

  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  const Size n       = freq_grid.size();
  const Size sz      = measurement_sensor.size();
  const Size nonzero = pol.nonzero_components();

  ARTS_USER_ERROR_IF(n < 2, "Must have a frequency grid")
  ARTS_USER_ERROR_IF(n != stds.size(), "Must have a standard deviation for each frequency point")
  ARTS_USER_ERROR_IF(stdr::any_of(stds, Cmp::le<0.0>()), "Standard deviation must be positive.\nstds := {:B,}", stds)
  ARTS_USER_ERROR_IF(nonzero == 0, "pol is 0")

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(freq_grid);
  auto p = std::make_shared<const SensorPosLosVector>(SensorPosLosVector{{.pos = pos, .los = los}});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n; i++) {
    try {
      StokvecMatrix w(1, n, pol);

      const gauss dist(freq_grid[i], stds[i]);
      for (Size j = 0; j < n; j++) { w[0, j] *= pdf(dist, freq_grid[j]); }

      measurement_sensor[i + sz] = {f, p, w};
      measurement_sensor[i + sz].normalize(pol);
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "{}", error)

  // Metadata: 1D sorted gridded field with frequency grid
  SortedGriddedField1 gf;
  gf.data_name     = "ffts";
  gf.gridname<0>() = "frequency";
  gf.grid<0>()     = freq_grid;
  gf.data.resize(n);
  gf.data = 0.0;
  measurement_sensor_meta.push_back(SensorMetaInfo{std::move(gf)});
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorAddSimpleGaussian(ArrayOfSensorObsel&    measurement_sensor,
                                         ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                         const AscendingGrid&   freq_grid,
                                         const Numeric&         std,
                                         const Vector3&         pos,
                                         const Vector2&         los,
                                         const Stokvec&         pol) {
  ARTS_TIME_REPORT

  const Vector stds(freq_grid.size(), std);
  measurement_sensorAddVectorGaussian(measurement_sensor, measurement_sensor_meta, freq_grid, stds, pos, los, pol);
}

namespace {
void measurement_sensorAddZenithResponse(ArrayOfSensorObsel&                         measurement_sensor,
                                         const std::shared_ptr<const AscendingGrid>& f,
                                         const Vector3&                              pos,
                                         const Vector2&                              los,
                                         const Stokvec&                              pol,
                                         const AscendingGrid&                        dzen_grid,
                                         const Vector&                               zen_weights) {
  const Size n_freq = f->size();
  const Size n_za   = dzen_grid.size();
  const Size sz     = measurement_sensor.size();

  ARTS_USER_ERROR_IF(n_za != zen_weights.size(), "dzen_grid and zen_weights must have same size");
  const Size nonzero = pol.nonzero_components();
  ARTS_USER_ERROR_IF(nonzero == 0, "pol is 0")

  measurement_sensor.resize(sz + n_freq);

  auto p_vec = SensorPosLosVector(n_za);
  for (Size i = 0; i < n_za; ++i) { p_vec[i] = SensorPosLos(pos, los + Vector2{dzen_grid[i], 0.0}); }
  auto p = std::make_shared<const SensorPosLosVector>(p_vec);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n_freq; i++) {
    try {
      sensor::SparseStokvecMatrix w(n_za, n_freq);
      for (Size j = 0; j < n_za; j++) w[j, i] = zen_weights[j] * pol;
      measurement_sensor[i + sz] = {f, p, std::move(w)};
      measurement_sensor[i + sz].normalize(pol);
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }
  ARTS_USER_ERROR_IF(error.size(), "{}", error)
}
}  // namespace

void measurement_sensorAddGaussianZenith(ArrayOfSensorObsel&    measurement_sensor,
                                         ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                         const AscendingGrid&   freq_grid,
                                         const Vector3&         pos,
                                         const Vector2&         los,
                                         const Stokvec&         pol,
                                         const AscendingGrid&   dzen_grid,
                                         const Numeric&         std_zen) try {
  ARTS_TIME_REPORT

  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  ARTS_USER_ERROR_IF(std_zen <= 0, "Standard deviation must be positive.");
  ARTS_USER_ERROR_IF(dzen_grid.size() == 0, "dzen_grid cannot be empty.");

  Vector      zen_weights(dzen_grid.size());
  const gauss dist(0.0, std_zen);
  for (Size i = 0; i < dzen_grid.size(); ++i) { zen_weights[i] = pdf(dist, dzen_grid[i]); }

  measurement_sensorAddZenithResponse(
      measurement_sensor, std::make_shared<const AscendingGrid>(freq_grid), pos, los, pol, dzen_grid, zen_weights);

  // Metadata: 1D sorted gridded field with frequency grid
  SortedGriddedField1 gf;
  gf.data_name     = "ffts";
  gf.gridname<0>() = "frequency";
  gf.grid<0>()     = freq_grid;
  gf.data.resize(freq_grid.size());
  gf.data = 0.0;
  measurement_sensor_meta.push_back(SensorMetaInfo{std::move(gf)});
}
ARTS_METHOD_ERROR_CATCH

namespace {
template <typename T, typename... Grids>
void measurement_sensorAddRawSensorTmpl(ArrayOfSensorObsel&                         measurement_sensor,
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

  ARTS_USER_ERROR_IF(not raw_sensor.ok() or raw_sensor.data.size() != NF * NZA * NAA * NALT * NLAT * NLON,
                     R"(The raw sensor weights data have an unexpected size or are not OK.

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

    ARTS_USER_ERROR_IF(lowest != -1 and static_cast<Index>(grids[i].i) <= lowest,
                       R"(The raw sensor input grids must be sorted as follows:

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

  Vector f_grid(NF);

  SensorPosLosVector poslos_grid(NZA * NAA * NALT * NLAT * NLON);
  auto               wposlos = poslos_grid.view_as(NZA, NAA, NALT, NLAT, NLON);

  StokvecMatrix w(NZA * NAA * NALT * NLAT * NLON, NF, Stokvec{0, 0, 0, 0});
  auto          wview = w.view_as(NZA, NAA, NALT, NLAT, NLON, NF);

  measurement_sensor.reserve(measurement_sensor.size() + freq_grid.size());

  for (auto f : freq_grid) {
    for (Size ifreq = 0; ifreq < NF; ifreq++) {
      f_grid[ifreq] = f + grids[idf].get(ifreq);

      for (Size iza = 0; iza < NZA; iza++) {
        for (Size iaa = 0; iaa < NAA; iaa++) {
          for (Size ialt = 0; ialt < NALT; ialt++) {
            for (Size ilat = 0; ilat < NLAT; ilat++) {
              for (Size ilon = 0; ilon < NLON; ilon++) {
                if (ifreq == 0) {
                  const Vector3 dpos{grids[idalt].get(ialt), grids[idlat].get(ilat), grids[idlon].get(ilon)};
                  const Vector2 dlos{grids[idza].get(iza), grids[idaa].get(iaa)};

                  wposlos[iza, iaa, ialt, ilat, ilon] = SensorPosLos(pos + dpos, los + dlos);
                }

                wview[iza, iaa, ialt, ilat, ilon, ifreq] = view[ifreq, iza, iaa, ialt, ilat, ilon];
              }
            }
          }
        }
      }
    }

    measurement_sensor.emplace_back(AscendingGrid{f_grid}, poslos_grid, w);
    if (normalize) measurement_sensor.back().normalize();
  }
}
}  // namespace

#define AddRawSensor(T, ...)                                                                                  \
  void measurement_sensorAddRawSensor(ArrayOfSensorObsel&                            measurement_sensor,      \
                                      ArrayOfSensorMetaInfo&                         measurement_sensor_meta, \
                                      const AscendingGrid&                           freq_grid,               \
                                      const Vector3&                                 pos,                     \
                                      const Vector2&                                 los,                     \
                                      const matpack::gridded_data_t<T, __VA_ARGS__>& raw_sensor,              \
                                      const Index&                                   normalize) try {         \
    ARTS_TIME_REPORT                                                                                          \
                                                                                                              \
    measurement_sensorAddRawSensorTmpl(                                                                       \
        measurement_sensor, freq_grid, pos, los, raw_sensor, static_cast<bool>(normalize));                   \
                                                                                                              \
    SortedGriddedField1 gf;                                                                                   \
    gf.data_name     = "raw";                                                                                 \
    gf.gridname<0>() = "frequency";                                                                           \
    gf.grid<0>()     = freq_grid;                                                                             \
    gf.data.resize(freq_grid.size());                                                                         \
    gf.data = 0.0;                                                                                            \
    measurement_sensor_meta.push_back(SensorMetaInfo{std::move(gf)});                                         \
  }                                                                                                           \
  ARTS_METHOD_ERROR_CATCH

AddRawSensor(Numeric, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);

#undef AddRawSensor

void measurement_sensorAddCamera(ArrayOfSensorObsel&    measurement_sensor,
                                 ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                 const AscendingGrid&   freq_grid,
                                 const Vector3&         pos,
                                 const Vector2&         los,
                                 const Stokvec&         pol,
                                 const Index&           n_h,
                                 const Index&           n_w,
                                 const Numeric&         ccd_h,
                                 const Numeric&         ccd_w,
                                 const Numeric&         focal_length,
                                 const Numeric&         aperture_diameter,
                                 const Numeric&         focus_distance) try {
  ARTS_TIME_REPORT

  const Size n_freq    = freq_grid.size();
  const Size n_pixels  = static_cast<Size>(n_h) * static_cast<Size>(n_w);
  const Size start_idx = measurement_sensor.size();

  ARTS_USER_ERROR_IF(n_freq == 0, "freq_grid must not be empty")
  ARTS_USER_ERROR_IF(n_h < 1, "n_h must be at least 1, got {}", n_h)
  ARTS_USER_ERROR_IF(n_w < 1, "n_w must be at least 1, got {}", n_w)
  ARTS_USER_ERROR_IF(ccd_h <= 0, "ccd_h must be positive, got {}", ccd_h)
  ARTS_USER_ERROR_IF(ccd_w <= 0, "ccd_w must be positive, got {}", ccd_w)
  ARTS_USER_ERROR_IF(focal_length <= 0, "focal_length must be positive, got {}", focal_length)
  ARTS_USER_ERROR_IF(aperture_diameter <= 0, "aperture_diameter must be positive, got {}", aperture_diameter)
  ARTS_USER_ERROR_IF(focus_distance <= 0, "focus_distance must be positive, got {}", focus_distance)
  ARTS_USER_ERROR_IF(pol.nonzero_components() == 0, "pol is 0")

  // Pixel pitch (distance between pixel centers)
  const Numeric pitch_h = ccd_h / static_cast<Numeric>(n_h);
  const Numeric pitch_w = ccd_w / static_cast<Numeric>(n_w);

  // Thin-lens adjusted image distance for the given focus distance.
  // Thin lens equation: 1/f = 1/d_o + 1/d_i
  //   => d_i = f * d_o / (d_o - f)
  ARTS_USER_ERROR_IF(focus_distance <= focal_length,
                     "focus_distance ({}) must be greater than focal_length ({})",
                     focus_distance,
                     focal_length)

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

  // Shared frequency grid
  auto f = std::make_shared<const AscendingGrid>(freq_grid);

  // Build a SINGLE shared poslos_grid containing all pixel directions.
  // collect_simulations() deduplicates by f_grid pointer and iterates
  // the (shared) poslos_grid.  All camera obsels must share the same
  // poslos_grid so that every pixel direction gets simulated.
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
  auto p = std::make_shared<const SensorPosLosVector>(std::move(all_poslos));

  // Build obsels.  Each pixel (ih,iw) at frequency ifreq has a sparse
  // weight matrix of shape (n_pixels, n_freq) with a single nonzero
  // entry at (pix_idx, ifreq).
  measurement_sensor.resize(start_idx + n_pixels * n_freq);

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Index ih = 0; ih < n_h; ++ih) {
    for (Index iw = 0; iw < n_w; ++iw) {
      try {
        const Size pix_idx = static_cast<Size>(ih) * static_cast<Size>(n_w) + static_cast<Size>(iw);

        for (Size ifreq = 0; ifreq < n_freq; ++ifreq) {
          sensor::SparseStokvecMatrix w(n_pixels, n_freq);
          w[pix_idx, ifreq]                                        = pol;
          measurement_sensor[start_idx + pix_idx * n_freq + ifreq] = {f, p, std::move(w)};
        }
      } catch (std::runtime_error& e) {
#pragma omp critical
        if (error.empty()) error = e.what();
      }
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "{}", error)

  // Build camera metadata: gridded_data_t<Numeric, AscendingGrid, AscendingGrid, Vector>
  // grid<0> = row angular offsets from boresight (degrees)
  // grid<1> = col angular offsets from boresight (degrees)
  // grid<2> = frequency / color
  sensor::CameraGriddedField gf;
  gf.data_name     = "camera";
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

  // Zero data — to be filled later by measurement_sensor_metaFromMeasurementVec
  gf.data.resize(n_h, n_w, n_freq);
  gf.data = 0.0;

  measurement_sensor_meta.push_back(SensorMetaInfo{std::move(gf)});
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensor_metaFromMeasurementVec(ArrayOfSensorMetaInfo& measurement_sensor_meta,
                                               const Vector&          measurement_vec) try {
  ARTS_TIME_REPORT

  // Compute the total element count across all meta entries.
  Size total = 0;
  for (const auto& meta : measurement_sensor_meta) total += meta.count();

  ARTS_USER_ERROR_IF(total != static_cast<Size>(measurement_vec.size()),
                     "Total meta count ({}) does not match measurement_vec size ({})",
                     total,
                     measurement_vec.size())

  Size start = 0;
  for (auto& meta : measurement_sensor_meta) {
    const Size n = meta.count();

    std::visit(
        [&](auto& gf) {
          // Flat copy: measurement_vec[start..start+n) -> gf.data elements
          for (Size i = 0; i < n; ++i) gf.data.elem_at(static_cast<Index>(i)) = measurement_vec[start + i];
        },
        meta.data);

    start += n;
  }
}
ARTS_METHOD_ERROR_CATCH
