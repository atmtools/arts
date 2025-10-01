#include <arts_omp.h>
#include <debug.h>
#include <obsel.h>
#include <rtepack.h>
#include <workspace.h>

#include <boost/math/distributions/normal.hpp>
#include <stdexcept>

void measurement_sensorInit(ArrayOfSensorObsel& measurement_sensor) {
  ARTS_TIME_REPORT

  measurement_sensor = {};
}

void measurement_sensorAddSimple(ArrayOfSensorObsel& measurement_sensor,
                                 const AscendingGrid& frequency_grid,
                                 const Vector3& pos,
                                 const Vector2& los,
                                 const Stokvec& pol) try {
  ARTS_TIME_REPORT

  const Index n = frequency_grid.size();
  const Size sz = measurement_sensor.size();

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  StokvecMatrix w(1, n, 0);
  for (Index i = 0; i < n; i++) {
    w[0, i]                    = pol;
    measurement_sensor[i + sz] = {f, p, w};
    w[0, i]                    = 0.0;
  }
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorAddVectorGaussian(ArrayOfSensorObsel& measurement_sensor,
                                         const AscendingGrid& frequency_grid,
                                         const Vector& stds,
                                         const Vector3& pos,
                                         const Vector2& los,
                                         const Stokvec& pol) try {
  ARTS_TIME_REPORT

  using gauss = boost::math::normal_distribution<Numeric>;
  using boost::math::pdf;

  const Size n       = frequency_grid.size();
  const Size sz      = measurement_sensor.size();
  const Size nonzero = pol.nonzero_components();

  ARTS_USER_ERROR_IF(n < 2, "Must have a frequency grid")
  ARTS_USER_ERROR_IF(n != stds.size(),
                     "Must have a standard deviation for each frequency point")
  ARTS_USER_ERROR_IF(stdr::any_of(stds, Cmp::le(0)),
                     "No negative standard deviation not allowed.\nstds := {:B,}", stds)
  ARTS_USER_ERROR_IF(nonzero == 0, "pol is 0")

  measurement_sensor.resize(sz + n);

  auto f = std::make_shared<const AscendingGrid>(frequency_grid);
  auto p = std::make_shared<const SensorPosLosVector>(
      SensorPosLosVector{{.pos = pos, .los = los}});

  String error;
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < n; i++) {
    try {
      StokvecMatrix w(1, n, pol);

      const gauss dist(frequency_grid[i], stds[i]);
      for (Size j = 0; j < n; j++) {
        w[0, j] *= pdf(dist, frequency_grid[j]);
      }

      measurement_sensor[i + sz] = {f, p, w};
      measurement_sensor[i + sz].normalize(pol);
    } catch (std::runtime_error& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "{}", error)
}
ARTS_METHOD_ERROR_CATCH

void measurement_sensorAddSimpleGaussian(ArrayOfSensorObsel& measurement_sensor,
                                         const AscendingGrid& frequency_grid,
                                         const Numeric& std,
                                         const Vector3& pos,
                                         const Vector2& los,
                                         const Stokvec& pol) {
  ARTS_TIME_REPORT

  const Vector stds(frequency_grid.size(), std);
  measurement_sensorAddVectorGaussian(
      measurement_sensor, frequency_grid, stds, pos, los, pol);
}

namespace {
std::span<SensorObsel> get_span(ArrayOfSensorObsel& measurement_sensor,
                                const Index& start,
                                const Index& end) {
  ARTS_USER_ERROR_IF(
      start > end,
      "Start index ({}) must not be greater than end index ({}).",
      start,
      end)

  ARTS_USER_ERROR_IF(
      start > static_cast<Index>(measurement_sensor.size()),
      "Start index ({}) must not be greater than the size of the measurement sensor ({}).",
      start,
      measurement_sensor.size())

  ARTS_USER_ERROR_IF(
      end > static_cast<Index>(measurement_sensor.size()),
      "End index ({}) must not be greater than the size of the measurement sensor ({}).",
      end,
      measurement_sensor.size())

  auto first = measurement_sensor.begin() + (start > 0) * start;
  auto last =
      (end < 0) ? measurement_sensor.end() : measurement_sensor.begin() + end;

  return {first, last};
}
}  // namespace

void measurement_sensorMakeExclusive(ArrayOfSensorObsel& measurement_sensor,
                                     const Index& start,
                                     const Index& end) try {
  ARTS_TIME_REPORT

  make_exclusive(get_span(measurement_sensor, start, end));
}
ARTS_METHOD_ERROR_CATCH

namespace {

template <typename T, typename... Grids>
void measurement_sensorAddRawSensorTmpl(
    ArrayOfSensorObsel& measurement_sensor,
    const AscendingGrid& frequency_grid,
    const Vector3& pos,
    const Vector2& los,
    const matpack::gridded_data_t<T, Grids...>& raw_sensor,
    const bool normalize)
  requires(sizeof...(Grids) <= 6 and
           (std::same_as<Grids, AscendingGrid> and ...) and
           (std::same_as<Numeric, T> or std::same_as<Stokvec, T>))
{
  constexpr Size idf   = 0;
  constexpr Size idza  = 1;
  constexpr Size idaa  = 2;
  constexpr Size idalt = 3;
  constexpr Size idlat = 4;
  constexpr Size idlon = 5;

  struct InternalGrid {
    const AscendingGrid* grid{nullptr};
    Size i{};

    [[nodiscard]] bool exists() const { return grid != nullptr; }
    [[nodiscard]] Size size() const { return exists() ? grid->size() : 1; }
    [[nodiscard]] Numeric get(Size i) const {
      return exists() ? (*grid)[i] : 0.0;
    }

    InternalGrid() = default;
    InternalGrid(const char* name,
                 const matpack::gridded_data_t<T, Grids...>& raw_sensor) {
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
  grids[idf]   = InternalGrid("df", raw_sensor);
  grids[idza]  = InternalGrid("dza", raw_sensor);
  grids[idaa]  = InternalGrid("daa", raw_sensor);
  grids[idalt] = InternalGrid("dalt", raw_sensor);
  grids[idlat] = InternalGrid("dlat", raw_sensor);
  grids[idlon] = InternalGrid("dlon", raw_sensor);

  const Size NF   = grids[idf].size();
  const Size NZA  = grids[idza].size();
  const Size NAA  = grids[idaa].size();
  const Size NALT = grids[idalt].size();
  const Size NLAT = grids[idlat].size();
  const Size NLON = grids[idlon].size();

  ARTS_USER_ERROR_IF(
      not raw_sensor.ok() or
          raw_sensor.data.size() != NF * NZA * NAA * NALT * NLAT * NLON,
      R"(The raw sensor weights data have an unexpected size or are not OK.

The internal method wants to view the data as a 6-dimensional
array of shape:  [{0:B}, {1:}, {2:}, {3:}, {4:}, {5:}]
where the dimensions are:
  Frequency perturbance     ("df"):   {0:}
  Zenith angle perturbance  ("dza"):  {1:}
  Azimuth angle perturbance ("daa"):  {2:}
  Altitude perturbance      ("dalt"): {3:}
  Latitude perturbance      ("dlat"): {4:}
  Longitude perturbance     ("dlon"): {5:}

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

    ARTS_USER_ERROR_IF(
        lowest != -1 and static_cast<Index>(grids[i].i) <= lowest,
        R"(The raw sensor input grids must be sorted as follows:

Frequency grid perturbance (named "df"),
Zenith angle perturbance   (named "dza"),
Azimuth angle perturbance  (named "daa"),
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
  auto wposlos = poslos_grid.view_as(NZA, NAA, NALT, NLAT, NLON);

  StokvecMatrix w(NZA * NAA * NALT * NLAT * NLON, NF, Stokvec{0, 0, 0, 0});
  auto wview = w.view_as(NZA, NAA, NALT, NLAT, NLON, NF);

  measurement_sensor.reserve(measurement_sensor.size() + frequency_grid.size());

  for (auto f : frequency_grid) {
    for (Size ifreq = 0; ifreq < NF; ifreq++) {
      f_grid[ifreq] = f + grids[idf].get(ifreq);

      for (Size iza = 0; iza < NZA; iza++) {
        for (Size iaa = 0; iaa < NAA; iaa++) {
          for (Size ialt = 0; ialt < NALT; ialt++) {
            for (Size ilat = 0; ilat < NLAT; ilat++) {
              for (Size ilon = 0; ilon < NLON; ilon++) {
                if (ifreq == 0) {
                  const Vector3 dpos{grids[idalt].get(ialt),
                                     grids[idlat].get(ilat),
                                     grids[idlon].get(ilon)};
                  const Vector2 dlos{grids[idza].get(iza),
                                     grids[idaa].get(iaa)};

                  wposlos[iza, iaa, ialt, ilat, ilon] =
                      SensorPosLos(pos + dpos, los + dlos);
                }

                wview[iza, iaa, ialt, ilat, ilon, ifreq] =
                    view[ifreq, iza, iaa, ialt, ilat, ilon];
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

#define AddRawSensor(T, ...)                                          \
  void measurement_sensorAddRawSensor(                                \
      ArrayOfSensorObsel& measurement_sensor,                         \
      const AscendingGrid& frequency_grid,                            \
      const Vector3& pos,                                             \
      const Vector2& los,                                             \
      const matpack::gridded_data_t<T, __VA_ARGS__>& raw_sensor,      \
      const Index& normalize) try {                                   \
    ARTS_TIME_REPORT                                                  \
                                                                      \
    measurement_sensorAddRawSensorTmpl(measurement_sensor,            \
                                       frequency_grid,                \
                                       pos,                           \
                                       los,                           \
                                       raw_sensor,                    \
                                       static_cast<bool>(normalize)); \
  }                                                                   \
  ARTS_METHOD_ERROR_CATCH

AddRawSensor(Numeric, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid);
AddRawSensor(Numeric, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(
    Numeric, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Numeric,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid);
AddRawSensor(Numeric,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(
    Stokvec, AscendingGrid, AscendingGrid, AscendingGrid, AscendingGrid);
AddRawSensor(Stokvec,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid);
AddRawSensor(Stokvec,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid,
             AscendingGrid);

#undef AddRawSensor
