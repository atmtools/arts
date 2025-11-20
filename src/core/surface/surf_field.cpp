#include "surf_field.h"

#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <debug.h>
#include <enumsSurfaceKey.h>
#include <lagrange_interp.h>

#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>

Surf::Point::Point()                                   = default;
Surf::Point::Point(const Point &)                      = default;
Surf::Point::Point(Point &&) noexcept                  = default;
Surf::Point &Surf::Point::operator=(const Point &)     = default;
Surf::Point &Surf::Point::operator=(Point &&) noexcept = default;
Surf::Data::Data()                                     = default;
Surf::Data::Data(const Data &)                         = default;
Surf::Data::Data(Data &&) noexcept                     = default;
Surf::Data &Surf::Data::operator=(const Data &)        = default;
Surf::Data &Surf::Data::operator=(Data &&) noexcept    = default;
Surf::Field::Field()                                   = default;
Surf::Field::Field(const Field &)                      = default;
Surf::Field::Field(Field &&) noexcept                  = default;
Surf::Field &Surf::Field::operator=(const Field &)     = default;
Surf::Field &Surf::Field::operator=(Field &&) noexcept = default;

Numeric &Surf::Point::operator[](SurfaceKey x) {
  switch (x) {
    case SurfaceKey::h: return elevation;
    case SurfaceKey::t: return temperature;
  }
  std::unreachable();
}

Numeric &Surf::Point::operator[](const SurfacePropertyTag &x) {
  return prop[x];
}

Numeric Surf::Point::operator[](SurfaceKey x) const {
  switch (x) {
    case SurfaceKey::h: return elevation;
    case SurfaceKey::t: return temperature;
  }
  std::unreachable();
}

Numeric Surf::Point::operator[](const SurfacePropertyTag &x) const {
  return prop.at(x);
}

namespace {
Vector3 to_xyz(Vector2 ell, Vector3 pos) {
  using Conversion::cosd;
  using Conversion::sind;
  using Math::pow2;

  const auto [a, b]          = ell;
  const auto [alt, lat, lon] = pos;

  const Numeric e = std::sqrt(1 - pow2(b / a));
  const Numeric N = a / std::sqrt(1 - pow2(e * sind(lat)));

  return {(N + alt) * cosd(lon) * cosd(lat),
          pos[1] = (N + alt) * sind(lon) * cosd(lat),
          pos[2] = (N * (1 - pow2(e)) + alt) * sind(lat)};
}

Vector2 from_xyz_dxyz(Vector3 xyz, Vector3 dxyz) {
  using Conversion::acosd;
  using Conversion::asind;
  using Conversion::atan2d;
  using Conversion::cosd;
  using Conversion::sind;
  using std::hypot;

  const auto [x, y, z]    = xyz;
  const auto [dx, dy, dz] = dxyz;

  const Numeric r          = hypot(z, y, x);
  const Vector3 rll        = {r, asind(z / r), atan2d(y, x)};
  const auto [_, lat, lon] = rll;

  const auto norm = hypot(dz, dy, dx);

  const Numeric slat = sind(lat);
  const Numeric clat = cosd(lat);
  const Numeric slon = sind(lon);
  const Numeric clon = cosd(lon);

  const Numeric dr = (clat * clon * dx + slat * dz + clat * slon * dy) / norm;
  const Numeric dlat =
      (-slat * clon * dx + clat * dz - slat * slon * dy) / (norm * r);
  const Numeric dlon = (-slon / clat * dx + clon / clat * dy) / norm / r;

  Vector2 los;
  los[0] = acosd(dr);
  los[1] = acosd(r * dlat / sind(los[0]));
  if (std::isnan(los[1])) {
    if (dlat >= 0) {
      los[1] = 0;
    } else {
      los[1] = 180;
    }
  } else if (dlon < 0) {
    los[1] = -los[1];
  }

  return los;
}

std::ostream &operator<<(std::ostream &os, const SurfaceKeyVal &key) {
  std::visit([&](auto &k) { os << k; }, key);
  return os;
}
}  // namespace

namespace Surf {
String Data::data_type() const {
  if (std::holds_alternative<GeodeticField2>(data)) return "GeodeticField2";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";

  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

[[nodiscard]] std::vector<SurfaceKeyVal> Point::keys() const {
  std::vector<SurfaceKeyVal> out;
  out.reserve(enumtyps::SurfaceKeyTypes.size());
  for (auto key : enumtyps::SurfaceKeyTypes) out.emplace_back(key);
  return out;
}

[[nodiscard]] Index Point::size() const { return nother(); }

Numeric &Point::operator[](const SurfaceKeyVal &x) {
  return std::visit(
      [this](auto &key) -> Numeric & {
        return const_cast<Point *>(this)->operator[](key);
      },
      x);
}

Numeric Point::operator[](const SurfaceKeyVal &x) const {
  return std::visit(
      [this](auto &key) -> Numeric { return this->operator[](key); }, x);
}

namespace {
Numeric numeric_interpolation(
    const GeodeticField2 &data,
    Numeric lat,
    Numeric lon,
    std::pair<InterpolationExtrapolation, InterpolationExtrapolation>
        lat_extrap,
    std::pair<InterpolationExtrapolation, InterpolationExtrapolation>
        lon_extrap) {
  if (not data.ok()) throw std::runtime_error("bad field");

  const Vector &lats = data.grid<0>();
  const Vector &lons = data.grid<1>();

  if (lat < lats.front()) {
    ARTS_USER_ERROR_IF(lat_extrap.first == InterpolationExtrapolation::None,
                       "No extrapolation allowed")
    if (lat_extrap.first == InterpolationExtrapolation::Zero) return 0.0;
    if (lat_extrap.first == InterpolationExtrapolation::Nearest)
      lat = lats.front();
  }

  if (lat > lats.back()) {
    ARTS_USER_ERROR_IF(lat_extrap.second == InterpolationExtrapolation::None,
                       "No extrapolation allowed")
    if (lat_extrap.second == InterpolationExtrapolation::Zero) return 0.0;
    if (lat_extrap.second == InterpolationExtrapolation::Nearest)
      lat = lats.back();
  }

  if (lat < lons.front()) {
    ARTS_USER_ERROR_IF(lon_extrap.first == InterpolationExtrapolation::None,
                       "No extrapolation allowed")
    if (lon_extrap.first == InterpolationExtrapolation::Zero) return 0.0;
    if (lon_extrap.first == InterpolationExtrapolation::Nearest)
      lat = lons.front();
  }

  if (lat > lons.back()) {
    ARTS_USER_ERROR_IF(lon_extrap.second == InterpolationExtrapolation::None,
                       "No extrapolation allowed")
    if (lon_extrap.second == InterpolationExtrapolation::Zero) return 0.0;
    if (lon_extrap.second == InterpolationExtrapolation::Nearest)
      lat = lons.back();
  }

  if (lats.size() == 1 and lons.size() == 1) {
    return data.data[0, 0];
  }

  if (lats.size() == 1) {
    using LatLag = lagrange_interp::lag_t<0>;
    using LonLag = lagrange_interp::lag_t<1, lagrange_interp::loncross>;

    return interp(data.data,
                  LatLag(lats, lat, lagrange_interp::ascending_grid_t{}),
                  LonLag(lons, lon, lagrange_interp::ascending_grid_t{}));
  }

  if (lons.size() == 1) {
    using LatLag = lagrange_interp::lag_t<1>;
    using LonLag = lagrange_interp::lag_t<0, lagrange_interp::loncross>;

    return interp(data.data,
                  LatLag(lats, lat, lagrange_interp::ascending_grid_t{}),
                  LonLag(lons, lon, lagrange_interp::ascending_grid_t{}));
  }

  {
    using LatLag = lagrange_interp::lag_t<1>;
    using LonLag = lagrange_interp::lag_t<1, lagrange_interp::loncross>;

    return interp(data.data,
                  LatLag(lats, lat, lagrange_interp::ascending_grid_t{}),
                  LonLag(lons, lon, lagrange_interp::ascending_grid_t{}));
  }
}

Numeric numeric_interpolation(
    const FunctionalData &f,
    Numeric lat,
    Numeric lon,
    std::pair<InterpolationExtrapolation, InterpolationExtrapolation>,
    std::pair<InterpolationExtrapolation, InterpolationExtrapolation>) {
  return f(lat, lon);
}

constexpr Numeric numeric_interpolation(
    Numeric x,
    Numeric,
    Numeric,
    std::pair<InterpolationExtrapolation, InterpolationExtrapolation>,
    std::pair<InterpolationExtrapolation, InterpolationExtrapolation>) {
  return x;
}

auto interpolation_function(Numeric lat, Numeric lon) {
  return [lat, lon](auto &data) {
    const auto call = [lat,
                       lon,
                       lat_extrap = std::pair{data.lat_low, data.lat_upp},
                       lon_extrap =
                           std::pair{data.lon_low, data.lon_upp}](auto &&x) {
      return numeric_interpolation(x, lat, lon, lat_extrap, lon_extrap);
    };

    return std::visit(call, data.data);
  };
}
}  // namespace

Data &Field::operator[](const SurfaceKey &key) { return other[key]; }

Data &Field::operator[](const SurfacePropertyTag &key) { return props[key]; }

Data &Field::operator[](const KeyVal &key) {
  return std::visit([this](auto &k) -> Data & { return this->operator[](k); },
                    key);
}

const Data &Field::operator[](const SurfaceKey &key) const {
  return other.at(key);
}

const Data &Field::operator[](const SurfacePropertyTag &key) const {
  return props.at(key);
}

const Data &Field::operator[](const KeyVal &key) const {
  return std::visit(
      [this](auto &k) -> const Data & { return this->operator[](k); }, key);
}

bool Field::contains(const SurfaceKey &key) const {
  return other.contains(key);
}

bool Field::contains(const SurfacePropertyTag &key) const {
  return props.contains(key);
}

bool Field::contains(const KeyVal &key) const {
  return std::visit([this](auto &k) -> bool { return this->contains(k); }, key);
}

Vector2 Field::normal(Numeric lat, Numeric lon, Numeric alt) const try {
  ARTS_USER_ERROR_IF(
      bad_ellipsoid(), "Ellipsoid must have positive axes: {:B,}", ellipsoid)

  constexpr Vector2 up{180, 0};

  if (not contains(SurfaceKey::h)) {
    return up;
  }

  const auto &z = other.at(SurfaceKey::h);

  if (std::holds_alternative<Numeric>(z.data)) {
    return up;
  }

  alt = std::isnan(alt) ? interpolation_function(lat, lon)(z) : alt;

  constexpr Numeric offset = 1e-3;
  const Numeric lat1 = (lat + offset > 90) ? (lat - offset) : (lat + offset),
                lon1 = lon, alt1 = interpolation_function(lat1, lon1)(z);
  const Numeric lat2 = lat, lon2 = lon + offset,
                alt2 = interpolation_function(lat2, lon2)(z);

  const Vector3 xyz  = to_xyz(ellipsoid, {alt, lat, lon});
  const Vector3 xyz1 = to_xyz(ellipsoid, {alt1, lat1, lon1});
  const Vector3 xyz2 = to_xyz(ellipsoid, {alt2, lat2, lon2});
  const Vector3 r1{xyz1[0] - xyz[0], xyz1[1] - xyz[1], xyz1[2] - xyz[2]};
  const Vector3 r2{xyz2[0] - xyz[0], xyz2[1] - xyz[1], xyz2[2] - xyz[2]};

  const Vector3 dxyz{r1[1] * r2[2] - r1[2] * r2[1],
                     r1[2] * r2[0] - r1[0] * r2[2],
                     r1[0] * r2[1] - r1[1] * r2[0]};

  return from_xyz_dxyz(xyz, dxyz);
} catch (std::exception &e) {
  throw std::runtime_error(std::format(
      "Cannot find a normal to the surface at position {} {}\nThe internal error reads: {}",
      lat,
      lon,
      std::string_view(e.what())));
}

Numeric Data::at(const Numeric lat, const Numeric lon) const {
  return interpolation_function(lat, lon)(*this);
}

Size Field::nprops() const { return props.size(); }
Size Field::nother() const { return other.size(); }
Size Field::size() const { return nprops() + nother(); }

std::vector<Field::KeyVal> Field::keys() const {
  std::vector<KeyVal> out;
  out.reserve(size());
  for (const auto &kv : other) out.emplace_back(kv.first);
  for (const auto &kv : props) out.emplace_back(kv.first);
  return out;
}

Point Field::at(Numeric lat, Numeric lon) const {
  Point out;

  for (auto &key : keys()) {
    out[key] = this->operator[](key).at(lat, lon);
  }

  out.normal = normal(lat, lon);

  return out;
}

Numeric Field::single_value(const KeyVal &key, Numeric lat, Numeric lon) const {
  ARTS_USER_ERROR_IF(
      not contains(key), "Surface field does not possess the key: {}", key)

  const auto interp = interpolation_function(lat, lon);

  return interp(this->operator[](key));
}

namespace {
constexpr std::pair<Numeric, Numeric> minmax(Numeric x) { return {x, x}; }

std::pair<Numeric, Numeric> minmax(const FunctionalData &) {
  ARTS_USER_ERROR("Cannot extract minmax from functional data");
  return {std::numeric_limits<Numeric>::lowest(),
          std::numeric_limits<Numeric>::max()};
}

std::pair<Numeric, Numeric> minmax(const GeodeticField2 &x) {
  return matpack::minmax(x.data);
}
}  // namespace

std::pair<Numeric, Numeric> Field::minmax_single_value(
    const KeyVal &key) const {
  ARTS_USER_ERROR_IF(
      not contains(key), "Surface field does not possess the key: {}", key)
  return std::visit([](auto &a) { return minmax(a); },
                    this->operator[](key).data);
}

bool Field::constant_value(const KeyVal &key) const {
  ARTS_USER_ERROR_IF(
      not contains(key), "Surface field does not possess the key: {}", key)

  return std::holds_alternative<Numeric>(this->operator[](key).data);
}

[[nodiscard]] ConstVectorView Data::flat_view() const {
  return std::visit(
      [](auto &X) -> ConstVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, GeodeticField2>)
          return X.data.view_as(X.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return ConstVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ConstVectorView{};
        assert(false);
      },
      data);
}

[[nodiscard]] VectorView Data::flat_view() {
  return std::visit(
      [](auto &X) -> VectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, GeodeticField2>)
          return X.data.view_as(X.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return VectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return VectorView{};
        assert(false);
      },
      data);
}

namespace {
std::array<std::pair<Index, Numeric>, 4> flat_weights_(const Numeric &,
                                                       const Numeric &,
                                                       const Numeric &) {
  constexpr auto v1 = std::pair<Index, Numeric>{0, 1.};
  constexpr auto v0 = std::pair<Index, Numeric>{0, 0.};
  return {v1, v0, v0, v0};
}

std::array<std::pair<Index, Numeric>, 4> flat_weights_(const FunctionalData &,
                                                       const Numeric &,
                                                       const Numeric &) {
  constexpr auto v0 = std::pair<Index, Numeric>{0, 0.};
  return {v0, v0, v0, v0};
}

std::array<std::pair<Index, Numeric>, 4> flat_weights_(const GeodeticField2 &v,
                                                       const Numeric &lat,
                                                       const Numeric &lon) {
  using LonLag = lagrange_interp::lag_t<1, lagrange_interp::loncross>;
  using LatLag = lagrange_interp::lag_t<1>;

  const auto slon = v.shape()[1];
  const bool d1   = v.shape()[0] == 1;
  const bool d2   = slon == 1;

  constexpr auto v1 = std::pair<Index, Numeric>{0, 1.};
  constexpr auto v0 = std::pair<Index, Numeric>{0, 0.};

  if (d1 and d2) {
    return {v1, v0, v0, v0};
  }

  if (d1) {
    const LonLag wlon(v.grid<1>(), lon, lagrange_interp::ascending_grid_t{});
    return {std::pair<Index, Numeric>{wlon.indx[0], wlon.data[0]},
            std::pair<Index, Numeric>{wlon.indx[1], wlon.data[1]},
            v0,
            v0};
  }

  if (d2) {
    const LatLag wlat(v.grid<0>(), lat, lagrange_interp::ascending_grid_t{});
    return {std::pair<Index, Numeric>{wlat.indx[0], wlat.data[0]},
            std::pair<Index, Numeric>{wlat.indx[1], wlat.data[1]},
            v0,
            v0};
  }

  const LatLag wlat(v.grid<0>(), lat, lagrange_interp::ascending_grid_t{});
  const LonLag wlon(v.grid<1>(), lon, lagrange_interp::ascending_grid_t{});
  auto sz = [slon](auto lat_pos, auto lon_pos) {
    return lat_pos * slon + lon_pos;
  };

  return {std::pair<Index, Numeric>{sz(wlat.indx[0], wlon.indx[0]),
                                    wlat.data[0] * wlon.data[0]},
          std::pair<Index, Numeric>{sz(wlat.indx[0], wlon.indx[1]),
                                    wlat.data[0] * wlon.data[1]},
          std::pair<Index, Numeric>{sz(wlat.indx[1], wlon.indx[0]),
                                    wlat.data[1] * wlon.data[0]},
          std::pair<Index, Numeric>{sz(wlat.indx[1], wlon.indx[1]),
                                    wlat.data[1] * wlon.data[1]}};
}
}  // namespace

//! Flat weights for the positions in an atmosphere
std::array<std::pair<Index, Numeric>, 4> Data::flat_weights(
    const Numeric &lat, const Numeric &lon) const {
  return std::visit([&](auto &v) { return flat_weights_(v, lat, lon); }, data);
}

bool Data::ok() const {
  if (std::holds_alternative<GeodeticField2>(data)) {
    auto &v = *std::get_if<GeodeticField2>(&data);
    return v.ok() and
           lagrange_interp::loncross::cycle(v.grid<1>().front()) ==
               v.grid<1>().front() and
           lagrange_interp::loncross::cycle(v.grid<1>().back()) ==
               v.grid<1>().back() and
           v.grid<0>().front() >= -90 and v.grid<0>().back() <= 90;
  }

  return true;
}

void Data::adjust_interpolation_extrapolation() {
  if (std::holds_alternative<GeodeticField2>(data)) {
    auto &field = std::get<GeodeticField2>(data);

    if (field.grid<0>().size() == 1) {
      lat_upp = InterpolationExtrapolation::Nearest;
      lat_low = InterpolationExtrapolation::Nearest;
    }

    if (field.grid<1>().size() == 1) {
      lon_upp = InterpolationExtrapolation::Nearest;
      lon_low = InterpolationExtrapolation::Nearest;
    }
  } else {
    lat_upp = InterpolationExtrapolation::Nearest;
    lat_low = InterpolationExtrapolation::Nearest;
    lon_upp = InterpolationExtrapolation::Nearest;
    lon_low = InterpolationExtrapolation::Nearest;
  }
}

Data::Data(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

Data::Data(GeodeticField2 x) : data(std::move(x)) {
  adjust_interpolation_extrapolation();
}

Data::Data(FunctionalData x) : data(std::move(x)) {
  adjust_interpolation_extrapolation();
}

Data &Data::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

Data &Data::operator=(GeodeticField2 x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

Data &Data::operator=(FunctionalData x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

[[nodiscard]] bool Field::bad_ellipsoid() const {
  return not(ellipsoid[1] > 0 and ellipsoid[0] >= ellipsoid[1]);
}

bool Point::contains(const SurfaceKeyVal &k) const {
  return std::visit([this](auto &key) -> bool { return this->contains(key); },
                    k);
}
}  // namespace Surf

bool operator==(const SurfaceKeyVal &lhs, SurfaceKey rhs) {
  auto *val = std::get_if<SurfaceKey>(&lhs);

  return val ? (*val == rhs) : false;
}

bool operator==(SurfaceKey lhs, const SurfaceKeyVal &rhs) { return rhs == lhs; }

bool operator==(const SurfaceKeyVal &lhs, const SurfacePropertyTag &rhs) {
  auto *val = std::get_if<SurfacePropertyTag>(&lhs);

  return val ? (*val == rhs) : false;
}

bool operator==(const SurfacePropertyTag &lhs, const SurfaceKeyVal &rhs) {
  return rhs == lhs;
}

static_assert(
    std::same_as<typename SurfaceField::KeyVal, SurfaceKeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");
