#include "surf.h"

#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <variant>

#include "arts_constexpr_math.h"
#include "arts_conversions.h"
#include "debug.h"
#include "interp.h"
#include "matpack_math.h"

template <>
const Surf::Data &FieldMap::
    Map<Surf::Data, SurfaceKey, SurfaceTypeTag, SurfacePropertyTag>::operator[](
        const KeyVal &k) const try {
  return std::visit(
      [this](auto &key) -> const Surf::Data & {
        return this->map<decltype(key)>().at(key);
      },
      k);
} catch (std::out_of_range &) {
  throw std::out_of_range(var_string("Key not found in map: \"", k, '\"'));
} catch (...) {
  throw;
}

template <>
Surf::Data &FieldMap::
    Map<Surf::Data, SurfaceKey, SurfaceTypeTag, SurfacePropertyTag>::operator[](
        const KeyVal &k) try {
  return std::visit(
      [this](auto &key) -> Surf::Data & {
        return const_cast<Map *>(this)->map<decltype(key)>()[key];
      },
      k);
}
ARTS_METHOD_ERROR_CATCH

template <>
bool FieldMap::Map<Surf::Data, SurfaceKey, SurfaceTypeTag, SurfacePropertyTag>::
    contains(const KeyVal &key) const {
  return std::visit(
      [this](auto &k) -> bool { return this->map<decltype(k)>().contains(k); },
      key);
}

namespace geodetic {
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
}  // namespace geodetic

std::ostream &operator<<(std::ostream &os, const SurfaceKeyVal &key) {
  std::visit([&](auto &k) { os << k; }, key);
  return os;
}

namespace Surf {
std::ostream &operator<<(std::ostream &os, const Point &surf) {
  os << "Elevation: " << surf.elevation << " m,\n";
  os << "Temperature: " << surf.temperature << " K\n,";
  os << "Wind Field: [u: " << surf.wind[0] << ", v: " << surf.wind[1]
     << ", w: " << surf.wind[2] << "] m/s,\n";
  os << "Normal: [za: " << surf.normal[0] << ", aa: " << surf.normal[1]
     << "] degrees,\n";

  if (surf.type.size()) {
    os << "Types: [";
    for (auto &type : surf.type)
      os << type.first << ": " << type.second << ", ";
    return os << "] ratio";
  }
  return os << "Unspecified surface type";
}

std::ostream &operator<<(std::ostream &os, const Field &surf) {
  const auto printer = [&](auto &d) {
    if constexpr (std::same_as<std::remove_cvref_t<decltype(d)>,
                               FunctionalData>)
      os << "Functional Data";
    else
      os << d;
  };

  bool first = true;

  const auto keys = surf.keys();
  for (auto &key : surf.keys()) {
    if (not first) os << '\n';
    first = false;

    std::visit(printer, key);
    os << ":\n";
    std::visit(printer, surf[key].data);
  }

  return os;
}

// Allow copy and move set implicitly from all types
Data &Data::operator=(const GriddedField2 &x) {
  data = x;
  return *this;
}
Data &Data::operator=(const Numeric &x) {
  data = x;
  return *this;
}
Data &Data::operator=(const FunctionalData &x) {
  data = x;
  return *this;
}
Data &Data::operator=(GriddedField2 &&x) {
  data = std::move(x);
  return *this;
}
Data &Data::operator=(FunctionalData &&x) {
  data = std::move(x);
  return *this;
}

String Data::data_type() const {
  if (std::holds_alternative<GriddedField2>(data)) return "GriddedField2";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";
  ARTS_ASSERT(false,
              "Cannot be reached, you have added a new type but not "
              "doen the plumbing...")
  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

[[nodiscard]] std::vector<SurfaceKeyVal> Point::keys() const {
  std::vector<SurfaceKeyVal> out;
  out.reserve(enumtyps::SurfaceKeyTypes.size() + type.size());
  for (auto key : enumtyps::SurfaceKeyTypes) out.emplace_back(key);
  for (auto &key : type) out.emplace_back(key.first);
  return out;
}

[[nodiscard]] Index Point::ntype() const {
  return static_cast<Index>(type.size());
}

[[nodiscard]] Index Point::size() const { return nother() + ntype(); }

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

namespace detail {
Numeric numeric_interpolation(
    const GriddedField2 &data,
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
    return data.data(0, 0);
  }

  if (lats.size() == 1) {
    using LatLag = my_interp::Lagrange<0>;
    using LonLag = my_interp::
        Lagrange<1, false, GridType::Cyclic, my_interp::cycle_m180_p180>;

    return interp(data.data, LatLag(0, lat, lats), LonLag(0, lon, lons));
  }

  if (lons.size() == 1) {
    using LatLag = my_interp::Lagrange<1>;
    using LonLag = my_interp::
        Lagrange<0, false, GridType::Cyclic, my_interp::cycle_m180_p180>;

    return interp(data.data, LatLag(0, lat, lats), LonLag(0, lon, lons));
  }

  {
    using LatLag = my_interp::Lagrange<1>;
    using LonLag = my_interp::
        Lagrange<1, false, GridType::Cyclic, my_interp::cycle_m180_p180>;

    return interp(data.data, LatLag(0, lat, lats), LonLag(0, lon, lons));
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
      return detail::numeric_interpolation(x, lat, lon, lat_extrap, lon_extrap);
    };

    return std::visit(call, data.data);
  };
}
}  // namespace detail

Vector2 Field::normal(Numeric lat, Numeric lon, Numeric alt) const try {
  constexpr Vector2 up{180, 0};

  if (not contains(SurfaceKey::h)) {
    return up;
  }

  const auto &z = this->operator[](SurfaceKey::h);

  if (std::holds_alternative<Numeric>(z.data)) {
    return up;
  }

  alt = std::isnan(alt) ? detail::interpolation_function(lat, lon)(z) : alt;

  constexpr Numeric offset = 1e-3;
  const Numeric lat1 = (lat + offset > 90) ? (lat - offset) : (lat + offset),
                lon1 = lon,
                alt1 = detail::interpolation_function(lat1, lon1)(z);
  const Numeric lat2 = lat, lon2 = lon + offset,
                alt2 = detail::interpolation_function(lat2, lon2)(z);

  const Vector3 xyz  = geodetic::to_xyz(ellipsoid, {alt, lat, lon});
  const Vector3 xyz1 = geodetic::to_xyz(ellipsoid, {alt1, lat1, lon1});
  const Vector3 xyz2 = geodetic::to_xyz(ellipsoid, {alt2, lat2, lon2});
  const Vector3 r1{xyz1[0] - xyz[0], xyz1[1] - xyz[1], xyz1[2] - xyz[2]};
  const Vector3 r2{xyz2[0] - xyz[0], xyz2[1] - xyz[1], xyz2[2] - xyz[2]};

  const Vector3 dxyz{r1[1] * r2[2] - r1[2] * r2[1],
                     r1[2] * r2[0] - r1[0] * r2[2],
                     r1[0] * r2[1] - r1[1] * r2[0]};

  return geodetic::from_xyz_dxyz(xyz, dxyz);
} catch (std::exception &e) {
  throw std::runtime_error(
      var_string("Cannot find a normal to the surface at position ",
                 lat,
                 ' ',
                 lon,
                 "\nThe internal error reads: ",
                 e.what()));
}

Point Field::at(Numeric lat, Numeric lon) const {
  Point out;

  const auto interp = detail::interpolation_function(lat, lon);

  for (auto &key : keys()) {
    out[key] = interp(this->operator[](key));
  }

  // Normalize the surface types
  const auto div = 1.0 / [&]() {
    Numeric d = 0.0;
    for (auto &x : out.type) d += x.second;
    return d;
  }();
  for (auto &a : out.type) {
    a.second *= div;
  }

  if (has(SurfaceKey::h)) out.normal = normal(lat, lon, out.elevation);

  return out;
}

Numeric Field::single_value(const KeyVal &key, Numeric lat, Numeric lon) const {
  ARTS_USER_ERROR_IF(
      not std::visit([this](auto &k) { return this->contains(k); }, key),
      "Surface field does not possess the key: ",
      key)

  const auto interp = detail::interpolation_function(lat, lon);

  return interp(this->operator[](key));
}

namespace detail {
constexpr std::pair<Numeric, Numeric> minmax(Numeric x) { return {x, x}; }

std::pair<Numeric, Numeric> minmax(const FunctionalData &) {
  ARTS_USER_ERROR("Cannot extract minmax from functional data");
  return {std::numeric_limits<Numeric>::lowest(),
          std::numeric_limits<Numeric>::max()};
}

std::pair<Numeric, Numeric> minmax(const GriddedField2 &x) {
  return ::minmax(x.data);
}
}  // namespace detail

std::pair<Numeric, Numeric> Field::minmax_single_value(
    const KeyVal &key) const {
  ARTS_USER_ERROR_IF(
      not std::visit([this](auto &k) { return this->contains(k); }, key),
      "Surface field does not possess the key: ",
      key)
  return std::visit([](auto &a) { return detail::minmax(a); },
                    this->operator[](key).data);
}

bool Field::constant_value(const KeyVal &key) const {
  ARTS_USER_ERROR_IF(
      not std::visit([this](auto &k) { return this->contains(k); }, key),
      "Surface field does not possess the key: ",
      key)

  return std::holds_alternative<Numeric>(this->operator[](key).data);
}

[[nodiscard]] ExhaustiveConstVectorView Data::flat_view() const {
  return std::visit(
      [](auto &X) -> ExhaustiveConstVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, GriddedField2>)
          return X.data.flat_view();
        else if constexpr (std::same_as<T, Numeric>)
          return ExhaustiveConstVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ExhaustiveConstVectorView{};
        ARTS_ASSERT(
            false,
            "Cannot be reached, you have added a new type but not done the plumbing...");
      },
      data);
}

[[nodiscard]] ExhaustiveVectorView Data::flat_view() {
  return std::visit(
      [](auto &X) -> ExhaustiveVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, GriddedField2>)
          return X.data.flat_view();
        else if constexpr (std::same_as<T, Numeric>)
          return ExhaustiveVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ExhaustiveVectorView{};
        ARTS_ASSERT(
            false,
            "Cannot be reached, you have added a new type but not done the plumbing...");
      },
      data);
}

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

std::array<std::pair<Index, Numeric>, 4> flat_weights_(const GriddedField2 &v,
                                                       const Numeric &lat,
                                                       const Numeric &lon) {
  using LonLag = my_interp::
      Lagrange<1, false, GridType::Cyclic, my_interp::cycle_m180_p180>;
  using LatLag = my_interp::Lagrange<1>;

  const auto slon = v.shape()[1];
  const bool d1   = v.shape()[0] == 1;
  const bool d2   = slon == 1;

  constexpr auto v1 = std::pair<Index, Numeric>{0, 1.};
  constexpr auto v0 = std::pair<Index, Numeric>{0, 0.};

  if (d1 and d2) {
    return {v1, v0, v0, v0};
  }

  if (d1) {
    const LonLag wlon(0, lon, v.grid<1>());
    return {std::pair<Index, Numeric>{wlon.pos, wlon.lx[0]},
            std::pair<Index, Numeric>{wlon.pos + 1, wlon.lx[1]},
            v0,
            v0};
  }

  if (d2) {
    const LatLag wlat(0, lat, v.grid<0>());
    return {std::pair<Index, Numeric>{wlat.pos, wlat.lx[0]},
            std::pair<Index, Numeric>{wlat.pos + 1, wlat.lx[1]},
            v0,
            v0};
  }

  const LatLag wlat(0, lat, v.grid<0>());
  const LonLag wlon(0, lon, v.grid<1>());
  auto sz = [slon](auto lat_pos, auto lon_pos) {
    return lat_pos * slon + lon_pos;
  };

  return {std::pair<Index, Numeric>{sz(wlat.pos, wlon.pos),
                                    wlat.lx[0] * wlon.lx[0]},
          std::pair<Index, Numeric>{sz(wlat.pos, wlon.pos + 1),
                                    wlat.lx[0] * wlon.lx[1]},
          std::pair<Index, Numeric>{sz(wlat.pos + 1, wlon.pos),
                                    wlat.lx[1] * wlon.lx[0]},
          std::pair<Index, Numeric>{sz(wlat.pos + 1, wlon.pos + 1),
                                    wlat.lx[1] * wlon.lx[1]}};
}

//! Flat weights for the positions in an atmosphere
std::array<std::pair<Index, Numeric>, 4> Data::flat_weights(
    const Numeric &lat, const Numeric &lon) const {
  return std::visit([&](auto &v) { return flat_weights_(v, lat, lon); }, data);
}
}  // namespace Surf

std::string std::formatter<SurfaceKeyVal>::to_string(
    const SurfaceKeyVal &v) const {
  std::string out;
  return std::visit(
      [fmt = tags.get_format_args()](const auto &val) {
        return std::vformat(fmt.c_str(), std::make_format_args(val));
      },
      v);
}

std::ostream &operator<<(std::ostream &os, const SurfaceTypeTag &ppt) {
  return os << ppt.name;
}

std::ostream &operator<<(std::ostream &os, const SurfacePropertyTag &ppt) {
  return os << ppt.name;
}
