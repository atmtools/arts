#include "surf.h"
#include "arts_constexpr_math.h"
#include "arts_conversions.h"
#include "debug.h"
#include "gridded_fields.h"
#include "interp.h"
#include <cmath>
#include <exception>
#include <stdexcept>
#include <variant>

namespace geodetic {
Vector3 to_xyz(Vector2 ell, Vector3 pos) {
  using Conversion::cosd;
  using Conversion::sind;
  using Math::pow2;

  const auto [a, b] = ell;
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

  const auto [x, y, z] = xyz;
  const auto [dx, dy, dz] = dxyz;

  const Numeric r = hypot(z, y, x);
  const Vector3 rll = {r, asind(z / r), atan2d(y, x)};
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
} // namespace geodetic

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
    if (not first)
      os << '\n';
    first = false;

    std::visit(printer, key);
    os << ":\n";
    std::visit(printer, surf[key].data);
  }

  return os;
}

String Data::data_type() const {
  if (std::holds_alternative<GriddedField2>(data)) return "GriddedField2";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";
  ARTS_ASSERT(false, "Cannot be reached, you have added a new type but not doen the plumbing...")
  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

  [[nodiscard]] std::vector<KeyVal> Point::keys() const {
    std::vector<KeyVal> out;
    out.reserve(enumtyps::KeyTypes.size() + type.size());
    for (auto key: enumtyps::KeyTypes) out.emplace_back(key);
    for (auto& key: type) out.emplace_back(key.first);
    return out;
  }

  [[nodiscard]] Index Point::ntype() const {return static_cast<Index>(type.size());}

  [[nodiscard]] Index Point::nelem() const {return nother() + ntype();}

Numeric &Point::operator[](const KeyVal &x) {
  return std::visit(
      [this](auto &key) -> Numeric & {
        return const_cast<Point *>(this)->operator[](key);
      },
      x);
}

Numeric Point::operator[](const KeyVal &x) const {
  return std::visit(
      [this](auto &key) -> Numeric { return this->operator[](key); }, x);
}

namespace detail {
Numeric
numeric_interpolation(const GriddedField2 &data, Numeric lat, Numeric lon,
                      std::pair<Extrapolation, Extrapolation> lat_extrap,
                      std::pair<Extrapolation, Extrapolation> lon_extrap) {
  const auto &lats = data.get_numeric_grid(0);
  const auto &lons = data.get_numeric_grid(1);

  if (lat < lats.front()) {
    ARTS_USER_ERROR_IF(lat_extrap.first == Extrapolation::None,
                       "No extrapolation allowed")
    if (lat_extrap.first == Extrapolation::Zero)
      return 0.0;
    if (lat_extrap.first == Extrapolation::Nearest)
      lat = lats.front();
  }

  if (lat > lats.back()) {
    ARTS_USER_ERROR_IF(lat_extrap.second == Extrapolation::None,
                       "No extrapolation allowed")
    if (lat_extrap.second == Extrapolation::Zero)
      return 0.0;
    if (lat_extrap.second == Extrapolation::Nearest)
      lat = lats.back();
  }

  if (lat < lons.front()) {
    ARTS_USER_ERROR_IF(lon_extrap.first == Extrapolation::None,
                       "No extrapolation allowed")
    if (lon_extrap.first == Extrapolation::Zero)
      return 0.0;
    if (lon_extrap.first == Extrapolation::Nearest)
      lat = lons.front();
  }

  if (lat > lons.back()) {
    ARTS_USER_ERROR_IF(lon_extrap.second == Extrapolation::None,
                       "No extrapolation allowed")
    if (lon_extrap.second == Extrapolation::Zero)
      return 0.0;
    if (lon_extrap.second == Extrapolation::Nearest)
      lat = lons.back();
  }

  if (lats.size() == 1 and lons.size() == 1) {
    return data.data(0, 0);
  }

  if (lats.size() == 1) {
    using LatLag = my_interp::Lagrange<0>;
    using LonLag = my_interp::Lagrange<1, false, my_interp::GridType::Cyclic,
                                       my_interp::cycle_m180_p180>;

    return interp(data.data, LatLag(0, lat, lats), LonLag(0, lon, lons));
  }

  if (lons.size() == 1) {
    using LatLag = my_interp::Lagrange<1>;
    using LonLag = my_interp::Lagrange<0, false, my_interp::GridType::Cyclic,
                                       my_interp::cycle_m180_p180>;

    return interp(data.data, LatLag(0, lat, lats), LonLag(0, lon, lons));
  }

  {
    using LatLag = my_interp::Lagrange<1>;
    using LonLag = my_interp::Lagrange<1, false, my_interp::GridType::Cyclic,
                                       my_interp::cycle_m180_p180>;

    return interp(data.data, LatLag(0, lat, lats), LonLag(0, lon, lons));
  }
}

Numeric numeric_interpolation(const FunctionalData &f, Numeric lat, Numeric lon,
                              std::pair<Extrapolation, Extrapolation>,
                              std::pair<Extrapolation, Extrapolation>) {
  return f(lat, lon);
}

constexpr Numeric
numeric_interpolation(Numeric x, Numeric, Numeric,
                      std::pair<Extrapolation, Extrapolation>,
                      std::pair<Extrapolation, Extrapolation>) {
  return x;
}

auto interpolation_function(Numeric lat, Numeric lon) {
  return [lat, lon](auto &data) {
    const auto call = [lat, lon,
                       lat_extrap = std::pair{data.lat_low, data.lat_upp},
                       lon_extrap =
                           std::pair{data.lon_low, data.lon_upp}](auto &&x) {
      return detail::numeric_interpolation(x, lat, lon, lat_extrap, lon_extrap);
    };

    return std::visit(call, data.data);
  };
}
} // namespace detail

Vector2 Field::normal(Vector2 ellipsoid, Numeric lat, Numeric lon,
                      Numeric alt) const try {
  constexpr Vector2 up{180, 0};

  if (not contains(Key::h)) {
    return up;
  }

  const auto &z = this->operator[](Key::h);

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

  const Vector3 xyz = geodetic::to_xyz(ellipsoid, {alt, lat, lon});
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
      var_string("Cannot find a normal to the surface at position ", lat, ' ',
                 lon, "\nThe internal error reads: ", e.what()));
}

Point Field::at(Numeric lat, Numeric lon, Vector2 ellipsoid) const {
  Point out;

  const auto interp = detail::interpolation_function(lat, lon);

  for (auto &key : keys()) {
    out[key] = interp(this->operator[](key));
  }

  // Normalize the surface types
  const auto div =
      1.0 / std::reduce(out.type.begin(), out.type.end(), Numeric{0.0},
                        [](auto &a, auto &x) { return a + x.second; });
  for (auto &a : out.type) {
    a.second *= div;
  }

  out.normal = normal(ellipsoid, lat, lon, out.elevation);

  return out;
}

Numeric Field::single_value(const KeyVal& key, Numeric lat, Numeric lon) const {
  ARTS_USER_ERROR_IF (not std::visit([this](auto& k){return this->contains(k);}, key), "Surface field has no elevation")

  const auto interp = detail::interpolation_function(lat, lon);
  
  return interp(this->operator[](key));
}
} // namespace Surf

std::ostream &operator<<(std::ostream &os, const SurfaceTypeTag &ppt) {
  return os << ppt.name;
}

std::ostream &operator<<(std::ostream &os, const SurfacePropertyTag &ppt) {
  return os << ppt.name;
}
