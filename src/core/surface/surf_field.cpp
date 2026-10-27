#include "surf_field.h"

#include <arts_constexpr_math.h>
#include <arts_conversions.h>
#include <debug.h>
#include <enumsSurfaceKey.h>
#include <functional_field_interp.h>

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

Numeric &Surf::Point::operator[](const SurfacePropertyTag &x) { return prop[x]; }

Numeric Surf::Point::operator[](SurfaceKey x) const {
  switch (x) {
    case SurfaceKey::h: return elevation;
    case SurfaceKey::t: return temperature;
  }
  std::unreachable();
}

Numeric Surf::Point::operator[](const SurfacePropertyTag &x) const { return prop.at(x); }

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

  const Numeric dr   = (clat * clon * dx + slat * dz + clat * slon * dy) / norm;
  const Numeric dlat = (-slat * clon * dx + clat * dz - slat * slon * dy) / (norm * r);
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

  arts::user_error("Cannot understand data type; is this a new type");
}

[[nodiscard]] std::vector<SurfaceKeyVal> Point::keys() const {
  std::vector<SurfaceKeyVal> out;
  out.reserve(enumtyps::SurfaceKeyTypes.size());
  for (auto key : enumtyps::SurfaceKeyTypes) out.emplace_back(key);
  return out;
}

[[nodiscard]] Index Point::size() const { return nother(); }

Numeric &Point::operator[](const SurfaceKeyVal &x) {
  return std::visit([this](auto &key) -> Numeric & { return const_cast<Point *>(this)->operator[](key); }, x);
}

Numeric Point::operator[](const SurfaceKeyVal &x) const {
  return std::visit([this](auto &key) -> Numeric { return this->operator[](key); }, x);
}

namespace {
std::array<field_interp::extrapolation, 2> extrapolation(const Data &d) {
  return {{{d.lat_low, d.lat_upp}, {d.lon_low, d.lon_upp}}};
}
}  // namespace

Data &Field::operator[](const SurfaceKey &key) { return other[key]; }

Data &Field::operator[](const SurfacePropertyTag &key) { return props[key]; }

Data &Field::operator[](const KeyVal &key) {
  return std::visit([this](auto &k) -> Data & { return this->operator[](k); }, key);
}

const Data &Field::operator[](const SurfaceKey &key) const { return other.at(key); }

const Data &Field::operator[](const SurfacePropertyTag &key) const { return props.at(key); }

const Data &Field::operator[](const KeyVal &key) const {
  return std::visit([this](auto &k) -> const Data & { return this->operator[](k); }, key);
}

bool Field::contains(const SurfaceKey &key) const { return other.contains(key); }

bool Field::contains(const SurfacePropertyTag &key) const { return props.contains(key); }

bool Field::contains(const KeyVal &key) const {
  return std::visit([this](auto &k) -> bool { return this->contains(k); }, key);
}

Vector2 Field::normal(Numeric lat, Numeric lon, Numeric alt) const try {
  if (bad_ellipsoid()) arts::user_error("Ellipsoid must have positive axes: {:B,}", ellipsoid);

  constexpr Vector2 up{180, 0};

  if (not contains(SurfaceKey::h)) { return up; }

  const auto &z = other.at(SurfaceKey::h);

  if (std::holds_alternative<Numeric>(z.data)) { return up; }

  alt = std::isnan(alt) ? z.at(lat, lon) : alt;

  constexpr Numeric offset = 1e-3;
  const Numeric     lat1 = (lat + offset > 90) ? (lat - offset) : (lat + offset), lon1 = lon, alt1 = z.at(lat1, lon1);
  const Numeric     lat2 = lat, lon2 = lon + offset, alt2 = z.at(lat2, lon2);

  const Vector3 xyz  = to_xyz(ellipsoid, {alt, lat, lon});
  const Vector3 xyz1 = to_xyz(ellipsoid, {alt1, lat1, lon1});
  const Vector3 xyz2 = to_xyz(ellipsoid, {alt2, lat2, lon2});
  const Vector3 r1{xyz1[0] - xyz[0], xyz1[1] - xyz[1], xyz1[2] - xyz[2]};
  const Vector3 r2{xyz2[0] - xyz[0], xyz2[1] - xyz[1], xyz2[2] - xyz[2]};

  const Vector3 dxyz{r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2], r1[0] * r2[1] - r1[1] * r2[0]};

  return from_xyz_dxyz(xyz, dxyz);
} catch (std::exception &e) {
  throw std::runtime_error(
      std::format("Cannot find a normal to the surface at position {} {}\n{}", lat, lon, indent(e.what())));
}

Numeric Data::at(const Numeric lat, const Numeric lon) const {
  return field_interp::at(data, extrapolation(*this), {lat, lon});
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

  static_assert(
      sizeof(SurfaceField) == sizeof(std::unordered_map<SurfaceKey, Data>) * 2 + 2 * sizeof(Numeric),
      "The loops below must be over all keys, a size change of SurfaceField indicates that the number of keys have changed");
  out.prop.reserve(props.size());

  for (auto &[key, data] : other) out[key] = data.at(lat, lon);
  for (auto &[key, data] : props) out[key] = data.at(lat, lon);

  out.normal = normal(lat, lon);

  return out;
}

Numeric Field::single_value(const KeyVal &key, Numeric lat, Numeric lon) const {
  if (not contains(key)) arts::user_error("Surface field does not possess the key: {}", key);

  return this->operator[](key).at(lat, lon);
}

namespace {
constexpr std::pair<Numeric, Numeric> minmax(Numeric x) { return {x, x}; }

std::pair<Numeric, Numeric> minmax(const FunctionalData &) {
  arts::user_error("Cannot extract minmax from functional data");
  return {std::numeric_limits<Numeric>::lowest(), std::numeric_limits<Numeric>::max()};
}

std::pair<Numeric, Numeric> minmax(const GeodeticField2 &x) { return matpack::minmax(x.data); }
}  // namespace

std::pair<Numeric, Numeric> Field::minmax_single_value(const KeyVal &key) const {
  if (not contains(key)) arts::user_error("Surface field does not possess the key: {}", key);
  return std::visit([](auto &a) { return minmax(a); }, this->operator[](key).data);
}

bool Field::constant_value(const KeyVal &key) const {
  if (not contains(key)) arts::user_error("Surface field does not possess the key: {}", key);

  return std::holds_alternative<Numeric>(this->operator[](key).data);
}

ConstVectorView Data::flat_view() const { return field_interp::flat_view(data); }

VectorView Data::flat_view() { return field_interp::flat_view(data); }

//! Flat weights for the positions on the surface
std::vector<std::pair<Index, Numeric>> Data::flat_weights(const Numeric &lat, const Numeric &lon) const {
  return field_interp::flat_weight(data, {lat, lon});
}

bool Data::ok() const { return field_interp::ok(data); }

void Data::adjust_interpolation_extrapolation() {
  using enum InterpolationExtrapolation;
  const auto [lat, lon] = field_interp::single_point_dims(data);
  if (lat) lat_low = lat_upp = Nearest;
  if (lon) lon_low = lon_upp = Nearest;
}

Data::Data(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

Data::Data(GeodeticField2 x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

Data::Data(FunctionalData x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

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

[[nodiscard]] bool Field::bad_ellipsoid() const { return not(ellipsoid[1] > 0 and ellipsoid[0] >= ellipsoid[1]); }

bool Point::contains(const SurfaceKeyVal &k) const {
  return std::visit([this](auto &key) -> bool { return this->contains(key); }, k);
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

bool operator==(const SurfacePropertyTag &lhs, const SurfaceKeyVal &rhs) { return rhs == lhs; }

static_assert(std::same_as<typename SurfaceField::KeyVal, SurfaceKeyVal>,
              "The order of arguments in the template of which Field inherits from is "
              "wrong.  KeyVal must be defined in the same way for this to work.");
