#include "subsurf_field.h"

#include <debug.h>
#include <enumsSubsurfaceKey.h>
#include <functional_field_interp.h>

#include <utility>

namespace Subsurface {
Point::Point()                             = default;
Point::Point(const Point &)                = default;
Point::Point(Point &&) noexcept            = default;
Point &Point::operator=(const Point &)     = default;
Point &Point::operator=(Point &&) noexcept = default;

Numeric Point::operator[](SubsurfaceKey x) const {
  switch (x) {
    case SubsurfaceKey::t:   return temperature;
    case SubsurfaceKey::rho: return density;
  }
  std::unreachable();
}

Numeric &Point::operator[](SubsurfaceKey x) {
  switch (x) {
    case SubsurfaceKey::t:   return temperature;
    case SubsurfaceKey::rho: return density;
  }
  std::unreachable();
}

Numeric Point::operator[](const SubsurfacePropertyTag &x) const { return props.at(x); }

Numeric &Point::operator[](const SubsurfacePropertyTag &x) { return props[x]; }

Numeric Point::operator[](const KeyVal &x) const {
  return std::visit([this](auto &key) { return this->operator[](key); }, x);
}

Numeric &Point::operator[](const KeyVal &x) {
  return std::visit([this](auto &key) -> Numeric & { return const_cast<Point *>(this)->operator[](key); }, x);
}

bool Point::contains(const KeyVal &x) const {
  return std::visit([this](auto &key) { return this->has(key); }, x);
}

Index Point::size() const {
  return nbasic();  // +...
}

Index Point::nbasic() const { return enumsize::SubsurfaceKeySize; }

std::vector<KeyVal> Point::keys() const {
  std::vector<KeyVal> out;
  out.reserve(size());
  for (auto &a : enumtyps::SubsurfaceKeyTypes) out.emplace_back(a);
  for (auto &a : props | stdv::keys) out.emplace_back(a);
  return out;
}

Data::Data()                            = default;
Data::Data(const Data &)                = default;
Data::Data(Data &&) noexcept            = default;
Data &Data::operator=(const Data &)     = default;
Data &Data::operator=(Data &&) noexcept = default;

void Data::adjust_interpolation_extrapolation() {
  using enum InterpolationExtrapolation;
  const auto [alt, lat, lon] = field_interp::single_point_dims(data);
  if (alt) alt_low = alt_upp = Nearest;
  if (lat) lat_low = lat_upp = Nearest;
  if (lon) lon_low = lon_upp = Nearest;
}

Data::Data(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

Data::Data(GeodeticField3 x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(GeodeticField3 x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

Data::Data(FunctionalData x) : data(std::move(x)) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(FunctionalData x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

String Data::data_type() const {
  if (std::holds_alternative<GeodeticField3>(data)) return "GeodeticField3";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";

  arts::user_error("Cannot understand data type; is this a new type");
}

namespace {
std::array<field_interp::extrapolation, 3> extrapolation(const Data &d) {
  return {{{d.alt_low, d.alt_upp}, {d.lat_low, d.lat_upp}, {d.lon_low, d.lon_upp}}};
}
}  // namespace

Numeric Data::at(const Numeric alt, const Numeric lat, const Numeric lon) const {
  return field_interp::at(data, extrapolation(*this), {alt, lat, lon});
}

Numeric Data::at(const Vector3 pos) const { return at(pos[0], pos[1], pos[2]); }

ConstVectorView Data::flat_view() const { return field_interp::flat_view(data); }

VectorView Data::flat_view() { return field_interp::flat_view(data); }

std::vector<std::pair<Index, Numeric>> Data::flat_weight(const Numeric alt,
                                                         const Numeric lat,
                                                         const Numeric lon) const {
  return field_interp::flat_weight(data, {alt, lat, lon});
}

std::vector<std::pair<Index, Numeric>> Data::flat_weight(const Vector3 pos) const {
  return flat_weight(pos[0], pos[1], pos[2]);
}

bool Data::ok() const { return field_interp::ok(data); }

Field::Field()                             = default;
Field::Field(const Field &)                = default;
Field::Field(Field &&) noexcept            = default;
Field &Field::operator=(const Field &)     = default;
Field &Field::operator=(Field &&) noexcept = default;

Index Field::nbasic() const { return other.size(); }
Index Field::size() const { return nbasic(); }

[[nodiscard]] std::vector<KeyVal> Field::keys() const {
  std::vector<KeyVal> out;
  out.reserve(size());
  for (const auto &[key, _] : other) out.emplace_back(key);
  for (const auto &[key, _] : props) out.emplace_back(key);
  return out;
}

Data &Field::operator[](const SubsurfaceKey &key) { return other[key]; }

Data &Field::operator[](const SubsurfacePropertyTag &key) { return props[key]; }

Data &Field::operator[](const KeyVal &key) {
  return std::visit([this](auto &k) -> Data & { return this->operator[](k); }, key);
}

const Data &Field::operator[](const SubsurfaceKey &key) const { return other.at(key); }

const Data &Field::operator[](const SubsurfacePropertyTag &key) const { return props.at(key); }

const Data &Field::operator[](const KeyVal &key) const {
  return std::visit([this](auto &k) -> const Data & { return this->operator[](k); }, key);
}

bool Field::contains(const SubsurfaceKey &key) const { return other.contains(key); }

bool Field::contains(const KeyVal &key) const {
  return std::visit([this](auto &k) { return this->contains(k); }, key);
}

Point Field::at(const Numeric alt, const Numeric lat, const Numeric lon) const try {
  if (alt < bottom_depth)
    arts::user_error(
        "Cannot get values below the deepest point of the subsurface, "
        "which is at: {} m.\nYour depth is: {} m.",
        bottom_depth,
        alt);

  Point out;

  static_assert(
      sizeof(SubsurfaceField) == sizeof(std::unordered_map<SubsurfaceKey, Data>) * 2 + sizeof(Numeric),
      "The loops below must be over all keys, a size change of SubsurfaceField indicates that the number of keys have changed");
  out.props.reserve(props.size());

  for (auto &[key, data] : other) out[key] = data.at(alt, lat, lon);
  for (auto &[key, data] : props) out[key] = data.at(alt, lat, lon);

  return out;
}
ARTS_METHOD_ERROR_CATCH

Point Field::at(const Vector3 pos) const { return at(pos[0], pos[1], pos[2]); }
}  // namespace Subsurface

std::string std::formatter<SubsurfacePoint>::to_string(const SubsurfacePoint &v) const {
  const std::string_view sep = tags.sep();

  std::string out = tags.vformat(
      R"("temperature": )"sv, v.temperature, sep, R"("density": )"sv, v.density, "\nproperties:\n"sv, v.props);

  return tags.bracket ? ("{" + out + "}") : out;
}

std::string std::formatter<SubsurfaceField>::to_string(const SubsurfaceField &v) const {
  std::string out;

  if (tags.short_str) {
    const std::string_view sep = " "sv;

    out = tags.vformat(R"("bottom_depth": )"sv,
                       v.bottom_depth,
                       sep,
                       R"("Other": )"sv,
                       v.other.size(),
                       sep,
                       R"("Prop": )"sv,
                       v.props.size());
  } else {
    const std::string_view sep = tags.sep();

    out = tags.vformat(R"("bottom_depth": )"sv, v.bottom_depth);

    if (not v.other.empty()) { out += tags.vformat(sep, R"("Other": )"sv, v.other); }

    if (not v.props.empty()) { out += tags.vformat(sep, R"("Prop": )"sv, v.props); }
  }

  return tags.bracket ? ("{" + out + "}") : out;
}
