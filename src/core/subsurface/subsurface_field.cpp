#include "subsurface.h"

#include "enumsSubsurfaceKey.h"

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

Numeric Point::operator[](const KeyVal &x) const {
  return std::visit([this](auto &key) { return this->operator[](key); }, x);
}

Numeric &Point::operator[](const KeyVal &x) {
  return std::visit(
      [this](auto &key) -> Numeric & {
        return const_cast<Point *>(this)->operator[](key);
      },
      x);
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
  return out;
}

Data::Data()                            = default;
Data::Data(const Data &)                = default;
Data::Data(Data &&) noexcept            = default;
Data &Data::operator=(const Data &)     = default;
Data &Data::operator=(Data &&) noexcept = default;

void Data::adjust_interpolation_extrapolation() {
  if (std::holds_alternative<CartesianSubsurfaceGriddedField3>(data)) {
    auto &field = std::get<CartesianSubsurfaceGriddedField3>(data);

    if (field.grid<0>().size() == 1) {
      alt_upp = InterpolationExtrapolation::Nearest;
      alt_low = InterpolationExtrapolation::Nearest;
    }

    if (field.grid<1>().size() == 1) {
      lat_upp = InterpolationExtrapolation::Nearest;
      lat_low = InterpolationExtrapolation::Nearest;
    }

    if (field.grid<2>().size() == 1) {
      lon_upp = InterpolationExtrapolation::Nearest;
      lon_low = InterpolationExtrapolation::Nearest;
    }
  } else {
    alt_upp = InterpolationExtrapolation::Nearest;
    alt_low = InterpolationExtrapolation::Nearest;
    lat_upp = InterpolationExtrapolation::Nearest;
    lat_low = InterpolationExtrapolation::Nearest;
    lon_upp = InterpolationExtrapolation::Nearest;
    lon_low = InterpolationExtrapolation::Nearest;
  }
}

Data::Data(Numeric x) : data(x) { adjust_interpolation_extrapolation(); }

Data &Data::operator=(Numeric x) {
  data = x;
  adjust_interpolation_extrapolation();
  return *this;
}

Data::Data(CartesianSubsurfaceGriddedField3 x) : data(std::move(x)) {
  adjust_interpolation_extrapolation();
}

Data &Data::operator=(CartesianSubsurfaceGriddedField3 x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

Data::Data(FunctionalData x) : data(std::move(x)) {
  adjust_interpolation_extrapolation();
}

Data &Data::operator=(FunctionalData x) {
  data = std::move(x);
  adjust_interpolation_extrapolation();
  return *this;
}

String Data::data_type() const {
  if (std::holds_alternative<CartesianSubsurfaceGriddedField3>(data))
    return "CartesianSubsurfaceGriddedField3";
  if (std::holds_alternative<Numeric>(data)) return "Numeric";
  if (std::holds_alternative<FunctionalData>(data)) return "FunctionalData";

  ARTS_USER_ERROR("Cannot understand data type; is this a new type")
}

namespace interp {
namespace {
using altlag1 = my_interp::Lagrange<1>;
using altlag0 = my_interp::Lagrange<0>;

using latlag1 = my_interp::Lagrange<1>;
using latlag0 = my_interp::Lagrange<0>;

using lonlag1 =
    my_interp::Lagrange<1, false, GridType::Cyclic, my_interp::cycle_m180_p180>;
using lonlag0 =
    my_interp::Lagrange<0, false, GridType::Cyclic, my_interp::cycle_m180_p180>;

using altlags = std::variant<altlag0, altlag1>;
using latlags = std::variant<latlag0, latlag1>;
using lonlags = std::variant<lonlag0, lonlag1>;

struct Limits {
  Numeric alt_low{std::numeric_limits<Numeric>::lowest()};
  Numeric alt_upp{std::numeric_limits<Numeric>::max()};
  Numeric lat_low{std::numeric_limits<Numeric>::lowest()};
  Numeric lat_upp{std::numeric_limits<Numeric>::max()};
  Numeric lon_low{std::numeric_limits<Numeric>::lowest()};
  Numeric lon_upp{std::numeric_limits<Numeric>::max()};
};

Limits find_limits(const Numeric &) { return {}; }

Limits find_limits(const FunctionalData &) { return {}; }

Limits find_limits(const CartesianSubsurfaceGriddedField3 &gf3) {
  return {.alt_low = gf3.grid<0>().back(),
          .alt_upp = gf3.grid<0>().front(),
          .lat_low = gf3.grid<1>().front(),
          .lat_upp = gf3.grid<1>().back(),
          .lon_low = gf3.grid<2>().front(),
          .lon_upp = gf3.grid<2>().back()};
}

struct ComputeLimit {
  InterpolationExtrapolation type{InterpolationExtrapolation::Linear};
  Numeric alt, lat, lon;
};

constexpr InterpolationExtrapolation combine(InterpolationExtrapolation a,
                                             InterpolationExtrapolation b) {
  using enum InterpolationExtrapolation;
  switch (a) {
    case None: return None;
    case Zero: {
      switch (b) {
        case None:    return None;
        case Zero:    return Zero;
        case Nearest: return Zero;
        case Linear:  return Zero;
      }
      std::unreachable();
    }
    case Nearest: {
      switch (b) {
        case None:    return None;
        case Zero:    return Zero;
        case Nearest: return Nearest;
        case Linear:  return Nearest;
      }
      std::unreachable();
    }
    case Linear: return b;
  }

  std::unreachable();
}

constexpr InterpolationExtrapolation combine(InterpolationExtrapolation a,
                                             InterpolationExtrapolation b,
                                             InterpolationExtrapolation c) {
  return combine(combine(a, b), c);
}

void select(InterpolationExtrapolation lowt,
            InterpolationExtrapolation uppt,
            Numeric lowv,
            Numeric uppv,
            Numeric v,
            Numeric &outv,
            InterpolationExtrapolation &outt) {
  if (v < lowv) {
    outt = lowt;
    if (outt == InterpolationExtrapolation::Nearest) v = lowv;
  } else if (uppv < v) {
    outt = uppt;
    if (outt == InterpolationExtrapolation::Nearest) v = uppv;
  }

  outv = v;
}

ComputeLimit find_limit(const Data &data,
                        const Limits &lim,
                        Numeric alt,
                        Numeric lat,
                        Numeric lon) {
  ComputeLimit out;
  InterpolationExtrapolation a{InterpolationExtrapolation::Linear},
      b{InterpolationExtrapolation::Linear},
      c{InterpolationExtrapolation::Linear};

  select(data.alt_low, data.alt_upp, lim.alt_low, lim.alt_upp, alt, out.alt, a);
  select(data.lat_low, data.lat_upp, lim.lat_low, lim.lat_upp, lat, out.lat, b);
  select(data.lon_low, data.lon_upp, lim.lon_low, lim.lon_upp, lon, out.lon, c);
  out.type = combine(a, b, c);

  return out;
}

Numeric get(const CartesianSubsurfaceGriddedField3 &gf3,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  if (not gf3.ok()) throw std::runtime_error("bad field");

  return std::visit(
      [&data = gf3.data](auto &&al, auto &&la, auto &&lo) {
        return my_interp::interp(data, al, la, lo);
      },
      gf3.grid<0>().size() == 1 ? altlags{gf3.lag<0, altlag0>(alt)}
                                : altlags{gf3.lag<0, altlag1>(alt)},
      gf3.grid<1>().size() == 1 ? latlags{gf3.lag<1, latlag0>(lat)}
                                : latlags{gf3.lag<1, latlag1>(lat)},
      gf3.grid<2>().size() == 1 ? lonlags{gf3.lag<2, lonlag0>(lon)}
                                : lonlags{gf3.lag<2, lonlag1>(lon)});
}

constexpr Numeric get(const Numeric num,
                      const Numeric,
                      const Numeric,
                      const Numeric) {
  return num;
}

Numeric get(const FunctionalData &fd,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  return fd(alt, lat, lon);
}

struct PositionalNumeric {
  const FieldData &data;
  const Numeric alt;
  const Numeric lat;
  const Numeric lon;

  operator Numeric() const {
    return std::visit([&](auto &d) { return get(d, alt, lat, lon); }, data);
  }
};

std::optional<Numeric> get_optional_limit(const Data &data,
                                          const Numeric alt,
                                          const Numeric lat,
                                          const Numeric lon) {
  const auto lim =
      find_limit(data,
                 std::visit([](auto &d) { return find_limits(d); }, data.data),
                 alt,
                 lat,
                 lon);

  ARTS_USER_ERROR_IF(
      lim.type == InterpolationExtrapolation::None,
      "Limit breached.  Position ({}, {}, {}) is out-of-bounds when no extrapolation is wanted",
      lim.alt,
      lim.lat,
      lim.lon)

  if (lim.type == InterpolationExtrapolation::Zero) return 0.0;

  if (lim.type == InterpolationExtrapolation::Nearest)
    return PositionalNumeric{
        .data = data.data, .alt = lim.alt, .lat = lim.lat, .lon = lim.lon};

  return std::nullopt;
}

std::array<std::pair<Index, Numeric>, 8> flat_weight_(
    const CartesianSubsurfaceGriddedField3 &gf3,
    const Numeric alt,
    const Numeric lat,
    const Numeric lon) {
  if (not gf3.ok()) throw std::runtime_error("bad field");

  const Index nalt = gf3.grid<0>().size();
  const Index nlat = gf3.grid<1>().size();
  const Index nlon = gf3.grid<2>().size();

  return std::visit(
      [NN = nlat * nlon, N = nlon]<typename ALT, typename LAT, typename LON>(
          const ALT &al, const LAT &la, const LON &lo) {
        const auto x = interpweights(al, la, lo);
        constexpr std::pair<Index, Numeric> v0{0, 0.0};
        std::array<std::pair<Index, Numeric>, 8> out{
            v0, v0, v0, v0, v0, v0, v0, v0};

        Index m = 0;
        for (Index i = 0; i < al.size(); i++) {
          for (Index j = 0; j < la.size(); j++) {
            for (Index k = 0; k < lo.size(); k++, ++m) {
              out[m] = {{(al.pos + i) * NN + (la.pos + j) * N + lo.pos + k},
                        x[i, j, k]};
            }
          }
        }
        return out;
      },
      nalt == 1 ? altlags{gf3.lag<0, altlag0>(alt)}
                : altlags{gf3.lag<0, altlag1>(alt)},
      nlat == 1 ? latlags{gf3.lag<1, latlag0>(lat)}
                : latlags{gf3.lag<1, latlag1>(lat)},
      nlon == 1 ? lonlags{gf3.lag<2, lonlag0>(lon)}
                : lonlags{gf3.lag<2, lonlag1>(lon)});
}
}  // namespace
}  // namespace interp

Numeric Data::at(const Numeric alt,
                 const Numeric lat,
                 const Numeric lon) const {
  return interp::get_optional_limit(*this, alt, lat, lon)
      .value_or(interp::PositionalNumeric{
          .data = data, .alt = alt, .lat = lat, .lon = lon});
}

Numeric Data::at(const Vector3 pos) const { return at(pos[0], pos[1], pos[2]); }

ConstVectorView Data::flat_view() const {
  return std::visit(
      [](auto &X) -> ConstVectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, CartesianSubsurfaceGriddedField3>)
          return X.data.view_as(X.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return ConstVectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return ConstVectorView{};
        else
          static_assert(
              RawDataType<T>,
              "Cannot be reached, you have added a new type but not done the plumbing...");
      },
      data);
}

VectorView Data::flat_view() {
  return std::visit(
      [](auto &X) -> VectorView {
        using T = std::remove_cvref_t<decltype(X)>;
        if constexpr (std::same_as<T, CartesianSubsurfaceGriddedField3>)
          return X.data.view_as(X.data.size());
        else if constexpr (std::same_as<T, Numeric>)
          return VectorView{X};
        else if constexpr (std::same_as<T, FunctionalData>)
          return VectorView{};
        else
          static_assert(
              RawDataType<T>,
              "Cannot be reached, you have added a new type but not done the plumbing...");
      },
      data);
}

namespace {
std::array<std::pair<Index, Numeric>, 8> flat_weight_(
    const CartesianSubsurfaceGriddedField3 &data,
    Numeric alt,
    Numeric lat,
    Numeric lon) {
  return interp::flat_weight_(data, alt, lat, lon);
}

std::array<std::pair<Index, Numeric>, 8> flat_weight_(Numeric,
                                                      Numeric,
                                                      Numeric,
                                                      Numeric) {
  constexpr std::pair<Index, Numeric> v0{0, 0.0};
  constexpr std::pair<Index, Numeric> v1{0, 1.0};
  return {v1, v0, v0, v0, v0, v0, v0, v0};
}

std::array<std::pair<Index, Numeric>, 8> flat_weight_(const FunctionalData &,
                                                      Numeric,
                                                      Numeric,
                                                      Numeric) {
  constexpr std::pair<Index, Numeric> v0{0, 0.0};
  return {v0, v0, v0, v0, v0, v0, v0, v0};
}
}  // namespace

std::array<std::pair<Index, Numeric>, 8> Data::flat_weight(
    const Numeric alt, const Numeric lat, const Numeric lon) const {
  return std::visit([&](auto &v) { return flat_weight_(v, alt, lat, lon); },
                    data);
}

std::array<std::pair<Index, Numeric>, 8> Data::flat_weight(
    const Vector3 pos) const {
  return flat_weight(pos[0], pos[1], pos[2]);
}

Field::Field()                             = default;
Field::Field(const Field &)                = default;
Field::Field(Field &&) noexcept            = default;
Field &Field::operator=(const Field &)     = default;
Field &Field::operator=(Field &&) noexcept = default;

const std::unordered_map<SubsurfaceKey, Data> &Field::basic() const {
  return map<SubsurfaceKey>();
}
std::unordered_map<SubsurfaceKey, Data> &Field::basic() {
  return map<SubsurfaceKey>();
}

Point Field::at(const Numeric alt, const Numeric lat, const Numeric lon) const
    try {
  ARTS_USER_ERROR_IF(
      alt > bottom_depth,
      "Cannot get values below the deepest point of the subsurface, which is at: {}"
      " m.\nYour depth is: {} m.",
      bottom_depth,
      alt)

  Point out;
  for (auto &&key : keys()) out[key] = operator[](key).at(alt, lat, lon);
  return out;
}
ARTS_METHOD_ERROR_CATCH

Point Field::at(const Vector3 pos) const { return at(pos[0], pos[1], pos[2]); }

Index Field::nbasic() const { return basic().size(); }
}  // namespace Subsurface

std::string std::formatter<SubsurfacePoint>::to_string(
    const SubsurfacePoint &v) const {
  const std::string_view sep = tags.sep(true);

  std::string out = tags.vformat(R"("temperature": )"sv,
                                 v.temperature,
                                 sep,
                                 R"("density": )"sv,
                                 v.density);

  return tags.bracket ? ("{" + out + "}") : out;
}

std::string std::formatter<SubsurfaceField>::to_string(
    const SubsurfaceField &v) const {
  std::string out;

  if (tags.short_str) {
    const std::string_view sep = tags.sep();

    out = tags.vformat(R"("bottom_depth": )"sv,
                       v.bottom_depth,
                       sep,
                       R"("Basic": )"sv,
                       v.basic().size());
  } else {
    const std::string_view sep = tags.sep(true);

    out = tags.vformat(R"("bottom_depth": )"sv,
                       v.bottom_depth,
                       sep,
                       R"("Basic": )"sv,
                       v.basic());
  }

  return tags.bracket ? ("{" + out + "}") : out;
}

std::string std::formatter<SubsurfaceKeyVal>::to_string(
    const SubsurfaceKeyVal &v) const {
  std::string out;
  return std::visit(
      [fmt = tags.get_format_args()](const auto &val) {
        return std::vformat(fmt.c_str(), std::make_format_args(val));
      },
      v);
}
