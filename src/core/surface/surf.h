#pragma once

#include <enumsInterpolationExtrapolation.h>
#include <enumsSurfaceKey.h>
#include <fieldmap.h>
#include <matpack.h>
#include <mystring.h>
#include <operators.h>

#include <limits>
#include <type_traits>
#include <unordered_map>
#include <variant>

struct SurfacePropertyTag {
  String name;

  auto operator<=>(const SurfacePropertyTag &x) const = default;
};

namespace std {
template <>
struct hash<SurfacePropertyTag> {
  std::size_t operator()(const SurfacePropertyTag &pp) const noexcept {
    return std::hash<String>{}(pp.name);
  }
};
}  // namespace std

using SurfaceKeyVal = std::variant<SurfaceKey, SurfacePropertyTag>;

namespace Surf {
template <typename T>
concept isSurfaceKey = std::same_as<std::remove_cvref_t<T>, SurfaceKey>;

template <typename T>
concept isSurfacePropertyTag =
    std::same_as<std::remove_cvref_t<T>, SurfacePropertyTag>;

template <typename T>
concept KeyType = isSurfaceKey<T> or isSurfacePropertyTag<T>;

struct Point {
  Numeric elevation{std::numeric_limits<Numeric>::quiet_NaN()};
  Numeric temperature{std::numeric_limits<Numeric>::quiet_NaN()};
  Vector2 normal{std::numeric_limits<Numeric>::quiet_NaN(),
                 std::numeric_limits<Numeric>::quiet_NaN()};
  std::unordered_map<SurfacePropertyTag, Numeric> prop;

  Point();
  Point(const Point &);
  Point(Point &&) noexcept;
  Point &operator=(const Point &);
  Point &operator=(Point &&) noexcept;

  Numeric &operator[](SurfaceKey x);
  Numeric &operator[](const SurfacePropertyTag &x);
  Numeric &operator[](const SurfaceKeyVal &x);

  Numeric operator[](SurfaceKey x) const;
  Numeric operator[](const SurfacePropertyTag &x) const;
  Numeric operator[](const SurfaceKeyVal &x) const;

  [[nodiscard]] std::vector<SurfaceKeyVal> keys() const;

  [[nodiscard]] Index size() const;
  [[nodiscard]] static constexpr Index nother() {
    return static_cast<Index>(enumtyps::SurfaceKeyTypes.size());
  }

  template <KeyType T>
  constexpr bool contains(T &&k) const {
    if constexpr (isSurfacePropertyTag<T>)
      return prop.contains(std::forward<T>(k));
    else if constexpr (isSurfaceKey<T>)
      return true;
    else
      return false;
  }

  template <KeyType... Ts>
  constexpr bool has(Ts &&...keys) const {
    return (contains<Ts>(std::forward<Ts>(keys)) and ...);
  }
};

using FunctionalData = NumericBinaryOperator;
using FieldData      = std::variant<GriddedField2, Numeric, FunctionalData>;

struct FunctionalDataAlwaysThrow {
  Numeric operator()(Numeric, Numeric) const {
    ARTS_USER_ERROR("You touched the field but did not set any data")
  }
};

//! Hold all atmospheric data
struct Data {
  FieldData data{FunctionalData{FunctionalDataAlwaysThrow{}}};
  InterpolationExtrapolation lat_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lat_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_low{InterpolationExtrapolation::None};

  // Standard
  Data();
  Data(const Data &);
  Data(Data &&) noexcept;
  Data &operator=(const Data &);
  Data &operator=(Data &&) noexcept;

  void adjust_interpolation_extrapolation();

  // Allow copy implicitly from all types
  explicit Data(Numeric x);
  Data &operator=(Numeric x);

  explicit Data(GriddedField2 x);
  Data &operator=(GriddedField2 x);

  explicit Data(FunctionalData x);
  Data &operator=(FunctionalData x);

  [[nodiscard]] String data_type() const;

  template <typename T>
  [[nodiscard]] T get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  template <typename T>
  [[nodiscard]] T get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  void rescale(Numeric);

  [[nodiscard]] ConstVectorView flat_view() const;

  [[nodiscard]] VectorView flat_view();

  //! Flat weights for the positions on the surface
  [[nodiscard]] std::array<std::pair<Index, Numeric>, 4> flat_weights(
      const Numeric &lat, const Numeric &lon) const;

  [[nodiscard]] Numeric at(const Numeric lat, const Numeric lon) const;
};

struct Field final : FieldMap::Map<Data, SurfaceKey, SurfacePropertyTag> {
  //! The ellipsoid used for the surface, in [a, b] in meters
  Vector2 ellipsoid;

  Field();
  Field(const Field &);
  Field(Field &&) noexcept;
  Field &operator=(const Field &);
  Field &operator=(Field &&) noexcept;

  /** Compute the values at a single point
   *
   * Note that this method uses the pos's alt to compute the normal only,
   * for field values, only the lat and lon of the pos is used.  It is
   * thus allowed that the pos is above or below the surface, though it
   * is adviced that this is checked at calling site
   *
   * @param lat The latitude
   * @param lon The longitude
   * @return Point values at the surface
   */
  [[nodiscard]] Point at(Numeric lat, Numeric lon) const;

  /** Get the normal of the surface at a given position
   *
   * @param lat The latitude
   * @param lon The longitude
   * @param alt The altitude of the point.  If NaN, it is computed first
   * @return Vector2 [za, aa]
   */
  [[nodiscard]] Vector2 normal(
      Numeric lat,
      Numeric lon,
      Numeric alt = std::numeric_limits<Numeric>::quiet_NaN()) const;

  [[nodiscard]] Numeric single_value(const KeyVal &key,
                                     Numeric lat,
                                     Numeric lon) const;

  [[nodiscard]] std::pair<Numeric, Numeric> minmax_single_value(
      const KeyVal &key) const;

  [[nodiscard]] bool constant_value(const KeyVal &key) const;
};

static_assert(
    std::same_as<typename Field::KeyVal, SurfaceKeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");
}  // namespace Surf

using SurfaceData         = Surf::Data;
using SurfacePoint        = Surf::Point;
using SurfaceField        = Surf::Field;
using ArrayOfSurfacePoint = Array<SurfacePoint>;

bool operator==(const SurfaceKeyVal &, SurfaceKey);
bool operator==(SurfaceKey, const SurfaceKeyVal &);
bool operator==(const SurfaceKeyVal &, const SurfacePropertyTag &);
bool operator==(const SurfacePropertyTag &, const SurfaceKeyVal &);
template <>
struct std::formatter<SurfacePropertyTag> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SurfacePropertyTag &v,
                              FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(ctx.out(), "{}{}{}", quote, v.name, quote);
  }
};

template <>
struct std::formatter<SurfaceKeyVal> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  [[nodiscard]] std::string to_string(const SurfaceKeyVal &v) const;

  template <class FmtContext>
  FmtContext::iterator format(const SurfaceKeyVal &v, FmtContext &ctx) const {
    tags.format(ctx, to_string(v));
    return ctx.out();
  }
};

template <>
struct std::formatter<Surf::FunctionalData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Surf::FunctionalData &,
                              FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(ctx.out(), "{}functional-data{}", quote, quote);
  }
};
template <>
struct std::formatter<Surf::Data> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Surf::Data &v, FmtContext &ctx) const {
    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.data, sep);
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.lat_upp, sep, v.lat_low, sep, v.lon_upp, sep, v.lon_low);
    tags.add_if_bracket(ctx, ']');
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<SurfaceField> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SurfaceField &v, FmtContext &ctx) const {
    tags.add_if_bracket(ctx, '{');
    std::format_to(ctx.out(), R"()");
    tags.format(ctx, R"("Ellipsoid": )"sv, v.ellipsoid);

    if (tags.short_str) {
      std::format_to(ctx.out(),
                     R"(, "SurfaceKey": {}, "SurfacePropertyTag": {})",
                     v.map<SurfaceKey>().size(),
                     v.map<SurfacePropertyTag>().size());
    } else {
      const std::string_view sep = tags.sep(true);
      tags.format(ctx,
                  sep,
                  R"("SurfaceKey": )"sv,
                  v.map<SurfaceKey>(),
                  sep,
                  R"("SurfacePropertyTag": )"sv,
                  v.map<SurfacePropertyTag>());
    }

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};

template <>
struct std::formatter<SurfacePoint> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SurfacePoint &v, FmtContext &ctx) const {
    tags.add_if_bracket(ctx, '{');
    std::format_to(ctx.out(),
                   R"("elevation": {}, "temperature": {}, "normal": {:B,})",
                   v.elevation,
                   v.temperature,
                   v.normal);

    if (tags.short_str) {
      std::format_to(ctx.out(), R"(, "SurfacePropertyTag": {})", v.prop.size());
    } else {
      const std::string_view sep = tags.sep(true);
      tags.format(ctx, sep, R"("SurfacePropertyTag": )"sv, v.prop);
    }

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};
