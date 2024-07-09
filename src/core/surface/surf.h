#pragma once

#include <matpack.h>

#include <limits>
#include <ostream>
#include <type_traits>
#include <unordered_map>
#include <variant>

#include "enums.h"
#include "fieldmap.h"
#include "mystring.h"

struct SurfaceTypeTag {
  String name;

  auto operator<=>(const SurfaceTypeTag &x) const = default;

  friend std::ostream &operator<<(std::ostream &, const SurfaceTypeTag &);
};

struct SurfacePropertyTag {
  String name;

  auto operator<=>(const SurfacePropertyTag &x) const = default;

  friend std::ostream &operator<<(std::ostream &, const SurfacePropertyTag &);
};

namespace std {
template <>
struct hash<SurfaceTypeTag> {
  std::size_t operator()(const SurfaceTypeTag &pp) const noexcept {
    return std::hash<String>{}(pp.name);
  }
};

template <>
struct hash<SurfacePropertyTag> {
  std::size_t operator()(const SurfacePropertyTag &pp) const noexcept {
    return std::hash<String>{}(pp.name);
  }
};
}  // namespace std

using SurfaceKeyVal =
    std::variant<SurfaceKey, SurfaceTypeTag, SurfacePropertyTag>;

std::ostream &operator<<(std::ostream &os, const SurfaceKeyVal &key);

namespace Surf {
template <typename T>
concept isSurfaceKey = std::same_as<std::remove_cvref_t<T>, SurfaceKey>;

template <typename T>
concept isSurfaceTypeTag = std::same_as<std::remove_cvref_t<T>, SurfaceTypeTag>;

template <typename T>
concept isSurfacePropertyTag =
    std::same_as<std::remove_cvref_t<T>, SurfacePropertyTag>;

template <typename T>
concept KeyType =
    isSurfaceKey<T> or isSurfaceTypeTag<T> or isSurfacePropertyTag<T>;

struct Point {
  Numeric elevation{std::numeric_limits<Numeric>::quiet_NaN()};
  Numeric temperature{std::numeric_limits<Numeric>::quiet_NaN()};
  Vector3 wind{std::numeric_limits<Numeric>::quiet_NaN(),
               std::numeric_limits<Numeric>::quiet_NaN(),
               std::numeric_limits<Numeric>::quiet_NaN()};
  Vector2 normal{std::numeric_limits<Numeric>::quiet_NaN(),
                 std::numeric_limits<Numeric>::quiet_NaN()};
  std::unordered_map<SurfaceTypeTag, Numeric> type;
  std::unordered_map<SurfacePropertyTag, Numeric> prop;

  template <KeyType Key>
  Numeric &operator[](Key &&x) {
    if constexpr (isSurfaceKey<Key>) {
      switch (std::forward<Key>(x)) {
        case SurfaceKey::h:
          return elevation;
        case SurfaceKey::t:
          return temperature;
        case SurfaceKey::wind_u:
          return wind[0];
        case SurfaceKey::wind_v:
          return wind[1];
        case SurfaceKey::wind_w:
          return wind[2];
      }
    } else if constexpr (isSurfaceTypeTag<Key>) {
      return type[std::forward<Key>(x)];
    } else if constexpr (isSurfacePropertyTag<Key>) {
      return prop[std::forward<Key>(x)];
    }

    return temperature;
  }

  template <KeyType Key>
  Numeric operator[](Key &&x) const {
    if constexpr (isSurfaceKey<Key>) {
      switch (std::forward<Key>(x)) {
        case SurfaceKey::h:
          return elevation;
        case SurfaceKey::t:
          return temperature;
        case SurfaceKey::wind_u:
          return wind[0];
        case SurfaceKey::wind_v:
          return wind[1];
        case SurfaceKey::wind_w:
          return wind[2];
      }
    } else if constexpr (isSurfaceTypeTag<Key>) {
      return type.at(std::forward<Key>(x));
    } else if constexpr (isSurfacePropertyTag<Key>) {
      return prop.at(std::forward<Key>(x));
    }

    return 0.0;
  }

  Numeric &operator[](const SurfaceKeyVal &x);

  Numeric operator[](const SurfaceKeyVal &x) const;

  [[nodiscard]] std::vector<SurfaceKeyVal> keys() const;

  [[nodiscard]] Index size() const;
  [[nodiscard]] Index ntype() const;
  [[nodiscard]] static constexpr Index nother() {
    return static_cast<Index>(enumtyps::SurfaceKeyTypes.size());
  }

  template <KeyType T>
  constexpr bool contains(T &&k) const {
    if constexpr (isSurfaceTypeTag<T>) return type.contains(std::forward<T>(k));
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

  friend std::ostream &operator<<(std::ostream &, const Point &);
};

using FunctionalData = std::function<Numeric(Numeric, Numeric)>;
using FieldData      = std::variant<GriddedField2, Numeric, FunctionalData>;

struct FunctionalDataAlwaysThrow {
  std::string error{"Undefined data"};
  Numeric operator()(Numeric, Numeric) const { ARTS_USER_ERROR(error) }
};

//! Hold all atmospheric data
struct Data {
  FieldData data{FunctionalData{FunctionalDataAlwaysThrow{
      "You touched the field but did not set any data"}}};
  InterpolationExtrapolation lat_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lat_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_low{InterpolationExtrapolation::None};

  // Standard
  Data()                        = default;
  Data(const Data &)            = default;
  Data(Data &&)                 = default;
  Data &operator=(const Data &) = default;
  Data &operator=(Data &&)      = default;

  // Allow copy and move construction implicitly from all types
  explicit Data(const GriddedField2 &x) : data(x) {}
  explicit Data(const Numeric &x) : data(x) {}
  explicit Data(const FunctionalData &x) : data(x) {}
  explicit Data(GriddedField2 &&x) : data(std::move(x)) {}
  explicit Data(FunctionalData &&x) : data(std::move(x)) {}

  [[nodiscard]] String data_type() const;

  // Allow copy and move set implicitly from all types
  Data &operator=(const GriddedField2 &x);
  Data &operator=(const Numeric &x);
  Data &operator=(const FunctionalData &x);
  Data &operator=(GriddedField2 &&x);
  Data &operator=(FunctionalData &&x);

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

  [[nodiscard]] ExhaustiveConstVectorView flat_view() const;

  [[nodiscard]] ExhaustiveVectorView flat_view();

  //! Flat weights for the positions on the surface
  [[nodiscard]] std::array<std::pair<Index, Numeric>, 4> flat_weights(
      const Numeric &lat, const Numeric &lon) const;
};

struct Field final
    : FieldMap::Map<Data, SurfaceKey, SurfaceTypeTag, SurfacePropertyTag> {
  //! The ellipsoid used for the surface, in [a, b] in meters
  Vector2 ellipsoid;

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

  friend std::ostream &operator<<(std::ostream &, const Field &);
};

static_assert(
    std::same_as<typename Field::KeyVal, SurfaceKeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");
}  // namespace Surf

using SurfacePoint = Surf::Point;
using SurfaceField = Surf::Field;

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
struct std::formatter<SurfaceTypeTag> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SurfaceTypeTag &v, FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(ctx.out(), "{}{}{}", quote, v.name, quote);
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
      std::format_to(
          ctx.out(),
          R"(, "SurfaceKey": {}, "SurfaceTypeTag": {}, "SurfacePropertyTag": {})",
          v.map<SurfaceKey>().size(),
          v.map<SurfaceTypeTag>().size(),
          v.map<SurfacePropertyTag>().size());
    } else {
      const std::string_view sep = tags.sep(true);
      tags.format(ctx,
                  sep,
                  R"("SurfaceKey": )"sv,
                  v.map<SurfaceKey>(),
                  sep,
                  R"("SurfaceTypeTag": )"sv,
                  v.map<SurfaceTypeTag>(),
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
                   R"("elevation": {}, "temperature": {}, "wind": )",
                   v.elevation,
                   v.temperature);
    tags.format(ctx, v.wind, R"(, "normal": )"sv, v.normal);

    if (tags.short_str) {
      std::format_to(ctx.out(),
                     R"(, "SurfaceTypeTag": {}, "SurfacePropertyTag": {})",
                     v.type.size(),
                     v.prop.size());
    } else {
      const std::string_view sep = tags.sep(true);
      tags.format(ctx,
                  sep,
                  R"("SurfaceTypeTag": )"sv,
                  v.type,
                  sep,
                  R"("SurfacePropertyTag": )"sv,
                  v.prop);
    }

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};
