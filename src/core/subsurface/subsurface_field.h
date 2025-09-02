#pragma once

#include <enumsInterpolationExtrapolation.h>
#include <enumsSubsurfaceKey.h>
#include <matpack.h>
#include <operators.h>

#include <unordered_map>

struct SubsurfacePropertyTag {
  String name;

  auto operator<=>(const SubsurfacePropertyTag &x) const = default;
};

namespace std {
template <>
struct hash<SubsurfacePropertyTag> {
  std::size_t operator()(const SubsurfacePropertyTag &pp) const noexcept {
    return std::hash<String>{}(pp.name);
  }
};
}  // namespace std

namespace Subsurface {
template <typename T>
concept isSubsurfaceKey = std::same_as<std::remove_cvref_t<T>, SubsurfaceKey>;

template <typename T>
concept isSubsurfacePropertyTag =
    std::same_as<std::remove_cvref_t<T>, SubsurfacePropertyTag>;

template <typename T>
concept KeyType = isSubsurfaceKey<T> or isSubsurfacePropertyTag<T>;

using KeyVal = std::variant<SubsurfaceKey, SubsurfacePropertyTag>;

struct Point {
  Numeric temperature{0};
  Numeric density{0};

  std::unordered_map<SubsurfacePropertyTag, Numeric> prop;

  Point();
  Point(const Point &);
  Point(Point &&) noexcept;
  Point &operator=(const Point &);
  Point &operator=(Point &&) noexcept;

  Numeric operator[](SubsurfaceKey x) const;
  Numeric &operator[](SubsurfaceKey x);

  Numeric operator[](const SubsurfacePropertyTag &x) const;
  Numeric &operator[](const SubsurfacePropertyTag &x);

  Numeric operator[](const KeyVal &) const;
  Numeric &operator[](const KeyVal &);

  template <KeyType T, KeyType... Ts, std::size_t N = sizeof...(Ts)>
  constexpr bool has(T &&key, Ts &&...keys) const {
    if constexpr (N > 0) {
      if constexpr (isSubsurfaceKey<T>) {
        return has(std::forward<Ts>(keys)...);
      } else if constexpr (isSubsurfacePropertyTag<T>) {
        return prop.contains(key) and has(std::forward<Ts>(keys)...);
      } else {
        static_assert(
            isSubsurfacePropertyTag<T> and not isSubsurfacePropertyTag<T>,
            "Unhandled type");
      }
    } else {
      if constexpr (isSubsurfaceKey<T>) {
        return true;
      } else if constexpr (isSubsurfacePropertyTag<T>) {
        return prop.contains(key);
      } else {
        static_assert(
            isSubsurfacePropertyTag<T> and not isSubsurfacePropertyTag<T>,
            "Unhandled type");
      }
    }
  }

  [[nodiscard]] bool contains(const KeyVal &) const;

  [[nodiscard]] Index size() const;
  [[nodiscard]] Index nbasic() const;

  [[nodiscard]] std::vector<KeyVal> keys() const;
};

using FunctionalData = NumericTernaryOperator;
using FieldData      = std::variant<GeodeticField3, Numeric, FunctionalData>;

template <typename T>
concept isGeodeticField3 =
    std::is_same_v<std::remove_cvref_t<T>, GeodeticField3>;

template <typename T>
concept isNumeric = std::is_same_v<std::remove_cvref_t<T>, Numeric>;

template <typename T>
concept isFunctionalDataType =
    std::is_same_v<std::remove_cvref_t<T>, FunctionalData>;

template <typename T>
concept RawDataType =
    isGeodeticField3<T> or isNumeric<T> or isFunctionalDataType<T>;

struct Data {
  FieldData data{FunctionalData{}};
  InterpolationExtrapolation alt_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation alt_low{InterpolationExtrapolation::None};
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

  explicit Data(GeodeticField3 x);
  Data &operator=(GeodeticField3 x);

  explicit Data(FunctionalData x);
  Data &operator=(FunctionalData x);

  [[nodiscard]] String data_type() const;

  template <RawDataType T>
  [[nodiscard]] const T &get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    if (out == nullptr) throw std::runtime_error(data_type());
    return *out;
  }

  template <RawDataType T>
  [[nodiscard]] T &get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    if (out == nullptr) throw std::runtime_error(data_type());
    return *out;
  }

  [[nodiscard]] Numeric at(const Numeric alt,
                           const Numeric lat,
                           const Numeric lon) const;

  [[nodiscard]] Numeric at(const Vector3 pos) const;

  [[nodiscard]] ConstVectorView flat_view() const;

  [[nodiscard]] VectorView flat_view();

  //! Flat weights for the positions in an atmosphere
  [[nodiscard]] std::array<std::pair<Index, Numeric>, 8> flat_weight(
      const Numeric alt, const Numeric lat, const Numeric lon) const;

  [[nodiscard]] std::array<std::pair<Index, Numeric>, 8> flat_weight(
      const Vector3 pos) const;

  [[nodiscard]] bool ok() const;
};

struct Field final {
  std::unordered_map<SubsurfaceKey, Data> other;
  std::unordered_map<SubsurfacePropertyTag, Data> prop;

  Numeric bottom_depth{std::numeric_limits<Numeric>::max()};

  Field();
  Field(const Field &);
  Field(Field &&) noexcept;
  Field &operator=(const Field &);
  Field &operator=(Field &&) noexcept;

  Data &operator[](const SubsurfaceKey &key);
  Data &operator[](const SubsurfacePropertyTag &key);
  Data &operator[](const KeyVal &key);

  const Data &operator[](const SubsurfaceKey &key) const;
  const Data &operator[](const SubsurfacePropertyTag &key) const;
  const Data &operator[](const KeyVal &key) const;

  [[nodiscard]] bool contains(const SubsurfaceKey &key) const;
  [[nodiscard]] bool contains(const KeyVal &key) const;

  //! Compute the values at a single point
  [[nodiscard]] Point at(const Numeric alt,
                         const Numeric lat,
                         const Numeric lon) const;

  //! Compute the values at a single point
  [[nodiscard]] Point at(const Vector3 pos) const;

  [[nodiscard]] Index size() const;
  [[nodiscard]] Index nbasic() const;

  [[nodiscard]] std::vector<KeyVal> keys() const;
};
}  // namespace Subsurface

using SubsurfaceField        = Subsurface::Field;
using SubsurfacePoint        = Subsurface::Point;
using SubsurfaceData         = Subsurface::Data;
using SubsurfaceKeyVal       = Subsurface::KeyVal;
using ArrayOfSubsurfacePoint = Array<SubsurfacePoint>;

template <>
struct std::formatter<SubsurfacePropertyTag> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SubsurfacePropertyTag &v,
                              FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return tags.format(ctx, quote, v.name, quote);
  }
};

template <>
struct xml_io_stream_name<SubsurfaceKeyVal> {
  static constexpr std::string_view name = "SubsurfaceKeyVal"sv;
};

template <>
struct std::formatter<SubsurfacePoint> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  [[nodiscard]] std::string to_string(const SubsurfacePoint &v) const;

  template <class FmtContext>
  FmtContext::iterator format(const SubsurfacePoint &v, FmtContext &ctx) const {
    return tags.format(ctx, to_string(v));
  }
};

template <>
struct std::formatter<SubsurfaceData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SubsurfaceData &v, FmtContext &ctx) const {
    const std::string_view sep = tags.sep();
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.data, sep);
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.alt_upp,
                sep,
                v.alt_low,
                sep,
                v.lat_upp,
                sep,
                v.lat_low,
                sep,
                v.lon_upp,
                sep,
                v.lon_low);
    tags.add_if_bracket(ctx, ']');
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<SubsurfaceField> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  [[nodiscard]] std::string to_string(const SubsurfaceField &v) const;

  template <class FmtContext>
  FmtContext::iterator format(const SubsurfaceField &v, FmtContext &ctx) const {
    return tags.format(ctx, to_string(v));
  }
};
