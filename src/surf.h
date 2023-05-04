#pragma once

#include "enums.h"
#include "fieldmap.h"
#include "gridded_fields.h"
#include "matpack_concepts.h"
#include "matpack_constexpr.h"
#include "mystring.h"

#include <concepts>
#include <limits>
#include <ostream>
#include <type_traits>
#include <variant>

using Vector2 = matpack::matpack_constant_data<Numeric, 2>;
using Vector3 = matpack::matpack_constant_data<Numeric, 3>;

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
template <> struct hash<SurfaceTypeTag> {
  std::size_t operator()(const SurfaceTypeTag &pp) const noexcept {
    return std::hash<String>{}(pp.name);
  }
};

template <> struct hash<SurfacePropertyTag> {
  std::size_t operator()(const SurfacePropertyTag &pp) const noexcept {
    return std::hash<String>{}(pp.name);
  }
};
} // namespace std

namespace Surf {
ENUMCLASS(Key, char, h, t, wind_u, wind_v, wind_w)

template <typename T>
concept isKey = std::same_as<std::remove_cvref_t<T>, Key>;

template <typename T>
concept isSurfaceTypeTag = std::same_as<std::remove_cvref_t<T>, SurfaceTypeTag>;

template <typename T>
concept isSurfacePropertyTag = std::same_as<std::remove_cvref_t<T>, SurfacePropertyTag>;

template <typename T>
concept KeyType = isKey<T> or isSurfaceTypeTag<T> or isSurfacePropertyTag<T>;

using KeyVal = std::variant<Key, SurfaceTypeTag, SurfacePropertyTag>;

struct Point {
  Numeric elevation;
  Numeric temperature;
  Vector3 wind;
  Vector2 normal;
  std::unordered_map<SurfaceTypeTag, Numeric> type;
  std::unordered_map<SurfacePropertyTag, Numeric> prop;

  template <KeyType Key> Numeric &operator[](Key &&x) {
    if constexpr (isKey<Key>) {
      switch (std::forward<Key>(x)) {
      case Surf::Key::h:
        return elevation;
      case Surf::Key::t:
        return temperature;
      case Surf::Key::wind_u:
        return wind[0];
      case Surf::Key::wind_v:
        return wind[1];
      case Surf::Key::wind_w:
        return wind[2];
      case Surf::Key::FINAL:
        return temperature;
      }
    } else if constexpr (isSurfaceTypeTag<Key>) {
      return type[std::forward<Key>(x)];
    } else if constexpr (isSurfacePropertyTag<Key>) {
      return prop[std::forward<Key>(x)];
    }
  }

  template <KeyType Key> Numeric operator[](Key &&x) const {
    if constexpr (isKey<Key>) {
      switch (std::forward<Key>(x)) {
      case Surf::Key::h:
        return elevation;
      case Surf::Key::t:
        return temperature;
      case Surf::Key::wind_u:
        return wind[0];
      case Surf::Key::wind_v:
        return wind[1];
      case Surf::Key::wind_w:
        return wind[2];
      case Surf::Key::FINAL:
        return std::numeric_limits<Numeric>::signaling_NaN();
      }
    } else if constexpr (isSurfaceTypeTag<Key>) {
      return type.at(std::forward<Key>(x));
    } else if constexpr (isSurfacePropertyTag<Key>) {
      return prop.at(std::forward<Key>(x));
    }
  }

  Numeric &operator[](const KeyVal &x);

  Numeric operator[](const KeyVal &x) const;

  [[nodiscard]] std::vector<KeyVal> keys() const;

  [[nodiscard]] Index nelem() const;
  [[nodiscard]] Index ntype() const;
  [[nodiscard]] static constexpr Index nother() {
    return static_cast<Index>(enumtyps::KeyTypes.size());
  }

  template <KeyType T> constexpr bool contains(T &&k) const {
    if constexpr (isSurfaceTypeTag<T>)
      return type.contains(std::forward<T>(k));
    if constexpr (isSurfacePropertyTag<T>)
      return prop.contains(std::forward<T>(k));
    else if constexpr (isKey<T>)
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
using FieldData = std::variant<GriddedField2, Numeric, FunctionalData>;

struct FunctionalDataAlwaysThrow {
  std::string error{"Undefined data"};
  Numeric operator()(Numeric, Numeric) const { ARTS_USER_ERROR(error) }
};

ENUMCLASS(Extrapolation, char, None, Zero, Nearest, Linear)

//! Hold all atmospheric data
struct Data {
  FieldData data{FunctionalData{FunctionalDataAlwaysThrow{
      "You touched the field but did not set any data"}}};
  Extrapolation lat_upp{Extrapolation::None};
  Extrapolation lat_low{Extrapolation::None};
  Extrapolation lon_upp{Extrapolation::None};
  Extrapolation lon_low{Extrapolation::None};

  // Standard
  Data() = default;
  Data(const Data &) = default;
  Data(Data &&) = default;
  Data &operator=(const Data &) = default;
  Data &operator=(Data &&) = default;

  // Allow copy and move construction implicitly from all types
  explicit Data(const GriddedField2 &x) : data(x) {}
  explicit Data(const Numeric &x) : data(x) {}
  explicit Data(const FunctionalData &x) : data(x) {}
  explicit Data(GriddedField2 &&x) : data(std::move(x)) {}
  explicit Data(FunctionalData &&x) : data(std::move(x)) {}

  [[nodiscard]] String data_type() const;

  // Allow copy and move set implicitly from all types
  Data &operator=(const GriddedField2 &x) {
    data = x;
    return *this;
  }
  Data &operator=(const Numeric &x) {
    data = x;
    return *this;
  }
  Data &operator=(const FunctionalData &x) {
    data = x;
    return *this;
  }
  Data &operator=(GriddedField2 &&x) {
    data = std::move(x);
    return *this;
  }
  Data &operator=(FunctionalData &&x) {
    data = std::move(x);
    return *this;
  }

  template <typename T> [[nodiscard]] T get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  template <typename T> [[nodiscard]] T get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  void rescale(Numeric);
};

struct Field final : FieldMap::Map<Data, Key, SurfaceTypeTag, SurfacePropertyTag> {
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
  [[nodiscard]] Point at(Numeric lat, Numeric lon, Vector2 ellipsoid) const;

  /** Get the normal of the surface at a given position
   *
   * @param lat The latitude
   * @param lon The longitude
   * @param alt The altitude of the point.  If NaN, it is computed first
   * @return Vector2 [za, aa]
   */
  [[nodiscard]] Vector2
  normal(Vector2 ellipsoid, Numeric lat, Numeric lon,
         Numeric alt = std::numeric_limits<Numeric>::quiet_NaN()) const;

  [[nodiscard]] Numeric single_value(const KeyVal& key, Numeric lat, Numeric lon) const;

  friend std::ostream &operator<<(std::ostream &, const Field &);
};

static_assert(
    std::same_as<typename Field::KeyVal, KeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");
} // namespace Surf

using SurfacePoint = Surf::Point;
using SurfaceField = Surf::Field;
