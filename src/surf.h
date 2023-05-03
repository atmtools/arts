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

namespace std {
template <> struct hash<SurfaceTypeTag> {
  std::size_t operator()(const SurfaceTypeTag &pp) const {
    return std::hash<String>{}(pp.name);
  }
};
} // namespace std

namespace Surf {
ENUMCLASS(Key, char, z, t, wind_u, wind_v, wind_w)

template <typename T>
concept isKey = std::same_as<std::remove_cvref_t<T>, Key>;

template <typename T>
concept isSurfaceTypeTag = std::same_as<std::remove_cvref_t<T>, SurfaceTypeTag>;

template <typename T>
concept KeyType = isKey<T> or isSurfaceTypeTag<T>;

using KeyVal = std::variant<Key, SurfaceTypeTag>;

struct Point {
  Numeric altitude;
  Numeric temperature;
  Vector3 wind;
  Vector2 normal;
  std::unordered_map<SurfaceTypeTag, Numeric> type;

  template <KeyType Key> Numeric &operator[](Key &&x) {
    if constexpr (isKey<Key>) {
      switch (std::forward<Key>(x)) {
      case Surf::Key::z:
        return altitude;
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
    }
  }

  template <KeyType Key> Numeric operator[](Key &&x) const {
    if constexpr (isKey<Key>) {
      switch (std::forward<Key>(x)) {
      case Surf::Key::z:
        return altitude;
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

  template <KeyType T, KeyType... Ts, std::size_t N = sizeof...(Ts)>
  constexpr bool has(T &&key, Ts &&...keys) const {
    const auto has_ = [](auto &x [[maybe_unused]], auto &&k [[maybe_unused]]) {
      if constexpr (isSurfaceTypeTag<T>)
        return x.specs.end() not_eq x.specs.find(std::forward<T>(k));
      else if constexpr (isKey<T>)
        return true;
      else return false;
    };

    if constexpr (N > 0)
      return has_(*this, std::forward<T>(key)) and
             has(std::forward<Ts>(keys)...);
    else
      return has_(*this, std::forward<T>(key));
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

struct Field final : FieldMap::Map<Data, Key, SurfaceTypeTag> {
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

  friend std::ostream &operator<<(std::ostream &, const Field &);
};

static_assert(
    std::same_as<typename Field::KeyVal, KeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");
} // namespace Surf

using SurfPoint = Surf::Point;
using SurfField = Surf::Field;
