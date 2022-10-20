#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <ostream>
#include <type_traits>
#include <utility>
#include <variant>

#include "artstime.h"
#include "debug.h"
#include "enums.h"
#include "gridded_fields.h"
#include "matpack.h"
#include "matpackI.h"
#include "matpackIV.h"
#include "species.h"
#include "species_tags.h"

namespace Atm {
ENUMCLASS(Key,
          char,
          temperature,
          pressure,
          wind_u,
          wind_v,
          wind_w,
          mag_u,
          mag_v,
          mag_w)

template <typename T>
concept isArrayOfSpeciesTag =
    std::is_same_v<std::remove_cvref_t<T>, ArrayOfSpeciesTag>;

template <typename T>
concept isKey = std::is_same_v<std::remove_cvref_t<T>, Key>;

template <typename T>
concept KeyType = isKey<T> or isArrayOfSpeciesTag<T>;

class Point {
  std::map<ArrayOfSpeciesTag, Numeric> species_content{};
  Numeric pressure{std::numeric_limits<Numeric>::min()};
  Numeric temperature{std::numeric_limits<Numeric>::min()};
  std::array<Numeric, 3> wind{0., 0., 0.};
  std::array<Numeric, 3> mag{0., 0., 0.};

  template <typename... Ts>
  void internal_set(KeyType auto&& lhs, auto&& rhs, Ts&&... ts) {
    set(std::forward<decltype(lhs)>(lhs), std::forward<decltype(rhs)>(rhs));
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

 public:
  template <typename... Ts>
  Point(Ts&&... ts) {
    static_assert((sizeof...(Ts) % 2) == 0, "Uneven number of inputs");
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

  Numeric operator[](KeyType auto&& x) const {
    using T = decltype(x);
    if constexpr (isArrayOfSpeciesTag<T>) {
      auto y = species_content.find(std::forward<T>(x));
      return y == species_content.end() ? 0 : y->second;
    } else {
      switch (std::forward<T>(x)) {
        case Key::temperature:
          return temperature;
        case Key::pressure:
          return pressure;
        case Key::wind_u:
          return wind[0];
        case Key::wind_v:
          return wind[1];
        case Key::wind_w:
          return wind[2];
        case Key::mag_u:
          return mag[0];
          break;
        case Key::mag_v:
          return mag[1];
          break;
        case Key::mag_w:
          return mag[2];
        case Key::FINAL: {
        }
      }
      return 0;
    }
  }

  [[nodiscard]] constexpr auto P() const { return pressure; }

  [[nodiscard]] constexpr auto T() const { return temperature; }

  [[nodiscard]] constexpr auto M() const { return mag; }

  [[nodiscard]] constexpr auto W() const { return wind; }

  void set(KeyType auto&& x, Numeric y) {
    using T = decltype(x);
    ARTS_USER_ERROR_IF(
        std::isnan(y) or std::isinf(y), "Bad input NaN or Inf: ", y)

    if constexpr (isArrayOfSpeciesTag<T>) {
      species_content[x] = y;
    } else {
      switch (std::forward<T>(x)) {
        case Key::temperature:
          ARTS_USER_ERROR_IF(y <= 0, "Bad temperature: ", y)
          temperature = y;
          break;
        case Key::pressure:
          ARTS_USER_ERROR_IF(y <= 0, "Bad pressure: ", y)
          pressure = y;
          break;
        case Key::wind_u:
          wind[0] = y;
          break;
        case Key::wind_v:
          wind[1] = y;
          break;
        case Key::wind_w:
          wind[2] = y;
          break;
        case Key::mag_u:
          mag[0] = y;
          break;
        case Key::mag_v:
          mag[1] = y;
          break;
        case Key::mag_w:
          mag[2] = y;
          break;
        case Key::FINAL:
          break;
      }
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const Point& atm);
};

//! All the field data; if these types grow too much we might want to reconsider...
using FunctionalData = std::function<Numeric(Time,Numeric,Numeric,Numeric)>;
using FieldData = std::variant<GriddedField4, Tensor4, Numeric, FunctionalData>;

template <typename T>
concept isGriddedField4 = std::is_same_v<std::remove_cvref_t<T>, GriddedField4>;

template <typename T>
concept isTensor4 = std::is_same_v<std::remove_cvref_t<T>, Tensor4>;

template <typename T>
concept isNumeric = std::is_same_v<std::remove_cvref_t<T>, Numeric>;

template <typename T>
concept isFunctionalDataType = std::is_same_v<std::remove_cvref_t<T>, FunctionalData>;

template <typename T>
concept RawDataType = isGriddedField4<T> or isNumeric<T> or isFunctionalDataType<T>;

template <typename T>
concept DataType = RawDataType<T> or isTensor4<T>;

class Field {
  std::map<Key, FieldData> other{};
  std::map<ArrayOfSpeciesTag, FieldData> specs{};

  //! The below only exist if regularized is true
  bool regularized{false};
  ArrayOfTime time{};
  Vector alt{}, lat{}, lon{};

  template <typename... Ts>
  void internal_set(KeyType auto&& lhs, RawDataType auto&& rhs, Ts&&... ts) {
    using T = decltype(lhs);
    if constexpr (isArrayOfSpeciesTag<T>) {
      specs[std::forward<T>(lhs)] = std::forward<decltype(rhs)>(rhs);
    } else {
      other[std::forward<T>(lhs)] = std::forward<decltype(rhs)>(rhs);
    }
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

 public:
  template <typename... Ts>
  Field(Ts&&... ts) {
    static_assert((sizeof...(Ts) % 2) == 0, "Uneven number of inputs");
    if constexpr (sizeof...(Ts)) internal_set(std::forward<Ts>(ts)...);
  }

  /*  // Make this unavailable?  It would make it easier to type erase in the future
  [[nodiscard]] const FieldData& get(KeyType auto&& x) const {
    using T = decltype(x);
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs.at(std::forward<T>(x));
    } else {
      return other.at(std::forward<T>(x));
    }
  }
  */

  [[nodiscard]] Shape<4> regularized_shape() const {
    return {time.nelem(), alt.nelem(), lat.nelem(), lon.nelem()};
  }

  void set(KeyType auto&& x, DataType auto&& y) {
    using T = decltype(x);
    using U = decltype(y);

    if constexpr (isTensor4<U>) {
      ARTS_USER_ERROR_IF(not regularized,
                         "Field needs to be regularized to set Tensor4 data")
      ARTS_USER_ERROR_IF(y.shape() not_eq regularized_shape(),
                         "The shape is wrong.  The input field is shape ",
                         y.shape(),
                         " but we have a regularized shape of ",
                         regularized_shape())
    } else {
      ARTS_USER_ERROR_IF(
          regularized,
          "Expects a regularized field (i.e., Tensor4-style gridded)")
    }

    if constexpr (isGriddedField4<U>) {
      ARTS_USER_ERROR_IF(
          "Time" not_eq y.get_grid_name(0) or
              "Altitude" not_eq y.get_grid_name(1) or
              "Latitude" not_eq y.get_grid_name(2) or
              "Longitude" not_eq y.get_grid_name(3),
          "The grids should be [Time x Altitude x Latitude x Longitude] but it is [",
          y.get_grid_name(0),
          " x ",
          y.get_grid_name(1),
          " x ",
          y.get_grid_name(2),
          " x ",
          y.get_grid_name(3), ']')
    }

    if constexpr (isArrayOfSpeciesTag<T>) {
      specs[std::forward<T>(x)] = std::move(y);
    } else {
      other[std::forward<T>(x)] = std::move(y);
    }
  }

  //! Regularizes the calculations so that all data is on a single grid
  Field& regularize(const ArrayOfTime& times,
                    const Vector& altitudes,
                    const Vector& latitudes,
                    const Vector& longitudes);

  //! Compute the values at a single point
  [[nodiscard]] Point at(Time time_point,
                         Numeric alt_point,
                         Numeric lat_point,
                         Numeric lon_point) const;

  friend std::ostream& operator<<(std::ostream& os, const Field& atm);
};
}  // namespace Atm

using AtmField = Atm::Field;
using AtmPoint = Atm::Point;