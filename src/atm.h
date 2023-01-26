#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>
#include <variant>

#include "arts_constants.h"
#include "debug.h"
#include "enums.h"
#include "gridded_fields.h"
#include "matpack.h"
#include "matpackI.h"
#include "matpackIV.h"
#include "quantum_numbers.h"
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
concept isQuantumIdentifier =
    std::is_same_v<std::remove_cvref_t<T>, QuantumIdentifier>;

template <typename T>
concept isKey = std::is_same_v<std::remove_cvref_t<T>, Key>;

template <typename T>
concept KeyType = isKey<T> or isArrayOfSpeciesTag<T> or isQuantumIdentifier<T>;

class Point {
  std::map<ArrayOfSpeciesTag, Numeric> species_content{};
  std::map<QuantumIdentifier, Numeric> nlte{};

  //FIXME: Remove
  std::shared_ptr<std::map<QuantumIdentifier, Numeric>> nlte_energy{};

public:
  Numeric pressure{std::numeric_limits<Numeric>::min()};
  Numeric temperature{std::numeric_limits<Numeric>::min()};
  std::array<Numeric, 3> wind{0., 0., 0.};
  std::array<Numeric, 3> mag{0., 0., 0.};

  template <KeyType T, typename U, typename... Ts>
  void internal_set(T&& lhs, U&& rhs, Ts&&... ts) {
    set(std::forward<T>(lhs), std::forward<U>(rhs));
    if constexpr (sizeof...(Ts) > 0) internal_set(std::forward<Ts>(ts)...);
  }

  template <KeyType T, typename... Ts>
  constexpr bool internal_has(T&& key, Ts&&... keys) const {
    bool has_this = has(std::forward<T>(key));
    if constexpr (sizeof...(Ts) > 0)
      return has_this and internal_has(std::forward<Ts>(keys)...);
    else
      return has_this;
  }

 public:
  template <typename... Ts, std::size_t N = sizeof...(Ts)>
  Point(Ts&&... ts) requires((N % 2) == 0) {
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  template<KeyType T>
  Numeric operator[](T&& x) const {
    if constexpr (isArrayOfSpeciesTag<T>) {
      auto y = species_content.find(std::forward<T>(x));
      return y not_eq species_content.end() ? 0 : y->second;
    } else if constexpr (isQuantumIdentifier<T>) {
      auto y = nlte.find(std::forward<T>(x));
      return y not_eq nlte.end() ? 0 : y->second;
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
        case Key::mag_v:
          return mag[1];
        case Key::mag_w:
          return mag[2];
        case Key::FINAL: {
        }
      }
      return temperature; // CANNOT REACH
    }
  }

  // FIXME: This data should be elsewhere???
  [[nodiscard]] Numeric energy_level(const QuantumIdentifier& x) const {
    ARTS_USER_ERROR_IF(
        not nlte_energy or nlte.size() not_eq nlte_energy->size(),
        "Incorrect NLTE energy term")
    auto y = nlte_energy->find(x);
    return y == nlte_energy->end() ? 0 : y->second;
  }

  // FIXME: This data should be elsewhere???
  void set_energy_level(const decltype(nlte_energy)& x) {
    ARTS_USER_ERROR_IF(x and x->size() and x->size() not_eq nlte.size(),
                       "Mismatching sizes")
    nlte_energy = x;
  }

  template<KeyType T>
  void set(T&& x, Numeric y) {
    if constexpr (isArrayOfSpeciesTag<T>) {
      species_content[x] = y;
    } else if constexpr (isQuantumIdentifier<T>) {
      nlte[x] = y;
    } else {
      switch (std::forward<T>(x)) {
        case Key::temperature:
          temperature = y;
          break;
        case Key::pressure:
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

  template<KeyType T>
  constexpr bool has(T&& key) const {
    if constexpr (isArrayOfSpeciesTag<T>)
      return species_content.end() not_eq
             species_content.find(std::forward<T>(key));
    else if constexpr (isKey<T>)
      return true;
    else if constexpr (isQuantumIdentifier<T>)
      return nlte.end() not_eq nlte.find(std::forward<T>(key));
  }

  template <typename... Ts>
  constexpr bool has_data(Ts&&... keys) const {
    if constexpr (sizeof...(Ts))
      return internal_has(std::forward<Ts>(keys)...);
    else
      return true;
  }

  [[nodiscard]] Numeric mean_mass() const;

  friend std::ostream& operator<<(std::ostream& os, const Point& atm);
};

//! All the field data; if these types grow too much we might want to reconsider...
using FunctionalData = std::function<Numeric(Numeric, Numeric, Numeric)>;
using FieldData = std::variant<GriddedField3, Tensor3, Numeric, FunctionalData>;

ENUMCLASS(Extrapolation, char, None, Zero, Nearest, Linear, Hydrostatic)

struct FunctionalDataAlwaysThrow {
  std::string error{"Undefined data"};
  Numeric operator()(Numeric, Numeric, Numeric) const { ARTS_USER_ERROR(error) }
};

//! Hold all atmospheric data
struct Data {
  FieldData data{FunctionalData{FunctionalDataAlwaysThrow{"You touched the field but did not set any data"}}};
  Extrapolation alt_upp{Extrapolation::None};
  Extrapolation alt_low{Extrapolation::None};
  Extrapolation lat_upp{Extrapolation::None};
  Extrapolation lat_low{Extrapolation::None};
  Extrapolation lon_upp{Extrapolation::None};
  Extrapolation lon_low{Extrapolation::None};

  [[nodiscard]] constexpr bool need_hydrostatic() const noexcept {
    return alt_low == Extrapolation::Hydrostatic or
        alt_upp == Extrapolation::Hydrostatic;
  }

  [[nodiscard]] std::tuple<Numeric, Numeric, Numeric> hydrostatic_coord(Numeric alt, Numeric lat, Numeric lon) const;
};

template <typename T>
concept isData = std::is_same_v<std::remove_cvref_t<T>, Data>;

template <typename T>
concept isGriddedField3 = std::is_same_v<std::remove_cvref_t<T>, GriddedField3>;

template <typename T>
concept isTensor3 = std::is_same_v<std::remove_cvref_t<T>, Tensor3>;

template <typename T>
concept isNumeric = std::is_same_v<std::remove_cvref_t<T>, Numeric>;

template <typename T>
concept isFunctionalDataType =
    std::is_same_v<std::remove_cvref_t<T>, FunctionalData>;

template <typename T>
concept RawDataType =
    isGriddedField3<T> or isNumeric<T> or isFunctionalDataType<T>;

template <typename T>
concept DataType = RawDataType<T> or isTensor3<T> or isData<T>;

class Field {
  std::map<Key, Data> other{};
  std::map<ArrayOfSpeciesTag, Data> specs{};
  std::map<QuantumIdentifier, Data> nlte{};

  //FIXME: Remove
  std::shared_ptr<std::map<QuantumIdentifier, Numeric>> nlte_energy{};

  //! The below only exist if regularized is true
  bool regularized{false};
  Vector alt{}, lat{}, lon{};

  //! The upper altitude limit of the atmosphere (the atmosphere INCLUDES this altitude)
  Numeric space_alt{std::numeric_limits<Numeric>::lowest()};

  template <KeyType T, RawDataType U, typename... Ts>
  void internal_set(T&& lhs, U&& rhs, Ts&&... ts) {
    set(std::forward<T>(lhs), std::forward<U>(rhs));
    if constexpr (sizeof...(Ts) > 0) internal_set(std::forward<Ts>(ts)...);
  }

  template <KeyType T, KeyType... Ts>
  bool internal_has(T&& key, Ts&&... keys) const {
    bool has_this = has(std::forward<T>(key));
    if constexpr (sizeof...(Ts))
      return has_this and internal_has(std::forward<Ts>(keys)...);
    else
      return has_this;
  }

  template<KeyType T>
  Data& get(T&& x) {
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs[std::forward<T>(x)];
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte[std::forward<T>(x)];
    } else {
      return other[std::forward<T>(x)];
    }
  }

  template<KeyType T>
  const Data& get(T&& x) const {
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs[std::forward<T>(x)];
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte[std::forward<T>(x)];
    } else {
      return other[std::forward<T>(x)];
    }
  }

 public:
  template <typename... Ts, std::size_t N = sizeof...(Ts)>
  Field(Ts&&... ts) requires((N % 2) == 0) {
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  void top_of_atmosphere(Numeric x) {space_alt = x;}
  [[nodiscard]] Numeric top_of_atmosphere() const {return space_alt;}

  [[nodiscard]] Shape<3> regularized_shape() const {
    return {alt.nelem(), lat.nelem(), lon.nelem()};
  }

  template<KeyType T, DataType U>
  void set(T &&x, U &&y) {
    if constexpr (isData<U>) {
      set_altitude_limits(x, y.alt_low, y.alt_upp);
      set_latitude_limits(x, y.lat_low, y.lat_upp);
      set_longitude_limits(x, y.lon_low, y.lon_upp);
      std::visit([this, x](auto&& d){this -> set(x, d);}, y.data);
    } else {
      if constexpr (isTensor3<U>) {
          ARTS_USER_ERROR_IF(
              not regularized,
              "Field needs to be regularized to set Tensor3 data")
          ARTS_USER_ERROR_IF(y.shape() not_eq regularized_shape(),
                             "The shape is wrong.  The input field is shape ",
                             y.shape(), " but we have a regularized shape of ",
                             regularized_shape())
      } else {
          ARTS_USER_ERROR_IF(
              regularized,
              "Expects a regularized field (i.e., Tensor3-style gridded)")
      }

      if constexpr (isGriddedField3<U>) {
          ARTS_USER_ERROR_IF("Pressure" not_eq y.get_grid_name(0) or
                                 "Latitude" not_eq y.get_grid_name(1) or
                                 "Longitude" not_eq y.get_grid_name(2),
                             "The grids should be [Pressure x Latitude x "
                             "Longitude] but it is [",
                             y.get_grid_name(0), " x ", y.get_grid_name(1),
                             " x ", y.get_grid_name(2), ']')
      }

      get(std::forward<T>(x)).data = std::move(y);
    }
  }

  template<KeyType T>
  Data& operator[](T&& x) { return get(x); }

  template<KeyType T>
  const Data& operator[](T&& x) const { return get(x); }

  template<KeyType T>
  void set_altitude_limits(T &&x, Extrapolation low,
                           Extrapolation upp) {
    auto &y = get(std::forward<T>(x));
    y.alt_low = low;
    y.alt_upp = upp;
  }

  template<KeyType T>
  void set_latitude_limits(T &&x, Extrapolation low,
                           Extrapolation upp) {
    auto &y = get(std::forward<T>(x));
    y.lat_low = low;
    y.lat_upp = upp;
  }

  template<KeyType T>
  void set_longitude_limits(T &&x, Extrapolation low,
                            Extrapolation upp) {
    auto &y = get(std::forward<T>(x));
    y.lon_low = low;
    y.lon_upp = upp;
  }

  // FIXME: This data should be elsewhere???
  [[nodiscard]] Numeric energy_level(const QuantumIdentifier &x) const {
    ARTS_USER_ERROR_IF(not nlte_energy or
                           nlte.size() not_eq nlte_energy->size(),
                       "Incorrect NLTE energy term")
    auto y = nlte_energy->find(x);
    return y == nlte_energy->end() ? 0 : y->second;
  }

  // FIXME: This data should be elsewhere???
  void set_energy_level(const QuantumIdentifier &x, Numeric y) {
    if (not nlte_energy)
      nlte_energy = std::make_shared<std::map<QuantumIdentifier, Numeric>>();
    nlte_energy->operator[](x) = y;
  }

  //! Regularizes the calculations so that all data is on a single grid
  Field &regularize(const Vector &, const Vector &, const Vector &);

  //! Compute the values at a single point
  [[nodiscard]] Point at(Numeric, Numeric, Numeric, const FunctionalData& g=FunctionalDataAlwaysThrow{}) const;

  template <KeyType T>
  bool has(T &&key) const {
    if constexpr (isArrayOfSpeciesTag<T>)
      return specs.end() not_eq specs.find(std::forward<T>(key));
    else if constexpr (isKey<T>)
      return other.end() not_eq other.find(std::forward<T>(key));
    else if constexpr (isQuantumIdentifier<T>)
      return nlte.end() not_eq nlte.find(std::forward<T>(key));
  }

  template <typename... Ts> bool has_data(Ts &&...keys) const {
    if constexpr (sizeof...(Ts))
      return internal_has(std::forward<Ts>(keys)...);
    else
      return true;
  }

  friend std::ostream &operator<<(std::ostream &os, const Field &atm);

private:
  //! Fits an atmospheric point without hydrostatic balance
  friend Point main_fitting(const Field& field, Numeric alt_point, Numeric lat_point, Numeric lon_point);
};

/** A wrapper to fix the input field to the expected format for Field
 * 
 * The input must contain all of the Pressure, Latitude, and Longitude grids
 *
 * Throws if anything goes wrong
 *
 * @param[in] gf A gridded field
 * @return GriddedField3 in the Field format
 */
GriddedField3 fix(const GriddedField3&);

/** A wrapper to fix the input field to the expected format for Field
 * 
 * The input must contain 2 of the Pressure, Latitude, and Longitude grids
 *
 * Throws if anything goes wrong
 *
 * @param[in] gf A gridded field
 * @return GriddedField3 in the Field format
 */
GriddedField3 fix(const GriddedField2&);

/** A wrapper to fix the input field to the expected format for Field
 * 
 * The input must contain 1 of the Pressure, Latitude, and Longitude grids
 *
 * Throws if anything goes wrong
 *
 * @param[in] gf A gridded field
 * @return GriddedField3 in the Field format
 */
GriddedField3 fix(const GriddedField1&);

class SimpleHydrostaticExpansion {
  Numeric x0;
  Numeric h0;
  Numeric H;

public:
  SimpleHydrostaticExpansion(Numeric x, Numeric h, Numeric T, Numeric mu,
                             Numeric g)
      : x0(x), h0(h), H(Constant::R * T / (mu * g)) {}

  Numeric operator()(Numeric alt) const {
    return x0 * std::exp(-(alt - h0) / H);
  }
};
} // namespace Atm

using AtmField = Atm::Field;
using AtmPoint = Atm::Point;
