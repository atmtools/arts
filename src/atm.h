#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

#include "arts_constants.h"
#include "compare.h"
#include "debug.h"
#include "enums.h"
#include "gridded_fields.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "quantum_numbers.h"
#include "species.h"
#include "species_tags.h"

namespace Atm {
ENUMCLASS(Key,
          char,
          t,
          p,
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

using KeyVal = std::variant<Key, ArrayOfSpeciesTag, QuantumIdentifier>;

template <typename T>
concept ListKeyType = requires (T a) {
  { a.size() } -> matpack::integral;
  { a[0] } -> KeyType;
};

template <typename T>
concept ListOfNumeric = requires (T a) {
  { matpack::mdshape(a) } -> std::same_as<std::array<Index, 1>>;
  { matpack::mdvalue(a, {Index{0}}) } -> std::same_as<Numeric>;
};

struct Point {
private:
  template <ListKeyType T, ListOfNumeric U, typename... Ts, size_t N = sizeof...(Ts)>
  void internal_set(T &&lhs, U &&rhs, Ts &&...ts) {
    ARTS_ASSERT(static_cast<Index>(lhs.size()) == matpack::mdshape(rhs)[0])
    for (decltype(lhs.size()) i = 0; i < lhs.size(); i++)
      this->operator[](lhs[i]) = mdvalue(rhs, {static_cast<Index>(i)});
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  template <KeyType T, typename U, typename... Ts, size_t N = sizeof...(Ts)>
  void internal_set(T &&lhs, U &&rhs, Ts &&...ts) {
    this->operator[](std::forward<T>(lhs)) = std::forward<U>(rhs);
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  std::unordered_map<ArrayOfSpeciesTag, Numeric> specs{};
  std::unordered_map<QuantumIdentifier, Numeric> nlte{};

public:
  Numeric pressure{0};
  Numeric temperature{0};
  std::array<Numeric, 3> wind{0., 0., 0.};
  std::array<Numeric, 3> mag{0., 0., 0.};

  template <typename... Ts, std::size_t N = sizeof...(Ts)>
  Point(Ts&&... ts) requires((N % 2) == 0) {
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  void reserve(Index nspec, Index nnlte) {
    specs.reserve(nspec);
    nlte.reserve(nnlte);
  }

  template<KeyType T>
  Numeric operator[](T&& x) const {
    if constexpr (isArrayOfSpeciesTag<T>) {
      auto y = specs.find(std::forward<T>(x));
      return y == specs.end() ? 0 : y->second;
    } else if constexpr (isQuantumIdentifier<T>) {
      auto y = nlte.find(std::forward<T>(x));
      return y == nlte.end() ? 0 : y->second;
    } else {
      switch (std::forward<T>(x)) {
        case Key::t:
          return temperature;
        case Key::p:
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
      ARTS_USER_ERROR("Cannot reach")
    }
  }

  template<KeyType T>
  Numeric& operator[](T&& x) {
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs[std::forward<T>(x)];
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte[std::forward<T>(x)];
    } else {
      switch (std::forward<T>(x)) {
        case Key::t:
          return temperature;
        case Key::p:
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

  Numeric operator[](Species::Species x) const noexcept;

  template <KeyType T, KeyType... Ts, std::size_t N = sizeof...(Ts)>
  constexpr bool has(T &&key, Ts &&...keys) const {
    const auto has_ = [](auto &x [[maybe_unused]],
                         auto &&k [[maybe_unused]]) {
      if constexpr (isArrayOfSpeciesTag<T>)
        return x.specs.end() not_eq
               x.specs.find(std::forward<T>(k));
      else if constexpr (isKey<T>)
        return true;
      else if constexpr (isQuantumIdentifier<T>)
        return x.nlte.end() not_eq x.nlte.find(std::forward<T>(k));
    };

    if constexpr (N > 0)
      return has_(*this, std::forward<T>(key)) and
             has(std::forward<Ts>(keys)...);
    else
      return has_(*this, std::forward<T>(key));
  }

  [[nodiscard]] Numeric mean_mass() const;

  [[nodiscard]] std::vector<KeyVal> keys() const;

  [[nodiscard]] Index nelem() const;
  [[nodiscard]] Index nspec() const;
  [[nodiscard]] Index nnlte() const;
  [[nodiscard]] static constexpr Index nother() {return static_cast<Index>(enumtyps::KeyTypes.size());}

  [[nodiscard]] constexpr bool zero_wind() const noexcept {
    return std::all_of(wind.begin(), wind.end(), Cmp::eq(0));
  }

  [[nodiscard]] constexpr bool zero_mag() const noexcept {
    return std::all_of(mag.begin(), mag.end(), Cmp::eq(0));
  }

  [[nodiscard]] bool is_lte() const noexcept;

  [[nodiscard]] std::pair<Numeric, Numeric>
  levels(const QuantumIdentifier &band) const {
    return {operator[](band.LowerLevel()), operator[](band.UpperLevel())};
  }

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

  // Standard
  Data() = default;
  Data(const Data&) = default;
  Data(Data&&) = default;
  Data& operator=(const Data&) = default;
  Data& operator=(Data&&) = default;

  // Allow copy and move construction implicitly from all types
  Data(const GriddedField3& x) : data(x) {}
  Data(const Tensor3& x) : data(x) {}
  Data(const Numeric& x) : data(x) {}
  Data(const FunctionalData& x) : data(x) {}
  Data(GriddedField3&& x) : data(std::move(x)) {}
  Data(Tensor3&& x) : data(std::move(x)) {}
  Data(FunctionalData&& x) : data(std::move(x)) {}

  // Allow copy and move set implicitly from all types
  Data& operator=(const GriddedField3& x) {data=x; return *this;}
  Data& operator=(const Tensor3& x) {data=x; return *this;}
  Data& operator=(const Numeric& x) {data=x; return *this;}
  Data& operator=(const FunctionalData& x) {data=x; return *this;}
  Data& operator=(GriddedField3&& x) {data=std::move(x); return *this;}
  Data& operator=(Tensor3&& x) {data=std::move(x); return *this;}
  Data& operator=(FunctionalData&& x) {data=std::move(x); return *this;}

  [[nodiscard]] constexpr bool need_hydrostatic() const noexcept {
    return alt_low == Extrapolation::Hydrostatic or
        alt_upp == Extrapolation::Hydrostatic;
  }

  [[nodiscard]] std::tuple<Numeric, Numeric, Numeric> hydrostatic_coord(Numeric alt, Numeric lat, Numeric lon) const;

  [[nodiscard]] String data_type() const;

  template <typename T>
  [[nodiscard]] T get() const {
    auto* out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  template <typename T>
  [[nodiscard]] T get() {
    auto* out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  void rescale(Numeric);
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

struct Field {
  template <KeyType T, RawDataType U, typename... Ts, std::size_t N = sizeof...(Ts)>
  void internal_set(T&& lhs, U&& rhs, Ts&&... ts) {
    this->operator[](std::forward<T>(lhs)) = std::forward<U>(rhs);
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  [[nodiscard]] Point internal_fitting(Numeric alt_point, Numeric lat_point, Numeric lon_point) const;

  std::unordered_map<Key, Data> other{};
  std::unordered_map<ArrayOfSpeciesTag, Data> specs{};
  std::unordered_map<QuantumIdentifier, Data> nlte{};
 
  //! Grid if regularized
  std::array<Vector, 3> grid{};

  //! The below only exist if regularized is true
  bool regularized{false};

  //! The upper altitude limit of the atmosphere (the atmosphere INCLUDES this altitude)
  Numeric top_of_atmosphere{std::numeric_limits<Numeric>::lowest()};

  template <typename... Ts, std::size_t N = sizeof...(Ts)>
  Field(Ts&&... ts) requires((N % 2) == 0) {
    if constexpr (N > 0) internal_set(std::forward<Ts>(ts)...);
  }

  [[nodiscard]] std::array<Index, 3> regularized_shape() const;

  template <KeyType T> Data &operator[](T &&x) {
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs[std::forward<T>(x)];
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte[std::forward<T>(x)];
    } else {
      return other[std::forward<T>(x)];
    }
  }

  template <KeyType T> const Data &operator[](T &&x) const {
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs.at(std::forward<T>(x));
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte.at(std::forward<T>(x));
    } else {
      return other.at(std::forward<T>(x));
    }
  }

  template <KeyType T> void erase_key(T &&x) {
    if constexpr (isArrayOfSpeciesTag<T>) {
      specs.erase(std::forward<T>(x));
    } else if constexpr (isQuantumIdentifier<T>) {
      nlte.erase(std::forward<T>(x));
    } else {
      other.erase(std::forward<T>(x));
    }
  }

  const Data &operator[](const KeyVal &) const;

  Data &operator[](const KeyVal &);

  //! Regularizes the calculations so that all data is on a single grid
  Field &regularize(const Vector &, const Vector &, const Vector &);

  //! Estimate the old atmospheric dim (ARTS-2 atmosphere_dim)
  [[deprecated, nodiscard]] Index old_atmosphere_dim_est() const;

  //! Checks that exactly dim number of dimensions have size greater than 1
  [[nodiscard]] bool regularized_atmosphere_dim(Index dim) const;

  //! Compute the values at a single point
  [[nodiscard]] Point at(Numeric, Numeric, Numeric, const FunctionalData& g=FunctionalDataAlwaysThrow{}) const;

  template <KeyType T, KeyType... Ts, std::size_t N = sizeof...(Ts)>
  bool has(T &&key, Ts &&...keys) const {
    const auto has_this_key = [this] (auto&& k) {
      if constexpr (isArrayOfSpeciesTag<T>)
        return specs.end() not_eq specs.find(std::forward<T>(k));
      else if constexpr (isKey<T>)
        return other.end() not_eq other.find(std::forward<T>(k));
      else if constexpr (isQuantumIdentifier<T>)
        return nlte.end() not_eq nlte.find(std::forward<T>(k));
    };

    if constexpr (N > 0)
      return has_this_key(std::forward<T>(key)) and has(std::forward<Ts>(keys)...);
    else
      return has_this_key(std::forward<T>(key));
  }

  [[nodiscard]] std::vector<KeyVal> keys() const;

  [[nodiscard]] Index nelem() const;
  [[nodiscard]] Index nspec() const;
  [[nodiscard]] Index nnlte() const;
  [[nodiscard]] Index nother() const;

  void throwing_check() const;

  [[nodiscard]] ArrayOfQuantumIdentifier nlte_keys() const;
  [[nodiscard]] ArrayOfArrayOfSpeciesTag spec_keys() const;

  friend std::ostream &operator<<(std::ostream &os, const Field &atm);
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

/** Extracts all the VMR from a regularized Field
 * 
 *  
 * @param[in] atm An atmosospheric field
 * @param[in] specs A list of atmospheric species
 * @return Tensor4 All gridded VMR
 */
Tensor4 extract_specs_content(const Field& atm, const ArrayOfArrayOfSpeciesTag& specs);
} // namespace Atm

using AtmField = Atm::Field;
using AtmPoint = Atm::Point;
using ArrayOfAtmPoint = Array<AtmPoint>;
