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
#include "fieldmap.h"
#include "gridded_fields.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "quantum_numbers.h"
#include "species.h"
#include "species_tags.h"

//! A type to name particulates (and let them be type-independent)
struct ParticulatePropertyTag {
  String name;

  auto operator<=>(const ParticulatePropertyTag &) const = default;

  friend std::ostream &operator<<(std::ostream &,
                                  const ParticulatePropertyTag &);
};

namespace std {
template <> struct hash<ParticulatePropertyTag> {
  std::size_t operator()(const ParticulatePropertyTag &pp) const {
    return std::hash<String>{}(pp.name);
  }
};
} // namespace std

namespace Atm {
ENUMCLASS(Key, char, t, p, wind_u, wind_v, wind_w, mag_u, mag_v, mag_w)

template <typename T>
concept isParticulatePropertyTag =
    std::is_same_v<std::remove_cvref_t<T>, ParticulatePropertyTag>;

template <typename T>
concept isArrayOfSpeciesTag =
    std::is_same_v<std::remove_cvref_t<T>, ArrayOfSpeciesTag>;

template <typename T>
concept isQuantumIdentifier =
    std::is_same_v<std::remove_cvref_t<T>, QuantumIdentifier>;

template <typename T>
concept isKey = std::is_same_v<std::remove_cvref_t<T>, Key>;

template <typename T>
concept KeyType = isKey<T> or isArrayOfSpeciesTag<T> or
                  isQuantumIdentifier<T> or isParticulatePropertyTag<T>;

using KeyVal = std::variant<Key, ArrayOfSpeciesTag, QuantumIdentifier,
                            ParticulatePropertyTag>;

template <typename T>
concept ListKeyType = requires(T a) {
  { a.size() } -> matpack::integral;
  { a[0] } -> KeyType;
};

template <typename T>
concept ListOfNumeric = requires(T a) {
  { matpack::mdshape(a) } -> std::same_as<std::array<Index, 1>>;
  { matpack::mdvalue(a, {Index{0}}) } -> std::same_as<Numeric>;
};

struct Point {
  std::unordered_map<ArrayOfSpeciesTag, Numeric> specs{};
  std::unordered_map<QuantumIdentifier, Numeric> nlte{};
  std::unordered_map<ParticulatePropertyTag, Numeric> partp{};

public:
  Numeric pressure{std::numeric_limits<Numeric>::quiet_NaN()};
  Numeric temperature{std::numeric_limits<Numeric>::quiet_NaN()};
  Vector3 wind{0, 0, 0};
  Vector3 mag{std::numeric_limits<Numeric>::quiet_NaN(),
              std::numeric_limits<Numeric>::quiet_NaN(),
              std::numeric_limits<Numeric>::quiet_NaN()};

  template <KeyType T> constexpr Numeric operator[](T &&x) const {
    if constexpr (isArrayOfSpeciesTag<T>) {
      auto y = specs.find(std::forward<T>(x));
      return y == specs.end() ? 0 : y->second;
    } else if constexpr (isParticulatePropertyTag<T>) {
      auto y = partp.find(std::forward<T>(x));
      return y == partp.end() ? 0 : y->second;
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

  template <KeyType T> constexpr Numeric &operator[](T &&x) {
    if constexpr (isArrayOfSpeciesTag<T>) {
      return specs[std::forward<T>(x)];
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte[std::forward<T>(x)];
    } else if constexpr (isParticulatePropertyTag<T>) {
      return partp[std::forward<T>(x)];
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

  Numeric operator[](const KeyVal &) const;
  Numeric &operator[](const KeyVal &);

  template <KeyType T, KeyType... Ts, std::size_t N = sizeof...(Ts)>
  constexpr bool has(T &&key, Ts &&...keys) const {
    const auto has_ = [](auto &x [[maybe_unused]], auto &&k [[maybe_unused]]) {
      if constexpr (isArrayOfSpeciesTag<T>)
        return x.specs.end() not_eq x.specs.find(std::forward<T>(k));
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
  [[nodiscard]] Index npart() const;
  [[nodiscard]] Index nnlte() const;
  [[nodiscard]] static constexpr Index nother() {
    return static_cast<Index>(enumtyps::KeyTypes.size());
  }

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

  void setZero();

  friend std::ostream &operator<<(std::ostream &os, const Point &atm);

  void set(const ArrayOfArrayOfSpeciesTag& sp, const ConstVectorView &x);
};

//! All the field data; if these types grow too much we might want to
//! reconsider...
using FunctionalData = std::function<Numeric(Numeric, Numeric, Numeric)>;
using FieldData = std::variant<GriddedField3, Numeric, FunctionalData>;

ENUMCLASS(Extrapolation, char, None, Zero, Nearest, Linear)

struct FunctionalDataAlwaysThrow {
  std::string error{"Undefined data"};
  Numeric operator()(Numeric, Numeric, Numeric) const { ARTS_USER_ERROR(error) }
};

template <typename T>
concept isGriddedField3 = std::is_same_v<std::remove_cvref_t<T>, GriddedField3>;

template <typename T>
concept isNumeric = std::is_same_v<std::remove_cvref_t<T>, Numeric>;

template <typename T>
concept isFunctionalDataType =
    std::is_same_v<std::remove_cvref_t<T>, FunctionalData>;

template <typename T>
concept RawDataType =
    isGriddedField3<T> or isNumeric<T> or isFunctionalDataType<T>;

//! Hold all atmospheric data
struct Data {
  FieldData data{FunctionalData{FunctionalDataAlwaysThrow{
      "You touched the field but did not set any data"}}};
  Extrapolation alt_upp{Extrapolation::None};
  Extrapolation alt_low{Extrapolation::None};
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
  explicit Data(const GriddedField3 &x) : data(x) {}
  explicit Data(const Numeric &x) : data(x) {}
  explicit Data(const FunctionalData &x) : data(x) {}
  explicit Data(GriddedField3 &&x) : data(std::move(x)) {}
  explicit Data(FunctionalData &&x) : data(std::move(x)) {}

  // Allow copy and move set implicitly from all types
  Data &operator=(const GriddedField3 &x) {
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
  Data &operator=(GriddedField3 &&x) {
    data = std::move(x);
    return *this;
  }
  Data &operator=(FunctionalData &&x) {
    data = std::move(x);
    return *this;
  }

  [[nodiscard]] String data_type() const;

  template <RawDataType T> [[nodiscard]] T get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  template <RawDataType T> [[nodiscard]] T get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  void rescale(Numeric);
};

template <typename T>
concept isData = std::is_same_v<std::remove_cvref_t<T>, Data>;

template <typename T>
concept DataType = RawDataType<T> or isData<T>;

struct Field final : FieldMap::Map<Data, Key, ArrayOfSpeciesTag,
                                   QuantumIdentifier, ParticulatePropertyTag> {
  //! The upper altitude limit of the atmosphere (the atmosphere INCLUDES this
  //! altitude)
  Numeric top_of_atmosphere{std::numeric_limits<Numeric>::lowest()};

  [[nodiscard]] const std::unordered_map<QuantumIdentifier, Data> &nlte() const;
  [[nodiscard]] const std::unordered_map<ArrayOfSpeciesTag, Data> &specs() const;
  [[nodiscard]] const std::unordered_map<Key, Data> &other() const;
  [[nodiscard]] const std::unordered_map<ParticulatePropertyTag, Data> &partp() const;

  [[nodiscard]] std::unordered_map<QuantumIdentifier, Data> &nlte();
  [[nodiscard]] std::unordered_map<ArrayOfSpeciesTag, Data> &specs();
  [[nodiscard]] std::unordered_map<Key, Data> &other();
  [[nodiscard]] std::unordered_map<ParticulatePropertyTag, Data> &partp();

  //! Compute the values at a single point in place
  void at(std::vector<Point> &out, const Vector &alt, const Vector &lat, const Vector &lon) const;

  //! Compute the values at a single point
  [[nodiscard]] std::vector<Point> at(const Vector &alt, const Vector &lat, const Vector &lon) const;

  [[nodiscard]] Index nspec() const;
  [[nodiscard]] Index npart() const;
  [[nodiscard]] Index nnlte() const;
  [[nodiscard]] Index nother() const;

  [[nodiscard]] ArrayOfQuantumIdentifier nlte_keys() const;

  friend std::ostream &operator<<(std::ostream &os, const Field &atm);
};

static_assert(
    std::same_as<typename Field::KeyVal, KeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");

/** A wrapper to fix the input field to the expected format for Field
 *
 * The input must contain all of the Pressure, Latitude, and Longitude grids
 *
 * Throws if anything goes wrong
 *
 * @param[in] gf A gridded field
 * @return GriddedField3 in the Field format
 */
GriddedField3 fix(const GriddedField3 &);

/** A wrapper to fix the input field to the expected format for Field
 *
 * The input must contain 2 of the Pressure, Latitude, and Longitude grids
 *
 * Throws if anything goes wrong
 *
 * @param[in] gf A gridded field
 * @return GriddedField3 in the Field format
 */
GriddedField3 fix(const GriddedField2 &);

/** A wrapper to fix the input field to the expected format for Field
 *
 * The input must contain 1 of the Pressure, Latitude, and Longitude grids
 *
 * Throws if anything goes wrong
 *
 * @param[in] gf A gridded field
 * @return GriddedField3 in the Field format
 */
GriddedField3 fix(const GriddedField1 &);
} // namespace Atm

using AtmField = Atm::Field;
using AtmPoint = Atm::Point;
using ArrayOfAtmPoint = Array<AtmPoint>;
