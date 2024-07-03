#pragma once

#include <matpack.h>
#include <quantum_numbers.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <exception>
#include <format>
#include <functional>
#include <iomanip>
#include <limits>
#include <ostream>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

#include "compare.h"
#include "debug.h"
#include "enums.h"
#include "fieldmap.h"
#include "isotopologues.h"
#include "species.h"

//! A type to name particulates (and let them be type-independent)
struct ParticulatePropertyTag {
  String name;

  auto operator<=>(const ParticulatePropertyTag &) const = default;

  friend std::ostream &operator<<(std::ostream &,
                                  const ParticulatePropertyTag &);
};

namespace std {
template <>
struct hash<ParticulatePropertyTag> {
  std::size_t operator()(const ParticulatePropertyTag &pp) const {
    return std::hash<String>{}(pp.name);
  }
};
}  // namespace std

namespace Atm {
template <typename T>
concept isParticulatePropertyTag =
    std::is_same_v<std::remove_cvref_t<T>, ParticulatePropertyTag>;

template <typename T>
concept isSpecies = std::is_same_v<std::remove_cvref_t<T>, SpeciesEnum>;

template <typename T>
concept isSpeciesIsotope =
    std::is_same_v<std::remove_cvref_t<T>, SpeciesIsotope>;

template <typename T>
concept isQuantumIdentifier =
    std::is_same_v<std::remove_cvref_t<T>, QuantumIdentifier>;

template <typename T>
concept isAtmKey = std::is_same_v<std::remove_cvref_t<T>, AtmKey>;

template <typename T>
concept KeyType = isAtmKey<T> or isSpecies<T> or isSpeciesIsotope<T> or
                  isQuantumIdentifier<T> or isParticulatePropertyTag<T>;

using KeyVal = std::variant<AtmKey,
                            SpeciesEnum,
                            SpeciesIsotope,
                            QuantumIdentifier,
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
  std::unordered_map<SpeciesEnum, Numeric> specs{};
  std::unordered_map<SpeciesIsotope, Numeric> isots{};
  std::unordered_map<QuantumIdentifier, Numeric> nlte{};
  std::unordered_map<ParticulatePropertyTag, Numeric> partp{};

  Numeric pressure{NAN};
  Numeric temperature{NAN};
  Vector3 wind{NAN, NAN, NAN};
  Vector3 mag{NAN, NAN, NAN};

  Point(const IsoRatioOption isots_key = IsoRatioOption::Builtin);
  Point(const Point &)            = default;
  Point(Point &&)                 = default;
  Point &operator=(const Point &) = default;
  Point &operator=(Point &&)      = default;

  template <KeyType T>
  constexpr Numeric operator[](T &&x) const try {
    if constexpr (isSpecies<T>) {
      return specs.at(std::forward<T>(x));
    } else if constexpr (isSpeciesIsotope<T>) {
      return isots.at(std::forward<T>(x));
    } else if constexpr (isParticulatePropertyTag<T>) {
      return partp.at(std::forward<T>(x));
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte.at(std::forward<T>(x));
    } else {
      switch (std::forward<T>(x)) {
        case AtmKey::t:
          return temperature;
        case AtmKey::p:
          return pressure;
        case AtmKey::wind_u:
          return wind[0];
        case AtmKey::wind_v:
          return wind[1];
        case AtmKey::wind_w:
          return wind[2];
        case AtmKey::mag_u:
          return mag[0];
        case AtmKey::mag_v:
          return mag[1];
        case AtmKey::mag_w:
          return mag[2];
      }
      ARTS_USER_ERROR("Cannot reach")
    }
  } catch (std::out_of_range &) {
    if constexpr (isSpecies<T>) {
      ARTS_USER_ERROR("Species VMR not found: ", std::quoted(toString<1>(x)))
    } else if constexpr (isSpeciesIsotope<T>) {
      ARTS_USER_ERROR("Isotopologue ration not found: ",
                      std::quoted(x.FullName()))
    } else if constexpr (isParticulatePropertyTag<T>) {
      ARTS_USER_ERROR("ParticulatePropertyTag value not found: ", x)
    } else if constexpr (isQuantumIdentifier<T>) {
      ARTS_USER_ERROR("Non-LTE level ratio not found: ", x)
    } else {
      ARTS_USER_ERROR("Key not found: \"", x, '"')
    }
  } catch (std::exception &) {
    throw;
  }

  template <KeyType T>
  constexpr Numeric &operator[](T &&x) {
    if constexpr (isSpecies<T>) {
      return specs[std::forward<T>(x)];
    } else if constexpr (isSpeciesIsotope<T>) {
      return isots[std::forward<T>(x)];
    } else if constexpr (isQuantumIdentifier<T>) {
      return nlte[std::forward<T>(x)];
    } else if constexpr (isParticulatePropertyTag<T>) {
      return partp[std::forward<T>(x)];
    } else {
      switch (std::forward<T>(x)) {
        case AtmKey::t:
          return temperature;
        case AtmKey::p:
          return pressure;
        case AtmKey::wind_u:
          return wind[0];
        case AtmKey::wind_v:
          return wind[1];
        case AtmKey::wind_w:
          return wind[2];
        case AtmKey::mag_u:
          return mag[0];
        case AtmKey::mag_v:
          return mag[1];
        case AtmKey::mag_w:
          return mag[2];
      }
      return temperature;  // CANNOT REACH
    }
  }

  Numeric operator[](const KeyVal &) const;
  Numeric &operator[](const KeyVal &);

  template <KeyType T, KeyType... Ts, std::size_t N = sizeof...(Ts)>
  constexpr bool has(T &&key, Ts &&...keys) const {
    const auto has_ = [](auto &x [[maybe_unused]], auto &&k [[maybe_unused]]) {
      if constexpr (isSpecies<T>)
        return x.specs.end() not_eq x.specs.find(std::forward<T>(k));
      else if constexpr (isSpeciesIsotope<T>)
        return x.isots.end() not_eq x.isots.find(std::forward<T>(k));
      else if constexpr (isAtmKey<T>)
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
  [[nodiscard]] Numeric mean_mass(SpeciesEnum) const;

  [[nodiscard]] std::vector<KeyVal> keys() const;

  [[nodiscard]] Index size() const;
  [[nodiscard]] Index nspec() const;
  [[nodiscard]] Index nisot() const;
  [[nodiscard]] Index npart() const;
  [[nodiscard]] Index nnlte() const;
  [[nodiscard]] static constexpr Index nother() {
    return static_cast<Index>(enumsize::AtmKeySize);
  }

  [[nodiscard]] constexpr bool zero_wind() const noexcept {
    return std::all_of(wind.begin(), wind.end(), Cmp::eq(0));
  }

  [[nodiscard]] constexpr bool zero_mag() const noexcept {
    return std::all_of(mag.begin(), mag.end(), Cmp::eq(0));
  }

  [[nodiscard]] bool is_lte() const noexcept;

  [[nodiscard]] std::pair<Numeric, Numeric> levels(
      const QuantumIdentifier &band) const {
    return {operator[](band.LowerLevel()), operator[](band.UpperLevel())};
  }

  void check_and_fix();

  friend std::ostream &operator<<(std::ostream &os, const Point &atm);
};

//! All the field data; if these types grow too much we might want to
//! reconsider...
using FunctionalData = std::function<Numeric(Numeric, Numeric, Numeric)>;
using FieldData      = std::variant<GriddedField3, Numeric, FunctionalData>;

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
  InterpolationExtrapolation alt_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation alt_low{InterpolationExtrapolation::None};
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
  explicit Data(const RawDataType auto &x) : data(x) {}
  explicit Data(RawDataType auto &&x) : data(std::move(x)) {}
  Data &operator=(const RawDataType auto &x) {
    data = x;
    return *this;
  }
  Data &operator=(RawDataType auto &&x) {
    data = std::move(x);
    return *this;
  }

  [[nodiscard]] String data_type() const;

  template <RawDataType T>
  [[nodiscard]] T get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  template <RawDataType T>
  [[nodiscard]] T get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr, "Does not contain correct type")
    return *out;
  }

  void rescale(Numeric);

  [[nodiscard]] Vector at(const Vector &alt,
                          const Vector &lat,
                          const Vector &lon) const;

  [[nodiscard]] Numeric at(const Numeric alt,
                           const Numeric lat,
                           const Numeric lon) const;

  [[nodiscard]] ExhaustiveConstVectorView flat_view() const;

  [[nodiscard]] ExhaustiveVectorView flat_view();

  //! Flat weights for the positions in an atmosphere
  [[nodiscard]] Array<std::array<std::pair<Index, Numeric>, 8>> flat_weights(
      const Vector &alt, const Vector &lat, const Vector &lon) const;
};

template <typename T>
concept isData = std::is_same_v<std::remove_cvref_t<T>, Data>;

template <typename T>
concept DataType = RawDataType<T> or isData<T>;

struct Field final : FieldMap::Map<Data,
                                   AtmKey,
                                   SpeciesEnum,
                                   SpeciesIsotope,
                                   QuantumIdentifier,
                                   ParticulatePropertyTag> {
  //! The upper altitude limit of the atmosphere (the atmosphere INCLUDES this
  //! altitude)
  Numeric top_of_atmosphere{std::numeric_limits<Numeric>::lowest()};

  Field(const IsoRatioOption isots_key = IsoRatioOption::Builtin);
  Field(const Field &)                = default;
  Field(Field &&) noexcept            = default;
  Field &operator=(const Field &)     = default;
  Field &operator=(Field &&) noexcept = default;

  [[nodiscard]] const std::unordered_map<QuantumIdentifier, Data> &nlte() const;
  [[nodiscard]] const std::unordered_map<SpeciesEnum, Data> &specs() const;
  [[nodiscard]] const std::unordered_map<SpeciesIsotope, Data> &isots() const;
  [[nodiscard]] const std::unordered_map<AtmKey, Data> &other() const;
  [[nodiscard]] const std::unordered_map<ParticulatePropertyTag, Data> &partp()
      const;

  [[nodiscard]] std::unordered_map<QuantumIdentifier, Data> &nlte();
  [[nodiscard]] std::unordered_map<SpeciesEnum, Data> &specs();
  [[nodiscard]] std::unordered_map<SpeciesIsotope, Data> &isots();
  [[nodiscard]] std::unordered_map<AtmKey, Data> &other();
  [[nodiscard]] std::unordered_map<ParticulatePropertyTag, Data> &partp();

  //! Compute the values at a single point in place
  void at(std::vector<Point> &out,
          const Vector &alt,
          const Vector &lat,
          const Vector &lon) const;

  //! Compute the values at a single point
  [[nodiscard]] std::vector<Point> at(const Vector &alt,
                                      const Vector &lat,
                                      const Vector &lon) const;

  //! Compute the values at a single point
  [[nodiscard]] Point at(const Numeric alt,
                         const Numeric lat,
                         const Numeric lon) const;

  //! Compute the values at a single point
  [[nodiscard]] Point at(const Vector3 pos) const;

  [[nodiscard]] Index nspec() const;
  [[nodiscard]] Index nisot() const;
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

std::ostream &operator<<(std::ostream &os, const Array<Point> &a);
}  // namespace Atm

using AtmKeyVal         = Atm::KeyVal;
using AtmField          = Atm::Field;
using AtmPoint          = Atm::Point;
using ArrayOfAtmPoint   = Array<AtmPoint>;
using AtmFunctionalData = Atm::FunctionalData;

bool operator==(const AtmKeyVal &, AtmKey);
bool operator==(AtmKey, const AtmKeyVal &);
bool operator==(const AtmKeyVal &, const SpeciesEnum &);
bool operator==(const SpeciesEnum &, const AtmKeyVal &);
bool operator==(const AtmKeyVal &, const QuantumIdentifier &);
bool operator==(const QuantumIdentifier &, const AtmKeyVal &);
bool operator==(const AtmKeyVal &, const ParticulatePropertyTag &);
bool operator==(const ParticulatePropertyTag &, const AtmKeyVal &);

std::ostream &operator<<(std::ostream &, const AtmKeyVal &);

template <>
struct std::formatter<ParticulatePropertyTag> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts> &...xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U> &x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const ParticulatePropertyTag &v,
                              FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(ctx.out(), "{}{}{}", quote, v.name, quote);
  }
};

template <>
struct std::formatter<AtmPoint> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts> &...xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U> &x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const AtmPoint &v, FmtContext &ctx) const {
    tags.add_if_bracket(ctx, '{');

    std::formatter<Vector3> vec3;
    std::format_to(ctx.out(),
                   R"("pressure": {}, "temperature": {}, "mag": )",
                   v.pressure,
                   v.temperature);
    vec3.format(v.mag, ctx);
    std::format_to(ctx.out(), R"("wind": )");
    vec3.format(v.wind, ctx);

    if (tags.short_str) {
      std::format_to(
          ctx.out(),
          R"(, "SpeciesEnum": {}, , "SpeciesIsotope": {}, "QuantumIdentifier": {}, "ParticulatePropertyTag": {})",
          v.specs.size(),
          v.isots.size(),
          v.nlte.size(),
          v.partp.size());

    } else {
      std::formatter<std::unordered_map<SpeciesEnum, Numeric>> specs{};
      std::formatter<std::unordered_map<SpeciesIsotope, Numeric>> isots{};
      std::formatter<std::unordered_map<QuantumIdentifier, Numeric>> nlte{};
      std::formatter<std::unordered_map<ParticulatePropertyTag, Numeric>>
          partp{};

      std::format_to(ctx.out(), R"(,
"SpeciesEnum": )");
      specs.format(v.specs, ctx);

      std::format_to(ctx.out(), R"(,
"SpeciesIsotope": )");
      isots.format(v.isots, ctx);

      std::format_to(ctx.out(), R"(,
"QuantumIdentifier": )");
      nlte.format(v.nlte, ctx);

      std::format_to(ctx.out(), R"(
"ParticulatePropertyTag": )");
      partp.format(v.partp, ctx);
    }

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};

template <>
struct std::formatter<Atm::FunctionalData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts> &...xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U> &x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Atm::FunctionalData &,
                              FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return std::format_to(ctx.out(), "{}functional-data{}", quote, quote);
  }
};

template <>
struct std::formatter<Atm::Data> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts> &...xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U> &x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Atm::Data &v, FmtContext &ctx) const {
    std::formatter<Atm::FieldData> data{};
    std::formatter<InterpolationExtrapolation> interp{};
    make_compat(interp, data);

    const std::string_view sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    data.format(v.data, ctx);
    std::format_to(ctx.out(), "{}", sep);
    tags.add_if_bracket(ctx, '[');

    interp.format(v.alt_upp, ctx);
    std::format_to(ctx.out(), "{}", sep);
    interp.format(v.alt_low, ctx);
    std::format_to(ctx.out(), "{}", sep);
    interp.format(v.lat_upp, ctx);
    std::format_to(ctx.out(), "{}", sep);
    interp.format(v.lat_low, ctx);
    std::format_to(ctx.out(), "{}", sep);
    interp.format(v.lon_upp, ctx);
    std::format_to(ctx.out(), "{}", sep);
    interp.format(v.lon_low, ctx);

    if (tags.bracket) std::format_to(ctx.out(), "]]");
    return ctx.out();
  }
};

template <>
struct std::formatter<AtmField> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts> &...xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U> &x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const AtmField &v, FmtContext &ctx) const {
    tags.add_if_bracket(ctx, '{');

    if (tags.short_str) {
      std::format_to(
          ctx.out(),
          R"("top_of_atmosphere": {}, "Base": {}, "SpeciesEnum": {}, "SpeciesIsotope": {}, "QuantumIdentifier": {}, "ParticulatePropertyTag": {})",
          v.top_of_atmosphere,
          v.other().size(),
          v.specs().size(),
          v.isots().size(),
          v.nlte().size(),
          v.partp().size());
    } else {
      std::formatter<std::unordered_map<QuantumIdentifier, Atm::Data>> nlte{};
      std::formatter<std::unordered_map<SpeciesEnum, Atm::Data>> specs{};
      std::formatter<std::unordered_map<SpeciesIsotope, Atm::Data>> isots{};
      std::formatter<std::unordered_map<AtmKey, Atm::Data>> other{};
      std::formatter<std::unordered_map<ParticulatePropertyTag, Atm::Data>>
          partp{};

      std::format_to(ctx.out(),
                     R"("top_of_atmosphere": {},
"AtmKey": )",
                     v.top_of_atmosphere);
      other.format(v.other(), ctx);

      std::format_to(ctx.out(), R"(,
"SpeciesEnum": )");
      specs.format(v.specs(), ctx);

      std::format_to(ctx.out(), R"(,
"SpeciesIsotope": )");
      isots.format(v.isots(), ctx);

      std::format_to(ctx.out(), R"(,
"QuantumIdentifier": )");
      nlte.format(v.nlte(), ctx);

      std::format_to(ctx.out(), R"(,
"ParticulatePropertyTag": )");
      partp.format(v.partp(), ctx);
    }

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};
