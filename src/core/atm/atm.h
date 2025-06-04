#pragma once

#include <compare.h>
#include <debug.h>
#include <enumsAtmKey.h>
#include <enumsInterpolationExtrapolation.h>
#include <enumsIsoRatioOption.h>
#include <fieldmap.h>
#include <isotopologues.h>
#include <matpack.h>
#include <operators.h>
#include <quantum_numbers.h>
#include <scattering/properties.h>
#include <species.h>

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>

AtmKey to_wind(const String &);
AtmKey to_mag(const String &);

namespace Atm {

template <typename T>
concept isScatteringSpeciesProperty =
    std::is_same_v<std::remove_cvref_t<T>, ScatteringSpeciesProperty>;

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
                  isQuantumIdentifier<T> or isScatteringSpeciesProperty<T>;

using KeyVal = std::variant<AtmKey,
                            SpeciesEnum,
                            SpeciesIsotope,
                            QuantumIdentifier,
                            ScatteringSpeciesProperty>;

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
  using SpeciesMap        = std::unordered_map<SpeciesEnum, Numeric>;
  using SpeciesIsotopeMap = std::unordered_map<SpeciesIsotope, Numeric>;
  using NlteMap           = std::unordered_map<QuantumIdentifier, Numeric>;
  using ScatteringSpeciesMap =
      std::unordered_map<ScatteringSpeciesProperty, Numeric>;

  SpeciesMap specs{};
  SpeciesIsotopeMap isots{};
  NlteMap nlte{};
  ScatteringSpeciesMap ssprops{};

  Numeric pressure{0};
  Numeric temperature{0};
  Vector3 wind{0, 0, 0};
  Vector3 mag{0, 0, 0};

  Point(const IsoRatioOption isots_key = IsoRatioOption::Builtin);
  Point(Numeric pressure_, Numeric temperature_) : pressure(pressure_), temperature(temperature_) {}
  Point(const Point &)            = default;
  Point(Point &&)                 = default;
  Point &operator=(const Point &) = default;
  Point &operator=(Point &&)      = default;

  Numeric operator[](SpeciesEnum x) const;
  Numeric operator[](const SpeciesIsotope &x) const;
  Numeric operator[](const QuantumIdentifier &x) const;
  Numeric operator[](const ScatteringSpeciesProperty &x) const;
  Numeric operator[](AtmKey x) const;

  Numeric &operator[](SpeciesEnum x);
  Numeric &operator[](const SpeciesIsotope &x);
  Numeric &operator[](const QuantumIdentifier &x);
  Numeric &operator[](const ScatteringSpeciesProperty &x);
  Numeric &operator[](AtmKey x);

  Numeric operator[](const KeyVal &) const;
  Numeric &operator[](const KeyVal &);

  [[nodiscard]] Numeric number_density() const;
  [[nodiscard]] Numeric number_density(const SpeciesEnum &spec) const;
  [[nodiscard]] Numeric number_density(const SpeciesIsotope &spec) const;

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
      else if constexpr (isScatteringSpeciesProperty<T>)
        return x.ssprops.end() not_eq x.ssprops.find(std::forward<T>(k));
    };

    if constexpr (N > 0)
      return has_(*this, std::forward<T>(key)) and
             has(std::forward<Ts>(keys)...);
    else
      return has_(*this, std::forward<T>(key));
  }

  [[nodiscard]] bool contains(const KeyVal &) const;

  [[nodiscard]] Numeric mean_mass() const;
  [[nodiscard]] Numeric mean_mass(SpeciesEnum) const;

  [[nodiscard]] std::vector<KeyVal> keys(bool keep_basic   = true,
                                         bool keep_specs   = true,
                                         bool keep_isots   = true,
                                         bool keep_nlte    = true,
                                         bool keep_ssprops = true) const;

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
      const QuantumIdentifier &band) const;

  void check_and_fix();
};

//! All the field data; if these types grow too much we might want to
//! reconsider...
using FunctionalData = NumericTernaryOperator;
using FieldData      = std::variant<GriddedField3, Numeric, FunctionalData>;

struct FunctionalDataAlwaysThrow {
  std::string error{"Undefined data"};
  Numeric operator()(Numeric, Numeric, Numeric) const {
    ARTS_USER_ERROR("{}", error)
  }
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

  void adjust_interpolation_extrapolation();

  // Allow copy implicitly from all types
  explicit Data(Numeric x);
  Data &operator=(Numeric x);

  explicit Data(GriddedField3 x);
  Data &operator=(GriddedField3 x);

  explicit Data(FunctionalData x);
  Data &operator=(FunctionalData x);

  [[nodiscard]] String data_type() const;

  template <RawDataType T>
  [[nodiscard]] const T &get() const {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr,
                       "Contains {}, not {}",
                       data_type(),
                       Data{T{}}.data_type())
    return *out;
  }

  template <RawDataType T>
  [[nodiscard]] T &get() {
    auto *out = std::get_if<std::remove_cvref_t<T>>(&data);
    ARTS_USER_ERROR_IF(out == nullptr,
                       "Contains {}, not {}",
                       data_type(),
                       Data{T{}}.data_type())
    return *out;
  }

  void rescale(Numeric x);

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
                                   ScatteringSpeciesProperty> {
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
  [[nodiscard]] const std::unordered_map<ScatteringSpeciesProperty, Data> &
  ssprops() const;

  [[nodiscard]] std::unordered_map<QuantumIdentifier, Data> &nlte();
  [[nodiscard]] std::unordered_map<SpeciesEnum, Data> &specs();
  [[nodiscard]] std::unordered_map<SpeciesIsotope, Data> &isots();
  [[nodiscard]] std::unordered_map<AtmKey, Data> &other();
  [[nodiscard]] std::unordered_map<ScatteringSpeciesProperty, Data> &ssprops();

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

  [[nodiscard]] Field gridded(const AscendingGrid &alt,
                              const AscendingGrid &lat,
                              const AscendingGrid &lon) const;
};

struct Xarr {
  Numeric toa{};
  std::vector<Field::KeyVal> keys{};
  AscendingGrid altitudes{};
  AscendingGrid latitudes{};
  AscendingGrid longitudes{};
  Tensor4 data{};

  Xarr(const Field &, std::vector<Field::KeyVal> keys = {});
};

static_assert(
    std::same_as<typename Field::KeyVal, KeyVal>,
    "The order of arguments in the template of which Field inherits from is "
    "wrong.  KeyVal must be defined in the same way for this to work.");

/** Returns true if the pressure is increasing
 * 
 * The pressure must be consistently increasing or decreasing.
 * Otherwise, the function will throw an error.
 * 
 * @param a 
 * @return true 
 * @return false 
 */
bool pressure_is_increasing(const std::span<const Point> &atm);

/** Create a 1D atmosphheric field from a profile
 * 
 * The pressure must be consistently increasing or decreasing.
 * 
 * The atmospheric profile is defined between two iterators.  These
 * must be contiguous in memory with begin < end.
 * 
 * @param atm A span of atmospheric points
 * @param altitudes The altitude grid
 * @param altitude_extrapolation The way that extrapolation should be done
 * @param top_of_atmosphere The top of the atmosphere altitude
 * @return Field 
 */
Field atm_from_profile(const std::span<const Point> &atm,
                       const AscendingGrid &altitudes,
                       const InterpolationExtrapolation &altitude_extrapolation,
                       const Numeric &top_of_atmosphere = NAN);

/** Extends the atmosphere by adding a new point at the given pressure
 * 
 * The new point is created by interpolating the neighbouring points. The altitude
 * grid for this extenstion is from the pressure profile of the atmosphere.
 * 
 * @param atm The profile to extend
 * @param extrapolation_type The type of extrapolation to use
 * @param new_pressure The new pressure to add
 * @param logarithmic If true, the pressure is treated as a logarithmic coordinate, otherwise it is treated as a linear coordinate
 */
void extend_in_pressure(std::vector<Point> &atm,
                        const Numeric &new_pressure,
                        const InterpolationExtrapolation extrapolation_type =
                            InterpolationExtrapolation::Nearest,
                        const bool logarithmic = true);
}  // namespace Atm

using AtmKeyVal         = Atm::KeyVal;
using AtmField          = Atm::Field;
using AtmPoint          = Atm::Point;
using ArrayOfAtmPoint   = Array<AtmPoint>;
using AtmFunctionalData = Atm::FunctionalData;
using AtmData           = Atm::Data;

bool operator==(const AtmKeyVal &, AtmKey);
bool operator==(AtmKey, const AtmKeyVal &);
bool operator==(const AtmKeyVal &, const SpeciesEnum &);
bool operator==(const SpeciesEnum &, const AtmKeyVal &);
bool operator==(const AtmKeyVal &, const QuantumIdentifier &);
bool operator==(const QuantumIdentifier &, const AtmKeyVal &);
bool operator==(const ScatteringSpeciesProperty &, const AtmKeyVal &);

template <>
struct std::formatter<AtmKeyVal> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  [[nodiscard]] std::string to_string(const AtmKeyVal &v) const;

  template <class FmtContext>
  FmtContext::iterator format(const AtmKeyVal &v, FmtContext &ctx) const {
    tags.format(ctx, to_string(v));
    return ctx.out();
  }
};

template <>
struct std::formatter<AtmPoint> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const AtmPoint &v, FmtContext &ctx) const {
    tags.add_if_bracket(ctx, '{');

    const std::string_view sep = tags.sep(true);

    tags.format(ctx,
                R"("pressure": )"sv,
                v.pressure,
                sep,
                R"("temperature": )"sv,
                v.temperature,
                sep,
                R"("mag" :)"sv,
                v.mag,
                sep,
                R"("wind": )"sv,
                v.wind,
                sep);

    if (tags.short_str) {
      tags.format(ctx,
                  R"("SpeciesEnum": )"sv,
                  v.specs.size(),
                  sep,
                  R"("SpeciesIsotope": )"sv,
                  v.isots.size(),
                  sep,
                  R"("QuantumIdentifier": )"sv,
                  v.nlte.size(),
                  sep,
                  R"("ScatteringSpeciesProperty": )"sv,
                  v.ssprops.size());

    } else {
      tags.format(ctx,
                  R"("SpeciesEnum": )"sv,
                  v.specs,
                  sep,
                  R"("SpeciesIsotope": )"sv,
                  v.isots,
                  sep,
                  R"("QuantumIdentifier": )"sv,
                  v.nlte,
                  sep,
                  R"("ScatteringSpeciesProperty": )"sv,
                  v.ssprops);
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

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Atm::Data &v, FmtContext &ctx) const {
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
struct std::formatter<AtmField> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const AtmField &v, FmtContext &ctx) const {
    tags.add_if_bracket(ctx, '{');

    if (tags.short_str) {
      const std::string_view sep = tags.sep();

      tags.format(ctx,
                  R"("top_of_atmosphere": )"sv,
                  v.top_of_atmosphere,
                  sep,
                  R"("Base": )"sv,
                  v.other().size(),
                  sep,
                  R"("SpeciesEnum": )"sv,
                  v.specs().size(),
                  sep,
                  R"("SpeciesIsotope": )"sv,
                  v.isots().size(),
                  sep,
                  R"("QuantumIdentifier": )"sv,
                  v.nlte().size(),
                  sep,
                  R"("ScatteringSpeciesProperty": )"sv,
                  v.ssprops().size());
    } else {
      const std::string_view sep = tags.sep(true);

      tags.format(ctx,
                  R"("top_of_atmosphere": )"sv,
                  v.top_of_atmosphere,
                  sep,
                  R"("Base": )"sv,
                  v.other(),
                  sep,
                  R"("SpeciesEnum": )"sv,
                  v.specs(),
                  sep,
                  R"("SpeciesIsotope": )"sv,
                  v.isots(),
                  sep,
                  R"("QuantumIdentifier": )"sv,
                  v.nlte(),
                  sep,
                  R"("ScatteringSpeciesProperty": )"sv,
                  v.ssprops());
    }

    tags.add_if_bracket(ctx, '}');
    return ctx.out();
  }
};

using SpeciesEnumVectors = std::unordered_map<SpeciesEnum, Vector>;
