#pragma once

#include <array.h>
#include <debug.h>
#include <enumsQuantumNumberType.h>
#include <isotopologues.h>
#include <nonstd.h>
#include <rational.h>

#include <algorithm>
#include <compare>
#include <cstddef>
#include <istream>
#include <limits>
#include <ostream>
#include <string_view>
#include <utility>
#include <vector>

constexpr Index quantum_number_error_value = -999'999'999;

namespace Quantum::Number {
//! Three tags, S: str, I: index, H: half-index
enum class QuantumNumberValueType : char { S, I, H };
std::ostream& operator<<(std::ostream&, QuantumNumberValueType);

//! Holds string values but can only hold sizeof(Index) long values
struct StringValue {
  static constexpr std::size_t N = 16;
  std::array<char, N> x{'\0'};  // First is \0 to be "empty"

  //! Returns the value in such a way that no \0 remains in the view
  [[nodiscard]] std::string_view val() const noexcept;

  //! Default initializer to a zero-first solution
  constexpr StringValue() = default;

  //! Set to expected value from a view of something at most N-char long
  explicit StringValue(std::string_view s);

  std::strong_ordering operator<=>(const StringValue& sv) const;
};

//! Holds integer values
struct IntegerValue {
  Index x{std::numeric_limits<Index>::lowest()};

  //! Returns the value as a rational
  [[nodiscard]] constexpr Rational val() const noexcept { return x; }

  constexpr std::strong_ordering operator<=>(const IntegerValue& i) const {
    return x <=> i.x;
  }
};

//! Holds half integer values, but only its denominator
struct HalfIntegerValue {
  Index x{std::numeric_limits<Index>::lowest()};

  //! Returns the value as a Rational, keeping that this is a half-integer
  [[nodiscard]] constexpr Rational val() const noexcept { return {x, 2}; }

  constexpr auto operator<=>(const HalfIntegerValue& h) const = default;
};

/** Common value type of a given quantum number type 
 * 
 * Guards against integer/half-integer comparisons failing for bad reasons
 * and allows IO to not need some information 
 *
 * Keep new entries to a 1-to-1 key-to-key match, and sort them the same as QuantumNumberType for ease of reading
 *
 * @param type A quantum number type
 * @return constexpr QuantumNumberValueType The common type
 */
QuantumNumberValueType common_value_type(QuantumNumberType type) noexcept;

/** Return a common type between a and b
 * 
 * If they are the same, returns the value
 * If either is H and the other is I, returns H
 * Otherwise returns FINAL as an error flag
 *
 * @param a A value type
 * @param b Another value type
 * @return constexpr QuantumNumberValueType with FINAL as error flag
 */
QuantumNumberValueType common_value_type(QuantumNumberValueType a,
                                         QuantumNumberValueType b) noexcept;

//! A union of the three type of values we need to consider
union ValueHolder {
  StringValue s;
  IntegerValue i;
  HalfIntegerValue h;

  ValueHolder(QuantumNumberValueType t) noexcept;
  ValueHolder(QuantumNumberType t) noexcept;

  ValueHolder(const ValueHolder&)                = default;
  ValueHolder(ValueHolder&&) noexcept            = default;
  ValueHolder& operator=(const ValueHolder&)     = default;
  ValueHolder& operator=(ValueHolder&&) noexcept = default;
};

/** A complete description of a value, its type and value
 * 
 * Intended to be returned from IO operations
 */
struct ValueDescription {
  QuantumNumberValueType type;
  ValueHolder val;

  ValueDescription(QuantumNumberValueType t) noexcept : type(t), val(t) {}
  ValueDescription(const ValueDescription&)                = default;
  ValueDescription(ValueDescription&&) noexcept            = default;
  ValueDescription& operator=(const ValueDescription&)     = default;
  ValueDescription& operator=(ValueDescription&&) noexcept = default;

  //! Debug output only
  friend std::ostream& operator<<(std::ostream& os, ValueDescription x);
};

constexpr std::strong_ordering cmp(std::strong_ordering&& a,
                                   std::strong_ordering&& b) {
  return a == std::strong_ordering::equal ? b : a;
}

/** The values of two levels
 * 
 * Its ValueDescription constructor ensures that we have two valid 
 * types and converts integer input in one to half-integer in case
 * there's a 'mismatch'
 */
struct TwoLevelValueHolder {
  ValueHolder upp;
  ValueHolder low;

  //! Constructor to ensure two ValueDescription have the same type, for IO purposes
  TwoLevelValueHolder(ValueDescription u,
                      ValueDescription l,
                      QuantumNumberType t);

  TwoLevelValueHolder(QuantumNumberType t) noexcept : upp(t), low(t) {}
  TwoLevelValueHolder(const TwoLevelValueHolder&)                = default;
  TwoLevelValueHolder(TwoLevelValueHolder&&) noexcept            = default;
  TwoLevelValueHolder& operator=(const TwoLevelValueHolder&)     = default;
  TwoLevelValueHolder& operator=(TwoLevelValueHolder&&) noexcept = default;

  [[nodiscard]] std::strong_ordering order(const TwoLevelValueHolder& tv,
                                           QuantumNumberValueType t) const;
};

/** Takes a rational and determine which type of quantum number it is,
 * returning this information or throwing a runtime error if there's
 * an error
 *
 * @param r_ A rational
 * @return constexpr ValueDescription 
 */
[[nodiscard]] ValueDescription value_holder(Rational r_);

[[nodiscard]] ValueDescription value_holder(std::string_view s);

/** Returns a rational if possible or RATIONAL_UNDEFINED otherwise
 * 
 * Spaces are not allowed and results in RATIONAL_UNDEFINED being returns
 *
 * Also only accepts rationals that match how quantum numbers are defined,
 * so half or full integers, in case a decimal string is given.  Note also
 * that only a single decimal value can be given
 *
 * If the rational is not a full integer or a half-integer, returns RATIONAL_UNDEFINED
 *
 * @param s Some view of a string
 * @return constexpr Rational 
 */
[[nodiscard]] Rational cast_qnrat(std::string_view s) noexcept;

/** Returns a value description for the quantum number
 * 
 * Note that several branches can throw as the input is assumed to be from
 * the user
 *
 * @param s Some view of a string
 * @return constexpr ValueDescription 
 */
[[nodiscard]]  ValueDescription value_holder(std::string_view s,
                                                      QuantumNumberType t) ;

//! Struct that converts to bool automatically but allows checking both energy levels matching status
struct LevelMatch {
  bool upp{true};
  bool low{true};

  //! Convert automatically to bool so exact matches are easy, non-explicit by design
  constexpr operator bool() const noexcept { return upp and low; }
};

/** Count all space-separated items in s
 * 
 * Example: "X 1 3123 1/3 1,,,,,2 " returns 5
 *
 * @param s Any set of characters
 * @return constexpr Index The number of space-separated items in s
 */
 Index count_items(std::string_view s) noexcept;

/** Strips spaces at the end of x before returning it
 * 
 * @param x any string view
 * @return constexpr std::string_view stripped
 */
 std::string_view rstrip(std::string_view x);

/** Strips spaces at the beginning x before returning it
 * 
 * @param x any string view
 * @return constexpr std::string_view stripped
 */
 std::string_view lstrip(std::string_view x);

/** Strips spaces at the beginning and end of x before returning it
 * 
 * @param x any string view
 * @return constexpr std::string_view stripped
 */
 std::string_view strip(std::string_view x);

/** Get a view of a number of space-separated items from the list
 * 
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=0, n=3 returns "X    1   3123"
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=2, n=3 returns "3123 1/3 1,,,,,2"
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=4, n=3 returns "1,,,,,2"
 * Example: "X    1   3123 1/3 1,,,,,2  " with i=5, n=3 returns ""
 *
 * @tparam n Number of items
 * @param s Any set of characters
 * @param i The first item from the original list in the list of n items
 * @param n The length of the list of items
 * @return constexpr std::string_view 
 */
std::string_view items(std::string_view s, std::size_t i, std::size_t n=1) noexcept;

//! A complete quantum number value with type information
struct Value {
  QuantumNumberType type;
  TwoLevelValueHolder qn;

   std::strong_ordering operator<=>(const Value& v) const;

   Value(QuantumNumberType t = QuantumNumberType::term)
      : type(t), qn(type) {}
  Value(const Value&)                = default;
  Value(Value&&) noexcept            = default;
  Value& operator=(const Value&)     = default;
  Value& operator=(Value&&) noexcept = default;

   Value(QuantumNumberType t, Rational upp_, Rational low_);

  //! Default constructor from some string of values
   Value(std::string_view s);

  //! Returns the upper quantum number rational if it exists or an undefined
  [[nodiscard]]  Rational upp() const noexcept;

  //! Returns the lower quantum number rational if it exists or an undefined
  [[nodiscard]] Rational low() const noexcept;

  //! Returns the upper quantum number string copy
  [[nodiscard]] String str_upp() const noexcept;

  //! Returns the lower quantum number string copy
  [[nodiscard]] String str_low() const noexcept;

  //! Legacy way to swap the values between two Values
  void swap_values(Value& x);

  //! Set level value
  void set(std::string_view s, bool upp);

  /** Returns a description of whether both levels match
   * 
   * Note that the LevelMatch type should automatically transform to bool
   * so there's no need for extra work if this is your target question
   *
   * @param other Another value
   * @return constexpr LevelMatch
   */
  [[nodiscard]] LevelMatch level_match(Value other) const noexcept;

  //! Standard output
  friend std::ostream& operator<<(std::ostream& os, Value x);

  //! Standard input
  friend std::istream& operator>>(std::istream& is, Value& x);

  bofstream& write(bofstream& bof) const;

  bifstream& read(bifstream& bif);

  [[nodiscard]] bool good() const;
};

//! Status of comparing two lists that are supposedly of some type
enum class CheckValue : char { Full, AinB, BinA, Miss };

//! Level-by-level version of CheckValue
struct CheckMatch {
  CheckValue upp{CheckValue::Full};
  CheckValue low{CheckValue::Full};

  //! Convert automatically to bool so exact matches are easy, non-explicit by design
  constexpr operator bool() const noexcept {
    return upp == CheckValue::Full and low == CheckValue::Full;
  }
};

//! Updates old by what a new check says it should be
 CheckValue update(CheckValue val, CheckValue res) noexcept ;

//! Updates old by what a new check says it should be
 CheckMatch update(CheckMatch val, CheckValue res) noexcept ;

//! Updates old by what a new check says it should be
 CheckMatch update(CheckMatch val, CheckMatch res) noexcept;

/** Checks if an array of types is sorted
 * 
 * @tparam N Number of types
 * @param types Array of types
 * @return true if it is sorted
 * @return false if it is not sorted
 */
template <size_t N>
constexpr bool is_sorted(
    const std::array<QuantumNumberType, N>& types) noexcept {
  for (size_t i = 1; i < N; i++)
    if (not(types[i - 1] < types[i])) return false;
  return true;
}

//! A list of many quantum numbers.  Should always remain sorted
struct ValueList {
  Array<Value> values;

  //! Internal sort function.  Should be called whenever new items are created
  void sort_by_type();

  //! Internal check function.  Remember to sort by type before calling this
  [[nodiscard]] bool has_unique_increasing_types() const;

  //! From text
  explicit ValueList(std::string_view s, bool legacy = false);

  //! From legacy text
  ValueList(std::string_view upp, std::string_view low);

  std::strong_ordering operator<=>(const ValueList& v) const;

  //! From values (resorted)
  explicit ValueList(Array<Value> values_);

  //! Empty
  ValueList() : values(0) {}

  //! For iterators
  Array<Value>::iterator begin() { return values.begin(); }
  Array<Value>::iterator end() { return values.end(); }
  Array<Value>::iterator cbegin() { return values.begin(); }
  Array<Value>::iterator cend() { return values.end(); }
  [[nodiscard]] Array<Value>::const_iterator begin() const {
    return values.begin();
  }
  [[nodiscard]] Array<Value>::const_iterator end() const {
    return values.end();
  }
  [[nodiscard]] Array<Value>::const_iterator cbegin() const {
    return values.cbegin();
  }
  [[nodiscard]] Array<Value>::const_iterator cend() const {
    return values.cend();
  }

  //! Should always be called before this object is handed to another user
  void finalize();

  //! Return number of quantum numbers
  [[nodiscard]] Index size() const ARTS_NOEXCEPT { return values.size(); }

  //! Finds whether two ValueList describe completely different sets of quantum numbers (e.g., local vs global)
  [[nodiscard]] bool perpendicular(const ValueList& that) const ARTS_NOEXCEPT;

  //! Returns whether all the Types are part of the list, the types must be sorted
  template <typename... Types>
  [[nodiscard]] bool has(Types... ts) const ARTS_NOEXCEPT {
    static_assert(sizeof...(Types) > 0);

    ARTS_ASSERT(is_sorted(std::array{QuantumNumberType(ts)...}))

    auto ptr = cbegin();
    auto end = cend();
    for (QuantumNumberType t : {QuantumNumberType(ts)...}) {
      ptr = std::find_if(ptr, end, [t](auto& x) { return x.type == t; });
      if (ptr == end) return false;
    }
    return true;
  }

  //! Returns the value of the Type (assumes it exist)
  const Value& operator[](QuantumNumberType t) const ARTS_NOEXCEPT;

  //! Legacy manipulation operator access
  Value& operator[](Index i) { return values.at(i); }

  //! Add for manipulation
  Value& add(QuantumNumberType t);

  //! Add for manipulation
  Value& add(Value v);

  //! Sets the value if it exists or adds it otherwise
  void set(Value v);

  //! Set a value in value list
  void set(Index i, std::string_view upp, std::string_view low);

  //! Returns upper and lower matching status
  [[nodiscard]] CheckMatch check_match(const ValueList& other) const
      ARTS_NOEXCEPT;

  //! ouptut stream if all values
  friend std::ostream& operator<<(std::ostream& os, const ValueList& vl);

  //! input stream must have pre-set size
  friend std::istream& operator>>(std::istream& is, ValueList& vl);

  //! Add a type without sorting (WARNING, many things might break if you don't sort in the end)
  void add_type_wo_sort(QuantumNumberType);

  [[nodiscard]] bool good() const;

  void reserve(Size);

  template <typename... Ts>
  auto& emplace_back(Ts&&... xs) {
    return values.emplace_back(std::forward<Ts>(xs)...);
  }
};

ValueList from_hitran(std::string_view upp, std::string_view low);

//! A logical struct for local quantum numbers
struct LocalState {
  ValueList val{};

  std::strong_ordering operator<=>(const LocalState& l) const;
  bool operator==(const LocalState& l) const ;
  bool operator!=(const LocalState& l) const;

  LocalState() = default;

  template <typename... Values>
  LocalState(Values... vals) : val(Array<Value>{Value(vals)...}) {}

  void set_unsorted_qns(const Array<QuantumNumberType>& vals);

  [[nodiscard]] String keys() const;

  [[nodiscard]] String values() const;

  [[nodiscard]] bool same_types_as(const LocalState& that) const;

  //! ouptut stream if all values
  friend std::ostream& operator<<(std::ostream& os, const LocalState& vl);

  //! input stream must have pre-set size
  friend std::istream& operator>>(std::istream& is, LocalState& vl);

  //! Test if there are bad quantum numbers (undefined ones)
  [[nodiscard]] bool good() const;
};

struct LevelTest {
  bool upp{true}, low{true};
};

//! A logical struct for global quantum numbers with species identifiers
struct GlobalState {
  static constexpr Index version = 1;  // Second version of quantum identifiers

  Index isotopologue_index{"Ar-8"_isot_index};
  ValueList val{};

  GlobalState() = default;

  explicit GlobalState(Index i, ValueList v = {})
      : isotopologue_index(i), val(std::move(v)) {}

  explicit GlobalState(const SpeciesIsotope& ir)
      : isotopologue_index(Species::find_species_index(ir)) {}

  explicit GlobalState(std::string_view s, Index v = version);

  [[nodiscard]] SpeciesIsotope Isotopologue() const noexcept;
  [[nodiscard]] SpeciesEnum Species() const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const GlobalState& gs);

  friend std::istream& operator>>(std::istream& is, GlobalState& gs);

  [[nodiscard]] GlobalState LowerLevel() const;
  [[nodiscard]] GlobalState UpperLevel() const;

  std::strong_ordering operator<=>(const GlobalState& g) const;
  bool operator==(const GlobalState& g) const;
  bool operator!=(const GlobalState& g) const;

  //! Checks whether all of the LHS is part of RHS
  [[nodiscard]] bool part_of(const GlobalState& other) const;

  //! Checks whether all of the RHS is part of LHS (uses reverse part_of)
  [[nodiscard]] bool may_be(const GlobalState& other) const;

  //! Checks whether all of the LHS is part of any of the RHS
  [[nodiscard]] LevelTest part_of(const GlobalState& g,
                                  const LocalState& l) const;

  //! Test if there are bad quantum numbers (undefined ones) or if the isotopologue is not a normal target
  [[nodiscard]] bool good() const;
};

//! StateMatchType operates so that a check less than a level should be 'better', bar None
enum class StateMatchType : char { Full, Level, Isotopologue, Species, None };

//! StateMatch, where you need to check the upp and low values when manipulating energy levels
struct StateMatch {
  StateMatchType type{StateMatchType::None};
  bool upp{false}, low{false};

  constexpr StateMatch() = default;

  StateMatch(const GlobalState& target,
             const LocalState& local,
             const GlobalState& global);

  StateMatch(const GlobalState& target, const GlobalState& key);

  //! It is of the desired type if it is less than the value, bar None
   bool operator==(StateMatchType x) const noexcept ;

   bool operator!=(StateMatchType x) const noexcept ;
};

//! VAMDC classes of quantum number cases
enum class VAMDC : char {
  asymcs,  // Schema for specifying the quantum numbers of closed-shell asymmetric top molecules
  asymos,  // Schema for specifying the quantum numbers of open-shell asymmetric top molecules
  dcs,  // Schema for specifying the quantum numbers of closed-shell, diatomic molecules
  hunda,  // Schema for specifying the quantum numbers for Hund's case (a) diatomic molecules
  hundb,  // Schema for specifying the quantum numbers for Hund's case (b) diatomic molecules
  lpcs,  // Schema for specifying the quantum numbers of closed-shell linear polyatomic molecules
  lpos,  // Schema for specifying the quantum numbers of open-shell linear polyatomic molecules
  ltcs,  // Schema for specifying the quantum numbers of closed-shell linear triatomic molecules
  ltos,  // Schema for specifying the quantum numbers of open-shell linear triatomic molecules
  nltcs,  // Schema for specifying the quantum numbers of closed-shell non-linear triatomic molecules
  nltos,  // Schema for specifying the quantum numbers of open-shell non-linear triatomic molecules
  sphcs,  // Schema for specifying the quantum numbers of closed-shell spherical top molecules
  sphos,  // Schema for specifying the quantum numbers of closed-shell spherical top molecules
  stcs  // Schema for specifying the quantum numbers of closed-shell, symmetric top molecules
};

/** Checks if a ValueList can belong to a given VAMDC type by
 * ensuring it cannot have some quantum numbers
 * 
 * @param l A list of values
 * @param type A type of VAMDC molecular model
 * @return true If it can belong to the VAMDC type
 * @return false If it cannot belong to the VAMDC type
 */
bool vamdcCheck(const ValueList& l, VAMDC type) ARTS_NOEXCEPT;

//! A default state of global quantum numbers
[[maybe_unused]] inline constexpr std::array global_types{
    QuantumNumberType::alpha,
    QuantumNumberType::config,
    QuantumNumberType::ElecStateLabel,
    QuantumNumberType::L,
    QuantumNumberType::Lambda,
    QuantumNumberType::Omega,
    QuantumNumberType::S,
    QuantumNumberType::Sigma,
    QuantumNumberType::SpinComponentLabel,
    QuantumNumberType::asSym,
    QuantumNumberType::elecInv,
    QuantumNumberType::elecRefl,
    QuantumNumberType::elecSym,
    QuantumNumberType::kronigParity,
    QuantumNumberType::l,
    QuantumNumberType::l1,
    QuantumNumberType::l10,
    QuantumNumberType::l11,
    QuantumNumberType::l12,
    QuantumNumberType::l2,
    QuantumNumberType::l3,
    QuantumNumberType::l4,
    QuantumNumberType::l5,
    QuantumNumberType::l6,
    QuantumNumberType::l7,
    QuantumNumberType::l8,
    QuantumNumberType::l9,
    QuantumNumberType::n,
    QuantumNumberType::parity,
    QuantumNumberType::r,
    QuantumNumberType::rotSym,
    QuantumNumberType::rovibSym,
    QuantumNumberType::sym,
    QuantumNumberType::tau,
    QuantumNumberType::term,
    QuantumNumberType::v,
    QuantumNumberType::v1,
    QuantumNumberType::v10,
    QuantumNumberType::v11,
    QuantumNumberType::v12,
    QuantumNumberType::v2,
    QuantumNumberType::v3,
    QuantumNumberType::v4,
    QuantumNumberType::v5,
    QuantumNumberType::v6,
    QuantumNumberType::v7,
    QuantumNumberType::v8,
    QuantumNumberType::v9,
    QuantumNumberType::vibInv,
    QuantumNumberType::vibRefl,
    QuantumNumberType::vibSym};

//! A default state of local quantum numbers
[[maybe_unused]] inline constexpr std::array local_types{QuantumNumberType::F,
                                                         QuantumNumberType::F1,
                                                         QuantumNumberType::F10,
                                                         QuantumNumberType::F11,
                                                         QuantumNumberType::F12,
                                                         QuantumNumberType::F2,
                                                         QuantumNumberType::F3,
                                                         QuantumNumberType::F4,
                                                         QuantumNumberType::F5,
                                                         QuantumNumberType::F6,
                                                         QuantumNumberType::F7,
                                                         QuantumNumberType::F8,
                                                         QuantumNumberType::F9,
                                                         QuantumNumberType::I,
                                                         QuantumNumberType::J,
                                                         QuantumNumberType::K,
                                                         QuantumNumberType::Ka,
                                                         QuantumNumberType::Kc,
                                                         QuantumNumberType::N};

std::ostream& operator<<(std::ostream& os, const Array<GlobalState>& a);
}  // namespace Quantum::Number

using QuantumNumberValue       = Quantum::Number::Value;
using QuantumNumberValueList   = Quantum::Number::ValueList;
using QuantumNumberLocalState  = Quantum::Number::LocalState;
using QuantumIdentifier        = Quantum::Number::GlobalState;
using ArrayOfQuantumIdentifier = Array<QuantumIdentifier>;

std::ostream& operator<<(std::ostream& os, const Array<QuantumNumberType>& a);

namespace std {
//! Allow Quantum::Number::TwoLevelValueHolder to be used in hashes
template <>
struct hash<Quantum::Number::TwoLevelValueHolder> {
  std::size_t operator()(const Quantum::Number::TwoLevelValueHolder& g) const {
    return std::hash<std::string_view>{}(g.upp.s.val()) ^
           (std::hash<std::string_view>{}(g.low.s.val()) << 1);
  }
};

//! Allow QuantumNumberValueList to be used in hashes
template <>
struct hash<QuantumNumberValueList> {
  std::size_t operator()(const QuantumNumberValueList& g) const {
    std::size_t out = 0;
    std::size_t i   = 1;
    for (auto& x : g) {
      out ^= (std::hash<QuantumNumberType>{}(x.type) ^
              (std::hash<Quantum::Number::TwoLevelValueHolder>{}(x.qn) << 1))
             << i;
      i++;
    }
    return out;
  }
};

//! Allow QuantumNumberLocalState to be used in hashes
struct LocalStateHash {
  std::size_t operator()(const QuantumNumberLocalState& g) const {
    return std::hash<QuantumNumberValueList>{}(g.val);
  }
};

//! Allow QuantumIdentifier to be used in hashes
template <>
struct hash<QuantumIdentifier> {
  std::size_t operator()(const QuantumIdentifier& g) const {
    return static_cast<std::size_t>(g.isotopologue_index) ^
           (std::hash<QuantumNumberValueList>{}(g.val) << 1);
  }
};
}  // namespace std

template <>
struct std::formatter<QuantumNumberValue> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumNumberValue& v,
                              FmtContext& ctx) const {
    const auto quote = tags.quote();
    const auto sep   = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx,
                v.type,
                sep,
                quote,
                v.str_upp(),
                quote,
                sep,
                quote,
                v.str_low(),
                quote);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<QuantumNumberValueList> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumNumberValueList& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.values);
  }
};

template <>
struct std::formatter<QuantumNumberLocalState> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumNumberLocalState& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.val);
  }
};

template <>
struct std::formatter<QuantumIdentifier> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const QuantumIdentifier& v,
                              FmtContext& ctx) const {
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.Isotopologue().FullName(), tags.sep(), v.val);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<Quantum::Number::QuantumNumberValueType> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Quantum::Number::QuantumNumberValueType& v,
                              FmtContext& ctx) const {
    const auto quote = tags.quote();

    switch (v) {
      case Quantum::Number::QuantumNumberValueType::I:
        return format_to(ctx.out(), "{}I{}", quote, quote);
      case Quantum::Number::QuantumNumberValueType::S:
        return format_to(ctx.out(), "{}S{}", quote, quote);
      case Quantum::Number::QuantumNumberValueType::H:
        return format_to(ctx.out(), "{}H{}", quote, quote);
    }

    return format_to(ctx.out(), "{}bad-value{}", quote, quote);
  }
};
