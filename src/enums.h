#ifndef enums_h
#define enums_h

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>

#include "debug.h"
#include "nonstd.h"

/** Checks if the enum number is good.
 * 
 * It is good as long as it is above-equal to 0 and less than FINAL
 * 
 * @param[in] x An enum class defined by ENUMCLASS or similar
 * @return Whether the enum is good
 */
template <typename EnumType>
constexpr bool good_enum(EnumType x) noexcept {
  return std::size_t(x) < std::size_t(EnumType::FINAL);
}

/** Internal string view array generator.
 * 
 * Automagically seperates a list by the commas in it,
 * removing spaces before creating the named variable
 * 
 * @param[in] strchars A list of characters
 * @return An array of views into the original strchars
 */
template <typename EnumType> constexpr
std::array<std::string_view, size_t(EnumType::FINAL)> enum_strarray(
  const std::string_view strchars) noexcept {
  std::array<std::string_view, size_t(EnumType::FINAL)> out;
  
  // Find the start
  std::string_view::size_type N0 = 0;
  
  // Set all the values
  for (auto& str: out) {
    // Find a comma but never look beyond the end of the string
    const std::string_view::size_type N1 =
      std::min(strchars.find(',', N0), strchars.size());
    
    // Set the string between start and the length of the string
    str = strchars.substr(N0, N1 - N0);
    
    // Remove spaces at the beginning and at the end of the string
    while (nonstd::isspace(str.front())) str.remove_prefix(1);
    while (nonstd::isspace(str.back())) str.remove_suffix(1);
    
    // Set the new start for the next iteration
    N0 = N1 + 1;
  }
  
  return out;
}

constexpr std::array<std::string_view, 0> enum_strarray() noexcept { return {}; }

/** A list of all enum types by index-conversion
 * 
 * Note that this assumes the enums are sorted from 0-FINAL
 * 
 * @return A list of all enum types
 */
template <typename EnumType> constexpr
std::array<EnumType, size_t(EnumType::FINAL)> enum_typarray() noexcept {
  std::array<EnumType, size_t(EnumType::FINAL)> out{};
  for (size_t i = 0; i < size_t(EnumType::FINAL); i++)
    out[i] = EnumType(i);
  return out;
}

/** Checks if the enum class type is good and otherwise throws
 * an error message composed by variadic input
 * 
 * @param[in] type The enum class type variable
 * @param[in] ...args A list of errors to print if type is bad
 */
template <typename EnumType, typename ... Messages> constexpr
void check_enum_error(EnumType type, Messages ... args) {
  ARTS_USER_ERROR_IF (not good_enum(type), args...)
}

/*! Enum style
 *
 * Generates a "enum class ENUMTYPE : TYPE"
 * with all the following arguments and terminated by FINAL
 *
 * Additionally, will fill a local namespace "enumstrs" with
 * a std::array<std::string_view, TYPE(ENUMTYPE::FINAL)> with all of
 * the names in the ENUMTYPE enum class (bar FINAL)
 *
 * Additionally, will fill a local namespace "enumtyps" with
 * a std::array<ENUMTYPE, TYPE(ENUMTYPE::FINAL)> with all of
 * the types in the ENUMTYPE enum class (bar FINAL)
 *
 * Additionally, will generate a constexpr "toString" function
 * that takes a ENUMTYPE object and returns either its partial
 * name or "BAD ENUMTYPE" as a string_view
 *
 * Additionally, will generate a constexpr "toENUMTYPE" function
 * that takes a std::string_view and returns a corresponding ENUMTYPE
 * object or ENUMTYPE::FINAL if the object is bad
 * 
 * Additionally, will generate intuitive std::istream& and
 * std::ostream& operators
 *
 * Use the "good_enum(ENUMTYPE)" template function to check
 * if the enum class object is any good
 * 
 * Use check_enum_error(ENUMTYPE, message...) to throw an error
 * composed by the variadic message if the enum value is bad
 *
 * \verbatim
 // Example:
 ENUMCLASS(Test, char, Value)
 
 // Generates (effectively):
 enum class Test : char {Value, FINAL};
 namespace enumstrs {
    constexpr std::array<std::string_view, 1> TestNames={"Value"};
 }
 namespace enumtyps {
    constexpr std::array<Test, 1> TestTypes={Test::Value};
 }
 constexpr std::string_view toString(Test x) noexcept;
 constexpr Test toTest(const std::string_view& x) noexcept;
 inline std::ostream &operator<<(std::ostream &os, const Test x);
 inline std::istream &operator>>(std::istream &is, Test &x);  // throws if not good_enum(x) at end
 \endverbatim
 */
#define ENUMCLASS(ENUMTYPE, TYPE, ...)                                    \
  enum class ENUMTYPE : TYPE { __VA_ARGS__, FINAL };                      \
                                                                          \
  namespace enumstrs {                                                    \
  constexpr auto ENUMTYPE##Names = enum_strarray<ENUMTYPE>(#__VA_ARGS__); \
  }                                                                       \
                                                                          \
  namespace enumtyps {                                                    \
  [[maybe_unused]]                                                        \
  constexpr auto ENUMTYPE##Types = enum_typarray<ENUMTYPE>();             \
  }                                                                       \
                                                                          \
  constexpr std::string_view toString(ENUMTYPE x) noexcept {              \
    if (good_enum(x))                                                     \
      return enumstrs::ENUMTYPE##Names[(TYPE)x];                          \
    return "BAD " #ENUMTYPE;                                              \
  }                                                                       \
                                                                          \
  constexpr ENUMTYPE to##ENUMTYPE(const std::string_view x) noexcept {    \
    for (TYPE i = 0; i < (TYPE)ENUMTYPE::FINAL; i++)                      \
      if (enumstrs::ENUMTYPE##Names[i] == x) return ENUMTYPE(i);          \
    return ENUMTYPE::FINAL;                                               \
  }                                                                       \
                                                                          \
  constexpr ENUMTYPE to##ENUMTYPE##OrThrow(const std::string_view x) {    \
    const ENUMTYPE out = to##ENUMTYPE(x);                                 \
    check_enum_error(out, "Cannot understand argument: \"", x, "\"\n"     \
                     "Valid " #ENUMTYPE " options are: ["                 \
                     #__VA_ARGS__ "]");                                   \
    return out;                                                           \
  }                                                                       \
                                                                          \
  inline std::ostream &operator<<(std::ostream &os, const ENUMTYPE x) {   \
    return os << toString(x);                                             \
  }                                                                       \
                                                                          \
  inline std::istream &operator>>(std::istream &is, ENUMTYPE &x) {        \
    std::string val;                                                      \
    is >> val;                                                            \
    x = to##ENUMTYPE##OrThrow(val);                                       \
    return is;                                                            \
  }

#define ENUMCLASS_EMPTY(ENUMTYPE, TYPE)                                        \
  enum class ENUMTYPE : TYPE { FINAL };                                        \
                                                                               \
  namespace enumstrs {                                                         \
  constexpr auto ENUMTYPE##Names = enum_strarray();                            \
  }                                                                            \
                                                                               \
  namespace enumtyps {                                                         \
  [[maybe_unused]] constexpr auto ENUMTYPE##Types = enum_typarray<ENUMTYPE>(); \
  }                                                                            \
                                                                               \
  constexpr std::string_view toString(ENUMTYPE x) noexcept {                   \
    if (good_enum(x)) return enumstrs::ENUMTYPE##Names[(TYPE)x];               \
    return "BAD " #ENUMTYPE;                                                   \
  }                                                                            \
                                                                               \
  constexpr ENUMTYPE to##ENUMTYPE(const std::string_view x) noexcept {         \
    for (TYPE i = 0; i < (TYPE)ENUMTYPE::FINAL; i++)                           \
      if (enumstrs::ENUMTYPE##Names[i] == x) return ENUMTYPE(i);               \
    return ENUMTYPE::FINAL;                                                    \
  }                                                                            \
                                                                               \
  constexpr ENUMTYPE to##ENUMTYPE##OrThrow(const std::string_view x) {         \
    ENUMTYPE out = to##ENUMTYPE(x);                                            \
    check_enum_error(out,                                                      \
                     "Cannot understand argument: \"",                         \
                     x,                                                        \
                     "\"\n"                                                    \
                     "Valid " #ENUMTYPE " options are: []");                   \
    return out;                                                                \
  }                                                                            \
                                                                               \
  inline std::ostream &operator<<(std::ostream &os, const ENUMTYPE x) {        \
    return os << toString(x);                                                  \
  }                                                                            \
                                                                               \
  inline std::istream &operator>>(std::istream &is, ENUMTYPE &x) {             \
    std::string val;                                                           \
    is >> val;                                                                 \
    x = to##ENUMTYPE##OrThrow(val);                                            \
    return is;                                                                 \
  }

#endif  // enums_h

