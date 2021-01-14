#ifndef enums_h
#define enums_h

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>

template <typename T>
constexpr bool good_enum(T x) {
  return long(x) < long(T::FINAL) and long(x) >= 0;
}

constexpr bool is_space_char(char x) {
  return x == ' '  or
         x == '\n' or
         x == '\r' or
         x == '\t' or
         x == '\f' or
         x == '\v';
}

template <typename T> constexpr
std::array<std::string_view, size_t(T::FINAL)> enum_strarray(
  const std::string_view strchars) {
  std::array<std::string_view, size_t(T::FINAL)> out;
  
  // Find the start
  std::string_view::size_type N0 = 0;
  
  // Set all the values
  for (auto& str: out) {
    // Find a comma but never look beyond the end of the string
    const std::string_view::size_type N1 = std::min(strchars.find(',', N0), strchars.size());
    
    // Set the string between start and the length of the string
    str = strchars.substr(N0, N1 - N0);
    
    // Remove spaces at the beginning and at the end of the string
    while (is_space_char(str.front())) str.remove_prefix(1);
    while (is_space_char(str.back())) str.remove_suffix(1);
    
    // Set the new start for the next iteration
    N0 = N1 + 1;
  }
  
  return out;
}

template <typename EnumType, typename ... Messages>
void EnumErrorQuery(EnumType type, Messages ... args) {
  if (not good_enum(type)) {
    std::ostringstream os;
    (os << ... << args);
    throw std::runtime_error(os.str());
  }
}

/* Enum style
 *
 * Generates a "enum class ENUMTYPE : long"
 * with all the following arguments and terminated by FINAL
 *
 * Additionally, will fill a local namespace "enumstrs" with
 * a std::array<std::string, long(ENUMTYPE::FINAL)> with all of
 * the names in the ENUMTYPE enum class
 *
 * Additionally, will generate a inlined "toString" function
 * that takes a ENUMTYPE object and returns either its partial
 * name or "BAD ENUMTYPE" as a string
 *
 * Additionally, will generate a inlined "toENUMTYPE" function
 * that takes a std::string and returns a corresponding ENUMTYPE
 * object or ENUMTYPE::FINAL if the object is bad
 *
 * Use the "good_enum(ENUMTYPE)" template function to check
 * if the enum classs object is any good
 *
 * Will be updated as soon as possible to ensure that all functions
 * that can be are turned into constexpr functions
 *
 * Example:
 * ENUMCLASS(Test, Value)
 *
 * Generates:
 * enum class Test : long {Value, FINAL};
 * namespace enumstrs {std::array<std::string, 1> TestNames={"Value"};};
 * std::string toString(Test x) noexcept;
 * Test toTest(const std::string& x) noexcept;
 */
#define ENUMCLASS(ENUMTYPE, TYPE, ...)                                    \
  enum class ENUMTYPE : TYPE { __VA_ARGS__, FINAL };                      \
                                                                          \
  namespace enumstrs {                                                    \
  constexpr auto ENUMTYPE##Names = enum_strarray<ENUMTYPE>(#__VA_ARGS__); \
  }                                                                       \
                                                                          \
  constexpr std::string_view toString(ENUMTYPE x) noexcept {              \
    if (good_enum(x))                                                     \
      return enumstrs::ENUMTYPE##Names[(TYPE)x];                          \
    else                                                                  \
      return "BAD " #ENUMTYPE;                                            \
  }                                                                       \
                                                                          \
  constexpr ENUMTYPE to##ENUMTYPE(const std::string_view x) noexcept {    \
    for (TYPE i = 0; i < (TYPE)ENUMTYPE::FINAL; i++)                      \
      if (enumstrs::ENUMTYPE##Names[i] == x) return ENUMTYPE(i);          \
    return ENUMTYPE::FINAL;                                               \
  }                                                                       \
                                                                          \
  inline std::ostream &operator<<(std::ostream &os, const ENUMTYPE x) {   \
    return os << toString(x);                                             \
  }                                                                       \
                                                                          \
  inline std::istream &operator>>(std::istream &is, ENUMTYPE &x) {        \
    std::string val;                                                      \
    is >> val;                                                            \
    x = to##ENUMTYPE(val);                                                \
    EnumErrorQuery(x, "Cannot understand value: ", val);                  \
    return is;                                                            \
  }

#endif  // enums_h

