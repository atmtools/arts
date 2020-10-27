#ifndef enums_h
#define enums_h

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <string>

template <typename T>
constexpr bool good_enum(T x) {
  return long(x) < long(T::FINAL) and long(x) >= 0;
}
template <typename T>
std::array<std::string, long(T::FINAL)> enum_strarray(
    const std::string &strchars) {
  std::array<std::string, long(T::FINAL)> out;
  std::istringstream x(strchars);
  for (long i = 0; i < long(T::FINAL); i++) {
    std::getline(x, out[i], ',');
    out[i].erase(std::remove(out[i].begin(), out[i].end(), ' '), out[i].end());
  }
  return out;
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
#define ENUMCLASS(ENUMTYPE, TYPE, ...)                                         \
  enum class ENUMTYPE : TYPE { __VA_ARGS__, FINAL };                           \
                                                                               \
  namespace enumstrs {                                                         \
  static const auto ENUMTYPE##Names = enum_strarray<ENUMTYPE>(#__VA_ARGS__);   \
  }                                                                            \
                                                                               \
  inline std::string toString(ENUMTYPE x) noexcept {                           \
    if (good_enum(x))                                                          \
      return enumstrs::ENUMTYPE##Names[(TYPE)x];                               \
    else                                                                       \
      return "BAD " #ENUMTYPE;                                                 \
  }                                                                            \
                                                                               \
  inline ENUMTYPE to##ENUMTYPE(const std::string &x) noexcept {                \
    for (TYPE i = 0; i < (TYPE)ENUMTYPE::FINAL; i++)                           \
      if (enumstrs::ENUMTYPE##Names[i] == x) return ENUMTYPE(i);               \
    return ENUMTYPE::FINAL;                                                    \
  }                                                                            \
                                                                               \
  inline std::ostream &operator<<(std::ostream &os, const ENUMTYPE x) {        \
    return os << toString(x);                                                  \
  }                                                                            \
                                                                               \
  inline std::istream &operator>>(std::istream &is, ENUMTYPE &x) {             \
    std::string val;                                                           \
    is >> val;                                                                 \
    x = to##ENUMTYPE(val);                                                     \
    if (x == ENUMTYPE::FINAL) throw std::runtime_error("Bad read " #ENUMTYPE); \
    return is;                                                                 \
  }

#endif  // enums_h

