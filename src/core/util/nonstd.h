#ifndef nonstd_h
#define nonstd_h

#include <cmath>
#include <utility>

/*! For std functions that are changed somehow */
namespace nonstd {
/*! abs(x) returns |x| using -x if x < 0 or x otherwise
 * 
 * Reason to re-implement: std::abs(x) is not officially constexpr
 * 
 * @param[in] x Any real value type
 * @return |x|
 */
template <class T>
constexpr T abs(T x) noexcept {
  return x < 0 ? -x : x;
}

/*! pow(x, v) returns x^v not using std::pow but using exp(v * log(x))
 * 
 * Reason to re-implement: factor 4 faster
 * 
 * @param[in] x Any real value type
 * @return x^v
 */
template <typename T, typename U>
constexpr auto pow(T&& v, U&& x) {
  return std::exp(std::forward<T>(x) * std::log(std::forward<U>(v)));
}

/*! Checks if the given character in 0123456789.
 * 
 * Reason to re-implement: std::isdigit(ch) is not officially constexpr
 * 
 * The int-interface is kept from the standard
 * 
 * @param[in] ch A character
 * @return int from a simple boolean.
 */
constexpr int isdigit(int ch) noexcept {
  return ch == '0' or ch == '1' or ch == '2' or ch == '3' or ch == '4' or
         ch == '5' or ch == '6' or ch == '7' or ch == '8' or ch == '9';
}

/** Returns true if x is a standard space-character
 * 
 * Reason to re-implement: std::isspace(ch) is not officially constexpr
 * 
 * @param[in] x a character
 * @return true if x is a space
 */
constexpr bool isspace(unsigned char ch) noexcept {
  return ch == ' ' or ch == '\n' or ch == '\r' or ch == '\t' or ch == '\f' or
         ch == '\v';
}

/** Returns true if x is a standard abc/ABC-character
 * 
 * @param[in] x a character
 * @return true if x is a space
 */
constexpr bool isabc(unsigned char ch) noexcept {
  return ch == 'a' or ch == 'b' or ch == 'c' or ch == 'd' or ch == 'e' or
         ch == 'f' or ch == 'g' or ch == 'h' or ch == 'i' or ch == 'j' or
         ch == 'k' or ch == 'l' or ch == 'm' or ch == 'n' or ch == 'o' or
         ch == 'p' or ch == 'q' or ch == 'r' or ch == 's' or ch == 't' or
         ch == 'u' or ch == 'v' or ch == 'w' or ch == 'x' or ch == 'y' or
         ch == 'z' or ch == 'A' or ch == 'B' or ch == 'C' or ch == 'D' or
         ch == 'E' or ch == 'F' or ch == 'G' or ch == 'H' or ch == 'I' or
         ch == 'J' or ch == 'K' or ch == 'L' or ch == 'M' or ch == 'N' or
         ch == 'O' or ch == 'P' or ch == 'Q' or ch == 'R' or ch == 'S' or
         ch == 'T' or ch == 'U' or ch == 'V' or ch == 'W' or ch == 'X' or
         ch == 'Y' or ch == 'Z';
}

/** Returns true if x is a standard bracket-character
 * 
 * @param[in] x a character
 * @return true if x is a space
 */
constexpr bool isbracket(unsigned char c) noexcept {
  return c == '(' or c == ')' or c == '[' or c == ']' or c == '{' or c == '}' or
         c == '<' or c == '>';
}

/*! Checks if the given value is nan
 * 
 * Reason to re-implement: std::isnan(d) is not officially constexpr
 * 
 * Only the three basic interfaces are implemented
 * 
 * @param[in] d A value
 * @return int from a simple boolean.
 */
constexpr bool isnan(double d) noexcept { return d not_eq d; }
constexpr bool isnan(long double d) noexcept { return d not_eq d; }
constexpr bool isnan(float d) noexcept { return d not_eq d; }
}  // namespace nonstd

#endif
