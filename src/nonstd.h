#ifndef nonstd_h
#define nonstd_h

/*! For std functions that are changed somehow */
namespace nonstd {
/*! abs(x) returns |x| using -x if x < 0 or x otherwise
 * 
 * Reason to re-implement: std::abs(x) is not officially constexpr
 * 
 * @param[in] x Any real value type
 * @return |x|
 */
template <class T> constexpr T abs(T x) noexcept {return x < 0 ? - x : x;}

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
  return ch == '0' or ch == '1' or ch == '2' or ch == '3' or
         ch == '4' or ch == '5' or ch == '6' or ch == '7' or
         ch == '8' or ch == '9';
}

/** Returns 1 if x is a standard space-character
 * 
 * Reason to re-implement: std::isspace(ch) is not officially constexpr
 * 
 * The int-interface is kept from the standard
 * 
 * @param[in] x a character
 * @return true if x is a space
 */
constexpr bool isspace(unsigned char ch) noexcept {
  return ch == ' ' or ch == '\n' or ch == '\r' or ch == '\t' or ch == '\f' or
         ch == '\v';
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
constexpr bool isnan(double d) noexcept {return d not_eq d;}
constexpr bool isnan(long double d) noexcept {return d not_eq d;}
constexpr bool isnan(float d) noexcept {return d not_eq d;}
} // namespace nonstd

#endif
