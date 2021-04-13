#ifndef nonstd_h
#define nonstd_h

/*! For std functions that are changed somehow */
namespace nonstd {
/*! abs(x) returns |x| using -x if x < 0 or x otherwise
 * 
 * Reason to re-implement:  std::abs(x) is not constexpr
 * 
 * @param[in] x Any real value type
 * @return |x|
 */
template <class T> constexpr T abs(T x) noexcept {return x < 0 ? - x : x;}
}  // nonstd

#endif
