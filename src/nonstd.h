#ifndef nonstd_h
#define nonstd_h

/*! For std functions that are not constexpr */
namespace nonstd {
/*! abs(x) returns |x| using -x if x < 0 or x otherwise */
template <class T> constexpr T abs(T x) noexcept {return x < 0 ? - x : x;}
}  // nonstd

#endif
