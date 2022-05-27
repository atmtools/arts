#ifndef python_interface_value_type_h
#define python_interface_value_type_h

#include <cmath>
#include <type_traits>
#include <utility>

#include <matpack.h>

namespace Python {
template <typename ValueType>
struct ValueHolder {
  ValueType val;
  constexpr operator ValueType&() noexcept { return val; }
  constexpr operator const ValueType&() const noexcept { return val; }
  constexpr ValueHolder& operator=(ValueType x) noexcept {using std::move; val = move(x); return *this; }

  friend std::ostream& operator<<(std::ostream& os, ValueHolder x) {return os << x.val;}

  constexpr bool operator<(ValueHolder x) const noexcept  {return val <  x.val;}
  constexpr bool operator<=(ValueHolder x) const noexcept {return val <= x.val;}
  constexpr bool operator>(ValueHolder x) const noexcept  {return val >  x.val;}
  constexpr bool operator>=(ValueHolder x) const noexcept {return val >= x.val;}
  constexpr bool operator==(ValueHolder x) const noexcept {return val == x.val;}
  constexpr bool operator!=(ValueHolder x) const noexcept {return val != x.val;}

  constexpr ValueType operator-() const noexcept {return - val;}
  constexpr ValueType operator+() const noexcept {return val;}
  
  constexpr ValueType operator+(ValueHolder x) const noexcept {return val + x.val;}
  constexpr ValueType operator-(ValueHolder x) const noexcept {return val - x.val;}
  constexpr ValueType operator*(ValueHolder x) const noexcept {return val * x.val;}
  constexpr ValueType operator/(ValueHolder x) const noexcept {return val / x.val;}
  
  constexpr ValueHolder& operator+=(ValueHolder x) noexcept { return operator=(operator+(x)); }
  constexpr ValueHolder& operator-=(ValueHolder x) noexcept { return operator=(operator-(x)); }
  constexpr ValueHolder& operator*=(ValueHolder x) noexcept { return operator=(operator*(x)); }
  constexpr ValueHolder& operator/=(ValueHolder x) noexcept { return operator=(operator/(x)); }
};

using Numeric_ = ValueHolder<Numeric>;
using Index_ = ValueHolder<Index>;

//! Return a possibly dangling reference to a value, be sure it is not dangling!
template <typename T> auto& as_ref(T& x) noexcept {
  static_assert(not std::is_const_v<T>, "Cannot return a constant");

  if constexpr (std::is_same_v<T, Numeric>) return reinterpret_cast<Numeric_&>(x);
  else if constexpr (std::is_same_v<T, Index>) return reinterpret_cast<Index_&>(x);
  else return x;
}
} // namespace Python

#endif  // python_interface_value_type_h
