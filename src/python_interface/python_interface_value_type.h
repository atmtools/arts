#ifndef python_interface_value_type_h
#define python_interface_value_type_h

#include <cmath>
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
  
  constexpr ValueType operator+(ValueHolder x) const noexcept {return val + x.val;}
  constexpr ValueType operator-(ValueHolder x) const noexcept {return val - x.val;}
  constexpr ValueType operator*(ValueHolder x) const noexcept {return val * x.val;}
  constexpr ValueType operator/(ValueHolder x) const noexcept {return val / x.val;}
  
  constexpr ValueHolder& operator+=(ValueHolder x) noexcept {val += x.val; return *this; }
  constexpr ValueHolder& operator-=(ValueHolder x) noexcept {val -= x.val; return *this; }
  constexpr ValueHolder& operator*=(ValueHolder x) noexcept {val *= x.val; return *this; }
  constexpr ValueHolder& operator/=(ValueHolder x) noexcept {val /= x.val; return *this; }

  constexpr ValueType pow(ValueHolder x) const noexcept {using std::pow; return ValueType(pow(val, Numeric(x.val)));}
};

using Numeric_ = ValueHolder<Numeric>;
using Index_ = ValueHolder<Index>;
}

#endif  // python_interface_value_type_h
