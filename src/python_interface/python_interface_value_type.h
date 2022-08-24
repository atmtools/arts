#ifndef python_interface_value_type_h
#define python_interface_value_type_h

#include <concepts>
#include <cmath>
#include <type_traits>
#include <utility>

#include <matpack.h>

namespace Python {
template <typename T>
concept Arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept PythonValueHolderType = requires(T a) {
  { a * a } -> Arithmetic;
  { 1 * a } -> Arithmetic;
  { 1. * a } -> Arithmetic;
  { a } -> std::totally_ordered;
  { a.val } -> std::totally_ordered;
};

template <typename type>
struct ValueHolder {
  type val;

  template <typename T> using common_type = std::common_type_t<type, T>;

  constexpr operator type&() noexcept { return val; }
  constexpr operator const type&() const noexcept { return val; }
  constexpr ValueHolder& operator=(type x) noexcept { val = x; return *this; }

  friend std::ostream& operator<<(std::ostream& os, ValueHolder x) {return os << x.val;}

  constexpr type operator-() const noexcept {return - val;}
  constexpr type operator+() const noexcept {return val;}

  template <typename T> friend constexpr auto operator<=>(ValueHolder a, PythonValueHolderType auto b) noexcept {return static_cast<common_type<T>>(a.val) <=> static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator+(ValueHolder a, PythonValueHolderType auto b) noexcept {return static_cast<common_type<T>>(a.val) + static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator-(ValueHolder a, PythonValueHolderType auto b) noexcept {return static_cast<common_type<T>>(a.val) - static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator*(ValueHolder a, PythonValueHolderType auto b) noexcept {return static_cast<common_type<T>>(a.val) * static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator/(ValueHolder a, PythonValueHolderType auto b) noexcept {return static_cast<common_type<T>>(a.val) / static_cast<common_type<T>>(b.val);}

  template <typename T> friend constexpr auto operator<=>(type a, ValueHolder b) noexcept {return a <=> b.val;}
  template <typename T> friend constexpr auto operator+(type a, ValueHolder b) noexcept {return a + b.val;}
  template <typename T> friend constexpr auto operator-(type a, ValueHolder b) noexcept {return a - b.val;}
  template <typename T> friend constexpr auto operator*(type a, ValueHolder b) noexcept {return a * b.val;}
  template <typename T> friend constexpr auto operator/(type a, ValueHolder b) noexcept {return a / b.val;}
  
  constexpr ValueHolder& operator+=(auto x) noexcept { return operator=(val + static_cast<type>(x)); }
  constexpr ValueHolder& operator-=(auto x) noexcept { return operator=(val - static_cast<type>(x)); }
  constexpr ValueHolder& operator*=(auto x) noexcept { return operator=(val * static_cast<type>(x)); }
  constexpr ValueHolder& operator/=(auto x) noexcept { return operator=(val / static_cast<type>(x)); }

  static_assert(PythonValueHolderType<ValueHolder>);
};

// Set the type and ensure they are correct
using Numeric_ = ValueHolder<Numeric>;
using Index_ = ValueHolder<Index>;

template <typename T> concept isNumeric = std::is_same_v<Numeric, std::remove_cvref_t<T>>;
template <typename T> concept isIndex = std::is_same_v<Numeric, std::remove_cvref_t<T>>;

//! Return a possibly dangling reference to a value, be sure it is not dangling!
template <typename T> auto& as_ref(T& x) noexcept {
  static_assert(not std::is_const_v<T>, "Cannot return a constant");

  if constexpr (isNumeric<T>) return reinterpret_cast<Numeric_&>(x);
  else if constexpr (isIndex<T>) return reinterpret_cast<Index_&>(x);
  else return x;
}
} // namespace Python

#endif  // python_interface_value_type_h
