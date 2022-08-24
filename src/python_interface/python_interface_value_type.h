#ifndef python_interface_value_type_h
#define python_interface_value_type_h

#include <concepts>
#include <cmath>
#include <type_traits>
#include <utility>

#include <matpack.h>

namespace Python {
template <typename type>
struct ValueHolder {
  type val;

  template <typename T> using common_type = std::common_type_t<type, T>;

  constexpr ValueHolder() = default;
  constexpr ValueHolder(const ValueHolder&) = default;
  constexpr ValueHolder(ValueHolder&&) = default;
  constexpr ValueHolder& operator=(const ValueHolder&) = default;
  constexpr ValueHolder& operator=(ValueHolder&&) = default;
  constexpr ValueHolder(const type& a) noexcept : val(a) {}
  constexpr ValueHolder(type&& a) noexcept : val(a) {}
  constexpr ValueHolder& operator=(const type& a) noexcept { val = a; return *this; }
  constexpr ValueHolder& operator=(type&& a) noexcept { val = a; return *this; }

  constexpr operator type&() noexcept { return val; }
  constexpr operator const type&() const noexcept { return val; }

  friend std::ostream& operator<<(std::ostream& os, ValueHolder x) {return os << x.val;}

  constexpr ValueHolder operator-() const noexcept {return - val;}
  constexpr ValueHolder operator+() const noexcept {return val;}

  template <typename T> friend constexpr auto operator<=>(ValueHolder a, ValueHolder<T> b) noexcept {return static_cast<common_type<T>>(a.val) <=> static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator+(ValueHolder a, ValueHolder<T> b) noexcept {return static_cast<common_type<T>>(a.val) + static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator-(ValueHolder a, ValueHolder<T> b) noexcept {return static_cast<common_type<T>>(a.val) - static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator*(ValueHolder a, ValueHolder<T> b) noexcept {return static_cast<common_type<T>>(a.val) * static_cast<common_type<T>>(b.val);}
  template <typename T> friend constexpr auto operator/(ValueHolder a, ValueHolder<T> b) noexcept {return static_cast<common_type<T>>(a.val) / static_cast<common_type<T>>(b.val);}

  friend constexpr auto operator<=>(std::floating_point auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) <=> static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator+(std::floating_point auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) + static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator-(std::floating_point auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) - static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator*(std::floating_point auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) * static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator/(std::floating_point auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) / static_cast<common_type<decltype(a)>>(b.val);}

  friend constexpr auto operator<=>(std::integral auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) <=> static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator+(std::integral auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) + static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator-(std::integral auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) - static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator*(std::integral auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) * static_cast<common_type<decltype(a)>>(b.val);}
  friend constexpr auto operator/(std::integral auto a, ValueHolder b) noexcept {return static_cast<common_type<decltype(a)>>(a) / static_cast<common_type<decltype(a)>>(b.val);}
  
  constexpr ValueHolder& operator+=(auto x) noexcept { return operator=(val + static_cast<type>(x)); }
  constexpr ValueHolder& operator-=(auto x) noexcept { return operator=(val - static_cast<type>(x)); }
  constexpr ValueHolder& operator*=(auto x) noexcept { return operator=(val * static_cast<type>(x)); }
  constexpr ValueHolder& operator/=(auto x) noexcept { return operator=(val / static_cast<type>(x)); }
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
