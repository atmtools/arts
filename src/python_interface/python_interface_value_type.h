#ifndef python_interface_value_type_h
#define python_interface_value_type_h

#include <matpack_concepts.h>

#include <cmath>
#include <concepts>
#include <memory>
#include <type_traits>
#include <utility>

namespace Python {
template <typename type>
struct ValueHolder {
  type* val;

  template <typename T>
  using common_type = std::common_type_t<type, T>;

  ValueHolder() : val(std::make_shared<type>()) {}
  ValueHolder(const ValueHolder& other) : val(new type{*other.val}) {}
  ValueHolder(ValueHolder&&) noexcept = default;
  ValueHolder& operator=(const ValueHolder& other) {
    *val = *other.val;
    return *this;
  }
  ValueHolder& operator=(ValueHolder&&) noexcept = default;
  ValueHolder(const type& a) : val(std::make_shared<type>(a)) {}
  ValueHolder(type&& a) noexcept : val(std::make_shared<type>(a)) {}
  ValueHolder& operator=(const type& a) noexcept {
    *val = a;
    return *this;
  }
  ValueHolder& operator=(type&& a) noexcept {
    *val = a;
    return *this;
  }

  operator type&() noexcept { return *val; }
  operator const type&() const noexcept { return *val; }
  operator std::shared_ptr<type>&() noexcept { return val; }
  operator const std::shared_ptr<type>&() const noexcept { return val; }

  friend std::ostream& operator<<(std::ostream& os, const ValueHolder& a) {
    return os << *a.val;
  }
};

// Set the type and ensure they are correct
using Numeric_ = ValueHolder<Numeric>;
using Index_ = ValueHolder<Index>;
}  // namespace Python

#endif  // python_interface_value_type_h
