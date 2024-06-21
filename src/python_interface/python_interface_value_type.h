#ifndef python_interface_value_type_h
#define python_interface_value_type_h

#include <matpack_concepts.h>
namespace Python {
template <typename type>
struct ValueHolder {
  type* val;

  template <typename T>
  using common_type = std::common_type_t<type, T>;

  ValueHolder() : val(new type{}) {}
  ValueHolder(const ValueHolder& other) : val(new type{*other.val}) {}
  ValueHolder& operator=(const ValueHolder& other) {
    *val = *other.val;
    return *this;
  }

  ValueHolder(type t) : val(new type{t}) {}
  ValueHolder& operator=(type a) noexcept {
    *val = a;
    return *this;
  }

  operator type&() noexcept { return *val; }
  operator const type&() const noexcept { return *val; }

  friend std::ostream& operator<<(std::ostream& os, const ValueHolder& a) {
    return os << *a.val;
  }
};

// Set the type and ensure they are correct
using Numeric_ = ValueHolder<Numeric>;
using Index_   = ValueHolder<Index>;
}  // namespace Python

#endif  // python_interface_value_type_h
