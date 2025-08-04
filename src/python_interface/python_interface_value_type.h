#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <mystring.h>

#include <memory>
namespace Python {
template <typename type>
struct ValueHolder {
  std::shared_ptr<type> val;

  template <typename T>
  using common_type = std::common_type_t<type, T>;

  ValueHolder(std::shared_ptr<type> x) : val(std::move(x)) {}
  ValueHolder() : val(new type{}) {}
  ValueHolder(const ValueHolder& other) : val(new type{*other.val}) {}
  ValueHolder& operator=(const ValueHolder& other) {
    *val = *other.val;
    return *this;
  }
  ValueHolder& operator=(std::shared_ptr<type> x) { val = std::move(x); }

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
}  // namespace Python

template <typename T>
struct std::formatter<Python::ValueHolder<T>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Python::ValueHolder<T>& v,
                              FmtContext& ctx) const {
    std::formatter<T> fmt{};

    if constexpr (std::same_as<T, String>) {
      if (tags.quoted) tags.format(ctx, '"');
    }

    fmt.format(*v.val, ctx);

    if constexpr (std::same_as<T, String>) {
      if (tags.quoted) tags.format(ctx, '"');
    }

    return ctx.out();
  }
};
