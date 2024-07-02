#pragma once

#include <matpack.h>

#include <functional>
#include <ostream>

#include "debug.h"

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif

template <typename R, typename... Args>
struct CustomOperator {
  using func_t = std::function<R(Args...)>;
  func_t f;

  friend std::ostream &operator<<(std::ostream &os, const CustomOperator &);

  R operator()(Args... args) const {
    if (f) return f(args...);
    ARTS_USER_ERROR("CustomOperator not set");
  }
};

#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

using NumericUnaryOperator = CustomOperator<Numeric, Numeric>;

using NumericTernaryOperator =
    CustomOperator<Numeric, Numeric, Numeric, Numeric>;

template <typename... WTs>
struct std::formatter<CustomOperator<WTs...>> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  template <typename... Ts>
  void make_compat(std::formatter<Ts> &...xs) const {
    tags.compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U> &x) {
    x._make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const CustomOperator<WTs...> &,
                              FmtContext &ctx) const {
    std::format_to(ctx, "?=f(?)");
    return ctx.out();
  }
};
