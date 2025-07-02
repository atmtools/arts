#pragma once

#include <matpack.h>
#include <xml.h>

#include <functional>
#include <iosfwd>

#include "debug.h"

#ifndef _MSC_VER
#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-template-friend"
#endif
#endif

template <typename R, typename... Args>
struct CustomOperator {
  using func_t = std::function<R(Args...)>;
  func_t f{[](Args...) -> R {
    throw std::runtime_error("CustomOperator not set");
    std::unreachable();
  }};

  friend std::ostream &operator<<(std::ostream &os, const CustomOperator &) {
    return os << "custom-operator";
  }

  R operator()(Args... args) const {
    if (f) return f(args...);
    ARTS_USER_ERROR("CustomOperator not set");
  }
};

#ifndef _MSC_VER
#if defined(__clang__)
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#endif

using NumericUnaryOperator  = CustomOperator<Numeric, Numeric>;
using NumericBinaryOperator = CustomOperator<Numeric, Numeric, Numeric>;
using NumericTernaryOperator =
    CustomOperator<Numeric, Numeric, Numeric, Numeric>;

template <typename... WTs>
struct std::formatter<CustomOperator<WTs...>> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const CustomOperator<WTs...> &,
                              FmtContext &ctx) const {
    return std::format_to(ctx.out(), "{0}functional-data{0}"sv, tags.quote());
  }
};

template <typename R, typename... Args>
struct xml_io_stream<CustomOperator<R, Args...>> {
  static constexpr std::string_view type_name = "CustomOperator"sv;

  static void write(std::ostream &os,
                    const CustomOperator<R, Args...> &x,
                    bofstream *pbofs      = nullptr,
                    std::string_view name = ""sv) {
    std::println(os, R"(<{0} name="{1}">)", type_name, name);

    xml_write_to_stream(os, x.f, pbofs);

    std::println(os, R"(</{0}>)", type_name);
  }

  static void read(std::istream &is,
                   CustomOperator<R, Args...> &x,
                   bifstream *pbifs = nullptr) {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    xml_read_from_stream(is, x.f, pbifs);

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
};
