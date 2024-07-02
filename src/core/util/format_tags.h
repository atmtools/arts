#pragma once

#include <configtypes.h>

#include <concepts>
#include <format>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <variant>

using namespace std::literals;

inline constexpr Index short_str_v_stp = 3;
inline constexpr Index short_str_v_cut = 8;

struct format_tags {
  bool names     = false;
  bool comma     = false;
  bool bracket   = false;
  bool short_str = false;

  template <typename T>
  constexpr void compat(std::formatter<T>& x) const {
    if constexpr (requires { x.inner_fmt().tags = *this; }) {
      x.inner_fmt().tags = *this;
    } else if constexpr (requires { x = *this; }) {
      x = *this;
    } else if constexpr (requires { x.tags = *this; }) {
      x.tags = *this;
    }
  }

  template <typename... Ts>
  constexpr void compat(std::formatter<Ts>&... x) const
    requires(sizeof...(Ts) > 1)
  {
    (compat(x), ...);
  }
};

constexpr std::format_parse_context::iterator parse_format_tags(
    format_tags& fmt, std::format_parse_context& ctx) {
  auto&& it = ctx.begin();

  while (it != ctx.end() and *it != '}') {
    if (*it == ',') {
      fmt.comma = true;
      ++it;
      continue;
    }

    if (*it == 'B') {
      fmt.bracket = true;
      ++it;
      continue;
    }

    if (*it == 's') {
      fmt.short_str = true;
      ++it;
      continue;
    }

    if (*it == 'N') {
      fmt.names = true;
      ++it;
      continue;
    }

    throw std::format_error("Invalid format args for arts-type");
  }

  return it;
}

template <typename T>
concept arts_inner_fmt =
    requires(std::formatter<T> x) { x.inner_fmt().tags; } and
    std::same_as<
        format_tags,
        std::remove_cvref_t<decltype(std::formatter<T>{}.inner_fmt().tags)>>;

template <typename T>
concept arts_formatter_compat =
    requires(std::formatter<T> x, std::formatter<T> y) {
      x.compat(y);
      x.make_compat(y);
    };

template <typename T>
concept arts_formattable =
    arts_inner_fmt<T> and arts_formatter_compat<T> and requires(T x) {
      std::format("{}", x);
      std::format("{:sNB,}", x);
    };

template <typename T>
concept arts_formattable_or_value_type =
    arts_formattable<T> or std::integral<T> or std::floating_point<T>;

template <arts_formattable_or_value_type... WTs>
struct std::formatter<std::variant<WTs...>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  template <typename T>
  constexpr void compat(const std::formatter<T>& x) {
    x.make_compat(*this);
  }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const std::variant<WTs...>& v,
                              FmtContext& ctx) const {
    std::visit(
        [*this, &ctx]<typename T>(const T& x) {
          std::formatter<T> fmt;
          make_compat(fmt);
          fmt.format(x, ctx);
        },
        v);
    return ctx.out();
  }
};

template <arts_formattable_or_value_type Key,
          arts_formattable_or_value_type Value>
struct std::formatter<std::unordered_map<Key, Value>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }

  [[nodiscard]] constexpr const auto& inner_fmt() const { return *this; }

  template <typename T>
  constexpr void compat(const std::formatter<T>& x) {
    x.make_compat(*this);
  }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
    tags.compat(xs...);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const std::unordered_map<Key, Value>& v,
                              FmtContext& ctx) const {
    std::formatter<Key> key{};
    std::formatter<Value> value{};

    make_compat(key, value);
    const std::string_view sep = tags.comma ? ", "sv : " "sv;

    if (tags.bracket) std::ranges::copy("{"sv, ctx.out());
    for (const auto& [k, val] : v) {
      key.format(k, ctx);
      std::ranges::copy(": "sv, ctx.out());
      value.format(val, ctx);
      std::ranges::copy(sep, ctx.out());
    }
    if (tags.bracket) std::ranges::copy("}"sv, ctx.out());

    return ctx.out();
  }
};

template <arts_formattable_or_value_type T, class Allocator>
struct std::formatter<std::vector<T, Allocator>> {
  std::formatter<T> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
    inner_fmt().make_compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const std::vector<T, Allocator>& v,
                              FmtContext& ctx) const {
    using std::ranges::views::take, std::ranges::views::drop;

    const auto n = v.size();

    bool first = true;

    if (inner_fmt().tags.bracket) std::ranges::copy("["sv, ctx.out());

    if (inner_fmt().tags.short_str and n > short_str_v_cut) {
      for (auto&& a : v | take(short_str_v_stp)) {
        if (not first and inner_fmt().tags.comma)
          std::ranges::copy(","sv, ctx.out());
        if (not first) std::ranges::copy("\n"sv, ctx.out());
        fmt.format(a, ctx);
        first = false;
      }

      if (inner_fmt().tags.comma) std::ranges::copy(","sv, ctx.out());
      std::ranges::copy("\n..."sv, ctx.out());

      for (auto&& a : v | drop(n - short_str_v_stp)) {
        if (not first and inner_fmt().tags.comma)
          std::ranges::copy(","sv, ctx.out());
        if (not first) std::ranges::copy("\n"sv, ctx.out());
        fmt.format(a, ctx);
        first = false;
      }
    } else {
      for (auto&& a : v) {
        if (not first and inner_fmt().tags.comma)
          std::ranges::copy(","sv, ctx.out());
        if (not first) std::ranges::copy("\n"sv, ctx.out());
        fmt.format(a, ctx);
        first = false;
      }
    }

    if (inner_fmt().tags.bracket) std::ranges::copy("]"sv, ctx.out());

    return ctx.out();
  }
};
