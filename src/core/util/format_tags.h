#pragma once

#include <format>

#include <configtypes.h>

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
