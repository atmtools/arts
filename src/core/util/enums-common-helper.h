#pragma once

#include <algorithm>
#include <array>
#include <istream>
#include <ostream>
#include <ranges>
#include <string>
#include <string_view>

#include "format_tags.h"

using namespace std::literals;

namespace stdr = std::ranges;

template <typename T> constexpr bool good_enum(T x) noexcept      = delete;
template <typename T> constexpr T    to(const std::string_view x) = delete;

template <typename T> struct enumdocs {
  static std::string_view           str() noexcept = delete;
  static constexpr std::string_view name           = "no-name"sv;
};

namespace enumstrs {
template <typename T, int N> struct enum_str_data;
}  // namespace enumstrs
