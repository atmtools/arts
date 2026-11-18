#pragma once

#include <array>
#include <cstdint>
#include <ranges>
#include <string>
#include <string_view>

#include "format_tags.h"

using namespace std::literals;

namespace stdr = std::ranges;

/* FNV-1a hash used by the generated to<Enum>() name lookups.  The
 * generated code switches on enumhash(x), so resolving a name costs
 * O(len) constexpr steps instead of a linear scan over all the enum's
 * names.  The hash only narrows the search; each case still verifies
 * the full string, so collisions cannot change the result. */
constexpr std::uint32_t enumhash(const std::string_view s) noexcept {
  std::uint32_t h = 2166136261u;
  for (const char c : s)
    h = (h ^ static_cast<std::uint32_t>(static_cast<unsigned char>(c))) *
        16777619u;
  return h;
}

template <typename T> constexpr bool good_enum(T x) noexcept      = delete;
template <typename T> constexpr T    to(const std::string_view x) = delete;

template <typename T> struct enumdocs {
  static std::string_view           str() noexcept = delete;
  static constexpr std::string_view name           = "no-name"sv;
};

namespace enumstrs {
template <typename T, int N> struct enum_str_data;
}  // namespace enumstrs
