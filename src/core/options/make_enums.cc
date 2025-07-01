#include <fstream>
#include <iostream>
#include <span>

#include "arts_options.h"

namespace {
void create_headers() {
  std::ofstream common("enums-common-helper.h");
  common << R"--(#pragma once

#include <algorithm>
#include <array>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>

#include <format_tags.h>

using namespace std::literals;

template <typename T> constexpr bool good_enum(T x) noexcept = delete;
template <typename T> constexpr T to(const std::string_view x) = delete;

template <typename T>
struct enumdocs {
  static std::string_view str() noexcept = delete;
  static constexpr std::string_view name = "no-name"sv;
};

namespace enumstrs {
    template <typename T, int N> struct enum_str_data;
}  // namespace enumstrs
)--";

  for (auto& opt : internal_options()) {
    std::ofstream os("enums" + opt.name + ".h");
    os << R"--(#pragma once

#include "enums-common-helper.h"
)--";

    os << opt.head() << '\n';
    os << opt.tail() << '\n';
  }
}

void create_header() try {
  std::ofstream os("enums.h");
  const auto& opts = internal_options();

  os << R"--(#pragma once

#include "enums-common-helper.h"

)--";

  for (const auto& opt : opts) {
    os << "#include \"enums" << opt.name << ".h\"\n";
  }
} catch (const std::exception& e) {
  std::cerr << "Error creating enum header: " << e.what() << '\n';
  throw;
}

template <typename T>
std::span<const T> partial_span(const std::vector<T>& vec,
                                std::size_t N,
                                std::size_t i) {
  if (i >= N) {
    throw std::out_of_range("Partial span exceeds vector size");
  }

  const std::size_t size = vec.size() / N;
  if (i == N - 1) {
    return std::span<const T>(vec.begin() + i * size, vec.end());
  }
  return std::span<const T>(vec.begin() + i * size, size);
}

void create_cc() {
  std::ofstream os("enums.cpp");

  os << "#include \"enums.h\"\n\n";

  for (const auto& opt : internal_options()) {
    os << opt.impl() << '\n';
  }
}
}  // namespace

int main() try {
  create_header();
  create_headers();
  create_cc();
} catch (const std::exception& e) {
  std::cerr << "Error: " << e.what() << '\n';
  return 1;
}
