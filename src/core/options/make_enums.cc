#include <fstream>
#include <iostream>

#include "arts_options.h"

void create_header(std::ostream& os) try {
  const auto& opts = internal_options();

  os << R"--(#pragma once

#include <array>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>

using namespace std::literals;

)--";

  for (const auto& opt : opts) {
    os << opt.head() << '\n';
  }

  os << R"--(template <typename T>
concept ValidArtsOption = false)--";

  for (const auto& opt : opts) {
    os << "\n  or std::same_as<T, " << opt.name << ">";
  }

  os << R"--(;

template <ValidArtsOption T> constexpr bool good_enum(T x) noexcept;
template <ValidArtsOption T> constexpr T to(const std::string_view x);

namespace enumstrs {
    template <ValidArtsOption T, int N> struct enum_str_data;
}  // namespace enumstrs

)--";

  for (const auto& opt : opts) {
    os << opt.tail() << '\n';
  }
} catch (const std::exception& e) {
  std::cerr << "Error creating enum header: " << e.what() << '\n';
  throw;
}

void create_cc(std::ostream& os) {
  const auto& opts = internal_options();

  os << "#include \"enums.h\"\n\n";

  for (const auto& opt : opts) {
    os << opt.impl() << '\n';
  }
}

int main() try {
  std::ofstream header("enums.h");
  std::ofstream cc("enums.cpp");

  create_header(header);
  create_cc(cc);

} catch (const std::exception& e) {
  std::cerr << "Error: " << e.what() << '\n';
  return 1;
}
