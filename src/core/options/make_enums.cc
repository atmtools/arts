#include <fstream>
#include <iostream>

#include "arts_options.h"

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

class bifstream;
class bofstream;

template <typename T> constexpr bool good_enum(T x) noexcept = delete;
template <typename T> constexpr T to(const std::string_view x) = delete;

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

void create_cc() {
  std::ofstream os("enums.cpp");

  const auto& opts = internal_options();

  os << "#include \"enums.h\"\n\n";

  for (const auto& opt : opts) {
    os << opt.impl() << '\n';
  }
}

int main() try {
  create_header();
  create_headers();
  create_cc();
} catch (const std::exception& e) {
  std::cerr << "Error: " << e.what() << '\n';
  return 1;
}
