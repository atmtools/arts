#include <cstdlib>
#include <filesystem>
#include <ostream>

#include "jpl_species_info.h"

namespace {
std::vector<JplSpeciesInfo> read(const std::filesystem::path& path) {
  std::vector<JplSpeciesInfo> out;

  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".xml") {
      xml_read_from_file(entry.path().string(), out.emplace_back());
    }
  }

  std::ranges::sort(out, {}, &JplSpeciesInfo::id);
  auto it = std::ranges::adjacent_find(out, {}, &JplSpeciesInfo::id);
  if (it != out.end()) { throw std::runtime_error(std::format("Duplicate ID found in JPL data: {}", it->id)); }

  return out;
}

void write_header(std::ostream& os, const std::vector<JplSpeciesInfo>& data) {
  std::println(os,
               R"(#pragma once

#include <array>

#include "jpl_species_info.h"

namespace Jpl {{
inline constexpr std::array<JplSpeciesInfo, {}> jpl_data{{)",
               data.size());

  for (const auto& [id, spec, has_qn, QT0, T0] : data) {
    std::println(os,
                 R"(   JplSpeciesInfo{{ .id={0}, .spec="{1}"_isot, .has_qn={2}, .QT0={3}, .T0={4} }}, )",
                 id,
                 spec,
                 has_qn,
                 QT0,
                 T0);
  }

  std::println(os, R"(}};
}}  // namespace Jpl)");
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) {
    std::println(stderr, "Usage: {} <arts-cat-data>", argv[0]);
    return EXIT_FAILURE;
  }

  const std::filesystem::path path = argv[1];
  auto                        data = read(path / "jpl/");

  std::ofstream h("auto_jpl_species_map.h", std::ios::out);

  write_header(h, data);

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "Error in main: {}", e.what());
  return EXIT_FAILURE;
}
