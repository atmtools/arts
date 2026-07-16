#include <enumsSpeciesEnum.h>
#include <nonstd.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <print>

#include "species_info.h"

namespace {
std::set<SpeciesIsotopologueInfo> read_split_species(
    const std::filesystem::path& path) {
  std::set<SpeciesIsotopologueInfo> species_set{};

  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".xml") {
      SpeciesIsotopologueInfo info;
      xml_read_from_file(entry.path().string(), info);

      if (species_set.count(info) > 0) {
        throw std::runtime_error(std::format(
            "Duplicate species isotopologue info found: {}", info.name()));
      }

      species_set.insert(info);
    }
  }

  return species_set;
}

void write_header(std::ostream& os,
                  const std::set<SpeciesIsotopologueInfo>& data) {
  std::print(os,
             R"(#pragma once

#include <isotopologues_class.h>
#include <species.h>

namespace Species {{
inline constexpr std::array Isotopologues{{)");

  std::size_t done{0};

  for (const auto& info : data) {
    const auto spec_index = static_cast<std::size_t>(info.species);

    while (done <= spec_index) {
      const auto spec = static_cast<SpeciesEnum>(done);

      std::println(os,
                   R"(
  Isotope{{.spec="{}"_spec}},)",
                   spec);
      done++;
    }

    if (nonstd::isdigit(info.code.front())) {
      std::println(
          os,
          R"(  Isotope{{.spec="{}"_spec, .isotname="{}"sv, .mass={}, .builtin_ratio={}, .gi={}}},)",
          info.species,
          info.code,
          info.mass,
          info.default_ratio,
          info.degeneracy);
    } else {
      std::println(os,
                   R"(  Isotope{{.spec="{}"_spec, .isotname="{}"sv}},)",
                   info.species,
                   info.code);
    }
  }

  SpeciesEnum specrest = static_cast<SpeciesEnum>(done);
  while (good_enum(specrest)) {
    const auto spec = static_cast<SpeciesEnum>(done);

    std::println(os,
                 R"(
  Isotope{{.spec="{}"_spec}},)",
                 spec);
    done++;
    specrest = static_cast<SpeciesEnum>(done);
  }

  std::println(os, R"(}};
}} // namespace Species)");
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) {
    std::println(
        stderr, "Usage: {} <arts-cat-data>", argv[0]);
    return EXIT_FAILURE;
  }

  const std::filesystem::path path = argv[1];
  auto species = read_split_species(path / "isotopologues/");

  if (species.empty()) {
    std::println(stderr,
                 "No species isotopologue info found in {}/{}",
                 argv[1],
                 "isotopologues/");
    return EXIT_FAILURE;
  }

  std::ofstream h("auto_isotopologues.h", std::ios::out);
  write_header(h, species);

  return EXIT_SUCCESS;
} catch (const std::exception& e) {
  std::println(stderr, "Error: {}", e.what());
  return EXIT_FAILURE;
} catch (...) {
  std::println(stderr, "An unknown error occurred.");
  return EXIT_FAILURE;
}
