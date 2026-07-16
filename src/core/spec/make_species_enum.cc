#include <filesystem>
#include <print>
#include <set>

#include "species_enum_info.h"

namespace {
std::set<SpeciesEnumInfo> read_split_species(const std::filesystem::path& path) {
  std::set<SpeciesEnumInfo> species_set{};

  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".xml") {
      SpeciesEnumInfo info;
      xml_read_from_file(entry.path().string(), info);

      if (species_set.count(info) > 0) {
        throw std::runtime_error(std::format("Duplicate species enum info found: {}", info.shortname));
      }

      species_set.insert(info);
    }
  }

  return species_set;
}

void write_header(std::ostream& os) {
  std::println(os,
               R"(#pragma once

#include <species_enum_info.h>

#include <vector>

const std::vector<SpeciesEnumInfo>& get_species_enum_info();
)");
}

void write_source(std::ostream& os, const std::set<SpeciesEnumInfo>& species_set) {
  std::println(os,
               R"(#include "auto_species_enum_info.h"

namespace {{
std::vector<SpeciesEnumInfo> get_species_enum_info_local() {{
  std::vector<SpeciesEnumInfo> info{{}};
  info.reserve({});
)",
               species_set.size());

  for (const auto& info : species_set) {
    std::println(os, R"(  info.emplace_back({}, "{}", "{}");)", info.enum_value, info.shortname, info.longname);
  }

  std::println(os, R"(
  return info;
}}
}}  // namespace

const std::vector<SpeciesEnumInfo>& get_species_enum_info() {{
  const static auto& info = get_species_enum_info_local();
  return info;
}})");
}

void write_isot_source(std::ostream& os, const std::set<SpeciesEnumInfo>& species_set) {
  std::println(os,
               R"(#include "isotopologues.h"

#include <debug.h>

namespace Species {{
ArrayOfSpeciesIsotope isotopologues(SpeciesEnum spec) {{
  switch (spec) {{
    case SpeciesEnum::Bath: break;)");

  for (const auto& info : species_set) {
    std::println(os,
                 R"(    case SpeciesEnum::{0}: {{
      static constexpr auto v = isotopologues<SpeciesEnum::{0}>();
      return {{v.begin(), v.end()}};
    }} break;)",
                 info.longname);
  }
  std::println(os,
               R"(    case SpeciesEnum::unused: break;
  }}

  ARTS_USER_ERROR("Cannot understand: {{}}", spec)
}}
}}  // namespace Species)");
}

void write_partfun_wrap_header(std::ostream& os, const std::set<SpeciesEnumInfo>& species_set) {
  std::println(os,
               R"(#include "partfun.h"

#include <auto_partfun.h>
#include <debug.h>

namespace PartitionFunctions {{
template <Derivatives d>
Numeric partfun_impl(Numeric T, const SpeciesIsotope& ir) {{
  switch (ir.spec) {{
    case SpeciesEnum::Bath: break;)");

  for (const auto& info : species_set) {
    std::println(os, R"(    case SpeciesEnum::{0}: return compute{0}<d>(T, ir.isotname);)", info.longname);
  }
  std::println(os,
               R"(    case SpeciesEnum::unused: break;
  }}

  ARTS_USER_ERROR("This is not a valid isotopologue: {{}}", ir)
}}

bool has_partfun(const SpeciesIsotope& ir) noexcept {{
  switch (ir.spec) {{
    case SpeciesEnum::Bath: break;)");
  for (const auto& info : species_set) {
    std::println(
        os, R"(    case SpeciesEnum::{0}: return std::ranges::binary_search(has{0}, ir.isotname);)", info.longname);
  }
  std::println(os,
               R"(    case SpeciesEnum::unused: break;
  }}

  return false;
}}
}}  // namespace PartitionFunctions)");
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) {
    std::println(stderr, "Usage: {} <arts-cat-data>", argv[0]);
    return EXIT_FAILURE;
  }

  const std::filesystem::path path        = argv[1];
  auto                        species_set = read_split_species(path / "species/");

  if (species_set.empty()) {
    std::println(stderr, "No species enum info found in {}/{}", argv[1], "species/");
    return EXIT_FAILURE;
  }

  auto h   = std::fstream("auto_species_enum_info.h", std::ios::out);
  auto cc  = std::fstream("auto_species_enum_info.cc", std::ios::out);
  auto cc2 = std::fstream("auto_isotopologues.cc", std::ios::out);
  auto h2  = std::fstream("auto_partition_functions.h", std::ios::out);

  write_source(cc, species_set);
  write_header(h);

  write_isot_source(cc2, species_set);
  write_partfun_wrap_header(h2, species_set);

} catch (const std::exception& e) {
  std::println(stderr, "Error: {}", e.what());
  return EXIT_FAILURE;
} catch (...) {
  std::println(stderr, "An unknown error occurred.");
  return EXIT_FAILURE;
}
