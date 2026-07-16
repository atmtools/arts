#include "species_info.h"

#include <double_imanip.h>

String SpeciesIsotopologueInfo::name() const { return std::format("{}-{}", species, code); }

void xml_io_stream<SpeciesIsotopologueInfo>::write(std::ostream&                  os,
                                                   const SpeciesIsotopologueInfo& v,
                                                   bofstream*,
                                                   std::string_view) {
  std::println(
      os, R"(<{0}> {1} {2} {3} {4} {5} </{0}>)", type_name, v.species, v.code, v.mass, v.default_ratio, v.degeneracy);
}

void xml_io_stream<SpeciesIsotopologueInfo>::read(std::istream& is, SpeciesIsotopologueInfo& v, bifstream*) try {
  std::string val{};
  is >> val;

  if (val != std::format("<{}>", type_name)) {
    throw std::runtime_error(std::format("Expected opening tag <{}>, but got '{}'", type_name, val));
  }

  is >> v.species >> v.code >> double_imanip{} >> v.mass >> v.default_ratio;
  is >> v.degeneracy;

  if (v.code.empty()) { throw std::runtime_error("Cannot read SpeciesIsotopologueInfo with empty species or code."); }

  is >> val;

  if (val != std::format("</{}>", type_name)) {
    throw std::runtime_error(std::format("Expected closing tag </{}>, but got '{}'", type_name, val));
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("Error reading {}: {}", type_name, e.what()));
}
