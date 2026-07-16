#include "species_enum_info.h"

void xml_io_stream<SpeciesEnumInfo>::write(std::ostream& os, const SpeciesEnumInfo& v, bofstream*, std::string_view) {
  if (v.shortname.empty() or v.longname.empty()) {
    throw std::runtime_error("Cannot write SpeciesEnumInfo with empty shortname or longname.");
  }

  if (v.shortname.find(' ') != std::string::npos) {
    throw std::runtime_error(std::format("Shortname cannot contain spaces: '{}'", v.shortname));
  }

  if (v.longname.find(' ') != std::string::npos) {
    throw std::runtime_error(std::format("Longname cannot contain spaces: '{}'", v.longname));
  }

  std::println(os, R"(<{0}> {3} {1} {2} </{0}>)"sv, type_name, v.shortname, v.longname, v.enum_value);
}

void xml_io_stream<SpeciesEnumInfo>::read(std::istream& is, SpeciesEnumInfo& v, bifstream*) try {
  std::string val{};
  is >> val;

  if (val != std::format("<{}>", type_name)) {
    throw std::runtime_error(std::format("Expected opening tag <{}>, but got '{}'", type_name, val));
  }

  is >> v.enum_value >> v.shortname >> v.longname;

  if (v.enum_value < 1) {
    throw std::runtime_error(std::format("Enum value must be positive, but got {}", v.enum_value));
  }

  is >> val;

  if (val != std::format("</{}>", type_name)) {
    throw std::runtime_error(std::format("Expected closing tag </{}>, but got '{}'", type_name, val));
  }
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("Error reading {}: {}", type_name, e.what()));
}
