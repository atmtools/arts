#include "properties.h"

#include <xml_io_base.h>
#include <xml_io_stream_core.h>

std::ostream& operator<<(std::ostream& os,
                                  const ScatteringSpeciesProperty& ssp) {
  return os << ssp.species_name << "_" << ssp.pproperty;
}

void xml_io_stream<ScatteringSpeciesProperty>::write(
    std::ostream& os,
    const ScatteringSpeciesProperty& x,
    bofstream* pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.species_name, pbofs, "Species name"sv);
  xml_write_to_stream(os, x.pproperty, pbofs, "The particulate property"sv);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<ScatteringSpeciesProperty>::read(
    std::istream& is, ScatteringSpeciesProperty& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.species_name, pbifs);
  xml_read_from_stream(is, x.pproperty, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
