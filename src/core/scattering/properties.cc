#include "properties.h"

#include <xml_io_base.h>
#include <xml_io_stream_core.h>

std::ostream& operator<<(std::ostream& os,
                         const ScatteringSpeciesProperty& ssp) {
  return os << ssp.species_name << "_" << ssp.pproperty;
}

ScatteringSpeciesProperty ScatteringSpeciesProperty::from_string(
    const std::string_view in) {
  std::string_view s = in;
  const auto pos     = s.find_last_of('_');
  if (pos == std::string_view::npos) {
    throw std::runtime_error(std::format(
        R"-x-(Bad input "{}"

Must be of the form "<species_name>_<particulate_property>"
)-x-",
        in));
  }

  const std::string_view species_name = s.substr(0, pos);
  const std::string_view pproperty    = s.substr(pos + 1);

  return {.species_name = std::string{species_name},
          .pproperty    = to<ParticulateProperty>(pproperty)};
}

void xml_io_stream<ScatteringSpeciesProperty>::write(
    std::ostream& os,
    const ScatteringSpeciesProperty& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.species_name, pbofs, "Species name"sv);
  xml_write_to_stream(os, x.pproperty, pbofs, "The particulate property"sv);

  tag.write_to_end_stream(os);
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
