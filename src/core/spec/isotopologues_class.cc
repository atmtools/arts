#include "isotopologues_class.h"

#include <debug.h>
#include <xml_io_base.h>

#include <cassert>

void xml_io_stream<SpeciesIsotope>::write(std::ostream&         os,
                                          const SpeciesIsotope& x,
                                          bofstream*,
                                          std::string_view name) try {
  XMLTag tag(type_name, "name", name, "isot", x.FullName());
  tag.write_to_stream(os);
  tag.write_to_end_stream(os);
}
ARTS_METHOD_ERROR_CATCH

void xml_io_stream<SpeciesIsotope>::read(std::istream& is, SpeciesIsotope& x, bifstream*) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String v;
  tag.get_attribute_value("isot", v);
  x = SpeciesIsotope::from_name(v);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
ARTS_METHOD_ERROR_CATCH
