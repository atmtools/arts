#include "jpl_species_info.h"

#include <isotopologues.h>

#include <algorithm>
#include <cassert>

#include "configtypes.h"
#include "xml_io_base.h"

void xml_io_stream<JplSpeciesInfo>::write(std::ostream&         os,
                                          const JplSpeciesInfo& x,
                                          bofstream*,
                                          std::string_view name) try {
  XMLTag tag(type_name,
             "name",
             name,
             "id",
             x.id,
             "spec",
             x.spec.FullName(),
             "has_qn",
             std::format("{}", x.has_qn),
             "QT0",
             x.QT0,
             "T0",
             x.T0);
  tag.write_to_stream(os);

  tag.write_to_end_stream(os);
}
ARTS_METHOD_ERROR_CATCH

void xml_io_stream<JplSpeciesInfo>::read(std::istream& is, JplSpeciesInfo& x, bifstream*) try {
  XMLTag tag, tag2;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  tag.get_attribute_value("id", x.id);

  String v;
  tag.get_attribute_value("spec", v);
  x.spec = SpeciesIsotope::from_name(v);

  tag.get_attribute_value("has_qn", v);
  std::ranges::transform(v, v.begin(), [](unsigned char c) { return std::tolower(c); });
  x.has_qn = v == "true";

  tag.get_attribute_value("QT0", x.QT0);
  tag.get_attribute_value("T0", x.T0);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
ARTS_METHOD_ERROR_CATCH
