#include "hitran_species_info.h"

#include <isotopologues.h>

#include <cassert>

void xml_io_stream<HitranSpeciesInfo>::write(std::ostream&            os,
                                             const HitranSpeciesInfo& x,
                                             bofstream*,
                                             std::string_view name) try {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  std::println(os, "{}", x);
  tag.write_to_end_stream(os);
}
ARTS_METHOD_ERROR_CATCH

void xml_io_stream<HitranSpeciesInfo>::read(std::istream& is, HitranSpeciesInfo& x, bifstream*) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String spec_name;
  is >> spec_name >> x.hitind >> x.hitchar >> x.ratio;
  x.spec = SpeciesIsotope::from_name(spec_name);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
ARTS_METHOD_ERROR_CATCH
