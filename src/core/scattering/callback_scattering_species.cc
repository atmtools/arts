#include "callback_scattering_species.h"

void xml_io_stream<CallbackScatteringSpecies>::write(std::ostream&                    os,
                                                     const CallbackScatteringSpecies& x,
                                                     bofstream*                       pbofs,
                                                     std::string_view                 name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);
  xml_write_to_stream(os, x.callback, pbofs);
  tag.write_to_end_stream(os);
}

void xml_io_stream<CallbackScatteringSpecies>::read(std::istream& is, CallbackScatteringSpecies& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);
  xml_read_from_stream(is, x.callback, pbifs);
  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
