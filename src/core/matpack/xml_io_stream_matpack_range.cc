#include "xml_io_stream_matpack_range.h"

void xml_io_stream<Range>::write(std::ostream &os,
                                 const Range &x,
                                 bofstream *pbofs,
                                 std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.nelem, pbofs);
  xml_write_to_stream(os, x.offset, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<Range>::read(std::istream &is, Range &x, bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.nelem, pbifs);
  xml_read_from_stream(is, x.offset, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
