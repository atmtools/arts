#include "xml_io_stream_matpack_range.h"

void xml_io_stream<Range>::write(std::ostream &os,
                                 const Range &x,
                                 bofstream *pbofs,
                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);
  xml_write_to_stream(os, x.nelem, pbofs);
  xml_write_to_stream(os, x.offset, pbofs);
  std::println(os, R"(/<{0}>)", type_name);
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
