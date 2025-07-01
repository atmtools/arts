#include "retrieval_target.h"

void xml_io_stream<PairOfBlockMatrix>::write(std::ostream &os,
                                             const PairOfBlockMatrix &x,
                                             bofstream *pbofs,
                                             std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.first, pbofs);
  xml_write_to_stream(os, x.second, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<PairOfBlockMatrix>::read(std::istream &is,
                                            PairOfBlockMatrix &x,
                                            bifstream *pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.first, pbifs);
  xml_read_from_stream(is, x.second, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
