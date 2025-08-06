#include "retrieval_target.h"

void xml_io_stream<PairOfBlockMatrix>::write(std::ostream &os,
                                             const PairOfBlockMatrix &x,
                                             bofstream *pbofs,
                                             std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.first, pbofs);
  xml_write_to_stream(os, x.second, pbofs);

  tag.write_to_end_stream(os);
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
