#include "xml_io_stream_array.h"

#include "xml_io_stream_core.h"

void xml_io_stream<ArrayOfString>::write(std::ostream& os,
                                         const ArrayOfString& x,
                                         bofstream* pbofs,
                                         std::string_view name) {
  xml_io_stream_array<String>::write(os, x, pbofs, name);
}

void xml_io_stream<ArrayOfString>::read(std::istream& is,
                                        ArrayOfString& x,
                                        bifstream* pbifs) {
  xml_io_stream_array<String>::read(is, x, pbifs);
}