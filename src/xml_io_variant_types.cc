#include <xml_io_variant_types.h>

//! Helper macro for Variant types
#define TMPL_XML_READ_WRITE_STREAM_VAR(T)                                     \
  void xml_read_from_stream(std::istream& is_xml, T& var, bifstream* pbifs) { \
    xml_read(is_xml, var, pbifs);                                             \
  }                                                                           \
                                                                              \
  void xml_write_to_stream(std::ostream& os_xml,                              \
                           const T& var,                                      \
                           bofstream* pbofs,                                  \
                           const String& name) {                              \
    xml_write(os_xml, var, pbofs, name);                                      \
  }


TMPL_XML_READ_WRITE_STREAM_VAR(Atm::FieldData)
TMPL_XML_READ_WRITE_STREAM_VAR(AtmKeyVal)
