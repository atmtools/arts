#include "debug.h"
#include "xml_io_partfun.h"

#include "xml_io_base.h"

//! Reads Data from XML input stream
/*!
  \param is_xml  XML Input stream
  \param data    PartitionFunctions::Data return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(std::istream& is_xml,
                          PartitionFunctions::Data& data,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  XMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("PartitionFunctionsData");

  String type;
  tag.get_attribute_value("type", type);
  data.type = PartitionFunctions::toTypeOrThrow(type);
  
  xml_read_from_stream(is_xml, data.data, pbifs, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/PartitionFunctionsData");
}

//! Writes PartitionFunctions::Data to XML output stream
/*!
 * \param os_xml  XML Output stream
 * \param data    PartitionFunctions::Data
 * \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 * \param name    Optional name attribute
 */
void xml_write_to_stream(std::ostream& os_xml,
                         const PartitionFunctions::Data& data,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  XMLTag open_tag(verbosity);
  XMLTag close_tag(verbosity);
  std::ostringstream v;
  
  open_tag.set_name("PartitionFunctionsData");
  if (name not_eq "") open_tag.add_attribute("name", name);
  open_tag.add_attribute("type", String(toString(data.type)));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  
  xml_write_to_stream(os_xml, data.data, pbofs, "Data", verbosity);
  
  close_tag.set_name("/PartitionFunctionsData");
  close_tag.write_to_stream(os_xml);
  
  os_xml << '\n';
}

namespace PartitionFunctions {
Data data_read_file(const std::filesystem::path& path) {
  Data out;
  
  xml_read_from_file_base(String(path.native()), out, Verbosity());
  
  return out;
}
} // namespace PartitionFunctions
