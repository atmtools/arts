#include "debug.h"
#include "xml_io_partfun.h"

#include "xml_io_base.h"

//! Reads Data from XML input stream
/*!
  \param is_xml  XML Input stream
  \param data    PartitionFunctions::Data return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
namespace PartitionFunctions {
void xml_read_from_stream(std::istream& is_xml,
                          Data& data,
                          bifstream* pbifs) {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name("PartitionFunctionsData");

  String type;
  tag.get_attribute_value("type", type);
  data.type = PartitionFunctions::toTypeOrThrow(type);
  
  xml_read_from_stream(is_xml, data.data, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_name("/PartitionFunctionsData");
}
} // namespace PartitionFunctions

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
                         const String& name) {
  XMLTag open_tag;
  XMLTag close_tag;
  std::ostringstream v;
  
  open_tag.set_name("PartitionFunctionsData");
  if (name not_eq "") open_tag.add_attribute("name", name);
  open_tag.add_attribute("type", String(toString(data.type)));
  open_tag.write_to_stream(os_xml);
  os_xml << '\n';
  
  xml_write_to_stream(os_xml, data.data, pbofs, "Data");
  
  close_tag.set_name("/PartitionFunctionsData");
  close_tag.write_to_stream(os_xml);
  
  os_xml << '\n';
}

namespace PartitionFunctions {
Data data_read_file(const std::filesystem::path& path) {
  Data out;
  
  xml_read_from_file_base(String(path.native()), out);
  
  return out;
}
} // namespace PartitionFunctions
