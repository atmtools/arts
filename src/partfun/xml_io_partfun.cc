#include "xml_io_partfun.h"

#include "debug.h"

//! Reads Data from XML input stream
/*!
  \param is_xml  XML Input stream
  \param data    PartitionFunctions::Data return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_io_stream<PartitionFunctionsData>::read(std::istream& is_xml,
                                                 PartitionFunctionsData& data,
                                                 bifstream* pbifs) {
  XMLTag tag;

  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  String type;
  tag.get_attribute_value("type", type);
  data.type = to<PartitionFunctionsType>(type);

  xml_io_stream<Matrix>::read(is_xml, data.data, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
}

//! Writes PartitionFunctions::Data to XML output stream
/*!
 * \param os_xml  XML Output stream
 * \param data    PartitionFunctions::Data
 * \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 * \param name    Optional name attribute
 */
void xml_io_stream<PartitionFunctionsData>::write(
    std::ostream& os_xml,
    const PartitionFunctionsData& data,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag open_tag;
  XMLTag close_tag;
  std::ostringstream v;

  open_tag.name = type_name;
  if (name not_eq "") open_tag.add_attribute("name", String{name});
  open_tag.add_attribute("type", String(toString(data.type)));
  open_tag.write_to_stream(os_xml);

  xml_io_stream<Matrix>::write(os_xml, data.data, pbofs, "Data");

  close_tag.name = type_name;
  close_tag.write_to_end_stream(os_xml);
}

namespace PartitionFunctions {
Data data_read_file(const std::filesystem::path& path) {
  Data out;

  xml_read_from_file_base(path.string(), out);

  return out;
}
}  // namespace PartitionFunctions
