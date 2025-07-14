////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#include "xml_io.h"

#include <bifstream.h>
#include <bofstream.h>
#include <config.h>
#include <debug.h>
#include <file.h>

#include <filesystem>

void xml_find_and_open_input_file(std::shared_ptr<std::istream>& ifs,
                                  const String& filename) {
  String xml_file = filename;
  find_xml_file(xml_file);

  // Open input stream:
  if (xml_file.substr(xml_file.length() - 3, 3) == ".gz")
#ifdef ENABLE_ZLIB
  {
    ifs = std::shared_ptr<std::istream>(new igzstream());
    xml_open_input_file(*(static_cast<igzstream*>(ifs.get())), xml_file);
  }
#else
  {
    throw std::runtime_error(
        "This arts version was compiled without zlib support.\n"
        "Thus zipped xml files cannot be read.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ifs = std::shared_ptr<std::istream>(new std::ifstream());
    xml_open_input_file(*(static_cast<std::ifstream*>(ifs.get())), xml_file);
  }

  // Read the file into memory first to significantly speed up
  // the parsing (13x to 18x faster).
  std::shared_ptr<std::stringstream> buffer(new std::stringstream());
  *buffer << ifs->rdbuf();
  ifs = buffer;
}

////////////////////////////////////////////////////////////////////////////
//   General XML functions (file header, start root tag, end root tag)
////////////////////////////////////////////////////////////////////////////

//! Throws XML parser runtime error
/*!
  This is used quite often inside the data parsing routines so it's a
  function for itself.

  \param tag        XMLTag
  \param str_error  Error description
*/
void xml_data_parse_error(XMLTag& tag, String str_error) {
  std::ostringstream os;
  os << "XML data parse error: Error reading ";
  tag.write_to_stream(os);
  os << str_error << "\n"
     << "Check syntax of XML file. A possible cause is that the file "
     << "contains NaN or Inf values.\n";
  throw std::runtime_error(os.str());
}

////////////////////////////////////////////////////////////////////////////
//   Default file name
////////////////////////////////////////////////////////////////////////////

//! Gives the default filename for the XML formats.
/*!
  The default name is only used if the filename is empty.

  \param filename filename
  \param varname variable name
*/
void filename_xml(const String& filename) {
  ARTS_USER_ERROR_IF(filename == "", "Must have filename")
}

//! Gives the default filename, with file index, for the XML formats.
/*!
  The default name is only used if the filename is empty.

  \param[out] filename   filename
  \param[in]  file_index Index appended to the filename
  \param[in]  digits     Width for padding with zeros
*/
void filename_xml_with_index(String& filename,
                             const Index& file_index,
                             const Index&) {
  ARTS_USER_ERROR_IF("" == filename, "Must have filename")
  filename = std::format("{}{}{}{}", filename, ".", file_index, ".xml");
}

String complete_basename(const String& basename) {
  if (basename.back() == '/') return basename;
  if (basename.back() == '.') return basename;
  if (std::filesystem::is_directory(basename)) return basename + "/";
  return basename + ".";
}
