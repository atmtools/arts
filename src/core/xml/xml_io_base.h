////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#pragma once

#include <config.h>
#include <enumsFileType.h>
#include <mystring.h>
#include <time_report.h>

#include <memory>
#include <sstream>

#include "xml_io_stream.h"

#ifdef ENABLE_ZLIB
#include <gzstream.h>
#endif

enum NumericType { NUMERIC_TYPE_FLOAT, NUMERIC_TYPE_DOUBLE };
enum EndianType { ENDIAN_TYPE_LITTLE, ENDIAN_TYPE_BIG };

////////////////////////////////////////////////////////////////////////////
//   XML parser classes
////////////////////////////////////////////////////////////////////////////

//! The ARTS XML tag class
/*!
  Handles reading, writing and constructing of XML tags.
*/
struct XMLTag {
  String name;                                /*!< Tag name */
  std::unordered_map<String, String> attribs; /*!< List of attributes */

  XMLTag() = default;

  void check_name(const std::string_view& expected_name);
  void check_end_name(const std::string_view& expected_name);

  void add_attribute(const std::string_view& aname,
                     const std::string_view& value);
  void add_attribute(const std::string_view& aname, const Index& value);
  void add_attribute(const std::string_view& aname, const Size& value);
  void add_attribute(const std::string_view& aname, const Numeric& value);

  template <typename... Args>
  void add_attribute(const std::string_view& name,
                     const auto& value,
                     const Args&... args)
    requires(sizeof...(Args) >= 2 and (sizeof...(Args) % 2 == 0))
  {
    add_attribute(name, value);
    add_attribute(args...);
  }

  template <typename... Args>
  XMLTag(std::string_view n, const Args&... args) : name(n), attribs{} {
    if (sizeof...(Args) != 0) add_attribute(args...);
  }

  void check_attribute(const std::string_view& aname,
                       const std::string_view& value);
  void check_attribute(const std::string_view& aname, const Index& value);
  void check_attribute(const std::string_view& aname, const Size& value);

  void get_attribute_value(const std::string_view& aname,
                           String& value,
                           std::string_view def = "");
  void get_attribute_value(const std::string_view& aname, Index& value);
  void get_attribute_value(const std::string_view& aname, Size& value);

  /** Returns value of attribute as type Numeric
   * 
   * Searches for the matching attribute and returns it value. If no
   * attribute with the given name exists, return value is set to
   * -1e99.
   * 
   * @param[in] aname Attribute name
   * @param[out] value Return value
   */
  void get_attribute_value(const std::string_view& aname, Numeric& value);

  void read_from_stream(std::istream& is);

  void write_to_stream(std::ostream& os);
  void write_to_end_stream(std::ostream& os);

  /** Returns if the attribute exists or not
   * 
   * @param[in] aname Attribute name
   * @return bool Does this attribute exist?
   */
  [[nodiscard]] bool has_attribute(const std::string_view& aname) const;
};

////////////////////////////////////////////////////////////////////////////
//   General XML handling routines
////////////////////////////////////////////////////////////////////////////

void xml_parse_error(const String& str_error);

void xml_data_parse_error(XMLTag& tag, const String& str_error);

void xml_set_stream_precision(std::ostream& os);

void parse_xml_tag_content_as_string(std::istream& is_xml, String& content);

void xml_parse_from_stream(std::istream&, ArrayOfString&, bifstream*, XMLTag&);

////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

//! Reads XML header and root tag
/*!
  Check whether XML file has correct version tag and reads arts root
  tag information.

  \param is     Input stream
  \param ftype  File type
  \param ntype  Numeric type
  \param etype  Endian type
*/
void xml_read_header_from_stream(std::istream& is,
                                 FileType& ftype,
                                 NumericType& ntype,
                                 EndianType& etype);

//! Reads closing root tag
/*!
  Checks whether XML file ends correctly with \</arts\>.

  \param is  Input stream
*/
void xml_read_footer_from_stream(std::istream& is);

//! Writes XML header and root tag
/*!
  \param os     Output stream
  \param ftype  File type
*/
void xml_write_header_to_stream(std::ostream& os, FileType ftype);

//! Write closing root tag
/*!
  \param os Output stream
*/
void xml_write_footer_to_stream(std::ostream& os);

//! Open file for XML input
/*!
  This function opens an XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(std::ifstream& ifs, const String& name);

//! Open file for XML output
/*!
  This function opens an XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(std::ofstream& file, const String& name);

#ifdef ENABLE_ZLIB

//! Open file for zipped XML input
/*!
  This function opens a zipped XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(igzstream& ifs, const String& name);

//! Open file for zipped XML output
/*!
  This function opens a zipped XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(ogzstream& file, const String& name);

#endif  // ENABLE_ZLIB

////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

//! Reads data from XML file
/*!
  This is a generic functions that is used to read the XML header and
  footer info and calls the overloaded functions to read the data.

  \param filename XML filename
  \param type Generic return value
*/
template <arts_xml_ioable T>
void xml_read_from_file_base(const String& filename, T& type) {
  // Open input stream:
  std::unique_ptr<std::istream> ifs;
  if (filename.size() > 2 && filename.substr(filename.length() - 3, 3) == ".gz")
#ifdef ENABLE_ZLIB
  {
    ifs = std::make_unique<igzstream>();
    xml_open_input_file(*static_cast<igzstream*>(ifs.get()), filename);
  }
#else
  {
    throw std::runtime_error(
        "This arts version was compiled without zlib support.\n"
        "Thus zipped xml files cannot be read.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ifs = std::make_unique<std::ifstream>();
    xml_open_input_file(*static_cast<std::ifstream*>(ifs.get()), filename);
  }

  // Read the file into memory first to significantly speed up
  // the parsing (13x to 18x faster).
  std::stringstream buffer;
  {
    ARTS_NAMED_TIME_REPORT("XmlBuffering")
    buffer << ifs->rdbuf();
  }
  // No need to check for error, because xml_open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the matrix from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the
  // filename.
  try {
    ARTS_NAMED_TIME_REPORT("XmlStreaming")
    FileType ftype;
    NumericType ntype;
    EndianType etype;

    xml_read_header_from_stream(buffer, ftype, ntype, etype);
    if (ftype == FileType::ascii) {
      xml_io_stream<T>::read(buffer, type, static_cast<bifstream*>(nullptr));
    } else {
      String bfilename = filename + ".bin";
      bifstream bifs(bfilename.c_str());
      xml_io_stream<T>::read(buffer, type, &bifs);
    }
    xml_read_footer_from_stream(buffer);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error reading file: " << filename << '\n' << e.what();
    throw std::runtime_error(os.str());
  }
}

//! Write data to XML file
/*!
  This is a generic functions that is used to write the XML header and
  footer info and calls the overloaded functions to write the data.

  \param filename   XML filename
  \param type       Generic input value
  \param no_clobber 0: Overwrite, 1: Use unique filename
  \param ftype      File type
*/
template <arts_xml_ioable T>
void xml_write_to_file_base(const String& filename,
                            const T& type,
                            const FileType ftype) {
  std::unique_ptr<std::ostream> ofs;

  if (ftype == FileType::zascii)
#ifdef ENABLE_ZLIB
  {
    ofs = std::make_unique<ogzstream>();
    xml_open_output_file(*static_cast<ogzstream*>(ofs.get()), filename);
  }
#else
  {
    throw std::runtime_error(
        "This arts version was compiled without zlib support.\n"
        "Thus zipped xml files cannot be written.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ofs = std::make_unique<std::ofstream>();
    xml_open_output_file(*static_cast<std::ofstream*>(ofs.get()), filename);
  }

  try {
    xml_write_header_to_stream(*ofs, ftype);
    if (ftype == FileType::ascii or ftype == FileType::zascii) {
      xml_io_stream<T>::write(*ofs, type, static_cast<bofstream*>(nullptr), "");
    } else {
      String bfilename = filename + ".bin";
      bofstream bofs(bfilename.c_str());
      xml_io_stream<T>::write(*ofs, type, &bofs, "");
    }

    xml_write_footer_to_stream(*ofs);
  } catch (const std::runtime_error& e) {
    std::ostringstream os;
    os << "Error writing file: " << filename << '\n' << e.what();
    throw std::runtime_error(os.str());
  }
}

template <arts_xml_ioable T>
void xml_read_from_stream(std::istream& i, T& v, bifstream* b = nullptr) try {
  xml_io_stream<T>::read(i, v, b);
} catch (const std::exception& e) {
  throw std::runtime_error(std::format(
      "Error streaming {}:\n{}", xml_io_stream<T>::type_name, e.what()));
}

template <arts_xml_ioable T>
void xml_write_to_stream(std::ostream& o,
                         const T& v,
                         bofstream* b       = nullptr,
                         std::string_view n = ""sv) {
  xml_io_stream<T>::write(o, v, b, n);
}
