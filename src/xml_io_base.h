////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#ifndef xml_io_base_h
#define xml_io_base_h

#include <memory>

#include "arts.h"
#include "xml_io_general_types.h"

#ifdef ENABLE_ZLIB
#include "gzstream.h"
#endif

enum FileType : Index {
  FILE_TYPE_ASCII = 0,
  FILE_TYPE_ZIPPED_ASCII = 1,
  FILE_TYPE_BINARY = 2
};

enum NumericType { NUMERIC_TYPE_FLOAT, NUMERIC_TYPE_DOUBLE };
enum EndianType { ENDIAN_TYPE_LITTLE, ENDIAN_TYPE_BIG };

////////////////////////////////////////////////////////////////////////////
//   XML parser classes
////////////////////////////////////////////////////////////////////////////

//! XML attribute class
/*!
  Holds the name and value of an XML attribute.
*/

class XMLAttribute {
 public:
  String name;  /*!< Attribute name */
  String value; /*!< Attribute value */
};

//! The ARTS XML tag class
/*!
  Handles reading, writing and constructing of XML tags.
*/
class XMLTag {
 public:
  XMLTag(const Verbosity& rverbosity) : verbosity(rverbosity){};

  String& get_name() { return name; }

  void check_name(const String& expected_name);

  void set_name(const String& new_name) { name = new_name; }

  void add_attribute(const String& aname, String value);

  void add_attribute(const String& aname, const Index& value);
  
  /** Adds value of attribute as type Numeric to tag
   * 
   * @param[in] aname Attribute name
   * @param[in] value Set value
   */
  void add_attribute(const String& aname, const Numeric& value);
  
  void check_attribute(const String& aname, const String& value);

  void get_attribute_value(const String& aname, String& value);
  
  void get_attribute_value(const String& aname, Index& value);
  
  /** Returns value of attribute as type Numeric
   * 
   * Searches for the matching attribute and returns it value. If no
   * attribute with the given name exists, return value is set to
   * -1e99.
   * 
   * @param[in] aname Attribute name
   * @param[out] value Return value
   */
  void get_attribute_value(const String& aname, Numeric& value);
  
  void read_from_stream(istream& is);

  void write_to_stream(ostream& os);
  
  /** Returns if the attribute exists or not
   * 
   * @param[in] aname Attribute name
   * @return bool Does this attribute exist?
   */
  [[nodiscard]] bool has_attribute(const String& aname) const;

 protected:
  String name;                 /*!< Tag name */
  Array<XMLAttribute> attribs; /*!< List of attributes */
  const Verbosity& verbosity;
};

////////////////////////////////////////////////////////////////////////////
//   General XML handling routines
////////////////////////////////////////////////////////////////////////////

void xml_parse_error(const String& str_error);

void xml_data_parse_error(XMLTag& tag, const String& str_error);

void xml_set_stream_precision(ostream& os);

void parse_xml_tag_content_as_string(std::istream& is_xml, String& content);

void xml_parse_from_stream(
    istream &, Vector &, bifstream *, XMLTag &, const Verbosity &verbosity);

void xml_parse_from_stream(
    istream &, ArrayOfString &, bifstream *, XMLTag &, const Verbosity &);

FileType string2filetype(const String& file_format);

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
void xml_read_header_from_stream(istream& is,
                                 FileType& ftype,
                                 NumericType& ntype,
                                 EndianType& etype,
                                 const Verbosity& verbosity);

//! Reads closing root tag
/*!
  Checks whether XML file ends correctly with \</arts\>.

  \param is  Input stream
*/
void xml_read_footer_from_stream(istream& is, const Verbosity& verbosity);

//! Writes XML header and root tag
/*!
  \param os     Output stream
  \param ftype  File type
*/
void xml_write_header_to_stream(ostream& os,
                                FileType ftype,
                                const Verbosity& verbosity);

//! Write closing root tag
/*!
  \param os Output stream
*/
void xml_write_footer_to_stream(ostream& os, const Verbosity& verbosity);

//! Open file for XML input
/*!
  This function opens an XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(ifstream& ifs,
                         const String& name,
                         const Verbosity& verbosity);

//! Open file for XML output
/*!
  This function opens an XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(ofstream& file, const String& name);

#ifdef ENABLE_ZLIB

//! Open file for zipped XML input
/*!
  This function opens a zipped XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(igzstream& ifs,
                         const String& name,
                         const Verbosity& verbosity);

//! Open file for zipped XML output
/*!
  This function opens a zipped XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(ogzstream& file, const String& name);

#endif // ENABLE_ZLIB

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
template <typename T>
void xml_read_from_file_base(const String& filename,
                        T& type,
                        const Verbosity& verbosity) {
  CREATE_OUT2;

  out2 << "  Reading " + filename + '\n';

  // Open input stream:
  std::unique_ptr<istream> ifs;
  if (filename.nelem() > 2 &&
      filename.substr(filename.length() - 3, 3) == ".gz")
#ifdef ENABLE_ZLIB
  {
    ifs = std::make_unique<igzstream>();
    xml_open_input_file(
        *static_cast<igzstream*>(ifs.get()), filename, verbosity);
  }
#else
  {
    throw runtime_error(
        "This arts version was compiled without zlib support.\n"
        "Thus zipped xml files cannot be read.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ifs = std::make_unique<ifstream>();
    xml_open_input_file(
        *static_cast<ifstream*>(ifs.get()), filename, verbosity);
  }

  // No need to check for error, because xml_open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the matrix from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the
  // filename.
  try {
    FileType ftype;
    NumericType ntype;
    EndianType etype;

    xml_read_header_from_stream(*ifs, ftype, ntype, etype, verbosity);
    if (ftype == FILE_TYPE_ASCII) {
      xml_read_from_stream(*ifs, type, NULL, verbosity);
    } else {
      String bfilename = filename + ".bin";
      bifstream bifs(bfilename.c_str());
      xml_read_from_stream(*ifs, type, &bifs, verbosity);
    }
    xml_read_footer_from_stream(*ifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading file: " << filename << '\n' << e.what();
    throw runtime_error(os.str());
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
template <typename T>
void xml_write_to_file_base(const String& filename,
                            const T& type,
                            const FileType ftype,
                            const Verbosity& verbosity) {
  CREATE_OUT2;

  std::unique_ptr<ostream> ofs;

  out2 << "  Writing " << filename << '\n';
  if (ftype == FILE_TYPE_ZIPPED_ASCII)
#ifdef ENABLE_ZLIB
  {
    ofs = std::make_unique<ogzstream>();
    xml_open_output_file(*static_cast<ogzstream*>(ofs.get()), filename);
  }
#else
  {
    throw runtime_error(
        "This arts version was compiled without zlib support.\n"
        "Thus zipped xml files cannot be written.");
  }
#endif /* ENABLE_ZLIB */
  else {
    ofs = std::make_unique<ofstream>();
    xml_open_output_file(*static_cast<ofstream*>(ofs.get()), filename);
  }

  try {
    xml_write_header_to_stream(*ofs, ftype, verbosity);
    if (ftype == FILE_TYPE_ASCII || ftype == FILE_TYPE_ZIPPED_ASCII) {
      xml_write_to_stream(*ofs, type, NULL, "", verbosity);
    } else {
      String bfilename = filename + ".bin";
      bofstream bofs(bfilename.c_str());
      xml_write_to_stream(*ofs, type, &bofs, "", verbosity);
    }

    xml_write_footer_to_stream(*ofs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error writing file: " << filename << '\n' << e.what();
    throw runtime_error(os.str());
  }
}

#endif
