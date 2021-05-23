/* Copyright (C) 2002-2012 Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#ifndef xml_io_h
#define xml_io_h

#include "arts.h"
#include "xml_io_types.h"

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
class ArtsXMLTag {
 public:
  ArtsXMLTag(const Verbosity& rverbosity) : verbosity(rverbosity){};

  String& get_name() { return name; }

  void check_name(const String& expected_name);

  void set_name(const String& new_name) { name = new_name; }

  void add_attribute(const String& aname, const String& value);

  void add_attribute(const String& aname, const Index& value);
  
  /** Adds value of attribute as type Numeric to tag
   * 
   * @param[in] aname Attribute name
   * @param[in] value Set value
   */
  void add_attribute(const String& aname, const Numeric& value);
  
  /** Adds value of attribute as type std::vector<QuantumNumberType> to tag
   * 
   * @param[in] aname Attribute name
   * @param[in] value Set value
   */
  void add_attribute(const String& aname, const std::vector<QuantumNumberType>& value);
  
  /** Adds value of attribute
   * 
   * @param[in] aname Attribute name
   * @param[in] value SpeciesTag(s) for all lines.  Basic initialization at self and bath
   * @param[in] self True if LineShape::self_broadening in list
   * @param[in] bath True if LineShape::bath_broadening in list
   */
  void add_attribute(const String& aname, const ArrayOfSpecies& value, const bool self, const bool bath);

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
  
  /** Returns value of attribute as type SpeciesTag
   * 
   * Searches for the matching attribute and returns it value. If no
   * attribute with the given name exists, it fails exceptionally.
   * 
   * @param[in] aname Attribute name
   * @param[out] value Return value
   */
  void get_attribute_value(const String& aname, SpeciesTag& value);
  
  /** Returns value of attribute as type ArrayOfSpeciesTag
   * 
   * Searches for the matching attribute and returns it value. If no
   * attribute with the given name exists, it fails exceptionally.
   * 
   * @param[in] aname Attribute name
   * @param[out] value SpeciesTag(s) for all lines.  Basic initialization at self and bath
   * @param[out] self True if LineShape::self_broadening in list
   * @param[out] bath True if LineShape::bath_broadening in list
   */
  void get_attribute_value(const String& aname, ArrayOfSpecies& value, bool& self, bool& bath);
  
  /** Returns value of attribute as type ArrayOfSpeciesTag
   * 
   * Searches for the matching attribute and returns it value
   * 
   * @param[in] aname Attribute name
   * @param[out] value Return value
   */
  void get_attribute_value(const String& aname, std::vector<QuantumNumberType>& value);
  
  /** Returns value of attribute as type ArrayOfSpeciesTag
   * 
   * Searches for the matching attribute and returns it value
   * 
   * @param[in] aname Attribute name
   * @param[in,out] value Return value
   */
  void get_attribute_value(const String& aname, QuantumNumbers& value);

  void read_from_stream(istream& is);

  void write_to_stream(ostream& os);

 private:
  String name;                 /*!< Tag name */
  Array<XMLAttribute> attribs; /*!< List of attributes */
  const Verbosity& verbosity;
};

////////////////////////////////////////////////////////////////////////////
//   General XML handling routines
////////////////////////////////////////////////////////////////////////////

void xml_parse_error(const String& str_error);

void xml_data_parse_error(ArtsXMLTag& tag, String str_error);

void xml_set_stream_precision(ostream& os);

void parse_xml_tag_content_as_string(std::istream& is_xml, String& content);

void xml_parse_from_stream(
    istream &, Vector &, bifstream *, ArtsXMLTag &, const Verbosity &verbosity);

void xml_parse_from_stream(
    istream &, ArrayOfString &, bifstream *, ArtsXMLTag &, const Verbosity &);

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_xml(String& filename, const String& varname);

void filename_xml_with_index(String& filename,
                             const Index& file_index,
                             const String& varname,
                             const Index& digits = 0);

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

/** Open plain or zipped xml file.
 *
 * Searches the include and data paths for the given filename.
 *
 * \param[out] ifs Pointer to input file stream
 * \param[in]  filename Input filename
 * \param[in]  verbosity Verbosity
 */
void xml_find_and_open_input_file(std::shared_ptr<istream>& ifs,
                                  const String& filename,
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
void xml_read_from_file(const String& filename,
                        T& type,
                        const Verbosity& verbosity) {
  CREATE_OUT2;

  String xml_file = filename;
  find_xml_file(xml_file, verbosity);
  out2 << "  Reading " + xml_file + '\n';

  // Open input stream:
  std::unique_ptr<istream> ifs;
  if (xml_file.nelem() > 2 &&
      xml_file.substr(xml_file.length() - 3, 3) == ".gz")
#ifdef ENABLE_ZLIB
  {
    ifs = std::make_unique<igzstream>();
    xml_open_input_file(
        *static_cast<igzstream*>(ifs.get()), xml_file, verbosity);
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
        *static_cast<ifstream*>(ifs.get()), xml_file, verbosity);
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
      String bfilename = xml_file + ".bin";
      bifstream bifs(bfilename.c_str());
      xml_read_from_stream(*ifs, type, &bifs, verbosity);
    }
    xml_read_footer_from_stream(*ifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading file: " << xml_file << '\n' << e.what();
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
void xml_write_to_file(const String& filename,
                       const T& type,
                       const FileType ftype,
                       const Index no_clobber,
                       const Verbosity& verbosity) {
  CREATE_OUT2;

  String efilename = add_basedir(filename);

  std::unique_ptr<ostream> ofs;

  if (no_clobber) make_filename_unique(efilename, ".xml");

  out2 << "  Writing " << efilename << '\n';
  if (ftype == FILE_TYPE_ZIPPED_ASCII)
#ifdef ENABLE_ZLIB
  {
    ofs = std::make_unique<ogzstream>();
    xml_open_output_file(*static_cast<ogzstream*>(ofs.get()), efilename);
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
    xml_open_output_file(*static_cast<ofstream*>(ofs.get()), efilename);
  }

  try {
    xml_write_header_to_stream(*ofs, ftype, verbosity);
    if (ftype == FILE_TYPE_ASCII || ftype == FILE_TYPE_ZIPPED_ASCII) {
      xml_write_to_stream(*ofs, type, NULL, "", verbosity);
    } else {
      String bfilename = efilename + ".bin";
      bofstream bofs(bfilename.c_str());
      xml_write_to_stream(*ofs, type, &bofs, "", verbosity);
    }

    xml_write_footer_to_stream(*ofs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error writing file: " << efilename << '\n' << e.what();
    throw runtime_error(os.str());
  }
}

#endif
