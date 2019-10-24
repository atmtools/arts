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
  \file   xml_io_private.h
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-11-06

  \brief This file contains private function declarations and
         template instantiation to handle XML data files.

*/

#ifndef xml_io_private_h
#define xml_io_private_h

#include <cfloat>
#include <memory>
#include <stdexcept>
#include "absorption.h"
#include "agenda_class.h"
#include "array.h"
#include "bifstream.h"
#include "bofstream.h"
#include "gas_abs_lookup.h"
#include "matpackVII.h"
#include "messages.h"
#include "optproperties.h"
#include "ppath.h"
#include "xml_io.h"

////////////////////////////////////////////////////////////////////////////
//   Functions to open and read XML files
////////////////////////////////////////////////////////////////////////////

void xml_open_output_file(ostream& file, const String& name);

void xml_open_input_file(ifstream& file,
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
  void add_attribute(const String& aname, const ArrayOfSpeciesTag& value, const bool self, const bool bath);

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
  void get_attribute_value(const String& aname, ArrayOfSpeciesTag& value, bool& self, bool& bath);
  
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

void xml_read_header_from_stream(istream& is,
                                 FileType& ftype,
                                 NumericType& ntype,
                                 EndianType& etype,
                                 const Verbosity& verbosity);

void xml_read_footer_from_stream(istream& is, const Verbosity& verbosity);

void xml_write_header_to_stream(ostream& os,
                                FileType ftype,
                                const Verbosity& verbosity);

void xml_write_footer_to_stream(ostream& os, const Verbosity& verbosity);

void xml_set_stream_precision(ostream& os);

void parse_xml_tag_content_as_string(std::istream& is_xml, String& content);

#endif /* xml_io_private_h */
