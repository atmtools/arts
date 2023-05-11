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
#include "xml_io_base.h"
#include "xml_io_arts_types.h"
#include <type_traits>

////////////////////////////////////////////////////////////////////////////
//   XML parser classes
////////////////////////////////////////////////////////////////////////////

//! The ARTS XML tag class
/*!
  Handles reading, writing and constructing of XML tags.
*/
class ArtsXMLTag : public XMLTag {
 public:
  ArtsXMLTag() {};

  using XMLTag::add_attribute;
  
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

  using XMLTag::get_attribute_value;

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
};

////////////////////////////////////////////////////////////////////////////
//   General XML handling routines
////////////////////////////////////////////////////////////////////////////

void xml_parse_from_stream(
    istream &, Vector &, bifstream *, ArtsXMLTag &);

void xml_parse_from_stream(
    istream &, ArrayOfString &, bifstream *, ArtsXMLTag &);

////////////////////////////////////////////////////////////////////////////
//   Generic IO routines for XML files
////////////////////////////////////////////////////////////////////////////

/** Open plain or zipped xml file.
 *
 * Searches the include and data paths for the given filename.
 *
 * \param[out] ifs Pointer to input file stream
 * \param[in]  filename Input filename
 */
void xml_find_and_open_input_file(std::shared_ptr<istream>& ifs,
                                  const String& filename);

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_xml(String& filename, const String& varname);

void filename_xml_with_index(String& filename,
                             const Index& file_index,
                             const String& varname,
                             const Index& digits = 0);

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
                        T& type) requires (std::same_as<T, std::remove_const_t<T>>) {
  String xml_file = filename;
  find_xml_file(xml_file);
  xml_read_from_file_base(xml_file, type);
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
                       const Index no_clobber) {
  String efilename{add_basedir(filename)};

  std::unique_ptr<ostream> ofs;

  if (no_clobber) efilename = make_filename_unique(efilename, ".xml");

  xml_write_to_file_base(efilename, type, ftype);
}

#endif
