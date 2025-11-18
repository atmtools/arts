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

#include <debug.h>
#include <file.h>
#include <time_report.h>

#include <type_traits>

#include "xml_io_base.h"

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
void xml_find_and_open_input_file(std::shared_ptr<std::istream>& ifs,
                                  const String& filename);

////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

void filename_xml(const String& filename);

void filename_xml_with_index(String& filename,
                             const Index& file_index,
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
template <arts_xml_ioable T>
String xml_read_from_file(const String& filename, T& type)
  requires(std::same_as<T, std::remove_const_t<T>>)
try {
  ARTS_TIME_REPORT

  String xml_file = filename;
  find_xml_file(xml_file);
  xml_read_from_file_base(xml_file, type);
  return xml_file;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading file {} containing {}:\n{}",
                  filename,
                  xml_io_stream_name_v<T>,
                  e.what()));
}

//! Extends data from XML file
/*!
  This is a generic functions that is used to read the XML header and
  footer info and calls the overloaded functions to extend the data.

  \param filename XML filename
  \param type Generic return value
*/
template <arts_xml_extendable T>
String xml_extend_from_file(const String& filename, T& type)
  requires(std::same_as<T, std::remove_const_t<T>>)
try {
  ARTS_TIME_REPORT

  String xml_file = filename;
  find_xml_file(xml_file);
  xml_extend_from_file_base(xml_file, type);
  return xml_file;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error extending file {} containing {}:\n{}",
                  filename,
                  xml_io_stream_name_v<T>,
                  e.what()));
}

//! Appends data from XML file
/*!
  This is a generic functions that is used to read the XML header and
  footer info and calls the overloaded functions to append the data.

  \param filename XML filename
  \param type Generic return value
*/
template <arts_xml_appendable T>
String xml_append_from_file(const String& filename, T& type)
  requires(std::same_as<T, std::remove_const_t<T>>)
try {
  ARTS_TIME_REPORT

  String xml_file = filename;
  find_xml_file(xml_file);
  xml_append_from_file_base(xml_file, type);
  return xml_file;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error appending file {} containing {}:\n{}",
                  filename,
                  xml_io_stream_name_v<T>,
                  e.what()));
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
String xml_write_to_file(const String& filename,
                         const T& type,
                         const FileType ftype,
                         const Index no_clobber) try {
  ARTS_TIME_REPORT

  String efilename{add_basedir(filename)};

  std::unique_ptr<std::ostream> ofs;

  if (no_clobber) efilename = make_filename_unique(efilename, ".xml");

  xml_write_to_file_base(efilename, type, ftype);

  return efilename;
}
ARTS_METHOD_ERROR_CATCH

String complete_basename(const String& basename);
