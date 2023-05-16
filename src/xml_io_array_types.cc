////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   xml_io_array_types.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2003-06-11

  \brief This file contains basic functions to handle XML data files.

*/

#include "arts.h"
#include "xml_io.h"

////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== ArrayOfAgenda ===========================================

//! Reads ArrayOfAgenda from XML input stream
/*!
  \param is_xml   XML Input stream
  \param aa       ArrayOfAgenda return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml _U_,
                          ArrayOfAgenda& aa _U_,
                          bifstream* pbifs _U_,
                          const Verbosity&) {
  throw runtime_error("Not supported.");
}

//! Writes ArrayOfAgenda to XML output stream
/*!
  \param os_xml   XML Output stream
  \param aa       ArrayOfAgenda
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml _U_,
                         const ArrayOfAgenda& aa _U_,
                         bofstream* pbofs _U_,
                         const String& name _U_,
                         const Verbosity&)

{
  throw runtime_error("ArrayOfAgendas can't be saved.");
}

//=== ArrayOfArrayOfSpeciesTag ================================================

//! Reads ArrayOfArrayOfSpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aastag  ArrayOfArrayOfSpeciesTag return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfSpeciesTag& aastag,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfSpeciesTag");

  tag.get_attribute_value("nelem", nelem);
  aastag.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) {
      xml_read_from_stream(is_xml, aastag[n], pbifs, verbosity);
    }
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfSpeciesTag: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfSpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aastag  ArrayOfArrayOfSpeciesTag
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfSpeciesTag& aastag,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfSpeciesTag");
  open_tag.add_attribute("nelem", aastag.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aastag.nelem(); n++) {
    xml_write_to_stream(os_xml, aastag[n], pbofs, "", verbosity);
  }

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfString ==========================================================

//! Reads ArrayOfArrayOfString from XML input stream
/*!
  \param is_xml   XML Input stream
  \param aastring  ArrayOfArrayOfString return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfString& aastring,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");

  tag.check_attribute("type", "ArrayOfString");

  tag.get_attribute_value("nelem", nelem);
  aastring.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aastring[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfString: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfString to XML output stream
/*!
  \param os_xml   XML Output stream
  \param aastring  ArrayOfArrayOfString
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfString& aastring,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfString");
  open_tag.add_attribute("nelem", aastring.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aastring.nelem(); n++)
    xml_write_to_stream(os_xml, aastring[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfPpath =========================================================

//! Reads ArrayOfPpath from XML input stream
/*!
  \param is_xml  XML Input stream
  \param appath  ArrayOfPpath return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfPpath& appath,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Ppath");

  tag.get_attribute_value("nelem", nelem);
  appath.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) {
      xml_read_from_stream(is_xml, appath[n], pbifs, verbosity);
    }
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfPpath: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfPpath to XML output stream
/*!
  \param os_xml  XML Output stream
  \param appath   ArrayOfPpath
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfPpath& appath,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Ppath");
  open_tag.add_attribute("nelem", appath.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < appath.nelem(); n++) {
    xml_write_to_stream(os_xml, appath[n], pbofs, "", verbosity);
  }

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfTensor3==================================================

//! Reads ArrayOfArrayOfTensor3 from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aatensor3  ArrayOfArrayOfTensor3 return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfTensor3& aatensor3,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfTensor3");

  tag.get_attribute_value("nelem", nelem);
  aatensor3.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aatensor3[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfTensor3: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfTensor3 to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aatensor3  ArrayOfArrayOfTensor3
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfTensor3& aatensor3,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfTensor3");
  open_tag.add_attribute("nelem", aatensor3.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aatensor3.nelem(); n++) {
    xml_write_to_stream(os_xml, aatensor3[n], pbofs, "", verbosity);
  }

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfTensor6==================================================

//! Reads ArrayOfArrayOfTensor6 from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aatensor6  ArrayOfArrayOfTensor6 return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfTensor6& aatensor6,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfTensor6");

  tag.get_attribute_value("nelem", nelem);
  aatensor6.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aatensor6[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfTensor6: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfTensor6 to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aatensor6  ArrayOfArrayOfTensor6
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfTensor6& aatensor6,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfTensor6");
  open_tag.add_attribute("nelem", aatensor6.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aatensor6.nelem(); n++)
    xml_write_to_stream(os_xml, aatensor6[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfGridPos =========================================================

//! Reads ArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param agpos   ArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfGridPos& agpos,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GridPos");

  tag.get_attribute_value("nelem", nelem);
  agpos.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agpos[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfGridPos: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param agpos   ArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfGridPos& agpos,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GridPos");
  open_tag.add_attribute("nelem", agpos.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agpos.nelem(); n++)
    xml_write_to_stream(os_xml, agpos[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfGridPos =====================================

//! Reads ArrayOfArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aagpos  ArrayOfArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfGridPos& aagpos,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfGridPos");

  tag.get_attribute_value("nelem", nelem);
  aagpos.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aagpos[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfGridPos: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aagpos  ArrayOfArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfGridPos& aagpos,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfGridPos");
  open_tag.add_attribute("nelem", aagpos.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagpos.nelem(); n++)
    xml_write_to_stream(os_xml, aagpos[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfArrayOfGridPos =====================================

//! Reads ArrayOfArrayOfArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aaagpos ArrayOfArrayOfArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfArrayOfGridPos& aaagpos,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfArrayOfGridPos");

  tag.get_attribute_value("nelem", nelem);
  aaagpos.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aaagpos[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfArrayOfGridPos: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aaagpos ArrayOfArrayOfArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfArrayOfGridPos& aaagpos,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfArrayOfGridPos");
  open_tag.add_attribute("nelem", aaagpos.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaagpos.nelem(); n++)
    xml_write_to_stream(os_xml, aaagpos[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfArrayOfArrayOfGridPos =====================================

//! Reads ArrayOfArrayOfArrayOfArrayOfGridPos from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aaaagpos ArrayOfArrayOfArrayOfArrayOfGridPos return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfArrayOfArrayOfGridPos& aaaagpos,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfArrayOfArrayOfGridPos");

  tag.get_attribute_value("nelem", nelem);
  aaaagpos.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aaaagpos[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfArrayOfArrayOfGridPos: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfArrayOfArrayOfGridPos to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aaaagpos   ArrayOfArrayOfArrayOfArrayOfGridPos
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfArrayOfArrayOfGridPos& aaaagpos,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfArrayOfArrayOfGridPos");
  open_tag.add_attribute("nelem", aaaagpos.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaaagpos.nelem(); n++)
    xml_write_to_stream(os_xml, aaaagpos[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfJacobianTarget ===========================================================

//! Reads ArrayOfJacobianTarget from XML input stream
/*!
  \param is_xml  XML Input stream
  \param ajt     ArrayOfJacobianTarget return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfJacobianTarget& ajt,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "JacobianTarget");

  tag.get_attribute_value("nelem", nelem);
  ajt.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, ajt[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfIndex: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfJacobianTarget to XML output stream
/*!
  \param os_xml  XML Output stream
  \param ajt     ArrayOfJacobianTarget
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfJacobianTarget& ajt,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "JacobianTarget");
  open_tag.add_attribute("nelem", ajt.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < ajt.nelem(); n++)
    xml_write_to_stream(os_xml, ajt[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfIndex ===========================================================

//! Reads ArrayOfIndex from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aindex  ArrayOfIndex return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfIndex& aindex,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Index");

  tag.get_attribute_value("nelem", nelem);
  aindex.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aindex[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfIndex: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfIndex to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aindex  ArrayOfIndex
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfIndex& aindex,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Index");
  open_tag.add_attribute("nelem", aindex.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aindex.nelem(); n++)
    xml_write_to_stream(os_xml, aindex[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfIndex =====================================================

//! Reads ArrayOfArrayOfIndex from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aaindex    ArrayOfArrayOfIndex return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfIndex& aaindex,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfIndex");

  tag.get_attribute_value("nelem", nelem);
  aaindex.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aaindex[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfIndex: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfIndex to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aaindex    ArrayOfArrayOfIndex
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfIndex& aaindex,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfIndex");
  open_tag.add_attribute("nelem", aaindex.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaindex.nelem(); n++)
    xml_write_to_stream(os_xml, aaindex[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfCIARecord ===========================================

//! Reads ArrayOfCIARecord from XML input stream
/*!
  \param is_xml   XML Input stream
  \param alspec   ArrayOfCIARecord return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfCIARecord& acr,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "CIARecord");

  tag.get_attribute_value("nelem", nelem);
  acr.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, acr[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfCIARecord: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfCIARecord to XML output stream
/*!
  \param os_xml   XML Output stream
  \param alspec   ArrayOfCIARecord
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfCIARecord& acr,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity)

{
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "CIARecord");
  open_tag.add_attribute("nelem", acr.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < acr.nelem(); n++)
    xml_write_to_stream(os_xml, acr[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfMatrix ==========================================================

//! Reads ArrayOfMatrix from XML input stream
/*!
  \param is_xml   XML Input stream
  \param amatrix  ArrayOfMatrix return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfMatrix& amatrix,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Matrix");

  tag.get_attribute_value("nelem", nelem);
  amatrix.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, amatrix[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfMatrix: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfMatrix to XML output stream
/*!
  \param os_xml   XML Output stream
  \param amatrix  ArrayOfMatrix
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfMatrix& amatrix,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Matrix");
  open_tag.add_attribute("nelem", amatrix.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < amatrix.nelem(); n++)
    xml_write_to_stream(os_xml, amatrix[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfMatrix====================================================

//! Reads ArrayOfArrayOfMatrix from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aamatrix   ArrayOfArrayOfMatrix return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfMatrix& aamatrix,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfMatrix");

  tag.get_attribute_value("nelem", nelem);
  aamatrix.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aamatrix[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfMatrix: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfMatrix to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aamatrix   ArrayOfArrayOfMatrix
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfMatrix& aamatrix,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfMatrix");
  open_tag.add_attribute("nelem", aamatrix.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aamatrix.nelem(); n++)
    xml_write_to_stream(os_xml, aamatrix[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfVector====================================================

//! Reads ArrayOfArrayOfVector from XML input stream
/*!
  \param is_xml     XML Input stream
  \param aaVector   ArrayOfArrayOfVector return value
  \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfVector& aavector,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfVector");

  tag.get_attribute_value("nelem", nelem);
  aavector.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aavector[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfVector: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfVector to XML output stream
/*!
  \param os_xml     XML Output stream
  \param aaVector   ArrayOfArrayOfVector
  \param pbofs      Pointer to binary file stream. NULL for ASCII output.
  \param name       Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfVector& aaVector,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfVector");
  open_tag.add_attribute("nelem", aaVector.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aaVector.nelem(); n++)
    xml_write_to_stream(os_xml, aaVector[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfSparse ==========================================================

//! Reads ArrayOfSparse from XML input stream
/*!
 \param is_xml   XML Input stream
 \param asparse  ArrayOfSparse return value
 \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfSparse& asparse,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Sparse");

  tag.get_attribute_value("nelem", nelem);
  asparse.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, asparse[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfSparse: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfSparse to XML output stream
/*!
 \param os_xml   XML Output stream
 \param asparse  ArrayOfSparse
 \param pbofs    Pointer to binary file stream. NULL for ASCII output.
 \param name     Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfSparse& asparse,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Sparse");
  open_tag.add_attribute("nelem", asparse.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < asparse.nelem(); n++)
    xml_write_to_stream(os_xml, asparse[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfQuantumIdentifier ================================================

//! Reads ArrayOfQuantumIdentifier from XML input stream
/*!
  \param is_xml  XML Input stream
  \param aqtag   ArrayOfQuantumIdentifier return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfQuantumIdentifier& aqtag,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "QuantumIdentifier");

  tag.get_attribute_value("nelem", nelem);
  aqtag.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) {
      aqtag[n] = QuantumIdentifier();
      xml_read_from_stream(is_xml, aqtag[n], pbifs, verbosity);
    }
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfQuantumIdentifier: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfQuantumIdentifier to XML output stream
/*!
  \param os_xml  XML Output stream
  \param aqtag   ArrayOfQuantumIdentifier
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfQuantumIdentifier& aqtag,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "QuantumIdentifier");
  open_tag.add_attribute("nelem", aqtag.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aqtag.nelem(); n++)
    xml_write_to_stream(os_xml, aqtag[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfRetrievalQuantity =======================================

//! Reads ArrayOfRetrievalQuantity from XML input stream
/*!
  \param is_xml    XML Input stream
  \param arq       ArrayOfRetrievalQuantity return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfRetrievalQuantity& arq,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "RetrievalQuantity");

  tag.get_attribute_value("nelem", nelem);
  arq.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, arq[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfRetrievalQuantity: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfRetrivalQuantity to XML output stream
/*!
  \param os_xml    XML Output stream
  \param arq       ArrayOfRetrievalQuantity
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfRetrievalQuantity& arq,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "RetrievalQuantity");
  open_tag.add_attribute("nelem", arq.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < arq.nelem(); n++)
    xml_write_to_stream(os_xml, arq[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfSpeciesTag ================================================

//! Reads ArrayOfSpeciesTag from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astag   ArrayOfSpeciesTag return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfSpeciesTag& astag,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "SpeciesTag");

  tag.get_attribute_value("nelem", nelem);
  astag.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, astag[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfSpeciesTag: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfSpeciesTag to XML output stream
/*!
  \param os_xml  XML Output stream
  \param astag   ArrayOfSpeciesTag
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfSpeciesTag& astag,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "SpeciesTag");
  open_tag.add_attribute("nelem", astag.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astag.nelem(); n++)
    xml_write_to_stream(os_xml, astag[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfSingleScatteringData===========================================

//! Reads ArrayOfSingleScatteringData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param assdata  ArrayOfSingleScatteringData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfSingleScatteringData& assdata,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "SingleScatteringData");

  tag.get_attribute_value("nelem", nelem);
  assdata.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, assdata[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfSingleScatteringData: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfSingleScatteringData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param assdata  ArrayOfSingleScatteringData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfSingleScatteringData& assdata,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "SingleScatteringData");
  open_tag.add_attribute("nelem", assdata.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < assdata.nelem(); n++)
    xml_write_to_stream(os_xml, assdata[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfSingleScatteringData===========================================

//! Reads ArrayOfArrayOfSingleScatteringData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param assdata  ArrayOfArrayOfSingleScatteringData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfSingleScatteringData& assdata,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfSingleScatteringData");

  tag.get_attribute_value("nelem", nelem);
  assdata.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, assdata[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfSingleScatteringData: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfSingleScatteringData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param assdata  ArrayOfArrayOfSingleScatteringData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfSingleScatteringData& assdata,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfSingleScatteringData");
  open_tag.add_attribute("nelem", assdata.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < assdata.nelem(); n++)
    xml_write_to_stream(os_xml, assdata[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfScatteringMetaData===========================================

//! Reads ArrayOfScatteringMetaData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param asmdata  ArrayOfScatteringMetaData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfScatteringMetaData& asmdata,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ScatteringMetaData");

  tag.get_attribute_value("nelem", nelem);
  asmdata.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, asmdata[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfScatteringMetaData: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfScatteringMetaData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param asmdata  ArrayOfScatteringMetaData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfScatteringMetaData& asmdata,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ScatteringMetaData");
  open_tag.add_attribute("nelem", asmdata.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < asmdata.nelem(); n++)
    xml_write_to_stream(os_xml, asmdata[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfScatteringMetaData===========================================

//! Reads ArrayOfArrayOfScatteringMetaData from XML input stream
/*!
  \param is_xml   XML Input stream
  \param asmdata  ArrayOfArrayOfScatteringMetaData return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfScatteringMetaData& asmdata,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfScatteringMetaData");

  tag.get_attribute_value("nelem", nelem);
  asmdata.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, asmdata[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfScatteringMetaData: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfScatteringMetaData to XML output stream
/*!
  \param os_xml   XML Output stream
  \param asmdata  ArrayOfArrayOfScatteringMetaData
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfScatteringMetaData& asmdata,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfScatteringMetaData");
  open_tag.add_attribute("nelem", asmdata.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < asmdata.nelem(); n++)
    xml_write_to_stream(os_xml, asmdata[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfGriddedField1 ===========================================

//! Reads ArrayOfGriddedField1 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGriddedField1 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfGriddedField1& agfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GriddedField1");

  tag.get_attribute_value("nelem", nelem);
  agfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfGriddedField1: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfGriddedField1 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGriddedField1
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfGriddedField1& agfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GriddedField1");
  open_tag.add_attribute("nelem", agfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem(); n++)
    xml_write_to_stream(os_xml, agfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfGriddedField2 ===========================================

//! Reads ArrayOfGriddedField2 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGriddedField2 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfGriddedField2& agfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GriddedField2");

  tag.get_attribute_value("nelem", nelem);
  agfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfGriddedField2: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfGriddedField2 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGriddedField2
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfGriddedField2& agfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GriddedField2");
  open_tag.add_attribute("nelem", agfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem(); n++)
    xml_write_to_stream(os_xml, agfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfGriddedField3 ===========================================

//! Reads ArrayOfGriddedField3 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGriddedField3 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfGriddedField3& agfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GriddedField3");

  tag.get_attribute_value("nelem", nelem);
  agfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfGriddedField3: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfGriddedField3 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGriddedField3
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfGriddedField3& agfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GriddedField3");
  open_tag.add_attribute("nelem", agfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem(); n++)
    xml_write_to_stream(os_xml, agfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfGriddedField1 ===========================================

//! Reads ArrayOfArrayOfGriddedField1 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param aagfield  ArrayOfArrayOfGriddedField1 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfGriddedField1& aagfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfGriddedField1");

  tag.get_attribute_value("nelem", nelem);
  aagfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aagfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfGriddedField1: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfGriddedField1 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param aagfield  ArrayOfArrayOfGriddedField1
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfGriddedField1& aagfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfGriddedField1");
  open_tag.add_attribute("nelem", aagfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagfield.nelem(); n++)
    xml_write_to_stream(os_xml, aagfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfGriddedField2 ===========================================

//! Reads ArrayOfArrayOfGriddedField2 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param aagfield  ArrayOfArrayOfGriddedField2 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfGriddedField2& aagfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfGriddedField2");

  tag.get_attribute_value("nelem", nelem);
  aagfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aagfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfGriddedField2: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfGriddedField2 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param aagfield  ArrayOfArrayOfGriddedField2
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfGriddedField2& aagfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfGriddedField2");
  open_tag.add_attribute("nelem", aagfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagfield.nelem(); n++)
    xml_write_to_stream(os_xml, aagfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfGriddedField3 ===========================================

//! Reads ArrayOfArrayOfGriddedField3 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param aagfield  ArrayOfArrayOfGriddedField3 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfGriddedField3& aagfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfGriddedField3");

  tag.get_attribute_value("nelem", nelem);
  aagfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aagfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfGriddedField3: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfGriddedField3 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param aagfield  ArrayOfArrayOfGriddedField3
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfGriddedField3& aagfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfGriddedField3");
  open_tag.add_attribute("nelem", aagfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aagfield.nelem(); n++)
    xml_write_to_stream(os_xml, aagfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfGriddedField4 ===========================================

//! Reads ArrayOfGriddedField4 from XML input stream
/*!
  \param is_xml   XML Input stream
  \param agfield  ArrayOfGriddedField4 return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfGriddedField4& agfield,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "GriddedField4");

  tag.get_attribute_value("nelem", nelem);
  agfield.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, agfield[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfGriddedField4: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfGriddedField4 to XML output stream
/*!
  \param os_xml   XML Output stream
  \param agfield  ArrayOfGriddedField4
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfGriddedField4& agfield,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "GriddedField4");
  open_tag.add_attribute("nelem", agfield.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < agfield.nelem(); n++)
    xml_write_to_stream(os_xml, agfield[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTelsemAtlas =========================================================

//! Reads ArrayOfTelsemAtlas from XML input stream
/*!
  \param is_xml      XML Input stream
  \param arr_telsem  ArrayOfTelsemAtlas return value
  \param pbifs       Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTelsemAtlas& arr_telsem,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "TelsemAtlas");

  tag.get_attribute_value("nelem", nelem);
  arr_telsem.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, arr_telsem[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTelsemAtlas: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTelsemAtlas to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor3  ArrayOfTelsemAtlas
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTelsemAtlas& arr_telsem,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "TelsemAtlas");
  open_tag.add_attribute("nelem", arr_telsem.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < arr_telsem.nelem(); n++)
    xml_write_to_stream(os_xml, arr_telsem[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTensor3=========================================================

//! Reads ArrayOfTensor3 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor3  ArrayOfTensor3 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTensor3& atensor3,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Tensor3");

  tag.get_attribute_value("nelem", nelem);
  atensor3.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, atensor3[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTensor3: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTensor3 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor3  ArrayOfTensor3
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTensor3& atensor3,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Tensor3");
  open_tag.add_attribute("nelem", atensor3.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor3.nelem(); n++)
    xml_write_to_stream(os_xml, atensor3[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTensor4=========================================================

//! Reads ArrayOfTensor4 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor4  ArrayOfTensor4 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTensor4& atensor4,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Tensor4");

  tag.get_attribute_value("nelem", nelem);
  atensor4.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, atensor4[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTensor4: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTensor4 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor4  ArrayOfTensor4
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTensor4& atensor4,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Tensor4");
  open_tag.add_attribute("nelem", atensor4.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor4.nelem(); n++)
    xml_write_to_stream(os_xml, atensor4[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTensor5=========================================================

//! Reads ArrayOfTensor5 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor5  ArrayOfTensor5 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTensor5& atensor5,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Tensor5");

  tag.get_attribute_value("nelem", nelem);
  atensor5.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, atensor5[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTensor5: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTensor5 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor5  ArrayOfTensor5
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTensor5& atensor5,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Tensor5");
  open_tag.add_attribute("nelem", atensor5.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor5.nelem(); n++)
    xml_write_to_stream(os_xml, atensor5[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTensor6=========================================================

//! Reads ArrayOfTensor6 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor6  ArrayOfTensor6 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTensor6& atensor6,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Tensor6");

  tag.get_attribute_value("nelem", nelem);
  atensor6.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, atensor6[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTensor6: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTensor6 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor6  ArrayOfTensor6
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTensor6& atensor6,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Tensor6");
  open_tag.add_attribute("nelem", atensor6.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor6.nelem(); n++)
    xml_write_to_stream(os_xml, atensor6[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTensor7=========================================================

//! Reads ArrayOfTensor7 from XML input stream
/*!
  \param is_xml    XML Input stream
  \param atensor7  ArrayOfTensor6 return value
  \param pbifs     Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTensor7& atensor7,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Tensor7");

  tag.get_attribute_value("nelem", nelem);
  atensor7.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, atensor7[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTensor7: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTensor6 to XML output stream
/*!
  \param os_xml    XML Output stream
  \param atensor7  ArrayOfTensor7
  \param pbofs     Pointer to binary file stream. NULL for ASCII output.
  \param name      Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTensor7& atensor7,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Tensor7");
  open_tag.add_attribute("nelem", atensor7.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atensor7.nelem(); n++)
    xml_write_to_stream(os_xml, atensor7[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfSun =========================================================

//! Reads ArrayOfSun from XML input stream
/*!
  \param is_xml  XML Input stream
  \param astar   astar return value
  \param pbifs   Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfSun& astar,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Sun");

  tag.get_attribute_value("nelem", nelem);
  astar.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++) {
      xml_read_from_stream(is_xml, astar[n], pbifs, verbosity);
    }
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfSun: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfSun to XML output stream
/*!
  \param os_xml  XML Output stream
  \param astar   ArrayOfSun
  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
  \param name    Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfSun& astar,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Sun");
  open_tag.add_attribute("nelem", astar.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astar.nelem(); n++) {
    xml_write_to_stream(os_xml, astar[n], pbofs, "", verbosity);
  }

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfString ==========================================================

//! Parse ArrayOfString from XML input stream
/*!
  \param is_xml   XML Input stream
  \param astring  ArrayOfString return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
  \param tag      XML tag object
*/
void xml_parse_from_stream(istream& is_xml,
                           ArrayOfString& astring,
                           bifstream* pbifs,
                           XMLTag& tag,
                           const Verbosity& verbosity) {
  Index nelem;

  tag.check_attribute("type", "String");

  tag.get_attribute_value("nelem", nelem);
  astring.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, astring[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfString: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }
}

//! Reads ArrayOfString from XML input stream
/*!
  \param is_xml   XML Input stream
  \param astring  ArrayOfString return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfString& astring,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  XMLTag tag(verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("Array");

  xml_parse_from_stream(is_xml, astring, pbifs, tag, verbosity);

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfString to XML output stream
/*!
  \param os_xml   XML Output stream
  \param astring  ArrayOfString
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfString& astring,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "String");
  open_tag.add_attribute("nelem", astring.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < astring.nelem(); n++)
    xml_write_to_stream(os_xml, astring[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfVector ==========================================================

//! Reads ArrayOfVector from XML input stream
/*!
  \param is_xml   XML Input stream
  \param avector  ArrayOfVector return value
  \param pbifs    Pointer to binary input stream. NULL in case of ASCII file.
*/
void xml_read_from_stream(istream& is_xml,
                          ArrayOfVector& avector,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Vector");

  tag.get_attribute_value("nelem", nelem);
  avector.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, avector[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfVector: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfVector to XML output stream
/*!
  \param os_xml   XML Output stream
  \param avector  ArrayOfVector
  \param pbofs    Pointer to binary file stream. NULL for ASCII output.
  \param name     Optional name attribute
*/
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfVector& avector,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Vector");
  open_tag.add_attribute("nelem", avector.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < avector.nelem(); n++)
    xml_write_to_stream(os_xml, avector[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTransmissionMatrix ======================================================

//! Reads ArrayOfTransmissionMatrix from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param apm        ArrayOfTransmissionMatrix return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTransmissionMatrix& atm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "TransmissionMatrix");

  tag.get_attribute_value("nelem", nelem);
  atm.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, atm[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTransmissionMatrix: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTransmissionMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param apm     ArrayOfTransmissionMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTransmissionMatrix& atm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "TransmissionMatrix");
  open_tag.add_attribute("nelem", atm.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < atm.nelem(); n++)
    xml_write_to_stream(os_xml, atm[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfTransmissionMatrix ======================================================

//! Reads ArrayOfArrayOfTransmissionMatrix from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param aatm       ArrayOfArrayOfTransmissionMatrix return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfTransmissionMatrix& aatm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfTransmissionMatrix");

  tag.get_attribute_value("nelem", nelem);
  aatm.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aatm[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfTransmissionMatrix: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfTransmissionMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param aatm    ArrayOfArrayOfTransmissionMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfTransmissionMatrix& aatm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfTransmissionMatrix");
  open_tag.add_attribute("nelem", aatm.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aatm.nelem(); n++)
    xml_write_to_stream(os_xml, aatm[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfRadiationVector ======================================================

//! Reads ArrayOfRadiationVector from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param apm        ArrayOfRadiationVector return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfRadiationVector& arv,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "RadiationVector");

  tag.get_attribute_value("nelem", nelem);
  arv.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, arv[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfRadiationVector: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfRadiationVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param arv     ArrayOfRadiationVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfRadiationVector& arv,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "RadiationVector");
  open_tag.add_attribute("nelem", arv.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < arv.nelem(); n++)
    xml_write_to_stream(os_xml, arv[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfRadiationVector ======================================================

//! Reads ArrayOfArrayOfRadiationVector from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param aaev       ArrayOfArrayOfRadiationVector return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfRadiationVector& aarv,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfRadiationVector");

  tag.get_attribute_value("nelem", nelem);
  aarv.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aarv[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfRadiationVector: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfRadiationVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param aarv    ArrayOfArrayOfRadiationVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfRadiationVector& aarv,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfRadiationVector");
  open_tag.add_attribute("nelem", aarv.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aarv.nelem(); n++)
    xml_write_to_stream(os_xml, aarv[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfPropagationMatrix ======================================================

//! Reads ArrayOfPropagationMatrix from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param apm        ArrayOfPropagationMatrix return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfPropagationMatrix& apm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "PropagationMatrix");

  tag.get_attribute_value("nelem", nelem);
  apm.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, apm[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfPropagationMatrix: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfPropagationMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param apm     ArrayOfPropagationMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfPropagationMatrix& apm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "PropagationMatrix");
  open_tag.add_attribute("nelem", apm.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < apm.nelem(); n++)
    xml_write_to_stream(os_xml, apm[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfAbsorptionLines ======================================================

//! Reads ArrayOfAbsorptionLines from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param aal        ArrayOfAbsorptionLines return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfAbsorptionLines& aal,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "AbsorptionLines");

  tag.get_attribute_value("nelem", nelem);
  aal.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aal[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfAbsorptionLines: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfAbsorptionLines to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param aal     ArrayOfAbsorptionLines
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfAbsorptionLines& aal,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "AbsorptionLines");
  open_tag.add_attribute("nelem", aal.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aal.nelem(); n++)
    xml_write_to_stream(os_xml, aal[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfAbsorptionLines ======================================================

//! Reads ArrayOfArrayOfAbsorptionLines from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param aal        ArrayOfArrayOfAbsorptionLines return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfAbsorptionLines& aal,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfAbsorptionLines");

  tag.get_attribute_value("nelem", nelem);
  aal.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aal[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfAbsorptionLines: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfAbsorptionLines to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param aal     ArrayOfArrayOfAbsorptionLines
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfAbsorptionLines& aal,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfAbsorptionLines");
  open_tag.add_attribute("nelem", aal.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aal.nelem(); n++)
    xml_write_to_stream(os_xml, aal[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfPropagationMatrix ======================================================

//! Reads ArrayOfArrayOfPropagationMatrix from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param aapm       ArrayOfArrayOfPropagationMatrix return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfPropagationMatrix& aapm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfPropagationMatrix");

  tag.get_attribute_value("nelem", nelem);
  aapm.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aapm[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfPropagationMatrix: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfPropagationMatrix to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param aapm    ArrayOfArrayOfPropagationMatrix
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfPropagationMatrix& aapm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfPropagationMatrix");
  open_tag.add_attribute("nelem", aapm.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aapm.nelem(); n++)
    xml_write_to_stream(os_xml, aapm[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfStokesVector ======================================================

//! Reads ArrayOfStokesVector from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param apm        ArrayOfStokesVector return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfStokesVector& apm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "StokesVector");

  tag.get_attribute_value("nelem", nelem);
  apm.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, apm[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfStokesVector: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfStokesVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param apm     ArrayOfStokesVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfStokesVector& apm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "StokesVector");
  open_tag.add_attribute("nelem", apm.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < apm.nelem(); n++)
    xml_write_to_stream(os_xml, apm[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfStokesVector ======================================================

//! Reads ArrayOfArrayOfStokesVector from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param aapm       ArrayOfArrayOfStokesVector return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfStokesVector& aapm,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfStokesVector");

  tag.get_attribute_value("nelem", nelem);
  aapm.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, aapm[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfStokesVector: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfStokesVector to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param aapm    ArrayOfArrayOfStokesVector
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfStokesVector& aapm,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfStokesVector");
  open_tag.add_attribute("nelem", aapm.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < aapm.nelem(); n++)
    xml_write_to_stream(os_xml, aapm[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfXsecRecord ======================================================

//! Reads ArrayOfXsecData from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param axd        ArrayOfXsecData return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfXsecRecord& axd,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "XsecRecord");

  tag.get_attribute_value("nelem", nelem);
  axd.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, axd[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfXsecRecord: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfXsecData to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param axd     ArrayOfXsecData
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfXsecRecord& axd,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "XsecRecord");
  open_tag.add_attribute("nelem", axd.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < axd.nelem(); n++)
    xml_write_to_stream(os_xml, axd[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfTime ======================================================

//! Reads ArrayOfTime from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param at         ArrayOfTime return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfTime& at,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "Time");
  
  tag.get_attribute_value("nelem", nelem);
  at.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, at[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfTime: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfTime to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param at      ArrayOfTime
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfTime& at,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "Time");
  open_tag.add_attribute("nelem", at.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < at.nelem(); n++)
    xml_write_to_stream(os_xml, at[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

//=== ArrayOfArrayOfTime ======================================================

//! Reads ArrayOfArrayOfTime from XML input stream
/*!
 * \param is_xml     XML Input stream
 * \param at         ArrayOfArrayOfTime return value
 * \param pbifs      Pointer to binary input stream. NULL in case of ASCII file.
 */
void xml_read_from_stream(istream& is_xml,
                          ArrayOfArrayOfTime& at,
                          bifstream* pbifs,
                          const Verbosity& verbosity) {
  ArtsXMLTag tag(verbosity);
  Index nelem;

  tag.read_from_stream(is_xml);
  tag.check_name("Array");
  tag.check_attribute("type", "ArrayOfTime");
  
  tag.get_attribute_value("nelem", nelem);
  at.resize(nelem);

  Index n;
  try {
    for (n = 0; n < nelem; n++)
      xml_read_from_stream(is_xml, at[n], pbifs, verbosity);
  } catch (const std::runtime_error& e) {
    ostringstream os;
    os << "Error reading ArrayOfArrayOfTime: "
       << "\n Element: " << n << "\n"
       << e.what();
    throw runtime_error(os.str());
  }

  tag.read_from_stream(is_xml);
  tag.check_name("/Array");
}

//! Writes ArrayOfArrayOfTime to XML output stream
/*!
 *  \param os_xml  XML Output stream
 *  \param at      ArrayOfArrayOfTime
 *  \param pbofs   Pointer to binary file stream. NULL for ASCII output.
 *  \param name    Optional name attribute
 */
void xml_write_to_stream(ostream& os_xml,
                         const ArrayOfArrayOfTime& at,
                         bofstream* pbofs,
                         const String& name,
                         const Verbosity& verbosity) {
  ArtsXMLTag open_tag(verbosity);
  ArtsXMLTag close_tag(verbosity);

  open_tag.set_name("Array");
  if (name.length()) open_tag.add_attribute("name", name);

  open_tag.add_attribute("type", "ArrayOfTime");
  open_tag.add_attribute("nelem", at.nelem());

  open_tag.write_to_stream(os_xml);
  os_xml << '\n';

  for (Index n = 0; n < at.nelem(); n++)
    xml_write_to_stream(os_xml, at[n], pbofs, "", verbosity);

  close_tag.set_name("/Array");
  close_tag.write_to_stream(os_xml);

  os_xml << '\n';
}

