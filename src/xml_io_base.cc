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
  \file   xml_io.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/

#include "xml_io_base.h"
#include "arts.h"
#include "bifstream.h"
#include "bofstream.h"
#include "file.h"


////////////////////////////////////////////////////////////////////////////
//   XMLTag implementation
////////////////////////////////////////////////////////////////////////////

//! Check tag name
/*!
  Checks whether the name of the tag is correct. Throws runtime
  error otherwise.

  \param expected_name Expected tag name
*/
void XMLTag::check_name(const String& expected_name) {
  if (name != expected_name)
    xml_parse_error("Tag <" + expected_name + "> expected but <" + name +
                    "> found.");
}

//! Adds a String attribute to tag
/*!

  \param aname Attribute name
  \param value Attribute value
*/
void XMLTag::add_attribute(const String& aname, const String& value) {
  XMLAttribute attr;

  attr.name = aname;
  attr.value = value;
  attribs.push_back(attr);
}

//! Adds an Index attribute to tag
/*!

  \param aname Attribute name
  \param value Attribute value
*/
void XMLTag::add_attribute(const String& aname, const Index& value) {
  ostringstream v;

  v << value;
  add_attribute(aname, v.str());
}

void XMLTag::add_attribute(const String& aname, const Numeric& value) {
  ostringstream v;
  
  v << value;
  add_attribute(aname, v.str());
}

//! Checks whether attribute has the expected value
/*!

  If the attribute has another value or is unknown an exception is
  thrown.

  \param aname Attribute name
  \param value Expected value
*/
void XMLTag::check_attribute(const String& aname, const String& value) {
  String actual_value;

  get_attribute_value(aname, actual_value);

  if (actual_value == "*not found*") {
    xml_parse_error("Required attribute " + aname + " does not exist");
  } else if (actual_value != value) {
    xml_parse_error("Attribute " + aname + " has value \"" + actual_value +
                    "\" but \"" + value + "\" was expected.");
  }
}


bool XMLTag::has_attribute(const String& aname) const {
  return std::any_of(attribs.cbegin(), attribs.cend(), [&](auto& attr){return attr.name == aname;});
}

//! Returns value of attribute as String
/*!
  Searches for the matching attribute and returns it value. If no
  attribute with the given name exists, return value is set to
  *not found*.

  \param aname Attribute name
  \param value Return value
*/
void XMLTag::get_attribute_value(const String& aname, String& value) {
  value = "";

  Array<XMLAttribute>::iterator it = attribs.begin();
  while (it != attribs.end()) {
    if (it->name == aname) {
      value = it->value;
      it = attribs.end();
    } else {
      it++;
    }
  }
}

//! Returns value of attribute as type Index
/*!
  Searches for the matching attribute and returns it value. If no
  attribute with the given name exists, return value is set to
  *not found*.

  \param aname Attribute name
  \param value Return value
*/
void XMLTag::get_attribute_value(const String& aname, Index& value) {
  String attribute_value;
  istringstream strstr("");

  get_attribute_value(aname, attribute_value);
  strstr.str(attribute_value);
  strstr >> value;
  if (strstr.fail()) {
    xml_parse_error("Error while parsing value of " + aname + " from <" + name +
                    ">");
  }
}

void XMLTag::get_attribute_value(const String& aname, Numeric& value) {
  String attribute_value;
  istringstream strstr("");
  
  get_attribute_value(aname, attribute_value);
  strstr.str(attribute_value);
  strstr >> value;
  if (strstr.fail()) {
    xml_parse_error("Error while parsing value of " + aname + " from <" + name +
    ">");
  }
}

//! Reads next XML tag
/*!
  Reads the name and attributes of the next XML tag from stream.

  \param is Input stream
*/
void XMLTag::read_from_stream(istream& is) {
  CREATE_OUT3;

  String token;
  stringbuf tag;
  istringstream sstr("");
  XMLAttribute attr;
  char ch = 0;

  attribs.clear();

  while (is.good() && isspace(is.peek())) {
    is.get();
  }

  is >> ch;

  if (ch != '<') {
    is >> token;
    token = ch + token;

    xml_parse_error("'<' expected but " + token + " found.");
  }

  is.get(tag, '>');

  // Hit EOF while looking for '>'
  if (is.bad() || is.eof()) {
    xml_parse_error("Unexpected end of file while looking for '>'");
  }

  if (is.get() != '>') {
    xml_parse_error("Closing > not found in tag: " + tag.str());
  }

  sstr.str(tag.str() + '>');
  out3 << "Read: " << sstr.str() << '\n';

  sstr >> name;

  if (name[name.length() - 1] == '>') {
    // Because closin > was found, the tag for sure has no
    // attributes, set token to ">" to skip reading of attributes
    name.erase(name.length() - 1, 1);
    token = ">";
  } else {
    // Tag may have attributes, so read next token
    sstr >> token;
  }
  out3 << "Name: " << name << '\n';

  //extract attributes
  while (token != ">") {
    String::size_type pos;

    pos = token.find("=", 0);
    if (pos == String::npos) {
      xml_parse_error("Syntax error in tag: " + tag.str());
    }

    attr.name = token.substr(0, pos);
    token.erase(0, pos + 1);

    if (token[0] != '\"') {
      xml_parse_error("Missing \" in tag: " + tag.str());
    }

    while ((pos = token.find("\"", 1)) == (String::size_type)String::npos &&
           token != ">") {
      String ntoken;
      sstr >> ntoken;
      if (!ntoken.length()) break;
      token += " " + ntoken;
    }

    if (pos == (String::size_type)String::npos) {
      xml_parse_error("Missing \" in tag: " + sstr.str());
    }

    if (pos == 1)
      attr.value = "";
    else
      attr.value = token.substr(1, pos - 1);

    attribs.push_back(attr);

    out3 << "Attr: " << attr.name << '\n';
    out3 << "Value: " << attr.value << '\n';

    if (token[token.length() - 1] == '>') {
      token = ">";
    } else {
      sstr >> token;
    }
  }

  out3 << '\n';

  // Skip comments
  if (name == "comment") {
    is.get(tag, '<');

    // Hit EOF while looking for '<'
    if (is.bad() || is.eof()) {
      xml_parse_error(
          "Unexpected end of file while looking for "
          "comment tag");
    }

    read_from_stream(is);
    check_name("/comment");
    read_from_stream(is);
  }
}

//! Write XML tag
/*!
  Puts the tag together and writes it to stream.

  \param os Output stream
*/
void XMLTag::write_to_stream(ostream& os) {
  os << "<" << name;

  Array<XMLAttribute>::iterator it = attribs.begin();

  while (it != attribs.end()) {
    os << ' ' << it->name << "=\"" << it->value << '\"';
    it++;
  }

  os << ">";
}

FileType string2filetype(const String& file_format) {
  if (file_format == "ascii") return FILE_TYPE_ASCII;
  if (file_format == "zascii") return FILE_TYPE_ZIPPED_ASCII;
  if (file_format == "binary") return FILE_TYPE_BINARY;

  throw std::runtime_error(
      "file_format contains illegal string. "
      "Valid values are:\n"
      "  ascii:  XML output\n"
      "  zascii: Zipped XML output\n"
      "  binary: XML + binary output");
}

////////////////////////////////////////////////////////////////////////////
//   Functions to open and read XML files
////////////////////////////////////////////////////////////////////////////

//! Open file for XML output
/*!
  This function opens an XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(ofstream& file, const String& name) {
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(ios::badbit | ios::failbit);

  // c_str explicitly converts to c String.
  try {
    file.open(name.c_str());
  } catch (const std::exception&) {
    ostringstream os;
    os << "Cannot open output file: " << name << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly. (In that case we would not
  // get here if there really was a problem, because of the exception
  // thrown by open().)
  if (!file) {
    ostringstream os;
    os << "Cannot open output file: " << name << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw runtime_error(os.str());
  }
}

#ifdef ENABLE_ZLIB

//! Open file for zipped XML output
/*!
  This function opens a zipped XML file for writing.

  \param file Output filestream
  \param name Filename
*/
void xml_open_output_file(ogzstream& file, const String& name) {
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(ios::badbit | ios::failbit);

  // c_str explicitly converts to c String.
  String nname = name;

  if (nname.nelem() < 3 || nname.substr(nname.length() - 3, 3) != ".gz") {
    nname += ".gz";
  }

  try {
    file.open(nname.c_str());
  } catch (const ios::failure&) {
    ostringstream os;
    os << "Cannot open output file: " << nname << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly. (In that case we would not
  // get here if there really was a problem, because of the exception
  // thrown by open().)
  if (!file) {
    ostringstream os;
    os << "Cannot open output file: " << nname << '\n'
       << "Maybe you don't have write access "
       << "to the directory or the file?";
    throw runtime_error(os.str());
  }
}

#endif /* ENABLE_ZLIB */

//! Open file for XML input
/*!
  This function opens an XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(ifstream& ifs,
                         const String& name,
                         const Verbosity& verbosity) {
  CREATE_OUT3;

  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  ifs.exceptions(ios::badbit);

  // c_str explicitly converts to c String.
  try {
    ifs.open(name.c_str());
  } catch (const ios::failure&) {
    ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!ifs) {
    ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw runtime_error(os.str());
  }

  out3 << "- Reading input file " << name << "\n";
}

#ifdef ENABLE_ZLIB

//! Open file for zipped XML input
/*!
  This function opens a zipped XML file for reading.

  \param ifs   Input filestream
  \param name  Filename
*/
void xml_open_input_file(igzstream& ifs,
                         const String& name,
                         const Verbosity& verbosity) {
  CREATE_OUT3;

  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  ifs.exceptions(ios::badbit);

  // c_str explicitly converts to c String.
  try {
    ifs.open(name.c_str());
  } catch (const ios::failure&) {
    ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw runtime_error(os.str());
  }

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!ifs) {
    ostringstream os;
    os << "Cannot open input file: " << name << '\n'
       << "Maybe the file does not exist?";
    throw runtime_error(os.str());
  }

  out3 << "- Reading input file " << name << "\n";
}

#endif /* ENABLE_ZLIB */

////////////////////////////////////////////////////////////////////////////
//   General XML functions (file header, start root tag, end root tag)
////////////////////////////////////////////////////////////////////////////

//! Throws XML parser runtime error
/*!
  This is used quite often inside the parsing routines so it's a
  function for itself.

  \param str_error Error description
*/
void xml_parse_error(const String& str_error) {
  ostringstream os;
  os << "XML parse error: " << str_error << '\n'
     << "Check syntax of XML file\n";
  throw runtime_error(os.str());
}

//! Throws XML parser runtime error
/*!
  This is used quite often inside the data parsing routines so it's a
  function for itself.

  \param tag        XMLTag
  \param str_error  Error description
*/
void xml_data_parse_error(XMLTag& tag, String str_error) {
  ostringstream os;
  os << "XML data parse error: Error reading ";
  tag.write_to_stream(os);
  os << str_error << "\n"
     << "Check syntax of XML file. A possible cause is that the file "
     << "contains NaN or Inf values.\n";
  throw runtime_error(os.str());
}

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
                                 const Verbosity& verbosity) {
  char str[6];
  stringbuf strbuf;
  XMLTag tag(verbosity);
  String strtype;

  while (!is.fail() && isspace(is.peek())) is.get();

  is.get(str, 6, ' ');

  if (string(str) != "<?xml") {
    xml_parse_error(
        "Input file is not a valid xml file "
        "(<?xml not found)");
  }

  is.get(strbuf, '>');
  is.get();

  if (is.fail()) {
    xml_parse_error("Input file is not a valid xml file");
  }

  tag.read_from_stream(is);
  tag.check_name("arts");

  // Check file format
  tag.get_attribute_value("format", strtype);
  if (strtype == "binary") {
    ftype = FILE_TYPE_BINARY;
  } else {
    ftype = FILE_TYPE_ASCII;
  }

  // Check endian type
  tag.get_attribute_value("endian_type", strtype);
  if (strtype == "little") {
    etype = ENDIAN_TYPE_LITTLE;
  }
  if (strtype == "big") {
    etype = ENDIAN_TYPE_BIG;
  }
  if (strtype == "") {
    /*      out1 << "  Warning: Endian type not specified in XML file, "
              <<    "assuming little endian (PC)\n";*/
    etype = ENDIAN_TYPE_LITTLE;
  } else {
    ostringstream os;
    os << "  Error: Unknown endian type \"" << strtype
       << "\" specified in XML file.\n";
    throw runtime_error(os.str());
  }

  // Check numeric type
  tag.get_attribute_value("numeric_type", strtype);
  if (strtype == "float") {
    ntype = NUMERIC_TYPE_FLOAT;
  } else if (strtype == "double") {
    ntype = NUMERIC_TYPE_DOUBLE;
  } else if (strtype == "") {
    /*      out1 << "  Warning: Numeric type not specified in XML file, "
              <<    "assuming double\n";*/
    ntype = NUMERIC_TYPE_DOUBLE;
  } else {
    ostringstream os;
    os << "  Error: Unknown numeric type \"" << strtype
       << "\" specified in XML file.\n";
    throw runtime_error(os.str());
  }
}

//! Reads closing root tag
/*!
  Checks whether XML file ends correctly with \</arts\>.

  \param is  Input stream
*/
void xml_read_footer_from_stream(istream& is, const Verbosity& verbosity) {
  XMLTag tag(verbosity);

  tag.read_from_stream(is);
  tag.check_name("/arts");
}

//! Writes XML header and root tag
/*!
  \param os     Output stream
  \param ftype  File type
*/
void xml_write_header_to_stream(ostream& os,
                                FileType ftype,
                                const Verbosity& verbosity) {
  XMLTag tag(verbosity);

  os << "<?xml version=\"1.0\"?>" << '\n';

  tag.set_name("arts");
  switch (ftype) {
    case FILE_TYPE_ASCII:
    case FILE_TYPE_ZIPPED_ASCII:
      tag.add_attribute("format", "ascii");
      break;
    case FILE_TYPE_BINARY:
      tag.add_attribute("format", "binary");
      break;
  }

  tag.add_attribute("version", "1");

  tag.write_to_stream(os);

  os << '\n';
}

//! Write closing root tag
/*!
  \param os Output stream
*/
void xml_write_footer_to_stream(ostream& os, const Verbosity& verbosity) {
  XMLTag tag(verbosity);

  tag.set_name("/arts");
  tag.write_to_stream(os);

  os << endl;
}

void xml_set_stream_precision(ostream& os) {
  // Determine the precision, depending on whether Numeric is double
  // or float:
  int precision;
#ifdef USE_FLOAT
  precision = FLT_DIG;
#else
#ifdef USE_DOUBLE
  precision = DBL_DIG;
#else
#error Numeric must be double or float
#endif
#endif

  os << setprecision(precision);
}

//! Get the content of an xml tag as a string
void parse_xml_tag_content_as_string(std::istream& is_xml, String& content) {
  char dummy;

  content = "";
  dummy = (char)is_xml.peek();
  while (is_xml && dummy != '<') {
    is_xml.get(dummy);
    content += dummy;
    dummy = (char)is_xml.peek();
  }

  if (!is_xml) throw std::runtime_error("Unexpected end of file.");
}
