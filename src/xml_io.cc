/* Copyright (C) 2002 Oliver Lemke <olemke@uni-bremen.de>

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
  \author Oliver Lemke <olemke@uni-bremen.de>
  \date   2002-05-10

  \brief This file contains basic functions to handle XML data files.

*/


#include "xml_io.h"
#include "xml_io_private.h"


////////////////////////////////////////////////////////////////////////////
//   XML parser classes
////////////////////////////////////////////////////////////////////////////

//! XML attribute class
/*!
  Holds the name and value of an XML attribute.
*/

class XMLAttribute
{
public:
  String name;                  /*!< Attribute name */
  String value;                 /*!< Attribute value */
};


//! The ARTS XML tag class
/*!
  Handles reading, writing and constructing of XML tags.
*/
class ArtsXMLTag
{
public:

  String&
  get_name () { return (name); }

  void
  check_name (const String& expected_name);

  void
  set_name (const String& new_name) { name = new_name; }

  void
  add_attribute (const String& aname, const String& value);

  void
  add_attribute (const String& aname, const Index& value);

  void
  check_attribute (const String& aname, const String& value);

  void
  get_attribute_value (const String& aname, String& value);

  void
  get_attribute_value (const String& aname, Index& value);

  void
  read_from_stream (istream& is);

  void
  write_to_stream (ostream& os);

private:
  String name;                  /*!< Tag name */
  Array<XMLAttribute> attribs;  /*!< List of attributes */
};


////////////////////////////////////////////////////////////////////////////
//   ArtsXMLTag implementation
////////////////////////////////////////////////////////////////////////////

//! Check tag name
/*!
  Checks whether the name of the tag is correct. Throws runtime
  error otherwise.
  
  \param is Input stream
  \param expected_name Expected tag name
*/
void
ArtsXMLTag::check_name (const String& expected_name)
{
  if (name != expected_name)
    {
      xml_parse_error ("Tag <" + expected_name + "> expected but <"
                       + name + "> found.");
    }
  
}


//! Adds a String attribute to tag
/*! 
  
  \param aname Attribute name
  \param value Attribute value
*/
void
ArtsXMLTag::add_attribute (const String& aname, const String& value)
{
  XMLAttribute attr;

  attr.name = aname;
  attr.value = value;
  attribs.push_back (attr);
}


//! Adds an Index attribute to tag
/*! 
  
  \param aname Attribute name
  \param value Attribute value
*/
void
ArtsXMLTag::add_attribute (const String& aname, const Index& value)
{
  ostringstream v;

  v << value;
  add_attribute (aname, v.str ());
}


//! Checks whether attribute has the expected value
/*! 

  If the attribute has another value or is unknown an exception is
  thrown.

  \param aname Attribute name
  \param value Expected value
*/
void
ArtsXMLTag::check_attribute (const String& aname, const String& value)
{
  String actual_value;
  
  get_attribute_value (aname, actual_value);

  if (actual_value == "*not found*")
    {
      xml_parse_error ("Required attribute " + aname
                       + " does not exist");
    }
  else if (actual_value != value)
    {
      xml_parse_error ("Attribute " + aname + " has value "
                       + actual_value + " but "
                       + value + " was expected.");
    }
}


//! Returns value of attribute as String
/*! 
  Searches for the matching attribute and returns it value. If no
  attribute with the given name exists, return value is set to
  *not found*.

  \param aname Attribute name
  \param value Return value
*/
void
ArtsXMLTag::get_attribute_value (const String& aname, String& value)
{
  value = "*not found*";

  Array<XMLAttribute>::iterator it = attribs.begin ();
  while (it != attribs.end ())
    {
      if (it->name == aname)
        {
          value = it->value;
          it = attribs.end ();
        }
      else
        {
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
void
ArtsXMLTag::get_attribute_value (const String& aname, Index& value)
{
  String attribute_value;
  istringstream strstr ("");

  get_attribute_value (aname, attribute_value);
  strstr.str (attribute_value);
  strstr >> value;
  if (strstr.fail ())
    {
      xml_parse_error ("Error while parsing value of " + aname
                       + " from <" + name + ">");
    }
}



//! Reads next XML tag
/*!
  Reads the name and attributes of the next XML tag from stream.
  
  \param is Input stream
*/
void
ArtsXMLTag::read_from_stream (istream& is)
{
  String        token;
  stringbuf     tag;
  istringstream sstr ("");
  XMLAttribute  attr;
  char          ch;

  attribs.clear ();

  while (is.good () && isspace (is.peek ()))
    {
      is.get ();
    }

  is >> ch;

  if (ch != '<')
    {
      is >> token;
      token = ch + token;

      xml_parse_error ("'<' expected but " + token + " found.");
    }

  is.get (tag, '>');

  // Hit EOF while looking for '>'
  if (is.bad () || is.eof ())
    {
      xml_parse_error ("Unexpected end of file while looking for '>'");
    }

  if (is.get () != '>')
    {
      xml_parse_error ("Closing > not found in tag: " + tag.str ());
    }


  sstr.str (tag.str () + '>');
  out2 << "Read: " << sstr.str () << '\n';

  sstr >> name;

  if (name [name.length () - 1] == '>')
    {
      // Because closin > was found, the tag for sure has no
      // attributes, set token to ">" to skip reading of attributes
      name.erase (name.length () - 1, 1);
      token = ">";
    }
  else
    {
      // Tag may have attributes, so read next token
      sstr >> token;
    }
  out2 << "Name: " << name << '\n';

  //extract attributes
  while (token != ">")
  {
    String::size_type pos;

    pos = token.find ("=", 0);
    if (pos == String::npos)
      {
        xml_parse_error ("Syntax error in tag: " + tag.str ());
      }

    attr.name = token.substr (0, pos);
    token.erase (0, pos + 1);

    if (token[0] != '\"')
      {
        xml_parse_error ("Missing \" in tag: " + tag.str ());
      }

    pos = token.find ("\"", 1);
    if (pos == (String::size_type)String::npos)
      {
        xml_parse_error ("Missing \" in tag: " + sstr.str ());
      }

    attr.value = token.substr (1, pos - 1);

    attribs.push_back (attr);

    out2 << "Attr: " << attr.name << '\n';
    out2 << "Value: " << attr.value << '\n';

    if (token[token.length () - 1] == '>')
      {
        token = ">";
      }
    else
      {
        sstr >> token;
      }
  }

  out2 << '\n';

  // Skip comments
  if (name == "comment")
    {
      is.get (tag, '<');

      // Hit EOF while looking for '<'
      if (is.bad () || is.eof ())
        {
          xml_parse_error ("Unexpected end of file while looking for "
                           "comment tag");
        }

      read_from_stream (is);
      check_name ("/comment");
      read_from_stream (is);
    }
}


//! Write XML tag
/*!
  Puts the tag together and writes it to stream.
  
  \param os Output stream
*/
void
ArtsXMLTag::write_to_stream (ostream& os)
{

  os << "<" << name;

  Array<XMLAttribute>::iterator it = attribs.begin ();
  
  while (it != attribs.end ())
    {
      os << ' ' << it->name
         << "=\"" << it->value << '\"';
      it++;
    }

  os << ">";
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
void
filename_xml (String&  filename,
              const String&  varname )
{
  if ("" == filename)
    {
      extern const String out_basename;
      filename = out_basename + "." + varname + ".xml";
    }
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
void
xml_open_output_file (ofstream& file, const String& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(ios::badbit |
                  ios::failbit);
  
  // c_str explicitly converts to c String.
  file.open(name.c_str() );

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly. (In that case we would not
  // get here if there really was a problem, because of the exception
  // thrown by open().)
  if (!file)
    {
      ostringstream os;
      os << "Cannot open output file: " << name << '\n'
         << "Maybe you don't have write access "
         << "to the directory or the file?";
      throw runtime_error(os.str());
    }
}


//! Open file for XML input
/*! 
  This function opens an XML file for reading.

  \param file Input filestream
  \param name Filename
*/
void
xml_open_input_file (ifstream& ifs, const String& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  ifs.exceptions (ios::badbit);

  // c_str explicitly converts to c String.
  ifs.open (name.c_str ());

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!ifs)
    {
      ostringstream os;
      os << "Cannot open input file: " << name << '\n'
         << "Maybe the file does not exist?";
      throw runtime_error (os.str ());
    }
}



////////////////////////////////////////////////////////////////////////////
//   General XML functions (file header, start root tag, end root tag)
////////////////////////////////////////////////////////////////////////////

//! Throws XML parser runtime error
/*!
  This is used quite often inside the parsing routines so it's a
  function for itself.
  
  \param str_error Error description
*/
void
xml_parse_error (const String& str_error)
{
  ostringstream os;
  os << "XML parse error: " << str_error << '\n'
     << "Check syntax of XML file\n";
  throw runtime_error (os.str ());
}


//! Reads XML header and root tag
/*! 
  Check whether XML file has correct version tag and reads arts root
  tag information.
  
  \param is Input stream
*/
void
xml_read_header_from_stream (istream& is)
{
  char str[6];
  stringbuf strbuf;
  ArtsXMLTag tag;

  is.get (str, 6);

  if (strcasecmp (str, "<?xml"))
    {
      xml_parse_error ("Input file is not a valid xml file "
                       "(<?xml not found)");
    }

  is.get (strbuf, '>');
  is.get ();

  if (is.fail ())
    {
      xml_parse_error ("Input file is not a valid xml file");
    }

  tag.read_from_stream (is);
  tag.check_name ("arts");
}


//! Reads closing root tag
/*! 
  Checks whether XML file ends correctly with </arts>.
  
  \param is Input stream
*/
void
xml_read_footer_from_stream (istream& is)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is);
  tag.check_name ("/arts");
}


//! Writes XML header and root tag
/*! 
  \param os Output stream
*/
void
xml_write_header_to_stream (ostream& os)
{
  ArtsXMLTag tag;

  os << "<?xml version=\"1.0\"?>"
     << '\n';

  tag.set_name ("arts");
  tag.add_attribute ("format", "ascii");
  tag.add_attribute ("version", "1");
  tag.write_to_stream (os);

  os << '\n';

}


//! Write closing root tag
/*! 
  \param os Output stream
*/
void
xml_write_footer_to_stream (ostream& os)
{
  ArtsXMLTag tag;
  
  tag.set_name ("/arts");
  tag.write_to_stream (os);

  os << endl;
}

void
xml_set_stream_precision (ostream &os)
{
  // Determine the precision, depending on whether Numeric is double
  // or float:  
  Index precision;
  switch (sizeof (Numeric))
    {
    case sizeof (float):
      precision = FLT_DIG;
      break;
    case sizeof (double):
      precision = DBL_DIG;
      break;
    default:
      out0 << "Numeric must be double or float\n";
      exit (1);
  }

  os << setprecision (precision);
}


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
template<typename T> void
xml_read_from_file (const String& filename,
                    T&            type)
{
  ifstream ifs;

  out2 << "  Reading " << filename << '\n';
  
  // Open input stream:
  xml_open_input_file (ifs, filename);

  // No need to check for error, because xml_open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the matrix from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the 
  // filename.
  try
    {
      xml_read_header_from_stream (ifs);
      xml_read_from_stream (ifs, type);
      xml_read_footer_from_stream (ifs);
    }
  catch (runtime_error e)
    {
      ostringstream os;
      os << "Error reading file: " << filename << '\n'
         << e.what ();
      throw runtime_error (os.str ());
    }
}


//! Write data to XML file
/*!
  This is a generic functions that is used to write the XML header and
  footer info and calls the overloaded functions to write the data.
  
  \param filename XML filename
  \param type Generic input value
*/
template<typename T> void
xml_write_to_file (const String& filename,
                   const T&      type)
{
  ofstream ofs;

  out2 << "  Writing " << filename << '\n';
  xml_open_output_file(ofs, filename);

  try
    {
      xml_write_header_to_stream (ofs);
      xml_write_to_stream (ofs, type);
      xml_write_footer_to_stream (ofs);
    }
  catch (runtime_error e)
    {
      ostringstream os;
      os << "Error writing file: " << filename << '\n'
         << e.what ();
      throw runtime_error (os.str ());
    }

}


////////////////////////////////////////////////////////////////////////////
//   Overloaded functions for reading/writing data from/to XML stream
////////////////////////////////////////////////////////////////////////////

//=== ArrayOfArrayOfSpeciesTag ================================================

//! Reads ArrayOfArrayOfSpeciesTag from XML input stream
/*!
  Checks whether the next tag in input stream is
  <Array type="ArrayOfSpeciesTag"> and if so, write the values to 'aastag'
  parameter.

  \param is     Input stream
  \param aastag ArrayOfArrayOfSpeciesTag return value
*/
void
xml_read_from_stream (istream&                  is,
                      ArrayOfArrayOfSpeciesTag& aastag)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfSpeciesTag");

  tag.get_attribute_value ("nelem", nelem);
  aastag.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, aastag[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfSpeciesTag: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfSpeciesTag to XML output stream
/*!
  \param os     Output stream
  \param aastag ArrayOfArrayOfSpeciesTag
*/
void
xml_write_to_stream (ostream&                        os,
                     const ArrayOfArrayOfSpeciesTag& aastag)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "ArrayOfSpeciesTag");
  open_tag.add_attribute ("nelem", aastag.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < aastag.nelem (); n++)
    {
      xml_write_to_stream (os, aastag[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfArrayOfTensor3==================================================

//! Reads ArrayOfArrayOfTensor3 from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="ArrayOfTensor3">
  and if so, write the values to 'aatensor3' parameter.

  \param is Input stream
  \param aatensor3 ArrayOfArrayOfTensor3 return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfArrayOfTensor3& aatensor3)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfTensor3");

  tag.get_attribute_value ("nelem", nelem);
  aatensor3.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, aatensor3[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfTensor3: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfTensor3 to XML output stream
/*!
  \param os Output stream
  \param aatensor3 ArrayOfArrayOfTensor3
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfArrayOfTensor3& aatensor3)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "ArrayOfTensor3");
  open_tag.add_attribute ("nelem", aatensor3.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < aatensor3.nelem (); n++)
    {
      xml_write_to_stream (os, aatensor3[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfArrayOfTensor6==================================================

//! Reads ArrayOfArrayOfTensor6 from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="ArrayOfTensor6">
  and if so, write the values to 'aatensor6' parameter.

  \param is Input stream
  \param aatensor6 ArrayOfArrayOfTensor6 return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfArrayOfTensor6& aatensor6)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "ArrayOfTensor6");

  tag.get_attribute_value ("nelem", nelem);
  aatensor6.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, aatensor6[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfArrayOfTensor6: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfArrayOfTensor6 to XML output stream
/*!
  \param os        Output stream
  \param aatensor6 ArrayOfArrayOfTensor6
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfArrayOfTensor6& aatensor6)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "ArrayOfTensor6");
  open_tag.add_attribute ("nelem", aatensor6.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < aatensor6.nelem (); n++)
    {
      xml_write_to_stream (os, aatensor6[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfIndex ==========================================================

//! Reads ArrayOfIndex from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="Index">
  and if so, write the values to 'aindex' parameter.

  \param is Input stream
  \param aindex ArrayOfIndex return value
*/
void
xml_read_from_stream (istream&      is,
                      ArrayOfIndex& aindex)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Index");

  tag.get_attribute_value ("nelem", nelem);
  aindex.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, aindex[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfIndex: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfIndex to XML output stream
/*!
  \param os Output stream
  \param aindex ArrayOfIndex
*/
void
xml_write_to_stream (ostream&            os,
                     const ArrayOfIndex& aindex)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "Index");
  open_tag.add_attribute ("nelem", aindex.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < aindex.nelem (); n++)
    {
      xml_write_to_stream (os, aindex[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== ArrayOfMatrix ==========================================================

//! Reads ArrayOfMatrix from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="Matrix">
  and if so, write the values to 'amatrix' parameter.

  \param is Input stream
  \param amatrix ArrayOfMatrix return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfMatrix& amatrix)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Matrix");

  tag.get_attribute_value ("nelem", nelem);
  amatrix.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, amatrix[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfMatrix: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfMatrix to XML output stream
/*!
  \param os Output stream
  \param amatrix ArrayOfMatrix
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfMatrix& amatrix)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "Matrix");
  open_tag.add_attribute ("nelem", amatrix.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < amatrix.nelem (); n++)
    {
      xml_write_to_stream (os, amatrix[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfSpeciesTag ================================================

//! Reads ArrayOfSpeciesTag from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="SpeciesTag">
  and if so, write the values to 'astag' parameter.

  \param is    Input stream
  \param astag ArrayOfSpeciesTag return value
*/
void
xml_read_from_stream (istream&           is,
                      ArrayOfSpeciesTag& astag)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "SpeciesTag");

  tag.get_attribute_value ("nelem", nelem);
  astag.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, astag[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfSpeciesTag: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfSpeciesTag to XML output stream
/*!
  \param os    Output stream
  \param astag ArrayOfSpeciesTag
*/
void
xml_write_to_stream (ostream&                 os,
                     const ArrayOfSpeciesTag& astag)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "SpeciesTag");
  open_tag.add_attribute ("nelem", astag.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < astag.nelem (); n++)
    {
      xml_write_to_stream (os, astag[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfTensor3=========================================================

//! Reads ArrayOfTensor3 from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="Tensor3">
  and if so, write the values to 'atensor3' parameter.

  \param is Input stream
  \param atensor3 ArrayOfTensor3 return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfTensor3& atensor3)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Tensor3");

  tag.get_attribute_value ("nelem", nelem);
  atensor3.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, atensor3[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfTensor3: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfTensor3 to XML output stream
/*!
  \param os Output stream
  \param atensor3 ArrayOfTensor3
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfTensor3& atensor3)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "Tensor3");
  open_tag.add_attribute ("nelem", atensor3.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < atensor3.nelem (); n++)
    {
      xml_write_to_stream (os, atensor3[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfTensor6=========================================================

//! Reads ArrayOfTensor6 from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="Tensor6">
  and if so, write the values to 'atensor6' parameter.

  \param is Input stream
  \param atensor6 ArrayOfTensor6 return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfTensor6& atensor6)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Tensor6");

  tag.get_attribute_value ("nelem", nelem);
  atensor6.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, atensor6[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfTensor6: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfTensor6 to XML output stream
/*!
  \param os Output stream
  \param atensor6 ArrayOfTensor6
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfTensor6& atensor6)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "Tensor6");
  open_tag.add_attribute ("nelem", atensor6.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < atensor6.nelem (); n++)
    {
      xml_write_to_stream (os, atensor6[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== ArrayOfString ==========================================================

//! Reads ArrayOfString from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="String">
  and if so, write the values to 'astring' parameter.

  \param is Input stream
  \param astring ArrayOfString return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfString& astring)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "String");

  tag.get_attribute_value ("nelem", nelem);
  astring.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, astring[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfString: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }

  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfString to XML output stream
/*!
  \param os Output stream
  \param astring ArrayOfString
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfString& astring)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "String");
  open_tag.add_attribute ("nelem", astring.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < astring.nelem (); n++)
    {
      xml_write_to_stream (os, astring[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== ArrayOfVector ==========================================================

//! Reads ArrayOfVector from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="Vector">
  and if so, write the values to 'amatrix' parameter.

  \param is Input stream
  \param avector ArrayOfVector return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfVector& avector)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  tag.check_attribute ("type", "Vector");

  tag.get_attribute_value ("nelem", nelem);
  avector.resize (nelem);

  Index n;
  try
    {
      for (n = 0; n < nelem; n++)
        {
          xml_read_from_stream (is, avector[n]);
        }
    } catch (runtime_error e) {
      ostringstream os;
      os << "Error reading ArrayOfVector: "
         << "\n Element: " << n
         << "\n" << e.what();
      throw runtime_error(os.str());
    }


  tag.read_from_stream (is);
  tag.check_name ("/Array");
}


//! Writes ArrayOfVector to XML output stream
/*!
  \param os Output stream
  \param amatrix ArrayOfVector
*/
void
xml_write_to_stream (ostream&             os,
                     const ArrayOfVector& avector)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Array");

  open_tag.add_attribute ("type", "Vector");
  open_tag.add_attribute ("nelem", avector.nelem ());

  open_tag.write_to_stream (os);
  os << '\n';

  for (Index n = 0; n < avector.nelem (); n++)
    {
      xml_write_to_stream (os, avector[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << '\n';
}


//=== GasAbsLookup ===========================================================

//! Reads GasAbsLookup from XML input stream
/*!
  Checks whether the next tag in input stream is <GasAbsLookup>
  and if so, verifies the order of the components and writes them
  to the 'gal' parameter.

  \param is  Input stream
  \param gal GasAbsLookup return value
*/
void
xml_read_from_stream (istream&      is,
                      GasAbsLookup& gal)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is);
  tag.check_name ("GasAbsLookup");

  xml_read_from_stream (is, gal.species);
  xml_read_from_stream (is, gal.nonlinear_species);
  xml_read_from_stream (is, gal.f_grid);
  xml_read_from_stream (is, gal.p_grid);
  xml_read_from_stream (is, gal.vmrs_ref);
  xml_read_from_stream (is, gal.t_ref);
  xml_read_from_stream (is, gal.t_pert);
  xml_read_from_stream (is, gal.nls_pert);
  xml_read_from_stream (is, gal.abs);

  tag.read_from_stream (is);
  tag.check_name ("/GasAbsLookup");
}


//! Writes GasAbsLookup to XML output stream
/*!
  \param os Output stream
  \param gas GasAbsLookup
*/
void
xml_write_to_stream (ostream&            os,
                     const GasAbsLookup& gal)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("GasAbsLookup");
  open_tag.write_to_stream (os);

  xml_write_to_stream (os, gal.species);
  xml_write_to_stream (os, gal.nonlinear_species);
  xml_write_to_stream (os, gal.f_grid);
  xml_write_to_stream (os, gal.p_grid);
  xml_write_to_stream (os, gal.vmrs_ref);
  xml_write_to_stream (os, gal.t_ref);
  xml_write_to_stream (os, gal.t_pert);
  xml_write_to_stream (os, gal.nls_pert);
  xml_write_to_stream (os, gal.abs);

  close_tag.set_name ("/GasAbsLookup");
  close_tag.write_to_stream (os);
  os << '\n';
}


//=== Index ==================================================================

//! Reads Index from XML input stream
/*!
  Checks whether the next tag in input stream is <Index> and if so,
  write the value to 'index' parameter.

  \param is Input stream
  \param index Index return value
*/
void
xml_read_from_stream (istream& is,
                      Index&   index)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is);
  tag.check_name ("Index");
  
  is >> index;
  if (is.fail ())
    {
      xml_parse_error ("Error while reading data");
    }

  out2 << "Data:" << index << '\n';

  tag.read_from_stream (is);
  tag.check_name ("/Index");
}


//! Writes Index to XML output stream
/*!
  \param os Output stream
  \param index Index value
*/
void
xml_write_to_stream (ostream&     os,
                     const Index& index)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Index");

  open_tag.write_to_stream (os);

  os << index;

  close_tag.set_name ("/Index");
  close_tag.write_to_stream (os);
  os << '\n';
}

//=== Matrix ==========================================================

//! Reads Matrix from XML input stream
/*!
  Checks whether the next tag in input stream is <Matrix> and if so,
  write the values to 'matrix' parameter.

  \param is Input stream
  \param matrix Matrix return value
*/
void
xml_read_from_stream (istream& is,
                      Matrix& matrix)
{
  ArtsXMLTag tag;
  Index nrows, ncols;

  tag.read_from_stream (is);
  tag.check_name ("Matrix");

  tag.get_attribute_value ("nrows", nrows);
  tag.get_attribute_value ("ncols", ncols);
  matrix.resize (nrows, ncols);

  for (Index r = 0; r < nrows; r++)
    {
      for (Index c = 0; c < ncols; c++)
        {
          is >> matrix (r, c);
          if (is.fail ())
            {
              ostringstream os;
              os << "Error reading Matrix:"
                 << "\n  Row   : " << r
                 << "\n  Column: " << c;
              xml_parse_error (os.str());
            }
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Matrix");
}


//! Writes Matrix to XML output stream
/*!
  \param os Output stream
  \param matrix Matrix
*/
void
xml_write_to_stream (ostream&     os,
                     const Matrix& matrix)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Matrix");
  open_tag.add_attribute ("nrows", matrix.nrows ());
  open_tag.add_attribute ("ncols", matrix.ncols ());

  open_tag.write_to_stream (os);
  os << '\n';

  xml_set_stream_precision (os);

  // Write the elements:
  for (Index r = 0; r < matrix.nrows (); ++r)
    {
      os << matrix (r,0);
      
      for (Index c = 1; c < matrix.ncols (); ++c)
        {
          os << " " << matrix (r,c);
        }
      
      os << '\n';
    }

  close_tag.set_name ("/Matrix");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== Numeric =========================================================

//! Reads Numeric from XML input stream
/*!
  Checks whether the next tag in input stream is <Numeric> and if so,
  write the value to 'numeric' parameter.

  \param is Input stream
  \param numeric Numeric return value
*/
void
xml_read_from_stream (istream& is,
                      Numeric& numeric)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is);
  tag.check_name ("Numeric");

  is >> numeric;
  if (is.fail ())
    {
      xml_parse_error ("Error while reading data");
    }

  out2 << "Data:" << numeric << '\n';

  tag.read_from_stream (is);
  tag.check_name ("/Numeric");
}


//! Writes Numeric to XML output stream
/*!
  \param os Output stream
  \param numeric Numeric value
*/
void
xml_write_to_stream (ostream&       os,
                     const Numeric& numeric)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Numeric");

  open_tag.write_to_stream (os);

  xml_set_stream_precision (os);

  os << numeric;

  close_tag.set_name ("/Numeric");
  close_tag.write_to_stream (os);
  os << '\n';
}


//=== SpeciesTag ================================================

//! Reads SpeciesTag from XML input stream
/*!
  Checks whether the next tag in input stream is <SpeciesTag>
  and if so, write the values to 'stag' parameter.

  \param is   Input stream
  \param stag SpeciesTag return value
*/
void
xml_read_from_stream (istream&    is,
                      SpeciesTag& stag)
{
  ArtsXMLTag tag;
  stringbuf  strbuf;
  char dummy;

  tag.read_from_stream (is);
  tag.check_name ("SpeciesTag");

  is.get (strbuf, '"');
  if (is.fail ())
    {
      xml_parse_error ("SpeciesTag must begin with \"");
    }

  is >> dummy;
  is.get (strbuf, '"');
  if (is.fail ())
    {
      xml_parse_error ("SpeciesTag must end with \"");
    }

  stag = SpeciesTag (strbuf.str ());

  tag.read_from_stream (is);
  tag.check_name ("/SpeciesTag");
}


//! Writes SpeciesTag to XML output stream
/*!
  \param os   Output stream
  \param stag SpeciesTag
*/
void
xml_write_to_stream (ostream&          os,
                     const SpeciesTag& stag)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("SpeciesTag");
  open_tag.write_to_stream (os);

  os << '\"' << stag.Name () << '\"';

  close_tag.set_name ("/SpeciesTag");
  close_tag.write_to_stream (os);
  os << '\n';
}


//=== String ===========================================================

//! Reads String from XML input stream
/*!
  Checks whether the next tag in input stream is <String> and if so,
  write the value to 'str' parameter.

  \param is Input stream
  \param str String return value
*/
void
xml_read_from_stream (istream& is,
                      String&  str)
{
  ArtsXMLTag tag;
  stringbuf strbuf;

  tag.read_from_stream (is);
  tag.check_name ("String");
  
  is.get (strbuf, '<');
  if (is.fail ())
    {
      xml_parse_error ("Error while reading data");
    }
  
  str = strbuf.str ();

  tag.read_from_stream (is);
  tag.check_name ("/String");
}


//! Writes String to XML output stream
/*!
  \param os Output stream
  \param str String value
*/
void
xml_write_to_stream (ostream&     os,
                     const String& str)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("String");

  open_tag.write_to_stream (os);

  os << str;

  close_tag.set_name ("/String");
  close_tag.write_to_stream (os);
  os << '\n';
}


//=== Tensor3 ================================================================

//! Reads Tensor3 from XML input stream
/*!
  Checks whether the next tag in input stream is <Tensor3> and if so,
  write the values to 'tensor' parameter.

  \param is Input stream
  \param tensor Tensor return value
*/
void
xml_read_from_stream (istream& is,
                      Tensor3& tensor)
{
  ArtsXMLTag tag;
  Index npages, nrows, ncols;

  tag.read_from_stream (is);
  tag.check_name ("Tensor3");

  tag.get_attribute_value ("npages", npages);
  tag.get_attribute_value ("nrows", nrows);
  tag.get_attribute_value ("ncols", ncols);
  tensor.resize (npages, nrows, ncols);

  for (Index p = 0; p < npages; p++)
    {
      for (Index r = 0; r < nrows; r++)
        {
          for (Index c = 0; c < ncols; c++)
            {
              is >> tensor (p, r, c);
              if (is.fail ())
                {
                  ostringstream os;
                  os << "Error reading Tensor3:"
                     << "\n  Page  : " << p
                     << "\n  Row   : " << r
                     << "\n  Column: " << c;
                  xml_parse_error (os.str());
                }
            }
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Tensor3");
}


//! Writes Tensor3 to XML output stream
/*!
  \param os Output stream
  \param tensor Tensor
*/
void
xml_write_to_stream (ostream&     os,
                     const Tensor3& tensor)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Tensor3");
  open_tag.add_attribute ("npages", tensor.npages ());
  open_tag.add_attribute ("nrows", tensor.nrows ());
  open_tag.add_attribute ("ncols", tensor.ncols ());

  open_tag.write_to_stream (os);
  os << '\n';

  xml_set_stream_precision (os);

  // Write the elements:
  for (Index p = 0; p < tensor.npages (); ++p)
    {
      for (Index r = 0; r < tensor.nrows (); ++r)
        {
          os << tensor (p, r, 0);
          for (Index c = 1; c < tensor.ncols (); ++c)
            {
              os << " " << tensor (p, r, c);
            }
          os << '\n';
        }
    }

  close_tag.set_name ("/Tensor3");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== Tensor4 =========================================================

//! Reads Tensor4 from XML input stream
/*!
  Checks whether the next tag in input stream is <Tensor4> and if so,
  write the values to 'tensor' parameter.

  \param is Input stream
  \param tensor Tensor return value
*/
void
xml_read_from_stream (istream& is,
                      Tensor4& tensor)
{
  ArtsXMLTag tag;
  Index nbooks, npages, nrows, ncols;

  tag.read_from_stream (is);
  tag.check_name ("Tensor4");

  tag.get_attribute_value ("nbooks", nbooks);
  tag.get_attribute_value ("npages", npages);
  tag.get_attribute_value ("nrows", nrows);
  tag.get_attribute_value ("ncols", ncols);
  tensor.resize (nbooks,npages, nrows, ncols);

  for (Index b = 0; b < nbooks; b++)
    {
      for (Index p = 0; p < npages; p++)
        {
          for (Index r = 0; r < nrows; r++)
            {
              for (Index c = 0; c < ncols; c++)
                {
                  is >> tensor (b, p, r, c);
                  if (is.fail ())
                    {
                      ostringstream os;
                      os << "Error reading Tensor4:"
                         << "\n  Book  : " << b
                         << "\n  Page  : " << p
                         << "\n  Row   : " << r
                         << "\n  Column: " << c;
                      xml_parse_error (os.str());
                    }
                }
            }
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Tensor4");
}


//! Writes Tensor4 to XML output stream
/*!
  \param os Output stream
  \param tensor Tensor
*/
void
xml_write_to_stream (ostream&     os,
                     const Tensor4& tensor)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Tensor4");
  open_tag.add_attribute ("nbooks", tensor.nbooks ());
  open_tag.add_attribute ("npages", tensor.npages ());
  open_tag.add_attribute ("nrows", tensor.nrows ());
  open_tag.add_attribute ("ncols", tensor.ncols ());

  open_tag.write_to_stream (os);
  os << '\n';

  xml_set_stream_precision (os);

  // Write the elements:
  for (Index b = 0; b < tensor.nbooks (); ++b)
    {
      for (Index p = 0; p < tensor.npages (); ++p)
        {
          for (Index r = 0; r < tensor.nrows (); ++r)
            {
              os << tensor (b, p, r, 0);
              for (Index c = 1; c < tensor.ncols (); ++c)
                {
                  os << " " << tensor (b, p, r, c);
                }
              os << '\n';
            }
        }
    }

  close_tag.set_name ("/Tensor4");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== Tensor5 =========================================================

//! Reads Tensor5 from XML input stream
/*!
  Checks whether the next tag in input stream is <Tensor5> and if so,
  write the values to 'tensor' parameter.

  \param is Input stream
  \param tensor Tensor return value
*/
void
xml_read_from_stream (istream& is,
                      Tensor5& tensor)
{
  ArtsXMLTag tag;
  Index nshelves, nbooks, npages, nrows, ncols;

  tag.read_from_stream (is);
  tag.check_name ("Tensor5");

  tag.get_attribute_value ("nshelves", nshelves);
  tag.get_attribute_value ("nbooks", nbooks);
  tag.get_attribute_value ("npages", npages);
  tag.get_attribute_value ("nrows", nrows);
  tag.get_attribute_value ("ncols", ncols);
  tensor.resize (nshelves, nbooks,npages, nrows, ncols);

  for (Index s = 0; s < nshelves; s++)
    {
      for (Index b = 0; b < nbooks; b++)
        {
          for (Index p = 0; p < npages; p++)
            {
              for (Index r = 0; r < nrows; r++)
                {
                  for (Index c = 0; c < ncols; c++)
                    {
                      is >> tensor (s, b, p, r, c);
                      if (is.fail ())
                        {
                          ostringstream os;
                          os << "Error reading Tensor5:"
                             << "\n  Shelf : " << s
                             << "\n  Book  : " << b
                             << "\n  Page  : " << p
                             << "\n  Row   : " << r
                             << "\n  Column: " << c;
                          xml_parse_error (os.str());
                        }
                    }
                }
            }
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Tensor5");
}


//! Writes Tensor5 to XML output stream
/*!
  \param os Output stream
  \param tensor Tensor
*/
void
xml_write_to_stream (ostream&     os,
                     const Tensor5& tensor)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Tensor5");
  open_tag.add_attribute ("nshelves", tensor.nshelves ());
  open_tag.add_attribute ("nbooks", tensor.nbooks ());
  open_tag.add_attribute ("npages", tensor.npages ());
  open_tag.add_attribute ("nrows", tensor.nrows ());
  open_tag.add_attribute ("ncols", tensor.ncols ());

  open_tag.write_to_stream (os);
  os << '\n';

  xml_set_stream_precision (os);

  // Write the elements:
  for (Index s = 0; s < tensor.nshelves (); ++s)
    {
      for (Index b = 0; b < tensor.nbooks (); ++b)
        {
          for (Index p = 0; p < tensor.npages (); ++p)
            {
              for (Index r = 0; r < tensor.nrows (); ++r)
                {
                  os << tensor (s, b, p, r, 0);
                  for (Index c = 1; c < tensor.ncols (); ++c)
                    {
                      os << " " << tensor (s, b, p, r, c);
                    }
                  os << '\n';
                }
            }
        }
    }

  close_tag.set_name ("/Tensor5");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== Tensor6 =========================================================

//! Reads Tensor6 from XML input stream
/*!
  Checks whether the next tag in input stream is <Tensor6> and if so,
  write the values to 'tensor' parameter.

  \param is Input stream
  \param tensor Tensor return value
*/
void
xml_read_from_stream (istream& is,
                      Tensor6& tensor)
{
  ArtsXMLTag tag;
  Index nvitrines, nshelves, nbooks, npages, nrows, ncols;

  tag.read_from_stream (is);
  tag.check_name ("Tensor6");

  tag.get_attribute_value ("nvitrines", nvitrines);
  tag.get_attribute_value ("nshelves", nshelves);
  tag.get_attribute_value ("nbooks", nbooks);
  tag.get_attribute_value ("npages", npages);
  tag.get_attribute_value ("nrows", nrows);
  tag.get_attribute_value ("ncols", ncols);
  tensor.resize (nvitrines, nshelves, nbooks,npages, nrows, ncols);

  for (Index v = 0; v < nvitrines; v++)
    {
      for (Index s = 0; s < nshelves; s++)
        {
          for (Index b = 0; b < nbooks; b++)
            {
              for (Index p = 0; p < npages; p++)
                {
                  for (Index r = 0; r < nrows; r++)
                    {
                      for (Index c = 0; c < ncols; c++)
                        {
                          is >> tensor (v, s, b, p, r, c);
                          if (is.fail ())
                            {
                              ostringstream os;
                              os << "Error reading Tensor6:"
                                 << "\n  Vitrine: " << v
                                 << "\n  Shelf  : " << s
                                 << "\n  Book   : " << b
                                 << "\n  Page   : " << p
                                 << "\n  Row    : " << r
                                 << "\n  Column : " << c;
                              xml_parse_error (os.str());
                            }
                        }
                    }
                }
            }
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Tensor6");
}


//! Writes Tensor6 to XML output stream
/*!
  \param os Output stream
  \param tensor Tensor
*/
void
xml_write_to_stream (ostream&     os,
                     const Tensor6& tensor)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Tensor6");
  open_tag.add_attribute ("nvitrines", tensor.nvitrines ());
  open_tag.add_attribute ("nshelves", tensor.nshelves ());
  open_tag.add_attribute ("nbooks", tensor.nbooks ());
  open_tag.add_attribute ("npages", tensor.npages ());
  open_tag.add_attribute ("nrows", tensor.nrows ());
  open_tag.add_attribute ("ncols", tensor.ncols ());

  open_tag.write_to_stream (os);
  os << '\n';

  xml_set_stream_precision (os);

  // Write the elements:
  for (Index v = 0; v < tensor.nvitrines (); ++v)
    {
      for (Index s = 0; s < tensor.nshelves (); ++s)
        {
          for (Index b = 0; b < tensor.nbooks (); ++b)
            {
              for (Index p = 0; p < tensor.npages (); ++p)
                {
                  for (Index r = 0; r < tensor.nrows (); ++r)
                    {
                      os << tensor (v, s, b, p, r, 0);
                      for (Index c = 1; c < tensor.ncols (); ++c)
                        {
                          os << " " << tensor (v, s, b, p, r, c);
                        }
                      os << '\n';
                    }
                }
            }
        }
    }

  close_tag.set_name ("/Tensor6");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== Tensor7 =========================================================

//! Reads Tensor7 from XML input stream
/*!
  Checks whether the next tag in input stream is <Tensor7> and if so,
  write the values to 'tensor' parameter.

  \param is Input stream
  \param tensor Tensor return value
*/
void
xml_read_from_stream (istream& is,
                      Tensor7& tensor)
{
  ArtsXMLTag tag;
  Index nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols;

  tag.read_from_stream (is);
  tag.check_name ("Tensor7");

  tag.get_attribute_value ("nlibraries", nlibraries);
  tag.get_attribute_value ("nvitrines", nvitrines);
  tag.get_attribute_value ("nshelves", nshelves);
  tag.get_attribute_value ("nbooks", nbooks);
  tag.get_attribute_value ("npages", npages);
  tag.get_attribute_value ("nrows", nrows);
  tag.get_attribute_value ("ncols", ncols);
  tensor.resize (nlibraries, nvitrines, nshelves, nbooks,npages, nrows, ncols);

  for (Index l = 0; l < nlibraries; l++)
    {
      for (Index v = 0; v < nvitrines; v++)
        {
          for (Index s = 0; s < nshelves; s++)
            {
              for (Index b = 0; b < nbooks; b++)
                {
                  for (Index p = 0; p < npages; p++)
                    {
                      for (Index r = 0; r < nrows; r++)
                        {
                          for (Index c = 0; c < ncols; c++)
                            {
                              is >> tensor (l, v, s, b, p, r, c);
                              if (is.fail ())
                                {
                                  xml_parse_error ("Error while reading data");
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Tensor7");
}


//! Writes Tensor7 to XML output stream
/*!
  \param os Output stream
  \param tensor Tensor
*/
void
xml_write_to_stream (ostream&     os,
                     const Tensor7& tensor)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Tensor7");
  open_tag.add_attribute ("nlibraries", tensor.nlibraries ());
  open_tag.add_attribute ("nvitrines", tensor.nvitrines ());
  open_tag.add_attribute ("nshelves", tensor.nshelves ());
  open_tag.add_attribute ("nbooks", tensor.nbooks ());
  open_tag.add_attribute ("npages", tensor.npages ());
  open_tag.add_attribute ("nrows", tensor.nrows ());
  open_tag.add_attribute ("ncols", tensor.ncols ());

  open_tag.write_to_stream (os);
  os << '\n';

  xml_set_stream_precision (os);

  // Write the elements:
  for (Index l = 0; l < tensor.nlibraries (); ++l)
    {
      for (Index v = 0; v < tensor.nvitrines (); ++v)
        {
          for (Index s = 0; s < tensor.nshelves (); ++s)
            {
              for (Index b = 0; b < tensor.nbooks (); ++b)
                {
                  for (Index p = 0; p < tensor.npages (); ++p)
                    {
                      for (Index r = 0; r < tensor.nrows (); ++r)
                        {
                          os << tensor (l, v, s, b, p, r, 0);
                          for (Index c = 1; c < tensor.ncols (); ++c)
                            {
                              os << " " << tensor (l, v, s, b, p, r, c);
                            }
                          os << '\n';
                        }
                    }
                }
            }
        }
    }

  close_tag.set_name ("/Tensor7");
  close_tag.write_to_stream (os);

  os << '\n';
}

//=== Vector ==========================================================

//! Reads Vector from XML input stream
/*!
  Checks whether the next tag in input stream is <Vector> and if so,
  write the values to 'vector' parameter.

  \param is Input stream
  \param vector Vector return value
*/
void
xml_read_from_stream (istream& is,
                      Vector& vector)
{
  ArtsXMLTag tag;
  Index nelem;

  tag.read_from_stream (is);
  tag.check_name ("Vector");
  
  tag.get_attribute_value ("nelem", nelem);
  vector.resize (nelem);
  
  for (Index n = 0; n < nelem; n++)
    {
      is >> vector[n];
      if (is.fail ())
        {
          ostringstream os;
          os << "Error reading Vector:"
             << "\n  Element: " << n;
          xml_parse_error (os.str());
        }
    }

  tag.read_from_stream (is);
  tag.check_name ("/Vector");
}


//! Writes Vector to XML output stream
/*!
  \param os Output stream
  \param vector Vector
*/
void
xml_write_to_stream (ostream&     os,
                     const Vector& vector)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;
  Index n = vector.nelem ();
  ostringstream v;

  // Convert nelem to string
  v << n;

  open_tag.set_name ("Vector");
  open_tag.add_attribute ("nelem", v.str ());

  open_tag.write_to_stream (os);

  os << '\n';
  for (Index i=0; i<n; ++i)
    os << vector[i] << '\n';

  close_tag.set_name ("/Vector");
  close_tag.write_to_stream (os);

  os << '\n';
}


////////////////////////////////////////////////////////////////////////////
//   Dummy funtion for groups for which 
//   IO function have not yet been implemented
////////////////////////////////////////////////////////////////////////////

// FIXME: These should be implemented, sooner or later...

void
xml_read_from_stream (istream& is,
                      Ppath& ppath)
{
  throw runtime_error("Method not implemented!");
}

void
xml_write_to_stream (ostream&     os,
                     const Ppath& ppath)
{
  throw runtime_error("Method not implemented!");
}

void
xml_read_from_stream (istream& is,
                      Agenda& agenda)
{
  throw runtime_error("Method not implemented!");
}

void
xml_write_to_stream (ostream&     os,
                     const Agenda& agenda)
{
  throw runtime_error("Method not implemented!");
}

void
xml_read_from_stream (istream& is,
                      GridPos& gp)
{
  throw runtime_error("Method not implemented!");
}

void
xml_write_to_stream (ostream&     os,
                     const GridPos& gp)
{
  throw runtime_error("Method not implemented!");
}


