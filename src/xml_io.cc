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
  \date   Fri May 10 09:43:57 2002
  
  \brief This file contains basic functions to handle XML data files.
  
*/


#include "xml_io.h"
#include <stdexcept>
#include "array.h"
#include "messages.h"



void
xml_parse_error (const String& str_error);

void
xml_read_header (istream& is);

void
xml_read_footer (istream& is);



class XMLAttribute
{
public:
  String name;
  String value;
};


class ArtsXMLTag
{
public:
  String&
  getName () { return (name); }

  void
  setName (const String& name);

  void
  addAttribute (const String& name, const String& value);

  void
  getAttributeValue (const String& name, String& value);

  void
  checkName (const String& expected_name);

  void
  readTag (istream& is);

  void
  writeStartTag (ostream& os);

  void
  writeEndTag (ostream& os);

private:
  String name;
  Array<XMLAttribute> attribs;
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
ArtsXMLTag::checkName (const String& expected_name)
{
  if (name != expected_name)
    {
      xml_parse_error ("Tag <" + expected_name + "> expected but <"
                       + getName () + "> found.");
    }
  
}


//! Reads next XML tag
/*!
  Reads the name and attributes of the next XML tag from stream.
  
  \param is Input stream
*/
void
ArtsXMLTag::readTag (istream& is)
{
  String       token;
  stringbuf    tag;
  istringstream sstr("");
  XMLAttribute attr;

  attribs.clear ();

  while (is.good () && isspace (is.peek ()))
    {
      is.get ();
    }

  is.get (tag, '>');

  // Hit EOF while looking for '>'
  if (is.bad () || is.eof ())
    {
      xml_parse_error ("Unexpected end of file while looking for '>'");
    }

  if (tag.str ()[0] != '<')
    {
      xml_parse_error ("Expected tag but found: "
                       + tag.str ());
    }

  if (is.get () != '>')
    {
      xml_parse_error ("Closing > not found in tag: " + tag.str ());
    }


  sstr.str (tag.str () + '>');
  cout << "Read: " << sstr.str () << endl;

  sstr >> name;

  // Remove opening < of tag
  name.erase (0, 1);

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
  cout << "Name: " << name << endl;

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
    if (pos == String::npos)
      {
        xml_parse_error ("Missing \" in tag: " + sstr.str ());
      }

    attr.value = token.substr (1, pos - 1);

    attribs.push_back (attr);

    cout << "Attr: " << attr.name << endl;
    cout << "Value: " << attr.value << endl;

    if (token[token.length () - 1] == '>')
      {
        token = ">";
      }
    else
      {
        sstr >> token;
      }
  }

  cout << endl;
}



////////////////////////////////////////////////////////////////////////////
//   Default file names
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
{}


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
xml_read_header (istream& is)
{
  String str;
  stringbuf strbuf;
  ArtsXMLTag tag;

  is.get (strbuf, '>');
  if (is.get () != '>' || strbuf.str () != "<?xml version=\"1.0\"?")
    {
      xml_parse_error ("First line does not contain: "
                       "<?xml version=\"1.0\"?>");
    }

  tag.readTag (is);
  tag.checkName ("arts");
}


//! Reads closing root tag
/*! 
  Checks whether XML file ends correctly with </arts>.
  
  \param is Input stream
*/
void
xml_read_footer (istream& is)
{
  ArtsXMLTag tag;

  tag.readTag (is);
  tag.checkName ("/arts");
}



////////////////////////////////////////////////////////////////////////////
//   Index IO routines for XML files
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
                          T&      type)
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
      xml_read_header (ifs);
      xml_read_from_stream (ifs, type);
      xml_read_footer (ifs);
    }
  catch (runtime_error e)
    {
      ostringstream os;
      os << "Error reading file: " << filename << '\n'
         << e.what ();
      throw runtime_error (os.str ());
    }
}



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

  tag.readTag (is);
  tag.checkName ("Index");
  
  is >> index;
  if (is.fail ())
    {
      xml_parse_error ("Error on reading data");
    }

  cout << "Data:" << index << endl;

  tag.readTag (is);
  tag.checkName ("/Index");
}


void
xml_write_index_to_file (const String& filename,
                         const Index&  index)
{}


void
xml_write_index_to_stream (ostream&     os,
                           const Index& index)
{}


////////////////////////////////////////////////////////////////////////////
//   Numeric IO routines for XML files
////////////////////////////////////////////////////////////////////////////

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

  tag.readTag (is);
  tag.checkName ("Numeric");
  
  is >> numeric;
  if (is.fail ())
    {
      xml_parse_error ("Error on reading data");
    }

  cout << "Data:" << numeric << endl;

  tag.readTag (is);
  tag.checkName ("/Numeric");
}



////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of template functions we need
////////////////////////////////////////////////////////////////////////////

template void xml_read_from_file<Index> (const String&, Index&);
template void xml_read_from_file<Numeric> (const String&, Numeric&);
