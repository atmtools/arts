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
#include <cfloat>
#include "array.h"
#include "messages.h"



void
xml_parse_error (const String& str_error);

void
xml_read_header_from_stream (istream& is);

void
xml_read_footer_from_stream (istream& is);

void
xml_write_header_to_stream (ostream& os);

void
xml_write_footer_to_stream (ostream& os);



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
  get_name () { return (name); }

  void
  check_name (const String& expected_name);

  void
  set_name (const String& new_name) { name = new_name; }

  void
  add_attribute (const String& aname, const String& value);

  void
  get_attribute_value (const String& aname, String& value);

  void
  read_from_stream (istream& is);

  void
  write_to_stream (ostream& os);

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
ArtsXMLTag::check_name (const String& expected_name)
{
  if (name != expected_name)
    {
      xml_parse_error ("Tag <" + expected_name + "> expected but <"
                       + name + "> found.");
    }
  
}


//! Adds an attribute to tag
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


//! Reads next XML tag
/*!
  Reads the name and attributes of the next XML tag from stream.
  
  \param is Input stream
*/
void
ArtsXMLTag::read_from_stream (istream& is)
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
      os << ' ' << (*it).name
         << "=\"" << (*it).value << '\"';
      it++;
    }

  os << ">";
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
  String str;
  stringbuf strbuf;
  ArtsXMLTag tag;

  is.get (strbuf, '>');
  if (is.get () != '>' || strbuf.str () != "<?xml version=\"1.0\"?")
    {
      xml_parse_error ("First line does not contain: "
                       "<?xml version=\"1.0\"?>");
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
     << endl;

  tag.set_name ("arts");
  tag.add_attribute ("format", "ascii");
  tag.add_attribute ("version", "1");
  tag.write_to_stream (os);

  os << endl;

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

//=== Index ===========================================================

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

  cout << "Data:" << index << endl;

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
  close_tag.set_name ("/Index");

  open_tag.write_to_stream (os);
  os << index;
  close_tag.write_to_stream (os);

  os << endl;
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
      xml_parse_error ("Error on reading data");
    }

  cout << "Data:" << numeric << endl;

  tag.read_from_stream (is);
  tag.check_name ("/Numeric");
}


//! Writes Numeric to XML output stream
/*!
  \param os Output stream
  \param numeric Numeric value
*/
void
xml_write_to_stream (ostream&     os,
                     const Numeric& numeric)
{
  ArtsXMLTag open_tag;
  ArtsXMLTag close_tag;

  open_tag.set_name ("Numeric");
  close_tag.set_name ("/Numeric");

  open_tag.write_to_stream (os);
  os << numeric;
  close_tag.write_to_stream (os);

  os << endl;
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

  tag.read_from_stream (is);
  tag.check_name ("Vector");
  
  //  is >> numeric;
  if (is.fail ())
    {
      xml_parse_error ("Error on reading data");
    }

  //  cout << "Data:" << numeric << endl;

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

  os << endl;
  for (Index i=0; i<n; ++i)
    os << vector[i] << endl;

  close_tag.set_name ("/Vector");
  close_tag.write_to_stream (os);

  os << endl;
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

  tag.read_from_stream (is);
  tag.check_name ("Matrix");
  
  //  is >> numeric;
  if (is.fail ())
    {
      xml_parse_error ("Error on reading data");
    }

  //  cout << "Data:" << numeric << endl;

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

  {
    ostringstream v;
    v << matrix.nrows ();
    open_tag.add_attribute ("nrows", v.str ());
  }

  {
    ostringstream v;
    v << matrix.ncols ();
    open_tag.add_attribute ("ncols", v.str ());
  }

  open_tag.write_to_stream (os);

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
  // Write the elements:
  for (Index r = 0; r < matrix.nrows (); ++r)
    {
      os << matrix (r,0);
      
      for (Index c = 1; c < matrix.ncols (); ++c)
        {
          os << " " << matrix (r,c);
        }
      
      os << endl;
    }

  os << endl;
  close_tag.set_name ("/Matrix");
  close_tag.write_to_stream (os);

  os << endl;
}

//=== ArrayOfMatrix ==========================================================

//! Reads ArrayOfMatrix from XML input stream
/*!
  Checks whether the next tag in input stream is <Array type="matrix">
  and if so, write the values to 'amatrix' parameter.

  \param is Input stream
  \param amatrix ArrayOfMatrix return value
*/
void
xml_read_from_stream (istream&       is,
                      ArrayOfMatrix& amatrix)
{
  ArtsXMLTag tag;

  tag.read_from_stream (is);
  tag.check_name ("Array");
  
  //  is >> numeric;
  if (is.fail ())
    {
      xml_parse_error ("Error on reading data");
    }

  //  cout << "Data:" << numeric << endl;

  tag.read_from_stream (is);
  tag.check_name ("/Matrix");
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
  {
    ostringstream v;
    v << amatrix.nelem ();
    open_tag.add_attribute ("nelem", v.str ());
  }

  open_tag.write_to_stream (os);

  for (Index n = 0; n < amatrix.nelem (); n++)
    {
      xml_write_to_stream (os, amatrix[n]);
    }

  close_tag.set_name ("/Array");
  close_tag.write_to_stream (os);

  os << endl;
}
////////////////////////////////////////////////////////////////////////////
//   Explicit instantiation of template functions we need
////////////////////////////////////////////////////////////////////////////

template void
xml_read_from_file<ArrayOfMatrix> (const String&, ArrayOfMatrix&);

template void
xml_read_from_file<Index> (const String&, Index&);

template void
xml_read_from_file<Matrix> (const String&, Matrix&);

template void
xml_read_from_file<Numeric> (const String&, Numeric&);

template void
xml_read_from_file<Vector> (const String&, Vector&);


template void
xml_write_to_file<ArrayOfMatrix> (const String&, const ArrayOfMatrix&);

template void
xml_write_to_file<Index> (const String&, const Index&);

template void
xml_write_to_file<Matrix> (const String&, const Matrix&);

template void
xml_write_to_file<Numeric> (const String&, const Numeric&);

template void
xml_write_to_file<Vector> (const String&, const Vector&);
