/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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
/**
   \file  file.cc

   This file contains basic functions to handle ASCII data files.

   \author Patrick Eriksson
   \date 2000-10-28 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <stdexcept>
#include <cmath>
#include <cfloat>
#include "arts.h"
#include "matpackI.h"
#include "array.h"
#include "messages.h"
#include "file.h"



////////////////////////////////////////////////////////////////////////////
//   Default file names
////////////////////////////////////////////////////////////////////////////

//// filename_ascii ////////////////////////////////////////////////////////
/**
   Gives the default file name for the ASCII formats.

   The default name is only used if the file name is empty.

   \param   filename Output:     file name
   \param    varname      variable name

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void filename_ascii(
              String&  filename,
        const String&  varname )
{
  if ( "" == filename )
  {
    extern const String out_basename;                       
    filename = out_basename+"."+varname+".aa";
  }
}



////////////////////////////////////////////////////////////////////////////
//   Functions to open and read ASCII files
////////////////////////////////////////////////////////////////////////////

//// open_output_file //////////////////////////////////////////////////////
/**
   Open a file for writing. If the file cannot be opened, the
   exception IOError is thrown. 
   @param     file File pointer 
   @param     name Name of the file to open
   @author    Stefan Buehler
   @version   1
   @exception ios_base::failure Could for example mean that the
                      directory is read only. */
void open_output_file(ofstream& file, const String& name)
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



//// open_input_file ///////////////////////////////////////////////////////
/**
   Open a file for reading. If the file cannot be opened, the
   exception IOError is thrown. 
   @param     file File pointer 
   @param     name Name of the file to open
   @author    Stefan Buehler
   @version   1
   @exception ios_base::failure Somehow the file cannot be opened. */
void open_input_file(ifstream& file, const String& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted.
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  file.exceptions(ios::badbit);

  // c_str explicitly converts to c String.
  file.open(name.c_str() );

  // See if the file is ok.
  // FIXME: This should not be necessary anymore in the future, when
  // g++ stream exceptions work properly.
  if (!file)
    {
      ostringstream os;
      os << "Cannot open input file: " << name << '\n'
         << "Maybe the file does not exist?";
      throw runtime_error(os.str());
    }
}



//// read_text_from_stream /////////////////////////////////////////////////
/**
   Read an ASCII stream and append the contents to the String array
   text.  TEXT IS NOT OVERWRITTEN, BUT APPENDED!

   @param text Output. The contents fo the file
   @param is Stream from which to read
   @exception IOError Some error occured during the read
   @version   1
   @author Stefan Buehler */
void read_text_from_stream(ArrayOfString& text, istream& is)
{
  String linebuffer;

  // Read as long as `is' is good.
  // Contary to what I understood from the book, the explicit check
  // for eof is necessary here, otherwise the last line is read twice
  // if it is not terminated by a newline character!
  while (is && is.good() && !is.eof())
    {
      // Read line from file into linebuffer:
      getline(is,linebuffer);
      
      // Append to end of text:
      text.push_back(linebuffer);
    }
  
  // Check for error:
  // FIXME: This should not be necessary anymore when stream
  // exceptions work properly.
  if ( !is.eof() ) {
    ostringstream os;
    os << "Read Error. Last line read:\n" << linebuffer;
    throw runtime_error(os.str());
  }

}



//// read_text_from_file ////////////////////////////////////////////////////
/**
   Reads an ASCII file and appends the contents to the String vector
   text. This uses the function @see read_text_from_stream. TEXT IS
   NOT OVERWRITTEN, BUT APPENDED!  

   @param text Output. The contents fo the file
   @param  filename Name of file to read
   @exception IOError
   @version   1
   @author Stefan Buehler */
void read_text_from_file(ArrayOfString& text, const String& name)
{
  ifstream ifs;

  // Open input stream:
  open_input_file(ifs, name);
  // No need to check for error, because open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the text from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the 
  // filename.
  try
    {
      read_text_from_stream(text,ifs);
    }
  catch (runtime_error x)
    {
      ostringstream os;
      os << "Error reading file: " << name << '\n'
         << x.what();
      throw runtime_error(os.str());
    }
}



//// replace_all //////////////////////////////////////////////////////////
/** 
    Replace all occurances of `what' in `s' with `with'.

    @param s Output. The String to act on.
    @param what The String to replace.
    @param with The replacement.
    
    @author Stefan Buehler */
void replace_all(String& s, const String& what, const String& with)
{
  Index j = s.find(what);
  while ( j != s.npos )
    {
      s.replace(j,1,with);
      j = s.find(what,j+with.size());
    }
}



////////////////////////////////////////////////////////////////////////////
//   Matrix/Vector IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

//// write_array_of_matrix_to_stream ///////////////////////////////////////
/** A helper function that writes an array of matrix to a stream. This
    is the generic output function for Vectors, Matrixs, and ArrayOf
    both. All these are converted first to ArrayOfMatrix, and then
    written by this function.

    @param os   Output. The stream to write to.
    @param am    The matrix to write.

    @author Stefan Buehler */
void write_array_of_matrix_to_stream(ostream& os,
                                     const ArrayOfMatrix& am)
{
  extern const String full_name;

  // Determine the precision, depending on whether Numeric is double
  // or float:  
  Index precision;
  switch (sizeof(Numeric)) {
  case sizeof(float)  : precision = FLT_DIG; break;
  case sizeof(double) : precision = DBL_DIG; break;
  default: out0 << "Numeric must be double or float\n"; exit(1);
  }

  os << "# Generated by "
     << full_name << ", "
     << __DATE__ << ", "
     << __TIME__ << "\n";

  // Number of array elements:
  const Index n = am.nelem();
  os << n << '\n';

  for (Index i=0; i<n; ++i)
    {
      // Number of elements:
      os << am[i].nrows() << ' ' << am[i].ncols() << '\n';

      os << setprecision(precision);
      // Write the elements:
      for (Index r=0; r<am[i].nrows(); ++r)
        {
          os << am[i](r,0);
      
          for (Index c=1; c<am[i].ncols(); ++c)
            {
              os << " " << am[i](r,c);
            }

          os << '\n';
        }
    }
}



//// write_array_of_matrix_to_file /////////////////////////////////////////
/** A helper function that writes an array of matrix to a file. Uses
    write_array_of_matrix_to_stream. 

    @param filename    The name of the file.
    @param am          The array of matrix to write.

    @author Stefan Buehler */
void write_array_of_matrix_to_file(const String& filename,
                                   const ArrayOfMatrix& am)
{
  ofstream of;

  out2 << "  Writing " << filename << '\n';
  open_output_file(of, filename);

  // Write the array of matrix to the stream:
  write_array_of_matrix_to_stream(of,am);
}

/** A helper function that skips lines containing comments. Comments
    are all lines that have a # at the beginning. (Maybe preceded by
    whitespace. 

    \param is The input stream.

    \author Stefan Buehler
    \date   2001-11-09
*/
void skip_comments(istream & is)
{
  bool comments=true;
  char c;
  String linebuffer;
  while (comments)
    {
      is >> ws;
      is.get(c);
      if ('#'==c)
          getline(is,linebuffer);
      else
        {
          is.unget();
          comments = false;
        }
    }
}

//// read_array_of_matrix_from_stream ///////////////////////////////////////
/** A helper function that reads an array of matrix from a stream.

    @param am   Output. The array of matrix to read.
    @param is   Output. The input stream.

    @author Stefan Buehler */
void read_array_of_matrix_from_stream(ArrayOfMatrix& am,
                                      istream& is)
{
  Index n;
  skip_comments(is);   // Skip comment lines.
  is >> n;
  //  cout << "n: " << n << "\n";
  am.resize(n);
  for( Index i=0; i<n; ++i )
    {
      //      cout << "- " << i << "\n";
      {
        Index nr,nc;
        skip_comments(is);   // Skip comment lines.
        is >> nr >> nc;
        //   cout << "nr: " << nr << "\n";
        //   cout << "nc: " << nc << "\n";
        am[i].resize(nr,nc);
        //        cout << "am[i]: " << am[i].nrows() << " / " << am[i].ncols() << "\n";
        for(Index ir=0; ir<nr; ++ir)
          for(Index ic=0; ic<nc; ++ic)
            {
              //                cout << "ir / ic = " << ir << " / " << ic << "\n";
              skip_comments(is);   // Skip comment lines.
              is >> am[i](ir,ic);
            }
      }
    }

  if ( is.fail() || is.bad() )
    throw runtime_error("Stream gave fail or bad.");

  is >> ws;

  if ( !is.eof() )
    throw runtime_error("Input finished, but end of stream not reached.");

  // Some output to the lowest priority stream:
  out3 << "  Dimensions:\n";
  out3 << "     "<< am.nelem() << "\n";
  for ( Index i=0; i<am.nelem(); ++i )
    out3 << "     " << am[i].nrows() << ", " << am[i].ncols() << "\n";
}



//// read_array_of_matrix_from_file /////////////////////////////////////////
/** A helper function that reads an array of matrix from a file. 
    Uses read_array_of_matrix_from_stream.

    @param am        Output. The array of matrix to read.
    @param filename  The name of the file to read.

    @author Stefan Buehler */
void read_array_of_matrix_from_file(ArrayOfMatrix& am,
                                    const String& filename)
{
  ifstream ifs;

  out2 << "  Reading " << filename << '\n';
  
  // Open input stream:
  open_input_file(ifs, filename);
  // No need to check for error, because open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the matrix from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the 
  // filename.
  try
    {
      read_array_of_matrix_from_stream(am,ifs);
    }
  catch (runtime_error x)
    {
      ostringstream os;
      os << "Error reading file: " << filename << '\n'
         << x.what();
      throw runtime_error(os.str());
    }
}



////////////////////////////////////////////////////////////////////////////
//   STRING IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

//// write_array_of_string_to_stream ///////////////////////////////////////
/**
   A helper function that writes an array of String to a stream. 

   \param   os Output:           the stream to write to
   \param    as           the String array to write

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void write_array_of_String_to_stream(
              ostream&         os,
        const ArrayOfString&   as )
{
  extern const String full_name;

  os << "# Generated by " << full_name << "\n";

  // Number of array elements:
  const Index n = as.nelem();
  os << n << '\n';

  for (Index i=0; i<n; ++i)
    os << as[i] << '\n';
}



//// write_array_of_String_to_file /////////////////////////////////////////
/**
   A help function that writes an array of String to a file.

   \param filename Output:    The name of the file.  
   \param  as          The array of String to write.  

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void write_array_of_String_to_file(
        const String&          filename,
        const ArrayOfString&   as )
{
  ofstream of;

  out2 << "  Writing " << filename << '\n';
  open_output_file(of, filename);

  // Write the array of String to the stream:
  write_array_of_String_to_stream(of,as);
}



//// read_array_of_String_from_stream ///////////////////////////////////////
/**
   A help function to read an array of String from a stream.

   \param as Output:   The array of String to read
   \param  is   The input stream

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void read_array_of_String_from_stream(
        ArrayOfString&   as,
        istream&         is )
{
  // First, skip all the lines that have a # at the beginning. (Maybe
  // preceded by whitespace.)
  bool comments=true;
  char c;
  String linebuffer;
  while (comments)
    {
      is >> ws;
      is.get(c);
      if ('#'==c)
          getline(is,linebuffer);

      else
        {
          is.unget();
          comments = false;
        }
    }

  // Read the Array of Strings from the input stream:
  {
    Index n;
    is >> n;
    //  cout << "n: " << n << "\n";
    as.resize(n);
    for( Index i=0; i<n; ++i )
      {
        //      cout << "- " << i << "\n";
        is >> as[i];
      }
  }

  if ( is.fail() || is.bad() )
    throw runtime_error("Stream gave fail or bad.");

  is >> ws;
  if ( !is.eof() )
    throw runtime_error("Input finished, but end of stream not reached.");

  // Some output to the lowest priority stream:
  out3 << "  Dimension: "
       << as.nelem() << ", ";
}



//// read_array_of_matrix_from_file /////////////////////////////////////////
/**
   A help function to read an array of String from a file.

   \param as Output:        The array of String to read
   \param  filename  The name of the file to read

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void read_array_of_String_from_file(
           ArrayOfString&   as,
     const String&          filename )
{
  ifstream ifs;

  out2 << "  Reading " << filename << '\n';
  
  // Open input stream:
  open_input_file(ifs, filename);
  // No need to check for error, because open_input_file throws a
  // runtime_error with an appropriate error message.

  // Read the matrix from the stream. Here we catch the exception,
  // because then we can issue a nicer error message that includes the 
  // filename.
  try
    {
      read_array_of_String_from_stream(as,ifs);
    }
  catch (runtime_error x)
    {
      ostringstream os;
      os << "Error reading file: " << filename << '\n'
         << x.what();
      throw runtime_error(os.str());
    }
}



