/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

   This file contains basic functions to handle ASCII and binary (HDF)
   data files.

   \author Patrick Eriksson
   \date 2000-10-28 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "messages.h"
#include "file.h"
#include <hdf5.h>



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
void open_output_file(ofstream& file, const string& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // FIXME: This does not yet work in  egcs-2.91.66, try again later.
  file.exceptions(ios::badbit |
		  ios::failbit);
  
  // c_str explicitly converts to c string.
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
void open_input_file(ifstream& file, const string& name)
{
  // Tell the stream that it should throw exceptions.
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  file.exceptions(ios::badbit |
		  ios::failbit);

  // c_str explicitly converts to c string.
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
   Read an ASCII stream and append the contents to the string array
   text.  TEXT IS NOT OVERWRITTEN, BUT APPENDED!

   @param text Output. The contents fo the file
   @param is Stream from which to read
   @exception IOError Some error occured during the read
   @version   1
   @author Stefan Buehler */
void read_text_from_stream(ARRAY<string>& text, istream& is)
{
  string linebuffer;

  // Read as long as `is' is good.
  // Contary to what I understood from the book, the explicit check
  // for eof is necessary here, otherwise the last line is read twice
  // if it is not terminated by a newline character!
  while (is && !is.eof())
    {
      // Read line from file into linebuffer:
      getline(is,linebuffer);

      //      cout << "lb:" << linebuffer << '\n';
      
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
   Reads an ASCII file and appends the contents to the string vector
   text. This uses the function @see read_text_from_stream. TEXT IS
   NOT OVERWRITTEN, BUT APPENDED!  

   @param text Output. The contents fo the file
   @param  filename Name of file to read
   @exception IOError
   @version   1
   @author Stefan Buehler */
void read_text_from_file(ARRAY<string>& text, const string& name)
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

    @param s Output. The string to act on.
    @param what The string to replace.
    @param with The replacement.
    
    @author Stefan Buehler */
void replace_all(string& s, const string& what, const string& with)
{
  string::size_type j = s.find(what);
  while ( j != string::npos )
    {
      //		cout << "j = " << j << '\n';
      s.replace(j,1,with);
      j = s.find(what,j+with.size());
    }
}



////////////////////////////////////////////////////////////////////////////
//   MATRIX/VECTOR IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

//// write_array_of_matrix_to_stream ///////////////////////////////////////
/** A helper function that writes an array of matrix to a stream. This
    is the generic output function for VECTORs, MATRIXs, and ARRAYof
    both. All these are converted first to ARRAYofMATRIX, and then
    written by this function.

    @param os   Output. The stream to write to.
    @param am    The matrix to write.

    @author Stefan Buehler */
void write_array_of_matrix_to_stream(ostream& os,
                                     const ARRAYofMATRIX& am)
{
  extern const string full_name;

  // Determine the precision, depending on whether Numeric is double
  // or float:  
  int precision;
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
  const size_t n = am.dim();
  os << n << '\n';

  for (size_t i=1; i<=n; ++i)
    {
      // Number of elements:
      os << am(i).dim(1) << ' ' << am(i).dim(2) << '\n';

      os << setprecision(precision);
      // Write the elements:
      for (size_t r=1; r<=am(i).dim(1); ++r)
        {
          os << am(i)(r,1);
      
          for (size_t c=2; c<=am(i).dim(2); ++c)
            {
              os << " " << am(i)(r,c);
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
void write_array_of_matrix_to_file(const string& filename,
                                   const ARRAYofMATRIX& am)
{
  ofstream of;

  out2 << "  Writing " << filename << '\n';
  open_output_file(of, filename);

  // Write the array of matrix to the stream:
  write_array_of_matrix_to_stream(of,am);
}



//// read_array_of_matrix_from_stream ///////////////////////////////////////
/** A helper function that reads an array of matrix from a stream.

    @param am   Output. The array of matrix to read.
    @param is   Output. The input stream.

    @author Stefan Buehler */
void read_array_of_matrix_from_stream(ARRAYofMATRIX& am,
                                      istream& is)
{
  // First, skip all the lines that have a # at the beginning. (Maybe
  // preceded by whitespace.)
  bool comments=true;
  char c;
  string linebuffer;
  while (comments)
    {
      is >> ws;
      is.get(c);
      if ('#'==c)
        {
          getline(is,linebuffer);
          //      cout << "C: " << linebuffer << endl;
        }
      else
        {
          is.unget();
          comments = false;
        }
    }

  // Read the array of matrix. The TNT package expects exactly this input
  // format.
  // For ARRAY: First one number indicating the dimension, then the
  // elements. 
  // For MATRIX: First two numbers indicating the dimensions, then the
  // elements. 
  is >> am;

  if ( is.fail() || is.bad() )
    throw runtime_error("Stream gave fail or bad.");

  is >> ws;

  if ( !is.eof() )
    throw runtime_error("Input finished, but end of stream not reached.");

  // Some output to the lowest priority stream:
  out3 << "  Dimensions: "
       << am.dim() << ", ";
  for ( size_t i=0; i<am.size(); ++i )
    {
       out3 << am[i].dim(1) << ", "
	    << am[i].dim(2) << "\n";
    }
}



//// read_array_of_matrix_from_file /////////////////////////////////////////
/** A helper function that reads an array of matrix from a file. 
    Uses read_array_of_matrix_from_stream.

    @param am        Output. The array of matrix to read.
    @param filename  The name of the file to read.

    @author Stefan Buehler */
void read_array_of_matrix_from_file(ARRAYofMATRIX& am,
                                    const string& filename)
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
//   Basic functions to handle HDF files
////////////////////////////////////////////////////////////////////////////

void binfile_open(
              hid_t&    fid,
        const string&   name )
{
  fid = H5Fcreate( name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  // Handle errors
}


void binfile_close(
              hid_t&    fid,
        const string&   name )
{
  herr_t  status;

  status = H5Fclose( fid );
  // Handle errors
}


void binfile_write_matrix(
        const hid_t&    fid,
        const MATRIX&   a )
{
  hid_t    space_id, set_id;
  hsize_t  dims[2];
  herr_t   status;

  // Create test data
  double   b[2][3];
  int      i,j;
  for ( i=0; i<2; i++ )
    for ( j=0; j<3; j++ )
      b[i][j] = 1.1;

  // Create data space
  //dims[0]  = a.dim(1);  
  //dims[1]  = a.dim(2);  
  dims[0]  = 2;  
  dims[1]  = 3;  
  space_id = H5Screate_simple( 2, dims, NULL );

  // Create the data set
  set_id = H5Dcreate( fid, "/matrix", H5T_IEEE_F64LE, space_id, H5P_DEFAULT);


  // Write the data
  status = H5Dwrite( set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                                                          H5P_DEFAULT, b );

  // Close
  status = H5Dclose( set_id );
  status = H5Sclose( space_id );
}

