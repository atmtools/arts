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

string get_numeric_type_string( )
{
  if ( sizeof(Numeric) == sizeof(float) )
    return "float";
  else if ( sizeof(Numeric) == sizeof(double) )
    return "double";
  else
    throw runtime_error("Numeric must be double or float.");
}


void binfile_open_out(
              hid_t&    fid,
        const string&   filename )
{
  fid = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if ( fid < 0 )
  {
    ostringstream os;
    os << "Cannot create output file: " << filename << '\n'
       << "Maybe you don't have write access to the directory or the file?";
    throw runtime_error(os.str());
  }
}


void binfile_open_in(
              hid_t&    fid,
        const string&   filename )
{
  fid = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
  if ( fid < 0 )
  {
    ostringstream os;
    os << "Cannot open input file: " << filename << '\n'
       << "Maybe you don't have read access to the directory or the file?";
    throw runtime_error(os.str());
  }
}


void binfile_close(
              hid_t&    fid,
        const string&   filename )
{
  herr_t  status;

  status = H5Fclose( fid );
  if ( status < 0 )
  {
    ostringstream os;
    os << "Cannot close file: " << filename;
    throw runtime_error(os.str());
  }
}


void binfile_get_dataids(
              hid_t&    set_id,
              hid_t&    space_id,
        const hid_t&    fid,
        const string&   filename,
        const string&   dataname )
{
  const string   s = "\\" + dataname;
 
  // Open data set
  set_id = H5Dopen( fid, s.c_str() );

  // Handle errors
  if ( set_id < 0 )
  {
    ostringstream os;
    os << "Cannot open set: "<<dataname<<" in file: "<<filename<<"\n"
       << "Maybe the file contains data of other type?";
    throw runtime_error(os.str());
  }

  // Open data space
  space_id = H5Dget_space( set_id );
  // Handle errors
  if ( space_id < 0 )
  {
    ostringstream os;
    os << "Cannot open data space in file: " << filename;
    throw runtime_error(os.str());
  }
}


void binfile_close_set(
              hid_t&    set_id,
        const string&   filename )
{
  herr_t status = H5Dclose( set_id );
  if ( status < 0 )
  {
    ostringstream os;
    os << "Cannot close data set in file: " << filename;
    throw runtime_error(os.str());
  }
}


void binfile_close_space(
              hid_t&    space_id,
        const string&   filename )
{
  herr_t status = H5Sclose( space_id );
  if ( status < 0 )
  {
    ostringstream os;
    os << "Cannot close data space in file: " << filename;
    throw runtime_error(os.str());
  }
}


void binfile_get_rank(
              size_t&   rank,
        const string&   filename,
        const hid_t&    set_id,
        const hid_t&    space_id )
{
  rank     = H5Sget_simple_extent_ndims( space_id );
}


void binfile_get_size(
              hsize_t*  dims,
        const string&   filename,
        const hid_t&    set_id,
        const hid_t&    space_id )
{
  hsize_t maxdims[5];
  H5Sget_simple_extent_dims( space_id, dims, maxdims );
}


void binfile_read_init(
              hid_t&    set_id, 
              size_t&   rank,
              hsize_t*  dims,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname,
        const size_t&   rank0,
        const hsize_t*  dims0 )
{
  hid_t  space_id;

  // Get identifier for data set and space 
  binfile_get_dataids( set_id, space_id, fid, filename, dataname );

  // Get and check data rank
  binfile_get_rank( rank, filename, set_id, space_id );
  if ( rank != rank0 )
  {
    ostringstream os;
    os << "The opened data set (" << dataname << ") has rank " << rank << ",\n"
       << "but a rank of " << rank0 << " was expected.";
    throw runtime_error(os.str());
  }

  // Get and check data size
  binfile_get_size( dims, filename, set_id, space_id );
  if ( dims0 != NULL )
  {
    for ( size_t i=0; i<rank; i++ )
    {
      if ( (dims0[i]>0) && (dims[i]!=dims0[i]) )
      {
	ostringstream os;
	os << "For dimension " << i+1 << "the data size is " << dims[i] <<"\n"
	   << "but the expected length is " << dims0[i];
	throw runtime_error(os.str());
      }
    }
  }

  // Close data space
  binfile_close_space( space_id, filename );
}


////////////////////////////////////////////////////////////////////////////
//   Core read and write functions for binary files
////////////////////////////////////////////////////////////////////////////

void binfile_write(
        const hid_t&    fid,
        const string&   filename,
        const string&   dataname,
        const string&   datatype,
        const size_t&   rank,
        const hsize_t*  dims,
        const void*     dpointer )
{
        hid_t    space_id, set_id;
        herr_t   status;
  const string   s = "\\" + dataname;

  // Create the data space
  space_id = H5Screate_simple( rank, dims, NULL );
  if ( space_id < 0 )
  {
    ostringstream os;
    os << "Cannot create data space in file: " << filename << '\n'
       << "Is the file opened for writing?";
    throw runtime_error(os.str());
  }

  // Create the data set and write data
  if ( datatype == "size_t" ) 
  {
    if ( sizeof(size_t) == sizeof(unsigned) )
    {
      set_id = H5Dcreate( fid, s.c_str(), H5T_STD_I32BE, space_id, 
                                                                H5P_DEFAULT );
      status = H5Dwrite( set_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );
      cout << "a\n";
    }
    else if ( sizeof(size_t) == sizeof(unsigned long) )
    { 
      set_id = H5Dcreate( fid, s.c_str(), H5T_STD_I64BE, space_id, 
                                                                H5P_DEFAULT );
      status = H5Dwrite( set_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );
      cout << "b\n";
    }
    else
      throw runtime_error("Unknown data size of size_t.");
  }

  else if ( datatype == "Numeric" ) 
  {
    // Determine if NUMERIC is float or double
    string numtype = get_numeric_type_string();

    if ( numtype == "float" )
    {
      set_id = H5Dcreate( fid, s.c_str(), H5T_IEEE_F32LE, space_id, 
                                                                 H5P_DEFAULT);
      status = H5Dwrite( set_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );
    }
    else if ( numtype == "double" ) 
    {
      set_id = H5Dcreate( fid, s.c_str(), H5T_IEEE_F64LE, space_id, 
                                                                 H5P_DEFAULT);
      status = H5Dwrite( set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );
    }
    else
      throw runtime_error("Numeric must be double or float.");
  }

  else
  {
    ostringstream os;
    os << "The data type " << datatype << " is not handled";
    throw runtime_error(os.str());
  }

  // Handle errors
  if ( set_id < 0 )
  {
    ostringstream os;
    os << "Cannot create data set in file: " << filename << '\n'
       << "Is the file opened for writing?";
    throw runtime_error(os.str());
  }
  //
  if ( status < 0 )
  {
    ostringstream os;
    os << "Cannot write data to file: " << filename;
    throw runtime_error(os.str());
  }

  // Close data set and space
  binfile_close_space( space_id, filename );
  binfile_close_set( set_id, filename );
}


void binfile_read(
              void*     dpointer,
              hid_t&    set_id,
        const string&   filename,
        const string&   datatype,
        const bool&     close_set )
{
  // Read the data
  herr_t status;
  if ( datatype == "size_t" ) 
    status = H5Dread( set_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );

  else if ( datatype == "Numeric" ) 
  {
    // Determine if NUMERIC is float or double
    string numtype = get_numeric_type_string();

    if ( numtype == "float" )
      status = H5Dread( set_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );
    else if ( numtype == "double" ) 
      status = H5Dread( set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                                                      H5P_DEFAULT, dpointer );
    else
      throw runtime_error("Numeric must be double or float.");
  }

  else
  {
    ostringstream os;
    os << "The data type " << datatype << " is not handled";
    throw runtime_error(os.str());
  }

  // Handle errors
  if ( status < 0 )
  {
    ostringstream os;
    os << "Cannot read data from file: " << filename << "\n"
       << "Maybe the file data have another data format than assumed";;
    throw runtime_error(os.str());
  }

  // Close data set?
  if ( close_set )
    binfile_close_set( set_id, filename );
}



////////////////////////////////////////////////////////////////////////////
//   Functions to read and write binary data for ARTS data types
////////////////////////////////////////////////////////////////////////////

void binfile_write_size_t(
        const string&   filename,
        const hid_t&    fid,
        const size_t&   x,
        const string&   dataname )
{
  size_t  a[1];
  a[0] = x;
  hsize_t  dims[1];
  dims[0] = 1;
  binfile_write( fid,  filename, dataname, "size_t", 1, dims, a );
}


void binfile_read_size_t(
              size_t&   x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname )
{
  // We expect rank 1 and a length of 1
  const size_t   rank0 = 1;
        hsize_t  dims0[rank0];  dims0[0] = 1;
  //-------------------------------------------
        hsize_t  dims[rank0];
        hid_t    set_id;
        size_t   rank, a[rank0];

  binfile_read_init( set_id, rank, dims, filename, fid, dataname, rank0,dims0);
  binfile_read( a, set_id, filename, "size_t", 1 );
  x = a[0];  
}


void binfile_write_numeric(
        const string&   filename,
        const hid_t&    fid,
        const Numeric&  x,
        const string&   dataname )
{
  Numeric  a[1];
  a[0] = x;
  hsize_t  dims[1];
  dims[0] = 1;
  binfile_write( fid,  filename, dataname, "Numeric", 1, dims, a );
}


void binfile_read_numeric(
              Numeric&  x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname )
{
  // We expect rank 1 and a length of 1
  const size_t   rank0 = 1;
        hsize_t  dims0[rank0];  dims0[0] = 1;
  //-------------------------------------------
        hsize_t  dims[rank0];
        hid_t    set_id;
        size_t   rank;
        Numeric  a[rank0];

  binfile_read_init( set_id, rank, dims, filename, fid, dataname, rank0,dims0);
  binfile_read( a, set_id, filename, "Numeric", 1 );
  x = a[0];  
}


void binfile_write_vector(
        const string&   filename,
        const hid_t&    fid,
        const VECTOR&   x,
        const string&   dataname )
{
  const size_t  n = x.dim();
        size_t  i;
       hsize_t  dims[1];
  dims[0] = n;

  Numeric a[n];
  for ( i=0; i<n; i++ )
    a[i] = x(i+1);
  binfile_write( fid,  filename, dataname, "Numeric", 1, dims, a );
}


void binfile_read_vector(
              VECTOR&   x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname )
{
  // We expect rank 1 (any length)
  const size_t   rank0 = 1;
  //-------------------------------------------
        hsize_t  dims[rank0];
        size_t   rank, n, i;
        hid_t    set_id;

  binfile_read_init( set_id, rank, dims, filename, fid, dataname, rank0, NULL);
  n = dims[0];
  x.newsize(n);
  Numeric a[n];
  binfile_read( a, set_id, filename, "Numeric", 1 );
  for ( i=0; i<n; i++ )
    x(i+1) = a[i];
}


void binfile_write_matrix(
        const string&   filename,
        const hid_t&    fid,
        const MATRIX&   x,
        const string&   dataname )
{
  const size_t  nrows = x.dim(1);
  const size_t  ncols = x.dim(2);
        size_t  i, j;
       hsize_t  dims[2];
  dims[0] = nrows;
  dims[1] = ncols;

  Numeric a[nrows][ncols];
  for ( i=0; i<nrows; i++ )
    for ( j=0; j<ncols; j++ )
      a[i][j] = x(i+1,j+1);
  binfile_write( fid,  filename, dataname, "Numeric", 2, dims, a );
}


void binfile_read_matrix(
              MATRIX&   x,
        const string&   filename,
        const hid_t&    fid,
        const string&   dataname )
{
  // We expect rank 2 (any size)
  const size_t   rank0 = 2;
  //-------------------------------------------
        hsize_t  dims[rank0];
        size_t   rank, nrows, ncols, i, j;
        hid_t    set_id;

  binfile_read_init( set_id, rank, dims, filename, fid, dataname, rank0, NULL);
  nrows = dims[0];
  ncols = dims[1];
  x.newsize(nrows,ncols);
  Numeric a[nrows][ncols];
  binfile_read( a, set_id, filename, "Numeric", 1 );
  for ( i=0; i<nrows; i++ )
    for ( j=0; j<ncols; j++ )
      x(i+1,j+1) = a[i][j];
}
       
