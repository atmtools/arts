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

   This file contains basic functions to handle ASCII and binary (HDF)
   data files.

   \author Patrick Eriksson
   \date 2000-10-28 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <hdf.h>
#include <stdexcept>
#include <math.h>
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

   \retval   filename     file name
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
    extern const String basename;                       
    filename = basename+"."+varname+".aa";
  }
}



//// filename_bin ////////////////////////////////////////////////////////////
/**
   Gives the default file name for the binary format.

   The default name is only used if the file name is empty.

   \retval   filename     file name
   \param    varname      variable name

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void filename_bin(
              String&  filename,
        const String&  varname )
{
  if ( "" == filename )
  {
    extern const String basename;                       
    filename = basename+"."+varname+".ab";
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
  // Badbit means that the entire stream is corrupted, failbit means
  // that the last operation has failed, but the stream is still
  // valid. We don't want either to happen!
  // On the other hand, end of file will not lead to an exception, you
  // have to check this manually!
  file.exceptions(ios::badbit |
		  ios::failbit);

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
  while (is && !is.eof())
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



//// read_array_of_matrix_from_stream ///////////////////////////////////////
/** A helper function that reads an array of matrix from a stream.

    @param am   Output. The array of matrix to read.
    @param is   Output. The input stream.

    @author Stefan Buehler */
void read_array_of_matrix_from_stream(ArrayOfMatrix& am,
                                      istream& is)
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

  // Read the Array of Matrices. 
  {
    Index n;
    is >> n;
    //  cout << "n: " << n << "\n";
    am.resize(n);
    for( Index i=0; i<n; ++i )
      {
	//      cout << "- " << i << "\n";
	{
	  Index nr,nc;
	  is >> nr >> nc;
	  //   cout << "nr: " << nr << "\n";
	  //   cout << "nc: " << nc << "\n";
	  am[i].resize(nr,nc);
	  //	  cout << "am[i]: " << am[i].nrows() << " / " << am[i].ncols() << "\n";
	  for(Index ir=0; ir<nr; ++ir)
	    for(Index ic=0; ic<nc; ++ic)
	      {
		//		cout << "ir / ic = " << ir << " / " << ic << "\n";
		is >> am[i](ir,ic);
	      }
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

   \retval   os           the stream to write to
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

   \retval filename    The name of the file.  
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

   \retval as   The array of String to read
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

   \retval as        The array of String to read
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



////////////////////////////////////////////////////////////////////////////
//   Help functions to handle HDF files
////////////////////////////////////////////////////////////////////////////

//// check_data_types //////////////////////////////////////////////////////
/**
   Checks that some data types have the expected length.

   The following data types are checked: Index, float and double.

   @exception runtime_error  Some data type has an unexpected length.

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void check_data_types()
{
  if ( sizeof(Index) != 4 )
    throw runtime_error("An Index is expected to be 4 bytes.");
  if ( sizeof(float) != 4 )
    throw runtime_error("A float is expected to be 4 bytes.");
  if ( sizeof(double) != 8 )
    throw runtime_error("A double is expected to be 8 bytes.");
  if ( sizeof(char) != 1 )
    throw runtime_error("A char is expected to be 1 byte.");
}



//// binfile_open_out //////////////////////////////////////////////////////
/**
   Opens a binary file for writing.

   \retval   fid          file identifier
   \param    filename     file name

   @exception runtime_error  The file cannot be opened.
   @exception runtime_error  Cannot init the VS interface.

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_open_out(
              int&      fid,
        const String&   filename )
{
  // Open the file for writing (deleting an old file with same name)
  fid = Hopen( filename.c_str(), DFACC_CREATE, 0 );
  if ( fid < 0 )
  {
    ostringstream os;
    os << "Cannot create output file: " << filename << '\n'
       << "Maybe you don't have write access to the directory or the file?";
    throw runtime_error(os.str());
  }
  out2 << "  Opened file " << filename << " for writing.\n";

  // Initialize the VS interface
  if ( Vstart( fid ) < 0 )
  {
    ostringstream os;
    os << "Cannot initialize the VS interafce in file: " << filename;
    throw runtime_error(os.str());
  }
}



//// binfile_open_in //////////////////////////////////////////////////////
/**
   Opens a binary file for reading.

   \retval   fid          file identifier
   \param    filename     file name

   @exception runtime_error  The file is not a HDF file.
   @exception runtime_error  The file cannot be opened.
   @exception runtime_error  Cannot init the VS interface.

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_open_in(
              int&      fid,
        const String&   filename )
{
  // Check if the file is HDF
  if ( Hishdf( filename.c_str() ) < 0 )
  {
    ostringstream os;
    os << "The file " << filename << " is not a HDF file.\n";
    throw runtime_error(os.str());
  }

  // Open the file for reading
  fid = Hopen( filename.c_str(), DFACC_READ, 0 );
  if ( fid < 0 )
  {
    ostringstream os;
    os << "Cannot open input file: " << filename << '\n'
       << "Maybe you don't have read access to the directory or the file?";
    throw runtime_error(os.str());
  }
  out2 << "  Opened file " << filename << " for reading.\n";

  // Initialize the VS interface
  if ( Vstart( fid ) < 0 )
  {
    ostringstream os;
    os << "Cannot initialize the VS interafce in file: " << filename;
    throw runtime_error(os.str());
  }
}



//// binfile_close ////////////////////////////////////////////////////////
/**
   Closes a binary file.

   \retval   fid          file identifier
   \param    filename     file name

   @exception runtime_error  Cannot close the VS interface.
   @exception runtime_error  The file cannot be closed.

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_close(
              int&      fid,
        const String&   filename )
{
  // Terminate access to the VS interface
  if ( Vend( fid ) < 0 )
  {
    ostringstream os;
    os << "Cannot terminate access to the VS interface in: " << filename;
    throw runtime_error(os.str());
  }

  // Close the file
  if ( Hclose( fid ) < 0 )
  {
    ostringstream os;
    os << "Cannot close file: " << filename;
    throw runtime_error(os.str());
  }

  out2 << "  Closed file " << filename << "\n";
}



//// binfile_write_size /////////////////////////////////////////////////////
/**
   Writes the size of the data as an attribute to the vdata.

   The size is written as [nrows,ncols]

   \param    filename     file name
   \param    dataname     name on data set
   \param    vdata_id     vdata identifier
   \param    nrows        number of data rows
   \param    ncols        number of data columns

   @exception runtime_error  Cannot create attribute.

   \author Patrick Eriksson              
   \date   2000-11-09
*/
void binfile_write_size( 
        const String&   filename, 
        const String&   dataname,
        const int&      vdata_id,
        const Index&   nrows,
        const Index&   ncols )
{
  Index v[2];
  v[0] = nrows;
  v[1] = ncols;

  if ( VSsetattr( vdata_id, _HDF_VDATA, "SIZE", DFNT_UINT32, 2, v ) < 0 )
  {
    ostringstream os;
    os << "Cannot write size data for " << dataname << " in file " << filename;
    throw runtime_error(os.str());
  }
}



//// binfile_read_init /////////////////////////////////////////////////////
/**
   Initilizes a binary data field for reading.

   The function returns an identifier for the data and the data size. 
   The size is also compared to the expected size.

   \retval   vdata_id     identifier for the Vdata
   \retval   nrows        number of data rows
   \retval   ncols        number of data columns
   \param    fid          file identifier
   \param    filename     file name
   \param    dataname     name on the data
   \param    storagetype  SCALAR, Vector, Matrix etc.
   \param    nrows0       expected number of data rows
   \param    ncols0       expected number of data columns

   @exception runtime_error  Cannot find or access the data

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_read_init(
              int&      vdata_id,
              Index&   nrows,
              Index&   ncols,
        const int&      fid,
        const String&   filename,
        const String&   dataname,
        const String&   storagetype,
        const Index&   nrows0,
        const Index&   ncols0 )
{
  // Find the Vdata in the file
  int  vdata_ref = VSfind( fid, dataname.c_str() );
  if ( vdata_ref <= 0 )
  {
    ostringstream os;
    os << "Cannot find the data " << dataname << " in file " <<filename<<"\n"
       << "Maybe the file contains data of other type";
    throw runtime_error(os.str());
  }

  // Attach the Vdata
  vdata_id = VSattach( fid, vdata_ref, "r" );
  if ( vdata_id <= 0 )
  {
    ostringstream os;
    os << "Cannot attach the data " << dataname << " in file " << filename;
    throw runtime_error(os.str());
  }

  // Get number of rows and columns
  Index  v[2];
  if ( VSgetattr( vdata_id, _HDF_VDATA, 0, v ) < 0 )
  {
    ostringstream os;
    os << "Cannot determine the size of " << dataname << "\n" 
       << "in file " << filename;
    throw runtime_error(os.str());
  }
  nrows = v[0];
  ncols = v[1];

  // Check if number of rows and columns are as expected
  if ( (nrows0>0) && (nrows!=nrows0) )
  {
    ostringstream os;
    os << nrows0 << " rows were expected, but the data have " <<nrows<<" rows";
    throw runtime_error(os.str());    
  }
  if ( (ncols0>0) && (ncols!=ncols0) )
  {
    ostringstream os;
    os << ncols0 << " columns were expected, but the data have " << ncols
       << " columns";
    throw runtime_error(os.str());    
  }

  // Set fields to read
  if ( (nrows>0) && (ncols>0) )
  {
    if ( VSsetfields( vdata_id, storagetype.c_str() ) < 0 )
    {
      cout << dataname << endl;
      ostringstream os;
      os << "Cannot find the field " << storagetype << " in file " << filename
         << "\n" << "Maybe the file contains data of other type";
      throw runtime_error(os.str());
    }
  }
}



//// binfile_read_end /////////////////////////////////////////////////////
/**
   Closes a binary data field for reading.

   \retval   vdata_id     identifier for the Vdata
   \param    filename     file name
   \param    dataname     name on the data

   @exception runtime_error  Cannot detach the data

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_read_end(
              int&      vdata_id,
        const String&   filename,
        const String&   dataname )
{
  if ( VSdetach( vdata_id ) < 0 )
  {
    ostringstream os;
    os << "Cannot detach the field " << dataname << " in file " << filename;
    throw runtime_error(os.str());
  }
}



//// binfile_get_datatype ///////////////////////////////////////////////////
/**
   Gets the data type of the file data.

   \retval  type_in_file   data type String (e.g. "DOUBLE")
   \param   vdata_id       identifier for the Vdata

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_get_datatype(
              String&   type_in_file,
        const int&      vdata_id )
{
  char   c[50];
  VSgetclass( vdata_id, c );
  type_in_file = c;
}



////////////////////////////////////////////////////////////////////////////
//   Core read and write functions for binary files
////////////////////////////////////////////////////////////////////////////

//// binfile_write ////////////////////////////////////////////////////////
/**
   Core function for writing to binary files.

   The data can be a scalar, a vector or a matrix holding Index or NUMERIC.   

   \param    fid          file identifier
   \param    filename     file name
   \param    dataname     name on data set
   \param    storagetype  type of data container (e.g. Vector or Matrix)
   \param    atomictype   basic data type (Index or NUMERIC)
   \param    nrows        number of data rows
   \param    ncols        number of data columns
   \param    dpointer     pointer to the start of the data to store

   @exception runtime_error  Unknown atomic data type.

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write(
        const int&      fid,
        const String&   filename,
        const String&   dataname,
        const String&   storagetype,
        const String&   atomictype,
        const Index&   nrows,
        const Index&   ncols,
        const uint8*    dpointer )
{ 
  // Check that data types have expected length
  check_data_types();

  out3 << "    Writing: " << dataname << "\n";

  // Create a new vdata
  int  vdata_id = VSattach( fid, -1, "w" );
  if ( vdata_id < 0 )
  {
    ostringstream os;
    os << "Cannot create a new vdata in file " << filename;
    throw runtime_error(os.str());
  }
 
  // Set name of the vdata
  if ( VSsetname( vdata_id, dataname.c_str() ) < 0 )
  {
    ostringstream os;
    os << "Cannot name the vdata " << dataname << " in file " << filename;
    throw runtime_error(os.str());
  }

  // Write data size
  binfile_write_size( filename, dataname, vdata_id, nrows, ncols );

  // Create the field
  int    status1, status2;
  //
  if ( atomictype == "INDEX" ) 
  {
    status1 = VSsetclass( vdata_id, "UINT"  );
    status2 = VSfdefine( vdata_id, storagetype.c_str(), DFNT_UINT32, 1);
  }     
  else if ( atomictype == "NUMERIC" ) 
  {
    if ( sizeof(Numeric) == 4 )
    {
      status1 = VSsetclass( vdata_id, "FLOAT"  );
      status2 = VSfdefine( vdata_id, storagetype.c_str(), DFNT_FLOAT32, 1); 
    }
    else
    {
      status1 = VSsetclass( vdata_id, "DOUBLE"  );
      status2 = VSfdefine( vdata_id, storagetype.c_str(), DFNT_FLOAT64, 1);
    }
  }

  else if ( atomictype == "CHAR" ) 
  {
    status1 = VSsetclass( vdata_id, "CHAR"  );
    status2 = VSfdefine( vdata_id, storagetype.c_str(), DFNT_CHAR, 1);     
  }
  else
  {
    ostringstream os;
    os << "The atomic data type " << atomictype << " is not handled";
    throw runtime_error(os.str());
  }

  // Handle error
  if ( status1 < 0 )
  {
    ostringstream os;
    os << "Cannot set class on " << dataname << " in file "<<filename;
    throw runtime_error(os.str());
  }
  if ( status2 < 0 )
  {
    ostringstream os;
    os << "Cannot create the field " << storagetype << " in file "<<filename;
    throw runtime_error(os.str());
  }

  // Finalize the definition of the field
  if ( VSsetfields( vdata_id, storagetype.c_str() ) < 0 )
  {
    ostringstream os;
    os << "Cannot set the field " << storagetype << " in file " << filename;
    throw runtime_error(os.str());
  }

  // Write data (if not empty)
  if ( (nrows>0) && (ncols>0) )
  {
    // Do actual writing
    const Index   nout = nrows*ncols;
    Index ndone = VSwrite( vdata_id, dpointer, nout, FULL_INTERLACE );
    if ( ndone != nout )
    {
      ostringstream os;
      os << "Could not write all data to field " << storagetype << "in file "
         << filename << "\nOut of memory?";
      throw runtime_error(os.str());
    }
  }

  // Detach the vdata
  if ( VSdetach( vdata_id ) < 0 )
  {
    ostringstream os;
    os << "Cannot detach the vdata " << dataname << " in file " << filename;
    throw runtime_error(os.str());
  }
}



//// binfile_read1 ///////////////////////////////////////////////////////////
/**
   Core function for reading data of index type from binary files.

   The data is read an index array.

   The reading is splitted to several files as type conversions can be needed.

   \retval   x            the index array to read
   \param    vdata_id     data identifier
   \param    nrows        number of values
   \param    filename     file name
   \param    dataname     name on data set

   @exception runtime_error  Unknown data type in file.

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read1(
              ArrayOfIndex&   x,
        const int&            vdata_id,
        const Index&         n,
        const String&         filename,
        const String&         dataname )
{
  // Check that data types have expected length
  check_data_types();

  out3 << "    Reading: " << dataname << "\n";

  // Get the data type of the data to read
  String type_in_file;
  binfile_get_datatype( type_in_file, vdata_id );

  // Reallocate x
  x.resize(n);

  if ( n > 0 )
  {
    // Do reading and copy data
    if ( type_in_file == "UINT" )
      //
      VSread( vdata_id, (uint8*)&x[0], n, FULL_INTERLACE );
      /*
    {
      Index  a[n];
      VSread( vdata_id, (uint8*)a, n, FULL_INTERLACE );
      for ( Index i=0; i<n; i++ )
	x[i] = a[i];
    }
      */

    else
    {
      ostringstream os;
      os << "Files with data type " << type_in_file << " are not handled";
      throw runtime_error(os.str());
    }
  }
}



//// binfile_read2 ///////////////////////////////////////////////////////////
/**
   Core function for reading data of numeric type from binary files.

   The data is read a matrix.

   The reading is splitted to several files as type conversions can be needed.

   \retval   x            the matrix to read
   \param    vdata_id     data identifier
   \param    nrows        number of data rows
   \param    ncols        number of data columns
   \param    filename     file name
   \param    dataname     name on data set

   @exception runtime_error  Unknown data type in file.

   \author Patrick Eriksson              
   \date   2000-11-01

   Adapted to 0-based indexing.
   \date   2001-01-08
   \author Stefan Buehler
*/
void binfile_read2(
              Matrix&   x,
        const int&      vdata_id,
        const Index&   nrows,
        const Index&   ncols,
        const String&   filename,
        const String&   dataname )
{
  // Check that data types have expected length
  check_data_types();

  out3 << "    Reading: " << dataname << "\n";

  // Get the data type of the data to read
  String type_in_file;
  binfile_get_datatype( type_in_file, vdata_id );

  // Reallocate x
  x.resize(nrows,ncols);

  if ( (nrows > 0) && (ncols > 0) )
  {
    // Do reading and copy data
    if ( type_in_file == "FLOAT" )

      if ( sizeof(Numeric) == 4 )
        //
	VSread( vdata_id, (uint8*)&x(0,0), nrows*ncols, FULL_INTERLACE );

      else
      {
	float *a = new float[nrows*ncols];
	Index i,j,j0;
	VSread( vdata_id, (uint8*)a, nrows*ncols, FULL_INTERLACE );
	for ( i=0; i<nrows; i++ )
	{
		j0 = i*ncols;
		for ( j=0; j<ncols; j++ )
		  x(i,j) = a[j0+j];
	}
	delete[] a;
      }
  
    else if ( type_in_file == "DOUBLE" )

      if ( sizeof(Numeric) == 8 )
        //
    	VSread( vdata_id, (uint8*)&x(0,0), nrows*ncols, FULL_INTERLACE );

      else
      {
	double *a = new double[nrows*ncols];
	Index i,j,j0;
	VSread( vdata_id, (uint8*)a, nrows*ncols, FULL_INTERLACE );
	for ( i=0; i<nrows; i++ )
	{
		j0 = i*ncols;
		for ( j=0; j<ncols; j++ )
		   x(i,j) = a[j0+j];
	}
	delete[] a;
      }
  
    else
    {
      ostringstream os;
      os << "Files with data type " << type_in_file << " are not handled";
      throw runtime_error(os.str());
    }
  }
}



//// binfile_read3 ///////////////////////////////////////////////////////////
/**
   Core function for reading text from binary files.

   The data is read a String.

   The reading is splitted to several files as type conversions can be needed.

   \retval   x            the String to read
   \param    vdata_id     data identifier
   \param    nrows        number of values
   \param    filename     file name
   \param    dataname     name on data set

   @exception runtime_error  Unknown data type in file.

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_read3(
              String&   x,
        const int&      vdata_id,
        const Index&   n,
        const String&   filename,
        const String&   dataname )
{
  // Check that data types have expected length
  check_data_types();

  out3 << "    Reading: " << dataname << "\n";

  // Get the data type of the data to read
  String type_in_file;
  binfile_get_datatype( type_in_file, vdata_id );

  // Reallocate x
  x.resize(n);

  if ( n > 0 )
  {
    // Do reading and copy data
    if ( type_in_file == "CHAR" )
      //
      VSread( vdata_id, (uint8*)&x[0], n, FULL_INTERLACE );

    /*
    {
      char  a[n];
      VSread( vdata_id, (uint8*)a, n, FULL_INTERLACE );
      for ( Index i=0; i<n; i++ )
	x[i] = a[i];
    }
    */
    else
    {
      ostringstream os;
      os << "Files with data type " << type_in_file << " are not handled";
      throw runtime_error(os.str());
    }
  }
}



////////////////////////////////////////////////////////////////////////////
//   Functions to read and write binary data for ARTS data types
////////////////////////////////////////////////////////////////////////////

//// binfile_write_index ///////////////////////////////////////////////////
/**
   Writes a value of type Index to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the data to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_index(
        const String&   filename,
        const int&      fid,
        const Index&   x,
        const String&   dataname )
{
  /*
  Index  a[1];
  a[0] = x;
  binfile_write( fid,  filename, dataname, "SCALAR", "INDEX", 1, 1, 
                                                                (uint8*)a );
  */
  binfile_write( fid,  filename, dataname, "SCALAR", "INDEX", 1, 1, 
                                                               (uint8*)&x );
}



//// binfile_read_index ///////////////////////////////////////////////////
/**
   Reads a value of type Index from a binary file.

   \param    x            the read index
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_index(
              Index&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname )
{
  int     vdata_id;
  Index  nrows, ncols;
  ArrayOfIndex  a;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "SCALAR", 1, 1 );
  binfile_read1( a, vdata_id, nrows, filename, dataname );
  x = a[0];
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_numeric ///////////////////////////////////////////////////
/**
   Writes a value of type NUMERIC to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the data to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_numeric(
        const String&   filename,
        const int&      fid,
        const Numeric&  x,
        const String&   dataname )
{
  //Numeric  a[1];
  //a[0] = x;
  binfile_write( fid,  filename, dataname, "SCALAR", "NUMERIC", 1, 1, 
                                                               (uint8*)&x );
}



//// binfile_read_numeric ///////////////////////////////////////////////////
/**
   Reads a value of type NUMERIC from a binary file.

   \param    x            the read numeric
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_numeric(
              Numeric&  x,
        const String&   filename,
        const int&      fid,
        const String&   dataname )
{
  int     vdata_id;
  Index  nrows, ncols;
  Matrix  a;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "SCALAR", 1, 1 );
  binfile_read2( a, vdata_id, nrows, ncols, filename, dataname );
  x = a(0,0);
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_vector ///////////////////////////////////////////////////
/**
   Writes a vector to a binary file.

   \param    fname        file name
   \param    fid          file identifier
   \param    x            the vector to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_vector(
        const String&   filename,
        const int&      fid,
        const Vector&   x,
        const String&   dataname )
{
  const Index  n = x.nelem();
  
  Numeric *a = new Numeric[n];
  for ( Index i=0; i<n; i++ )
    {
      a[i] = x[i];
      //      cout << "a[i] = " << a[i] << "\n";
    }
  
  binfile_write( fid,  filename, dataname, "VECTOR", "NUMERIC", n, 1, 
                                                               (uint8*)a );

  delete a;
}



//// binfile_read_vector ///////////////////////////////////////////////////
/**
   Reads a vector from a binary file.

   \param    x            the read vector
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_vector(
              Vector&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname )
{
  int     vdata_id;
  Index  nrows, ncols;
  Matrix  a;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "VECTOR", 0, 1 );
  binfile_read2( a, vdata_id, nrows, ncols, filename, dataname );
  x.resize(nrows);
  for ( Index i=0; i<nrows; i++ )
    x[i] = a(i,0);
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_matrix ///////////////////////////////////////////////////
/**
   Writes a matrix to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the matrix to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01

   Adapted to 0-based indexing.
   \date   2001-01-08
   \author Stefan Buehler
*/
void binfile_write_matrix(
        const String&   filename,
        const int&      fid,
        const Matrix&   x,
        const String&   dataname )
{
  const Index  nrows = x.nrows();
  const Index  ncols = x.ncols();

  Numeric *a = new Numeric[nrows*ncols];
  for ( Index r=0; r<nrows; r++ )
    for ( Index c=0; c<ncols; c++ )
      a[r*ncols+c] = x(r,c);

  binfile_write( fid,  filename, dataname, "MATRIX", "NUMERIC", nrows, ncols, 
		 (uint8*)a );
  delete a;
}



//// binfile_read_matrix ///////////////////////////////////////////////////
/**
   Reads a matrix from a binary file.

   \param    x            the read matrix
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_matrix(
              Matrix&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname )
{
  int     vdata_id;
  Index  nrows, ncols;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "MATRIX", 0, 0 );
  binfile_read2( x, vdata_id, nrows, ncols, filename, dataname );
  binfile_read_end( vdata_id, filename, dataname );
}
       


//// binfile_write_indexarray //////////////////////////////////////////////
/**
   Writes an ArrayOfIndex to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_indexarray(
        const String&         filename,
        const int&            fid,
        const ArrayOfIndex&   x,
        const String&         dataname )
{
  const Index  n = x.nelem();

  /*
  Index a[n];
  for ( Index i=0; i<n; i++ )
    a[i] = x[i];
  binfile_write( fid,  filename, dataname, "ARRAY", "INDEX", n, 1, 
                                                                (uint8*)a );
  */
  binfile_write( fid,  filename, dataname, "ARRAY", "INDEX", n, 1, 
                                                               (uint8*)&x[0] );
}



//// binfile_read_indexarray ///////////////////////////////////////////////
/**
   Reads a ArrayOfIndex from a binary file.

   \param    x            the read index array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_indexarray(
              ArrayOfIndex&   x,
        const String&         filename,
        const int&            fid,
        const String&         dataname )
{
  int     vdata_id;
  Index  nrows, ncols;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "ARRAY", 0, 1 );
  binfile_read1( x, vdata_id, nrows, filename, dataname );
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_vectorarray //////////////////////////////////////////////
/**
   Writes an ArrayOfVector to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_vectorarray(
        const String&          filename,
        const int&             fid,
        const ArrayOfVector&   x,
        const String&          dataname )
{
  const Index  n = x.nelem();

  // Write number of vectors
  binfile_write_index( filename, fid, n, "N_"+dataname );  

  // Write each matrix seperately
  for (Index i=0; i<n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_write_vector( filename, fid, x[i], os.str() );    
  }
}



//// binfile_read_vectorarray ///////////////////////////////////////////////
/**
   Reads a ArrayOfVector from a binary file.

   \param    x            the read vector array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_vectorarray(
              ArrayOfVector&   x,
        const String&          filename,
        const int&             fid,
        const String&          dataname )
{
  // Read the number of vectors
  Index  n;
  binfile_read_index( n, filename, fid, "N_"+dataname );  

  // Reallocate x and read vectors
  x.resize(n);
  for (Index i=0; i<n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_read_vector( x[i], filename, fid, os.str() );    
  }
}



//// binfile_write_matrixarray //////////////////////////////////////////////
/**
   Writes an ArrayOfMatrix to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_matrixarray(
        const String&          filename,
        const int&             fid,
        const ArrayOfMatrix&   x,
        const String&          dataname )
{
  const Index  n = x.nelem();

  // Write number of matrices
  binfile_write_index( filename, fid, n, "N_"+dataname );  

  // Write each matrix seperately
  for ( Index i=0; i<n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_write_matrix( filename, fid, x[i], os.str() );    
  }
}



//// binfile_read_matrixarray ///////////////////////////////////////////////
/**
   Reads a ArrayOfMatrix from a binary file.

   \param    x            the read matrix array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_matrixarray(
              ArrayOfMatrix&   x,
        const String&          filename,
        const int&             fid,
        const String&          dataname )
{
  // Read the number of matrices
  Index  n;
  binfile_read_index( n, filename, fid, "N_"+dataname );  

  // Reallocate x and read matrices
  x.resize(n);
  for (Index i=0; i<n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_read_matrix( x[i], filename, fid, os.str() );    
  }
}



//// binfile_write_String ///////////////////////////////////////////////////
/**
   Writes a String to a binary file.

   \param    fname        file name
   \param    fid          file identifier
   \param    s            the String to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_String(
        const String&   filename,
        const int&      fid,
        const String&   s,
        const String&   dataname )
{
  const Index  n = s.length();

  binfile_write( fid,  filename, dataname, "STRING", "CHAR", n, 1, 
                                                           (uint8*)s.c_str() );
}



//// binfile_read_String ////////////////////////////////////////////////////
/**
   Reads a String from a binary file.

   \param    x            the read String
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_String(
              String&   x,
        const String&   filename,
        const int&      fid,
        const String&   dataname )
{
  int     vdata_id;
  Index  nrows, ncols;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "STRING", 0, 1 );
  binfile_read3( x, vdata_id, nrows, filename, dataname );
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_Stringarray //////////////////////////////////////////////
/**
   Writes an ArrayOfSTRING to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_Stringarray(
        const String&          filename,
        const int&             fid,
        const ArrayOfString&   x,
        const String&          dataname )
{
  const Index  n = x.nelem();

  // Write number of matrices
  binfile_write_index( filename, fid, n, "N_"+dataname );  

  // Write each String seperately
  for ( Index i=0; i<n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_write_String( filename, fid, x[i], os.str() );    
  }
}



//// binfile_read_Stringarray ///////////////////////////////////////////////
/**
   Reads a ArrayOfSTRING from a binary file.

   \param    x            the read String array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_Stringarray(
              ArrayOfString&   x,
        const String&          filename,
        const int&             fid,
        const String&          dataname )
{
  // Read the number of matrices
  Index  n;
  binfile_read_index( n, filename, fid, "N_"+dataname );  

  // Reallocate x and read matrices
  x.resize(n);
  for (Index i=0; i<n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_read_String( x[i], filename, fid, os.str() );    
  }
}
