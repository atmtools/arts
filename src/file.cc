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
#include <hdf.h>



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
          getline(is,linebuffer);

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
//   STRING IO routines for ASCII files
////////////////////////////////////////////////////////////////////////////

//// write_array_of_string_to_stream ///////////////////////////////////////
/**
   A helper function that writes an array of string to a stream. 

   \retval   os           the stream to write to
   \param    as           the string array to write

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void write_array_of_string_to_stream(
              ostream&         os,
        const ARRAYofstring&   as )
{
  extern const string full_name;

  os << "# Generated by "
     << full_name << ", "
     << __DATE__ << ", "
     << __TIME__ << "\n";

  // Number of array elements:
  const size_t n = as.dim();
  os << n << '\n';

  for (size_t i=1; i<=n; ++i)
    os << as(i) << '\n';
}



//// write_array_of_string_to_file /////////////////////////////////////////
/**
   A help function that writes an array of string to a file.

   \retval filename    The name of the file.  
   \param  as          The array of string to write.  

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void write_array_of_string_to_file(
        const string&          filename,
        const ARRAYofstring&   as )
{
  ofstream of;

  out2 << "  Writing " << filename << '\n';
  open_output_file(of, filename);

  // Write the array of string to the stream:
  write_array_of_string_to_stream(of,as);
}



//// read_array_of_string_from_stream ///////////////////////////////////////
/**
   A help function to read an array of string from a stream.

   \retval as   The array of string to read
   \param  is   The input stream

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void read_array_of_string_from_stream(
        ARRAYofstring&   as,
        istream&         is )
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
          getline(is,linebuffer);

      else
        {
          is.unget();
          comments = false;
        }
    }

  is >> as;

  if ( is.fail() || is.bad() )
    throw runtime_error("Stream gave fail or bad.");

  is >> ws;
  if ( !is.eof() )
    throw runtime_error("Input finished, but end of stream not reached.");

  // Some output to the lowest priority stream:
  out3 << "  Dimension: "
       << as.dim() << ", ";
}



//// read_array_of_matrix_from_file /////////////////////////////////////////
/**
   A help function to read an array of string from a file.

   \retval as        The array of string to read
   \param  filename  The name of the file to read

   \author Patrick Eriksson              
   \date   2000-11-04
*/
void read_array_of_string_from_file(
           ARRAYofstring&   as,
     const string&          filename )
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
      read_array_of_string_from_stream(as,ifs);
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

   The following data types are checked: INDEX, float and double.

   @exception runtime_error  Some data type has an unexpected length.

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void check_data_types()
{
  if ( sizeof(size_t) != 4 )
    throw runtime_error("An INDEX is expected to be 4 bytes.");
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
        const string&   filename )
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
        const string&   filename )
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
        const string&   filename )
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
        const string&   filename, 
        const string&   dataname,
        const int&      vdata_id,
        const size_t&   nrows,
        const size_t&   ncols )
{
  size_t v[2];
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
   \param    storagetype  SCALAR, VECTOR, MATRIX etc.
   \param    nrows0       expected number of data rows
   \param    ncols0       expected number of data columns

   @exception runtime_error  Cannot find or access the data

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_read_init(
              int&      vdata_id,
              size_t&   nrows,
              size_t&   ncols,
        const int&      fid,
        const string&   filename,
        const string&   dataname,
        const string&   storagetype,
        const size_t&   nrows0,
        const size_t&   ncols0 )
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
  size_t  v[2];
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
        const string&   filename,
        const string&   dataname )
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

   \retval  type_in_file   data type string (e.g. "DOUBLE")
   \param   vdata_id       identifier for the Vdata

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_get_datatype(
              string&   type_in_file,
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

   The data can be a scalar, a vector or a matrix holding INDEX or NUMERIC.   

   \param    fid          file identifier
   \param    filename     file name
   \param    dataname     name on data set
   \param    storagetype  type of data container (e.g. VECTOR or MATRIX)
   \param    atomictype   basic data type (INDEX or NUMERIC)
   \param    nrows        number of data rows
   \param    ncols        number of data columns
   \param    dpointer     pointer to the start of the data to store

   @exception runtime_error  Unknown atomic data type.

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write(
        const int&      fid,
        const string&   filename,
        const string&   dataname,
        const string&   storagetype,
        const string&   atomictype,
        const size_t&   nrows,
        const size_t&   ncols,
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
    const size_t   nout = nrows*ncols;
    size_t ndone = VSwrite( vdata_id, dpointer, nout, FULL_INTERLACE );
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
              ARRAYofsizet&   x,
        const int&            vdata_id,
        const size_t&         n,
        const string&         filename,
        const string&         dataname )
{
  // Check that data types have expected length
  check_data_types();

  out3 << "    Reading: " << dataname << "\n";

  // Get the data type of the data to read
  string type_in_file;
  binfile_get_datatype( type_in_file, vdata_id );

  // Reallocate x
  x.newsize(n);

  if ( n > 0 )
  {
    // Do reading and copy data
    if ( type_in_file == "UINT" )
    {
      size_t  a[n];
      VSread( vdata_id, (uint8*)a, n, FULL_INTERLACE );
      for ( size_t i=0; i<n; i++ )
	x[i] = a[i];
    }
  
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
*/
void binfile_read2(
              MATRIX&   x,
        const int&      vdata_id,
        const size_t&   nrows,
        const size_t&   ncols,
        const string&   filename,
        const string&   dataname )
{
  // Check that data types have expected length
  check_data_types();

  out3 << "    Reading: " << dataname << "\n";

  // Get the data type of the data to read
  string type_in_file;
  binfile_get_datatype( type_in_file, vdata_id );

  // Reallocate x
  x.newsize(nrows,ncols);

  if ( (nrows > 0) && (ncols > 0) )
  {
    // Do reading and copy data
    if ( type_in_file == "FLOAT" )
    {
      float  a[nrows*ncols];
      size_t i,j,j0;
      VSread( vdata_id, (uint8*)a, nrows*ncols, FULL_INTERLACE );
      for ( i=1; i<=nrows; i++ )
      {
	j0 = (i-1)*ncols-1;
	for ( j=1; j<=ncols; j++ )
	  x(i,j) = a[j0+j];
      }
    }
  
    else if ( type_in_file == "DOUBLE" )
    {
      double a[nrows*ncols];
      size_t i,j,j0;
      VSread( vdata_id, (uint8*)a, nrows*ncols, FULL_INTERLACE );
      for ( i=1; i<=nrows; i++ )
      {
	j0 = (i-1)*ncols-1;
	for ( j=1; j<=ncols; j++ )
	   x(i,j) = a[j0+j];
      }
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

   The data is read a string.

   The reading is splitted to several files as type conversions can be needed.

   \retval   x            the string to read
   \param    vdata_id     data identifier
   \param    nrows        number of values
   \param    filename     file name
   \param    dataname     name on data set

   @exception runtime_error  Unknown data type in file.

   \author Patrick Eriksson              
   \date   2000-11-07
*/
void binfile_read3(
              string&   x,
        const int&      vdata_id,
        const size_t&   n,
        const string&   filename,
        const string&   dataname )
{
  // Check that data types have expected length
  check_data_types();

  out3 << "    Reading: " << dataname << "\n";

  // Get the data type of the data to read
  string type_in_file;
  binfile_get_datatype( type_in_file, vdata_id );

  // Reallocate x
  x.resize(n);

  if ( n > 0 )
  {
    // Do reading and copy data
    if ( type_in_file == "CHAR" )
    {
      char  a[n];
      VSread( vdata_id, (uint8*)a, n, FULL_INTERLACE );
      for ( size_t i=0; i<n; i++ )
	x[i] = a[i];
    }
  
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
   Writes a value of type INDEX to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the data to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_index(
        const string&   filename,
        const int&      fid,
        const size_t&   x,
        const string&   dataname )
{
  size_t  a[1];
  a[0] = x;
  binfile_write( fid,  filename, dataname, "SCALAR", "INDEX", 1, 1, 
                                                                (uint8*)a );
}



//// binfile_read_index ///////////////////////////////////////////////////
/**
   Reads a value of type INDEX from a binary file.

   \param    x            the read index
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_index(
              size_t&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname )
{
  int     vdata_id;
  size_t  nrows, ncols;
  ARRAYofsizet  a;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "SCALAR", 1, 1 );
  binfile_read1( a, vdata_id, nrows, filename, dataname );
  x = a(1);
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
        const string&   filename,
        const int&      fid,
        const Numeric&  x,
        const string&   dataname )
{
  Numeric  a[1];
  a[0] = x;
  binfile_write( fid,  filename, dataname, "SCALAR", "NUMERIC", 1, 1, 
                                                                (uint8*)a );
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
        const string&   filename,
        const int&      fid,
        const string&   dataname )
{
  int     vdata_id;
  size_t  nrows, ncols;
  MATRIX  a;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "SCALAR", 1, 1 );
  binfile_read2( a, vdata_id, nrows, ncols, filename, dataname );
  x = a(1,1);
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
        const string&   filename,
        const int&      fid,
        const VECTOR&   x,
        const string&   dataname )
{
  const size_t  n = x.dim();

  Numeric a[n];
  for ( size_t i=0; i<n; i++ )
    a[i] = x(i+1);
  binfile_write( fid,  filename, dataname, "VECTOR", "NUMERIC", n, 1, 
                                                                (uint8*)a );
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
              VECTOR&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname )
{
  int     vdata_id;
  size_t  nrows, ncols;
  MATRIX  a;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "VECTOR", 0, 1 );
  binfile_read2( a, vdata_id, nrows, ncols, filename, dataname );
  x.newsize(nrows);
  for ( size_t i=1; i<=nrows; i++ )
    x(i) = a(i,1);
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
*/
void binfile_write_matrix(
        const string&   filename,
        const int&      fid,
        const MATRIX&   x,
        const string&   dataname )
{
  const size_t  nrows = x.dim(1);
  const size_t  ncols = x.dim(2);
        size_t  i, j;

  Numeric a[nrows*ncols];
  size_t j0;
  for ( i=1; i<=nrows; i++ )
  {
    j0 = (i-1)*ncols-1;
    for ( j=1; j<=ncols; j++ )
      a[j0+j] = x(i,j);
  }
  binfile_write( fid,  filename, dataname, "MATRIX", "NUMERIC", nrows, ncols, 
                                                                (uint8*)a );
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
              MATRIX&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname )
{
  int     vdata_id;
  size_t  nrows, ncols;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "MATRIX", 0, 0 );
  binfile_read2( x, vdata_id, nrows, ncols, filename, dataname );
  binfile_read_end( vdata_id, filename, dataname );
}
       


//// binfile_write_indexarray //////////////////////////////////////////////
/**
   Writes an ARRAYofINDEX to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_indexarray(
        const string&         filename,
        const int&            fid,
        const ARRAYofsizet&   x,
        const string&         dataname )
{
  const size_t  n = x.dim();

  size_t a[n];
  for ( size_t i=0; i<n; i++ )
    a[i] = x(i+1);
  binfile_write( fid,  filename, dataname, "ARRAY", "INDEX", n, 1, 
                                                                (uint8*)a );
}



//// binfile_read_indexarray ///////////////////////////////////////////////
/**
   Reads a ARRAYofINDEX from a binary file.

   \param    x            the read index array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_indexarray(
              ARRAYofsizet&   x,
        const string&         filename,
        const int&            fid,
        const string&         dataname )
{
  int     vdata_id;
  size_t  nrows, ncols;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "ARRAY", 0, 1 );
  binfile_read1( x, vdata_id, nrows, filename, dataname );
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_vectorarray //////////////////////////////////////////////
/**
   Writes an ARRAYofVECTOR to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_vectorarray(
        const string&          filename,
        const int&             fid,
        const ARRAYofVECTOR&   x,
        const string&          dataname )
{
  const size_t  n = x.dim();

  // Write number of vectors
  binfile_write_index( filename, fid, n, "N_"+dataname );  

  // Write each matrix seperately
  for (size_t i=1; i<=n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_write_vector( filename, fid, x(i), os.str() );    
  }
}



//// binfile_read_vectorarray ///////////////////////////////////////////////
/**
   Reads a ARRAYofVECTOR from a binary file.

   \param    x            the read vector array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_vectorarray(
              ARRAYofVECTOR&   x,
        const string&          filename,
        const int&             fid,
        const string&          dataname )
{
  // Read the number of vectors
  size_t  n;
  binfile_read_index( n, filename, fid, "N_"+dataname );  

  // Reallocate x and read vectors
  x.newsize(n);
  for (size_t i=1; i<=n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_read_vector( x(i), filename, fid, os.str() );    
  }
}



//// binfile_write_matrixarray //////////////////////////////////////////////
/**
   Writes an ARRAYofMATRIX to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_matrixarray(
        const string&          filename,
        const int&             fid,
        const ARRAYofMATRIX&   x,
        const string&          dataname )
{
  const size_t  n = x.dim();

  // Write number of matrices
  binfile_write_index( filename, fid, n, "N_"+dataname );  

  // Write each matrix seperately
  for ( size_t i=1; i<=n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_write_matrix( filename, fid, x(i), os.str() );    
  }
}



//// binfile_read_matrixarray ///////////////////////////////////////////////
/**
   Reads a ARRAYofMATRIX from a binary file.

   \param    x            the read matrix array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_matrixarray(
              ARRAYofMATRIX&   x,
        const string&          filename,
        const int&             fid,
        const string&          dataname )
{
  // Read the number of matrices
  size_t  n;
  binfile_read_index( n, filename, fid, "N_"+dataname );  

  // Reallocate x and read matrices
  x.newsize(n);
  for (size_t i=1; i<=n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_read_matrix( x(i), filename, fid, os.str() );    
  }
}



//// binfile_write_string ///////////////////////////////////////////////////
/**
   Writes a string to a binary file.

   \param    fname        file name
   \param    fid          file identifier
   \param    s            the string to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_string(
        const string&   filename,
        const int&      fid,
        const string&   s,
        const string&   dataname )
{
  const size_t  n = s.length();

  binfile_write( fid,  filename, dataname, "STRING", "CHAR", n, 1, 
                                                           (uint8*)s.c_str() );
}



//// binfile_read_string ////////////////////////////////////////////////////
/**
   Reads a string from a binary file.

   \param    x            the read string
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_string(
              string&   x,
        const string&   filename,
        const int&      fid,
        const string&   dataname )
{
  int     vdata_id;
  size_t  nrows, ncols;

  binfile_read_init( vdata_id, nrows, ncols, fid, filename, dataname, 
                                                           "STRING", 0, 1 );
  binfile_read3( x, vdata_id, nrows, filename, dataname );
  binfile_read_end( vdata_id, filename, dataname );
}



//// binfile_write_stringarray //////////////////////////////////////////////
/**
   Writes an ARRAYofSTRING to a binary file.

   \param    filename     file name
   \param    fid          file identifier
   \param    x            the array to store
   \param    dataname     name to give the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_write_stringarray(
        const string&          filename,
        const int&             fid,
        const ARRAYofstring&   x,
        const string&          dataname )
{
  const size_t  n = x.dim();

  // Write number of matrices
  binfile_write_index( filename, fid, n, "N_"+dataname );  

  // Write each string seperately
  for ( size_t i=1; i<=n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_write_string( filename, fid, x(i), os.str() );    
  }
}



//// binfile_read_stringarray ///////////////////////////////////////////////
/**
   Reads a ARRAYofSTRING from a binary file.

   \param    x            the read string array
   \param    filename     file name
   \param    fid          file identifier
   \param    dataname     name of the data set

   \author Patrick Eriksson              
   \date   2000-11-01
*/
void binfile_read_stringarray(
              ARRAYofstring&   x,
        const string&          filename,
        const int&             fid,
        const string&          dataname )
{
  // Read the number of matrices
  size_t  n;
  binfile_read_index( n, filename, fid, "N_"+dataname );  

  // Reallocate x and read matrices
  x.newsize(n);
  for (size_t i=1; i<=n; i++ )
  {
    ostringstream os;
    os << dataname << i;
    binfile_read_string( x(i), filename, fid, os.str() );    
  }
}
