/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>
                            Patrick Eriksson <patrick@rss.chalmers.se>
                            
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
   \file   m_hdf.cc

   This file contains all functions to read and write the binary HDF files. 

   \author Patrick Eriksson
   \date 2001-09-26 
*/

////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"

// This has to be below the include for arts.h, sinde HDF_SUPPORT is
// defined in config.h, which is loaded by arts.h.
#ifdef HDF_SUPPORT

#include <hdf.h>
#include <cmath>
#include "math_funcs.h"
#include "messages.h"
#include "auto_md.h"
#include "make_array.h"


////////////////////////////////////////////////////////////////////////////
//   Some help functions
////////////////////////////////////////////////////////////////////////////

//// filename_bin ////////////////////////////////////////////////////////////
/**
   Gives the default file name for the binary format.

   The default name is only used if the file name is empty.

   \param   filename Output:     file name
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
    extern const String out_basename;                       
    filename = out_basename+"."+varname+".ab";
  }
}



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

   \param   fid Output:          file identifier
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

   \param   fid Output:          file identifier
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

   \param   fid Output:          file identifier
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

   \param   vdata_id Output:     identifier for the Vdata
   \param   nrows Output:        number of data rows
   \param   ncols Output:        number of data columns
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

   \param   vdata_id Output:     identifier for the Vdata
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

   \param  type_in_file Output:   data type String (e.g. "DOUBLE")
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

   \param   x Output:            the index array to read
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

   \param   x Output:            the matrix to read
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

   \param   x Output:            the String to read
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
    a[i] = x[i];
  
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




////////////////////////////////////////////////////////////////////////////
//   The workspace methods (sorted after workspace variable)
////////////////////////////////////////////////////////////////////////////

//=== Index ============================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
// This function shall be modified to handle Index
void IndexWriteBinary(
        const Index&      v,
        const String&   var_name,
        const String&   f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_index( filename, fid, v, "Index" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
// This function shall be modified to handle Index
void IndexReadBinary(
              Index&      v,
        const String&   var_name,
        const String&   f )
{
  Index vtemp;  //To be removed
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_index( vtemp, filename, fid, "Index" );
  v = (int) vtemp;
  binfile_close( fid, filename );
}



//=== NUMERIC ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void NumericWriteBinary(
        const Numeric&  v,
        const String&   var_name,
        const String&   f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_numeric( filename, fid, v, "NUMERIC" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void NumericReadBinary(
              Numeric&  v,
        const String&   var_name,
        const String&   f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_numeric( v, filename, fid, "NUMERIC" );
  binfile_close( fid, filename );
}



//=== Vector ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void VectorWriteBinary(
        const Vector&  v,
        const String&  var_name,
        const String&  f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vector( filename, fid, v, "VECTOR" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void VectorReadBinary(
              Vector&  v,
        const String&  var_name,
        const String&  f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vector( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
}



//=== Matrix ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void MatrixWriteBinary(
        const Matrix&  v,
        const String&  var_name,
        const String&  f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrix( filename, fid, v, "MATRIX" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void MatrixReadBinary(
              Matrix&  v,
        const String&  var_name,
        const String&  f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrix( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== ArrayOfIndex =====================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfIndexWriteBinary(
        const ArrayOfIndex&  v,
        const String&        var_name,
        const String&        f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_indexarray( filename, fid, v, "INDEXARRAY" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfIndexReadBinary(
              ArrayOfIndex&  v,
        const String&        var_name,
        const String&        f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_indexarray( v, filename, fid, "INDEXARRAY" );
  binfile_close( fid, filename );
}



//=== ArrayOfVector ====================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfVectorWriteBinary(
        const ArrayOfVector&  v,
        const String&         var_name,
        const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, v, "VECTOR" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfVectorReadBinary(
              ArrayOfVector&  v,
        const String&         var_name,
        const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vectorarray( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
}



//=== ArrayOfMatrix ====================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfMatrixWriteBinary(
        const ArrayOfMatrix&  v,
        const String&         var_name,
        const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrixarray( filename, fid, v, "MATRIX" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfMatrixReadBinary(
              ArrayOfMatrix&  v,
        const String&         var_name,
        const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrixarray( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== STRING ===============================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void StringWriteBinary(
        const String&  v,
        const String&  var_name,
        const String&  f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_String( filename, fid, v, "STRING" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void StringReadBinary(
              String&  v,
        const String&  var_name,
        const String&  f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_String( v, filename, fid, "STRING" );
  binfile_close( fid, filename );
}



//=== ArrayOfString ====================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfStringWriteBinary(
        const ArrayOfString&  v,
        const String&         var_name,
        const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_Stringarray( filename, fid, v, "STRING" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-11-02
*/
void ArrayOfStringReadBinary(
              ArrayOfString&  v,
        const String&         var_name,
        const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_Stringarray( v, filename, fid, "STRING" );
  binfile_close( fid, filename );
}




#endif // HDF_SUPPORT
