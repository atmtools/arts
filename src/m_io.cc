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
   \file   m_io.cc

   This file contains IO workspace methods. 

   The functions are of two types:
     1. Functions with overall influence on ARTS
     2. Reading and writing to/from files.
     3. Creation by workspace method keywords and generic input

   \author Patrick Eriksson
   \date 2000-11-01 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "atm_funcs.h"          
#include "file.h"
#include "math_funcs.h"
#include "messages.h"
#include "md.h"
#include <hdf.h>



//**************************************************************************
//
//  1. Overall ARTS functions
//
//**************************************************************************

void Exit()
{
  out1 << "  Forced exit.\n";
  exit(1);
}



//**************************************************************************
//
//  2. File functions
//
//**************************************************************************

////////////////////////////////////////////////////////////////////////////
//   File workspace methods (sorted after workspace variable)
////////////////////////////////////////////////////////////////////////////

// The ASCII functions are created by Stefan Buehler around 00????
// The binary functions are created by Patrick Eriksson around 001102

//=== INDEX ============================================================

// This function shall be modified to handle INDEX
void IndexWriteAscii(
        const int&      v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ARRAYofMATRIX and write to file
  ARRAYofMATRIX am(1);
  resize( am[0], 1, 1 );
  setto( am[0], (Numeric) v);
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle INDEX
void IndexReadAscii(
	      int&      v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ARRAYofMATRIX and move to Numeric
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.size()!=1) || (am[0].nrows()!=1) || (am[0].ncols()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single value.";
    throw runtime_error(os.str());
  }

  Numeric a = am[0][0][0];

  if ( (a-floor(a)) != 0 )
    throw runtime_error("The value in the file is not an integer.");
  if ( a < 0 )
    throw runtime_error("The value in the file is negative.");

  v = (int) a;
}


// This function shall be modified to handle INDEX
void IndexWriteBinary(
	const int&      v,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_index( filename, fid, v, "INDEX" );
  binfile_close( fid, filename );
}


// This function shall be modified to handle INDEX
void IndexReadBinary(
	      int&      v,
        const string&   var_name,
	const string&   f )
{
  size_t vtemp;  //To be removed
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_index( vtemp, filename, fid, "INDEX" );
  v = (int) vtemp;
  binfile_close( fid, filename );
}



//=== NUMERIC ==========================================================

void NumericWriteAscii(
        const Numeric&  v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ARRAYofMATRIX and write to file
  ARRAYofMATRIX am(1);
  resize( am[0], 1, 1 );
  setto( am[0], v);
  write_array_of_matrix_to_file(filename,am);
}


void NumericReadAscii(
	      Numeric&  v,
        const string&   v_name,
	const string&   f )
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ARRAYofMATRIX and move to Numeric
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.size()!=1) || (am[0].nrows()!=1) || (am[0].ncols()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single numeric value";
    throw runtime_error(os.str());
  }

  v = am[0][0][0];
}



void NumericWriteBinary(
        const Numeric&  v,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_numeric( filename, fid, v, "NUMERIC" );
  binfile_close( fid, filename );
}



void NumericReadBinary(
	      Numeric&  v,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_numeric( v, filename, fid, "NUMERIC" );
  binfile_close( fid, filename );
}



//=== VECTOR ==========================================================

void VectorWriteAscii(// WS Output:
		       const VECTOR& v,
		       // WS Variable Names:
		       const string& v_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Convert the vector to a matrix:
  //  to_matrix(m,v);
  MATRIX m(v.size(),1);
  copy( v, columns(m)[0] );

  // Convert the matrix to an array of matrix:
  ARRAYofMATRIX am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void VectorReadAscii(// WS Generic Output:
			VECTOR& v,
			// WS Generic Output Names:
			const string& v_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read an array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to a matrix.
  if ( 1 != am.size() )
   throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                       "but the dimension of the array is not 1.");
  MATRIX m(am[0]);

  // Convert the matrix to a vector:
  to_vector(v,m);
}



void VectorWriteBinary(
        const VECTOR&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vector( filename, fid, v, "VECTOR" );
  binfile_close( fid, filename );
}



void VectorReadBinary(
	      VECTOR&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vector( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
}



//=== MATRIX ==========================================================

void MatrixWriteAscii(// WS Generic Input:
		       const MATRIX& m,
		       // WS Generic Input Names:
		       const string& m_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, m_name );

  // Convert the matrix to an array of matrix:
  ARRAYofMATRIX am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void MatrixReadAscii(// WS Generic Output:
			MATRIX& m,
			// WS Generic Output Names:
			const string& m_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, m_name );

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  //  cout << "am.size(): " << am.size() << "\n";

  // Convert the array of matrix to a matrix.
  if ( 1 != am.size() )
   throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                       "but the dimension of the array is not 1.");

  m = am[0];
}



void MatrixWriteBinary(
        const MATRIX&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrix( filename, fid, v, "MATRIX" );
  binfile_close( fid, filename );
}



void MatrixReadBinary(
	      MATRIX&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrix( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== ARRAYofINDEX =====================================================

// This function shall be modified to handle ARRAYofINDEX
void ArrayOfIndexWriteAscii(
        const ARRAYofsizet&   v,
        const string&         v_name,
	const string&         f )
{
  string filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ARRAYofMATRIX and write to file
  const size_t  n = v.size();
  ARRAYofMATRIX am(1);
  resize(am[0],n,1);
  for ( size_t i=0; i<n; i++ )
    am[0][i][0] = (Numeric) v[i];
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle INDEX
void ArrayOfIndexReadAscii(
	      ARRAYofsizet&   v,
        const string&         v_name,
	const string&         f )
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ARRAYofMATRIX 
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.size()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains more than one vector.";
    throw runtime_error(os.str());
  }

  // Move values to v using a vector
  VECTOR x;
  to_vector( x, am[0] );
  const size_t  n = x.size();
  resize(v,n);
  for ( size_t i=0; i<n; i++ )
  {  
    if ( (x[i]-floor(x[i])) != 0 )
      throw runtime_error("A value in the file is not an integer.");
    if ( x[i] < 0 )
      throw runtime_error("A value in the file is negative.");
    v[i] = (size_t) x[i];
  }
}



void ArrayOfIndexWriteBinary(
        const ARRAYofsizet&  v,
        const string&        var_name,
	const string&        f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_indexarray( filename, fid, v, "INDEXARRAY" );
  binfile_close( fid, filename );
}



void ArrayOfIndexReadBinary(
	      ARRAYofsizet&  v,
        const string&        var_name,
	const string&        f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_indexarray( v, filename, fid, "INDEXARRAY" );
  binfile_close( fid, filename );
}



//=== ARRAYofVECTOR ====================================================

void ArrayOfVectorWriteAscii(// WS Output:
			      const ARRAYofVECTOR& av,
			      // WS Variable Names:
			      const string& av_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, av_name );

  // Convert the array of vector to an array of matrix:
  ARRAYofMATRIX am(av.size());
  for (size_t i=0; i<av.size(); ++i)
    {
      //      to_matrix(am[i],av[i]);
      resize( am[i], av[i].size(), 1 );
      copy( av[i], columns(am[i])[0] ); 
    }

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void ArrayOfVectorReadAscii(// WS Generic Output:
			       ARRAYofVECTOR& av,
			       // WS Generic Output Names:
			       const string& av_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, av_name );

  // Read an array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to an array of vector.
  resize(av,am.size());
  for (size_t i=0; i<am.size(); ++i)
    {
      to_vector(av[i],am[i]);
    }
}



void ArrayOfVectorWriteBinary(
        const ARRAYofVECTOR&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, v, "VECTOR" );
  binfile_close( fid, filename );
}



void ArrayOfVectorReadBinary(
	      ARRAYofVECTOR&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vectorarray( v, filename, fid, "VECTOR" );
  binfile_close( fid, filename );
}



//=== ARRAYofMATRIX ====================================================

void ArrayOfMatrixWriteAscii(// WS Generic Input:
			      const ARRAYofMATRIX& am,
			      // WS Generic Input Names:
			      const string& am_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, am_name );

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void ArrayOfMatrixReadAscii(// WS Generic Output:
			       ARRAYofMATRIX& am,
			       // WS Generic Output Names:
			       const string& am_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, am_name );

  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);
}



void ArrayOfMatrixWriteBinary(
        const ARRAYofMATRIX&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrixarray( filename, fid, v, "MATRIX" );
  binfile_close( fid, filename );
}



void ArrayOfMatrixReadBinary(
	      ARRAYofMATRIX&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrixarray( v, filename, fid, "MATRIX" );
  binfile_close( fid, filename );
}



//=== STRING ===============================================================

void StringWriteAscii( // WS Generic Input:
		       const string& s,
		       // WS Generic Input Names:
		       const string& s_name,
		       // Control Parameters:
		       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, s_name );

  // Convert the string to an array of string:
  ARRAYofstring as(1);
  as[0] = s;

  // Write the array of matrix to the file.
  write_array_of_string_to_file(filename,as);
}



void StringReadAscii(   // WS Generic Output:
			string& s,
			// WS Generic Output Names:
			const string& s_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, s_name );

  // Read the array of matrix from the file:
  ARRAYofstring as;
  read_array_of_string_from_file(as,filename);

  // Convert the array of string to a string.
  if ( 1 != as.size() )
   throw runtime_error("You tried to convert an array of string to a string,\n"
                       "but the dimension of the array is not 1.");

  s = as[0];
}



void StringWriteBinary(
        const string&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_string( filename, fid, v, "STRING" );
  binfile_close( fid, filename );
}



void StringReadBinary(
	      string&  v,
        const string&  var_name,
	const string&  f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_string( v, filename, fid, "STRING" );
  binfile_close( fid, filename );
}



//=== ARRAYofstring ====================================================

void ArrayOfStringWriteAscii( // WS Generic Input:
			      const ARRAYofstring& as,
			      // WS Generic Input Names:
			      const string& as_name,
			      // Control Parameters:
			      const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, as_name );

  // Write the array of matrix to the file.
  write_array_of_string_to_file(filename,as);
}



void ArrayOfStringReadAscii(   // WS Generic Output:
			       ARRAYofstring& as,
			       // WS Generic Output Names:
			       const string& as_name,
			       // Control Parameters:
			       const string& f)
{
  string filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, as_name );

  // Read the array of matrix from the file:
  read_array_of_string_from_file(as,filename);
}



void ArrayOfStringWriteBinary(
        const ARRAYofstring&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_stringarray( filename, fid, v, "STRING" );
  binfile_close( fid, filename );
}



void ArrayOfStringReadBinary(
	      ARRAYofstring&  v,
        const string&         var_name,
	const string&         f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_stringarray( v, filename, fid, "STRING" );
  binfile_close( fid, filename );
}



//=== MAYBESPARSE ====================================================

void HmatrixReadAscii(// WS Generic Output:
			Hmatrix& h,
			// WS Generic Output Names:
			const string& h_name,
			// Control Parameters:
			const string& f)
{
  string filename = f;
 
  // Create default filename if empty  
  filename_ascii( filename, h_name );

  // Read the array of matrix from the file:
  ARRAYofMATRIX am;
  read_array_of_matrix_from_file(am,filename);

  h.issparse = 0;
  h.full     = am[0];
  // Remember to clear the other field
}



//=== LOS ==================================================================

void LosWriteBinary(
        const LOS&      los,
        const string&   var_name,
	const string&   f )
{
  int    fid;
  string filename = f;
  filename_bin( filename, var_name );

  // Handle ground by adding 1 and handle as INDEX
  const size_t   n = los.ground.size();
  ARRAYofsizet   ground(n);
  for ( size_t i=0; i<n; i++ )
    ground[i] = size_t( los.ground[i] + 1 );

  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, los.p, "LOS.P" );
  binfile_write_vector( filename, fid, los.l_step, "LOS.L_STEP" );
  binfile_write_indexarray( filename, fid, ground, "LOS.GROUND" );
  binfile_write_indexarray( filename, fid, los.start, "LOS.START" );
  binfile_write_indexarray( filename, fid, los.stop, "LOS.STOP" );
  binfile_close( fid, filename );
}



void LosReadBinary(
	      LOS&      los,
        const string&   var_name,
	const string&   f )
{
  int   fid;
  string filename = f;
  filename_bin( filename, var_name );

  // Ground is stored as 1-based
  ARRAYofsizet   ground;

  binfile_open_in( fid, filename );
  binfile_read_vectorarray( los.p, filename, fid, "LOS.P" );
  binfile_read_vector( los.l_step, filename, fid, "LOS.L_STEP" );
  binfile_read_indexarray( ground, filename, fid, "LOS.GROUND" );
  binfile_read_indexarray( los.start, filename, fid, "LOS.START" );
  binfile_read_indexarray( los.stop, filename, fid, "LOS.STOP" );
  binfile_close( fid, filename );

  const size_t   n = ground.size();
  resize(los.ground,n);
  for ( size_t i=0; i<n; i++ )
    los.ground[i] = int ( ground[i] - 1 );
}





//**************************************************************************
//
//   3. Creation by workspace method keywords
//
//**************************************************************************

//=== INDEX ============================================================

void IntSet(// WS Generic Output:
            int& x,
            // WS Generic Output Names:
            const string& x_name,
            // Control Parameters:
            const int& value)
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



//=== NUMERIC ==========================================================

void NumericSet(// WS Generic Output:
                Numeric& x,
                // WS Generic Output Names:
                const string& x_name,
                // Control Parameters:
                const Numeric& value)
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



//=== VECTOR ==========================================================

void VectorSet(           VECTOR&  x, 
                    const string&  x_name,
                    const int&     n,
                    const Numeric& value )
{
  resize(x,n);
  setto(x,value);
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
}



void VectorSetLengthFromVector(
              VECTOR&  x, 
        const string&  x_name,
        const VECTOR&  z,
        const string&  z_name,
        const Numeric& value )
{
  const size_t  n = z.size();
  resize(x,n);
  setto(x,value);
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
}



void VectorLinSpace(      VECTOR&  x, 
                    const string&  x_name,
                    const Numeric& start,
                    const Numeric& stop,
                    const Numeric& step )
{
  x = linspace(start,stop,step);
  out2 << "  Creating " << x_name << " as linearly spaced vector\n";
  out3 << "         length: " << x.size() << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.size() > 1 )
  {
    out3 << "      step size: " << x[1]-x[0] << "\n";
    out3 << "     last value: " << x[x.size()-1] << "\n";
  }
}



void VectorNLinSpace(     VECTOR&  x, 
                    const string&  x_name,
                    const Numeric& start,
                    const Numeric& stop,
                    const int& n )
{
  x = nlinspace(start,stop,n);
  out2 << "  Creating " << x_name << " as linearly spaced vector\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.size() > 1 )
  {
    out3 << "      step size: " << x[1]-x[0] << "\n";
    out3 << "     last value: " << x[x.size()-1] << "\n";
  }
}



void VectorNLogSpace(     VECTOR&  x, 
                    const string&  x_name,
                    const Numeric& start,
                    const Numeric& stop,
                    const int&     n )
{
  x = nlogspace(start,stop,n);
  out2 << "  Creating " << x_name << " as logarithmically spaced vector\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.size() > 1 )
    out3 << "     last value: " << x[x.size()-1] << "\n";
}



void VectorCopy(
                      VECTOR&   y2,
                const string&   name_y2,
                const VECTOR&   y1,
                const string&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  resize( y2, y1.size() );
  copy( y1, y2 );
}



void VectorCopyFromArrayOfVector(
                      VECTOR&          v,
                const string&          v_name,
                const ARRAYofVECTOR&   va,
                const string&          va_name,
                const int&             i )
{
  if ( i < 0 )
    throw runtime_error("The index must be >= 0.");
  if ( size_t(i) >= va.size() )
  {
    ostringstream os;
    os << "The vector array has only " << va.size() << "elements.";
    throw runtime_error(os.str());
  }

  out2 << "  Copies " << v_name << " from vector " << i << " in " << va_name
                                                                   << "\n";
  resize(v,va[i].size());
  copy(va[i],v);
}



void VectorPlanck(
                    VECTOR&   y,
              const string&   y_name,
              const VECTOR&   f,
              const string&   f_name,
              const Numeric&  t )
{
  if ( t > 0 )
  {
    resize( y, f.size() );
    planck( y, f, t );
    out2<<"  Setting " << y_name << " to blackbody radiation for "<<t<<" K.\n";
  }
  else
    throw runtime_error("The temperature must be > 0.");
}



void VectorCalcLog10(
                    VECTOR&   out,
              const string&   out_name,
              const VECTOR&   in,
              const string&   in_name )
{
  out2<<"  " << out_name << " = log10( " << in_name << " )\n";

  out = transf( in, log10 );
}



void VectorRandUniform(
                    VECTOR&   y,
              const string&   y_name,
              const Numeric&  x_low,
              const Numeric&  x_high,
              const int&      n )
{
  out2<<"  Filling " << y_name << " with uniform random data.\n";
  resize( y, n );
  rand_uniform( y, x_low, x_high );
}



void VectorRandGaussian(
                    VECTOR&   y,
              const string&   y_name,
              const Numeric&  stddev,
              const int&      n )
{
  out2<<"  Filling " << y_name << " with Gaussian random data.\n";
  resize( y, n );
  rand_gaussian( y, stddev );
}



//=== MATRIX ==========================================================

void MatrixCopy(
              MATRIX&   y2,
        const string&   name_y2,
        const MATRIX&   y1,
        const string&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  resize( y2, y1.nrows(), y1.ncols() );
  copy( y1, y2 );
}



void MatrixFillWithVector(
              MATRIX&   m,
        const string&   name_m,
        const VECTOR&   y,
        const string&   name_y,
        const int&      n )
{
  out2 << "  Creates" << name_m << " by copying " << name_y << n << "times.\n";
  const size_t nrows = y.size();
  resize( m, nrows, n );
  for ( INDEX row=0; row<nrows; row++ ) 
    mtl::set(m[row],y[row]);
}



//=== STRING ===============================================================

void StringSet(           string&  s, 
                    const string&  s_name,
                    const string&  s2 )
{
  s = s2;
  out3 << "  Setting " << s_name << " to " << s2 << "\n"; 
}



//=== ARRAYofSTRING ========================================================

void ArrayOfStringSet(    
              ARRAYofstring&  sa, 
        const string&         sa_name,
        const ARRAYofstring&  sa2 )
{
  resize(sa,sa2.size());
  copy(sa2,sa);
  out3 << "  Setting " << sa_name << "\n"; 
}



