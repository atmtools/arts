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
#include "auto_md.h"
#include <hdf.h>



//**************************************************************************
//
//  1. Overall ARTS functions
//
//**************************************************************************

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void Exit()
{
  out1 << "  Forced exit.\n";
  exit(1);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-14
*/
void RandSetSeed( )
{
  srand( (unsigned int) time( NULL ) );
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-05-15
*/
void Test( )
{
  // This function can be used to test stuff.
  Matrix z(5,3);
  Vector z0(5,1.0);
  SYMMETRIC s;

  SymmetricDiagonal( s, "xxx", 5, 2 );

  rand_data_gaussian( z, z0, s );
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

//=== Index ============================================================

// This function shall be modified to handle Index
void IndexWriteAscii(
        const int&      v,
        const String&   v_name,
	const String&   f )
{
  String filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ArrayofMatrix and write to file
  ArrayofMatrix am(1);
  resize( am[0], 1, 1 );
  setto( am[0], (Numeric) v);
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle Index
void IndexReadAscii(
	      int&      v,
        const String&   v_name,
	const String&   f )
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ArrayofMatrix and move to Numeric
  ArrayofMatrix am;
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
// This function shall be modified to handle Index
void IndexWriteBinary(
	const int&      v,
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
   \date   2000-?-?
*/
// This function shall be modified to handle Index
void IndexReadBinary(
	      int&      v,
        const String&   var_name,
	const String&   f )
{
  size_t vtemp;  //To be removed
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_index( vtemp, filename, fid, "Index" );
  v = (int) vtemp;
  binfile_close( fid, filename );
}



//=== NUMERIC ==========================================================

void NumericWriteAscii(
        const Numeric&  v,
        const String&   v_name,
	const String&   f )
{
  String filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ArrayofMatrix and write to file
  ArrayofMatrix am(1);
  resize( am[0], 1, 1 );
  setto( am[0], v);
  write_array_of_matrix_to_file(filename,am);
}


void NumericReadAscii(
	      Numeric&  v,
        const String&   v_name,
	const String&   f )
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ArrayofMatrix and move to Numeric
  ArrayofMatrix am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.size()!=1) || (am[0].nrows()!=1) || (am[0].ncols()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single numeric value";
    throw runtime_error(os.str());
  }

  v = am[0][0][0];
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
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
   \date   2000-?-?
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

void VectorWriteAscii(// WS Output:
		       const Vector& v,
		       // WS Variable Names:
		       const String& v_name,
		       // Control Parameters:
		       const String& f)
{
  String filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Convert the vector to a matrix:
  //  to_matrix(m,v);
  Matrix m(v.size(),1);
  copy( v, columns(m)[0] );

  // Convert the matrix to an array of matrix:
  ArrayofMatrix am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void VectorReadAscii(// WS Generic Output:
			Vector& v,
			// WS Generic Output Names:
			const String& v_name,
			// Control Parameters:
			const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read an array of matrix from the file:
  ArrayofMatrix am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to a matrix.
  if ( 1 != am.size() )
   throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                       "but the dimension of the array is not 1.");
  Matrix m(am[0]);

  // Convert the matrix to a vector:
  to_vector(v,m);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
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
  binfile_write_vector( filename, fid, v, "Vector" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
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
  binfile_read_vector( v, filename, fid, "Vector" );
  binfile_close( fid, filename );
}



//=== Matrix ==========================================================

void MatrixWriteAscii(// WS Generic Input:
		       const Matrix& m,
		       // WS Generic Input Names:
		       const String& m_name,
		       // Control Parameters:
		       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, m_name );

  // Convert the matrix to an array of matrix:
  ArrayofMatrix am(1,m);

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void MatrixReadAscii(// WS Generic Output:
			Matrix& m,
			// WS Generic Output Names:
			const String& m_name,
			// Control Parameters:
			const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, m_name );

  // Read the array of matrix from the file:
  ArrayofMatrix am;
  read_array_of_matrix_from_file(am,filename);

  //  cout << "am.size(): " << am.size() << "\n";

  // Convert the array of matrix to a matrix.
  if ( 1 != am.size() )
   throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                       "but the dimension of the array is not 1.");

  m = am[0];
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
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
  binfile_write_matrix( filename, fid, v, "Matrix" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
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
  binfile_read_matrix( v, filename, fid, "Matrix" );
  binfile_close( fid, filename );
}



//=== ArrayofIndex =====================================================

// This function shall be modified to handle ArrayofIndex
void ArrayOfIndexWriteAscii(
        const Arrayofsizet&   v,
        const String&         v_name,
	const String&         f )
{
  String filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ArrayofMatrix and write to file
  const size_t  n = v.size();
  ArrayofMatrix am(1);
  resize(am[0],n,1);
  for ( size_t i=0; i<n; i++ )
    am[0][i][0] = (Numeric) v[i];
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle Index
void ArrayOfIndexReadAscii(
	      Arrayofsizet&   v,
        const String&         v_name,
	const String&         f )
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ArrayofMatrix 
  ArrayofMatrix am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.size()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains more than one vector.";
    throw runtime_error(os.str());
  }

  // Move values to v using a vector
  Vector x;
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfIndexWriteBinary(
        const Arrayofsizet&  v,
        const String&        var_name,
	const String&        f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_indexarray( filename, fid, v, "IndexArray" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfIndexReadBinary(
	      Arrayofsizet&  v,
        const String&        var_name,
	const String&        f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_indexarray( v, filename, fid, "IndexArray" );
  binfile_close( fid, filename );
}



//=== ArrayofVector ====================================================

void ArrayOfVectorWriteAscii(// WS Output:
			      const ArrayofVector& av,
			      // WS Variable Names:
			      const String& av_name,
			      // Control Parameters:
			      const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, av_name );

  // Convert the array of vector to an array of matrix:
  ArrayofMatrix am(av.size());
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
			       ArrayofVector& av,
			       // WS Generic Output Names:
			       const String& av_name,
			       // Control Parameters:
			       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, av_name );

  // Read an array of matrix from the file:
  ArrayofMatrix am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to an array of vector.
  resize(av,am.size());
  for (size_t i=0; i<am.size(); ++i)
    {
      to_vector(av[i],am[i]);
    }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfVectorWriteBinary(
        const ArrayofVector&  v,
        const String&         var_name,
	const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, v, "Vector" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfVectorReadBinary(
	      ArrayofVector&  v,
        const String&         var_name,
	const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_vectorarray( v, filename, fid, "Vector" );
  binfile_close( fid, filename );
}



//=== ArrayofMatrix ====================================================

void ArrayOfMatrixWriteAscii(// WS Generic Input:
			      const ArrayofMatrix& am,
			      // WS Generic Input Names:
			      const String& am_name,
			      // Control Parameters:
			      const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, am_name );

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void ArrayOfMatrixReadAscii(// WS Generic Output:
			       ArrayofMatrix& am,
			       // WS Generic Output Names:
			       const String& am_name,
			       // Control Parameters:
			       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, am_name );

  // Read the array of matrix from the file:
  read_array_of_matrix_from_file(am,filename);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfMatrixWriteBinary(
        const ArrayofMatrix&  v,
        const String&         var_name,
	const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_out( fid, filename );
  binfile_write_matrixarray( filename, fid, v, "Matrix" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfMatrixReadBinary(
	      ArrayofMatrix&  v,
        const String&         var_name,
	const String&         f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );
  binfile_open_in( fid, filename );
  binfile_read_matrixarray( v, filename, fid, "Matrix" );
  binfile_close( fid, filename );
}



//=== STRING ===============================================================

void StringWriteAscii( // WS Generic Input:
		       const String& s,
		       // WS Generic Input Names:
		       const String& s_name,
		       // Control Parameters:
		       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, s_name );

  // Convert the String to an array of String:
  ArrayofString as(1);
  as[0] = s;

  // Write the array of matrix to the file.
  write_array_of_String_to_file(filename,as);
}



void StringReadAscii(   // WS Generic Output:
			String& s,
			// WS Generic Output Names:
			const String& s_name,
			// Control Parameters:
			const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, s_name );

  // Read the array of matrix from the file:
  ArrayofString as;
  read_array_of_String_from_file(as,filename);

  // Convert the array of String to a String.
  if ( 1 != as.size() )
   throw runtime_error("You tried to convert an array of String to a String,\n"
                       "but the dimension of the array is not 1.");

  s = as[0];
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
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
   \date   2000-?-?
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



//=== ArrayofString ====================================================

void ArrayOfStringWriteAscii( // WS Generic Input:
			      const ArrayofString& as,
			      // WS Generic Input Names:
			      const String& as_name,
			      // Control Parameters:
			      const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, as_name );

  // Write the array of matrix to the file.
  write_array_of_String_to_file(filename,as);
}



void ArrayOfStringReadAscii(   // WS Generic Output:
			       ArrayofString& as,
			       // WS Generic Output Names:
			       const String& as_name,
			       // Control Parameters:
			       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, as_name );

  // Read the array of matrix from the file:
  read_array_of_String_from_file(as,filename);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfStringWriteBinary(
        const ArrayofString&  v,
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
   \date   2000-?-?
*/
void ArrayOfStringReadBinary(
	      ArrayofString&  v,
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



//=== SYMMETRIC ==========================================================

void SymmetricWriteAscii(// WS Generic Input:
		       const SYMMETRIC& m,
		       // WS Generic Input Names:
		       const String& m_name,
		       // Control Parameters:
		       const String& f)
{
  const Index   n = m.nrows();
  
  Matrix   x(n,n);

  copy( m, x );

  MatrixWriteAscii( x, m_name, f );
}



void SymmetricReadAscii(// WS Generic Output:
			SYMMETRIC& m,
			// WS Generic Output Names:
			const String& m_name,
			// Control Parameters:
			const String& f)
{
  Matrix   x;

  MatrixReadAscii( x, m_name, f );

  const Index   n = x.nrows();
  
  assert( n == x.ncols() );

  copy( x, m );
}



void SymmetricWriteBinary(// WS Generic Input:
		       const SYMMETRIC& m,
		       // WS Generic Input Names:
		       const String& m_name,
		       // Control Parameters:
		       const String& f)
{
  const Index   n = m.nrows();
  
  Matrix   x(n,n);

  copy( m, x );

  MatrixWriteBinary( x, m_name, f );
}



void SymmetricReadBinary(// WS Generic Output:
			SYMMETRIC& m,
			// WS Generic Output Names:
			const String& m_name,
			// Control Parameters:
			const String& f)
{
  Matrix   x;

  MatrixReadBinary( x, m_name, f );

  const Index   n = x.nrows();
  
  assert( n == x.ncols() );

  copy( x, m );
}



//=== MAYBESPARSE ====================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void HmatrixReadAscii(// WS Generic Output:
			Hmatrix& h,
			// WS Generic Output Names:
			const String& h_name,
			// Control Parameters:
			const String& f)
{
  String filename = f;
 
  // Create default filename if empty  
  filename_ascii( filename, h_name );

  // Read the array of matrix from the file:
  ArrayofMatrix am;
  read_array_of_matrix_from_file(am,filename);

  h.issparse = 0;
  h.full     = am[0];
  // Remember to clear the other field
}



//=== LOS ==================================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void LosWriteBinary(
        const LOS&      los,
        const String&   var_name,
	const String&   f )
{
  int    fid;
  String filename = f;
  filename_bin( filename, var_name );

  binfile_open_out( fid, filename );
  binfile_write_vectorarray( filename, fid, los.p,      "LOS.P" );
  binfile_write_vectorarray( filename, fid, los.psi,    "LOS.PSI" );
  binfile_write_vectorarray( filename, fid, los.z,      "LOS.Z" );
  binfile_write_vector(      filename, fid, los.l_step, "LOS.L_STEP" );
  binfile_write_indexarray(  filename, fid, los.ground, "LOS.GROUND" );
  binfile_write_indexarray(  filename, fid, los.start,  "LOS.START" );
  binfile_write_indexarray(  filename, fid, los.stop,   "LOS.STOP" );
  binfile_close( fid, filename );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void LosReadBinary(
	      LOS&      los,
        const String&   var_name,
	const String&   f )
{
  int   fid;
  String filename = f;
  filename_bin( filename, var_name );

  binfile_open_in( fid, filename );
  binfile_read_vectorarray( los.p,      filename, fid, "LOS.P" );
  binfile_read_vectorarray( los.psi,    filename, fid, "LOS.PSI" );
  binfile_read_vectorarray( los.z,      filename, fid, "LOS.Z" );
  binfile_read_vector(      los.l_step, filename, fid, "LOS.L_STEP" );
  binfile_read_indexarray(  los.ground, filename, fid, "LOS.GROUND" );
  binfile_read_indexarray(  los.start,  filename, fid, "LOS.START" );
  binfile_read_indexarray(  los.stop,   filename, fid, "LOS.STOP" );
  binfile_close( fid, filename );
}





//**************************************************************************
//
//   3. Creation by workspace method keywords
//
//**************************************************************************

//=== Index ============================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void IntSet(// WS Generic Output:
            int& x,
            // WS Generic Output Names:
            const String& x_name,
            // Control Parameters:
            const int& value)
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



//=== NUMERIC ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void NumericSet(// WS Generic Output:
                Numeric& x,
                // WS Generic Output Names:
                const String& x_name,
                // Control Parameters:
                const Numeric& value)
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-03-28
*/
void NumericCopyFirstOfVector(
                // WS Generic Output:
                      Numeric&  x,
                // WS Generic Output Names:
                const String&   x_name,
                // Control Parameters:
                const Vector&   v,
                const String&   v_name )
{
  x = v[0];
  out3 << "  Setting " << x_name << " to the first value of " << v_name << ".\n";
}



//=== Vector ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorSet(           Vector&  x, 
                    const String&  x_name,
                    const int&     n,
                    const Numeric& value )
{
  resize(x,n);
  setto(x,value);
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorSetLengthFromVector(
              Vector&  x, 
        const String&  x_name,
        const Vector&  z,
        const String&  z_name,
        const Numeric& value )
{
  const size_t  n = z.size();
  resize(x,n);
  setto(x,value);
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorLinSpace(      Vector&  x, 
                    const String&  x_name,
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorNLinSpace(     Vector&  x, 
                    const String&  x_name,
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorNLogSpace(     Vector&  x, 
                    const String&  x_name,
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
                      Vector&   y2,
                const String&   name_y2,
                const Vector&   y1,
                const String&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  resize( y2, y1.size() );
  copy( y1, y2 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-06-12
*/void VectorFlip(
                      Vector&   y2,
                const String&   name_y2,
                const Vector&   y1,
                const String&   name_y1 )
{
  out2 << "  Flips " << name_y2 << " to create " << name_y1 << "\n";

  Index n = y1.size();

  Vector dum( n );
  for ( Index i=0; i<n; i++ )
    dum[n-1-i] = y1[i];

  resize( y2, n );
  copy( dum, y2 );
}



void VectorCopyFromArrayOfVector(
                      Vector&          v,
                const String&          v_name,
                const ArrayofVector&   va,
                const String&          va_name,
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorPlanck(
                    Vector&   y,
              const String&   y_name,
              const Vector&   f,
              const String&   f_name,
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorCalcLog10(
                    Vector&   out,
              const String&   out_name,
              const Vector&   in,
              const String&   in_name )
{
  out2<<"  " << out_name << " = log10( " << in_name << " )\n";

  resize( out, in.size() );
  out = transf( in, log10 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-17
*/
void VectorAdd(
                    Vector&   out,
              const String&   out_name,
              const Vector&   in,
              const String&   in_name,
              const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " + " << value << "\n";

  const size_t  n = in.size();

  // Note that in and out can be the same vector
  Vector dummy(n);
  copy( in, dummy );
  resize( out, n );
  add ( dummy, Vector(n,value), out );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-17
*/
void VectorScale(
                    Vector&   out,
              const String&   out_name,
              const Vector&   in,
              const String&   in_name,
              const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  const size_t  n = in.size();

  // Note that in and out can be the same vector
  Vector dummy(n);
  copy( in, dummy );
  resize( out, n );
  for ( Index i=0; i<n; i++ )
    out[i] = dummy[i] * value;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorRandUniform(
                    Vector&   y,
              const String&   y_name,
              const Numeric&  x_low,
              const Numeric&  x_high,
              const int&      n )
{
  out2<<"  Filling " << y_name << " with uniform random data.\n";
  resize( y, n );
  rand_uniform( y, x_low, x_high );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorRandGaussian(
                    Vector&   y,
              const String&   y_name,
              const Numeric&  stddev,
              const int&      n )
{
  out2<<"  Filling " << y_name << " with Gaussian random data.\n";
  resize( y, n );
  rand_gaussian( y, stddev );
}



//=== Matrix ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-02-21
*/
void MatrixSet(           Matrix&  x, 
                    const String&  x_name,
                    const int&     nrows,
                    const int&     ncols,
                    const Numeric& value )
{
  resize( x, nrows, ncols );
  setto( x, value );
  out2 << "  Creating " << x_name << " as a constant matrix\n"; 
  out3 << "          nrows : " << nrows << "\n";
  out3 << "          ncols : " << ncols << "\n";
  out3 << "          value : " << value << "\n";
}



void MatrixCopy(
              Matrix&   y2,
        const String&   name_y2,
        const Matrix&   y1,
        const String&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  resize( y2, y1.nrows(), y1.ncols() );
  copy( y1, y2 );
}



void MatrixFillWithVector(
              Matrix&   m,
        const String&   name_m,
        const Vector&   y,
        const String&   name_y,
        const int&      n )
{
  out2 << "  Creates" << name_m << " by copying " << name_y << n << "times.\n";
  const size_t nrows = y.size();
  resize( m, nrows, n );
  for ( Index row=0; row<nrows; row++ ) 
    mtl::set(m[row],y[row]);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-17
*/
void MatrixScale(
                    Matrix&   out,
              const String&   out_name,
              const Matrix&   in,
              const String&   in_name,
              const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  const size_t  n1 = in.nrows();
  const size_t  n2 = in.ncols();

  // Note that in and out can be the same matrix
  Matrix dummy(n1,n2);
  copy( in, dummy );
  resize( out, n1, n2 );
  size_t i,j;
  for ( i=0; i<n1; i++ )
  {
    for ( j=0; j<n2; j++ )
      out[i][j] = dummy[i][j] * value;
  }  
}



//=== STRING ===============================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void StringSet(           String&  s, 
                    const String&  s_name,
                    const String&  s2 )
{
  s = s2;
  out3 << "  Setting " << s_name << " to " << s2 << "\n"; 
}



//=== ArrayofSTRING ========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfStringSet(    
              ArrayofString&  sa, 
        const String&         sa_name,
        const ArrayofString&  sa2 )
{
  resize(sa,sa2.size());
  copy(sa2,sa);
  out3 << "  Setting " << sa_name << "\n"; 
}



//=== SYMMETRIC ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-02-21
*/
void SymmetricDiagonal(
                          SYMMETRIC&  x, 
                    const String&     x_name,
                    const int&        nrows,
                    const Numeric&    value )
{
  resize( x, nrows, nrows );
  for ( Index i=0; i<Index(nrows); i++ )
    x[i][i] = value;

  out2 << "  Creating " << x_name << " as a diagonal (symmetric) matrix\n"; 
  out3 << "          nrows : " << nrows << "\n";
  out3 << "          value : " << value << "\n";
}

void SymmetricCopy(
           SYMMETRIC&   y2,
        const String&   name_y2,
     const SYMMETRIC&   y1,
        const String&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  resize( y2, y1.nrows(), y1.ncols() );
  copy( y1, y2 );
}
