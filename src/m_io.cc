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
/*!
  \file   m_io.cc
  \brief  Workspace functions for reading and writing data from/to files.

  The file includes now the ASCII functions from ARTS-1, but these will be
  replaced with XML-type functions, both for ASCII and binary.  

  The functions are sorted in alphabetical order.

  \author Patrick Eriksson
  \date 2002-05-08 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "arts.h"
#include "file.h"
#include "math_funcs.h"
#include "messages.h"
#include "auto_md.h"
#include "make_array.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////


//**************************************************************************
//
//  2. ASCII file functions
//
//**************************************************************************

////////////////////////////////////////////////////////////////////////////
//   File workspace methods (sorted after workspace variable)
////////////////////////////////////////////////////////////////////////////

// The ASCII functions are created by Stefan Buehler


//=== Index ============================================================

// This function shall be modified to handle Index
void IndexWriteAscii(
        const Index&      v,
        const String&   v_name,
	const String&   f )
{
  String filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ArrayOfMatrix and write to file
  ArrayOfMatrix am(1);
  am[0].resize( 1, 1 );
  am[0] = static_cast<Numeric>(v); // Matpack can set all elements like this.
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle Index
void IndexReadAscii(
	      Index&      v,
        const String&   v_name,
	const String&   f )
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ArrayOfMatrix and move to Numeric
  ArrayOfMatrix am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.nelem()!=1) || (am[0].nrows()!=1) || (am[0].ncols()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single value.";
    throw runtime_error(os.str());
  }

  Numeric a = am[0](0,0);

  if ( (a-floor(a)) != 0 )
    throw runtime_error("The value in the file is not an integer.");
  if ( a < 0 )
    throw runtime_error("The value in the file is negative.");

  v = (int) a;
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

  // Store the value in an ArrayOfMatrix and write to file
  ArrayOfMatrix am(1);
  am[0].resize( 1, 1 );
  am[0] = v;			// Matpack can set all elements like this.
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

  // Read the value as ArrayOfMatrix and move to Numeric
  ArrayOfMatrix am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.nelem()!=1) || (am[0].nrows()!=1) || (am[0].ncols()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains not a single numeric value";
    throw runtime_error(os.str());
  }

  v = am[0](0,0);
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

  // Convert the vector to a 1-column matrix:
  Matrix m(static_cast<const Vector>(v));
  // The static_cast to const Vector here is necessary to suppress a
  // warning message from the compiler about different possible
  // conversion paths. 

  // Convert the matrix to an array of matrix:
  MakeArray<Matrix> am(m);
  // MakeArray is the special kind of array with explicit
  // initialization. Here, the generated array has only 1 element,
  // which is initialized from Matrix m.

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
  ArrayOfMatrix am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to a matrix.
  if ( 1 != am.nelem() )
   throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                       "but the dimension of the array is not 1.");
  Matrix m(am[0]);

  // Check that this really is a 1-column matrix:
  if ( 1 != m.ncols() )
    throw runtime_error("You tried to convert a matrix to a vector,\n"
			"but it has more than one column.");

  // Convert the 1-column matrix to a vector:
  v.resize(m.nrows());
  v = m(Range(joker),0);	
  // (The m(Range(joker),0) picks out first column of m, = operator copies
  // to v.)
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
  MakeArray<Matrix> am(m);
  // MakeArray is the special kind of array with explicit
  // initialization. Here, the generated array has only 1 element,
  // which is initialized from Matrix m.

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
  ArrayOfMatrix am;
  read_array_of_matrix_from_file(am,filename);

  //  cout << "am.nelem(): " << am.nelem() << "\n";

  // Convert the array of matrix to a matrix.
  if ( 1 != am.nelem() )
   throw runtime_error("You tried to convert an array of matrix to a matrix,\n"
                       "but the dimension of the array is not 1.");

  m.resize( am[0].nrows(), am[0].ncols() );
  m = am[0];
}



//=== ArrayOfIndex =====================================================

// This function shall be modified to handle ArrayOfIndex
void ArrayOfIndexWriteAscii(
        const ArrayOfIndex&   v,
        const String&         v_name,
	const String&         f )
{
  String filename = f;

  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Store the value in an ArrayOfMatrix and write to file
  const Index  n = v.nelem();
  ArrayOfMatrix am(1);
  am[0].resize(n,1);
  for ( Index i=0; i<n; i++ )
    am[0](i,0) = static_cast<Numeric>(v[i]);
  write_array_of_matrix_to_file(filename,am);
}


// This function shall be modified to handle Index
void ArrayOfIndexReadAscii(
			   ArrayOfIndex&   v,
			   const String&   v_name,
			   const String&   f )
{
  // FIXME: This function is crap. Put the whole ASCII file stuff
  // should be changed in the future, so I leave it for now.
  
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, v_name );

  // Read the value as ArrayOfMatrix 
  ArrayOfMatrix am;
  read_array_of_matrix_from_file(am,filename);
  if ( (am.nelem()!=1) )
  {
    ostringstream os;
    os << "The file " << filename << " contains more than one vector.";
    throw runtime_error(os.str());
  }
  
  Matrix m(am[0]);		// Create Matrix and initialize from
				// first element of am.

  // Check that this really is a 1-column matrix:
  if ( 1 != m.ncols() )
    throw runtime_error("You tried to convert a matrix to a vector,\n"
			"but it has more than one column.");

  // Convert the 1-column matrix to a vector:
  Vector x(m(Range(joker),0));
  // (The m(Range(joker),0) picks out first column of m. Matpack can
  // initialize the new vector x from that.

  const Index  n = x.nelem();
  v.resize(n);
  for ( Index i=0; i<n; i++ )
  {  
    if ( (x[i]-floor(x[i])) != 0 )
      throw runtime_error("A value in the file is not an integer.");
    if ( x[i] < 0 )
      throw runtime_error("A value in the file is negative.");
    v[i] = (Index) x[i];
  }
}



//=== ArrayOfVector ====================================================

void ArrayOfVectorWriteAscii(// WS Output:
			      const ArrayOfVector& av,
			      // WS Variable Names:
			      const String& av_name,
			      // Control Parameters:
			      const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, av_name );

  // Convert the array of vector to an array of matrix:
  ArrayOfMatrix am(av.nelem());
  for (Index i=0; i<av.nelem(); ++i)
    {
      //      to_matrix(am[i],av[i]);
      am[i].resize( av[i].nelem(), 1 );
      am[i] = av[i];		// Matpack can copy the content of a
				// Vector to a 1-column Matrix.
    }

  // Write the array of matrix to the file.
  write_array_of_matrix_to_file(filename,am);
}



void ArrayOfVectorReadAscii(// WS Generic Output:
			       ArrayOfVector& av,
			       // WS Generic Output Names:
			       const String& av_name,
			       // Control Parameters:
			       const String& f)
{
  String filename = f;
  
  // Create default filename if empty  
  filename_ascii( filename, av_name );

  // Read an array of matrix from the file:
  ArrayOfMatrix am;
  read_array_of_matrix_from_file(am,filename);

  // Convert the array of matrix to an array of vector.
  av.resize(am.nelem());
  for (Index i=0; i<am.nelem(); ++i)
    {
      // Check that this really is a 1-column matrix:
      if ( 1 != am[i].ncols() )
	throw runtime_error("You tried to convert a matrix to a vector,\n"
			    "but it has more than one column.");

      // Convert the 1-column matrix to a vector:
      av[i].resize(am[i].nrows());
      av[i] = am[i](Range(joker),0);	
      // (The am[i](Range(joker),0) picks out first column of am[i], = operator copies
      // to v.)
    }
}



//=== ArrayOfMatrix ====================================================

void ArrayOfMatrixWriteAscii(// WS Generic Input:
			      const ArrayOfMatrix& am,
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
			       ArrayOfMatrix& am,
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
  ArrayOfString as(1);
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
  ArrayOfString as;
  read_array_of_String_from_file(as,filename);

  // Convert the array of String to a String.
  if ( 1 != as.nelem() )
   throw runtime_error("You tried to convert an array of String to a String,\n"
                       "but the dimension of the array is not 1.");

  s = as[0];
}



//=== ArrayOfString ====================================================

void ArrayOfStringWriteAscii( // WS Generic Input:
			      const ArrayOfString& as,
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
			       ArrayOfString& as,
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



