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
   \file   m_io.cc

   This file contains IO workspace methods. 

   The functions are of two types:
     1. Functions with overall influence on ARTS
     2. Reading and writing to/from ASCII files.
     3. Creation by workspace method keywords and generic input

   \author Patrick Eriksson
   \date 2000-11-01 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "arts.h"
#include "atm_funcs.h"          
#include "file.h"
#include "math_funcs.h"
#include "messages.h"
#include "auto_md.h"
#include "make_array.h"



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
  exit(0);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-05-15
*/
void Test( )
{
  // This function can be used to test stuff.

  //  assert_bool( 2, "XCVF" );
  check_length_ncol( Vector(6), "XCVF", Matrix(5,5), "ERTFD" );
}




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
void IndexSet(// WS Generic Output:
            Index& x,
            // WS Generic Output Names:
            const String& x_name,
            // Control Parameters:
            const Index& value)
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
                    const Index&     n,
                    const Numeric& value )
{
  x.resize(n);
  x = value;			// Matpack can set all elements like this.
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
  const Index  n = z.nelem();
  x.resize(n);
  x = value;			// Matpack can set all elements like this.
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
  linspace(x,start,stop,step);
  out2 << "  Creating " << x_name << " as linearly spaced vector\n";
  out3 << "         length: " << x.nelem() << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.nelem() > 1 )
  {
    out3 << "      step size: " << x[1]-x[0] << "\n";
    out3 << "     last value: " << x[x.nelem()-1] << "\n";
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
                    const Index& n )
{
  if ( n<2 )
    throw runtime_error("The number of points must be > 1."); 
  nlinspace(x,start,stop,n);
  out2 << "  Creating " << x_name << " as linearly spaced vector\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.nelem() > 1 )
  {
    out3 << "      step size: " << x[1]-x[0] << "\n";
    out3 << "     last value: " << x[x.nelem()-1] << "\n";
  }
}



void VectorPressuresForLinAltitudes(
                     // WS Generic Output:
                     Vector&         p,
                     // WS Generic Output Names:
                     const String&   p_name,
                     // WS Input:
                     const Vector&   p_abs,
                     const Vector&   z_abs,
                     // Control Parameters:
                     const Numeric&  delta_z,
                     const Numeric&  p_start,
                     const Numeric&  p_stop)

{
  Vector p_lim(2), z_lim(2);
  p_lim[0] = p_start; 
  p_lim[1] = p_stop; 
  
  interpp(z_lim,p_abs,z_abs,p_lim);

  Vector z;
  linspace(z,z_lim[0],z_lim[1],delta_z);
  p.resize( z.nelem());
  z2p(p, z_abs, p_abs, z);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorNLogSpace( Vector&  x, 
		      const String&  x_name,
		      const Numeric& start,
		      const Numeric& stop,
		      const Index&     n )
{
  if ( n<2 )
    throw runtime_error("The number of points must be > 1."); 
  if ( (start<=0) || (stop<=0) )
    throw runtime_error("Only positive numbers are allowed."); 

  x.resize(n);
  x = nlogspace(start,stop,n);
  out2 << "  Creating " << x_name << " as logarithmically spaced vector\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.nelem() > 1 )
    out3 << "     last value: " << x[x.nelem()-1] << "\n";
}



void VectorCopy(
                      Vector&   y2,
                const String&   name_y2,
                const Vector&   y1,
                const String&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2.resize( y1.nelem() );
  y2 = y1;			// Matpack can copy the contents of
				// vectors like this. The dimensions
				// must be the same! 
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-06-12
*/
void VectorFlip(
                      Vector&   y2,
                const String&   name_y2,
                const Vector&   y1,
                const String&   name_y1 )
{
  out2 << "  Flips " << name_y2 << " to create " << name_y1 << "\n";

  Index n = y1.nelem();

  Vector dum( n );
  for ( Index i=0; i<n; i++ )
    dum[n-1-i] = y1[i];

  y2.resize( n );
  y2 = dum;			// Matpack can copy the contents of
				// vectors like this. The dimensions
				// must be the same! 
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
    y.resize( f.nelem() );
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

  out.resize( in.nelem() );
  transform( out, log10, in );	// out = log10(in)
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

  // Note that in and out can be the same vector
  if (&out==&in)
    {
      // Out and in are the same. Just add the scalar value.
      out += value;		// With Matpack you can add a scalar
				// to all elements of a vector like
				// this. 
    }
  else
    {
      // Out and in are different. We first have to copy in to out,
      // then add the scalar value.

      out.resize( in.nelem() );
      out = in;			// Matpack can copy the contents of
				// vectors like this. The dimensions
				// must be the same! 
      out += value;
    }
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

  // Note that in and out can be the same vector
  if (&out==&in)
    {
      // Out and in are the same. Just multiply by the scalar value.
      out *= value;		// With Matpack you can add a scalar
				// to all elements of a vector like
				// this. 
    }
  else
    {
      // Out and in are different. We first have to copy in to out,
      // then multiply by the scalar value.

      out.resize( in.nelem() );
      out = in;			// Matpack can copy the contents of
				// vectors like this. The dimensions
				// must be the same! 
      out *= value;
    }
}

/**
   Compute y = M*x.

   Works also if y and x are the same Vector.

   For more information see the the online help (arts -d
   FUNCTION_NAME).

   \author Stefan Buehler
   \date   2001-10-02
*/
void VectorMatrixMultiply(// WS Generic Output:
                          Vector& y,
                          // WS Generic Output Names:
                          const String& y_name,
                          // WS Generic Input:
                          const Matrix& M,
                          const Vector& x,
                          // WS Generic Input Names:
                          const String& M_name,
                          const String& x_name)
{
  // Check that dimensions are right, x must match columns of M:
  check_length_ncol( x, x_name, M, M_name );

  // Temporary for the result:
  Vector dummy( M.nrows() );

  mult( dummy, M, x );

  // Copy result to y:

  y.resize( dummy.nelem() );

  y = dummy;
}


//=== Matrix ==========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-02-21
*/
void MatrixSet(           Matrix&  x, 
                    const String&  x_name,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric& value )
{
  x.resize( nrows, ncols );
  x = value;			// Matpack can set all elements like this.
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
  y2.resize( y1.nrows(), y1.ncols() );
  y2 = y1;			// Matpack can copy the contents of
				// matrices like this. The dimensions
				// must be the same! 
}



void MatrixFillWithVector(
              Matrix&   m,
        const String&   name_m,
        const Vector&   y,
        const String&   name_y,
        const Index&      n )
{
  out2 << "  Creates" << name_m << " by copying " << name_y << n << "times.\n";
  m.resize( y.nelem(), n );
  for ( Index i=0; i<n; ++i ) 
    m(Range(joker),i) = y;	// Copy content of vector y to this
				// column of Matrix m.
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

  // Note that in and out can be the same matrix
  if (&out==&in)
    {
      // Out and in are the same. Just multiply by the scalar value.
      out *= value;		// With Matpack you can multiply a scalar
				// to all elements of a matrix like
				// this. 
    }
  else
    {
      // Out and in are different. We first have to copy in to out,
      // then multiply by the scalar value.

      out.resize( in.nrows(), in.ncols() );
      out = in;			// Matpack can copy the contents of
				// matrices like this. The dimensions
				// must be the same! 
      out *= value;
    }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-02-21
*/
void MatrixDiagonal(
		    Matrix&           x, 
                    const String&     x_name,
                    const Index&      nrows,
                    const Numeric&    value )
{
  x.resize( nrows, nrows );
  for ( Index i=0; i<Index(nrows); i++ )
    x(i,i) = value;

  out2 << "  Creating " << x_name << " as a diagonal matrix\n"; 
  out3 << "          nrows : " << nrows << "\n";
  out3 << "          value : " << value << "\n";
}

/**
   Compute Y = M*X.

   Works also if Y and X are the same Matrix.

   For more information see the the online help (arts -d
   FUNCTION_NAME).

   \author Stefan Buehler
   \date   2001-10-02
*/
void MatrixMatrixMultiply(// WS Generic Output:
                          Matrix& Y,
                          // WS Generic Output Names:
                          const String& Y_name,
                          // WS Generic Input:
                          const Matrix& M,
                          const Matrix& X,
                          // WS Generic Input Names:
                          const String& M_name,
                          const String& X_name)
{
  // Check that dimensions are right, M.ncols() must match X.nrows():
  check_ncol_nrow( M, M_name, X, X_name );

  // Temporary for the result:
  Matrix dummy( M.nrows(), X.ncols() );

  mult( dummy, M, X );

  // Copy result to Y:

  Y.resize( dummy.nrows(), dummy.ncols() );

  Y = dummy;
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



//=== ArrayOfSTRING ========================================================

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void ArrayOfStringSet(    
              ArrayOfString&  sa, 
        const String&         sa_name,
        const ArrayOfString&  sa2 )
{
  sa.resize(sa2.nelem());
  sa = sa2;			// Arrays can be copied like this.
  out3 << "  Setting " << sa_name << "\n"; 
}

