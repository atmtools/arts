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
  \file   m_basic_types.cc
  \brief  Workspace functions for straightforward operations on variables 
          of basic types.

  This file includes workspace functions for variables of basic types, 
  such as Matrix and ArrayOfIndex. The functions are mainly of two types:
  1. Initiation of variables by keyword arguments, such as *StringSet*.
  2. Basic math, such as *MatrixVectorMultiply*.

  The functions are sorted in alphabetical order.

  \author Patrick Eriksson
  \date 2002-05-08 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "make_array.h"
#include "matpackI.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void IndexSet(    Index&    x,
            const String&   x_name,
            const Index&    value )
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void NumericSet(      Numeric&   x,
                const String&    x_name,
                const Numeric&   value )
{
  x = value;
  out3 << "  Setting " << x_name << " to " << value << ".\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author ???
   \date   ????-?-?
*/
void VectorCopy(      Vector&   y2,
                const String&   name_y2,
                const Vector&   y1,
                const String&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2.resize( y1.nelem() );
  y2 = y1;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorLinSpace(      Vector&    x, 
                    const String&    x_name,
                    const Numeric&   start,
                    const Numeric&   stop,
                    const Numeric&   step )
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
void VectorNLinSpace(     Vector&    x, 
                    const String&    x_name,
                    const Numeric&   start,
                    const Numeric&   stop,
                    const Index&     n )
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorNLogSpace(       Vector&    x, 
		      const String&    x_name,
		      const Numeric&   start,
		      const Numeric&   stop,
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



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorSet(           Vector&    x, 
                    const String&    x_name,
                    const Index&     n,
                    const Numeric&   value )
{
  x.resize(n);
  x = value;		
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void VectorSetTakingLengthFromVector(
              Vector&  	 x, 
        const String&  	 x_name,
        const Vector&  	 z,
        const String&  	 z_name,
        const Numeric& 	 value )
{
  const Index  n = z.nelem();
  x.resize(n);
  x = value;		
  out2 << "  Creating " << x_name << " as a constant vector\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
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

