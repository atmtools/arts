/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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
#include "math_funcs.h"
#include "matpackI.h"
#include "messages.h"



////////////////////////////////////////////////////////////////////////////
//   The functions
////////////////////////////////////////////////////////////////////////////

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
  sa = sa2;
  out3 << "  Setting " << sa_name << "\n"; 
}



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

   \author ???
   \date   ????-??-??
*/
void MatrixCopy(
              Matrix&   y2,
        const String&   name_y2,
        const Matrix&   y1,
        const String&   name_y1 )
{
  out2 << "  " << name_y2 << " = " << name_y1 << "\n";
  y2.resize( y1.nrows(), y1.ncols() );
  y2 = y1;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-11
*/
void MatrixFillWithVector(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v,
        // WS Generic Input Names:
        const String&   v_name,
        // Control Parameters:
        const Index&    nrows,
        const Index&    ncols )
{
  if( nrows==0 && ncols==0 )
    throw runtime_error("0nly one of the keyword arguments can be 0."); 

  const Index nv = v.nelem();

  if( ncols == 0 )
    {
      out2 << "  Creates the matrix " << m_name << " by putting in " << v_name
           << " as rows.\n";
      out3 << "          nrows : " << nrows << "\n";
      out3 << "          ncols : " << nv << "\n";
      m.resize( nrows, nv );
      for( Index i=0; i<nrows; i++ )
        m(i,Range(joker)) = v;
    }
  else if( nrows == 0 )
    {
      out2 << "  Creates the matrix " << m_name << " by putting in " << v_name
           << " as columns.\n";
      out3 << "          nrows : " << nv << "\n";
      out3 << "          ncols : " << ncols << "\n";
      m.resize( nv, ncols );
      for( Index i=0; i<ncols; i++ )
        m(Range(joker),i) = v;
    }
  else 
    throw runtime_error(
                    "The size argument for either rows or columns must be 0.");
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
      out *= value;  
    }
  else
    {
      // Out and in are different. We first have to copy in to out,
      // then multiply by the scalar value.
      out.resize( in.nrows(), in.ncols() );
      out = in;	
      out *= value;
    }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-02-21
*/
void MatrixSet(           Matrix&    x, 
                    const String&    x_name,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( nrows, ncols );
  x = value;
  out2 << "  Creating " << x_name << " as a constant matrix\n"; 
  out3 << "          nrows : " << nrows << "\n";
  out3 << "          ncols : " << ncols << "\n";
  out3 << "          value : " << value << "\n";
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



void Tensor3FillWithVector(
        // WS Generic Output:
              Tensor3&   t,
        // WS Generic Output Names:
        const String&   t_name,
        // WS Generic Input:
        const Vector&   v,
        // WS Generic Input Names:
        const String&   v_name,
        // Control Parameters:
        const Index& npages,
        const Index& nrows,
        const Index& ncols )
{
  if( ( (npages==0) + (nrows==0) + (ncols==0) ) > 1 )
    throw runtime_error("0nly one of the keyword arguments can be 0."); 

  const Index nv = v.nelem();

  if( ncols == 0 )
    {
      out2 << "  Creates the tensor " << t_name << " by putting in " << v_name
           << "\n  perpendicular to the column dimension.\n";
      out3 << "          npages : " << npages << "\n";
      out3 << "          nrows  : " << nrows << "\n";
      out3 << "          ncols  : " << nv << "\n";
      t.resize( npages, nrows, nv );
      for( Index i=0; i<npages; i++ )
	{
	  for( Index j=0; j<nrows; j++ )
	    t(i,j,Range(joker)) = v;
	}
    }
  else if( nrows == 0 )
    {
      out2 << "  Creates the tensor " << t_name << " by putting in " << v_name
           << "\n  perpendicular to the row dimension.\n";
      out3 << "          npages : " << npages << "\n";
      out3 << "          nrows  : " << nv << "\n";
      out3 << "          ncols  : " << ncols << "\n";
      t.resize( npages, nv, ncols );
      for( Index i=0; i<npages; i++ )
	{
	  for( Index j=0; j<ncols; j++ )
	    t(i,Range(joker),j) = v;
	}
    }
  else if( npages == 0 )
    {
      out2 << "  Creates the tensor " << t_name << " by putting in " << v_name
           << "\n  perpendicular to the page dimension.\n";
      out3 << "          npages : " << nv << "\n";
      out3 << "          nrows  : " << nrows << "\n";
      out3 << "          ncols  : " << ncols << "\n";
      t.resize( nv, nrows, ncols );
      for( Index i=0; i<nrows; i++ )
	{
	  for( Index j=0; j<ncols; j++ )
	    t(Range(joker),i,j) = v;
	}
    }
  else 
    throw runtime_error(
             "The size argument for either pages, rows or columns must be 0.");
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-17
*/
void VectorAddScalar(
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
      out += value;	
    }
  else
    {
      // Out and in are different. We first have to copy in to out,
      // then add the scalar value.
      out.resize( in.nelem() );
      out = in;
      out += value;
    }
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
  out2 << "  Creating " << x_name << " as linearly spaced vector.\n";
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
  out2 << "  Creating " << x_name << " as linearly spaced vector.\n";
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
  out2 << "  Creating " << x_name << " as logarithmically spaced vector.\n";
  out3 << "         length: " << n << "\n";
  out3 << "    first value: " << x[0] << "\n";
  if ( x.nelem() > 1 )
    out3 << "     last value: " << x[x.nelem()-1] << "\n";
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
      out *= value;
    }
  else
    {
      // Out and in are different. We first have to copy in to out,
      // then multiply by the scalar value.
      out.resize( in.nelem() );
      out = in;
      out *= value;
    }
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
  out2 << "  Creating " << x_name << " as a constant vector.\n"; 
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
  out2 << "  Creating " << x_name << " as a constant vector.\n"; 
  out3 << "         length : " << n << "\n";
  out3 << "          value : " << value << "\n";
}

