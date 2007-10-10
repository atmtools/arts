/* Copyright (C) 2002-2007
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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



/*===========================================================================
  ===  File description
  ===========================================================================*/

/*!
  \file   m_basic_types.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions for straightforward operations on variables 
          of basic types.

  This file includes workspace functions for variables of basic types, 
  such as Matrix and ArrayOfIndex. The functions are mainly of two types: <br>
  1. Initiation of variables by keyword arguments, such as *StringSet*. <br>
  2. Basic math, such as *MatrixVectorMultiply*.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "arts.h"
#include "array.h"
#include "matpackI.h"
#include "matpackIII.h"
#include "matpackIV.h"
#include "matpackV.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "mystring.h"
#include "exceptions.h"
#include "make_array.h"
#include "math_funcs.h"
#include "messages.h"
#include "logic.h"
#include "sorting.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfMatrixSet(// WS Generic Output:
                      ArrayOfMatrix&  aom,
                      // WS Generic Output Names:
                      const String&   aom_name,
                      // WS Generic Input:
                      const Matrix&   m,
                      // WS Generic Input Names:
                      const String&   m_name,
                      // Control Parameters:
                      const Index&    element )
{
  // Check input index, if larger than number of elements in
  // the array, return error message
  if (element>aom.nelem()) {
    ostringstream os;
    os << "The element index "<<element<<" is too large, there are only "
       << aom.nelem() <<" elements in "<<aom_name<<".\n";
    throw runtime_error(os.str());
  } else if (element<0 || element==aom.nelem()) {
    aom.push_back(m);
    out2 << "  Appending Matrix " << m_name << " to ArrayOfMatrix "
       << aom_name << ".\n";
  } else {
    aom[element] = m;
    out2 << "  Setting element "<<element<<" in ArrayOfMatrix "
         << aom_name<<" to Matrix "<<m_name<<".\n";
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfStringSet(
              ArrayOfString&  sa,
        const String&         sa_name,
        const ArrayOfString&  sa2 )
{
  sa.resize(sa2.nelem());
  sa = sa2;
  out2 << "  Setting " << sa_name << " to the given string array.\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FlagOff(    
            Index&    x,
      const String&   x_name )
{
  x = 0;
  out2 << "  " << x_name << " = 0\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FlagOn(      
             Index&    x,
       const String&   x_name )
{
  x = 1;
  out2 << "  " << x_name << " = 1\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IndexSet(    Index&    x,
            const String&   x_name,
            const Index&    value )
{
  x = value;
  out3 << "  " << x_name << " = " << value << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexSetFromArrayOfMatrixLength(// WS Generic Output:
                                     Index& n,
                                     // WS Generic Output Names:
                                     const String& n_name,
                                     // WS Generic Input:
                                     const ArrayOfMatrix& am,
                                     // WS Generic Input Names:
                                     const String& am_name)
{
  n = am.nelem();
  out3 << "  Setting " << n_name << " = " << n
       << " from length of " << am_name << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void IndexStep(    
            Index&     xout,       
      const String&    outname,
      const Index&     xin,
      const String&    inname )
{
  xout = xin + 1;
  out3 << "  " << outname << " = " << inname << " + 1\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
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
  out3 << "  Multiplying " << M_name << " with " << X_name
       << " and storing the result in " << Y_name << ".\n";

  // Check that dimensions are right, M.ncols() must match X.nrows():
  if (M.ncols()!=X.nrows())
    {
      ostringstream os;
      os << "Matrix dimensions must be consistent!\n"
         << M_name << ".ncols() = " << M.ncols() << "\n"
         << X_name << ".nrows() = " << X.nrows();
      throw runtime_error( os.str() );
    }

  // Temporary for the result:
  Matrix dummy( M.nrows(), X.ncols() );

  mult( dummy, M, X );

  // Copy result to Y:

  Y.resize( dummy.nrows(), dummy.ncols() );

  Y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix1ColFromVector(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v,
        // WS Generic Input Names:
        const String&   v_name )
{
  const Index nv = v.nelem();

  out2 << "  Creates the matrix " << m_name << " by putting in " << v_name
       << " as only column.\n";
  m.resize( nv, 1 );
  m( joker,0 ) = v;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix2ColFromVectors(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v1,
        const Vector&   v2,
        // WS Generic Input Names:
        const String&   v1_name,
        const String&   v2_name )
{
  const Index nv = v1.nelem();

  if( v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  out2 << "  Creates the matrix " << m_name << " by putting in\n  "
       << v1_name << " as first column and\n  "
       << v2_name << " as second column.\n";
  m.resize( nv, 2 );
  m( joker,0 ) = v1;
  m( joker,1 ) = v2;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix3ColFromVectors(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v1,
    const Vector&   v2,
        const Vector&   v3,
    // WS Generic Input Names:
        const String&   v1_name,
        const String&   v2_name,
        const String&   v3_name )
{
  const Index nv = v1.nelem();

  if( v3.nelem() != nv || v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  out2 << "  Creates the matrix " << m_name << " by putting in\n  "
       << v1_name << " as first column,\n  "
       << v2_name << " as second column and\n  "
       << v3_name << " as third column.\n";
  m.resize( nv, 3 );
  m( joker,0 ) = v1;
  m( joker,1 ) = v2;
  m( joker,2 ) = v3;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix1RowFromVector(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v,
        // WS Generic Input Names:
        const String&   v_name )
{
  const Index nv = v.nelem();

  out2 << "  Creates the matrix " << m_name << " by putting in " << v_name
       << " as only row.\n";
  m.resize( 1, nv );
  m( 0, joker ) = v;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix2RowFromVectors(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v1,
        const Vector&   v2,
    // WS Generic Input Names:
        const String&   v1_name,
        const String&   v2_name )
{
  const Index nv = v1.nelem();

  if( v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  out2 << "  Creates the matrix " << m_name << " by putting in\n  "
       << v1_name << " as first row and\n  "
       << v2_name << " as second row.\n";
  m.resize( 2, nv );
  m( 0, joker ) = v1;
  m( 1, joker ) = v2;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix3RowFromVectors(
        // WS Generic Output:
              Matrix&   m,
        // WS Generic Output Names:
        const String&   m_name,
        // WS Generic Input:
        const Vector&   v1,
    const Vector&   v2,
        const Vector&   v3,
    // WS Generic Input Names:
        const String&   v1_name,
        const String&   v2_name,
        const String&   v3_name )
{
  const Index nv = v1.nelem();

  if( v3.nelem() != nv || v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  out2 << "  Creates the matrix " << m_name << " by putting in\n  "
       << v1_name << " as first row,\n  "
       << v2_name << " as second row and\n  "
       << v3_name << " as third row.\n";
  m.resize( 3, nv );
  m( 0, joker ) = v1;
  m( 1, joker ) = v2;
  m( 2, joker ) = v3;

}


/* Workspace method: Doxygen documentation will be auto-generated */
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


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixSet(           Matrix&    x, 
                    const String&    x_name,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( nrows, ncols );
  x = value;
  out2 << "  " << x_name << " = " << value << "\n"; 
  out3 << "             nrows : " << nrows << "\n";
  out3 << "             ncols : " << ncols << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericSet(      Numeric&   x,
                const String&    x_name,
                const Numeric&   value )
{
  x = value;
  out3 << "  " << x_name << " = " << value << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void StringSet(           String&  s, 
                    const String&  s_name,
                    const String&  s2 )
{
  s = s2;
  out3 << "  " << s_name << " = " << s2 << "\n"; 
}


/* Workspace method: Doxygen documentation will be auto-generated */
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
      out3 << "            npages : " << npages << "\n";
      out3 << "            nrows  : " << nrows << "\n";
      out3 << "            ncols  : " << nv << "\n";
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
      out3 << "            npages : " << npages << "\n";
      out3 << "            nrows  : " << nv << "\n";
      out3 << "            ncols  : " << ncols << "\n";
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
      out3 << "            npages : " << nv << "\n";
      out3 << "            nrows  : " << nrows << "\n";
      out3 << "            ncols  : " << ncols << "\n";
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


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3Scale(        Tensor3&  out,
                    const String&   out_name,
                    const Tensor3&  in,
                    const String&   in_name,
                    const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.npages(), in.nrows(), in.ncols() );
    out = in;
    out *= value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3Set(          Tensor3&   x, 
                    const String&    x_name,
                    const Index&     npages,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( npages, nrows, ncols );
  x = value;
  out2 << "  " << x_name << " = " << value  << "\n";
  out3 << "            npages : " << npages << "\n";
  out3 << "             nrows : " << nrows  << "\n";
  out3 << "             ncols : " << ncols  << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4Scale(        Tensor4&  out,
                    const String&   out_name,
                    const Tensor4&  in,
                    const String&   in_name,
                    const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.nbooks(), in.npages(), in.nrows(), in.ncols() );
    out = in;
    out *= value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4Set(          Tensor4&   x, 
                    const String&    x_name,
                    const Index&     nbooks,      
                    const Index&     npages,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( nbooks, npages, nrows, ncols );
  x = value;
  out2 << "  " << x_name << " = " << value  << "\n";
  out3 << "            nbooks : " << nbooks << "\n";
  out3 << "            npages : " << npages << "\n";
  out3 << "             nrows : " << nrows  << "\n";
  out3 << "             ncols : " << ncols  << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor5Scale(        Tensor5&  out,
                    const String&   out_name,
                    const Tensor5&  in,
                    const String&   in_name,
                    const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.nshelves(), in.nbooks(), in.npages(),
      in.nrows(), in.ncols() );
    out = in;
    out *= value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor5Set(          Tensor5&   x, 
                    const String&    x_name,
                    const Index&     nshelves,     
                    const Index&     nbooks,      
                    const Index&     npages,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( nshelves, nbooks, npages, nrows, ncols );
  x = value;
  out2 << "  " << x_name << " = " << value    << "\n";
  out3 << "          nshelves : " << nshelves << "\n";
  out3 << "            nbooks : " << nbooks   << "\n";
  out3 << "            npages : " << npages   << "\n";
  out3 << "             nrows : " << nrows    << "\n";
  out3 << "             ncols : " << ncols    << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor6Scale(        Tensor6&  out,
                    const String&   out_name,
                    const Tensor6&  in,
                    const String&   in_name,
                    const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.nvitrines(), in.nshelves(), in.nbooks(),
      in.npages(), in.nrows(), in.ncols() );
    out = in;
    out *= value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor6Set(          Tensor6&   x, 
                    const String&    x_name,
                    const Index&     nvitrines,   
                    const Index&     nshelves,     
                    const Index&     nbooks,      
                    const Index&     npages,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( nvitrines, nshelves, nbooks, npages, nrows, ncols );
  x = value;
  out2 << "  " << x_name << " = " << value     << "\n";
  out3 << "         nvitrines : " << nvitrines << "\n";    
  out3 << "          nshelves : " << nshelves  << "\n";
  out3 << "            nbooks : " << nbooks    << "\n";
  out3 << "            npages : " << npages    << "\n";
  out3 << "             nrows : " << nrows     << "\n";
  out3 << "             ncols : " << ncols     << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor7Scale(        Tensor7&  out,
                    const String&   out_name,
                    const Tensor7&  in,
                    const String&   in_name,
                    const Numeric&  value )
{
  out2<<"  " << out_name << " = " << in_name << " * " << value << "\n";

  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out *= value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.nlibraries(), in.nvitrines(), in.nshelves(),
      in.nbooks(), in.npages(), in.nrows(), in.ncols() );
    out = in;
    out *= value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor7Set(          Tensor7&   x, 
                    const String&    x_name,
                    const Index&     nlibraries,          
                    const Index&     nvitrines,   
                    const Index&     nshelves,     
                    const Index&     nbooks,      
                    const Index&     npages,
                    const Index&     nrows,
                    const Index&     ncols,
                    const Numeric&   value )
{
  x.resize( nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols );
  x = value;
  out2 << "  " << x_name << " = " << value      << "\n";
  out3 << "        nlibraries : " << nlibraries << "\n";
  out3 << "         nvitrines : " << nvitrines  << "\n";
  out3 << "          nshelves : " << nshelves   << "\n";
  out3 << "            nbooks : " << nbooks     << "\n";
  out3 << "            npages : " << npages     << "\n";
  out3 << "             nrows : " << nrows      << "\n";
  out3 << "             ncols : " << ncols      << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
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


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorInsertGridPoints(// WS Generic Output:
                            Vector& og,                  // Output grid
                            // WS Generic Output Names:
                            const String& og_name,
                            // WS Generic Input:
                            const Vector& ingrid,        // Input grid 
                            const Vector& points,        // Points to insert
                            // WS Generic Input Names:
                            const String& ingrid_name,
                            const String& points_name)
{

  out2 << "  Inserting points from *" << points_name
       << "* into grid in *" << ingrid_name
       << "*,\n"
       << "  and returning the result in *" << og_name << "*.\n";
    

  // First make duplikates of the input vectors, in case one of them
  // happens to be identical to the output vector. Also, we can fool
  // around with these, if we want.
  Vector ig(ingrid);
  Vector p(points);

  // Check how the input grid is sorted. If the grid is sorted in
  // descending order, we simply turn it around. (But don't
  // forget to turn it back at the end!) 
  Index ascending;              // 1=ascending, 0=descending
  if (is_increasing(ig))
    {
      ascending = 1;
    }
  else if (is_decreasing(ig))
    {
      ascending = 0;

      // Turn grid round.

      // Copy ig to dummy vector in reverse order:
      const Vector dummy = ig[Range(ig.nelem()-1,ig.nelem(),-1)];

      // Copy dummy back to ig vector:
      ig = dummy;
    }
  else
    {
      ostringstream os;
      os << "The Vector *" << ingrid_name <<  "* must be either\n"
         << "strictly increasing or strictly decreasing,\n"
         << "but this is not the case.\n";
      os << "The vector contains:\n"
         << ig;
      throw runtime_error( os.str() );
    }

  // Sort also the vector of points to insert in increasing order:
  {
    ArrayOfIndex si;                   // Sorted indices
    get_sorted_indexes (si, p);        // Get sorted p indices
    const Vector dummy = p;            // Copy p to dummy
    // Copy back dummy to p in right order:
    for (Index j = 0; j < p.nelem(); j++)
      p[j] = dummy[si[j]];
  }

  // The idea is to step through both ig and p, and build up the
  // output in a temporary array.
  Array<Numeric> x;
  Index iig=0, ip=0;            // indices to ig and p
  Index sk=0;                   // skip count
  while ( iig<ig.nelem() && ip<p.nelem() )
    {
      if ( p[ip]<ig[iig] )
        {
          x.push_back(p[ip]);
          ++ip;
        }
      else if ( p[ip]>ig[iig] )
        {
          x.push_back(ig[iig]);
          ++iig;
        }
      else
        {
          out3 << "  Skipping point " << p[ip] << ", which is already "
               << "in the original grid.\n";
          ++ip;
          ++sk;
        }
    }
  
  out2 << "  " << sk << " points skipped.\n";

  // Add remaining points of either p or ig, depending on which is
  // longer:
  if ( ip==p.nelem() )
    {
      // p has reached its end.
      while ( iig<ig.nelem() )
        {
          x.push_back(ig[iig]);
          ++iig;
        }
    }
  else if ( iig==ig.nelem() )
    {
      // ig has reached its end
      while ( ip<p.nelem() )
        {
          x.push_back(p[ip]);
          ++ip;
        }
    }
  else
    {
      // We should never be here.
      assert(false);
      arts_exit();
    }

  // Ok, x should now contain the new grid.

  og.resize(x.nelem());

  // Copy to result vector, turn around if necessary.
  if (ascending)
    for ( Index i=0; i<x.nelem(); ++i )
      og[i] = x[i];               // Just copy.
  else
    for ( Index i=0; i<x.nelem(); ++i )
      og[i] = x[x.nelem()-1-i];   // Copy in reverse order.

}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorLinSpace(      Vector&    x, 
                    const String&    x_name,
                    const Numeric&   start,
                    const Numeric&   stop,
                    const Numeric&   step )
{
  linspace(x,start,stop,step);
  out2 << "  Creating " << x_name << " as linearly spaced vector.\n";
  out3 << "        length : " << x.nelem() << "\n";
  out3 << "   first value : " << x[0] << "\n";
  if ( x.nelem() > 1 )
  {
    out3 << "          step size : " << x[1]-x[0] << "\n";
    out3 << "         last value : " << x[x.nelem()-1] << "\n";
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
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
  out3 << "  Multiplying " << M_name << " with " << x_name
       << " and storing the result in " << y_name << ".\n";
  
  // Check that dimensions are right, x must match columns of M:
  if (M.ncols()!=x.nelem())
    {
      ostringstream os;
      os << "Matrix and vector dimensions must be consistent!\n"
         << M_name << ".ncols() = " << M.ncols() << "\n"
         << x_name << ".nelem() = " << x.nelem();
      throw runtime_error( os.str() );
    }

  // Temporary for the result:
  Vector dummy( M.nrows() );

  mult(dummy,M,x);

  y.resize(dummy.nelem());

  y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
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
  out3 << "            length : " << n << "\n";
  out3 << "       first value : " << x[0] << "\n";
  if ( x.nelem() > 1 )
    {
      out3 << "         step size : " << x[1]-x[0] << "\n";
      out3 << "        last value : " << x[x.nelem()-1] << "\n";
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
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

  nlogspace(x,start,stop,n);
  out2 << "  Creating " << x_name << " as logarithmically spaced vector.\n";
  out3 << "            length : " << n << "\n";
  out3 << "       first value : " << x[0] << "\n";
  if ( x.nelem() > 1 )
    out3 << "        last value : " << x[x.nelem()-1] << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
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


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSet(           Vector&    x, 
                    const String&    x_name,
                    const Index&     n,
                    const Numeric&   value )
{
  x.resize(n);
  x = value;            
  out2 << "  Creating " << x_name << " as a constant vector.\n"; 
  out3 << "            length : " << n << "\n";
  out3 << "             value : " << value << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSetExplicitly( Vector&       x, 
                          const String& x_name,
                          const Vector& values )
{
  x = values;
  
  out2 << "  Creating " << x_name << ".\n"; 
  out3 << "  " << x_name << " = " << x << "\n";
}

