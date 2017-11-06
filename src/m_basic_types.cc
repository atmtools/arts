/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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

#include <cmath>
#include "arts.h"
#include "array.h"
#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"
#include "matpackIV.h"
#include "matpackV.h"
#include "matpackVI.h"
#include "matpackVII.h"
#include "mystring.h"
#include "exceptions.h"
#include "math_funcs.h"
#include "messages.h"
#include "lin_alg.h"
#include "logic.h"
#include "sorting.h"
#include "gridded_fields.h"
#include "optproperties.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfIndexSet(ArrayOfIndex& aoi, 
                     const ArrayOfIndex& values,
                     const Verbosity&)
{
  aoi = values;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfIndexSetConstant(ArrayOfIndex& aoi,
                     const Index&   nelem,
                     const Index&   value,
                     const Verbosity&)
{
  aoi.resize(nelem);
  for( Index i=0; i<nelem; i++ )
    aoi[i] = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfIndexLinSpace(ArrayOfIndex&  x,
                          const Index&   start,
                          const Index&   stop,
                          const Index&   step,
                          const Verbosity& verbosity)
{
    CREATE_OUT2;
    CREATE_OUT3;

    Index n = (Index) floor( (stop-start)/step ) + 1;
    if (n < 1) n = 1;

    x.resize(n);

    for ( Index i = 0; i < n; i++ )
        x[i] = start + i*step;

    out2 << "  Creating a linearly spaced ArrayOfIndex.\n";
    out3 << "        length : " << x.nelem() << "\n";
    out3 << "   first value : " << x[0] << "\n";

    if ( x.nelem() > 1 )
    {
        out3 << "     step size : " << x[1]-x[0] << "\n";
        out3 << "    last value : " << x[x.nelem()-1] << "\n";
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfStringSet(ArrayOfString&  sa,
                      const ArrayOfString&  sa2,
                      const Verbosity&)
{
  sa.resize(sa2.nelem());
  sa = sa2;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FlagOff(Index& x, const Verbosity&)
{
  x = 0;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void FlagOn(Index& x, const Verbosity&)
{
  x = 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IndexAdd(Index&   out,
              const Index&   in,
              const Index&   value,
              const Verbosity&)
{
    out = value + in;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IndexSet(Index& x,
              const Index& value,
              const Verbosity&)
{
  x = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IndexStepDown(Index& xout,       
               const Index& xin,
               const Verbosity&)
{
  xout = xin - 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void IndexStepUp(Index& xout,       
               const Index& xin,
               const Verbosity&)
{
  xout = xin + 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixAddScalar(Matrix&   out,
                     const Matrix&   in,
                     const Numeric&  value,
                     const Verbosity&)
{
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
      out.resize( in.nrows(), in.ncols() );
      out = in;
      out += value;
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixCopySparse(Matrix&   out,
                     const Sparse&   in,
                     const Verbosity&)
{
  // There is probably a more efficient way to do this
  out.resize( in.nrows(), in.ncols() );
  for( Index r=0; r<in.nrows(); r++ )
    {
      for( Index c=0; c<in.ncols(); c++ )
        {
          out(r,c) = in(r,c);
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixExtractFromTensor3(
      // WS Generic Output:
      Matrix&          m,
      // WS Input:
      // WS Generic Input:
      const Tensor3&    t3,
      const Index&     index,
      // Control Parameters:
      const String& direction,
      const Verbosity&)
{
  if (direction=="page")
    {
      if( index >= t3.npages() )
        {
          ostringstream os;
          os << "The index " << index 
             << " is outside the page range of the Matrix.";
          throw runtime_error( os.str() );

        }

      m.resize( t3.nrows(), t3.ncols() );
      m = t3( index, joker, joker );
    }
  else if (direction=="row")
    {
      if( index >= t3.nrows() )
        {
          ostringstream os;
          os << "The index " << index 
             << " is outside the row range of the Matrix.";
          throw runtime_error( os.str() );

        }

      m.resize( t3.npages(), t3.ncols() );
      m = t3( joker, index, joker );
    }
  else if (direction=="column")
    {
      if( index >= t3.ncols() )
        {
          ostringstream os;
          os << "The index " << index 
             << " is outside the column range of the Matrix.";
          throw runtime_error( os.str() );

        }

      m.resize( t3.npages(), t3.nrows() );
      m = t3( joker, joker, index );
    }
  else
    {
      ostringstream os;
      os << "Keyword *direction* must be either *page* or *row* or *column*,"
         << "but you gave: " << direction << ".";
      throw runtime_error( os.str() );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixMatrixMultiply(// WS Generic Output:
                          Matrix& Y,
                          // WS Generic Input:
                          const Matrix& M,
                          const Matrix& X,
                          const Verbosity&)
{
  // Check that dimensions are right, M.ncols() must match X.nrows():
  if (M.ncols()!=X.nrows())
    {
      ostringstream os;
      os << "Matrix dimensions must be consistent!\n"
         << "Matrix1.ncols() = " << M.ncols() << "\n"
         << "Matrix2.nrows() = " << X.nrows();
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
void MatrixVectorMultiply(// WS Generic Output:
                          Vector& Y,
                          // WS Generic Input:
                          const Matrix& M,
                          const Vector& X,
                          const Verbosity&)
{
  // Check that dimensions are right, M.ncols() must match X.nrows():
  if (M.ncols()!=X.nelem())
    {
      ostringstream os;
      os << "Matrix and vector dimensions must be consistent!\n"
         << "Matrix.ncols() = " << M.ncols() << "\n"
         << "Vector.nelem() = " << X.nelem();
      throw runtime_error( os.str() );
    }

  // Temporary for the result:
  Vector dummy( M.nrows() );

  mult( dummy, M, X );

  // Copy result to Y:
  Y = dummy;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix1ColFromVector(// WS Generic Output:
                          Matrix&   m,
                          // WS Generic Input:
                          const Vector&   v,
                          const Verbosity&)
{
  const Index nv = v.nelem();

  m.resize( nv, 1 );
  m( joker,0 ) = v;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix2ColFromVectors(// WS Generic Output:
                           Matrix&   m,
                           // WS Generic Input:
                           const Vector&   v1,
                           const Vector&   v2,
                           const Verbosity&)
{
  const Index nv = v1.nelem();

  if( v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  m.resize( nv, 2 );
  m( joker,0 ) = v1;
  m( joker,1 ) = v2;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix3ColFromVectors(// WS Generic Output:
                           Matrix&   m,
                           // WS Generic Input:
                           const Vector&   v1,
                           const Vector&   v2,
                           const Vector&   v3,
                           const Verbosity&)
{
  const Index nv = v1.nelem();

  if( v3.nelem() != nv || v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  m.resize( nv, 3 );
  m( joker,0 ) = v1;
  m( joker,1 ) = v2;
  m( joker,2 ) = v3;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix1RowFromVector(// WS Generic Output:
                          Matrix&   m,
                          // WS Generic Input:
                          const Vector&   v,
                          const Verbosity&)
{
  const Index nv = v.nelem();

  m.resize( 1, nv );
  m( 0, joker ) = v;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix2RowFromVectors(// WS Generic Output:
                           Matrix&   m,
                           // WS Generic Input:
                           const Vector&   v1,
                           const Vector&   v2,
                           const Verbosity&)
{
  const Index nv = v1.nelem();

  if( v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  m.resize( 2, nv );
  m( 0, joker ) = v1;
  m( 1, joker ) = v2;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void Matrix3RowFromVectors(// WS Generic Output:
                           Matrix&   m,
                           // WS Generic Input:
                           const Vector&   v1,
                           const Vector&   v2,
                           const Vector&   v3,
                           const Verbosity&)
{
  const Index nv = v1.nelem();

  if( v3.nelem() != nv || v2.nelem() != nv )
    throw runtime_error("Vectors must be of the same size.");
  
  m.resize( 3, nv );
  m( 0, joker ) = v1;
  m( 1, joker ) = v2;
  m( 2, joker ) = v3;

}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixIdentity(Matrix&   out,
                 const Index&   n,
                 const Numeric&  value,
                 const Verbosity&)
{
  out.resize(n,n);
  id_mat( out );
  if( value != 1 )
    { out *= value; }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixScale(Matrix&   out,
                 const Matrix&   in,
                 const Numeric&  value,
                 const Verbosity&)
{
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
void MatrixSet(Matrix& x, 
               const Matrix& values,
               const Verbosity&)
{
  x = values;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void MatrixSetConstant(Matrix&    x, 
                       const Index&     nrows,
                       const Index&     ncols,
                       const Numeric&   value,
                       const Verbosity&)
{
  x.resize( nrows, ncols );
  x = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericAdd(Numeric&   out,
                const Numeric&   in,
                const Numeric&   value,
                const Verbosity&)
{
  out = value + in;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericInvScale(Numeric&   out,
                     const Numeric&   in,
                     const Numeric&   value,
                     const Verbosity&)
{
  out = in / value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericScale(Numeric&   out,
                  const Numeric&   in,
                  const Numeric&   value,
                  const Verbosity&)
{
  out = value * in;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void NumericSet(Numeric&   x,
                const Numeric&   value,
                const Verbosity&)
{
  x = value;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void SparseSparseMultiply(// WS Generic Output:
                          Sparse& Y,
                          // WS Generic Input:
                          const Sparse& M,
                          const Sparse& X,
                          const Verbosity&)
{
  // Check that dimensions are right, M.ncols() must match X.nrows():
  if (M.ncols()!=X.nrows())
    {
      ostringstream os;
      os << "Matrix dimensions must be consistent!\n"
         << "Matrix1.ncols() = " << M.ncols() << "\n"
         << "Matrix2.nrows() = " << X.nrows();
      throw runtime_error( os.str() );
    }

  // Temporary for the result:
  Sparse dummy( M.nrows(), X.ncols() );

  mult( dummy, M, X );

  // Copy result to Y:
  Y = dummy;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void SparseMatrixIdentity( Sparse        &X,
                           const Index   &n,
                           const Numeric &value,
                           const Verbosity &)
{
  X.resize(n,n);
  id_mat( X );

  if( value != 1.0 )
      X *= value;
}

template<typename MatrixType>
void DiagonalMatrix(MatrixType&           X,
                    const Vector&     diag,
                    const Verbosity& /*v*/);

/* Workspace method: Doxygen documentation will be auto-generated */
template<>
void DiagonalMatrix(Matrix&           X,
                    const Vector&     diag,
                    const Verbosity& /*v*/)
{
    Index n = diag.nelem();
    X.resize(n,n);

    for (Index i = 0; i < n; ++i) {
        X(i, i) = diag[i];
    }
}
template void DiagonalMatrix(Matrix&, const Vector&, const Verbosity&);

/* Workspace method: Doxygen documentation will be auto-generated */
template<>
void DiagonalMatrix(Sparse&           X,
                    const Vector&     diag,
                    const Verbosity& /*v*/)
{
    Index n = diag.nelem();
    X.resize(n,n);

    ArrayOfIndex indices(n);

    for (Index i = 0; i < n; ++i) {
        indices[i] = i;
    }

    X.insert_elements(n, indices, indices, diag);
}
template void DiagonalMatrix(Sparse&, const Vector&, const Verbosity&);

/* Workspace method: Doxygen documentation will be auto-generated */
void StringSet(String&  s, 
               const String&  s2,
               const Verbosity&)
{
  s = s2;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3AddScalar(Tensor3&  out,
                      const Tensor3&  in,
                      const Numeric&  value,
                      const Verbosity&)
{
  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out += value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.npages(), in.nrows(), in.ncols() );
    out = in;
    out += value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor3Scale(Tensor3&  out,
                  const Tensor3&  in,
                  const Numeric&  value,
                  const Verbosity&)
{
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
void Tensor3SetConstant(Tensor3&   x, 
                        const Index&     npages,
                        const Index&     nrows,
                        const Index&     ncols,
                        const Numeric&   value,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  x.resize( npages, nrows, ncols );
  x = value;
  
  out2 << "  Tensor3 = " << value  << "\n";
  out3 << "            npages : " << npages << "\n";
  out3 << "             nrows : " << nrows  << "\n";
  out3 << "             ncols : " << ncols  << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4AddScalar(Tensor4&  out,
                      const Tensor4&  in,
                      const Numeric&  value,
                      const Verbosity&)
{
  // Note that in and out can be the same vector
  if (&out==&in) {
    // Out and in are the same. Just multiply by the scalar value.
    out += value;
  } else {
    // Out and in are different. We first have to copy in to out,
    // then multiply by the scalar value.
    out.resize( in.nbooks(), in.npages(), in.nrows(), in.ncols() );
    out = in;
    out += value;
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor4Scale(Tensor4&  out,
                  const Tensor4&  in,
                  const Numeric&  value,
                  const Verbosity&)
{
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
void Tensor4SetConstant(Tensor4&   x, 
                        const Index&     nbooks,      
                        const Index&     npages,
                        const Index&     nrows,
                        const Index&     ncols,
                        const Numeric&   value,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  x.resize( nbooks, npages, nrows, ncols );
  x = value;
  
  out2 << "  Tensor4 = " << value  << "\n";
  out3 << "            nbooks : " << nbooks << "\n";
  out3 << "            npages : " << npages << "\n";
  out3 << "             nrows : " << nrows  << "\n";
  out3 << "             ncols : " << ncols  << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor5Scale(Tensor5&  out,
                  const Tensor5&  in,
                  const Numeric&  value,
                  const Verbosity&)
{
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
void Tensor5SetConstant(Tensor5&   x, 
                        const Index&     nshelves,     
                        const Index&     nbooks,      
                        const Index&     npages,
                        const Index&     nrows,
                        const Index&     ncols,
                        const Numeric&   value,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  x.resize( nshelves, nbooks, npages, nrows, ncols );
  x = value;
  
  out2 << "  Tensor5 = " << value    << "\n";
  out3 << "          nshelves : " << nshelves << "\n";
  out3 << "            nbooks : " << nbooks   << "\n";
  out3 << "            npages : " << npages   << "\n";
  out3 << "             nrows : " << nrows    << "\n";
  out3 << "             ncols : " << ncols    << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor6Scale(Tensor6&  out,
                  const Tensor6&  in,
                  const Numeric&  value,
                  const Verbosity&)
{
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
void Tensor6SetConstant(Tensor6&   x, 
                        const Index&     nvitrines,   
                        const Index&     nshelves,     
                        const Index&     nbooks,      
                        const Index&     npages,
                        const Index&     nrows,
                        const Index&     ncols,
                        const Numeric&   value,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  x.resize( nvitrines, nshelves, nbooks, npages, nrows, ncols );
  x = value;
  
  out2 << "  Tensor6 = " << value     << "\n";
  out3 << "         nvitrines : " << nvitrines << "\n";    
  out3 << "          nshelves : " << nshelves  << "\n";
  out3 << "            nbooks : " << nbooks    << "\n";
  out3 << "            npages : " << npages    << "\n";
  out3 << "             nrows : " << nrows     << "\n";
  out3 << "             ncols : " << ncols     << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Tensor7Scale(Tensor7&  out,
                  const Tensor7&  in,
                  const Numeric&  value,
                  const Verbosity&)
{
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
void Tensor7SetConstant(Tensor7&   x, 
                        const Index&     nlibraries,          
                        const Index&     nvitrines,   
                        const Index&     nshelves,     
                        const Index&     nbooks,      
                        const Index&     npages,
                        const Index&     nrows,
                        const Index&     ncols,
                        const Numeric&   value,
                        const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  x.resize( nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols );
  x = value;
  
  out2 << "  Tensor7 = " << value      << "\n";
  out3 << "        nlibraries : " << nlibraries << "\n";
  out3 << "         nvitrines : " << nvitrines  << "\n";
  out3 << "          nshelves : " << nshelves   << "\n";
  out3 << "            nbooks : " << nbooks     << "\n";
  out3 << "            npages : " << npages     << "\n";
  out3 << "             nrows : " << nrows      << "\n";
  out3 << "             ncols : " << ncols      << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorAddScalar(Vector&   out,
                     const Vector&   in,
                     const Numeric&  value,
                     const Verbosity&)
{
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
void VectorAddVector(Vector&   c,
                     const Vector&   a,
                     const Vector&   b,
                     const Verbosity&)
{
  // If anything is changed in the method, implement the same change in
  // VectorSubtractVector

  // b has length 1. Here we easily can avoid just adding 0.
  if( b.nelem() == 1 )
    {
      // a and c are the same WSV
      if( &c==&a )
        {
          if( b[0] !=0 )
            { 
              c += b[0]; 
            }
        }
      else
        {
          c = a;
          if( b[0] !=0 )
            { 
              c += b[0];
            }
        }
    }

  // b is a vector
  else if( b.nelem() == a.nelem() )
    {
      // a and c are the same WSV
      if( &c==&a )
        { 
          c += b; 
        }
      else
        {
          c = a;
          c += b;
        }
    }

  else
    throw runtime_error( "The vector *b* must have length 1 or match *a* in length." );
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSubtractVector(Vector&   c,
                     const Vector&   a,
                     const Vector&   b,
                     const Verbosity&)
{
  // If anything is changed in the method, implement the same change in
  // VectorAddVector

  // b has length 1. Here we easily can avoid just adding 0.
  if( b.nelem() == 1 )
    {
      // a and c are the same WSV
      if( &c==&a )
        {
          if( b[0] !=0 )
            { 
              c -= b[0]; 
            }
        }
      else
        {
          c = a;
          if( b[0] !=0 )
            { 
              c -= b[0];
            }
        }
    }

  // b is a vector
  else if( b.nelem() == a.nelem() )
    {
      // a and c are the same WSV
      if( &c==&a )
        { 
          c -= b; 
        }
      else
        {
          c = a;
          c -= b;
        }
    }

  else
    throw runtime_error( "The vector *b* must have length 1 or match *a* in length." );
}



/* Workspace method: Doxygen documentation will be auto-generated */
void VectorCrop(       Vector&   out,
                 const Vector&   in,
                 const Numeric&  min_value,
                 const Numeric&  max_value,
                 const Verbosity&)
{
  const Index nin  = in.nelem();

  Index nout = 0;
  //
  for( Index i=0; i<nin; i++ )
    { 
      if( in[i] >= min_value  &&  in[i] <= max_value )
        { nout += 1; }
    }

  // Make copy if in-vector, as it also can be the out one
  Vector c(in);
  
  out.resize( nout );

  nout = 0;
  //
  for( Index i=0; i<nin; i++ )
    { 
      if( c[i] >= min_value  &&  c[i] <= max_value )
        { 
          out[nout] = c[i];
          nout += 1; 
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated 

   2004-09-15 Patrick Eriksson

   Added keyword to control if row or column is extracted.

   2007-07-24 Stefan Buehler */
void VectorExtractFromMatrix(
      // WS Generic Output:
      Vector&          v,
      // WS Input:
      // WS Generic Input:
      const Matrix&    m,
      const Index&     index,
      // Control Parameters:
      const String& direction,
      const Verbosity&)
{
  if (direction=="row")
    {
      if( index >= m.nrows() )
        {
          ostringstream os;
          os << "The index " << index 
             << " is outside the row range of the Matrix.";
          throw runtime_error( os.str() );

        }

      v.resize( m.ncols() );
      v = m( index, joker );
    }
  else if (direction=="column")
    {
      if( index >= m.ncols() )
        {
          ostringstream os;
          os << "The index " << index 
             << " is outside the column range of the Matrix.";
          throw runtime_error( os.str() );

        }

      v.resize( m.nrows() );
      v = m( joker, index );
    }
  else
    {
      ostringstream os;
      os << "Keyword *direction* must be either *row* or *column*,"
         << "but you gave: " << direction << ".";
      throw runtime_error( os.str() );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorFlip(Vector&   out,
                const Vector&   in,
                const Verbosity&)
{
  const Index n = in.nelem();

  // Note that in and out can be the same vector
  if (&out==&in)
    {
      // Out and in are the same. A copy is needed
      const Vector v = in;
      for( Index i=0; i<n; i++ )
        out[i] = v[n-1-i];
    }
  else
    {
      // Out and in are different. 
      out.resize( n );
      for( Index i=0; i<n; i++ )
        out[i] = in[n-1-i];
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorInsertGridPoints(// WS Generic Output:
                            Vector& og,                  // Output grid
                            // WS Generic Input:
                            const Vector& ingrid,        // Input grid 
                            const Vector& points,        // Points to insert
                            const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  // First make duplicates of the input vectors, in case one of them
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
      os << "The input Vector must be either\n"
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
void VectorLinSpace(Vector&    x, 
                    const Numeric&   start,
                    const Numeric&   stop,
                    const Numeric&   step,
                    const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  linspace(x,start,stop,step);
  
  out2 << "  Creating a linearly spaced vector.\n";
  out3 << "        length : " << x.nelem() << "\n";
  out3 << "   first value : " << x[0] << "\n";
  
  if ( x.nelem() > 1 )
  {
    out3 << "          step size : " << x[1]-x[0] << "\n";
    out3 << "         last value : " << x[x.nelem()-1] << "\n";
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorLogSpace(Vector&    x, 
                    const Numeric&   start,
                    const Numeric&   stop,
                    const Numeric&   step,
                    const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  linspace(x,log(start),log(stop),step);
  transform(x,exp,x);
  
  out2 << "  Creating a logarithmically spaced vector.\n";
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
                          // WS Generic Input:
                          const Matrix& M,
                          const Vector& x,
                          const Verbosity&)
{
  // Check that dimensions are right, x must match columns of M:
  if (M.ncols()!=x.nelem())
    {
      ostringstream os;
      os << "Matrix and vector dimensions must be consistent!\n"
         << "Matrix.ncols() = " << M.ncols() << "\n"
         << "Vector.nelem() = " << x.nelem();
      throw runtime_error( os.str() );
    }

  // Temporary for the result:
  Vector dummy( M.nrows() );

  mult(dummy,M,x);

  y.resize(dummy.nelem());

  y = dummy;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorNLinSpace(Vector&    x, 
                     const Index&     n,
                     const Numeric&   start,
                     const Numeric&   stop,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  if ( n<2 ) 
    throw runtime_error("The number of points must be > 1."); 
  nlinspace(x,start,stop,n);
  
  out2 << "  Creating a linearly spaced vector.\n";
  out3 << "            length : " << n << "\n";
  out3 << "       first value : " << x[0] << "\n";
  
  if ( x.nelem() > 1 )
    {
      out3 << "         step size : " << x[1]-x[0] << "\n";
      out3 << "        last value : " << x[x.nelem()-1] << "\n";
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorNLogSpace(Vector&    x, 
                     const Index&     n,
                     const Numeric&   start,
                     const Numeric&   stop,
                     const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  if ( n<2 )
    throw runtime_error("The number of points must be > 1."); 
  if ( (start<=0) || (stop<=0) )
    throw runtime_error("Only positive numbers are allowed."); 

  nlogspace(x,start,stop,n);
  
  out2 << "  Creating a logarithmically spaced vector.\n";
  out3 << "            length : " << n << "\n";
  out3 << "       first value : " << x[0] << "\n";
  
  if ( x.nelem() > 1 )
    out3 << "        last value : " << x[x.nelem()-1] << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorReshapeMatrix(
      Vector&          v,
      const Matrix&    m,
      const String& direction,
      const Verbosity&)
{
  const Index nrows = m.nrows();
  const Index ncols = m.ncols();

  v.resize( nrows*ncols );

  Index iv = 0;

  if (direction=="column")
    {
      for( Index col=0; col<ncols; col++ )
        {
          for( Index row=0; row<nrows; row++ )
            { 
              v[iv] = m(row,col);
              iv++;
            }
        }
    }
  else if (direction=="row")
    {
      for( Index row=0; row<nrows; row++ )
        {
          for( Index col=0; col<ncols; col++ )
            { 
              v[iv] = m(row,col);
              iv++;
            }
        }
    }
  else
    {
      ostringstream os;
      os << "Keyword *direction* must be either *row* or *column*,"
         << "but you gave: " << direction << ".";
      throw runtime_error( os.str() );
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorScale(Vector&   out,
                 const Vector&   in,
                 const Numeric&  value,
                 const Verbosity&)
{
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
void VectorSetConstant(Vector&    x, 
                       const Index&     n,
                       const Numeric&   value,
                       const Verbosity& verbosity)
{
  CREATE_OUT2;
  CREATE_OUT3;
  
  x.resize(n);
  x = value;            
  
  out2 << "  Creating a constant vector.\n"; 
  out3 << "            length : " << n << "\n";
  out3 << "             value : " << value << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorSet(Vector&       x, 
               const Vector& values,
               const Verbosity&)
{
  x = values;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorVectorMultiply(// WS Generic Output:
                          Vector& y,
                          // WS Generic Input:
                          const Vector& x1,
                          const Vector& x2,
                          const Verbosity&)
{
  // Check that dimensions are right, x1 must match x2:
  if (x1.nelem()!=x2.nelem())
    {
      ostringstream os;
      os << "Both vectors have to have identical dimensions!\n"
         << "Vector1.nelem() = " << x1.nelem() << "\n"
         << "Vector2.nelem() = " << x2.nelem();
      throw runtime_error( os.str() );
    }

  Vector dummy;
  dummy.resize(x1.nelem());

  for ( Index i=0; i<x1.nelem(); i++ )
    dummy[i] = x1[i]*x2[i];

  y = dummy;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Numeric&   var1,
             const Numeric&   var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
  Numeric maxdiff = var1-var2;

  if( isnan(var1)  ||  isnan(var2) )
    {
      if( isnan(var1)  &&  isnan(var2) )
        { maxdiff = 0; }
      else if( isnan(var1) )
        {
          ostringstream os;
          os << "Nan found in " << var1name << ", but there is no "
             << "NaN at same position in " << var2name << ".\nThis "
             << "is not allowed.";
          throw runtime_error(os.str());
        }
      else 
        {
          ostringstream os;
          os << "Nan found in " << var2name << ", but there is no "
             << "NaN at same position in " << var1name << ".\nThis "
             << "is not allowed.";
          throw runtime_error(os.str());
        }
    }      

  if( abs(maxdiff) > maxabsdiff )
  {
    ostringstream os;
    os << var1name << "-" << var2name << " FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsdiff << endl
       << "but the value deviates with:  " << maxdiff << endl;
    throw runtime_error(os.str());
  }
  
  CREATE_OUT2;
  out2 << "   " << var1name << "-" << var2name
       << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Vector&    var1,
             const Vector&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
  const Index n = var1.nelem();

  if( var2.nelem() != n )
  {
    ostringstream os;
    os << var1name << " (" << n << ") and "
       << var2name << " (" << var2.nelem() << ") do not have the same size.";
    throw runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;
  for( Index i=0; i <n; i++ )
    {
      Numeric diff = var1[i] - var2[i];

      if( isnan(var1[i])  ||  isnan(var2[i]) )
        {
          if( isnan(var1[i])  &&  isnan(var2[i]) )
            { diff = 0; }
          else if( isnan(var1[i]) )
            {
              ostringstream os;
              os << "Nan found in " << var1name << ", but there is no "
                 << "NaN at same position in " << var2name << ".\nThis "
                 << "is not allowed.";
              throw runtime_error(os.str());
            }
          else 
            {
              ostringstream os;
              os << "Nan found in " << var2name << ", but there is no "
                 << "NaN at same position in " << var1name << ".\nThis "
                 << "is not allowed.";
              throw runtime_error(os.str());
            }
        }      

      if( abs(diff) > abs(maxdiff) )
      { maxdiff = diff; }
    }
    
  if( isnan(maxdiff)  ||  abs(maxdiff) > maxabsdiff )
  {
    ostringstream os;
    os << var1name << "-" << var2name << " FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsdiff << endl
       << "but the vectors deviate with: " << maxdiff << endl;
    throw runtime_error(os.str());
  }
  
  CREATE_OUT2;
  out2 << "   " << var1name << "-" << var2name
       << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Matrix&    var1,
             const Matrix&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
  const Index nrows = var1.nrows();
  const Index ncols = var1.ncols();

  if( var2.nrows() != nrows  ||  var2.ncols() != ncols )
  {
    ostringstream os;
    os << var1name << " (" << nrows << "," << ncols << ") and "
       << var2name << " (" << var2.nrows() << "," << var2.ncols()
       << ") do not have the same size.";
    throw runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for( Index r=0; r<nrows; r++ )
    { 
      for( Index c=0; c<ncols; c++ )
        {
          Numeric diff = var1(r,c) - var2(r,c);

          if( isnan(var1(r,c))  ||  isnan(var2(r,c)) )
            {
              if( isnan(var1(r,c))  &&  isnan(var2(r,c)) )
                { diff = 0; }
              else if( isnan(var1(r,c)) )
                {
                  ostringstream os;
                  os << "Nan found in " << var1name << ", but there is no "
                     << "NaN at same position in " << var2name << ".\nThis "
                     << "is not allowed.";
                  throw runtime_error(os.str());
                }
              else 
                {
                  ostringstream os;
                  os << "Nan found in " << var2name << ", but there is no "
                     << "NaN at same position in " << var1name << ".\nThis "
                     << "is not allowed.";
                  throw runtime_error(os.str());
                }
            }      

          if( abs(diff) > abs(maxdiff) )
            { maxdiff = diff; }
        }
    }

  if( abs(maxdiff) > maxabsdiff )
    {
      ostringstream os;
      os << var1name << "-" << var2name << " FAILED!\n";
      if (error_message.length()) os << error_message << "\n";
      os << "Max allowed deviation set to : " << maxabsdiff << endl
         << "but the matrices deviate with: " << maxdiff << endl;
      throw runtime_error(os.str());
    }

  CREATE_OUT2;
  out2 << "   " << var1name << "-" << var2name
       << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor3&   var1,
             const Tensor3&   var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    const Index ncols = var1.ncols();
    const Index nrows = var1.nrows();
    const Index npages = var1.npages();
    
    if(var2.ncols() != ncols   ||
       var2.nrows() != nrows  ||
       var2.npages() != npages )
    {
      ostringstream os;
      os << var1name << " and " << var2name << " do not have the same size.";
      throw runtime_error(os.str());
    }
    
    Numeric maxdiff = 0.0;
    
    for( Index c=0; c<ncols; c++ )
        for( Index r=0; r<nrows; r++ )
            for( Index p=0; p<npages; p++ )
              {
                Numeric diff = var1(p,r,c) - var2(p,r,c);

                if( isnan(var1(p,r,c))  ||  isnan(var2(p,r,c)) )
                  {
                    if( isnan(var1(p,r,c))  &&  isnan(var2(p,r,c)) )
                      { diff = 0; }
                    else if( isnan(var1(p,r,c)) )
                      {
                        ostringstream os;
                        os << "Nan found in " << var1name << ", but there is no "
                           << "NaN at same position in " << var2name << ".\nThis "
                           << "is not allowed.";
                        throw runtime_error(os.str());
                      }
                    else 
                      {
                        ostringstream os;
                        os << "Nan found in " << var2name << ", but there is no "
                           << "NaN at same position in " << var1name << ".\nThis "
                           << "is not allowed.";
                        throw runtime_error(os.str());
                      }
                  }      

                if( abs(diff) > abs(maxdiff) )
                  { maxdiff = diff; }
              }
    
    if( abs(maxdiff) > maxabsdiff )
      {
        ostringstream os;
        os << var1name << "-" << var2name << " FAILED!\n";
        if (error_message.length()) os << error_message << "\n";
        os << "Max allowed deviation set to : " << maxabsdiff << endl
           << "but the tensors deviate with: " << maxdiff << endl;
        throw runtime_error(os.str());
      }
    
    CREATE_OUT2;
    out2 << "   " << var1name << "-" << var2name
         << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor4&   var1,
             const Tensor4&   var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    const Index ncols = var1.ncols();
    const Index nrows = var1.nrows();
    const Index npages = var1.npages();
    const Index nbooks = var1.nbooks();
    
    if(var2.ncols() != ncols   ||
       var2.nrows() != nrows  ||
       var2.npages() != npages ||
       var2.nbooks() != nbooks )
    {
      ostringstream os;
      os << var1name << " and " << var2name << " do not have the same size.";
      throw runtime_error(os.str());
    }
    
    Numeric maxdiff = 0.0;
    
    for( Index c=0; c<ncols; c++ )
        for( Index r=0; r<nrows; r++ )
            for( Index p=0; p<npages; p++ )
                for( Index b=0; b<nbooks; b++ )
              {
                Numeric diff = var1(b,p,r,c) - var2(b,p,r,c);

                if( isnan(var1(b,p,r,c))  ||  isnan(var2(b,p,r,c)) )
                  {
                    if( isnan(var1(b,p,r,c))  &&  isnan(var2(b,p,r,c)) )
                      { diff = 0; }
                    else if( isnan(var1(b,p,r,c)) )
                      {
                        ostringstream os;
                        os << "Nan found in " << var1name << ", but there is no "
                           << "NaN at same position in " << var2name << ".\nThis "
                           << "is not allowed.";
                        throw runtime_error(os.str());
                      }
                    else 
                      {
                        ostringstream os;
                        os << "Nan found in " << var2name << ", but there is no "
                           << "NaN at same position in " << var1name << ".\nThis "
                           << "is not allowed.";
                        throw runtime_error(os.str());
                      }
                  }      

                if( abs(diff) > abs(maxdiff) )
                  { maxdiff = diff; }
              }
    
    if( abs(maxdiff) > maxabsdiff )
      {
        ostringstream os;
        os << var1name << "-" << var2name << " FAILED!\n";
        if (error_message.length()) os << error_message << "\n";
        os << "Max allowed deviation set to : " << maxabsdiff << endl
           << "but the tensors deviate with: " << maxdiff << endl;
        throw runtime_error(os.str());
      }
    
    CREATE_OUT2;
    out2 << "   " << var1name << "-" << var2name
         << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor5&   var1,
             const Tensor5&   var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    const Index ncols = var1.ncols();
    const Index nrows = var1.nrows();
    const Index npages = var1.npages();
    const Index nbooks = var1.nbooks();
    const Index nshelves = var1.nshelves();

    if(var2.ncols() != ncols   ||
       var2.nrows() != nrows  ||
       var2.npages() != npages   ||
       var2.nbooks() != nbooks   ||
       var2.nshelves() != nshelves )
    {
        ostringstream os;
        os << var1name << " and " << var2name << " do not have the same size.";
        throw runtime_error(os.str());
    }

    Numeric maxdiff = 0.0;

    for( Index c=0; c<ncols; c++ )
        for( Index r=0; r<nrows; r++ )
            for( Index p=0; p<npages; p++ )
                for( Index b=0; b<nbooks; b++ )
                    for( Index s=0; s<nshelves; s++ )
                    {
                        Numeric diff = var1(s,b,p,r,c) - var2(s,b,p,r,c);

                        if( isnan(var1(s,b,p,r,c))  ||  isnan(var2(s,b,p,r,c)) )
                        {
                            if( isnan(var1(s,b,p,r,c))  &&  isnan(var2(s,b,p,r,c)) )
                            { diff = 0; }
                            else if( isnan(var1(s,b,p,r,c)) )
                            {
                                ostringstream os;
                                os << "Nan found in " << var1name << ", but there is no "
                                << "NaN at same position in " << var2name << ".\nThis "
                                << "is not allowed.";
                                throw runtime_error(os.str());
                            }
                            else
                            {
                                ostringstream os;
                                os << "Nan found in " << var2name << ", but there is no "
                                << "NaN at same position in " << var1name << ".\nThis "
                                << "is not allowed.";
                                throw runtime_error(os.str());
                            }
                        }

                        if( abs(diff) > abs(maxdiff) )
                        { maxdiff = diff; }
                    }

    if( abs(maxdiff) > maxabsdiff )
    {
        ostringstream os;
        os << var1name << "-" << var2name << " FAILED!\n";
        if (error_message.length()) os << error_message << "\n";
        os << "Max allowed deviation set to : " << maxabsdiff << endl
        << "but the tensors deviate with: " << maxdiff << endl;
        throw runtime_error(os.str());
    }
    
    CREATE_OUT2;
    out2 << "   " << var1name << "-" << var2name
    << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Tensor7&   var1,
             const Tensor7&   var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    const Index ncols = var1.ncols();
    const Index nrows = var1.nrows();
    const Index npages = var1.npages();
    const Index nbooks = var1.nbooks();
    const Index nshelves = var1.nshelves();
    const Index nvitrines = var1.nvitrines();
    const Index nlibraries = var1.nlibraries();
    
    if(var2.ncols() != ncols   ||
       var2.nrows() != nrows  ||
       var2.npages() != npages   ||
       var2.nbooks() != nbooks   ||
       var2.nshelves() != nshelves   ||
       var2.nvitrines() != nvitrines   ||
       var2.nlibraries() != nlibraries       )
    {
      ostringstream os;
      os << var1name << " and " << var2name << " do not have the same size.";
      throw runtime_error(os.str());
    }
    
    Numeric maxdiff = 0.0;
    
    for( Index c=0; c<ncols; c++ )
        for( Index r=0; r<nrows; r++ )
            for( Index p=0; p<npages; p++ )
                for( Index b=0; b<nbooks; b++ )
                    for( Index s=0; s<nshelves; s++ )
                        for( Index v=0; v<nvitrines; v++ )
                            for( Index l=0; l<nlibraries; l++ )
            {
              Numeric diff = var1(l,v,s,b,p,r,c) - var2(l,v,s,b,p,r,c);

              if( isnan(var1(l,v,s,b,p,r,c))  ||  isnan(var2(l,v,s,b,p,r,c)) )
                {
                  if( isnan(var1(l,v,s,b,p,r,c))  &&  isnan(var2(l,v,s,b,p,r,c)) )
                    { diff = 0; }
                  else if( isnan(var1(l,v,s,b,p,r,c)) )
                    {
                      ostringstream os;
                      os << "Nan found in " << var1name << ", but there is no "
                         << "NaN at same position in " << var2name << ".\nThis "
                         << "is not allowed.";
                      throw runtime_error(os.str());
                    }
                  else 
                    {
                      ostringstream os;
                      os << "Nan found in " << var2name << ", but there is no "
                         << "NaN at same position in " << var1name << ".\nThis "
                         << "is not allowed.";
                      throw runtime_error(os.str());
                    }
                }      

              if( abs(diff) > abs(maxdiff) )
                { maxdiff = diff; }
            }

    if( abs(maxdiff) > maxabsdiff )
      {
        ostringstream os;
        os << var1name << "-" << var2name << " FAILED!\n";
        if (error_message.length()) os << error_message << "\n";
        os << "Max allowed deviation set to : " << maxabsdiff << endl
           << "but the tensors deviate with: " << maxdiff << endl;
        throw runtime_error(os.str());
      }
    
    CREATE_OUT2;
    out2 << "   " << var1name << "-" << var2name
         << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const ArrayOfVector&    var1,
             const ArrayOfVector&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    if( var1.nelem() != var2.nelem() )
    {
        ostringstream os;
        os << "The two arrays do not have the same size." << endl
           << var1name << " nelem: " << var1.nelem() << endl
           << var2name << " nelem: " << var2.nelem() << endl;
        throw runtime_error(os.str());
    }

    bool failed = false;
    ostringstream fail_msg;
    for (Index i = 0; i < var1.nelem(); i++)
    {
        try
        {
            ostringstream vn1, vn2;
            vn1 << var1name << "[" << i << "]";
            vn2 << var2name << "[" << i << "]";
            Compare(var1[i], var2[i], maxabsdiff, error_message,
                    vn1.str(), vn2.str(), "", "", verbosity);
        }
        catch (runtime_error e)
        {
            failed = true;
            fail_msg << endl << e.what() << endl
                     << "Mismatch at array index: " << i << endl;
        }
    }

    if (failed)
        throw runtime_error(fail_msg.str());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const ArrayOfMatrix&    var1,
             const ArrayOfMatrix&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
  if( var1.nelem() != var2.nelem() )
  {
    ostringstream os;
    os << "The two arrays do not have the same size." << endl
    << var1name << " nelem: " << var1.nelem() << endl
    << var2name << " nelem: " << var2.nelem() << endl;
    throw runtime_error(os.str());
  }
  
  bool failed = false;
  ostringstream fail_msg;
  for (Index i = 0; i < var1.nelem(); i++)
  {
    try
    {
      ostringstream vn1, vn2;
      vn1 << var1name << "[" << i << "]";
      vn2 << var2name << "[" << i << "]";
      Compare(var1[i], var2[i], maxabsdiff, error_message,
              vn1.str(), vn2.str(), "", "", verbosity);
    }
    catch (runtime_error e)
    {
      failed = true;
      fail_msg << endl << e.what() << endl
      << "Mismatch at array index: " << i << endl;
    }
  }
  
  if (failed)
    throw runtime_error(fail_msg.str());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const GriddedField3&    var1,
             const GriddedField3&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    for (Index i = 0; i < var1.get_dim(); i++)
    {
        if( var1.get_grid_size(i) != var2.get_grid_size(i) )
        {
            ostringstream os;
            os << var1name << " and "
            << var2name << " grid " << i << " do not have the same size: "
            << var1.get_grid_size(i) << " != " << var2.get_grid_size(i);
            throw runtime_error(os.str());
        }
        if( var1.get_grid_name(i) != var2.get_grid_name(i) )
        {
            ostringstream os;
            os << var1name << " and "
            << var2name << " grid " << i << " do not have the same name: "
            << var1.get_grid_name(i) << " != " << var2.get_grid_name(i);
            throw runtime_error(os.str());
        }
    }

    Compare(var1.data, var2.data, maxabsdiff, error_message,
            var1name, var2name, "", "", verbosity);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const Sparse&    var1,
             const Sparse&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
  const Index nrows = var1.nrows();
  const Index ncols = var1.ncols();

  if( var2.nrows() != nrows  ||  var2.ncols() != ncols )
  {
    ostringstream os;
    os << var1name << " (" << nrows << "," << ncols << ") and "
       << var2name << " (" << var2.nrows() << "," << var2.ncols()
       << ") do not have the same size.";
    throw runtime_error(os.str());
  }

  Numeric maxdiff = 0.0;

  for( Index r=0; r<nrows; r++ )
    { 
      for( Index c=0; c<ncols; c++ )
        {
          Numeric diff = var1(r,c) - var2(r,c);

          if( isnan(var1(r,c))  ||  isnan(var2(r,c)) )
            {
              if( isnan(var1(r,c))  &&  isnan(var2(r,c)) )
                { diff = 0; }
              else if( isnan(var1(r,c)) )
                {
                  ostringstream os;
                  os << "Nan found in " << var1name << ", but there is no "
                     << "NaN at same position in " << var2name << ".\nThis "
                     << "is not allowed.";
                  throw runtime_error(os.str());
                }
              else 
                {
                  ostringstream os;
                  os << "Nan found in " << var2name << ", but there is no "
                     << "NaN at same position in " << var1name << ".\nThis "
                     << "is not allowed.";
                  throw runtime_error(os.str());
                }
            }      

          if( abs(diff) > abs(maxdiff) )
            { maxdiff = diff; }
        }
    }

  if( abs(maxdiff) > maxabsdiff )
    {
      ostringstream os;
      os << var1name << "-" << var2name << " FAILED!\n";
      if (error_message.length()) os << error_message << "\n";
      os << "Max allowed deviation set to : " << maxabsdiff << endl
         << "but the matrices deviate with: " << maxdiff << endl;
      throw runtime_error(os.str());
    }

  CREATE_OUT2;
  out2 << "   " << var1name << "-" << var2name
       << " OK (maximum difference = " << maxdiff << ").\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void Compare(const SingleScatteringData&    var1,
             const SingleScatteringData&    var2,
             const Numeric&   maxabsdiff,
             const String&    error_message,
             const String&    var1name,
             const String&    var2name,
             const String&,
             const String&,
             const Verbosity& verbosity)
{
    if (var1.ptype != var2.ptype)
    {
        std::ostringstream os;
        os << "The particle types don't match: " << std::endl
        << var1name << " = " << PTypeToString(var1.ptype) << ", "
        << var2name << " = " << PTypeToString(var2.ptype) << std::endl;
        throw std::runtime_error(os.str());
    }
    Compare(var1.f_grid, var2.f_grid, maxabsdiff, error_message,
            var1name+".f_grid", var2name+".f_grid", "", "", verbosity);
    Compare(var1.T_grid, var2.T_grid, maxabsdiff, error_message,
            var1name+".T_grid", var2name+".T_grid", "", "", verbosity);
    Compare(var1.za_grid, var2.za_grid, maxabsdiff, error_message,
            var1name+".za_grid", var2name+".za_grid", "", "", verbosity);
    Compare(var1.aa_grid, var2.aa_grid, maxabsdiff, error_message,
            var1name+".aa_grid", var2name+".aa_grid", "", "", verbosity);
    Compare(var1.pha_mat_data, var2.pha_mat_data, maxabsdiff, error_message,
            var1name+".pha_mat_data", var2name+".pha_mat_data", "", "", verbosity);
    Compare(var1.ext_mat_data, var2.ext_mat_data, maxabsdiff, error_message,
            var1name+".ext_mat_data", var2name+".ext_mat_data", "", "", verbosity);
    Compare(var1.abs_vec_data, var2.abs_vec_data, maxabsdiff, error_message,
            var1name+".abs_vec_data", var2name+".abs_vec_data", "", "", verbosity);
}



/* Workspace method: Doxygen documentation will be auto-generated */
void CompareRelative(const Vector&    var1,
                     const Vector&    var2,
                     const Numeric&   maxabsreldiff,
                     const String&    error_message,
                     const String&    var1name,
                     const String&    var2name,
                     const String&,
                     const String&,
                     const Verbosity& verbosity)
{
  const Index n = var1.nelem();

  if( var2.nelem() != n )
  {
    ostringstream os;
    os << var1name << " (" << n << ") and "
       << var2name << " (" << var2.nelem() << ") do not have the same size.";
    throw runtime_error(os.str());
  }

  Numeric maxreldiff = 0.0;
  for(Index i = 0; i <n; i++)
  {
    Numeric diff;
    
    if(isnan(var1[i])  ||  isnan(var2[i]))
    {
      if(isnan(var1[i])  &&  isnan(var2[i]))
        diff = 0.0;
      else if( isnan(var1[i]) )
      {
        ostringstream os;
        os << "Nan found in " << var1name << ", but there is no "
            << "NaN at same position in " << var2name << ".\nThis "
            << "is not allowed.";
        throw runtime_error(os.str());
      }
      else 
      {
        ostringstream os;
        os << "Nan found in " << var2name << ", but there is no "
            << "NaN at same position in " << var1name << ".\nThis "
            << "is not allowed.";
        throw runtime_error(os.str());
      }
    }
    else if(var1[i] == 0.0)
      if(var2[i] == 0.0)
        diff = 0.0;
      else
        diff = 1.0;
    else 
      diff = (var1[i] - var2[i])/var1[i];

    if( abs(diff) > abs(maxreldiff) )
      maxreldiff = diff;
  }
    
  if(isnan(maxreldiff)  ||  abs(maxreldiff) > maxabsreldiff)
  {
    ostringstream os;
    os << var1name << "-" << var2name << " FAILED!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsreldiff*100.0 << "%" << endl
       << "but the vectors deviate with: " << maxreldiff*100.0 << "%" << endl;
    throw runtime_error(os.str());
  }
  
  CREATE_OUT2;
  out2 << "   " << var1name << "-" << var2name
       << " OK (maximum difference = " << maxreldiff*100.0 << "%" << ").\n";
}
