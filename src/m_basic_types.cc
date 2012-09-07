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
#include "make_array.h"
#include "math_funcs.h"
#include "messages.h"
#include "logic.h"
#include "sorting.h"


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
void MatrixCompare(const Matrix&    matrix1,
                   const Matrix&    matrix2,
                   const Numeric&   maxabsdiff,
                   const String&    error_message,
                   const Verbosity& verbosity)
{
  const Index nrows = matrix1.nrows();
  const Index ncols = matrix1.ncols();

  if( matrix2.nrows() != nrows  ||  matrix2.ncols() != ncols )
    throw runtime_error( "The two matrices do not have the same size." );

  Numeric maxdiff = 0.0;

  for( Index r=0; r<nrows; r++ )
    { 
      for( Index c=0; c<ncols; c++ )
        {
          const Numeric diff = abs( matrix1(r,c) - matrix2(r,c) );
          if( diff > maxdiff )
            { maxdiff = diff; }
        }
    }

  if( maxdiff > maxabsdiff )
    {
      ostringstream os;
      os << "Checked failed!\n";
      if (error_message.length()) os << error_message << "\n";
      os << "Max allowed deviation set to : " << maxabsdiff << endl
         << "but the matrices deviate with: " << maxdiff << endl;
      throw runtime_error(os.str());
    }

  CREATE_OUT2;
  out2 << "   Check OK (maximum difference = " << maxdiff << ").\n";
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
void NumericCompare(const Numeric&   n1,
                    const Numeric&   n2,
                    const Numeric&   maxabsdiff,
                    const String&    error_message,
                    const Verbosity& verbosity)
{
  const Numeric maxdiff = abs(n1-n2);
  if( maxdiff > maxabsdiff )
  {
    ostringstream os;
    os << "Checked failed!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsdiff << endl
       << "but the value deviates with:  " << maxdiff << endl;
    throw runtime_error(os.str());
  }
  
  CREATE_OUT2;
  out2 << "   Check OK (maximum difference = " << maxdiff << ").\n";
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
void Tensor7Compare(const Tensor7&  aa,
                   const Tensor7&   bb,
                   const Numeric&   maxabsdiff,
                   const String&    error_message,
                   const Verbosity& verbosity)
{
    const Index ncols = aa.ncols();
    const Index nrows = aa.nrows();
    const Index npages = aa.npages();
    const Index nbooks = aa.nbooks();
    const Index nshelves = aa.nshelves();
    const Index nvitrines = aa.nvitrines();
    const Index nlibraries = aa.nlibraries();
    
    if(bb.ncols() != ncols   ||
       bb.nrows() != nrows  ||
       bb.npages() != npages   ||
       bb.nbooks() != nbooks   ||
       bb.nshelves() != nshelves   ||
       bb.nvitrines() != nvitrines   ||
       bb.nlibraries() != nlibraries       )
        throw runtime_error( "The two tensors do not have the same size." );
    
    Numeric maxdiff = 0.0;
    
    for( Index c=0; c<ncols; c++ )
        for( Index r=0; r<nrows; r++ )
            for( Index p=0; p<npages; p++ )
                for( Index b=0; b<nbooks; b++ )
                    for( Index s=0; s<nshelves; s++ )
                        for( Index v=0; v<nvitrines; v++ )
                            for( Index l=0; l<nlibraries; l++ )
            {
              const Numeric diff = abs( aa(l,v,s,b,p,r,c) - bb(l,v,s,b,p,r,c) );
              if( diff > maxdiff )
                { maxdiff = diff; }
            }

    
    if( maxdiff > maxabsdiff )
      {
        ostringstream os;
        os << "Checked failed!\n";
        if (error_message.length()) os << error_message << "\n";
        os << "Max allowed deviation set to : " << maxabsdiff << endl
        << "but the tensors deviate with: " << maxdiff << endl;
        throw runtime_error(os.str());
      }
    
    CREATE_OUT2;
    out2 << "   Check OK (maximum difference = " << maxdiff << ").\n";
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
void VectorCompare(const Vector&    vector1,
                   const Vector&    vector2,
                   const Numeric&   maxabsdiff,
                   const String&    error_message,
                   const Verbosity& verbosity)
{
  const Index n = vector1.nelem();

  if( vector2.nelem() != n )
    throw runtime_error( "The two vectors do not have the same size." );

  Numeric maxdiff = 0.0;
  for( Index i=0; i <n; i++ )
    {
      const Numeric diff = abs( vector1[i] - vector2[i] );
      if( diff > maxdiff )
      { maxdiff = diff; }
    }
    
  if( maxdiff > maxabsdiff )
  {
    ostringstream os;
    os << "Checked failed!\n";
    if (error_message.length()) os << error_message << "\n";
    os << "Max allowed deviation set to: " << maxabsdiff << endl
       << "but the vectors deviate with: " << maxdiff << endl;
    throw runtime_error(os.str());
  }
  
  CREATE_OUT2;
  out2 << "   Check OK (maximum difference = " << maxdiff << ").\n";
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

