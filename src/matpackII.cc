/* Copyright (C) 2001-2012
   Stefan Buehler <sbuehler@ltu.se>
   Mattias Ekstroem <ekstrom@rss.chalmers.se>

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

/*!
  \file   matpackII.cc
  \author Stefan Buehler <sbuehler@ltu.se>
  \author Mattias Ekstroem
  \author Oliver Lemke
  \author Simon Pfreundschuh <simonpf@chalmers.se> (Eigen interface)
  \date   Wed Jan 06 10:17:15 2016

  \brief  Implementation of sparse matrices.

  Notes:

  There are two different ways to index:
  S.rw(3,4) = 1;                // Read and write
  cout << S.ro(3,4);            // Read only

  This distinction is necessary, because rw() creates elements if they
  don't already exist.

  The normal index operator "()" correspondes to ro, so "S(3,4)" is
  the same as S.ro(3,4).
*/

// #include <vector>
#include <algorithm>
#include <set>
#include <iostream>             // For debugging.
#include <cmath>
#include <iterator>
#include "matpackII.h"
#include <Eigen/Core>

using std::vector;
using std::setw;
using std::cout;
using std::endl;


// Simple member Functions
// ----------------

//! Returns true if variable size is zero.
bool Sparse::empty() const
{
    return ( matrix.rows() == 0 || matrix.cols() == 0);
}

//! Returns the number of rows.
Index Sparse::nrows() const
{
    return matrix.rows();
}

//! Returns the number of columns.
Index Sparse::ncols() const
{
    return matrix.cols();
}

//! Returns the number of nonzero elements.
Index Sparse::nnz() const
{
    return matrix.nonZeros();
}

// Index Operators
// ---------------

//! Read and write index operator.
/*!
  This has to correctly handle two cases:

  1. The data element exists. In this case the operator acts similar
     to the const index operator in that the element is returned.

  2. The element does not exist. In this case it is created.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.
*/
Numeric& Sparse::rw(Index r, Index c)
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<nrows() );
  assert( c<ncols() );

  return matrix.coeffRef( (int) r, (int) c );

}

//! Plain index operator.
/*!
  This is the same as the .ro index operator.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices.

*/
Numeric Sparse::operator() (Index r, Index c) const
{
    return matrix.coeff( (int) r, (int) c );
}

//! Read only index operator.
/*!
  This has to correctly handle two cases:

  1. The data element exists. In this case the element is returned.

  2. The element does not exist. In this case the value 0 is returned.

  \param r Row index.
  \param c Column index.

  \return The data element with these indices, or zero.

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003 */
Numeric Sparse::ro (Index r, Index c) const
{
  // Check if indices are valid:
  assert( 0<=r );
  assert( 0<=c );
  assert( r<nrows() );
  assert( c<ncols() );

  return matrix.coeff( (int) r,(int) c );
}


// Constructors
// ------------

//! Default constructor.
/*!
  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
Sparse::Sparse() :
    matrix()
{
  // Nothing to do here
}

//! Constructor setting size.
/*!
  \param r Row dimension of new sparse matrix.
  \param c Column dimension of new sparse matrix.
*/
Sparse::Sparse(Index r, Index c) :
    matrix( (int) r,(int) c )
{
  // Nothing to do here.
}

//! Convert to dense matrix.
/*!
  Converts a given sparse matrix to a dense matrix. Intended mainly
  for testing purposes.

  \return The dense representation of the given sparse matrix.
*/
Sparse::operator Matrix()
{
    Index m,n;
    m = nrows(); n = ncols();
    Matrix A( m, n );
    A = 0.0;

    for (int i=0; i < m; i++)
    {
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it( matrix, i );
        for (; it; ++it)
        {
            A( it.row(), it.col() ) = it.value();
        }
    }
    return A;
}

void Sparse::list_elements( Vector &values,
                            ArrayOfIndex &row_indices,
                            ArrayOfIndex &column_indices ) const
{
    Index m, n;
    m = nrows(); n = ncols();
    Matrix A( m, n );
    A = 0.0;

    values.resize( nnz() );
    row_indices.resize( nnz() );
    column_indices.resize( nnz() );

    Index j = 0;
    for (int i=0; i < m; i++)
    {
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it( matrix, i );
        for (it; it; ++it)
        {
            values[j] = it.value();
            row_indices[j] = it.row();
            column_indices[j] = it.col();
            j++;
        }
    }
}




//! Insert row function
/*!
  Inserts a Vector as row of elements at the given position.

  The row index must agree with the size of the matrix. This
  function can not be used to expand the matrix.
  Only non-zero values will be stored. If the destination row
  already exist it will be overwritten.

  \param r Where to insert the row
  \param v Vector to be inserted.
*/
void Sparse::insert_row(Index r, Vector v)
{

  // Check if the row index and the Vector length are valid
  assert( 0<=r );
  assert( r<nrows() );
  assert( v.nelem()==ncols() );

  for ( int i = 0; i < ncols(); i++ )
  {
      if (v[i] != 0)
          matrix.coeffRef( (int) r,i) = v[i];
  }

}

//! Make Identity matrix
/*!
  This functions sets the given square matrix to the identity matrix.
*/
void Sparse::identity()
{

    assert( ncols() == nrows() );

    matrix.setIdentity();
}

//! Resize function.
/*!
  If the size is already correct this function does nothing.

  All data is lost after resizing! The new Sparse is not initialised
  so it will be empty.

  \param r New row dimension.
  \param c New column dimension.
*/
void Sparse::resize(Index r, Index c)
{
    assert( 0<=r );
    assert( 0<=c );

    matrix.resize( (int) r, (int) c );
}



//! Output operator for Sparse.
/*!
  \param os Output stream.
  \param v Sparse matrix to print.

  \return Output stream.
*/
std::ostream& operator<<(std::ostream& os, const Sparse& M)
{
    os << M.matrix;
    return os;
}

// General matrix functions

//! Absolute value of sparse matrix elements
/*!
  Computes the absolute values of the elements in sparse matrix B.

  The output matrix A must have been initialized with the correct size.

  \param A Output: Absolute value matrix.
  \param B Original matrix.

  \author Mattias Ekstrom
  \date   2005-03-21
*/
void abs(       Sparse& A,
          const Sparse& B )
{
  // Check dimensions
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == B.ncols() );

  A.matrix = B.matrix.cwiseAbs();
}


//! Sparse matrix - Vector multiplication.
/*!
  This calculates the product

  y = M*x, where M is sparse.

  Output comes first!

  Dimensions of y, M, and x must match. No memory reallocation takes
  place, only the data is copied.

  \param y Output: The multiplication result.
  \param M Matrix for multiplication (sparse).
  \param x Vector for multiplication.
*/
void mult( VectorView y,
           const Sparse& M,
           ConstVectorView x )
{
  // Check dimensions:
  assert( y.nelem() == M.nrows() );
  assert( M.ncols() == x.nelem() );

  // Typedefs for Eigen interface
  typedef Eigen::Matrix< Numeric, Eigen::Dynamic, 1, Eigen::ColMajor>
      EigenColumnVector;
  typedef Eigen::Stride< 1, Eigen::Dynamic > Stride;
  typedef Eigen::Map< EigenColumnVector, 0, Stride > ColumnMap;

  Numeric *data;
  data = x.mdata + x.mrange.get_start();
  ColumnMap x_map( data, x.nelem(), Stride( 1, x.mrange.get_stride() ) );
  data = y.mdata + y.mrange.get_start();
  ColumnMap y_map( data, y.nelem(), Stride( 1, y.mrange.get_stride() ) );

  y_map = M.matrix * x_map;

}

//! SparseMatrix - Matrix multiplication.
/*!
  Calculates the matrix product:

  A = B*C, where B is sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix (full).
  \param B First matrix to multiply (sparse).
  \param C Second matrix to multiply (full).

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
void mult( MatrixView A,
           const Sparse& B,
           const ConstMatrixView& C )
{
  // Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Typedefs for Eigen interface
  typedef Eigen::Matrix< Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      EigenMatrix;
  typedef Eigen::Stride< Eigen::Dynamic, Eigen::Dynamic > Stride;
  typedef Eigen::Map< EigenMatrix, 0, Stride > MatrixMap;

  Index row_stride, column_stride;
  row_stride = C.mrr.get_stride();
  column_stride = C.mcr.get_stride();

  Numeric *data;
  data = C.mdata + C.mrr.get_start() + C.mcr.get_start();
  Stride c_stride( row_stride, column_stride );
  MatrixMap C_map( data, C.nrows(), C.ncols(), c_stride );

  row_stride = A.mrr.get_stride();
  column_stride = A.mcr.get_stride();
  data = A.mdata + A.mrr.get_start() + A.mcr.get_start();
  Stride a_stride( row_stride, column_stride );
  MatrixMap A_map( data, A.nrows(), A.ncols(), a_stride );

  A_map = B.matrix * C_map;

}

//! Matrix - SparseMatrix multiplication.
/*!
  Calculates the matrix product:

  A = B*C, where C is sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix (full).
  \param B First matrix to multiply (sparse).
  \param C Second matrix to multiply (full).

  \author Stefan Buehler <sbuehler@ltu.se>
  \date   Tue Jul 15 15:05:40 2003
*/
void mult( MatrixView A,
           const ConstMatrixView& B,
           const Sparse& C )
{
  // Check dimensions:
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  // Typedefs for Eigen interface.
  typedef Eigen::Matrix< Numeric, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      EigenMatrix;
  typedef Eigen::Stride< Eigen::Dynamic, Eigen::Dynamic > Stride;
  typedef Eigen::Map< EigenMatrix, 0, Stride > MatrixMap;

  Index row_stride, column_stride;
  row_stride = B.mrr.get_stride();
  column_stride = B.mcr.get_stride();

  Numeric *data;
  data = B.mdata + B.mrr.get_start() + B.mcr.get_start();
  Stride b_stride( row_stride, column_stride );
  MatrixMap B_map( data, B.nrows(), B.ncols(), b_stride );

  row_stride = A.mrr.get_stride();
  column_stride = A.mcr.get_stride();
  data = A.mdata + A.mrr.get_start() + A.mcr.get_start();
  Stride a_stride( row_stride, column_stride );
  MatrixMap A_map( data, A.nrows(), A.ncols(), a_stride );

  A_map = B_map * C.matrix;

}

//! Transpose of sparse matrix
/*!
  Computes the transpose of the sparse matrix B.

  The output matrix A must have been initialized with the correct size.

  \param A Output: Transposed matrix.
  \param B Original matrix.

  \author Mattias Ekstroem
  \date   2003-08-05
*/
void transpose( Sparse& A,
                const Sparse& B )
{
  // Check dimensions
  assert( A.nrows() == B.ncols() );
  assert( A.ncols() == B.nrows() );

  A.matrix = B.matrix.transpose();
}


//! Sparse - Sparse multiplication.
/*!
  Calculates A = B*C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of A, B, and C must match. No memory reallocation takes
  place, only the data is copied.

  \param A Output: Result matrix.
  \param B First product matrix.
  \param C Second product matrix.

  \author Mattias Ekstroem
  \date   2003-08-06
*/
void mult( Sparse& A,
           const Sparse& B,
           const Sparse& C )
{
  // Check dimensions and make sure that A is empty
  assert( A.nrows() == B.nrows() );
  assert( A.ncols() == C.ncols() );
  assert( B.ncols() == C.nrows() );

  A.matrix = B.matrix * C.matrix;
}


//! Sparse - Sparse addition.
/*!
  Calculates A = B+C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of B, and C must match.  A will be resized.

  \param A Output: Result matrix.
  \param B First summand matrix.
  \param C Second summand matrix.

  \author Oliver Lemke
  \date   2009-09-03
*/
void add( Sparse& A,
          const Sparse& B,
          const Sparse& C )
{
  // Check dimensions
  assert( B.ncols() == C.ncols() );
  assert( B.nrows() == C.nrows() );

  A.resize( B.nrows(), B.ncols() );
  A.matrix = B.matrix + C.matrix;
}


//! Sparse - Sparse subtraction.
/*!
  Calculates A = B-C, where result A is sparse, and B and C are also sparse.

  Output comes first!

  Dimensions of B, and C must match.  A will be resized.

  \param A Output: Result matrix.
  \param B First subtrahend matrix.
  \param C Second subtrahend matrix.

  \author Oliver Lemke
  \date   2009-09-03
*/
void sub( Sparse& A,
          const Sparse& B,
          const Sparse& C )
{
  A.resize( B.nrows(), B.ncols() );

  // Check dimensions
  assert( B.ncols() == C.ncols() );
  assert( B.nrows() == C.nrows() );

  A.matrix = B.matrix - C.matrix;
}

