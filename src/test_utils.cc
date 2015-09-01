/* Copyright (C) 2001-2012 Stefan Buehler <sbuehler@ltu.se>

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
  \file   test_utils.cc
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Sun May  3 21:56:19 2015
*/

#include "arts.h"
#include <cmath>
#include "lin_alg.h"
#include "test_utils.h"

using std::abs;

//! Fill matrix with random values.
/*!

Fills the given matrix with random values of type Numeric in the range
[0, range], if positive is set to true, or [-range, range], if positive
is set to false.

  \param[in,out] A The matrix to be filled.
  \param[in] range The range of the values to draw the values from.
  \param[in] positive If true the matrix is filled with values from the interval
                  [0,range], otherwise the values are taken from the interval
                  [-range, range].
*/
void random_fill_matrix( MatrixView A,
                         Numeric range,
                         bool positive )
{
    Index m = A.nrows();
    Index n = A.ncols();

    Rand<Numeric> rand( positive ? 0 : - range , range );

    for (Index i = 0; i<m; i++)
    {
        for (Index j = 0; j<n; j++)
        {

            A( i, j ) = (Numeric) rand();

        }
    }
}

void random_fill_matrix_symmetric( MatrixView A,
                                   Numeric range,
                                   bool positive )
{
    random_fill_matrix( A, range, positive);
    Matrix M( A );
    A += transpose( M );
}

//! Fill vector with random values.
/*!

  Fills the given vector with random values of type Numeric drawn
  from the range [0, range], if positive is set to true, or from the
  range [-range, range], if positive == false.

  \param[in,out] v The vector to be filled.
  \param[in] range The range from which the values are taken.
  \param[in] positive If true, the values are taken from the interval [0, range],
                      otherwise from the range [-range, range].
*/
void random_fill_vector( VectorView v,
                         Numeric range,
                         bool positive )
{

    Index n = v.nelem();

    Rand<Numeric> rand( positive ? 0 : - range , range );

    for (Index i = 0; i < n; i++)
    {
        v[i] = rand();
    }
}

//! Pick random random submatrix of size m times n.
/*!
  Randomly chooses a submatrix of the given matrix A and returns the
  corresponding MatrixView.

  \param[in] A The matrix to choose the submatrix from.
  \param[in] m Number of rows of the submatrix.
  \param[in] n Number of columns of the submatrix.

  \return ConstMatrixView corresponding to a randomly chosen m-by-n submatrix.
*/
MatrixView random_submatrix( MatrixView A,
                             int m,
                             int n )
{
    Index m0( A.nrows() ), n0( A.ncols() );
    assert( (0 <= m) && (m <= m0) );
    assert( (0 <= n) && (m <= n0) );

    Rand<Index> rand_m( 0, (m0 - m - 1) ), rand_n( 0, (n0 - n - 1) );
    Index m1, n1;
    m1 = rand_m();
    n1 = rand_n();

    Range r1( m1, m ), r2( n1, n );
    return A( r1, r2 );
}

//! Generate random sub-range of the range [0, n-1].
/*!
  Generate random Range object such that 0 <= extent <= n
  and 0 <= start < n - extent.

  \param n The range [0, n-1] to pick the sub-range from.
  \return Random sub-range.
*/
Range random_range( Index n )
{
    Rand<Index> extent( 1, n );
    Index e = extent();

    Index s = 0;
    if ( 0 <= (n - e - 1) )
    {
        Rand<Index> start( 0, n - e - 1 );
        s = start();
    }
    return Range( s, e );
}

//! Maximum element-wise error of two vectors.
/*!
  If relative is true, the maximum element-wise error is computed.
  Otherwise the absolute error is computed.

  \param[in] v1 The first vector.
  \param[in] v2 The reference vector used to normalize the relative error.
  \param[in] relative If true the relative error is computed, otherwise the absolute
                  error is computed.

  \return The maximum relative or absolute element-wise error.
*/
Numeric max_error( ConstVectorView v1,
                   ConstVectorView v2,
                   bool relative )
{

    Index n = min( v1.nelem(), v2.nelem() );

    Numeric max = 0.0, err = 0.0;

    for ( Index i = 0; i < n; i++ )
    {

        err = 0.0;

        if ( relative )
        {

            if ( v2[i] != 0.0 )
            {
                err = abs( ( v2[i] - v1[i] ) / v2[i] );
            }

        } else {

            err = abs( v2[i] - v2[i] );

        }

        if (err > max)
        {
            max = err;
        }
    }
    return err;
}

//! Maximum element-wise error of two matrices.
/*!
  If relative is true, the maximum element-wise error is computed.
  Otherwise the absolute error is computed.

  \param[in] A1 The first matrix.
  \param[in] A2 The reference matrix used to normalize the relative error.
  \param[in] relative If true the relative error is computed, otherwise the absolute
                  error is computed.

  \return The maximum relative or absolute element-wise error.
*/
Numeric max_error( ConstMatrixView A1,
                   ConstMatrixView A2,
                   bool relative )
{

    Index m = min( A1.nrows(), A2.nrows() );
    Index n = min( A1.ncols(), A2.ncols() );

    Numeric max = 0.0, err = 0.0;

    for ( Index i = 0; i < m; i++ )
    {

        for ( Index j = 0; j < n; j++ )
        {

            err = 0.0;

            if ( relative )
            {

                if ( A2(i,j) != 0.0 )
                {
                    err = abs( ( A2(i,j) - A1(i,j) ) / A2(i,j) );
                }

            } else {
                err = A2(i, j) - A1(i, j);
            }

            if (err > max)
            {
                max = err;
            }
        }
    }

    return max;
}
