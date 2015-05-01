/*!
  \file   oem.cc
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri Apr 17 16:39:25 2015

  \brief Optimal estimation method for retrieval.
*/

#include "arts.h"
#include "lin_alg.h"
#include "oem.h"

//! MAP solution for the inversion of a linear model.
/*!

Solves the inverse of a linear forward model by computing the MAP
solution as described in Rodgers, Inverse Methods for Atmospheric Sounding,
p. 67. In particular, the m-form (eq. (4.6)) is used, i.e. the inversion is
performed for a m times m matrix.
For the execution a n times m matrix, 2 m times matrices and two vectors with
m and n elements respectively are allocated. Otherwise no mermory is reallocates.
The given Matrix and Vector views may not overlap.

  \param[out] x The optimal, estimated state vector consisting of n elements.
  \param[in] y The measurement vector consisting of m elements.
  \param[in] xa The mean a priori state vector
  \param[in] K The weighting function (m,n)-matrix
  \param[in] Se The measurement error covariance (m,m)-matrix
  \param[in] Sa The a priori covariance (n,n)-matrix
*/
void linear_oem( VectorView x,
		 ConstVectorView y,
		 ConstVectorView xa,
		 ConstMatrixView K,
		 ConstMatrixView Se,
		 ConstMatrixView Sa,
                 bool mform )
{

    Index m = K.nrows();
    Index n = K.ncols();

    // Check dimensions for consistency.
    assert( x.nelem() == n );
    assert( xa.nelem() == n );
    assert( y.nelem() == m );

    assert( (Se.ncols() == m) && (Se.nrows() == m) );
    assert( (Sa.ncols() == n) && (Sa.nrows() == n) );



    if ( !mform )
    {

	// n-form (eq. (4.4)).
	Matrix SeInv(m,m);
	Matrix tmp_nm(n,m), tmp_nm2(n,m), tmp_nn(n,n), tmp_nn2(n,n);
	ArrayOfIndex indx(n);
	Vector tmp_n(n), tmp_m(m);

	id_mat(tmp_nn2); // tmp_nn2 = I
	inv( SeInv, Se );

	mult( tmp_nm, transpose(K), SeInv );
	mult( tmp_nm2, Sa, tmp_nm ); // tmp_nm2 = S_a K^T S_e^(-1)
	mult( tmp_nn, tmp_nm2, K);
	tmp_nn2 += tmp_nn;

	mult( tmp_n, tmp_nm2, y);
	tmp_n += xa;

	// Use LU decomposition instead of inversion to save lots of time.
	ludcmp( tmp_nn, indx, tmp_nn2 );
	lubacksub( x, tmp_nn, tmp_n, indx );

    } else {

        // m-form (eq. (4.6)).
	Matrix tmp_nm(n,m);
	Matrix tmp_mm(m,m);
	Matrix Inv(m,m);
	Vector tmp_m(m), tmp_n(n);

	mult( tmp_nm, Sa, transpose(K) ); // tmp_nm = S_a * K^T
	mult( tmp_mm, K, tmp_nm);
	tmp_mm += Se;
	inv( Inv, tmp_mm);

	mult( tmp_m, K, xa);
	tmp_m *= (-1.0);
	tmp_m += y;          // tmp_m = y - K*x_a

	mult( x, Inv, tmp_m);
	mult( tmp_m, tmp_nm, x );

	x = tmp_m;
	x += xa;
    }
}
