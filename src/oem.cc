/*!
  \file   oem.cc
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri Apr 17 16:39:25 2015

  \brief Optimal inversion methods for linear and non-linear inverse problems.
*/

#include "arts.h"
#include "lin_alg.h"
#include "oem.h"
#include "stdlib.h"
#include "math.h"

//! Linear OEM.
/*!

Computes the inverse of a linear forward model by computing the MAP
solution as described in Rodgers, Inverse Methods for Atmospheric Sounding,
p. 67. In particular, the m-form (eq. (4.6)) is used, i.e. the inversion is
performed for a m times m matrix.

For the execution a n times m matrix, 2 m times matrices and two vectors with
m and n elements respectively are allocated. The given Matrix and Vector views
may not overlap.

  \param[out] x The optimal, estimated state vector consisting of n elements.
  \param[in] y The measurement vector consisting of m elements.
  \param[in] xa The mean a priori state vector
  \param[in] K The weighting function (m,n)-matrix
  \param[in] Se The measurement error covariance (m,m)-matrix
  \param[in] Sa The a priori covariance (n,n)-matrix
*/
void oem_linear( VectorView x,
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
        Vector tmp_m(m), tmp_n(n);

        mult( tmp_nm, Sa, transpose(K) ); // tmp_nm = S_a * K^T
        mult( tmp_mm, K, tmp_nm);
        tmp_mm += Se;

        mult( tmp_m, K, xa);
        tmp_m *= (-1.0);
        tmp_m += y;          // tmp_m = y - K*x_a

        solve( x, tmp_mm, tmp_m);
        mult( tmp_m, tmp_nm, x );

        x = tmp_m;
        x += xa;
    }
}

//! Non-linear OEM using Gauss-Newton method.
/*!
  Computes the optimal nonlinear inverse of a given forward model using the
  Gauss-Newton method as given in eq. (5.8) in Rodgers book. The forward model
  is given by an object implementing the interface described by the ForwardModel
  class. Convergence is determined using equation (5.29). The given
  tolerance is scaled by the problem size n. If the method doesn't converge it
  abords after the given maximum number of iterations.

  During the execution, space for up to 6 additional matrices and vectors is
  allocated. The given Matrix and Vector views should not overlap.

  \param[out] x The optimal inverse x.
  \param[in] y The measurement vector containing m measurements.
  \param[in] xa The size-n a-priori state vector.
  \param[in] K The ForwardModel representing the forward model to invert.
  \param[in] Se The measurement error covariance (m,m)-matrix
  \param[in] Sa The a priori covariance (n,n)-matrix
  \param[in] tol The convergence criterium before scaling by the problem size.
  \param maxIter Tha maximum number of iterations in case of no convergence.

  \return True if the method has converged, false otherwise.
*/
bool oem_gauss_newton( VectorView x,
                       ConstVectorView y,
                       ConstVectorView xa,
                       ForwardModel &K,
                       ConstMatrixView Se,
                       ConstMatrixView Sa,
                       Numeric tol,
                       Index maxIter )
{

    Index n(x.nelem()), m(y.nelem());
    Numeric di2;
    Matrix Ki(m,n), KiTSeInv(n,m), KiTSeInvKi(n,n), SeInv(m,m), SaInv(n,n);
    Vector xi(n), yi(m), tm(m), tn1(n), tn2(n), dx(n);

    Index iter = 0;
    bool converged = false;

    // Compute matrix inverses.
    inv(SeInv, Se);
    inv(SaInv, Sa);

    // Set the starting vector.
    x = xa;

    while ( (!converged) && (iter < maxIter) )
    {
        // Compute Jacobian and y_i.
        K.evaluate_jacobian( yi, Ki, x);

        mult( KiTSeInv, transpose(Ki), SeInv );
        mult( KiTSeInvKi, KiTSeInv, Ki );
        KiTSeInvKi += SaInv;

        // tm = K_i(x_i - x_a)

        tn1 = x;
        tn1 -= xa;
        mult( tn2, SaInv, tn1 );

        tm = y;
        tm -= yi;
        // tn1 = K_i^T * S_e^{-1} * (y - F(x_i))
        mult( tn1, KiTSeInv, tm );

        // This vector is used to test for convergence later.
        // See eqn. (5.31).
        tn1 -= tn2;

        solve( dx, KiTSeInvKi, tn1 );
        x += dx;

        di2 = dx * tn1;
        if ( di2 <= tol * (Numeric) n )
        {
            converged = true;
        }


        // Increase convergence counter and check for convergence.
        iter++;

    }

    return converged;

}

//! Non-linear OEM using Gauss-Newton method.
/*!
  Computes the optimal nonlinear inverse of a given forward model using the
  Gauss-Newton method as given in eq. (5.10) in Rodgers book. The forward model
  is given by an object implementing the interface described by the ForwardModel
  class. Convergence is determined using equation (5.33). The given
  tolerance is scaled by the problem size m. If the method doesn't converge it
  abords after the given maximum number of iterations.

  During the execution, space for up to 6 additional matrices and vectors is
  allocated. The given Matrix and Vector views should not overlap.

  \param[out] x The optimal inverse x.
  \param[in] y The measurement vector containing m measurements.
  \param[in] xa The size-n a-priori state vector.
  \param[in] K The ForwardModel representing the forward model to invert.
  \param[in] Se The measurement error covariance (m,m)-matrix
  \param[in] Sa The a priori covariance (n,n)-matrix
  \param[in] tol The convergence criterium before scaling by the problem size.
  \param maxIter Tha maximum number of iterations in case of no convergence.
*/
bool oem_gauss_newton_m_form( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel &K,
                              ConstMatrixView Se,
                              ConstMatrixView Sa,
                              Numeric tol,
                              Index maxIter )
{
    Index n( x.nelem() ), m( y.nelem() );
    Numeric di2;
    Matrix Ki(m,n), SaKiT(m, n), KiSaKiT(m, m);
    Vector xi(n), yi(m), tm(m), tm2(m), tn(n), yi_old(m);

    Index iter = 0;
    bool converged = false;

    // Set the starting vector.
    x = xa;

    while ((!converged) && (iter < maxIter))
    {

        // Compute Jacobian and y_i.
        K.evaluate_jacobian( yi, Ki, x );
        // If not in the first iteration, check for
        // convergence here.
        if (iter > 0)
        {
            yi_old -= yi;

            solve( tm, Se, yi_old );
            mult( tm2, KiSaKiT, tm );
            // TODO: Optimize using LU decomp.
            solve( tm, Se, tm2);
            di2 = yi_old * tm;
            if ( fabs( di2 ) < tol * (Numeric) n )
            {
                converged = true;
                break;
            }
        }

        mult( SaKiT, Sa, transpose(Ki) );
        mult( KiSaKiT, Ki, SaKiT );
        KiSaKiT += Se;

        // tm = K_i(x_i - x_a)
        tn = x;
        tn -= xa;
        mult( tm, Ki, tn );

        tm -= yi;
        tm += y;

        solve( tm2, KiSaKiT, tm );
        mult( x, SaKiT, tm2 );
        x += xa;

        // Increase iteration counter and store yi.
        iter++;
        yi_old = yi;
    }
    return converged;
}

//! Non-linear OEM using Gauss-Newton method.
/*!
  Computes the optimal nonlinear inverse of a given forward model using the
  Gauss-Newton method as given in eq. (5.9) in Rodgers book. The forward model
  is given by an object implementing the interface described by the ForwardModel
  class. Convergence is determined using equation (5.29). The given
  tolerance is scaled by the problem size m. If the method doesn't converge it
  abords after the given maximum number of iterations.

  During the execution, space for up to 6 additional matrices and vectors is
  allocated. The given Matrix and Vector views should not overlap.

  \param[out] x The optimal inverse x.
  \param[in] y The measurement vector containing m measurements.
  \param[in] xa The size-n a-priori state vector.
  \param[in] K The ForwardModel representing the forward model to invert.
  \param[in] Se The measurement error covariance (m,m)-matrix
  \param[in] Sa The a priori covariance (n,n)-matrix
  \param[in] tol The convergence criterium before scaling by the problem size.
  \param maxIter Tha maximum number of iterations in case of no convergence.

  \return True if the method has converged, false otherwise.
*/
bool oem_gauss_newton_n_form( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel &K,
                              ConstMatrixView Se,
                              ConstMatrixView Sa,
                              Numeric tol,
                              Index maxIter )
{
    Index n( x.nelem() ), m( y.nelem() );
    Numeric di2;
    Matrix Ki(m,n), KiTSeInvKi(n, n), KiTSeInv(m,n), SeInv(m,m), SaInv(n,n);
    Vector xi(n), dx(n), dx_old(n), yi(m), tm(m), tn(n), sInvDx(n);

    Index iter = 0;
    bool converged = false;

    // Required matrix inverses.
    inv(SeInv, Se);
    inv(SaInv, Sa);

    // Set the starting vector.
    x = xa;
    dx_old = 0;

    while ((!converged) && (iter < maxIter))
    {

        // Compute Jacobian and y_i.
        K.evaluate_jacobian( yi, Ki, x );

        mult( KiTSeInv, transpose( Ki ), SeInv );
        mult( KiTSeInvKi, KiTSeInv, Ki );
        KiTSeInvKi += SaInv;

        // tm = K_i(x_i - x_a)
        tn = x;
        tn -= xa;
        mult( tm, Ki, tn );

        tm -= yi;
        tm += y;

        mult( tn, KiTSeInv, tm );
        solve( dx, KiTSeInvKi, tn );

        x = xa;
        x += dx;

        // Check for convergence.
        tm = y;
        tm -= yi;
        mult( sInvDx, KiTSeInv, tm );
        mult( tn, SaInv, dx_old );
        sInvDx -= tn;

        dx -= dx_old;
        di2 = dx * sInvDx;

        if ( di2 <= tol * (Numeric) m )
        {
            converged = true;
        }

        // Increase iteration counter and store yi.
        iter++;
        dx_old = dx;
    }
    return converged;
}

//! Non-linear OEM using Levenberg-Marquardt method.
/*!
  Inverts a given non-linear forward model using the Levenberg-Marquardt method.
  The implementation follows Eq. (5.36) in Rodger's book.

  Communication with the forward model is performed in the same way as in
  oem_gauss_newton().

  The start value for gamma is given by the parameter gamma_start. If a new x
  value is found, gamma is decreased by a factor gamma_scale_dec. If the cost
  is increased, gamma is increased by a factor gamma_scale_inc. If gamma falls
  below gamma_threshold, it is set to zero. If no lower cost can be obtained
  with gamma = 0, gamma is set to 1. If gamma becomes larger than gamma_max and
  the cost can not be reduced, the iteration is stopped.
  During the execution, space for up to 6 additional matrices and vectors is
  allocated. The given Matrix and Vector views should not overlap.

  \param[out] x The optimal estimator of the state vector.
  \param[in] y The measurement vector.
  \param xa The a priori state vector.
  \param K A forward model object implementing the FowardModel class.
  \param Se The covariance matrix of the measurement error.
  \param Sa The covariance matrix of the a prioi error.
  \param tol The convergence criterion before scaling by the problem size.
  \param maxIter The maximum number of iteration before abortion.
  \param gamma_start The start value of the gamma factor.
  \param gamma_scale_dec The factor to decrease gamma by if the cost function
  was descreased.
  \param gamma_scale_inc The factor to increase gamma by if the cost function
  could not be decreased.
  \param gamma_max The maximum gamma value. If gamma == gamma_max and the cost
  function can not be decreased, the iteration is aborted.
  \param gamma_threshold The minimum value that gamma can take on before it is
  set to zero.

  \return True if the method has converged, false otherwise.
*/
bool oem_levenberg_marquardt( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel &K,
                              ConstMatrixView Se,
                              ConstMatrixView Sa,
                              Numeric tol,
                              Index maxIter,
                              Numeric gamma_start,
                              Numeric gamma_scale_dec,
                              Numeric gamma_scale_inc,
                              Numeric gamma_max,
                              Numeric gamma_threshold )
{

    Index n( x.nelem() ), m( y.nelem() );
    Numeric di2, cost, cost_old, gamma;

    Vector xnew(n), yi(m);
    Vector tm(m), tn1(n), tn2(n), tn3(n), dx(n);
    Matrix Ki(m,n);
    Matrix SeInv(m,m), SaInv(n,n);
    Matrix KiTSeInv(n,m), KiTSeInvKi(n,n), KiTSeInvKi2(n,n);

    // Compute frequently used Matrix inverses.
    inv(SeInv, Se);
    inv(SaInv, Sa);

    // Set starting vector.
    x = xa;

    gamma = gamma_start;
    cost_old = 0.0;

    Index iter = 0;
    bool converged = false;

    while ( (!converged) && (iter < maxIter) )
    {
        // Compute Jacobian and y_i.
        K.evaluate_jacobian( yi, Ki, x);

        mult( KiTSeInv, transpose(Ki), SeInv );
        mult( KiTSeInvKi, KiTSeInv, Ki );

        tn1 = x;
        tn1 -= xa;
        mult( tn2, SaInv, tn1 );

        tm = y;
        tm -= yi;
        mult( tn1, KiTSeInv, tm );

        tn1 -= tn2;

        // Compute old_cost for first iteration.
        if (iter == 0)
        {
            mult( yi, SeInv, tm );
            cost_old = tm * yi;
        }

        bool found_x = false;
        while ( !found_x )
        {
            KiTSeInvKi2 = SaInv;
            KiTSeInvKi2 *= ( 1.0 + gamma );
            KiTSeInvKi2 += KiTSeInvKi;

            // This vector is used to test for convergence later.
            // See eqn. (5.31).

            solve( dx, KiTSeInvKi2, tn1 );
            xnew = x;
            xnew += dx;

            // Evaluate cost function.

            K.evaluate( yi, xnew );

            tm = y;
            tm -= yi;
            mult( yi, SeInv, tm );
            cost = tm*yi;

            tn2 = xnew;
            tn2 -= xa;
            mult( tn3, SaInv, tn2 );
            cost += tn3 * tn2;

            // If cost has decreased, keep new x and
            // scale gamma.
            if (cost < cost_old)
            {

                if ( gamma >= (gamma_threshold * gamma_scale_dec))
                    gamma /= gamma_scale_dec;
                else
                    gamma = 0;

                x = xnew;
                cost_old = cost;
                found_x = true;

            }
            // Else try to increase gamma and look for a
            // new x.
            else
            {
                if ( gamma < gamma_threshold )
                           gamma = gamma_threshold;
                else
                {
                    if ( gamma < gamma_max )
                    {
                        gamma *= gamma_scale_inc;
                        if (gamma > gamma_max)
                            gamma = gamma_max;
                    }
                    // Gamma exceeds maximum. Abort.
                    else
                    {
                        iter = maxIter;
                        break;
                    }
                }
            }
        }

        // Increase iteration counter and check for convergence.
        iter++;

        di2 = dx * tn1;
        if ( di2 <= tol * (Numeric) n )
        {
            converged = true;
        }

    }

    return converged;
}
