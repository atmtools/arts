/*!
  \file   oem.cc
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri Apr 17 16:39:25 2015

  \brief Optimal inversion methods for linear and non-linear inverse problems.
*/

#include "arts.h"
#include <iostream>
#include "lin_alg.h"
#include "logic.h"
#include "math.h"
//#include "oem.h"
#include "stdlib.h"
#include "arts_omp.h"

using std::ostream;
using std::endl;
using std::setw;
using std::scientific;
using std::fixed;


//------------------------------------------------------------------------------------
//
//   Functions for displaying progress of inversions
//
//------------------------------------------------------------------------------------

void separator( ostream& stream,
                Index length )
{
    for (Index i = 0; i < length; i++)
        stream << "-";
    stream << endl;
}

//! Initial log message, linear
/*!

  \param[in] stream Stream to print message to.
*/
void log_init_li( ostream& stream )
{
    stream << endl;
    stream << "Starting linear OEM inversion:" << endl << endl;
    separator( stream, 52 );
    stream << setw(6) << "Step";
    stream << setw(15) << "     Chi^2     ";
    stream << setw(15) << "     Chi^2_x   ";
    stream << setw(15) << "     Chi^2_y   ";
    stream << endl;
    separator( stream, 52 );
}

//! Step log message, linear
/*!

  \param[in] stream Stream to print message to.
  \param[in] step_number Current step number.
  \param[in] cost Current value of cost function.
  \param[in] cost_x x-component of current cost function value.
  \param[in] cost_y y-component of current cost function value.
*/
void log_step_li( ostream& stream,
                  Index step_number,
                  Numeric cost,
                  Numeric cost_x,
                  Numeric cost_y )
{
    stream << setw(5) << step_number << " ";
    stream << scientific;
    stream << setw(15) << cost;
    stream << setw(15) << cost_x;
    stream << setw(15) << cost_y;
    stream << endl;
}

//! Final log message, linear
/*!

  \param[in] stream Stream to print message to.
  \param[in] converged converged flag.
  \param[in] cost Value of cost function.
  \param[in] cost_x x-component of cost function.
  \param[in] cost_y y-component of cost function.
  \param[in] iter Final step number.
  \param[in] max_iter  Maximum step number.
*/
void log_finalize_li( ostream& stream )
{
    separator( stream, 52 );
    stream << endl;
}

//! Initial log message, Gauss-Newton
/*!

  \param[in] stream Stream to print message to.
  \param[in] tol Tolerance.
  \param[in] max_iter Maximum Iterations.
*/
void log_init_gn( ostream& stream,
                  Numeric tol,
                  Index max_iter )
{
    stream << endl;
    stream << "Starting OEM inversion: " << endl << endl;
    stream << "\tMethod: Gauss-Newton" << endl;
    stream << "\tStop criterion: " << tol << endl;
    stream << "\tMax. iterations: " << max_iter << endl << endl;
    separator( stream, 67 );
    stream << setw(6) << "Step";
    stream << setw(15) << "     Chi^2     ";
    stream << setw(15) << "     Chi^2_x   ";
    stream << setw(15) << "     Chi^2_y   ";
    stream << setw(15) << "     d_i^2     ";
    stream << endl;
    separator( stream, 67 );
}

//! Step log message, Gauss-Newton
/*!

  \param[in] stream Stream to print message to.
  \param[in] step_number Current step number.
  \param[in] cost Current value of cost function.
  \param[in] cost_x x-component of current cost function value.
  \param[in] cost_y y-component of current cost function value.
  \param[in] di2 Convergence criterion.
*/
void log_step_gn( ostream& stream,
                  Index step_number,
                  Numeric cost,
                  Numeric cost_x,
                  Numeric cost_y,
                  Numeric di2 )
{
    stream << setw(5) << step_number << " ";
    stream << scientific;
    stream << setw(15) << cost;
    stream << setw(15) << cost_x;
    stream << setw(15) << cost_y;
    if (di2 < 0)
        stream << setw(15) << "";
    else
        stream << setw(15) << di2;
    stream << endl;
}

//! Final log message, Gauss-Newton
/*!

  \param[in] stream Stream to print message to.
  \param[in] converged converged flag.
  \param[in] cost Value of cost function.
  \param[in] cost_x x-component of cost function.
  \param[in] cost_y y-component of cost function.
  \param[in] iter Final step number.
  \param[in] max_iter  Maximum step number.
*/
void log_finalize_gn( ostream& stream,
                      bool converged,
                      Index iter,
                      Index max_iter )
{
    separator( stream, 67 );
    stream << endl;
    if ( converged )
        stream << "\tConverged: YES";
    else if ( iter == max_iter )
    {
        stream << "\tConverged: NO" << endl;
        stream << "\tMaximum no. of iterations was reached.";
    }
    stream << endl << endl;
}

//! Initial log message, Levenberg-Marquardt
/*!

  \param[in] stream Stream to print message to.
  \param[in] tol Tolerance.
  \param[in] max_iter Maximum Iterations.
*/
void log_init_lm( ostream& stream,
                  Numeric tol,
                  Index max_iter )
{
    stream << endl;
    stream << "Starting OEM inversion: " << endl  << endl;
    stream << "     Method: Levenberg-Marquardt " << endl;
    stream << "     Stop criterion: " << tol << endl;
    stream << "     Max. iterations: " << max_iter << endl << endl;
    separator( stream, 75 );
    stream << setw(6) << "Step";
    stream << setw(15) << "     Chi^2      ";
    stream << setw(15) << "     Chi^2_x    ";
    stream << setw(15) << "     Chi^2_y    ";
    stream << setw(9) << "  gamma  ";
    stream << setw(15) << "     d_i^2     ";
    stream << endl;
    separator( stream, 75 );
}

void log_gamma_step_lm( ostream& stream,
                        Numeric cost,
                        Numeric cost_x,
                        Numeric cost_y,
                        Numeric gamma )
{
    stream << setw(6) << "";
    stream << scientific;
    stream << setw(15) << cost;
    stream << setw(15) << cost_x;
    stream << setw(15) << cost_y;
    stream.unsetf(ios_base::floatfield);
    stream << setw(9) << gamma;
    stream << endl;
}

//! Step log message, Levenberg-Marquardt
/*!

  \param[in] stream Stream to print message to.
  \param[in] step_number Current step number.
  \param[in] cost Current value of cost function.
  \param[in] cost_x x-component of current cost function value.
  \param[in] cost_y y-component of current cost function value.
  \param[in] gamma Current gamma value.
  \param[in] di2 Convergence criterion.
*/
void log_step_lm( ostream& stream,
                  Index step_number,
                  Numeric cost,
                  Numeric cost_x,
                  Numeric cost_y,
                  Numeric gamma,
                  Numeric di2 )
{
    stream << setw(6) << step_number;
    stream << scientific;
    stream << setw(15) << cost;
    stream << setw(15) << cost_x;
    stream << setw(15) << cost_y;
    stream.unsetf(ios_base::floatfield);
    stream << setw(9) << gamma;
    stream << scientific;
    stream << setw(15) << di2;
    stream << endl;
}

//! Final log message, Levenberg-Marquardt
/*!

  \param[in] stream Stream to print message to.
  \param[in] converged converged flag.
  \param[in] cost Value of cost function.
  \param[in] cost_x x-component of cost function.
  \param[in] cost_y y-component of cost function.
  \param[in] gamma Final gamma value.
  \param[in] gamma_max Maximum gamma value.
  \param[in] iter Final step number.
  \param[in] max_iter  Maximum step number.
*/
void log_finalize_lm( ostream& stream,
                      bool converged,
                      Numeric cost,
                      Numeric cost_x,
                      Numeric cost_y,
                      Numeric gamma,
                      Numeric gamma_max,
                      Index iter,
                      Index max_iter )
{
    separator( stream, 75 );
    stream << endl << "Finished Levenberg-Marquardt iterations:" << endl;
    if ( converged )
        stream << "\t Converged: YES" << endl;
    else if ( gamma == gamma_max )
    {
        stream << "\t Converged: NO" << endl;
        stream << "\t Iteration aborted because gamma > gamma_max.";
    }
    else if ( iter == max_iter )
    {
        stream << "\t Converged: NO" << endl;
        stream << "\t Iteration aborted because maximum no. of iterations was reached.";
    }
    stream << endl;

    stream << "\t Chi^2: " << cost << endl;
    stream << "\t Chi^2_x: " << cost_x << endl;
    stream << "\t Chi^2_y: " << cost_y << endl;
    stream << endl << endl << endl;
}

//------------------------------------------------------------------------------------
//
//   Calculation of the cost function
//
//------------------------------------------------------------------------------------

//! Calculation of y-part of cost function
/*!
  This version is suitable if no term is already at hand.

  \param[out] cost_y    Const function value
  \param[in]  y         Measurement vector.
  \param[in]  yf        Fitted measurement.
  \param[in]  SeInv     Inverse of relevenaty covariance matrix
  \param[in]  normfac   Normalisation factor. The cost is scaled with 1/normfac

  \author Patrick Eriksson
  \date   2015-10-05
*/
void oem_cost_y( Numeric& cost_y,
                 ConstVectorView y,
                 ConstVectorView yf,
                 ConstMatrixView SeInv,
                 const Numeric&  normfac )
{
  Vector dy = y; dy -= yf;
  Vector tmp( y.nelem() );
  mult( tmp, SeInv, dy );
  cost_y = dy * tmp;
  cost_y /= normfac;
}

//! Calculation of x-part of cost function
/*!
  This version is suitable if no term is already at hand.

  \param[out] cost_x    Const function value
  \param[in]  x         State vector.
  \param[in]  xa        A priori state.
  \param[in]  SxInv     Inverse of relevenaty covariance matrix
  \param[in]  normfac   Normalisation factor. The cost is scaled with 1/normfac

  \author Patrick Eriksson
  \date   2015-10-05
*/
void oem_cost_x( Numeric& cost_x,
                 ConstVectorView x,
                 ConstVectorView xa,
                 ConstMatrixView SxInv,
                 const Numeric&  normfac )
{

  Vector dx = x; dx -= xa;

  Vector tmp( x.nelem() );
  mult( tmp, SxInv, dx );
  cost_x = dx * tmp;
  cost_x /= normfac;
}

//------------------------------------------------------------------------------------
//
//  Algebra Helper Functions
//
//------------------------------------------------------------------------------------

//! Multiply matrix element-wise by outer product of vector.
/*!

  Multiplies the matrix B by the outer product b' * b of b. This is used to
  scale the a priori covariance matrix to avoid numerical problems. The Matrices
  A and B can be the same but shouldn't overlap in other ways.

  \param[out] A The matrix B divided by b'*b.
  \param[in] B  The matrix B to divide.
  \param[in] b The vector used to form the outer product.
*/
void mult_outer( MatrixView A,
                 ConstMatrixView B,
                 ConstVectorView b )
{
    Index n;

    n = b.nelem();

    // A,B must be n x n.
    assert( is_size( A, n, n ) );
    assert( is_size( B, n, n ) );

    for ( Index i = 0; i < n; i++)
    {
        for ( Index j = 0; j < n; j++)
        {
            A(i,j) = B(i,j) * b[i] * b[j];
        }
    }
}

//! Scale columns.
/*!
  Scales the columns of a matrix B by the elements of a vector b.

  \param A The scaled matrix
  \param B The matrix to scale
  \param b The vector containing the scaling factors
*/
void scale_columns( MatrixView A,
                    ConstMatrixView B,
                    ConstVectorView b )
{
    Index m,n;

    m = A.nrows();
    n = A.ncols();

    // A,B must be n x n.
    assert( is_size( B, m, n ) );
    assert( is_size( b, n ) );

    for ( Index i = 0; i < n; i++)
    {
        for ( Index j = 0; j < n; j++)
        {
            A(i,j) = B(i,j) * b[j];
        }
    }
}

//! Scale rows.
/*!
  Scales the rows of a matrix B by the elements of a vector b.

  \param A The scaled matrix
  \param B The matrix to scale
  \param b The vector containing the scaling factors
*/
void scale_rows( MatrixView A,
                 ConstMatrixView B,
                 ConstVectorView b )
{
    Index m,n;

    m = A.nrows();
    n = A.ncols();

    // A,B must be n x n.
    assert( is_size( B, m, n ) );
    assert( is_size( b, n ) );

    for ( Index i = 0; i < n; i++)
    {
        for ( Index j = 0; j < n; j++)
        {
            A(i,j) = B(i,j) * b[i];
        }
    }
}

//------------------------------------------------------------------------------------
//
//  Linear OEM Class
//
//------------------------------------------------------------------------------------

//! Constructor
/*!
  Construct a LinearOEM instace for the computation of the optimal estimator of
  a linear forward model described by a Jacobian J, the inverse of the measurement
  space covariance matrix SeInv, an a priori state vector and the inverse of the
  state space covariance matrix SxInv.

  The constructor allocates the memory that is necessary for the computation of
  the optimal estimator without the computation of the gain matrix. That is,
  space for a n x m and a n x n matrix is allocated as well as space for a
  vector of length n and an integer array of length n.

  \param J The Jacobian
  \param SeInv The measurement space covariance matrix
  \param xa The a priori vector
  \param SxInv The state space covariance matrix
*/
LinearOEM::LinearOEM( ConstMatrixView J_,
                      ConstMatrixView SeInv_,
                      ConstVectorView xa_,
                      ConstMatrixView SxInv_ )
    : J(J_), SeInv(SeInv_), SxInv(SxInv_), x_norm()
{
    n = J.ncols();
    m = J.nrows();

    assert( is_size( SeInv_, m, m ) );
    assert( is_size( SxInv_, n, n ) );
    assert( is_size( xa_, n ) );

    xa = xa_;

    matrices_set = false;
    gain_set = false;
    x_norm_set = false;
    form = NFORM;

    // Allocate memory for matrices.
    LU = Matrix( n, n );
    indx = ArrayOfIndex( n );
    tmp_nn_1 = Matrix( n, n );
    tmp_nm_1 = Matrix( n, m );
    tmp_n_1 = Vector( n );

}

//! Set normalization vector.
/*!
  Sets the normalization vector that is used to scale the state space
  covariance matrix Sx. Avoids numerical problems due to different scaling of
  variables.

  \param x_norm_ The normalization vector.
*/
void LinearOEM::set_x_norm( ConstVectorView x_norm_ )
{
    x_norm = x_norm_;
    x_norm_set = true;
    matrices_set = false;
    gain_set = false;

}

//! Return normalization vector.
/*!
  Returns view of the current normalization vector.

  \return ConstVectorView of the current normalization vector.
*/
ConstVectorView LinearOEM::get_x_norm()
{
    return x_norm;
}

//! Compute optimal estimator.
/*!

  Compute optimal estimator x for a given measurement vector y using the n-form
  formulation of the linear OEM problem without computing the gain matrix.

  This reduced form uses only matrix-vector multiplication and solution of a
  simple linear system of size n x n and is therefore faster than the
  computation that uses the gain matrix.

  If the computation has already been performed once for another measurement
  vector, intermediate results can be reused and the computation only requires
  the solution of an already QR decomposed linear system by backsubstitution.

  \param x The optimal estimator.
  \param y The measurement vector
  \param y0 The offset vector at the linearization point.

  \return Error code indicating success or failure of the computation.
*/
Index LinearOEM::compute( Vector &x,
                         ConstVectorView y,
                         ConstVectorView y0 )
{

    // Set up timers.
    Index t1, t2, t3;
    t1 = timer.add_timer( "State vector computation" );
    timer.mark( t1 );
    t2 = timer.add_timer( "Matrix multiplication");
    t3 = timer.add_timer( "Linear system solving");

    if (!(matrices_set || gain_set))
    {

        timer.mark( t2 );
        mult( tmp_nm_1, transpose(J), SeInv );
        mult( tmp_nn_1, tmp_nm_1, J);
        timer.mark( t2 );

        tmp_nn_1 += SxInv;

        if ( x_norm_set )
        {
            mult_outer( tmp_nn_1, tmp_nn_1, x_norm );
            scale_rows( tmp_nm_1, tmp_nm_1, x_norm );
        }

        timer.mark( t3 );
        ludcmp( LU, indx, tmp_nn_1 );
        timer.mark( t3 );

        matrices_set = true;
    }

    tmp_m_1 = y;
    tmp_m_1 -= y0;

    if (!gain_set)
    {
        mult( tmp_n_1, tmp_nm_1, tmp_m_1 );

        timer.mark( t3 );
        lubacksub( x, LU, tmp_n_1, indx);
        timer.mark( t3 );

        if (x_norm_set)
            x *= x_norm;
    }
    else
    {
        mult( x, G, tmp_m_1 );
    }

    x += xa;

    timer.mark( t1 );

    return 0;
}

//! Compute optimal estimator.
/*!

  Compute optimal estimator x for a given measurement vector y using the n-form
  formulation of the linear OEM problem with computing the gain matrix.

  The computation of the Gain matrix requires the allocation of an additional
  n-times-n and a n-times-m matrices.

  If the computation has already been performed once for another measurement
  vector, intermediate results can be reused and the computation only requires
  the solution of an already QR decomposed linear system by backsubstitution.

  \param x The optimal estimator.
  \param y The measurement vector
  \param y0 The offset vector at the linearization point.

  \return Error code indicating success or failure of the computation.
*/
Index LinearOEM::compute( Vector &x,
                          MatrixView G_,
                          ConstVectorView y,
                          ConstVectorView y0 )
{

    if (!gain_set)
        compute_gain_matrix();

    compute( x, y, y0 );
    G_ = G;

    return 0;
}

void LinearOEM::compute_gain_matrix()
{

    // Setup timers.
    Index t1, t2, t3;
    t1 = timer.add_timer( "Gain Matrix Computation" );
    t2 = timer.add_timer( "Matrix Matrix Mult.");
    t3 = timer.add_timer( "Linear System" );
    timer.mark( t1 );

    // Assure that G has the right size.
    G.resize( n, m );
    tmp_nn_2.resize( n, n );

    if (!(matrices_set))
    {
        timer.mark( t2 );
        mult( tmp_nm_1, transpose(J), SeInv );
        mult( tmp_nn_1, tmp_nm_1, J);
        timer.mark( t2 );

        tmp_nn_1 += SxInv;

        if ( x_norm_set )
        {
            mult_outer( tmp_nn_1, tmp_nn_1, x_norm );
            scale_rows( tmp_nm_1, tmp_nm_1, x_norm );
        }

        timer.mark( t3 );
        ludcmp( LU, indx, tmp_nn_1 );
        timer.mark( t3 );

        matrices_set = true;
    }

    // Invert J^T S_e^{-1} J using the already computed LU decomposition.
    tmp_n_1 = 0.0;

    timer.mark( t3 );
    for ( Index i = 0; i < n; i++ )
    {
        tmp_n_1[i] = 1.0;
        lubacksub( tmp_nn_2(joker,i), LU, tmp_n_1, indx);
        tmp_n_1[i] = 0.0;
    }
    timer.mark( t3 );

    timer.mark( t2 );
    mult( G, tmp_nn_2, tmp_nm_1 );
    timer.mark( t2 );

    if ( x_norm_set )
    {
        scale_rows( G, G, x_norm );
        matrices_set = false;
    }
    gain_set = true;

    timer.mark( t1 );
}

//! Compute fit
/*!
  Compute fitted measurement vector from a given forward model and estimated
  state vector x.

  \param[out] y The fitted measurement vector y.
  \param[in] x The estimated state vector x.
  \param[in] F The ForwardModel instance representing the forward model.

  \return Error code indicating the success of the operation.
*/
Index LinearOEM::compute_fit( Vector &yf,
                            const Vector &x,
                            ForwardModel &F )
{
    try
    {
        F.evaluate( yf, x );
    }
    catch (int e)
    {
        err = 9;
        return err;
    }

    return 0;
}

//! Compute fit
/*!
  Compute fitted measurement vector and costs from a given forward model and
  estimated state vector x.

  \param[out] y The fitted measurement vector y.
  \param[out] cost_x The cost part corresponding to the state vector x.
  \param[out] cost_y The cost part corresponding to the measurement vector y.
  \param[in] x The estimated state vector x.
  \param[in] F The ForwardModel instance representing the forward model.

  \return Error code indicating the success of the operation.
*/
Index LinearOEM::compute_fit( Vector &yf,
                            Numeric &cost_x,
                            Numeric &cost_y,
                            Vector &x,
                            ConstVectorView y,
                            ForwardModel &F )
{

    try
    {
        F.evaluate( yf, x );
    }
    catch (int e)
    {
        return 9;
    }

    oem_cost_y( cost_y, y, yf, SeInv, (Numeric) m);
    oem_cost_x( cost_x, x, xa, SxInv, (Numeric) m);

    return 0;
}

//------------------------------------------------------------------------------------
//
//  Linear OEM
//
//------------------------------------------------------------------------------------

//! Linear OEM, n-form.
/*!

Computes the inverse of a linear forward model by computing the MAP
solution as described in Rodgers, Inverse Methods for Atmospheric Sounding,
p. 67. This function uses the n-form (eq. (4.3)) which requires the solution
of a linear system of equations of size n-times-n.

Requires the inverses of the covariance matrices for the state and measurement
vector to be provided as arguments.

The returned Index value indicates if the computation was successful or if
an error occured. The following encoding for errors is used:

    0 - Success
    9 - ERROR: Evaluation of forward model failed.

  \return Convergence status, see *oem_diagnostics*
  \param[out] x The optimal, estimated state vector consisting of n elements.
  \param[out] G The gain matrix.
  \param[in] xa The a priori state vector
  \param[in] yf The value of the forward model at a priori.
  \param[in] y The measurement vector consisting of m elements.
  \param[in] K The weighting function (m,n)-matrix
  \param[in] SeInv The inverse of the measurement error covariance (m,m)-matrix
  \param[in] SxInv The inverse of the a priori covariance (n,n)-matrix

  \return Error code indicating success or failure of the method.
*/
Index oem_linear_nform( Vector& x,
                        Matrix& G,
                        Matrix& J,
                        Vector& yf,
                        Numeric& cost_y,
                        Numeric& cost_x,
                        ForwardModel &F,
                        ConstVectorView xa,
                        ConstVectorView x_norm,
                        ConstVectorView y,
                        ConstMatrixView SeInv,
                        ConstMatrixView SxInv,
                        const Numeric& cost_start,
                        const bool& verbose )
{
    const Index m = J.nrows();
    const Index n = J.ncols();

    // Check dimensions for consistency.
    assert( xa.nelem() == n );
    assert( x_norm.nelem() == 0  || x_norm.nelem() == n );
    assert( y.nelem() == m );
    assert( (SeInv.ncols() == m) && (SeInv.nrows() == m) );
    assert( (SxInv.ncols() == n) && (SxInv.nrows() == n) );

    LinearOEM oem( J, SeInv, xa, SxInv );

    if ( x_norm.nelem() == n )
        oem.set_x_norm( x_norm );

    // Initialize log output.
    if (verbose)
      {
        log_init_li( cout );
        log_step_li( cout, 0, cost_start, 0, cost_start );
      }

    oem.compute( x, G, y, yf );
    oem.compute_fit( yf, cost_x, cost_y, x, y, F );

    // Finalize log output.
    if (verbose)
      {
        log_step_li( cout, 1, cost_y+cost_x, cost_x, cost_y );
        log_finalize_li( cout );
      }

    // Return convergence status
    return oem.get_error();
}

//------------------------------------------------------------------------------------
//
// Non-linear OEM Class
//
//------------------------------------------------------------------------------------

//! Create NonLinearOEM object.
/*!
  Set up a non-linear OEM computation. Given the inverse measurement
  and state covariance matrices SeInv, SxInv, the a priori
  state vector xa and a ForwardModel instance representing the forward model,
  this function allocates the memory required for the computation and
  sets the default values for the computation.

  Mermory for the following matrices and vectors is allocated:

  - 3 [n,m] matrix
  - 2 [n,n] matrices
  - 2 n-element vectors
  - 3 n-element vectors

  The default values for the computation are:

  max_iter = 100
  tol = 10e-3
  gamma_start = 2.0
  gamma_decrease = 2.0
  gamma_increase = 3.0
  gamma_max = 100.0
  gamma_threshold = 1.0

  For a description of the functionality of the parameters see the documentation
  of the gauss_newton() and levenberg_marquardt() class methods.

  \param[in] SeInv_ The inverse of the measurement covariance matrix
  \param[in] xa_ The a priori state vector
  \param[in] SxInv_  The inverse of the state covariance matrix
  \param[in] F_ The forward model
  \param[in] method_ The computational method to use for the inversion. Either
  GAUSS_NEWTON or LEVENBERG_MARQUARDT.
*/
template<typename Se_t, typename Sx_t>
NonLinearOEM<Se_t, Sx_t>::NonLinearOEM( const Se_t      &SeInv_,
                                        ConstVectorView  xa_,
                                        const Sx_t      &SxInv_,
                                        ForwardModel    &F_,
                                        OEMMethod        method_ )
    : SeInv(SeInv_), SxInv(SxInv_), xa(xa_), F(F_), method(method_)
{

    m = SeInv_.nrows();
    n = SxInv_.nrows();

    assert( xa_.nelem() == n );
    assert( (SeInv_.ncols() == m) && (SeInv_.nrows() == m) );
    assert( (SxInv_.ncols() == n) && (SxInv_.nrows() == n) );

    xa = xa_;

    // Allocate memory for internal matrices and vectors.
    J = Matrix(m,n);
    G = Matrix(n,m);
    tmp_nm_1 = Matrix(n,m);
    tmp_nn_1 = Matrix(n,n);
    tmp_nn_2 = Matrix(n,n);

    tmp_m_1 = Vector(m);
    tmp_n_1 = Vector(n);
    tmp_n_2 = Vector(n);
    dx = Vector(n);
    yi = Vector(m);

    // Initialize internal state variables.
    matrices_set = false;
    gain_set = false;
    x_norm_set = false;
    conv = false;

    // Default iteration parameters.
    tol = 1e-5;
    max_iter = 100;

    ga_start = 4.0;
    ga_max = 100.0;
    ga_decrease = 2.0;
    ga_increase = 3.0;
    ga_threshold = 1.0;

}

//! Set normalization vector.
/*!
  Set the vector used to normalize the state covariance matrix to improve
  numerical robustness.

  \param x_norm_ The normalizatoin vector.
*/
template<typename Se_t, typename Sx_t>
void NonLinearOEM<Se_t, Sx_t>::set_x_norm( ConstVectorView x_norm_ )
{
    x_norm = x_norm_;
    x_norm_set = true;
}

//! Return current normalization vector.
/*!
  \return The current normalization vector
*/
template<typename Se_t, typename Sx_t>
ConstVectorView NonLinearOEM<Se_t, Sx_t>::get_x_norm()
{
    return x_norm;
}

//! Perform OEM computation.
/*!
  Computes the Bayesian optimal estimator for the state vector x, given a
  measurement vector y. The computation is performed using the n-form of the
  chosen computational method (Gauss-Newton or Levenberg-Marquardt) as given
  in equation (5.8) by Rodgers.

  The iteration aborts when the maximum number of iterations max_iter is reached,
  the criterion (5.29) in Rodgers book is satisfied or when the maximum gamma
  value is reached lofs.com/(Levenberg-Marquardt).

  \param x The Bayesian optimal estimator for the measurement vector y
  \param y The measurement vector
  \param verbose If true, print log to standard out.
*/
template<typename Se_t, typename Sx_t>
Index NonLinearOEM<Se_t, Sx_t>::compute( Vector &x,
                                         ConstVectorView y,
                                         bool verbose )
{

    if (method == GAUSS_NEWTON)
    {
        gauss_newton( x, y, verbose );
    }
    else if (method == LEVENBERG_MARQUARDT)
    {
        levenberg_marquardt( x, y, verbose );
    }

    if (verbose)
        cout << timer << endl;

    return err;
}

//! Perform OEM computation.
/*!
  Performs the same computations as the standard compute method and also
  computes the gain matrix.

  \param x The Bayesian optimal estimator for the measurement vector y
  \param G_ The gain matrix
  \param y The measurement vector
  \param verbose If true, print log to standard out.
*/
template<typename Se_t, typename Sx_t>
Index NonLinearOEM<Se_t, Sx_t>::compute( Vector          &x,
                                         MatrixView       G_,
                                         ConstVectorView  y,
                                         bool             verbose )
{

    if (method == GAUSS_NEWTON)
    {
        gauss_newton( x, y, verbose );
    }
    else if (method == LEVENBERG_MARQUARDT)
    {
        levenberg_marquardt( x, y, verbose );
    }

    compute_gain_matrix( x );
    G_ = G;

    if (verbose)
        cout << timer << endl;

    return err;
}

//! Gauss-Newton method.
/*!
  Solve the given non-linear inverse problem using Gauss-Newton method. The
  algorithm is a straight-forward implementation of the method described in
  section 5.3 in Rodgers' book. The iteration terminates when the maximum
  number of iterations is reached or when the convergence criterion given
  by Rodgers in equation (5.29) is satisfied up to the given tolerance.

  \param x The optimal estimator for the state vector.
  \param y The measurement vector.
  \param verbose If true, print log output to standard out.
*/
template<typename Se_t, typename Sx_t>
void NonLinearOEM<Se_t, Sx_t>::gauss_newton( Vector          &x,
                                             ConstVectorView  y,
                                             bool             verbose )
{

    Index t1, t2, t3, t4;
    t1 = timer.add_timer( "Gauss-Newton Iteration" );
    t2 = timer.add_timer( "Matrix Matrix Mult." );
    t3 = timer.add_timer( "Linear System" );
    t4 = timer.add_timer( "Jacobian Evaluation" );
    timer.mark( t1 );

    Numeric di2 = -1.0;

    cost_x = 0.0; cost_y = 0.0;
    iter = 0;
    conv = false;
    err = 1;

    // Initialize log output.
    if (verbose)
    {
          log_init_gn( cout, tol, max_iter );
    }

    // Start stop watch.

    // Set the starting vector.
    x = xa;

    while ( (!conv) && (iter < max_iter) )
    {
        // Compute Jacobian and y_i.

        timer.mark( t1 );
        timer.mark( t4 );
        try
        {
            F.evaluate_jacobian( yi, J, x);
        }
        catch (int e)
        {
            err = 9;
            return void();
        }
        timer.mark( t1 );
        timer.mark( t4 );


        timer.mark( t2 );
        mult( tmp_nm_1, transpose(J), SeInv );
        mult( tmp_nn_1, tmp_nm_1, J );
        timer.mark( t2 );

        tmp_nn_1 += SxInv;

        if (x_norm_set)
        {
            mult_outer( tmp_nn_1, tmp_nn_1, x_norm );
            scale_rows( tmp_nm_1, tmp_nm_1, x_norm );
        }

        // tm1 = K_i(x_i - x_a)

        tmp_n_1 = x;
        tmp_n_1 -= xa;
        mult( tmp_n_2, SxInv, tmp_n_1 );

        if (x_norm_set)
            tmp_n_2 *= x_norm;

        tmp_m_1 = y;
        tmp_m_1 -= yi;
        // tmp_n_1 = K_i^T * S_e^{-1} * (y - F(x_i))
        mult( tmp_n_1, tmp_nm_1, tmp_m_1 );

        // This vector is used to test for convergence later.
        // See eqn. (5.31).
        tmp_n_1 -= tmp_n_2;

        // Compute cost function and convergence criterion.
        oem_cost_x( cost_x, x, xa, SxInv, (Numeric) m);
        oem_cost_y( cost_y, y, yi, SeInv, (Numeric) m);

        // Print log.
        if ( verbose )
            log_step_gn( cout, iter, cost_x + cost_y, cost_x, cost_y, di2 );

        timer.mark( t3 );
        solve( dx, tmp_nn_1, tmp_n_1 );
        timer.mark( t3 );

        if (x_norm_set)
            dx *= x_norm;
        x += dx;
        iter++;

        // Compute convergence criterion.
        di2 = dx * tmp_n_1;

        // If converged, stop iteration.
        if ( di2 <= tol * (Numeric) n )
        {
            conv = true;
            err = 0;
        }
    }

    if ( iter == max_iter )
        err = 1;

    // Stop stop watch.

    // Finalize log output.
    if ( verbose )
        log_finalize_gn( cout, conv, iter, max_iter );

}

//! Levenberg-Marquardt method.
/*!
  Solve the given non-linear inverse problem using Levenberg-Marquardt method.
  The algorithm is a straight-forward implementation of the method described in
  section 5.3 in Rodgers' book. The iteration terminates when the maximum
  number of iterations is reached, the convergence criterion given
  by Rodgers in equation (5.29) is satisfied up to the given tolerance or when
  the given maximum gamma value is reached.

  The Levenberg-Marquardt method uses a scaling gamma factor in the computation
  of the current iteration step dx. If the current gamma leads to a new x vector

      x_{i+1} = x_{i} + dx

  that decreases the cost function, then the new x-value x_{i+1} is accepted and
  gamma is decrease by a factor gamma_decrease. If the value of the cost
  function for x_{i+1} has increased, gamma is increased by the factor
  gamma_increase and a new x_{i+1} is computed.

  If gamma decreases below gamma_threshold, gamma is set to 0.0. If the current
  gamma is 0.0 and is increased, it is set to 1.0. If gamma increases to more
  than gamma_max, it is set to gamma_max. If the current gamma is gamma_max and
  no new value for x can be found that decreases the cost function, the
  iteration is aborted.

  \param x The optimal estimator for the state vector.
  \param y The measurement vector.
  \param verbose If true, print log output to standard out.
*/
template<typename Se_t, typename Sx_t>
void NonLinearOEM<Se_t, Sx_t>::levenberg_marquardt( Vector          &x,
                                                    ConstVectorView  y,
                                                    bool             verbose )
{
    // Setup timers.
    Index t1, t2, t3, t4;
    t1 = timer.add_timer( "Levenberg-Marquardt Iteration" );
    t2 = timer.add_timer( "Matrix Matrix Mult." );
    t3 = timer.add_timer( "Linear System" );
    t4 = timer.add_timer( "Jacobian Evaluation" );
    timer.mark( t1 );

    Numeric cost_old, cost, di2;
    di2 = -1.0;

    if ( verbose )
        log_init_lm( cout, tol, max_iter );

    // Set starting vector.
    x = xa;

    conv = false;
    iter = 0;
    cost_x = 0.0;
    cost_y = 0.0;
    Numeric gamma = ga_start;

    while ( (!conv) && (iter < max_iter) && (gamma <= ga_max) )
    {
        // Compute Jacobian and y_i.
        timer.mark( t1 );
        timer.mark( t4 );
        try
        {
            F.evaluate_jacobian( yi, J, x);
        }
        catch (int e)
        {
            err = 9;
            return void();
        }
        timer.mark( t1 );
        timer.mark( t4 );

        timer.mark(t2);
        mult( tmp_nm_1, transpose(J), SeInv );
        mult( tmp_nn_1, tmp_nm_1, J );
        timer.mark(t2);

        tmp_n_1 = x;
        tmp_n_1 -= xa;
        mult( tmp_n_2, SxInv, tmp_n_1 );

        tmp_m_1 = y;
        tmp_m_1 -= yi;
        mult( tmp_n_1, tmp_nm_1, tmp_m_1 );

        tmp_n_1 -= tmp_n_2;
        tmp_n_2 = tmp_n_1;
        if (x_norm_set)
            tmp_n_1 *= x_norm;

        // Compute costs and print log message.
        oem_cost_x( cost_x, x, xa, SxInv, (Numeric) m);
        oem_cost_y( cost_y, y, yi, SeInv, (Numeric) m);
        cost_old = cost_y + cost_x;

        if ( verbose )
            log_step_gn( cout, iter, cost_x + cost_y, cost_x, cost_y, di2 );

        bool found_x = false;
        while ( !found_x )
        {
            tmp_nn_2 = SxInv;
            tmp_nn_2 *= ( 1.0 + gamma );
            tmp_nn_2 += tmp_nn_1;

            if ( x_norm_set )
            {
                mult_outer( tmp_nn_2, tmp_nn_2, x_norm );
            }

            // This vector is used to test for convergence later.
            // See eqn. (5.31).

            timer.mark( t3 );
            solve( dx, tmp_nn_2, tmp_n_1 );
            timer.mark( t3 );

            if (x_norm_set)
                dx *= x_norm;
            xnew = x;
            xnew += dx;

            // Evaluate cost function.


            timer.mark(t1);
            timer.mark(t4);
            try
            {
                F.evaluate( yi, xnew );
            }
            catch ( int e )
            {
                err = 9;
                return void();
            }
            timer.mark(t1);
            timer.mark(t4);


            oem_cost_x( cost_x, xnew, xa, SxInv, (Numeric) m);
            oem_cost_y( cost_y, y, yi, SeInv, (Numeric) m);
            cost = cost_x + cost_y;

            // If cost has decreased, keep new x and
            // scale gamma.
            if (cost < cost_old)
            {

                if ( gamma >= (ga_threshold * ga_decrease))
                    gamma /= ga_decrease;
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
                if ( gamma < ga_threshold )
                           gamma = ga_threshold;
                else
                {
                    if ( gamma < ga_max )
                    {
                        gamma *= ga_increase;
                        if (gamma > ga_max)
                            gamma = ga_max;
                    }
                    // Gamma exceeds maximum. Abort.
                    else
                    {
                        gamma = ga_max + 1.0;
                        break;
                    }
                }
            }
        }

        // Increase iteration counter and check for convergence.
        iter++;

        di2 = dx * tmp_n_2;
        di2 /= (Numeric) n;

        if ( di2 <= tol )
        {
            conv = true;
            err = 0;
       }

    }

    if ( iter == max_iter )
        err = 1;
    if ( gamma >= ga_max )
        err = 2;

    timer.mark( t1 );

    // Final log output.
    if ( verbose )
    {
        cost = cost_x + cost_y;
        log_finalize_lm( cout, conv, cost, cost_x, cost_y,
                         gamma, ga_max, iter, max_iter );
    }
}

//! Compute fit.
/*!
  Compute fit for a given estimated state vector.

  \param yf The fitted measurement vector
  \param x The estimated state vector

  \return Error code indicating success or failure of the computation of
  the fit.
*/
template<typename Se_t, typename Sx_t>
Index NonLinearOEM<Se_t, Sx_t>::compute_fit( Vector                &yf,
                                                      const Vector &x )
{
    try
    {
        F.evaluate( yf, x );
    }
    catch (int e)
    {
        err = 9;
        return err;
    }
    return 0;

}

//! Compute fit and evaluate cost function.
/*!
  Compute the fitted measurement vector for a given state vector and
  evaluate the cost function.

  \param yf The fitted measurement vector
  \param cost_x_ The x-part of the cost function value
  \param cost_y_  The y-part of the cost function value
  \param x The estimated state vector
  \param y The original measurement vector

  \return Error code indicating success or failure of the computation of
  the fit.
*/
template<typename Se_t, typename Sx_t>
Index NonLinearOEM<Se_t, Sx_t>::compute_fit( Vector &yf,
                                             Numeric &cost_x_,
                                             Numeric &cost_y_,
                                             const Vector &x,
                                             ConstVectorView y )
{

    try
    {
        F.evaluate( yf, x );
    }
    catch (int e)
    {
        err = 9;
        return err;
    }

    oem_cost_y( cost_y_, y, yf, SeInv, (Numeric) m);
    oem_cost_x( cost_x_, x, xa, SxInv, (Numeric) m);

    return 0;

}

//! Compute gain matrix.
/*!
  Compute the gain matrix for a given state vector x.

  \param x The state vector
*/
template<typename Se_t, typename Sx_t>
void NonLinearOEM<Se_t, Sx_t>::compute_gain_matrix( Vector& x )
{
    Index t1, t2, t3, t4;
    t1 = timer.add_timer( "Gain Matrix Computation" );
    t2 = timer.add_timer( "Matrix Matrix Mult." );
    t3 = timer.add_timer( "Linear System" );
    t4 = timer.add_timer( "Jacobian Evaluation" );
    timer.mark( t1 );

    timer.mark( t4 );
    try
    {
        F.evaluate_jacobian( tmp_m_1, J, x);
    }
    catch (int e)
    {
        err = 9;
        return void();
    }
    timer.mark( t4 );

    timer.mark( t2 );
    mult( tmp_nm_1, transpose(J), SeInv );
    mult( tmp_nn_1, tmp_nm_1, J );
    timer.mark( t2 );

    tmp_nn_1 += SxInv;

    if (x_norm_set)
    {
        mult_outer( tmp_nn_1, tmp_nn_1, x_norm );
        scale_columns( tmp_nm_1, tmp_nm_1, x_norm );
    }

    timer.mark( t3 );
    inv( tmp_nn_2, tmp_nn_1 );
    timer.mark( t3 );

    timer.mark( t2 );
    mult( G, tmp_nn_2, tmp_nm_1 );
    timer.mark( t2 );

    timer.mark( t1 );
}

//! Gauss-Newton non-linear OEM using precomputed inverses, n-form.
/*!
  Computes the optimal nonlinear inverse of a given forward model using the
  Gauss-Newton method as given in eq. (5.8) in Rodgers book. The forward model
  is given by an object implementing the interface described by the ForwardModel
  class. Convergence is determined using equation (5.29). The given
  tolerance is scaled by the problem size n. If the method doesn't converge it
  abords after the given maximum number of iterations.

  During execution two additional n-times-m and one n-times-n matrix is
  allocated. In addition to that, space for 4 length-n vectors and two length-m
  vectors is allocated. The given Matrix and Vector views may not overlap.

  The returned Index value indicates if the computation was successful or if
  an error occured. The following encoding for errors is used:

      0 - Success
      1 - ERROR: Computation didn't converge. Maximum of iterations reached.
      2 - ERROR: Computation didn't converge. Maximum gamma value reached.
      9 - ERROR: Evaluation of forward model failed.

  \param[out] x The optimal inverse x.
  \param[out] G The gain matrix corresponding to x.
  \param[out] J The jacobian as computed in the last LM iteration.
  \param[out] yf The fitted measurement vector.
  \param[out] cost_y The cost related to the measurement vector.
  \param[out] cost_x The cost related to the state vector.
  \param[out] iter The number of iterations.
  \param[in] F The ForwardModel representing the forward model to invert.
  \param[in] xa The size-n a-priori state vector.
  \param[in] x_norm Normalization vector for the state covariance matrix.
  \param[in] y The length-m, input measurement vector.
  \param[in] SeInv The inverse of the measurement error covariance (m,m)-matrix
  \param[in] SxInv The inverse of the a priori covariance (n,n)-matrix
  \param[in] max_iter Tha maximum number of iterations in case of no convergence.
  \param[in] tol The convergence criterium before scaling by the problem size.
  \param[in] verbose If true, log message are printed to stdout.

  \return Error code indicating success or failure of the method.
*/
template <typename Se_t, typename Sx_t>
Index oem_gauss_newton( Vector& x,
                        Matrix& G,
                        Matrix& J,
                        Vector& yf,
                        Numeric& cost_y,
                        Numeric& cost_x,
                        Index& iter,
                        ForwardModel &F,
                        ConstVectorView xa,
                        ConstVectorView x_norm,
                        ConstVectorView y,
                        const Se_t &SeInv,
                        const Sx_t &SxInv,
                        const Index max_iter,
                        const Numeric tol,
                        bool verbose )
{
    const Index m = J.nrows();
    const Index n = J.ncols();

    // Check dimensions for consistency.
    assert( xa.nelem() == n );
    assert( x_norm.nelem() == 0  || x_norm.nelem() == n );
    assert( y.nelem() == m );
    assert( (SeInv.ncols() == m) && (SeInv.nrows() == m) );
    assert( (SxInv.ncols() == n) && (SxInv.nrows() == n) );

    NonLinearOEM<Se_t, Sx_t>
        oem( SeInv, xa, SxInv, F, GAUSS_NEWTON);
    oem.tolerance( tol );
    oem.maximum_iterations( max_iter );

    if ( x_norm.nelem() == n )
        oem.set_x_norm( x_norm );

    oem.compute( x, G, y, verbose );
    oem.compute_fit( yf, cost_x, cost_y, x, y );

    J = oem.get_jacobian();


    iter = oem.iterations();

    return oem.get_error();
}



//! Non-linear OEM using Levenberg-Marquardt method.
/*!
  Inverts a given non-linear forward model using the Levenberg-Marquardt method.
  The implementation follows Eq. (5.36) in Rodger's book.

  Communication with the forward model is performed in the same way as in
  oem_gauss_newton().

  The start value for gamma is given by the parameter gamma_start. If a new x
  value is found, gamma is decreased by a factor gamma_scale_dec. If the ost
  is increased, gamma is increased by a factor gamma_scale_inc. If gamma falls
  below gamma_threshold, it is set to zero. If no lower cost can be obtained
  with gamma = 0, gamma is set to 1. If gamma becomes larger than gamma_max and
  the cost can not be reduced, the iteration is stopped.

  During the execution, space for two n-times-n and one n-times-m matrices is
  allocated as well as space for 5 length-n vectors and two length-m vectors.

  The returned Index value indicates if the computation was successful or if
  an error occured. The following encoding for errors is used:

      0 - Success
      1 - ERROR: Computation didn't converge. Maximum of iterations reached.
      2 - ERROR: Computation didn't converge. Maximum gamma value reached.
      9 - ERROR: Evaluation of forward model failed.

  \param[out] x The optimal estimator of the state vector.
  \param[out] yf The fitted state vector as computed in the second-last LM
  iteration.
  \param[out] G The gain matrix.
  \param[out] J The jacobian as computed in the second-last LM iteration.
  \param[in] y The size-m input measurement vector.
  \param[in] xa The size-n a priori state vector.
  \param[in] K A forward model object implementing the FowardModel class.
  \param[in] Seinv The inverse of the covariance matrix of the measurement error.
  \param[in] Sxinv The inverse of the covariance matrix of the a priori error.
  \param[in] tol The normalized convergence criterion.
  \param[in] max_iter The maximum number of iterations before abortion.
  \param[in] gamma_start The start value of the gamma factor.
  \param[in] gamma_scale_dec The factor to decrease gamma by if the cost function
  was descreased.
  \param[in] gamma_scale_inc The factor to increase gamma by if the cost function
  could not be decreased.
  \param[in] gamma_max The maximum gamma value. If gamma == gamma_max and the cost
  function can not be decreased, the iteration is aborted.
  \param[in] gamma_threshold The minimum value that gamma can take on before it is
  set to zero.
  \param[in] verbose If true, log messages are printed to stdout.

  \return Error code indicating success or failure of the method.
*/
template <typename Se_t, typename Sx_t>
Index oem_levenberg_marquardt( Vector& x,
                               Matrix& G,
                               Matrix& J,
                               Vector& yf,
                               Numeric& cost_y,
                               Numeric& cost_x,
                               Index& iter,
                               ForwardModel &F,
                               ConstVectorView xa,
                               ConstVectorView x_norm,
                               ConstVectorView y,
                               const Se_t &SeInv,
                               const Sx_t &SxInv,
                               Index max_iter,
                               Numeric tol,
                               Numeric gamma_start,
                               Numeric gamma_decrease,
                               Numeric gamma_increase,
                               Numeric gamma_max,
                               Numeric gamma_threshold,
                               bool verbose )
{

    const Index m = J.nrows();
    const Index n = J.ncols();

    // Check dimensions for consistency.
    assert( xa.nelem() == n );
    assert( x_norm.nelem() == 0  || x_norm.nelem() == n );
    assert( y.nelem() == m );
    assert( (SeInv.ncols() == m) && (SeInv.nrows() == m) );
    assert( (SxInv.ncols() == n) && (SxInv.nrows() == n) );

    NonLinearOEM<Se_t, Sx_t>
        oem( SeInv, xa, SxInv, F, LEVENBERG_MARQUARDT);
    oem.tolerance( tol );
    oem.maximum_iterations( max_iter );
    oem.gamma_start( gamma_start );
    oem.gamma_decrease( gamma_decrease );
    oem.gamma_increase( gamma_increase );
    oem.gamma_max( gamma_max );
    oem.gamma_threshold( gamma_threshold );

    if ( x_norm.nelem() == n )
        oem.set_x_norm( x_norm );

    oem.compute( x, G, y, verbose );
    oem.compute_fit( yf, cost_x, cost_y, x, y );

    J = oem.get_jacobian();

    iter = oem.iterations();

    return oem.get_error();

}

//------------------------------------------------------------------------------------
//
//   OEM versions so far not used
//
//------------------------------------------------------------------------------------


//! Linear OEM, m-form.
/*!

Computes the inverse of a linear forward model by computing the MAP
solution as described in Rodgers, Inverse Methods for Atmospheric Sounding,
p. 67. This function uses the m-form ( Eq. (4.6) ) which requires the solution
of a linear system of equations of size m-times-m.

For the execution 1 n-times-m matrices, 2 m-times-m matrices and a vector with m
elements are allocated. The given Matrix and Vector views may not overlap.

  \param[out] x The optimal, estimated state vector consisting of n elements.
  \param[in] y The measurement vector consisting of m elements.
  \param[in,out] yf On input yf should contain the value of the forward model
  at the linearization point. On output yf should contain the fitted measurement
  vector.
  \param[in] xa The mean a priori state vector.
  \param[in] K The weighting function (m,n)-matrix.
  \param[in] Se The measurement error covariance (m,m)-matrix.
  \param[in] Sx The a priori covariance (n,n)-matrix.
  \param[out] G The gain matrix.
*/
void oem_linear_mform( VectorView x,
                       MatrixView G,
                       ConstVectorView xa,
                       ConstVectorView yf,
                       ConstVectorView y,
                       ConstMatrixView K,
                       ConstMatrixView Se,
                       ConstMatrixView Sx )
{
    Index m = K.nrows();
    Index n = K.ncols();

    // Check dimensions for consistency.
    assert( x.nelem() == n );
    assert( xa.nelem() == n );
    assert( y.nelem() == m );

    assert( (Se.ncols() == m) && (Se.nrows() == m) );
    assert( (Sx.ncols() == n) && (Sx.nrows() == n) );

    // m-form (eq. (4.6)).
    Matrix tmp_nm(n,m);
    Matrix tmp_mm(m,m), tmp_mm2(m,m);
    Vector tmp_m(m);

    mult( tmp_nm, Sx, transpose(K) ); // tmp_nm = S_a * K^T
    mult( tmp_mm, K, tmp_nm);
    tmp_mm += Se;

    // Compute gain matrix.
    inv( tmp_mm2, tmp_mm );
    mult( G, tmp_nm, tmp_mm2 );

    tmp_m = y;
    tmp_m -= yf;
    mult( x, G, tmp_m );

    x += xa;
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
  \param[in] Sx The a priori covariance (n,n)-matrix
  \param[in] tol The convergence criterium before scaling by the problem size.
  \param max_iter Tha maximum number of iterations in case of no convergence.
*/
bool oem_gauss_newton_m_form( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel &K,
                              ConstMatrixView Se,
                              ConstMatrixView Sx,
                              Numeric tol,
                              Index max_iter )
{
    Index n( x.nelem() ), m( y.nelem() );
    Numeric di2;
    Matrix Ki(m,n), SxKiT(m, n), KiSxKiT(m, m);
    Vector xi(n), yi(m), tm(m), tm2(m), tn(n), yi_old(m);

    Index iter = 0;
    bool converged = false;

    // Set the starting vector.
    x = xa;

    while ((!converged) && (iter < max_iter))
    {

        // Compute Jacobian and y_i.
        K.evaluate_jacobian( yi, Ki, x );
        // If not in the first iteration, check for
        // convergence here.
        if (iter > 0)
        {
            yi_old -= yi;

            solve( tm, Se, yi_old );
            mult( tm2, KiSxKiT, tm );
            // TODO: Optimize using LU decomp.
            solve( tm, Se, tm2);
            di2 = yi_old * tm;
            if ( fabs( di2 ) < tol * (Numeric) n )
            {
                converged = true;
                break;
            }
        }

        mult( SxKiT, Sx, transpose(Ki) );
        mult( KiSxKiT, Ki, SxKiT );
        KiSxKiT += Se;

        // tm = K_i(x_i - x_a)
        tn = x;
        tn -= xa;
        mult( tm, Ki, tn );

        tm -= yi;
        tm += y;

        solve( tm2, KiSxKiT, tm );
        mult( x, SxKiT, tm2 );
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
  \param[in] Sx The a priori covariance (n,n)-matrix
  \param[in] tol The convergence criterium before scaling by the problem size.
  \param max_iter Tha maximum number of iterations in case of no convergence.

  \return True if the method has converged, false otherwise.
*/
bool oem_gauss_newton_n_form( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel &K,
                              ConstMatrixView Se,
                              ConstMatrixView Sx,
                              Numeric tol,
                              Index max_iter )
{
    Index n( x.nelem() ), m( y.nelem() );
    Numeric di2;
    Matrix Ki(m,n), KiTSeInvKi(n, n), KiTSeInv(m,n), SeInv(m,m), SxInv(n,n);
    Vector xi(n), dx(n), dx_old(n), yi(m), tm(m), tn(n), sInvDx(n);

    Index iter = 0;
    bool converged = false;

    // Required matrix inverses.
    inv(SeInv, Se);
    inv(SxInv, Sx);

    // Set the starting vector.
    x = xa;
    dx_old = 0;

    while ((!converged) && (iter < max_iter))
    {

        // Compute Jacobian and y_i.
        K.evaluate_jacobian( yi, Ki, x );
        mult( KiTSeInv, transpose( Ki ), SeInv );
        mult( KiTSeInvKi, KiTSeInv, Ki );
        KiTSeInvKi += SxInv;

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
        mult( tn, SxInv, dx_old );
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


