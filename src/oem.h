/*!
  \file   oem.h
  \author Simon Pfreundschuh <simonpf@chalmers.se>
  \date   Fri Apr 17 16:17:54 2015

  \brief Optimal estimation method for retrieval.
*/

#ifndef oem_h
#define oem_h

#include "logic.h"
#include "matpackI.h"
#include "timings.h"
#include <string>
#include <vector>

//! The Forward Model Class
/*!

  Abstract class to provide a communication interface between
  non-linear oem methods and the forward model.

*/
class ForwardModel
{
public:
//! Return a linearization, evaluate the forward model at the given point xi
//  and write the results into J and yi, respectively.

    virtual void evaluate_jacobian (  VectorView &yi,
                                      MatrixView &J,
                                      const ConstVectorView &xi ) = 0;
    virtual void evaluate ( VectorView &yi,
                            const ConstVectorView &xi ) = 0;
};

//! OEM Form Enum
/*!

  Enum representing the formulation of the OEM equations. The equations
  for the computation of the optimal estimator can be formulated in two ways.
  When the n-form is used, a linear system of size (n,n) has to be solved. When
  the m-form is used the linear system to be solved has size (m,m). The choice
  should therefore be made depending on which of the parameters m,n is smaller.

 */
enum OEMForm { NFORM, MFORM };

//! Linear OEM Class
/*!
  Class that represents a linear OEM computation using the n-form for the
  computation of the optimal estimator. Given a linear forward model
  model described by a Jacobian J, an inverse measurement covariance Matrix
  SeInv and an inverse a priori covariance Matrix SxInv the linear OEM class
  can be used to compute the optimal estimator and the gain matrix. The object
  stores intermediate results which can be used to considerable subsequent
  computations.
 */
class LinearOEM
{
public:

    // Constructor
    LinearOEM( ConstMatrixView J,
               ConstMatrixView SeInv,
               ConstVectorView xa,
               ConstMatrixView SxInv );

    // Get and set normalization vector.
    void set_x_norm( ConstVectorView );
    ConstVectorView get_x_norm();

   //! Return gain matrix.
   /*!
    If not already computed, compute the gain matrix and return a matrix view of
    it.

     \return ConstMatrixView of the gain matrix of the linear model.
    */
    ConstMatrixView get_G()
    {
        if (gain_set)
            return G;
        else
        {
            compute_gain_matrix();
            return G;
        }
    }

    //! Reset intermediate results.
    /*!
     Resets the internally stored result and forces their recomputation on the
     next call of compute(...).
    */
    void reset()
    {
        matrices_set = false;
        gain_set = false;
        x_norm_set = false;
    }

    //! Get error code.
    /*!
      \return The internal error code.
    */
    Index get_error()
    {
        return err;
    }

    // Compute optimal estimator, simple method.
    Index compute( Vector &x,
                   ConstVectorView y,
                   ConstVectorView y0 );

    // Compute optimal estimator using the Gain matrix.
    Index compute( Vector &x,
                   MatrixView G,
                   ConstVectorView y,
                   ConstVectorView y0 );

    // Compute the fit of a given estimator.
    Index compute_fit( Vector &yf,
                       const Vector &x,
                       ForwardModel &F );

    // Compute fit and cost of a given estimator.
    Index compute_fit( Vector &yf,
                       Numeric &cost_x,
                       Numeric &cost_y,
                       Vector &x,
                       ConstVectorView y,
                       ForwardModel &F );
private:

    // Hide default constructor.
    LinearOEM();

    void compute_gain_matrix();

    // Model parameters.
    Index n,m;

    bool matrices_set, gain_set, x_norm_set;
    OEMForm form;

    ConstMatrixView J, SeInv, SxInv;
    Matrix G;

    Timings timer;

    // Internal matrices and vectors needed for the computations.
    Matrix tmp_nn_1, tmp_nn_2, tmp_nm_1, tmp_mn_1, LU;
    Vector xa, tmp_m_1, tmp_n_1, x_norm;
    Index err;

    ArrayOfIndex indx;
};

//! OEM Method Enum
/*!
  OEM method enum representing the two methods for solving non-linear inverse
  problems, Gauss-Newton and Levenberg-Marquardt.

*/
enum OEMMethod { GAUSS_NEWTON, LEVENBERG_MARQUARDT  };

//! Non-Linear OEM Class
/*!

  Class to represent non-linear OEM computations. Given a forward model
  described by the inverses of the measurement and state covariance matrices
  SeInv and SxInv, the a priori vector xa and a ForwardModel instance F, the
  class can computes the Bayesian optimal estimator using either the
  Gauss-Newton method or Levenber-Marquardt method, as described in Rodgers
  book. The form used is the n-form given in formula (5.8).

  The NonLinearOEM object contains references to the matrices, vectors and the
  forward model defining the problem, internal matrices, that are required
  during the computation and the computation state.

*/
class NonLinearOEM
{

public:

    // Constructor
    NonLinearOEM( ConstMatrixView SeInv,
                  ConstVectorView xa,
                  ConstMatrixView SxInv,
                  ForwardModel &F,
                  OEMMethod method );

    // Get and set normalization vector.
    void set_x_norm( ConstVectorView );
    ConstVectorView get_x_norm();

   //! Set iteration maximum.

   /*! Sets the number of iterations that are performed before the computation
     is aborted. Default is 100.

     \param max The new maximum number of iterations.
   */
    void maximum_iterations( Index max )
    {
        max_iter = max;
    }

    //! Get iteration maximum.
    /*!
      \return The currently set iteration maximum.
    */
    Index maximum_iterations( )
    {
        return max_iter;
    }

    //! Set convergence criterion.
    /*!
      Convergence is determined using equation (5.29) in Rodger's book. Note
      that the provided tolerance value is scaled by n before comparing to
      d_i^2.

      \param tol_ The new convergence criterion.
    */
    void tolerance( Numeric tol_)
    {
        tol = tol_;
    }

    //! Get current convergence criterion.
    /*!
      \return The current convergence criterion.
    */
    Numeric tolerance( )
    {
        return tol;
    }

    //! Get number of iterations.
    /*!
       Returns the number of iterations performed in the latest computation
       before abortion or the convergence criterion was met.

      \return Number of iterations of the latest computation.
    */
    Index iterations()
    {
        return iter;
    }

    //! Set gamma start.
    /*!
      The set initial gamma value.

      Note: Only affects Levenberg-Marquardt method.

      \param[in] ga The new gamma_start value.
    */
    void gamma_start( Numeric ga )
    {
        ga_start = ga;
    }

    //! Set gamma scaling factor.
    /*!
      Set the factor that is used to INCREASE gamma if no new x-value can be
      found that decreases the cost function.

      Note: Only affects Levenberg-Marquardt method.

      \param[in] ga The new gamma_increase value.
    */
    void gamma_increase( Numeric ga )
    {
        ga_increase = ga;
    }

    //! Set gamma scaling factor.
    /*!
      Set the factor that is used to DECREASE gamma if a new x-value was found
      found that decreases the cost function.

      Note: Only affects Levenberg-Marquardt method.

      \param[in] ga The new gamma_decrease value.
    */
    void gamma_decrease( Numeric ga )
    {
        ga_decrease = ga;
    }

    //! Set maximum gamma value.
    /*!
      Gamma is increased until it reaches gamma_max. If no new x-value can be
      found with a current gamma of gamma_max, the iteration aborts.

      Note: Only affects Levenberg-Marquardt method.

      \param[in] ga The new gamma_max value.
    */
    void gamma_max( Numeric ga )
    {
        ga_max = ga;
    }

    //! Set gamma threshold value.
    /*!
      If gamma is decreased until it reaches the threshold, it is set to zero.
      If the current gamma is zero and no new x-value that decreases the cost
      function can be found gamma is set to gamma_thresh.

      Note: Only affects Levenberg-Marquardt method.

      \param[in] ga The new gamma_threshold value.
    */
    void gamma_threshold( Numeric ga )
    {
        ga_threshold = ga;
    }

    //! Get error code.
    /*!
      \return The internal error code.
    */
    Index get_error()
    {
        return err;
    }

    //! Get Jacobian.
    /*!
      \return ConstMatrixView of the Jacobian.
    */
    ConstMatrixView get_jacobian()
    {
        return J;
    }

    // Perform OEM calculation.
    Index compute( Vector &x,
                   ConstVectorView y,
                   bool verbose );

    // Perform OEM calculation and compute gain matrix.
    Index compute( Vector &x,
                   MatrixView G_,
                   ConstVectorView y,
                   bool verbose );

    // Comute the fitted measurement vector.
    Index compute_fit( Vector &yf,
                       const Vector &x );

    // Compute fitted measurement vector and evaluate cost function.
    Index compute_fit( Vector &yf,
                       Numeric &cost_x,
                       Numeric &cost_y,
                       const Vector &x,
                       ConstVectorView y );

   //! Return convergence status.
   /*!
     \return The convergence status
   */
    bool converged()
    {
        return conv;
    }

private:

    // Hide standard constructor.
    NonLinearOEM();

    // Internal member functions.
    void compute_gain_matrix( Vector& x );
    void gauss_newton( Vector &x,
                       ConstVectorView y,
                       bool verbose );
    void levenberg_marquardt( Vector &x,
                              ConstVectorView y,
                              bool verbose );

    // References to model data.
    ConstMatrixView SeInv, SxInv;
    ConstVectorView xa;
    ForwardModel &F;

    // Internal state variables.
    OEMMethod method;
    bool matrices_set, gain_set, x_norm_set, conv;
    Index m, n, iter, max_iter, err;
    Numeric ga_max, ga_start, ga_threshold, ga_decrease, ga_increase;

    Numeric tol, cost_x, cost_y;

    Timings timer;

    // Internal matrices for intermediate results.
    Matrix G, J, tmp_nm_1, tmp_nn_1, tmp_nn_2;
    Vector dx, yi, xnew, x_norm, tmp_m_1, tmp_m_2, tmp_n_1, tmp_n_2;

};


void oem_cost_y( Numeric& cost_y,
                 ConstVectorView y,
                 ConstVectorView yf,
                 ConstMatrixView SeInv,
                 const Numeric&  normfac );

// Optimal estimation method for linear models, n-form.
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
                        const bool& verbose );

// Optimal estimation method for linear models, n-form.
Index oem_linear_mform( Vector& x,
                        Matrix& G,
                        ConstVectorView xa,
                        ConstVectorView y,
                        ConstVectorView yf,
                        ConstMatrixView J,
                        ConstMatrixView SeInv,
                        ConstMatrixView SxInv );

// Optimal estimation for non-linear models using Gauss-Newton method.
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
                        ConstMatrixView SeInv,
                        ConstMatrixView SxInv,
                        const Index max_iter,
                        const Numeric tol,
                        bool verbose );

Index oem_levenberg_marquardt( Vector &x,
                               Matrix &G,
                               Matrix &J,
                               Vector &yf,
                               Numeric &cost_y,
                               Numeric &cost_x,
                               Index &iter,
                               ForwardModel &F,
                               ConstVectorView xa,
                               ConstVectorView x_norm,
                               ConstVectorView y,
                               ConstMatrixView SeInv,
                               ConstMatrixView SxInv,
                               Index max_iter,
                               Numeric tol,
                               Numeric gamma_start,
                               Numeric gamma_decrease,
                               Numeric gamma_increase,
                               Numeric gamma_max,
                               Numeric gamma_threshold,
                               bool verbose );

#endif // oem_h
