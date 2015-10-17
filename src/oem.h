/*!
  \file   oem.h
  \author simon <simonpf@chalmers.se>
  \date   Fri Apr 17 16:17:54 2015

  \brief Optimal estimation method for retrieval.
*/

#ifndef oem_h
#define oem_h

#include "logic.h"
#include "matpackI.h"

//! The Forward Model Class
/*!

  Abstract class to provide a communication interface between
  non-linear oem methods and the forward model.

*/
class ForwardModel
{
public:
//! Return a linearization, evaluate the forward model at the given point xi
//  and write the results into Ki and yi, respectively.

    virtual void evaluate_jacobian (  VectorView &yi,
                                      MatrixView &Ki,
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

    // Compute optimal estimator, simple method.
    int compute( Vector &x,
                 ConstVectorView y,
                 ConstVectorView y0 );

    // Compute optimal estimator using the Gain matrix.
    int compute( Vector &x,
                 MatrixView G,
                 ConstVectorView y,
                 ConstVectorView y0 );

    // Compute the fit of a given estimator.
    int compute_fit( Vector &yf,
                     const Vector &x,
                     ForwardModel &F );

    // Compute fit and cost of a given estimator.
    int compute_fit( Vector &yf,
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

    // Internal matrices and vectors needed for the computations.
    Matrix tmp_nn_1, tmp_nn_2, tmp_nm_1, tmp_mn_1, LU;
    Vector xa, tmp_m_1, tmp_n_1, x_norm;

    ArrayOfIndex indx;
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
                        Index& used_iter,
                        ForwardModel &F,
                        ConstVectorView xa,
                        ConstVectorView x_norm,
                        ConstVectorView y,
                        ConstMatrixView SeInv,
                        ConstMatrixView SxInv,
                        const Numeric& cost_start,
                        const Index maxIter,
                        const Numeric tol,
                        const bool& verbose );

bool oem_levenberg_marquardt( Vector& x,
                              Vector& yf,
                              Matrix& G,
                              Matrix& J,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ConstMatrixView SeInv,
                              ConstMatrixView SxInv,
                              ForwardModel &K,
                              Numeric tol,
                              Index maxIter,
                              Numeric gamma_start,
                              Numeric gamma_scale_dec,
                              Numeric gamma_scale_inc,
                              Numeric gamma_max,
                              Numeric gamma_threshold,
                              bool verbose );

#endif // oem_h
