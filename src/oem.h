/*!
  \file   oem.h
  \author simon <simonpf@chalmers.se>
  \date   Fri Apr 17 16:17:54 2015

  \brief Optimal estimation method for retrieval.
*/

#ifndef oem_h
#define oem_h

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

// Optimal estimation method for linear models.
void oem_linear_nform( Vector& x,
                       Matrix& G,
                       ConstVectorView xa,
                       ConstVectorView yf,
                       ConstVectorView y,
                       ConstMatrixView K,
                       ConstMatrixView SeInv,
                       ConstMatrixView SaInv );

// Optimal estimation method for linear models.
void oem_linear_mform( VectorView x,
                       MatrixView G,
                       ConstVectorView xa,
                       ConstVectorView yf,
                       ConstVectorView y,
                       ConstMatrixView K,
                       ConstMatrixView Se,
                       ConstMatrixView Sa );

// Optimal estimation for non-linear models using Gauss-Newton method.
bool oem_gauss_newton( Vector& x,
                       Vector& yf,
                       Matrix& G,
                       Matrix& J,
                       ConstVectorView y,
                       ConstVectorView xa,
                       ConstMatrixView SeInv,
                       ConstMatrixView SaInv,
                       ForwardModel& K,
                       Numeric tol,
                       Index maxiter,
                       bool verbose );

// Optimal estimation for non-linear models using Gauss-Newton method.
bool oem_gauss_newton_n_form( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel& K,
                              ConstMatrixView Se,
                              ConstMatrixView Sa,
                              Numeric tol,
                              Index maxiter );

// Optimal estimation for non-linear models using Gauss-Newton method.
bool oem_gauss_newton_m_form( VectorView x,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ForwardModel& K,
                              ConstMatrixView Se,
                              ConstMatrixView Sa,
                              Numeric tol,
                              Index maxiter );

// Optimal estimation for non-linear models using Levenberg-Marquardt method.
bool oem_levenberg_marquardt( Vector &x,
                              Vector &yf,
                              Matrix &G,
                              Matrix &J,
                              ConstVectorView y,
                              ConstVectorView xa,
                              ConstMatrixView SeInv,
                              ConstMatrixView SaInv,
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
