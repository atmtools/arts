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


void oem_cost_y( Numeric& cost_y, 
                 ConstVectorView y, 
                 ConstVectorView yf,
                 ConstMatrixView SeInv,
                 const Numeric&  normfac );

// Optimal estimation method for linear models.
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

// Optimal estimation for non-linear models using Levenberg-Marquardt method.
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
