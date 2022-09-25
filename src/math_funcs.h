/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
   \file   math_funcs.h
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   \date 2000-09-18 

   Contains the decleration of the functions in math_funcs.cc.
*/

#ifndef math_funcs_h
#define math_funcs_h

#include "array.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_math.h"

Numeric fac(const Index n);

Index integer_div(const Index& x, const Index& y);

Numeric last(ConstVectorView x);

Index last(const ArrayOfIndex& x);

void linspace(Vector& x,
              const Numeric start,
              const Numeric stop,
              const Numeric step);

void nlinspace(Vector& x,
               const Numeric start,
               const Numeric stop,
               const Index n);

void nlinspace(VectorView x,
               const Numeric start,
               const Numeric stop,
               const Index n);

void nlogspace(Vector& x,
               const Numeric start,
               const Numeric stop,
               const Index n);

Numeric trapz(ConstVectorView x,
              ConstVectorView y);

void cumsum(VectorView csum,
            ConstVectorView x);

Numeric AngIntegrate_trapezoid(ConstMatrixView Integrand,
                               ConstVectorView za_grid,
                               ConstVectorView aa_grid);

Numeric AngIntegrate_trapezoid(ConstVectorView Integrand,
                               ConstVectorView za_grid);

Numeric AngIntegrate_trapezoid_opti(ConstMatrixView Integrand,
                                    ConstVectorView za_grid,
                                    ConstVectorView aa_grid,
                                    ConstVectorView grid_stepsize);

Numeric sign(const Numeric& x);

//! min_geq
/*! 
    Compares two values and returns the smallest positive >= *limit. That is,
    values below *limit* are ignored. If none of the values is >= *limit, 
    *limit-1 is returned.

    \param n1     Numeric value 1
    \param n2     Numeric value 2
    \param limit  Limit for considering n1 and n2

    \return  The smallest above or limit-1

    \author Patrick Eriksson
    \date   2020-08-12
*/
Numeric min_geq(const Numeric n1,
                const Numeric n2,
                const Numeric limit);

void mgd(VectorView psd,
         const Vector& x,
         const Numeric& n0,
         const Numeric& mu,
         const Numeric& la,
         const Numeric& ga);

void mgd_with_derivatives(VectorView psd,
                          MatrixView jac_data,
                          const Vector& x,
                          const Numeric& n0,
                          const Numeric& mu,
                          const Numeric& la,
                          const Numeric& ga,
                          const bool& do_n0_jac,
                          const bool& do_mu_jac,
                          const bool& do_la_jac,
                          const bool& do_ga_jac);

/**! Shape functions for normalized PSD.
 *
 * This function implements the shape function F(X, alpha, beta) from
 * as proposed by Delanoe et al. in "Normalized particle size distribution
 * for remote sensing application".
 *
 * @param[OUT] psd On return contains the values of F corresponding to
 *                 the values in x.
 * @param[OUT] jac_data On return contains the first derivative of F w.r.t
 *                      x evaluated at the values in x.
 * @param[IN] x The values at which to evaluate the shape functions and
 *              and derivatives.
 * @param[IN] alpha The alpha parameter of the shape function.
 * @param[IN] beta  The beta parameter of the shape function.
 */
void delanoe_shape_with_derivative(VectorView psd,
                                   MatrixView jac_data,
                                   const Vector& x,
                                   const Numeric& alpha,
                                   const Numeric& beta);

Numeric mod_gamma_dist(
    Numeric x, Numeric N0, Numeric Lambda, Numeric mu, Numeric gamma);

void unitl(Vector& x);

void flat(VectorView x, ConstMatrixView X);
void flat(VectorView x, ConstTensor3View X);

void reshape(MatrixView X, ConstVectorView x);
void reshape(Tensor3View X, ConstVectorView x);

void calculate_weights_linear(Vector& x, Vector& w, const Index nph);

/** Calculates trapezoidal integration weights for arbitray grid
 *
 * @param w Vector integration weights
 * @param x Vector evaluation points
 */
void calculate_int_weights_arbitrary_grid(Vector& w, const Vector& x);


/** Checks for negative values */
template <typename MatpackType>
constexpr bool any_negative(const MatpackType& var) noexcept {
  if (var.empty()) return false;
  if (min(var) < 0) return true;
  return false;
}

/** Computes std::pow(-1, x) without std::pow
 * 
 * @param x Index
 * @return constexpr Index 
 */
constexpr Index pow_negative_one(Index x) noexcept { return (x % 2) ? -1 : 1; }


#endif  // math_funcs_h
