/* Copyright (C) 2002-2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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

#include "matpackI.h"
#include "matpackIII.h"

Numeric fac(const Index n);

Index integer_div( const Index& x, const Index& y );

Numeric last( ConstVectorView x );

Index last( const ArrayOfIndex& x );

void linspace(                      
              Vector&     x,           
              const Numeric  start,    
              const Numeric  stop,        
              const Numeric  step );

void nlinspace(         
              Vector&     x, 
              const Numeric     start,     
              const Numeric     stop,        
              const Index       n );

void nlinspace(         
              VectorView        x, 
              const Numeric     start,     
              const Numeric     stop,        
              const Index       n );

void nlogspace(         
               Vector&     x, 
               const Numeric     start,     
               const Numeric     stop,        
               const Index         n );

Numeric AngIntegrate_trapezoid(ConstMatrixView Integrand,
                               ConstVectorView za_grid,
                               ConstVectorView aa_grid);

Numeric AngIntegrate_trapezoid(ConstVectorView Integrand,
                               ConstVectorView za_grid);

Numeric AngIntegrate_trapezoid_opti(ConstMatrixView Integrand,
                                    ConstVectorView za_grid,
                                    ConstVectorView aa_grid,
                                    ConstVectorView grid_stepsize);

Numeric sign( const Numeric& x );

void mgd(
          VectorView  psd,
    const Vector&     x,
    const Numeric&    n0,
    const Numeric&    mu,
    const Numeric&    la,
    const Numeric&    ga );

void mgd_with_derivatives(
          VectorView  psd,
          MatrixView  jac_data,
    const Vector&     x,
    const Numeric&    n0,
    const Numeric&    mu,
    const Numeric&    la,
    const Numeric&    ga,
    const bool&       do_n0_jac,
    const bool&       do_mu_jac,
    const bool&       do_la_jac,
    const bool&       do_ga_jac );

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
void delanoe_shape_with_derivative(
    VectorView  psd,
    MatrixView  jac_data,
    const Vector&     x,
    const Numeric&    alpha,
    const Numeric&    beta);

Numeric mod_gamma_dist(Numeric x,
                       Numeric N0,
                       Numeric Lambda,
                       Numeric mu,
                       Numeric gamma);

void unitl( Vector& x );

void flat( VectorView x, ConstMatrixView X );
void flat( VectorView x, ConstTensor3View X );

void reshape( MatrixView X, ConstVectorView x );
void reshape( Tensor3View X, ConstVectorView x );

#endif  // math_funcs_h
