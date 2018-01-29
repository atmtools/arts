/* Copyright (C) 2017
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * \file   propagationmatrix.h
 * \brief  Stuff related to the propagation matrix.
 * 
 * The reason is that the naive approach to keep the full propagation matrix is memory intensive and slow
 * 
 * \author Richard Larsson
 * \date   2017-06-23
 */

#include "propagationmatrix.h"
#include "complex.h"
#include "arts_omp.h"
#include "lin_alg.h"


void compute_transmission_matrix(Tensor3View T,
                                 const Numeric& r,
                                 const PropagationMatrix& upper_level,
                                 const PropagationMatrix& lower_level,
                                 const Index iz,
                                 const Index ia)
{
  const Index mstokes_dim = upper_level.StokesDimensions();
  const Index mfreqs = upper_level.NumberOfFrequencies();
  assert(r > 0.);
  assert(T.npages() == mfreqs);
  assert(T.nrows() == mstokes_dim);
  assert(T.ncols() == mstokes_dim);
  assert(lower_level.StokesDimensions() == mstokes_dim);
  assert(lower_level.NumberOfFrequencies() == mfreqs);
  
  if(mstokes_dim == 1)
  {
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      T(i, 0, 0) = exp(-0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]));
    }
  }
  else if(mstokes_dim == 2)
  {
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      MatrixView F = T(i, joker, joker);
      
      const Numeric a = -0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (upper_level.K12(iz, ia)[i] + lower_level.K12(iz, ia)[i]);
                    
      const Numeric exp_a = exp(a);
      
      if(b == 0.)
      {
        F(0, 1) = F(1, 0) = 0.;
        F(0, 0) = F(1, 1) = exp_a;
        continue;
      }
      
      const Numeric C0 = (b*cosh(b) - a*sinh(b))/b;
      const Numeric C1 = sinh(b)/b;  
      
      F(0, 0) = F(1, 1) = C0 + C1 * a;
      F(0, 1) = F(1, 0) = C1 * b;
      
      F *= exp_a;
    }
  }
  else if(mstokes_dim == 3)
  {
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      MatrixView F = T(i, joker, joker);
      
      const Numeric a = -0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (upper_level.K12(iz, ia)[i] + lower_level.K12(iz, ia)[i]), 
                    c = -0.5 * r * (upper_level.K13(iz, ia)[i] + lower_level.K13(iz, ia)[i]), 
                    u = -0.5 * r * (upper_level.K23(iz, ia)[i] + lower_level.K23(iz, ia)[i]);
                    
      const Numeric exp_a = exp(a);
      
      if(b == 0. and c == 0. and u == 0.)
      {
        F = 0.;
        F(0, 0) = F(1, 1) = F(2, 2) = exp_a;
        continue;
      }
      
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      
      const Numeric x = sqrt(b2 + c2 - u2), x2 = x * x, inv_x2 = 1.0/x2;
      const Numeric sinh_x = sinh(x), cosh_x = cosh(x);
      
      const Numeric C0 = (a2 * (cosh_x - 1) - a * x * sinh_x + x2) * inv_x2; // approaches (1-a)*exp_a for low x
      const Numeric C1 = (2 * a * (1 - cosh_x) + x * sinh_x) * inv_x2;  // approaches (exp_a) for low_x
      const Numeric C2 = (cosh_x - 1) * inv_x2; // Approaches infinity for low x
      
      F(0, 0) = F(1, 1) = F(2, 2) = C0 + C1 * a;
      F(0, 0) += C2 * (a2 + b2 + c2);
      F(1, 1) += C2 * (a2 + b2 - u2);
      F(2, 2) += C2 * (a2 + c2 - u2);
      
      F(0, 1) = F(1, 0) = C1 * b;
      F(0, 1) += C2 * (2*a*b - c*u);
      F(1, 0) += C2 * (2*a*b + c*u);
      
      F(0, 2) = F(2, 0) = C1 * c;
      F(0, 2) += C2 * (2*a*c + b*u);
      F(2, 0) += C2 * (2*a*c - b*u);
      
      F(1, 2) =  C1 * u + C2 * (2*a*u + b*c);
      F(2, 1) = -C1 * u - C2 * (2*a*u - b*c);
      
      F *= exp_a;
    }
  }
  else if(mstokes_dim == 4)
  {
    static const Numeric sqrt_05 = sqrt(0.5);
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      MatrixView F = T(i, joker, joker);
      
      const Numeric a = -0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (upper_level.K12(iz, ia)[i] + lower_level.K12(iz, ia)[i]), 
                    c = -0.5 * r * (upper_level.K13(iz, ia)[i] + lower_level.K13(iz, ia)[i]), 
                    d = -0.5 * r * (upper_level.K14(iz, ia)[i] + lower_level.K14(iz, ia)[i]), 
                    u = -0.5 * r * (upper_level.K23(iz, ia)[i] + lower_level.K23(iz, ia)[i]), 
                    v = -0.5 * r * (upper_level.K24(iz, ia)[i] + lower_level.K24(iz, ia)[i]), 
                    w = -0.5 * r * (upper_level.K34(iz, ia)[i] + lower_level.K34(iz, ia)[i]);
                    
      const Numeric exp_a = exp(a);
      
      if(b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.)
      {
        F = 0.;
        F(0, 0) = F(1, 1) = F(2, 2) = F(3, 3) = exp_a;
        continue;
      }
      
      const Numeric b2 = b * b, c2 = c * c,
                    d2 = d * d, u2 = u * u,
                    v2 = v * v, w2 = w * w;
      
      const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
      
      Numeric Const1;
      Const1  = b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2)
              + c2 * (c2 * 0.5 +      d2 - u2 + v2 - w2)
              + d2 * (d2 * 0.5 +           u2 - v2 - w2)
              + u2 * (u2 * 0.5 +                v2 + w2)
              + v2 * (v2 * 0.5 +                     w2)
              + 4 * (b * d * u * w - b * c * v * w - c * d * u * v);
      Const1 *= 2;
      Const1 += w2 * w2;
      
      if(Const1 > 0.0)
        Const1 = sqrt(Const1);
      else
        Const1 = 0.0;
      
      const Complex sqrt_BpA = sqrt(Complex(Const2 + Const1, 0.0));
      const Complex sqrt_BmA = sqrt(Complex(Const2 - Const1, 0.0));
      const Numeric x = sqrt_BpA.real() * sqrt_05;
      const Numeric y = sqrt_BmA.imag() * sqrt_05;
      const Numeric x2 = x * x;
      const Numeric y2 = y * y;
      const Numeric cos_y = cos(y);
      const Numeric sin_y = sin(y);
      const Numeric cosh_x = cosh(x);
      const Numeric sinh_x = sinh(x);
      const Numeric x2y2 = x2 + y2;
      const Numeric inv_x2y2 = 1.0 / x2y2;
      
      Numeric C0, C1, C2, C3;
      Numeric inv_y = 0.0, inv_x = 0.0;  // Init'd to remove warnings
      
      // X and Y cannot both be zero
      if(x == 0.0)
      {
        inv_y = 1.0 / y;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (1.0 - cos_y) * inv_x2y2;
        C3 = (1.0 - sin_y*inv_y) * inv_x2y2;
      }
      else if(y == 0.0)
      {
        inv_x = 1.0 / x;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (cosh_x - 1.0) * inv_x2y2;
        C3 = (sinh_x*inv_x - 1.0) * inv_x2y2;
      }
      else
      {
        inv_x = 1.0 / x;
        inv_y = 1.0 / y;
        
        C0 = (cos_y*x2 + cosh_x*y2) * inv_x2y2;
        C1 = (sin_y*x2*inv_y + sinh_x*y2*inv_x) * inv_x2y2;
        C2 = (cosh_x - cos_y) * inv_x2y2;
        C3 = (sinh_x*inv_x - sin_y*inv_y) * inv_x2y2;
      }
      
      // Diagonal Elements
      F(0, 0) = F(1, 1) = F(2, 2) = F(3, 3) = C0;
      F(0, 0) += C2 * (b2 + c2 + d2);
      F(1, 1) += C2 * (b2 - u2 - v2);
      F(2, 2) += C2 * (c2 - u2 - w2);
      F(3, 3) += C2 * (d2 - v2 - w2);
      
      // Linear main-axis polarization
      F(0, 1) = F(1, 0) = C1 * b;
      F(0, 1) += C2 * (-c *  u -  d *  v) + C3 * ( b * ( b2 + c2 + d2) - u * ( b *  u -  d *  w) - v * ( b *  v +  c *  w));
      F(1, 0) += C2 * ( c * u + d * v) + C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) + d * (b * d + u * w));
      
      // Linear off-axis polarization
      F(0, 2) = F(2, 0) = C1 * c;
      F(0, 2) += C2 * ( b * u - d * w) + C3 * (c * (b2 + c2 + d2)  - u * (c * u + d * v) - w * (b * v + c * w));
      F(2, 0) += C2 * (-b * u + d * w) + C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) + d * (c * d - u * v));
      
      // Circular polarization
      F(0, 3) = F(3, 0) = C1 * d;
      F(0, 3) += C2 * ( b * v + c * w) + C3 * (d * (b2 + c2 + d2)  - v * (c * u + d * v) + w * (b * u - d * w));
      F(3, 0) += C2 * (-b * v - c * w) + C3 * (b * (b * d + u * w) + c * (c * d - u * v) - d * (-d2 + v2 + w2));
      
      // Circular polarization rotation
      F(1, 2) = F(2, 1) = C2 * (b * c - v * w);
      F(1, 2) +=  C1 * u + C3 * ( c * (c * u + d * v) - u * (-b2 + u2 + v2) - w * (b * d + u * w));
      F(2, 1) += -C1 * u + C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) - v * (c * d - u * v));
      
      // Linear off-axis polarization rotation
      F(1, 3) = F(3, 1) = C2 * (b * d + u * w);
      F(1, 3) +=  C1 * v + C3 * ( d * (c * u + d * v) - v * (-b2 + u2 + v2) + w * (b * c - v * w));
      F(3, 1) += -C1 * v + C3 * (-b * (b * v + c * w) - u * (c * d - u * v) + v * (-d2 + v2 + w2));
      
      // Linear main-axis polarization rotation
      F(2, 3) = F(3, 2) = C2 * (c * d - u * v);
      F(2, 3) +=  C1 * w + C3 * (-d * (b * u - d * w) + v * (b * c - v * w) - w * (-c2 + u2 + w2));
      F(3, 2) += -C1 * w + C3 * (-c * (b * v + c * w) + u * (b * d + u * w) + w * (-d2 + v2 + w2));
      
      F *= exp_a;
    }
  }
}


void compute_transmission_matrix_from_averaged_matrix_at_frequency(MatrixView T, 
                                                                   const Numeric& r, 
                                                                   const PropagationMatrix& averaged_propagation_matrix,
                                                                   const Index iv,
                                                                   const Index iz,
                                                                   const Index ia)
{
  static const Numeric sqrt_05 = sqrt(0.5);
  const Index mstokes_dim = averaged_propagation_matrix.StokesDimensions();
  assert(r > 0.);
  assert(T.nrows() == mstokes_dim);
  assert(T.ncols() == mstokes_dim);
  assert(iv < averaged_propagation_matrix.NumberOfFrequencies());
  
  if(mstokes_dim == 1)
  {
    T(0, 0) = exp(-r * averaged_propagation_matrix.Kjj(iz, ia)[iv]);
  }
  else if(mstokes_dim == 2)
  { 
    const Numeric a = -r * averaged_propagation_matrix.Kjj(iz, ia)[iv], 
                  b = -r * averaged_propagation_matrix.K12(iz, ia)[iv];
                  
    const Numeric exp_a = exp(a);
    
    if(b == 0.)
    {
      T(0, 1) = T(1, 0) = 0.;
      T(0, 0) = T(1, 1) = exp_a;
    }
    else
    {
      const Numeric C0 = (b*cosh(b) - a*sinh(b))/b;
      const Numeric C1 = sinh(b)/b;  
      
      T(0, 0) = T(1, 1) = C0 + C1 * a;
      T(0, 1) = T(1, 0) = C1 * b;
      
      T *= exp_a;
    }
  }
  else if(mstokes_dim == 3)
  {
    const Numeric a = -r * averaged_propagation_matrix.Kjj(iz, ia)[iv], 
                  b = -r * averaged_propagation_matrix.K12(iz, ia)[iv], 
                  c = -r * averaged_propagation_matrix.K13(iz, ia)[iv], 
                  u = -r * averaged_propagation_matrix.K23(iz, ia)[iv];
                  
    const Numeric exp_a = exp(a);
    
    if(b == 0. and c == 0. and u == 0.)
    {
      T = 0.;
      T(0, 0) = T(1, 1) = T(2, 2) = exp_a;
    }
    else
    {
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      
      const Numeric x = sqrt(b2 + c2 - u2), x2 = x * x, inv_x2 = 1.0/x2;
      const Numeric sinh_x = sinh(x), cosh_x = cosh(x);
      
      const Numeric C0 = (a2 * (cosh_x - 1) - a * x * sinh_x + x2) * inv_x2; // approaches (1-a)*exp_a for low x
      const Numeric C1 = (2 * a * (1 - cosh_x) + x * sinh_x) * inv_x2;  // approaches (exp_a) for low_x
      const Numeric C2 = (cosh_x - 1) * inv_x2; // Approaches infinity for low x
      
      T(0, 0) = T(1, 1) = T(2, 2) = C0 + C1 * a;
      T(0, 0) += C2 * (a2 + b2 + c2);
      T(1, 1) += C2 * (a2 + b2 - u2);
      T(2, 2) += C2 * (a2 + c2 - u2);
      
      T(0, 1) = T(1, 0) = C1 * b;
      T(0, 1) += C2 * (2*a*b - c*u);
      T(1, 0) += C2 * (2*a*b + c*u);
      
      T(0, 2) = T(2, 0) = C1 * c;
      T(0, 2) += C2 * (2*a*c + b*u);
      T(2, 0) += C2 * (2*a*c - b*u);
      
      T(1, 2) =  C1 * u + C2 * (2*a*u + b*c);
      T(2, 1) = -C1 * u - C2 * (2*a*u - b*c);
      
      T *= exp_a;
    }
  }
  else if(mstokes_dim == 4)
  {
    const Numeric a = -r * averaged_propagation_matrix.Kjj(iz, ia)[iv], 
                  b = -r * averaged_propagation_matrix.K12(iz, ia)[iv], 
                  c = -r * averaged_propagation_matrix.K13(iz, ia)[iv], 
                  d = -r * averaged_propagation_matrix.K14(iz, ia)[iv], 
                  u = -r * averaged_propagation_matrix.K23(iz, ia)[iv], 
                  v = -r * averaged_propagation_matrix.K24(iz, ia)[iv], 
                  w = -r * averaged_propagation_matrix.K34(iz, ia)[iv];
                  
    const Numeric exp_a = exp(a);
    
    if(b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.)
    {
      T = 0.;
      T(0, 0) = T(1, 1) = T(2, 2) = T(3, 3) = exp_a;
    }
    else
    {
      const Numeric b2 = b * b, c2 = c * c,
                    d2 = d * d, u2 = u * u,
                    v2 = v * v, w2 = w * w;
      
      const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
      
      Numeric Const1;
      Const1  = b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2);
      Const1 += c2 * (c2 * 0.5 + d2 - u2 + v2 - w2);
      Const1 += d2 * (d2 * 0.5 + u2 - v2 - w2);
      Const1 += u2 * (u2 * 0.5 + v2 + w2);
      Const1 += v2 * (v2 * 0.5 + w2);
      Const1 *= 2;
      Const1 += 8 * (b * d * u * w - b * c * v * w - c * d * u * v);
      Const1 += w2 * w2;
      
      if(Const1 > 0.0)
        Const1 = sqrt(Const1);
      else
        Const1 = 0.0;
      
      const Complex sqrt_BpA = sqrt(Complex(Const2 + Const1, 0.0));
      const Complex sqrt_BmA = sqrt(Complex(Const2 - Const1, 0.0));
      const Numeric x = sqrt_BpA.real() * sqrt_05;
      const Numeric y = sqrt_BmA.imag() * sqrt_05;
      const Numeric x2 = x * x;
      const Numeric y2 = y * y;
      const Numeric cos_y = cos(y);
      const Numeric sin_y = sin(y);
      const Numeric cosh_x = cosh(x);
      const Numeric sinh_x = sinh(x);
      const Numeric x2y2 = x2 + y2;
      const Numeric inv_x2y2 = 1.0 / x2y2;
      
      Numeric C0, C1, C2, C3;
      Numeric inv_y = 0.0, inv_x = 0.0;  // Init'd to remove warnings
      
      // X and Y cannot both be zero
      if(x == 0.0)
      {
        inv_y = 1.0 / y;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (1.0 - cos_y) * inv_x2y2;
        C3 = (1.0 - sin_y*inv_y) * inv_x2y2;
      }
      else if(y == 0.0)
      {
        inv_x = 1.0 / x;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (cosh_x - 1.0) * inv_x2y2;
        C3 = (sinh_x*inv_x - 1.0) * inv_x2y2;
      }
      else
      {
        inv_x = 1.0 / x;
        inv_y = 1.0 / y;
        
        C0 = (cos_y*x2 + cosh_x*y2) * inv_x2y2;
        C1 = (sin_y*x2*inv_y + sinh_x*y2*inv_x) * inv_x2y2;
        C2 = (cosh_x - cos_y) * inv_x2y2;
        C3 = (sinh_x*inv_x - sin_y*inv_y) * inv_x2y2;
      }
      
      // Diagonal Elements
      T(0, 0) = T(1, 1) = T(2, 2) = T(3, 3) = C0;
      T(0, 0) += C2 * (b2 + c2 + d2);
      T(1, 1) += C2 * (b2 - u2 - v2);
      T(2, 2) += C2 * (c2 - u2 - w2);
      T(3, 3) += C2 * (d2 - v2 - w2);
      
      // Linear main-axis polarization
      T(0, 1) = T(1, 0) = C1 * b;
      T(0, 1) += C2 * (-c *  u -  d *  v) + C3 * ( b * ( b2 + c2 + d2) - u * ( b *  u -  d *  w) - v * ( b *  v +  c *  w));
      T(1, 0) += C2 * ( c * u + d * v) + C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) + d * (b * d + u * w));
      
      // Linear off-axis polarization
      T(0, 2) = T(2, 0) = C1 * c;
      T(0, 2) += C2 * ( b * u - d * w) + C3 * (c * (b2 + c2 + d2)  - u * (c * u + d * v) - w * (b * v + c * w));
      T(2, 0) += C2 * (-b * u + d * w) + C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) + d * (c * d - u * v));
      
      // Circular polarization
      T(0, 3) = T(3, 0) = C1 * d;
      T(0, 3) += C2 * ( b * v + c * w) + C3 * (d * (b2 + c2 + d2)  - v * (c * u + d * v) + w * (b * u - d * w));
      T(3, 0) += C2 * (-b * v - c * w) + C3 * (b * (b * d + u * w) + c * (c * d - u * v) - d * (-d2 + v2 + w2));
      
      // Circular polarization rotation
      T(1, 2) = T(2, 1) = C2 * (b * c - v * w);
      T(1, 2) +=  C1 * u + C3 * ( c * (c * u + d * v) - u * (-b2 + u2 + v2) - w * (b * d + u * w));
      T(2, 1) += -C1 * u + C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) - v * (c * d - u * v));
      
      // Linear off-axis polarization rotation
      T(1, 3) = T(3, 1) = C2 * (b * d + u * w);
      T(1, 3) +=  C1 * v + C3 * ( d * (c * u + d * v) - v * (-b2 + u2 + v2) + w * (b * c - v * w));
      T(3, 1) += -C1 * v + C3 * (-b * (b * v + c * w) - u * (c * d - u * v) + v * (-d2 + v2 + w2));
      
      // Linear main-axis polarization rotation
      T(2, 3) = T(3, 2) = C2 * (c * d - u * v);
      T(2, 3) +=  C1 * w + C3 * (-d * (b * u - d * w) + v * (b * c - v * w) - w * (-c2 + u2 + w2));
      T(3, 2) += -C1 * w + C3 * (-c * (b * v + c * w) + u * (b * d + u * w) + w * (-d2 + v2 + w2));
      
      T *= exp_a;
    }
  }
}


void compute_transmission_matrix_and_derivative(Tensor3View T, 
                                                Tensor4View dT_dx_upper_level, 
                                                Tensor4View dT_dx_lower_level, 
                                                const Numeric& r, 
                                                const PropagationMatrix& upper_level,
                                                const PropagationMatrix& lower_level,
                                                const ArrayOfPropagationMatrix& dupper_level_dx, 
                                                const ArrayOfPropagationMatrix& dlower_level_dx,
                                                const Index iz,
                                                const Index ia)
{
  
  const Index mstokes_dim = upper_level.StokesDimensions();
  const Index mfreqs = upper_level.NumberOfFrequencies();
  const Index nppd = dupper_level_dx.nelem();
  
  DEBUG_ONLY(
    assert(r > 0.);
    assert(T.npages() == mfreqs);
    assert(T.nrows() == mstokes_dim);
    assert(T.ncols() == mstokes_dim);
    assert(lower_level.StokesDimensions() == mstokes_dim);
    assert(lower_level.NumberOfFrequencies() == mfreqs);
    
    assert(nppd == dlower_level_dx.nelem());
    assert(nppd == dT_dx_upper_level.nbooks());
    assert(nppd == dT_dx_lower_level.nbooks());
    assert(mfreqs == dT_dx_upper_level.npages());
    assert(mfreqs == dT_dx_lower_level.npages());
    assert(mstokes_dim == dT_dx_upper_level.nrows());
    assert(mstokes_dim == dT_dx_lower_level.nrows());
    assert(mstokes_dim == dT_dx_upper_level.ncols());
    assert(mstokes_dim == dT_dx_lower_level.ncols());
    
    for(const PropagationMatrix& d_dx : dupper_level_dx)
    {
      if(d_dx.StokesDimensions() == 1 and d_dx.NumberOfFrequencies() == 0) // Catch cases with outer-perturbed jacobians
        continue;
      assert(d_dx.StokesDimensions() == mstokes_dim);
      assert(d_dx.NumberOfFrequencies() == mfreqs or d_dx.NumberOfFrequencies() == 0);
    }
    
    for(const PropagationMatrix& d_dx : dlower_level_dx)
    {
      if(d_dx.StokesDimensions() == 1 and d_dx.NumberOfFrequencies() == 0) // Catch cases with outer-perturbed jacobians
        continue;
      assert(d_dx.StokesDimensions() == mstokes_dim);
      assert(d_dx.NumberOfFrequencies() == mfreqs or d_dx.NumberOfFrequencies() == 0);
    }
  );
  
  if(mstokes_dim == 1)
  {
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      T(i, 0, 0) = exp(-0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]));
      for(Index j = 0; j < nppd; j++)
      {
        if(dupper_level_dx[j].NumberOfFrequencies())
          dT_dx_upper_level(j, i, 0, 0) = -0.5 * r * T(i, 0, 0) * dupper_level_dx[j].Kjj(iz, ia)[i];
        
        if(dlower_level_dx[j].NumberOfFrequencies())
          dT_dx_lower_level(j, i, 0, 0) = -0.5 * r * T(i, 0, 0) * dlower_level_dx[j].Kjj(iz, ia)[i];
      }
    }
  }
  else if(mstokes_dim == 2)
  {
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      MatrixView F = T(i, joker, joker);
      
      const Numeric a = -0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (upper_level.K12(iz, ia)[i] + lower_level.K12(iz, ia)[i]);
                    
      const Numeric exp_a = exp(a);
      
      if(b == 0.)
      {
        F(0, 1) = F(1, 0) = 0.;
        F(0, 0) = F(1, 1) = exp_a;
        for(Index j = 0; j < nppd; j++)
        {
          if(dupper_level_dx[j].NumberOfFrequencies())
          {
            dupper_level_dx[j].RightMultiplyAtPosition(dT_dx_upper_level(j, i, joker, joker), F, i, iz, ia);
            dT_dx_upper_level(j, i, joker, joker) *= -0.5 * r;
          }
          
          if(dlower_level_dx[j].NumberOfFrequencies())
          {
            dlower_level_dx[j].RightMultiplyAtPosition(dT_dx_lower_level(j, i, joker, joker), F, i ,iz, ia);
            dT_dx_lower_level(j, i, joker, joker) *= -0.5 * r;
          }
        }
        continue;
      }
      
      const Numeric C0 = cosh(b) - a*sinh(b)/b;
      const Numeric C1 = sinh(b)/b;  
      
      F(0, 0) = F(1, 1) = C0 + C1 * a;
      F(0, 1) = F(1, 0) = C1 * b;
      
      F *= exp_a;
      
      for(Index j = 0; j < nppd; j++)
      {
        if(not dlower_level_dx[j].NumberOfFrequencies())
          continue;
        MatrixView dF = dT_dx_lower_level(j, i, joker, joker);
        
        const Numeric da = -0.5 * r * dlower_level_dx[j].Kjj()[i],
                      db = -0.5 * r * dlower_level_dx[j].K12()[i];
        
        const Numeric dC0 = -a*cosh(b)*db/b + a*sinh(b)*db/b/b + sinh(b)*db - sinh(b)*da/b;
        const Numeric dC1 = (cosh(b) - C1)*db/b;
        
        dF(0, 0) = dF(1, 1) = dC0 + C1 * da + dC1 * a + F(0, 0) * da;
        dF(0, 1) = dF(1, 0) = C1 * db + dC1 * b + F(0, 1) * da;
      }
      
      for(Index j = 0; j < nppd; j++)
      {
        if(not dupper_level_dx[j].NumberOfFrequencies())
          continue;
        
        MatrixView dF = dT_dx_upper_level(j, i, joker, joker);
        
        const Numeric da = -0.5 * r * dupper_level_dx[j].Kjj(iz, ia)[i],
                      db = -0.5 * r * dupper_level_dx[j].K12(iz, ia)[i];
        
        const Numeric dC0 = -a*cosh(b)*db/b + a*sinh(b)*db/b/b + sinh(b)*db - sinh(b)*da/b;
        const Numeric dC1 = (cosh(b) - C1)*db/b;
        
        dF(0, 0) = dF(1, 1) = dC0 + C1 * da + dC1 * a + F(0, 0) * da;
        dF(0, 1) = dF(1, 0) = C1 * db + dC1 * b + F(0, 1) * da;
      }
    }
  }
  else if(mstokes_dim == 3)
  {
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      MatrixView F = T(i, joker, joker);
      
      const Numeric a = -0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (upper_level.K12(iz, ia)[i] + lower_level.K12(iz, ia)[i]), 
                    c = -0.5 * r * (upper_level.K13(iz, ia)[i] + lower_level.K13(iz, ia)[i]), 
                    u = -0.5 * r * (upper_level.K23(iz, ia)[i] + lower_level.K23(iz, ia)[i]);
                    
      const Numeric exp_a = exp(a);
      
      if(b == 0. and c == 0. and u == 0.)
      {
        F = 0.;
        F(0, 0) = F(1, 1) = F(2, 2) = exp_a;
        for(Index j = 0; j < nppd; j++)
        {
          if(dupper_level_dx[j].NumberOfFrequencies())
          {
            dupper_level_dx[j].RightMultiplyAtPosition(dT_dx_upper_level(j, i, joker, joker), F, i, iz, ia);
            dT_dx_upper_level(j, i, joker, joker) *= -0.5 * r;
          }
          
          if(dlower_level_dx[j].NumberOfFrequencies())
          {
            dlower_level_dx[j].RightMultiplyAtPosition(dT_dx_lower_level(j, i, joker, joker), F, i, iz, ia);
            dT_dx_lower_level(j, i, joker, joker) *= -0.5 * r;
          }
        }
        continue;
      }
      
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      
      const Numeric x = sqrt(b2 + c2 - u2), x2 = x * x, inv_x2 = 1.0/x2;
      const Numeric sinh_x = sinh(x), cosh_x = cosh(x);
      
      const Numeric C0 = (a2 * (cosh_x - 1) - a * x * sinh_x) * inv_x2 + 1; // approaches (1-a)*exp_a for low x
      const Numeric C1 = (2 * a * (1 - cosh_x) + x * sinh_x) * inv_x2;  // approaches (exp_a) for low_x
      const Numeric C2 = (cosh_x - 1) * inv_x2; // Approaches infinity for low x
      
      F(0, 0) = F(1, 1) = F(2, 2) = C0 + C1 * a;
      F(0, 0) += C2 * (a2 + b2 + c2);
      F(1, 1) += C2 * (a2 + b2 - u2);
      F(2, 2) += C2 * (a2 + c2 - u2);
      
      F(0, 1) = F(1, 0) = C1 * b;
      F(0, 1) += C2 * (2*a*b - c*u);
      F(1, 0) += C2 * (2*a*b + c*u);
      
      F(0, 2) = F(2, 0) = C1 * c;
      F(0, 2) += C2 * (2*a*c + b*u);
      F(2, 0) += C2 * (2*a*c - b*u);
      
      F(1, 2) =  C1 * u + C2 * (2*a*u + b*c);
      F(2, 1) = -C1 * u - C2 * (2*a*u - b*c);
      
      F *= exp_a;
      
      for(Index j = 0; j < nppd; j++)
      {
        if(not dlower_level_dx[j].NumberOfFrequencies())
          continue;
        
        MatrixView dF = dT_dx_lower_level(j, i, joker, joker);
        
        const Numeric da = -0.5 * r * dlower_level_dx[j].Kjj(iz, ia)[i],
                      db = -0.5 * r * dlower_level_dx[j].K12(iz, ia)[i],
                      dc = -0.5 * r * dlower_level_dx[j].K13(iz, ia)[i],
                      du = -0.5 * r * dlower_level_dx[j].K23(iz, ia)[i];
        
        const Numeric da2 = 2 * a * da, db2 = 2 * b * db, dc2 = 2 * c * dc, du2 = 2 * u * du;
        
        const Numeric dx  = (db2 + dc2 -du2)/x/2;
        
        const Numeric dC0 = -2 * (C0 - 1) * dx/x + (da2 * (cosh_x - 1) + a2 * sinh_x*dx - a * b * cosh_x * dx - a * sinh_x * dx - b * sinh_x * da) * inv_x2;
        const Numeric dC1 = -2 * C1 * dx / x + (2 * da * (1 - cosh_x) - 2 * a * sinh_x * dx + x*cosh_x*dx + sinh_x*dx) * inv_x2;
        const Numeric dC2 = (sinh_x/x - 2*C2)*dx/x;
        
        dF(0, 0) = dF(1, 1) = dF(2, 2) = dC0 + dC1 * a + C1 * da;
        dF(0, 0) += dC2 * (a2 + b2 + c2) + C2 * (da2 + db2 + dc2) + F(0, 0) * da;
        dF(1, 1) += dC2 * (a2 + b2 - u2) + C2 * (da2 + db2 - du2) + F(1, 1) * da;
        dF(2, 2) += dC2 * (a2 + c2 - u2) + C2 * (da2 + dc2 - du2) + F(2, 2) * da;
        
        dF(0, 1) = dF(1, 0) = dC1 * b + C1 * db;
        dF(0, 1) += dC2 * (2*a*b - c*u) + C2 * (2*da*b + 2*a*db - dc*u - c*du) + F(0, 1) * da;
        dF(1, 0) += dC2 * (2*a*b + c*u) + C2 * (2*da*b + 2*a*db + dc*u + c*du) + F(1, 0) * da;
        
        dF(0, 2) = dF(2, 0) = dC1 * c + C1 * dc;
        dF(0, 2) += dC2 * (2*a*c + b*u) + C2 * (2*da*c + 2*a*dc + db*u + b*du) + F(0, 2) * da;
        dF(2, 0) += dC2 * (2*a*c - b*u) + C2 * (2*da*c + 2*a*dc - db*u - b*du) + F(2, 0) * da;
        
        dF(1, 2) =  dC1 * u + C1 * du + dC2 * (2*a*u + b*c) + C2 * (2*da*u + 2*a*du + db*c + b*dc) + F(1, 2) * da;
        dF(2, 1) = -dC1 * u - C1 * du - dC2 * (2*a*u - b*c) - C2 * (2*da*u + 2*a*du - db*c - b*dc) + F(2, 1) * da;
      }
      
      for(Index j = 0; j < nppd; j++)
      {
        if(not dupper_level_dx[j].NumberOfFrequencies())
          continue;
        
        MatrixView dF = dT_dx_upper_level(j, i, joker, joker);
        
        const Numeric da = -0.5 * r * dupper_level_dx[j].Kjj(iz, ia)[i],
                      db = -0.5 * r * dupper_level_dx[j].K12(iz, ia)[i],
                      dc = -0.5 * r * dupper_level_dx[j].K13(iz, ia)[i],
                      du = -0.5 * r * dupper_level_dx[j].K23(iz, ia)[i];
        
        const Numeric da2 = 2 * a * da, db2 = 2 * b * db, dc2 = 2 * c * dc, du2 = 2 * u * du;
        
        const Numeric dx  = (db2 + dc2 -du2)/x/2;
        
        const Numeric dC0 = -2 * (C0 - 1) * dx/x + (da2 * (cosh_x - 1) + a2 * sinh_x*dx - a * b * cosh_x * dx - a * sinh_x * dx - b * sinh_x * da) * inv_x2;
        const Numeric dC1 = -2 * C1 * dx / x + (2 * da * (1 - cosh_x) - 2 * a * sinh_x * dx + x*cosh_x*dx + sinh_x*dx) * inv_x2;
        const Numeric dC2 = (sinh_x/x - 2*C2)*dx/x;
        
        dF(0, 0) = dF(1, 1) = dF(2, 2) = dC0 + dC1 * a + C1 * da;
        dF(0, 0) += dC2 * (a2 + b2 + c2) + C2 * (da2 + db2 + dc2) + F(0, 0) * da;
        dF(1, 1) += dC2 * (a2 + b2 - u2) + C2 * (da2 + db2 - du2) + F(1, 1) * da;
        dF(2, 2) += dC2 * (a2 + c2 - u2) + C2 * (da2 + dc2 - du2) + F(2, 2) * da;
        
        dF(0, 1) = dF(1, 0) = dC1 * b + C1 * db;
        dF(0, 1) += dC2 * (2*a*b - c*u) + C2 * (2*da*b + 2*a*db - dc*u - c*du) + F(0, 1) * da;
        dF(1, 0) += dC2 * (2*a*b + c*u) + C2 * (2*da*b + 2*a*db + dc*u + c*du) + F(1, 0) * da;
        
        dF(0, 2) = dF(2, 0) = dC1 * c + C1 * dc;
        dF(0, 2) += dC2 * (2*a*c + b*u) + C2 * (2*da*c + 2*a*dc + db*u + b*du) + F(0, 2) * da;
        dF(2, 0) += dC2 * (2*a*c - b*u) + C2 * (2*da*c + 2*a*dc - db*u - b*du) + F(2, 0) * da;
        
        dF(1, 2) =  dC1 * u + C1 * du + dC2 * (2*a*u + b*c) + C2 * (2*da*u + 2*a*du + db*c + b*dc) + F(1, 2) * da;
        dF(2, 1) = -dC1 * u - C1 * du - dC2 * (2*a*u - b*c) - C2 * (2*da*u + 2*a*du - db*c - b*dc) + F(2, 1) * da;
      }
    }
  }
  else if(mstokes_dim == 4)
  {
    static const Numeric sqrt_05 = sqrt(0.5);
    
    #pragma omp parallel for \
    if (!arts_omp_in_parallel() and mfreqs >= arts_omp_get_max_threads())
    for(Index i = 0; i < mfreqs; i++)
    {
      MatrixView F = T(i, joker, joker);
      
      
      const Numeric a = -0.5 * r * (upper_level.Kjj(iz, ia)[i] + lower_level.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (upper_level.K12(iz, ia)[i] + lower_level.K12(iz, ia)[i]), 
                    c = -0.5 * r * (upper_level.K13(iz, ia)[i] + lower_level.K13(iz, ia)[i]), 
                    d = -0.5 * r * (upper_level.K14(iz, ia)[i] + lower_level.K14(iz, ia)[i]), 
                    u = -0.5 * r * (upper_level.K23(iz, ia)[i] + lower_level.K23(iz, ia)[i]), 
                    v = -0.5 * r * (upper_level.K24(iz, ia)[i] + lower_level.K24(iz, ia)[i]), 
                    w = -0.5 * r * (upper_level.K34(iz, ia)[i] + lower_level.K34(iz, ia)[i]);
      
      const Numeric exp_a = exp(a);

      if(b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.)
      {
        F = 0.;
        F(0, 0) = F(1, 1) = F(2, 2) = F(3, 3) = exp_a;
        for(Index j = 0; j < nppd; j++)
        {
          if(dupper_level_dx[j].NumberOfFrequencies())
          {
            dupper_level_dx[j].RightMultiplyAtPosition(dT_dx_upper_level(j, i, joker, joker), F, i, iz, ia);
            dT_dx_upper_level(j, i, joker, joker) *= -0.5 * r;
          }
          
          if(dlower_level_dx[j].NumberOfFrequencies())
          {
            dlower_level_dx[j].RightMultiplyAtPosition(dT_dx_lower_level(j, i, joker, joker), F, i, iz, ia);
            dT_dx_lower_level(j, i, joker, joker) *= -0.5 * r;
          }
        }
        continue;
      }
      
      const Numeric b2 = b * b, c2 = c * c,
                    d2 = d * d, u2 = u * u,
                    v2 = v * v, w2 = w * w;
      
      const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
      
      Numeric Const1;
      Const1  = b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2);
      Const1 += c2 * (c2 * 0.5 + d2 - u2 + v2 - w2);
      Const1 += d2 * (d2 * 0.5 + u2 - v2 - w2);
      Const1 += u2 * (u2 * 0.5 + v2 + w2);
      Const1 += v2 * (v2 * 0.5 + w2);
      Const1 *= 2;
      Const1 += 8 * (b * d * u * w - b * c * v * w - c * d * u * v);
      Const1 += w2 * w2;
      
      if(Const1 > 0.0)
        Const1 = sqrt(Const1);
      else
        Const1 = 0.0;
      
      
      const Complex sqrt_BpA = sqrt(Complex(Const2 + Const1, 0.0));
      const Complex sqrt_BmA = sqrt(Complex(Const2 - Const1, 0.0));
      const Numeric x = sqrt_BpA.real() * sqrt_05;
      const Numeric y = sqrt_BmA.imag() * sqrt_05;
      const Numeric x2 = x * x;
      const Numeric y2 = y * y;
      const Numeric cos_y = cos(y);
      const Numeric sin_y = sin(y);
      const Numeric cosh_x = cosh(x);
      const Numeric sinh_x = sinh(x);
      const Numeric x2y2 = x2 + y2;
      const Numeric inv_x2y2 = 1.0 / x2y2;
      
      Numeric C0, C1, C2, C3;
      Numeric inv_y=0.0, inv_x=0.0;  // Init'd to remove warnings
      
      // X and Y cannot both be zero
      if(x == 0.0)
      {
        inv_y = 1.0 / y;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (1.0 - cos_y) * inv_x2y2;
        C3 = (1.0 - sin_y*inv_y) * inv_x2y2;
      }
      else if(y == 0.0)
      {
        inv_x = 1.0 / x;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (cosh_x - 1.0) * inv_x2y2;
        C3 = (sinh_x*inv_x -1.0) * inv_x2y2;
      }
      else
      {
        inv_x = 1.0 / x;
        inv_y = 1.0 / y;
        
        C0 = (cos_y*x2 + cosh_x*y2) * inv_x2y2;
        C1 = (sin_y*x2*inv_y + sinh_x*y2*inv_x) * inv_x2y2;
        C2 = (cosh_x - cos_y) * inv_x2y2;
        C3 = (sinh_x*inv_x - sin_y*inv_y) * inv_x2y2;
      }
      
      F(0, 0) = F(1, 1) = F(2, 2) = F(3, 3) = C0;
      F(0, 0) += C2 * (b2 + c2 + d2);
      F(1, 1) += C2 * (b2 - u2 - v2);
      F(2, 2) += C2 * (c2 - u2 - w2);
      F(3, 3) += C2 * (d2 - v2 - w2);
      
      F(0, 1) = F(1, 0) = C1 * b;
      F(0, 1) += C2 * (-c *  u -  d *  v) + C3 * ( b * ( b2 + c2 + d2) - u * ( b *  u -  d *  w) - v * ( b *  v +  c *  w));
      F(1, 0) += C2 * ( c * u + d * v) + C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) + d * (b * d + u * w));
      
      
      F(0, 2) = F(2, 0) = C1 * c;
      F(0, 2) += C2 * ( b * u - d * w) + C3 * (c * (b2 + c2 + d2)  - u * (c * u + d * v) - w * (b * v + c * w));
      F(2, 0) += C2 * (-b * u + d * w) + C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) + d * (c * d - u * v));
      
      
      F(0, 3) = F(3, 0) = C1 * d;
      F(0, 3) += C2 * ( b * v + c * w) + C3 * (d * (b2 + c2 + d2)  - v * (c * u + d * v) + w * (b * u - d * w));
      F(3, 0) += C2 * (-b * v - c * w) + C3 * (b * (b * d + u * w) + c * (c * d - u * v) - d * (-d2 + v2 + w2));
      
      
      F(1, 2) = F(2, 1) = C2 * (b * c - v * w);
      F(1, 2) +=  C1 * u + C3 * ( c * (c * u + d * v) - u * (-b2 + u2 + v2) - w * (b * d + u * w));
      F(2, 1) += -C1 * u + C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) - v * (c * d - u * v));
      
      
      F(1, 3) = F(3, 1) = C2 * (b * d + u * w);
      F(1, 3) +=  C1 * v + C3 * ( d * (c * u + d * v) - v * (-b2 + u2 + v2) + w * (b * c - v * w));
      F(3, 1) += -C1 * v + C3 * (-b * (b * v + c * w) - u * (c * d - u * v) + v * (-d2 + v2 + w2));
      
      
      F(2, 3) = F(3, 2) = C2 * (c * d - u * v);
      F(2, 3) +=  C1 * w + C3 * (-d * (b * u - d * w) + v * (b * c - v * w) - w * (-c2 + u2 + w2));
      F(3, 2) += -C1 * w + C3 * (-c * (b * v + c * w) + u * (b * d + u * w) + w * (-d2 + v2 + w2));
      
      F *= exp_a;
      
      if(nppd)
      {
        const Numeric inv_x2 = inv_x * inv_x;
        const Numeric inv_y2 = inv_y * inv_y;
        
        for(Index j = 0; j < nppd; j++)
        {
          if(not dupper_level_dx[j].NumberOfFrequencies())
            continue;
          
          MatrixView dF = dT_dx_upper_level(j, i, joker, joker);
          
          const Numeric da = -0.5 * r * dupper_level_dx[j].Kjj(iz, ia)[i],
                        db = -0.5 * r * dupper_level_dx[j].K12(iz, ia)[i],
                        dc = -0.5 * r * dupper_level_dx[j].K13(iz, ia)[i],
                        dd = -0.5 * r * dupper_level_dx[j].K14(iz, ia)[i],
                        du = -0.5 * r * dupper_level_dx[j].K23(iz, ia)[i],
                        dv = -0.5 * r * dupper_level_dx[j].K24(iz, ia)[i],
                        dw = -0.5 * r * dupper_level_dx[j].K34(iz, ia)[i];
          
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c,
                        dd2 = 2 * dd * d, du2 = 2 * du * u,
                        dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          
          Numeric dConst1;
          if(Const1 > 0.)
          {
            dConst1 = db2 * ( b2 * 0.5 +  c2 +  d2 -  u2 -  v2 +  w2);
            dConst1 += b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2);
            
            dConst1 += dc2 * (c2 * 0.5 +  d2 -  u2 +  v2 -  w2);
            dConst1 += c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2);
            
            dConst1 += dd2 * (d2 * 0.5 +  u2 -  v2 -  w2);
            dConst1 += d2 * (dd2 * 0.5 + du2 - dv2 - dw2);
            
            dConst1 += du2 * (u2 * 0.5 +  v2 +  w2);
            dConst1 += u2 * (du2 * 0.5 + dv2 + dw2);
            
            dConst1 += dv2 * (v2 * 0.5 +  w2);
            dConst1 += v2 * (dv2 * 0.5 + dw2);
            
            dConst1 += 4 * ((db *  d *  u *  w - db *  c *  v *  w - dc *  d *  u *  v +
            b * dd *  u *  w -  b * dc *  v *  w -  c * dd *  u *  v +
            b *  d * du *  w -  b *  c * dv *  w -  c *  d * du *  v +
            b *  d *  u * dw -  b *  c *  v * dw -  c *  d *  u * dv));
            dConst1 += dw2 * w2;
            dConst1 /= Const1;
          }
          else
            dConst1 = 0.0;
          
          Numeric dC0, dC1, dC2, dC3;
          if(x == 0.0)
          {
            const Numeric dy = (0.5 * (dConst2 - dConst1)/sqrt_BmA).imag() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*y*dy * C2 * inv_x2y2 + dy*sin_y * inv_x2y2;
            dC3 = -2*y*dy * C3 * inv_x2y2 + (dy*sin_y*inv_y2 - cos_y*dy*inv_y) * inv_x2y2;;
          }
          else if(y == 0.0)
          {
            const Numeric dx = (0.5 * (dConst2 + dConst1)/sqrt_BpA).real() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*x*dx * C2 * inv_x2y2 + dx*sinh_x * inv_x2y2;
            dC3 = -2*x*dx * C3 * inv_x2y2 + (cosh_x*dx*inv_x - dx*sinh_x*inv_x2) * inv_x2y2;
          }
          else
          { 
            const Numeric dx = (0.5 * (dConst2 + dConst1)/sqrt_BpA).real() * sqrt_05;
            const Numeric dy = (0.5 * (dConst2 - dConst1)/sqrt_BmA).imag() * sqrt_05;
            const Numeric dy2 = 2 * y * dy;
            const Numeric dx2 = 2 * x * dx;
            const Numeric dx2dy2 = dx2 + dy2;
            
            dC0 = -dx2dy2 * C0 * inv_x2y2 + 
            (2*cos_y*dx*x + 2*cosh_x*dy*y + dx*sinh_x*y2 - dy*sin_y*x2) * inv_x2y2;
            
            dC1 = -dx2dy2 * C1 * inv_x2y2 + 
            (cos_y*dy*x2*inv_y + dx2*sin_y*inv_y - dy*sin_y*x2*inv_y2 - 
            dx*sinh_x*y2*inv_x2 + cosh_x*dx*y2*inv_x + dy2*sinh_x*inv_x) * inv_x2y2;
            
            dC2 =  -dx2dy2 * C2 * inv_x2y2 + (dx*sinh_x + dy*sin_y) * inv_x2y2;
            
            dC3 =  -dx2dy2 * C3 * inv_x2y2 + 
            (dy*sin_y*inv_y2 - cos_y*dy*inv_y + cosh_x*dx*inv_x - dx*sinh_x*inv_x2) * inv_x2y2;
          }
          
          dF(0, 0) = dF(1, 1) = dF(2, 2) = dF(3, 3) = dC0;
          dF(0, 0) += dC2 * ( b2 +  c2 +  d2) + C2 * (db2 + dc2 + dd2);
          dF(1, 1) += dC2 * ( b2 -  u2 -  v2) + C2 * (db2 - du2 - dv2);
          dF(2, 2) += dC2 * ( c2 -  u2 -  w2) + C2 * (dc2 - du2 - dw2);
          dF(3, 3) += dC2 * ( d2 -  v2 -  w2) + C2 * (dd2 - dv2 - dw2);
          
          
          dF(0, 1) = dF(1, 0) = db * C1 + b * dC1;
          
          dF(0, 1) += dC2 * (- c *  u -  d *  v) 
          +   C2 * (-dc *  u - dd *  v 
          -           c * du -  d * dv)
          +  dC3 * ( b * ( b2 +  c2 +  d2) -  u * ( b *  u -  d *  w) -  v * ( b *  v +  c *  w))
          +   C3 * (db * ( b2 +  c2 +  d2) - du * ( b *  u -  d *  w) - dv * ( b *  v +  c *  w)
          +          b * (db2 + dc2 + dd2) -  u * (db *  u - dd *  w) -  v * (db *  v + dc *  w)
          -                                   u * ( b * du -  d * dw) -  v * ( b * dv +  c * dw));
          dF(1, 0) += dC2 * ( c *  u +  d *  v) 
          +   C2 * (dc *  u + dd *  v
          +          c * du +  d * dv)
          +  dC3 * (- b * (- b2 +  u2 +  v2) +  c * ( b *  c -  v *  w) +  d * ( b *  d +  u *  w))
          +   C3 * (-db * (- b2 +  u2 +  v2) + dc * ( b *  c -  v *  w) + dd * ( b *  d +  u *  w)
          -           b * (-db2 + du2 + dv2) +  c * (db *  c - dv *  w) +  d * (db *  d + du *  w)
          +                                     c * ( b * dc -  v * dw) +  d * ( b * dd +  u * dw));
          
          
          dF(0, 2) = dF(2, 0) = dC1 * c + C1 * dc;
          
          dF(0, 2) += dC2 * ( b *  u -  d *  w)
          +   C2 * (db *  u - dd *  w
          +          b * du -  d * dw)
          +  dC3 * ( c * ( b2 +  c2 +  d2) -  u * ( c *  u +  d *  v) -  w * ( b *  v +  c *  w))
          +   C3 * (dc * ( b2 +  c2 +  d2) - du * ( c *  u +  d *  v) - dw * ( b *  v +  c *  w)
          +          c * (db2 + dc2 + dd2) -  u * (dc *  u + dd *  v) -  w * (db *  v + dc *  w)
          -                                   u * ( c * du +  d * dv) -  w * ( b * dv +  c * dw));
          dF(2, 0) += dC2 * (- b *  u +  d *  w)
          +   C2 * (-db *  u + dd *  w
          -           b * du +  d * dw)
          +  dC3 * ( b * ( b *  c -  v *  w) -  c * (- c2 +  u2 +  w2) +  d * ( c *  d -  u *  v))
          +   C3 * (db * ( b *  c -  v *  w) - dc * (- c2 +  u2 +  w2) + dd * ( c *  d -  u *  v)
          +          b * (db *  c - dv *  w) -  c * (-dc2 + du2 + dw2) +  d * (dc *  d - du *  v)
          +          b * ( b * dc -  v * dw)                           +  d * ( c * dd -  u * dv));
          
          
          dF(0, 3) = dF(3, 0) = dC1 * d + C1 * dd;
          
          dF(0, 3) += dC2 * ( b *  v +  c *  w)
          +   C2 * (db *  v + dc *  w
          +          b * dv +  c * dw)
          +  dC3 * ( d * ( b2 +  c2 +  d2) -  v * ( c *  u +  d *  v) +  w * ( b *  u -  d *  w))
          +   C3 * (dd * ( b2 +  c2 +  d2) - dv * ( c *  u +  d *  v) + dw * ( b *  u -  d *  w)
          +          d * (db2 + dc2 + dd2) -  v * (dc *  u + dd *  v) +  w * (db *  u - dd *  w)
          -                                   v * ( c * du +  d * dv) +  w * ( b * du -  d * dw));
          dF(3, 0) += dC2 * (- b *  v -  c *  w)
          +   C2 * (-db *  v - dc *  w
          -           b * dv -  c * dw)
          +  dC3 * ( b * ( b *  d +  u *  w) +  c * ( c *  d -  u *  v) -  d * (- d2 +  v2 +  w2))
          +   C3 * (db * ( b *  d +  u *  w) + dc * ( c *  d -  u *  v) - dd * (- d2 +  v2 +  w2)
          +          b * (db *  d + du *  w) +  c * (dc *  d - du *  v) -  d * (-dd2 + dv2 + dw2)
          +          b * ( b * dd +  u * dw) +  c * ( c * dd -  u * dv)                          );
          
          
          dF(1, 2) = dF(2, 1) = dC2 * (b * c - v * w) + C2 * (db * c + b * dc - dv * w - v * dw);
          
          dF(1, 2) += dC1 *  u 
          +   C1 * du
          +  dC3 * ( c * ( c *  u +  d *  v) -  u * (- b2 +  u2 +  v2) -  w * ( b *  d +  u *  w))
          +   C3 * (dc * ( c *  u +  d *  v) - du * (- b2 +  u2 +  v2) - dw * ( b *  d +  u *  w)
          +          c * (dc *  u + dd *  v) -  u * (-db2 + du2 + dv2) -  w * (db *  d + du *  w)
          +          c * ( c * du +  d * dv)                           -  w * ( b * dd +  u * dw));
          dF(2, 1) += -dC1 *  u
          -    C1 * du
          +   dC3 * (- b * ( b *  u -  d *  w) +  u * (- c2 +  u2 +  w2) -  v * ( c *  d -  u *  v))
          +    C3 * (-db * ( b *  u -  d *  w) + du * (- c2 +  u2 +  w2) - dv * ( c *  d -  u *  v)
          -            b * (db *  u - dd *  w) +  u * (-dc2 + du2 + dw2) -  v * (dc *  d - du *  v)
          -            b * ( b * du -  d * dw)                           -  v * ( c * dd -  u * dv));
          
          
          dF(1, 3) = dF(3, 1) = dC2 * (b * d + u * w) + C2 * (db * d + b * dd + du * w + u * dw);
          
          dF(1, 3) += dC1 *  v
          +   C1 * dv
          +  dC3 * ( d * ( c *  u +  d *  v) -  v * (- b2 +  u2 +  v2) +  w * ( b *  c -  v *  w))
          +   C3 * (dd * ( c *  u +  d *  v) - dv * (- b2 +  u2 +  v2) + dw * ( b *  c -  v *  w)
          +          d * (dc *  u + dd *  v) -  v * (-db2 + du2 + dv2) +  w * (db *  c - dv *  w)
          +          d * ( c * du +  d * dv)                           +  w * ( b * dc -  v * dw));
          dF(3, 1) += -dC1 *  v
          -    C1 * dv
          +   dC3 * (- b * ( b *  v +  c *  w) -  u * ( c *  d -  u *  v) +  v * (- d2 +  v2 +  w2))
          +    C3 * (-db * ( b *  v +  c *  w) - du * ( c *  d -  u *  v) + dv * (- d2 +  v2 +  w2)
          -            b * (db *  v + dc *  w) -  u * (dc *  d - du *  v) +  v * (-dd2 + dv2 + dw2)
          -            b * ( b * dv +  c * dw) -  u * ( c * dd -  u * dv)                          );
          
          
          dF(2, 3) = dF(3, 2) = dC2 * (c * d - u * v) + C2 * (dc * d + c * dd - du * v - u * dv);
          
          dF(2, 3) += dC1 *  w
          +   C1 * dw
          +  dC3 * (- d * ( b *  u -  d *  w) +  v * ( b *  c -  v *  w) -  w * (- c2 +  u2 +  w2))
          +   C3 * (-dd * ( b *  u -  d *  w) + dv * ( b *  c -  v *  w) - dw * (- c2 +  u2 +  w2)
          -           d * (db *  u - dd *  w) +  v * (db *  c - dv *  w) -  w * (-dc2 + du2 + dw2)
          -           d * ( b * du -  d * dw) +  v * ( b * dc -  v * dw)                          );
          dF(3, 2) += -dC1 *  w
          -    C1 * dw
          +   dC3 * (- c * ( b *  v +  c *  w) +  u * ( b *  d +  u *  w) +  w * (- d2 +  v2 +  w2))
          +    C3 * (-dc * ( b *  v +  c *  w) + du * ( b *  d +  u *  w) + dw * (- d2 +  v2 +  w2)
          -            c * (db *  v + dc *  w) +  u * (db *  d + du *  w) +  w * (-dd2 + dv2 + dw2)
          -            c * ( b * dv +  c * dw) +  u * ( b * dd +  u * dw)                          );
          
          dF *= exp_a;
          
          // Finalize derivation by the chain rule
          dF(0, 0) += F(0, 0) * da;
          dF(0, 1) += F(0, 1) * da;
          dF(0, 2) += F(0, 2) * da;
          dF(0, 3) += F(0, 3) * da;
          dF(1, 0) += F(1, 0) * da;
          dF(1, 1) += F(1, 1) * da;
          dF(1, 2) += F(1, 2) * da;
          dF(1, 3) += F(1, 3) * da;
          dF(2, 0) += F(2, 0) * da;
          dF(2, 1) += F(2, 1) * da;
          dF(2, 2) += F(2, 2) * da;
          dF(2, 3) += F(2, 3) * da;
          dF(3, 0) += F(3, 0) * da;
          dF(3, 1) += F(3, 1) * da;
          dF(3, 2) += F(3, 2) * da;
          dF(3, 3) += F(3, 3) * da; 
        }
        
        for(Index j = 0; j < nppd; j++)
        {
          if(not dlower_level_dx[j].NumberOfFrequencies())
            continue;
          
          MatrixView dF = dT_dx_lower_level(j, i, joker, joker);
          
          const Numeric da = -0.5 * r * dlower_level_dx[j].Kjj(iz, ia)[i],
                        db = -0.5 * r * dlower_level_dx[j].K12(iz, ia)[i],
                        dc = -0.5 * r * dlower_level_dx[j].K13(iz, ia)[i],
                        dd = -0.5 * r * dlower_level_dx[j].K14(iz, ia)[i],
                        du = -0.5 * r * dlower_level_dx[j].K23(iz, ia)[i],
                        dv = -0.5 * r * dlower_level_dx[j].K24(iz, ia)[i],
                        dw = -0.5 * r * dlower_level_dx[j].K34(iz, ia)[i];
          
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c,
                        dd2 = 2 * dd * d, du2 = 2 * du * u,
                        dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          
          Numeric dConst1;
          if(Const1 > 0.)
          {
            dConst1 = db2 * ( b2 * 0.5 +  c2 +  d2 -  u2 -  v2 +  w2);
            dConst1 += b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2);
            
            dConst1 += dc2 * (c2 * 0.5 +  d2 -  u2 +  v2 -  w2);
            dConst1 += c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2);
            
            dConst1 += dd2 * (d2 * 0.5 +  u2 -  v2 -  w2);
            dConst1 += d2 * (dd2 * 0.5 + du2 - dv2 - dw2);
            
            dConst1 += du2 * (u2 * 0.5 +  v2 +  w2);
            dConst1 += u2 * (du2 * 0.5 + dv2 + dw2);
            
            dConst1 += dv2 * (v2 * 0.5 +  w2);
            dConst1 += v2 * (dv2 * 0.5 + dw2);
            
            dConst1 += 4 * ((db *  d *  u *  w - db *  c *  v *  w - dc *  d *  u *  v +
            b * dd *  u *  w -  b * dc *  v *  w -  c * dd *  u *  v +
            b *  d * du *  w -  b *  c * dv *  w -  c *  d * du *  v +
            b *  d *  u * dw -  b *  c *  v * dw -  c *  d *  u * dv));
            dConst1 += dw2 * w2;
            dConst1 /= Const1;
          }
          else
            dConst1 = 0.0;
          
          Numeric dC0, dC1, dC2, dC3;
          if(x == 0.0)
          {
            const Numeric dy = (0.5 * (dConst2 - dConst1)/sqrt_BmA).imag() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*y*dy * C2 * inv_x2y2 + dy*sin_y * inv_x2y2;
            dC3 = -2*y*dy * C3 * inv_x2y2 + (dy*sin_y*inv_y2 - cos_y*dy*inv_y) * inv_x2y2;;
          }
          else if(y == 0.0)
          {
            const Numeric dx = (0.5 * (dConst2 + dConst1)/sqrt_BpA).real() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*x*dx * C2 * inv_x2y2 + dx*sinh_x * inv_x2y2;
            dC3 = -2*x*dx * C3 * inv_x2y2 + (cosh_x*dx*inv_x - dx*sinh_x*inv_x2) * inv_x2y2;
          }
          else
          { 
            const Numeric dx = (0.5 * (dConst2 + dConst1)/sqrt_BpA).real() * sqrt_05;
            const Numeric dy = (0.5 * (dConst2 - dConst1)/sqrt_BmA).imag() * sqrt_05;
            const Numeric dy2 = 2 * y * dy;
            const Numeric dx2 = 2 * x * dx;
            const Numeric dx2dy2 = dx2 + dy2;
            
            dC0 = -dx2dy2 * C0 * inv_x2y2 + 
            (2*cos_y*dx*x + 2*cosh_x*dy*y + dx*sinh_x*y2 - dy*sin_y*x2) * inv_x2y2;
            
            dC1 = -dx2dy2 * C1 * inv_x2y2 + 
            (cos_y*dy*x2*inv_y + dx2*sin_y*inv_y - dy*sin_y*x2*inv_y2 - 
            dx*sinh_x*y2*inv_x2 + cosh_x*dx*y2*inv_x + dy2*sinh_x*inv_x) * inv_x2y2;
            
            dC2 =  -dx2dy2 * C2 * inv_x2y2 + (dx*sinh_x + dy*sin_y) * inv_x2y2;
            
            dC3 =  -dx2dy2 * C3 * inv_x2y2 + 
            (dy*sin_y*inv_y2 - cos_y*dy*inv_y + cosh_x*dx*inv_x - dx*sinh_x*inv_x2) * inv_x2y2;
          }
          
          dF(0, 0) = dF(1, 1) = dF(2, 2) = dF(3, 3) = dC0;
          dF(0, 0) += dC2 * ( b2 +  c2 +  d2) + C2 * (db2 + dc2 + dd2);
          dF(1, 1) += dC2 * ( b2 -  u2 -  v2) + C2 * (db2 - du2 - dv2);
          dF(2, 2) += dC2 * ( c2 -  u2 -  w2) + C2 * (dc2 - du2 - dw2);
          dF(3, 3) += dC2 * ( d2 -  v2 -  w2) + C2 * (dd2 - dv2 - dw2);
          
          
          dF(0, 1) = dF(1, 0) = db * C1 + b * dC1;
          
          dF(0, 1) += dC2 * (- c *  u -  d *  v) 
          +   C2 * (-dc *  u - dd *  v 
          -           c * du -  d * dv)
          +  dC3 * ( b * ( b2 +  c2 +  d2) -  u * ( b *  u -  d *  w) -  v * ( b *  v +  c *  w))
          +   C3 * (db * ( b2 +  c2 +  d2) - du * ( b *  u -  d *  w) - dv * ( b *  v +  c *  w)
          +          b * (db2 + dc2 + dd2) -  u * (db *  u - dd *  w) -  v * (db *  v + dc *  w)
          -                                   u * ( b * du -  d * dw) -  v * ( b * dv +  c * dw));
          dF(1, 0) += dC2 * ( c *  u +  d *  v) 
          +   C2 * (dc *  u + dd *  v
          +          c * du +  d * dv)
          +  dC3 * (- b * (- b2 +  u2 +  v2) +  c * ( b *  c -  v *  w) +  d * ( b *  d +  u *  w))
          +   C3 * (-db * (- b2 +  u2 +  v2) + dc * ( b *  c -  v *  w) + dd * ( b *  d +  u *  w)
          -           b * (-db2 + du2 + dv2) +  c * (db *  c - dv *  w) +  d * (db *  d + du *  w)
          +                                     c * ( b * dc -  v * dw) +  d * ( b * dd +  u * dw));
          
          
          dF(0, 2) = dF(2, 0) = dC1 * c + C1 * dc;
          
          dF(0, 2) += dC2 * ( b *  u -  d *  w)
          +   C2 * (db *  u - dd *  w
          +          b * du -  d * dw)
          +  dC3 * ( c * ( b2 +  c2 +  d2) -  u * ( c *  u +  d *  v) -  w * ( b *  v +  c *  w))
          +   C3 * (dc * ( b2 +  c2 +  d2) - du * ( c *  u +  d *  v) - dw * ( b *  v +  c *  w)
          +          c * (db2 + dc2 + dd2) -  u * (dc *  u + dd *  v) -  w * (db *  v + dc *  w)
          -                                   u * ( c * du +  d * dv) -  w * ( b * dv +  c * dw));
          dF(2, 0) += dC2 * (- b *  u +  d *  w)
          +   C2 * (-db *  u + dd *  w
          -           b * du +  d * dw)
          +  dC3 * ( b * ( b *  c -  v *  w) -  c * (- c2 +  u2 +  w2) +  d * ( c *  d -  u *  v))
          +   C3 * (db * ( b *  c -  v *  w) - dc * (- c2 +  u2 +  w2) + dd * ( c *  d -  u *  v)
          +          b * (db *  c - dv *  w) -  c * (-dc2 + du2 + dw2) +  d * (dc *  d - du *  v)
          +          b * ( b * dc -  v * dw)                           +  d * ( c * dd -  u * dv));
          
          
          dF(0, 3) = dF(3, 0) = dC1 * d + C1 * dd;
          
          dF(0, 3) += dC2 * ( b *  v +  c *  w)
          +   C2 * (db *  v + dc *  w
          +          b * dv +  c * dw)
          +  dC3 * ( d * ( b2 +  c2 +  d2) -  v * ( c *  u +  d *  v) +  w * ( b *  u -  d *  w))
          +   C3 * (dd * ( b2 +  c2 +  d2) - dv * ( c *  u +  d *  v) + dw * ( b *  u -  d *  w)
          +          d * (db2 + dc2 + dd2) -  v * (dc *  u + dd *  v) +  w * (db *  u - dd *  w)
          -                                   v * ( c * du +  d * dv) +  w * ( b * du -  d * dw));
          dF(3, 0) += dC2 * (- b *  v -  c *  w)
          +   C2 * (-db *  v - dc *  w
          -           b * dv -  c * dw)
          +  dC3 * ( b * ( b *  d +  u *  w) +  c * ( c *  d -  u *  v) -  d * (- d2 +  v2 +  w2))
          +   C3 * (db * ( b *  d +  u *  w) + dc * ( c *  d -  u *  v) - dd * (- d2 +  v2 +  w2)
          +          b * (db *  d + du *  w) +  c * (dc *  d - du *  v) -  d * (-dd2 + dv2 + dw2)
          +          b * ( b * dd +  u * dw) +  c * ( c * dd -  u * dv)                          );
          
          
          dF(1, 2) = dF(2, 1) = dC2 * (b * c - v * w) + C2 * (db * c + b * dc - dv * w - v * dw);
          
          dF(1, 2) += dC1 *  u 
          +   C1 * du
          +  dC3 * ( c * ( c *  u +  d *  v) -  u * (- b2 +  u2 +  v2) -  w * ( b *  d +  u *  w))
          +   C3 * (dc * ( c *  u +  d *  v) - du * (- b2 +  u2 +  v2) - dw * ( b *  d +  u *  w)
          +          c * (dc *  u + dd *  v) -  u * (-db2 + du2 + dv2) -  w * (db *  d + du *  w)
          +          c * ( c * du +  d * dv)                           -  w * ( b * dd +  u * dw));
          dF(2, 1) += -dC1 *  u
          -    C1 * du
          +   dC3 * (- b * ( b *  u -  d *  w) +  u * (- c2 +  u2 +  w2) -  v * ( c *  d -  u *  v))
          +    C3 * (-db * ( b *  u -  d *  w) + du * (- c2 +  u2 +  w2) - dv * ( c *  d -  u *  v)
          -            b * (db *  u - dd *  w) +  u * (-dc2 + du2 + dw2) -  v * (dc *  d - du *  v)
          -            b * ( b * du -  d * dw)                           -  v * ( c * dd -  u * dv));
          
          
          dF(1, 3) = dF(3, 1) = dC2 * (b * d + u * w) + C2 * (db * d + b * dd + du * w + u * dw);
          
          dF(1, 3) += dC1 *  v
          +   C1 * dv
          +  dC3 * ( d * ( c *  u +  d *  v) -  v * (- b2 +  u2 +  v2) +  w * ( b *  c -  v *  w))
          +   C3 * (dd * ( c *  u +  d *  v) - dv * (- b2 +  u2 +  v2) + dw * ( b *  c -  v *  w)
          +          d * (dc *  u + dd *  v) -  v * (-db2 + du2 + dv2) +  w * (db *  c - dv *  w)
          +          d * ( c * du +  d * dv)                           +  w * ( b * dc -  v * dw));
          dF(3, 1) += -dC1 *  v
          -    C1 * dv
          +   dC3 * (- b * ( b *  v +  c *  w) -  u * ( c *  d -  u *  v) +  v * (- d2 +  v2 +  w2))
          +    C3 * (-db * ( b *  v +  c *  w) - du * ( c *  d -  u *  v) + dv * (- d2 +  v2 +  w2)
          -            b * (db *  v + dc *  w) -  u * (dc *  d - du *  v) +  v * (-dd2 + dv2 + dw2)
          -            b * ( b * dv +  c * dw) -  u * ( c * dd -  u * dv)                          );
          
          
          dF(2, 3) = dF(3, 2) = dC2 * (c * d - u * v) + C2 * (dc * d + c * dd - du * v - u * dv);
          
          dF(2, 3) += dC1 *  w
          +   C1 * dw
          +  dC3 * (- d * ( b *  u -  d *  w) +  v * ( b *  c -  v *  w) -  w * (- c2 +  u2 +  w2))
          +   C3 * (-dd * ( b *  u -  d *  w) + dv * ( b *  c -  v *  w) - dw * (- c2 +  u2 +  w2)
          -           d * (db *  u - dd *  w) +  v * (db *  c - dv *  w) -  w * (-dc2 + du2 + dw2)
          -           d * ( b * du -  d * dw) +  v * ( b * dc -  v * dw)                          );
          dF(3, 2) += -dC1 *  w
          -    C1 * dw
          +   dC3 * (- c * ( b *  v +  c *  w) +  u * ( b *  d +  u *  w) +  w * (- d2 +  v2 +  w2))
          +    C3 * (-dc * ( b *  v +  c *  w) + du * ( b *  d +  u *  w) + dw * (- d2 +  v2 +  w2)
          -            c * (db *  v + dc *  w) +  u * (db *  d + du *  w) +  w * (-dd2 + dv2 + dw2)
          -            c * ( b * dv +  c * dw) +  u * ( b * dd +  u * dw)                          );
          
          dF *= exp_a;
          
          // Finalize derivation by the chain rule
          dF(0, 0) += F(0, 0) * da;
          dF(0, 1) += F(0, 1) * da;
          dF(0, 2) += F(0, 2) * da;
          dF(0, 3) += F(0, 3) * da;
          dF(1, 0) += F(1, 0) * da;
          dF(1, 1) += F(1, 1) * da;
          dF(1, 2) += F(1, 2) * da;
          dF(1, 3) += F(1, 3) * da;
          dF(2, 0) += F(2, 0) * da;
          dF(2, 1) += F(2, 1) * da;
          dF(2, 2) += F(2, 2) * da;
          dF(2, 3) += F(2, 3) * da;
          dF(3, 0) += F(3, 0) * da;
          dF(3, 1) += F(3, 1) * da;
          dF(3, 2) += F(3, 2) * da;
          dF(3, 3) += F(3, 3) * da; 
        }
      }
    }
  }
}


void PropagationMatrix::AddZeemanPiComponent(ConstVectorView attenuation, 
                                             ConstVectorView phase, 
                                             const Numeric& nd, 
                                             const Numeric& theta, 
                                             const Numeric& eta,
                                             ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, C2E=cos(2*eta), S2E=sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ST2    * attenuation[iv] * extra[iv];
      K12()[iv] -= nd * ST2C2E * attenuation[iv] * extra[iv];
      K13()[iv] -= nd * ST2S2E * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * ST2S2E * phase[iv] * extra[iv];
        K34()[iv] += 2 * nd * ST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ST2    * attenuation[iv];
      K12()[iv] -= nd * ST2C2E * attenuation[iv];
      K13()[iv] -= nd * ST2S2E * attenuation[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * ST2S2E * phase[iv];
        K34()[iv] += 2 * nd * ST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanPiComponentThetaDerivative(ConstVectorView attenuation, 
                                                            ConstVectorView phase, 
                                                            const Numeric& nd, 
                                                            const Numeric& theta, 
                                                            const Numeric& eta,
                                                            ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), dST=cos(theta), dST2=2*ST*dST, dST2C2E = dST2*cos(2*eta), dST2S2E = dST2*sin(2*eta);
  
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * dST2    * attenuation[iv] * extra[iv];
      K12()[iv] -= nd * dST2C2E * attenuation[iv] * extra[iv];
      K13()[iv] -= nd * dST2S2E * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * dST2S2E * phase[iv] * extra[iv];
        K34()[iv] += 2 * nd * dST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * dST2    * attenuation[iv];
      K12()[iv] -= nd * dST2C2E * attenuation[iv];
      K13()[iv] -= nd * dST2S2E * attenuation[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * dST2S2E * phase[iv];
        K34()[iv] += 2 * nd * dST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanPiComponentEtaDerivative(ConstVectorView attenuation, 
                                                          ConstVectorView phase, 
                                                          const Numeric& nd, 
                                                          const Numeric& theta, 
                                                          const Numeric& eta,
                                                          ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, C2E=cos(2*eta), S2E=sin(2*eta), dST2C2E=-2*ST2*S2E, dST2S2E=2*ST2*C2E;
  
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      K12()[iv] -= nd * dST2C2E * attenuation[iv] * extra[iv];
      K13()[iv] -= nd * dST2S2E * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * dST2S2E * phase[iv] * extra[iv];
        K34()[iv] += 2 * nd * dST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      K12()[iv] -= nd * dST2C2E * attenuation[iv];
      K13()[iv] -= nd * dST2S2E * attenuation[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * dST2S2E * phase[iv];
        K34()[iv] += 2 * nd * dST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanPiComponentDerivative(ConstVectorView attenuation,
                                                       ConstVectorView dattenuation, 
                                                       ConstVectorView phase, 
                                                       ConstVectorView dphase, 
                                                       const Numeric& nd, 
                                                       const Numeric& theta, 
                                                       const Numeric& dtheta, 
                                                       const Numeric& eta,
                                                       const Numeric& deta,
                                                       ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, C2E=cos(2*eta), S2E=sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
  const Numeric dST = cos(theta)*dtheta, dST2 = 2*ST*dST, dC2E = -2*S2E*deta, dS2E = 2*C2E*deta,
  dST2C2E = dST2*C2E + ST2*dC2E, dST2S2E = dST2*S2E + ST2*dS2E;
  
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ((ST2    + dST2)    * attenuation[iv] + ST2    * dattenuation[iv]) * extra[iv];
      K12()[iv] -= nd * ((ST2C2E + dST2C2E) * attenuation[iv] + ST2C2E * dattenuation[iv]) * extra[iv];
      K13()[iv] -= nd * ((ST2S2E + dST2S2E) * attenuation[iv] + ST2S2E * dattenuation[iv]) * extra[iv];
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * ((ST2S2E + dST2S2E) * phase[iv] + ST2S2E * dphase[iv]) * extra[iv];
        K34()[iv] += 2 * nd * ((ST2C2E + dST2C2E) * phase[iv] + ST2C2E * dphase[iv]) * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ((ST2    + dST2)    * attenuation[iv] + ST2    * dattenuation[iv]);
      K12()[iv] -= nd * ((ST2C2E + dST2C2E) * attenuation[iv] + ST2C2E * dattenuation[iv]);
      K13()[iv] -= nd * ((ST2S2E + dST2S2E) * attenuation[iv] + ST2S2E * dattenuation[iv]);
      if(not mvectortype)
      {
        K24()[iv] -= 2 * nd * ((ST2S2E + dST2S2E) * phase[iv] + ST2S2E * dphase[iv]);
        K34()[iv] += 2 * nd * ((ST2C2E + dST2C2E) * phase[iv] + ST2C2E * dphase[iv]);
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaMinusComponent(ConstVectorView attenuation, 
                                                     ConstVectorView phase, 
                                                     const Numeric& nd, 
                                                     const Numeric& theta, 
                                                     const Numeric& eta,
                                                     ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, CT=cos(theta), ONEplusCT2=1+CT*CT, C2E=cos(2*eta), S2E=sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ONEplusCT2 * attenuation[iv] * extra[iv];
      K12()[iv] += nd * ST2C2E     * attenuation[iv] * extra[iv];
      K13()[iv] += nd * ST2S2E     * attenuation[iv] * extra[iv];
      K14()[iv] += 2 * nd * CT     * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K23()[iv] += 4 * nd * CT     * phase[iv] * extra[iv];
        K24()[iv] += 2 * nd * ST2S2E * phase[iv] * extra[iv];
        K34()[iv] -= 2 * nd * ST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ONEplusCT2 * attenuation[iv];
      K12()[iv] += nd * ST2C2E     * attenuation[iv];
      K13()[iv] += nd * ST2S2E     * attenuation[iv];
      K14()[iv] += 2 * nd * CT     * attenuation[iv];
      if(not mvectortype)
      {
        K23()[iv] += 4 * nd * CT     * phase[iv];
        K24()[iv] += 2 * nd * ST2S2E * phase[iv];
        K34()[iv] -= 2 * nd * ST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaMinusComponentEtaDerivative(ConstVectorView attenuation, 
                                                                  ConstVectorView phase, 
                                                                  const Numeric& nd,
                                                                  const Numeric& theta,
                                                                  const Numeric& eta,
                                                                  ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, C2E=cos(2*eta), S2E=sin(2*eta), dST2C2E=-2*ST2*S2E, dST2S2E=2*ST2*C2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      K12()[iv] += nd * dST2C2E     * attenuation[iv] * extra[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K24()[iv] += 2 * nd * dST2S2E * phase[iv] * extra[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      K12()[iv] += nd * dST2C2E     * attenuation[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv];
      if(not mvectortype)
      {
        K24()[iv] += 2 * nd * dST2S2E * phase[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaMinusComponentThetaDerivative(ConstVectorView attenuation,
                                                                    ConstVectorView phase,
                                                                    const Numeric& nd,
                                                                    const Numeric& theta,
                                                                    const Numeric& eta,
                                                                    ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), CT=cos(theta), C2E=cos(2*eta), S2E=sin(2*eta);
  const Numeric dST = CT, dST2 = 2*ST*dST, dCT = -ST, dST2C2E = dST2*C2E, dST2S2E = dST2*S2E, dONEplusCT2 = 2*CT*dCT;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * dONEplusCT2 * attenuation[iv] * extra[iv];
      K12()[iv] += nd * dST2C2E     * attenuation[iv] * extra[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv] * extra[iv];
      K14()[iv] += 2 * nd * dCT     * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K23()[iv] += 4 * nd * dCT     * phase[iv] * extra[iv];
        K24()[iv] += 2 * nd * dST2S2E * phase[iv] * extra[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * dONEplusCT2 * attenuation[iv];
      K12()[iv] += nd * dST2C2E     * attenuation[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv];
      K14()[iv] += 2 * nd * dCT     * attenuation[iv];
      if(not mvectortype)
      {
        K23()[iv] += 4 * nd * dCT     * phase[iv];
        K24()[iv] += 2 * nd * dST2S2E * phase[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaMinusComponentDerivative(ConstVectorView attenuation,
                                                               ConstVectorView dattenuation, 
                                                               ConstVectorView phase, 
                                                               ConstVectorView dphase, 
                                                               const Numeric& nd, 
                                                               const Numeric& theta, 
                                                               const Numeric& dtheta, 
                                                               const Numeric& eta,
                                                               const Numeric& deta,
                                                               ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, CT=cos(theta), ONEplusCT2=1+CT*CT, C2E=cos(2*eta), S2E=sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
  const Numeric dST = CT*dtheta, dST2 = 2*ST*dST, dCT = -ST*dtheta, dONEplusCT2 = 2*CT*dCT, dC2E = -2*S2E*deta, dS2E = 2*C2E*deta,
  dST2C2E = dST2*C2E + ST2*dC2E, dST2S2E = dST2*S2E + ST2*dS2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ((ONEplusCT2 + dONEplusCT2) * attenuation[iv] + ONEplusCT2 * dattenuation[iv]) * extra[iv];
      K12()[iv] += nd * ((ST2C2E     + dST2C2E)     * attenuation[iv] + ST2C2E     * dattenuation[iv]) * extra[iv];
      K13()[iv] += nd * ((ST2S2E     + dST2S2E)     * attenuation[iv] + ST2S2E     * dattenuation[iv]) * extra[iv];
      K14()[iv] += 2 * nd * ((CT     + dCT)         * attenuation[iv] + CT         * dattenuation[iv]) * extra[iv];
      if(not mvectortype)
      {
        K23()[iv] += 4 * nd * ((CT     + dCT)     * phase[iv] + CT     * dphase[iv]) * extra[iv];
        K24()[iv] += 2 * nd * ((ST2S2E + dST2S2E) * phase[iv] + ST2S2E * dphase[iv]) * extra[iv];
        K34()[iv] -= 2 * nd * ((ST2C2E + dST2C2E) * phase[iv] + ST2C2E * dphase[iv]) * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ((ONEplusCT2 + dONEplusCT2) * attenuation[iv] + ONEplusCT2 * dattenuation[iv]);
      K12()[iv] += nd * ((ST2C2E     + dST2C2E)     * attenuation[iv] + ST2C2E     * dattenuation[iv]);
      K13()[iv] += nd * ((ST2S2E     + dST2S2E)     * attenuation[iv] + ST2S2E     * dattenuation[iv]);
      K14()[iv] += 2 * nd * ((CT     + dCT)         * attenuation[iv] + CT         * dattenuation[iv]);
      if(not mvectortype)
      {
        K23()[iv] += 4 * nd * ((CT     + dCT)     * phase[iv] + CT     * dphase[iv]);
        K24()[iv] += 2 * nd * ((ST2S2E + dST2S2E) * phase[iv] + ST2S2E * dphase[iv]);
        K34()[iv] -= 2 * nd * ((ST2C2E + dST2C2E) * phase[iv] + ST2C2E * dphase[iv]);
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaPlusComponent(ConstVectorView attenuation,
                                                    ConstVectorView phase,
                                                    const Numeric& nd,
                                                    const Numeric& theta, 
                                                    const Numeric& eta,
                                                    ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, CT=cos(theta), ONEplusCT2=1+CT*CT, C2E=cos(2*eta), S2E=sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ONEplusCT2 * attenuation[iv] * extra[iv];
      K12()[iv] += nd * ST2C2E     * attenuation[iv] * extra[iv];
      K13()[iv] += nd * ST2S2E     * attenuation[iv] * extra[iv];
      K14()[iv] -= 2 * nd * CT     * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K23()[iv] -= 4 * nd * CT     * phase[iv] * extra[iv];
        K24()[iv] += 2 * nd * ST2S2E * phase[iv] * extra[iv];
        K34()[iv] -= 2 * nd * ST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ONEplusCT2 * attenuation[iv];
      K12()[iv] += nd * ST2C2E     * attenuation[iv];
      K13()[iv] += nd * ST2S2E     * attenuation[iv];
      K14()[iv] -= 2 * nd * CT     * attenuation[iv];
      if(not mvectortype)
      {
        K23()[iv] -= 4 * nd * CT     * phase[iv];
        K24()[iv] += 2 * nd * ST2S2E * phase[iv];
        K34()[iv] -= 2 * nd * ST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaPlusComponentThetaDerivative(ConstVectorView attenuation,
                                                                   ConstVectorView phase,
                                                                   const Numeric& nd,
                                                                   const Numeric& theta,
                                                                   const Numeric& eta,
                                                                   ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), CT=cos(theta), C2E=cos(2*eta), S2E=sin(2*eta);
  const Numeric dST=CT, dST2 = 2*ST*dST, dCT=-ST, dONEplusCT2=2*CT*dCT, dST2C2E=dST2*C2E, dST2S2E=dST2*S2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * dONEplusCT2 * attenuation[iv] * extra[iv];
      K12()[iv] += nd * dST2C2E     * attenuation[iv] * extra[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv] * extra[iv];
      K14()[iv] -= 2 * nd * dCT     * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K23()[iv] -= 4 * nd * dCT     * phase[iv] * extra[iv];
        K24()[iv] += 2 * nd * dST2S2E * phase[iv] * extra[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * dONEplusCT2 * attenuation[iv];
      K12()[iv] += nd * dST2C2E     * attenuation[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv];
      K14()[iv] -= 2 * nd * dCT     * attenuation[iv];
      if(not mvectortype)
      {
        K23()[iv] -= 4 * nd * dCT     * phase[iv];
        K24()[iv] += 2 * nd * dST2S2E * phase[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaPlusComponentEtaDerivative(ConstVectorView attenuation,
                                                                 ConstVectorView phase,
                                                                 const Numeric& nd,
                                                                 const Numeric& theta,
                                                                 const Numeric& eta,
                                                                 ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, C2E=cos(2*eta), S2E=sin(2*eta), dST2C2E=-2*ST2*S2E, dST2S2E=2*ST2*C2E;
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      K12()[iv] += nd * dST2C2E     * attenuation[iv] * extra[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv] * extra[iv];
      if(not mvectortype)
      {
        K24()[iv] += 2 * nd * dST2S2E * phase[iv] * extra[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv] * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      K12()[iv] += nd * dST2C2E     * attenuation[iv];
      K13()[iv] += nd * dST2S2E     * attenuation[iv];
      if(not mvectortype)
      {
        K24()[iv] += 2 * nd * dST2S2E * phase[iv];
        K34()[iv] -= 2 * nd * dST2C2E * phase[iv];
      }
    }
  }
}


void PropagationMatrix::AddZeemanSigmaPlusComponentDerivative(ConstVectorView attenuation,
                                                              ConstVectorView dattenuation, 
                                                              ConstVectorView phase, 
                                                              ConstVectorView dphase, 
                                                              const Numeric& nd, 
                                                              const Numeric& theta, 
                                                              const Numeric& dtheta, 
                                                              const Numeric& eta,
                                                              const Numeric& deta,
                                                              ConstVectorView extra)
{
  assert(mstokes_dim == 4);
  assert(mza == 1);
  assert(maa == 1);
  assert(mfreqs == attenuation.nelem() and mfreqs == phase.nelem());
  assert(extra.nelem() == 0 or extra.nelem() == mfreqs);
  
  const Numeric ST=sin(theta), ST2=ST*ST, CT=cos(theta), ONEplusCT2=1+CT*CT, C2E=cos(2*eta), S2E=sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
  const Numeric dST = CT*dtheta, dST2 = 2*ST*dST, dCT = -ST*dtheta, dONEplusCT2 = 2*CT*dCT, dC2E = -2*S2E*deta, dS2E = 2*C2E*deta,
  dST2C2E = dST2*C2E + ST2*dC2E, dST2S2E = dST2*S2E + ST2*dS2E;
  
  
  if(extra.nelem())
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ((ONEplusCT2 + dONEplusCT2) * attenuation[iv] + ONEplusCT2 * dattenuation[iv]) * extra[iv];
      K12()[iv] += nd * ((ST2C2E     + dST2C2E)     * attenuation[iv] + ST2C2E     * dattenuation[iv]) * extra[iv];
      K13()[iv] += nd * ((ST2S2E     + dST2S2E)     * attenuation[iv] + ST2S2E     * dattenuation[iv]) * extra[iv];
      K14()[iv] -= 2 * nd * ((CT     + dCT)         * attenuation[iv] + CT         * dattenuation[iv]) * extra[iv];
      if(not mvectortype)
      {
        K23()[iv] -= 4 * nd * ((CT     + dCT)     * phase[iv] + CT     * dphase[iv]) * extra[iv];
        K24()[iv] += 2 * nd * ((ST2S2E + dST2S2E) * phase[iv] + ST2S2E * dphase[iv]) * extra[iv];
        K34()[iv] -= 2 * nd * ((ST2C2E + dST2C2E) * phase[iv] + ST2C2E * dphase[iv]) * extra[iv];
      }
    }
  }
  else
  {
    for(Index iv = 0; iv < mfreqs; iv++)
    {
      Kjj()[iv] += nd * ((ONEplusCT2 + dONEplusCT2) * attenuation[iv] + ONEplusCT2 * dattenuation[iv]);
      K12()[iv] += nd * ((ST2C2E     + dST2C2E)     * attenuation[iv] + ST2C2E     * dattenuation[iv]);
      K13()[iv] += nd * ((ST2S2E     + dST2S2E)     * attenuation[iv] + ST2S2E     * dattenuation[iv]);
      K14()[iv] -= 2 * nd * ((CT     + dCT)         * attenuation[iv] + CT         * dattenuation[iv]);
      if(not mvectortype)
      {
        K23()[iv] -= 4 * nd * ((CT     + dCT)     * phase[iv] + CT     * dphase[iv]);
        K24()[iv] += 2 * nd * ((ST2S2E + dST2S2E) * phase[iv] + ST2S2E * dphase[iv]);
        K34()[iv] -= 2 * nd * ((ST2C2E + dST2C2E) * phase[iv] + ST2C2E * dphase[iv]);
      }
    }
  }
}


void PropagationMatrix::MatrixInverseAtPosition(MatrixView ret, const Index iv, const Index iz, const Index ia) const
{
  assert(ret.ncols() == mstokes_dim and ret.nrows() == mstokes_dim);
  assert(iv < mfreqs);
  assert(iz < mza);
  assert(ia < maa);
  
  switch(mstokes_dim)
  {
    case 1: 
      assert(Kjj(iz, ia)[iv] not_eq 0); 
      ret(0, 0) = 1.0/Kjj(iz, ia)[iv]; 
      break;
    case 2: 
    {
      const Numeric a = Kjj(iz, ia)[iv], a2 = a*a, 
                    b = K12(iz, ia)[iv], b2 = b*b;
      
      const Numeric f = a2 - b2;
      
      // condition to be invertible
      assert(f not_eq 0);
      
      const Numeric div = 1.0/f;
      
      ret(1, 1) = ret(0, 0) =  Kjj(iz, ia)[iv]*div; 
      ret(1, 0) = ret(0, 1) = -K12(iz, ia)[iv]*div;
    }
      break;
    case 3:
    {
      const Numeric a = Kjj(iz, ia)[iv], a2 = a*a, 
                    b = K12(iz, ia)[iv], b2 = b*b,
                    c = K13(iz, ia)[iv], c2 = c*c,
                    u = K23(iz, ia)[iv], u2 = u*u;
      
      const Numeric f = a*(a2 - b2 - c2 + u2);
      
      // condition to be invertible
      assert(f not_eq 0);
      
      const Numeric div = 1.0/f;
      
      ret(0, 0) = (a2 + u2)*div;
      ret(0, 1) = -(a*b + c*u)*div;
      ret(0, 2) = (-a*c + b*u)*div;
      
      ret(1, 0) = (-a*b + c*u)*div;
      ret(1, 1) = (a2 - c2)*div;
      ret(1, 2) = (-a*u + b*c)*div;
      
      ret(2, 0) = -(a*c + b*u)*div;
      ret(2, 1) = (a*u + b*c)*div;
      ret(2, 2) = (a2 - b2)*div;
    }
      break;
    case 4:
    {
      const Numeric a = Kjj(iz, ia)[iv], a2 = a*a, 
                    b = K12(iz, ia)[iv], b2 = b*b, 
                    c = K13(iz, ia)[iv], c2 = c*c,
                    u = K23(iz, ia)[iv], u2 = u*u,
                    d = K14(iz, ia)[iv], d2 = d*d, 
                    v = K24(iz, ia)[iv], v2 = v*v, 
                    w = K34(iz, ia)[iv], w2 = w*w;
      
      const Numeric f = a2*a2 - a2*b2 - a2*c2 - a2*d2 + a2*u2 + a2*v2 + a2*w2 - b2*w2 + 2*b*c*v*w - 2*b*d*u*w - c2*v2 + 2*c*d*u*v - d2*u2;
      
      // condition to be invertible
      assert(f not_eq 0);
      
      const Numeric div = 1.0/f;
      
      ret(0, 0) = a*(a2 + u2 + v2 + w2)*div;
      ret(0, 1) = (-a2*b - a*c*u - a*d*v - b*w2 + c*v*w - d*u*w)*div;
      ret(0, 2) = (-a2*c + a*b*u - a*d*w + b*v*w - c*v2 + d*u*v)*div;
      ret(0, 3) = (-a2*d + a*b*v + a*c*w - b*u*w + c*u*v - d*u2)*div;
      
      ret(1, 0) = (-a2*b + a*c*u + a*d*v - b*w2 + c*v*w - d*u*w)*div;
      ret(1, 1) = a*(a2 - c2 - d2 + w2)*div;
      ret(1, 2) = (-a2*u + a*b*c - a*v*w + b*d*w - c*d*v + d2*u)*div;
      ret(1, 3) = (-a2*v + a*b*d + a*u*w - b*c*w + c2*v - c*d*u)*div;
      
      ret(2, 0) = (-a2*c - a*b*u + a*d*w + b*v*w - c*v2 + d*u*v)*div;
      ret(2, 1) = (a2*u + a*b*c - a*v*w - b*d*w + c*d*v - d2*u)*div;
      ret(2, 2) = a*(a2 - b2 - d2 + v2)*div;
      ret(2, 3) = (-a2*w + a*c*d - a*u*v + b2*w - b*c*v + b*d*u)*div;
      
      
      ret(3, 0) = (-a2*d - a*b*v - a*c*w - b*u*w + c*u*v - d*u2)*div;
      ret(3, 1) = (a2*v + a*b*d + a*u*w + b*c*w - c2*v + c*d*u)*div;
      ret(3, 2) = (a2*w + a*c*d - a*u*v - b2*w + b*c*v - b*d*u)*div;
      ret(3, 3) = a*(a2 - b2 - c2 + u2)*div;
    }
      break;
    default:
      throw std::runtime_error("Strange stokes dimensions");
      break;
  }
}


void PropagationMatrix::AddAverageAtPosition(ConstMatrixView mat1, ConstMatrixView mat2, const Index iv, const Index iz, const Index ia)
{
  switch(mstokes_dim)
  {
    case 4: mdata(ia ,iz, iv, 3) += (mat1(3, 0) + mat2(3, 0)) * 0.5; mdata(ia ,iz, iv, 5) += (mat1(1, 3) + mat2(1, 3)) * 0.5; mdata(ia ,iz, iv, 6) += (mat1(2, 3) + mat2(2, 3)) * 0.5; /* FALLTHROUGH */
    case 3: mdata(ia ,iz, iv, 2) += (mat1(2, 0) + mat2(2, 0)) * 0.5; mdata(ia ,iz, iv, mstokes_dim) += (mat1(1, 2) + mat2(1, 2)) * 0.5; /* FALLTHROUGH */
    case 2: mdata(ia ,iz, iv, 1) += (mat1(1, 0) + mat2(1, 0)) * 0.5; /* FALLTHROUGH */
    case 1: mdata(ia ,iz, iv, 0) += (mat1(0, 0) + mat2(0, 0)) * 0.5; /* FALLTHROUGH */
  }
}


void PropagationMatrix::MultiplyAndAdd(const Numeric x, const PropagationMatrix& y)
{
  assert(mstokes_dim == y.mstokes_dim);
  assert(mfreqs == y.mfreqs);
  
  for(Index i = 0; i < maa; i++)
    for(Index j = 0; j < mza; j++)
      for(Index k = 0; k < mfreqs; k++)
        for(Index l = 0; l < NumberOfNeededVectors(); l++)
          mdata(i, j, k, l) += x * y.mdata(i, j, k, l);
}


bool PropagationMatrix::FittingShape(ConstMatrixView x) const
{
  Index nelem = x.nrows();
  if(mstokes_dim == nelem)
  {
    if(x.ncols() == nelem)
    {
      switch(mstokes_dim)
      {
        case 1: return true; break;
        case 2: if(x(0, 0) == x(1, 1)) return true; break;
        case 3: if(x(0, 0) == x(1, 1) and x(0, 0) == x(2, 2) and x(0, 1) == x(1, 0) and x(2, 0) == x(0, 2) and x(1, 2) == -x(2, 1)) return true; break;
        case 4: if(x(0, 0) == x(1, 1) and x(0, 0) == x(2, 2) and x(0, 0) == x(3, 3) and x(0, 1) == x(1, 0) and x(2, 0) == x(0, 2) and x(3, 0) == x(0, 3) and x(1, 2) == -x(2, 1) and x(1, 3) == -x(3, 1) and x(3, 2) == -x(2, 3)) return true; break;
        default: throw std::runtime_error("Stokes dimension does not agree with accepted values"); break;
      }
    }
  }
  return false;
}


void PropagationMatrix::GetTensor3(Tensor3View tensor3, Index iz, Index ia)
{
  assert(not mvectortype);
  switch(mstokes_dim)
  {
    case 4: tensor3(joker, 3, 3) = mdata(ia, iz, joker, 0); tensor3(joker, 3, 2) = mdata(ia, iz, joker, 6); 
            tensor3(joker, 3, 1) = mdata(ia, iz, joker, 5); tensor3(joker, 0, 3) = mdata(ia, iz, joker, 3); 
            tensor3(joker, 3, 0) = mdata(ia, iz, joker, 3); tensor3(joker, 2, 3) = mdata(ia, iz, joker, 6); 
            tensor3(joker, 1, 3) = mdata(ia, iz, joker, 5); 
            tensor3(joker, 3, 2) *= -1; tensor3(joker, 3, 1) *= -1; /* FALLTHROUGH */
    case 3: tensor3(joker, 2, 2) = mdata(ia, iz, joker, 0); tensor3(joker, 2, 1) = mdata(ia, iz, joker, mstokes_dim); 
            tensor3(joker, 2, 0) = mdata(ia, iz, joker, 2); tensor3(joker, 0, 2) = mdata(ia, iz, joker, 2); 
            tensor3(joker, 1, 2) = mdata(ia, iz, joker, mstokes_dim); tensor3(joker, 2, 1) *= -1; /* FALLTHROUGH */
    case 2: tensor3(joker, 1, 1) = mdata(ia, iz, joker, 0); 
            tensor3(joker, 1, 0) = mdata(ia, iz, joker, 1); tensor3(joker, 0, 1) = mdata(ia, iz, joker, 1); /* FALLTHROUGH */
    case 1: tensor3(joker, 0, 0) = mdata(ia, iz, joker, 0); break;
    default: throw std::runtime_error("Stokes dimension does not agree with accepted values"); break;
  }
}


void PropagationMatrix::CalculationCase(ArrayOfCaseOfPropagationMatrix& cases) const
{
  assert(mfreqs == cases.nelem());
  for(Index i = 0; i < maa; i++)
  for(Index j = 0; j < mza; j++)
  for(Index k = 0; k < mfreqs; k++)
  {
    if(mdata(i,j,k, joker).sum() == mdata(i,j,k, 0))
      cases[k] = CaseOfPropagationMatrix::Diagonal;
    else
      cases[k] = CaseOfPropagationMatrix::FullDimensional;
  }
}


void PropagationMatrix::CalculationCaseMaximize(ArrayOfCaseOfPropagationMatrix& cases) const
{
  assert(mfreqs == cases.nelem());
  for(Index i = 0; i < maa; i++)
  for(Index j = 0; j < mza; j++)
  for(Index k = 0; k < mfreqs; k++)
  {
    if(mdata(i,j,k, joker).sum() not_eq mdata(i,j,k, 0) and cases[k] == CaseOfPropagationMatrix::Diagonal)
      cases[k] = CaseOfPropagationMatrix::FullDimensional;
  }
}


void PropagationMatrix::LeftMultiplyAtPosition(MatrixView out, ConstMatrixView in, const Index iv, const Index iz, const Index ia) const
{
  assert(not mvectortype);
  assert(out.ncols() == mstokes_dim and out.nrows() == mstokes_dim);
  assert(in.ncols() == mstokes_dim and in.nrows() == mstokes_dim);
  switch(mstokes_dim)
  {
    case 1: 
      out(0, 0) = Kjj(iz, ia)[iv] * in(0, 0); 
      break;
    case 2:
    {
      const Numeric a = Kjj(iz, ia)[iv], b = K12(iz, ia)[iv], 
      m11 = in(0, 0), m12 = in(0, 1),
      m21 = in(1, 0), m22 = in(1, 1);
      out(0, 0) = a*m11 + b*m21;
      out(0, 1) = a*m12 + b*m22;
      out(1, 0) = a*m21 + b*m11;
      out(1, 1) = a*m22 + b*m12;
    }
    break;
    case 3:
    {
      const Numeric a = Kjj(iz, ia)[iv], b = K12(iz, ia)[iv],
      c = K13(iz, ia)[iv], u = K23(iz, ia)[iv], 
      m11 = in(0, 0), m12 = in(0, 1), m13 = in(0, 2),
      m21 = in(1, 0), m22 = in(1, 1), m23 = in(1, 2),
      m31 = in(2, 0), m32 = in(2, 1), m33 = in(2, 2);
      out(0, 0) = a*m11 + b*m21 + c*m31;
      out(0, 1) = a*m12 + b*m22 + c*m32;
      out(0, 2) = a*m13 + b*m23 + c*m33;
      out(1, 0) = a*m21 + b*m11 + m31*u; 
      out(1, 1) = a*m22 + b*m12 + m32*u;
      out(1, 2) = a*m23 + b*m13 + m33*u;
      out(2, 0) = a*m31 + c*m11 - m21*u;
      out(2, 1) = a*m32 + c*m12 - m22*u;
      out(2, 2) = a*m33 + c*m13 - m23*u;
    }
    break;
    case 4:
    {
      const Numeric a = Kjj(iz, ia)[iv], b = K12(iz, ia)[iv],
      c = K13(iz, ia)[iv], u = K23(iz, ia)[iv], 
      d = K14(iz, ia)[iv], v = K24(iz, ia)[iv], w = K34(iz, ia)[iv],
      m11 = in(0, 0), m12 = in(0, 1), m13 = in(0, 2), m14 = in(0, 3),
      m21 = in(1, 0), m22 = in(1, 1), m23 = in(1, 2), m24 = in(1, 3),
      m31 = in(2, 0), m32 = in(2, 1), m33 = in(2, 2), m34 = in(2, 3),
      m41 = in(3, 0), m42 = in(3, 1), m43 = in(3, 2), m44 = in(3, 3);
      out(0, 0) = a*m11 + b*m21 + c*m31 + d*m41;
      out(0, 1) = a*m12 + b*m22 + c*m32 + d*m42;
      out(0, 2) = a*m13 + b*m23 + c*m33 + d*m43;
      out(0, 3) = a*m14 + b*m24 + c*m34 + d*m44;
      out(1, 0) = a*m21 + b*m11 + m31*u + m41*v;
      out(1, 1) = a*m22 + b*m12 + m32*u + m42*v;
      out(1, 2) = a*m23 + b*m13 + m33*u + m43*v;
      out(1, 3) = a*m24 + b*m14 + m34*u + m44*v;
      out(2, 0) = a*m31 + c*m11 - m21*u + m41*w;
      out(2, 1) = a*m32 + c*m12 - m22*u + m42*w;
      out(2, 2) = a*m33 + c*m13 - m23*u + m43*w;
      out(2, 3) = a*m34 + c*m14 - m24*u + m44*w;
      out(3, 0) = a*m41 + d*m11 - m21*v - m31*w;
      out(3, 1) = a*m42 + d*m12 - m22*v - m32*w;
      out(3, 2) = a*m43 + d*m13 - m23*v - m33*w;
      out(3, 3) = a*m44 + d*m14 - m24*v - m34*w;
    }
  }
}


void PropagationMatrix::RightMultiplyAtPosition(MatrixView out, ConstMatrixView in, const Index iv, const Index iz, const Index ia) const
{
  assert(not mvectortype);
  assert(out.ncols() == mstokes_dim and out.nrows() == mstokes_dim);
  assert(in.ncols() == mstokes_dim and in.nrows() == mstokes_dim);
  
  switch(mstokes_dim)
  {
    case 1: 
      out(0, 0) = in(0, 0) * Kjj(iz, ia)[iv];
      break;
    case 2:
    {
      const Numeric a = Kjj(iz, ia)[iv], b = K12(iz, ia)[iv], 
      m11 = in(0, 0), m12 = in(0, 1),
      m21 = in(1, 0), m22 = in(1, 1);
      out(0, 0) = a*m11 + b*m12;
      out(0, 1) = a*m12 + b*m11;
      out(1, 0) = a*m21 + b*m22;
      out(1, 1) = a*m22 + b*m21;
    }               
    break;
    case 3:
    {
      const Numeric a = Kjj(iz, ia)[iv], b = K12(iz, ia)[iv],
      c = K13(iz, ia)[iv], u = K23(iz, ia)[iv], 
      m11 = in(0, 0), m12 = in(0, 1), m13 = in(0, 2),
      m21 = in(1, 0), m22 = in(1, 1), m23 = in(1, 2),
      m31 = in(2, 0), m32 = in(2, 1), m33 = in(2, 2);
      out(0, 0) = a*m11 + b*m12 + c*m13;
      out(0, 1) = a*m12 + b*m11 - m13*u;
      out(0, 2) = a*m13 + c*m11 + m12*u;
      out(1, 0) = a*m21 + b*m22 + c*m23;
      out(1, 1) = a*m22 + b*m21 - m23*u;
      out(1, 2) = a*m23 + c*m21 + m22*u;
      out(2, 0) = a*m31 + b*m32 + c*m33;
      out(2, 1) = a*m32 + b*m31 - m33*u;
      out(2, 2) = a*m33 + c*m31 + m32*u;
    }
    break;
    case 4:
    {
      const Numeric a = Kjj(iz, ia)[iv], b = K12(iz, ia)[iv],
      c = K13(iz, ia)[iv], u = K23(iz, ia)[iv], 
      d = K14(iz, ia)[iv], v = K24(iz, ia)[iv], w = K34(iz, ia)[iv],
      m11 = in(0, 0), m12 = in(0, 1), m13 = in(0, 2), m14 = in(0, 3),
      m21 = in(1, 0), m22 = in(1, 1), m23 = in(1, 2), m24 = in(1, 3),
      m31 = in(2, 0), m32 = in(2, 1), m33 = in(2, 2), m34 = in(2, 3),
      m41 = in(3, 0), m42 = in(3, 1), m43 = in(3, 2), m44 = in(3, 3);
      out(0, 0) = a*m11 + b*m12 + c*m13 + d*m14;
      out(0, 1) = a*m12 + b*m11 - m13*u - m14*v;
      out(0, 2) = a*m13 + c*m11 + m12*u - m14*w;
      out(0, 3) = a*m14 + d*m11 + m12*v + m13*w;
      out(1, 0) = a*m21 + b*m22 + c*m23 + d*m24;
      out(1, 1) = a*m22 + b*m21 - m23*u - m24*v;
      out(1, 2) = a*m23 + c*m21 + m22*u - m24*w;
      out(1, 3) = a*m24 + d*m21 + m22*v + m23*w;
      out(2, 0) = a*m31 + b*m32 + c*m33 + d*m34;
      out(2, 1) = a*m32 + b*m31 - m33*u - m34*v;
      out(2, 2) = a*m33 + c*m31 + m32*u - m34*w;
      out(2, 3) = a*m34 + d*m31 + m32*v + m33*w;
      out(3, 0) = a*m41 + b*m42 + c*m43 + d*m44;
      out(3, 1) = a*m42 + b*m41 - m43*u - m44*v;
      out(3, 2) = a*m43 + c*m41 + m42*u - m44*w;
      out(3, 3) = a*m44 + d*m41 + m42*v + m43*w;
    }
  }
}


void PropagationMatrix::RemoveAtPosition(ConstMatrixView x, const Index iv, const Index iz, const Index ia)
{ 
  switch(mstokes_dim) 
  { 
    case 4: mdata(ia, iz, iv, 5) -= x(1, 3); mdata(ia, iz, iv, 6) -= x(2, 3); mdata(ia, iz, iv, 3) -= x(0, 3); /* FALLTHROUGH */
    case 3: mdata(ia, iz, iv, 2) -= x(0, 2); mdata(ia, iz, iv, mstokes_dim) -= x(1, 2); /* FALLTHROUGH */
    case 2: mdata(ia, iz, iv, 1) -= x(0, 1); /* FALLTHROUGH */
    case 1: mdata(ia, iz, iv, 0) -= x(0, 0); /* FALLTHROUGH */
  }
}


void PropagationMatrix::AddAtPosition(ConstMatrixView x, const Index iv, const Index iz, const Index ia)
{ 
  switch(mstokes_dim) 
  { 
    case 4: mdata(ia, iz, iv, 5) += x(1, 3); mdata(ia, iz, iv, 6) += x(2, 3); mdata(ia, iz, iv, 3) += x(0, 3); /* FALLTHROUGH */
    case 3: mdata(ia, iz, iv, 2) += x(0, 2); mdata(ia, iz, iv, mstokes_dim) += x(1, 2); /* FALLTHROUGH */
    case 2: mdata(ia, iz, iv, 1) += x(0, 1); /* FALLTHROUGH */
    case 1: mdata(ia, iz, iv, 0) += x(0, 0); /* FALLTHROUGH */
  }
}


void PropagationMatrix::MultiplyAtPosition(ConstMatrixView x, const Index iv, const Index iz, const Index ia)
{ 
  switch(mstokes_dim) 
  { 
    case 4: mdata(ia, iz, iv, 5) *= x(1, 3); mdata(ia, iz, iv, 6) *= x(2, 3); mdata(ia, iz, iv, 3) *= x(0, 3); /* FALLTHROUGH */
    case 3: mdata(ia, iz, iv, 2) *= x(0, 2); mdata(ia, iz, iv, mstokes_dim) *= x(1, 2); /* FALLTHROUGH */
    case 2: mdata(ia, iz, iv, 1) *= x(0, 1); /* FALLTHROUGH */
    case 1: mdata(ia, iz, iv, 0) *= x(0, 0); /* FALLTHROUGH */
  }
}


void PropagationMatrix::DivideAtPosition(ConstMatrixView x, const Index iv, const Index iz, const Index ia)
{ 
  switch(mstokes_dim) 
  { 
    case 4: mdata(ia, iz, iv, 5) /= x(1, 3); mdata(ia, iz, iv, 6) /= x(2, 3); mdata(ia, iz, iv, 3) /= x(0, 3); /* FALLTHROUGH */
    case 3: mdata(ia, iz, iv, 2) /= x(0, 2); mdata(ia, iz, iv, mstokes_dim) /= x(1, 2); /* FALLTHROUGH */
    case 2: mdata(ia, iz, iv, 1) /= x(0, 1); /* FALLTHROUGH */
    case 1: mdata(ia, iz, iv, 0) /= x(0, 0); /* FALLTHROUGH */
  }
}


void PropagationMatrix::SetAtPosition(ConstMatrixView x, const Index iv, const Index iz, const Index ia)
{ 
  switch(mstokes_dim) 
  { 
    case 4: mdata(ia, iz, iv, 5)  = x(1, 3); mdata(ia, iz, iv, 6)  = x(2, 3); mdata(ia, iz, iv, 3)  = x(0, 3); /* FALLTHROUGH */
    case 3: mdata(ia, iz, iv, 2)  = x(0, 2); mdata(ia, iz, iv, mstokes_dim)  = x(1, 2); /* FALLTHROUGH */
    case 2: mdata(ia, iz, iv, 1)  = x(0, 1); /* FALLTHROUGH */
    case 1: mdata(ia, iz, iv, 0)  = x(0, 0); /* FALLTHROUGH */
  }
}


void PropagationMatrix::MatrixAtPosition(MatrixView ret, const Index iv, const Index iz, const Index ia) const
{
  switch(mstokes_dim)
  {
    case 4: ret(3, 3) = mdata(ia, iz, iv, 0); ret(3, 1) = -mdata(ia, iz, iv, 5); ret(1, 3) = mdata(ia, iz, iv, 5); ret(3, 2) = -mdata(ia, iz, iv, 6); ret(2, 3) = mdata(ia, iz, iv, 6); ret(0, 3) = ret(3, 0) = mdata(ia, iz, iv, 3); /* FALLTHROUGH */
    case 3: ret(2, 2) = mdata(ia, iz, iv, 0); ret(2, 1) = -mdata(ia, iz, iv, 3); ret(1, 2) = mdata(ia, iz, iv, 3); ret(2, 0) = ret(0, 2) = mdata(ia, iz, iv, 2); /* FALLTHROUGH */
    case 2: ret(1, 1) = mdata(ia, iz, iv, 0); ret(1, 0) = ret(0, 1) = mdata(ia, iz, iv, 1); /* FALLTHROUGH */
    case 1: ret(0, 0) = mdata(ia, iz, iv, 0); /* FALLTHROUGH */
  }
}


Numeric PropagationMatrix::operator()(const Index iv, const Index is1, const Index is2, const Index iz, const Index ia) const
{
  assert(is1 < mstokes_dim and is2 < mstokes_dim);
  assert(iv < mfreqs);
  assert(iz < mza);
  assert(ia < maa);
  
  switch(is1)
  {
    case 0:
      switch(is2)
      {
        case 0: return  mdata(ia, iz, iv, 0); break;
        case 1: return  mdata(ia, iz, iv, 1); break;
        case 2: return  mdata(ia, iz, iv, 2); break;
        case 3: return  mdata(ia, iz, iv, 3); break;
        default: throw std::runtime_error("out of bounds");
      }
      break;
    case 1:
      assert(not mvectortype);
      switch(is2)
      {
        case 0: return  mdata(ia, iz, iv, 1); break;
        case 1: return  mdata(ia, iz, iv, 0); break;
        case 2: return  mdata(ia, iz, iv, mstokes_dim); break;
        case 3: return  mdata(ia, iz, iv, 5); break;
        default: throw std::runtime_error("out of bounds");
      }
    case 2:
      assert(not mvectortype);
      switch(is2)
      {
        case 0: return  mdata(ia, iz, iv, 2); break;
        case 1: return -mdata(ia, iz, iv, mstokes_dim); break; 
        case 2: return  mdata(ia, iz, iv, 0); break;
        case 3: return  mdata(ia, iz, iv, 6); break;
        default: throw std::runtime_error("out of bounds");
      }
      break;
    case 3:
      assert(not mvectortype);
      switch(is2)
      {
        case 0: return  mdata(ia, iz, iv, 3); break;
        case 1: return -mdata(ia, iz, iv, 5); break;
        case 2: return -mdata(ia, iz, iv, 6); break;
        case 3: return  mdata(ia, iz, iv, 0); break;
        default: throw std::runtime_error("out of bounds");
      }
      break;
        default: throw std::runtime_error("out of bounds");
  }
}

// Needs to be implemented in this file!!!
std::ostream& operator<<(std::ostream& os, const PropagationMatrix& pm)
{
  os << pm.GetData() << "\n";
  return os;
}


std::ostream& operator<<(std::ostream& os, const ArrayOfPropagationMatrix& apm)
{
  for(auto& pm : apm)
    os << pm;
  return os;
}


std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfPropagationMatrix& aapm)
{
  for(auto& apm : aapm)
    os << apm;
  return os;
}

// Needs to be implemented in this file!!!
std::ostream& operator<<(std::ostream& os, const StokesVector& sv)
{
  os << sv.GetData() << "\n";
  return os;
}


std::ostream& operator<<(std::ostream& os, const ArrayOfStokesVector& asv)
{
  for(auto& sv : asv)
    os << sv;
  return os;
}


std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfStokesVector& aasv)
{
  for(auto& asv : aasv)
    os << asv;
  return os;
}


void get_diydx_replacement(MatrixView diydx_this,
                           MatrixView diydx_next,
                           ConstMatrixView iy,
                           ConstMatrixView sibi,
                           const StokesVector& nlte_this,
                           const StokesVector& nlte_next,
                           const StokesVector& dnltedx_this,
                           const StokesVector& dnltedx_next,
                           const PropagationMatrix& K_this,
                           const PropagationMatrix& K_next,
                           const PropagationMatrix& dKdx_this,
                           const PropagationMatrix& dKdx_next,
                           ConstTensor3View T_this,
                           ConstTensor3View dTdx_this,
                           ConstTensor3View dTdx_next,
                           ConstTensor3View PiT_this,
                           ConstTensor3View PiT_next,
                           const Numeric& temperature_this,
                           const Numeric& temperature_next,
                           const Numeric& dt,
                           ConstVectorView dBdx_this,
                           ConstVectorView dBdx_next,
                           const Numeric& r,
                           const bool& do_Bsource,
                           const bool& do_HSE,
                           const bool& do_nonlte )
{
  const Index stokes_dim = K_this.StokesDimensions();
  const Index frequency_dim = K_next.NumberOfFrequencies();
  ArrayOfCaseOfPropagationMatrix cases(frequency_dim);
  K_this.CalculationCase(cases);
  K_next.CalculationCaseMaximize(cases);
  
  #pragma omp parallel for \
  if (!arts_omp_in_parallel() and frequency_dim >= arts_omp_get_max_threads())
  for(Index i = 0; i < frequency_dim; i++)
  {
    if(cases[i] == CaseOfPropagationMatrix::Diagonal)
    {
      /*
        * Solves ds/dx = T [ dTi/dx * (si-bi-Ki^(-1)ji) + 
        *        (1 - Ti) * ( dbi/dx + K^(-2)dKi/dxji + K^(-1)dji/dx ) ]
        * 
        * In extmat_case 1, dTi/dx = T*dKi/dx.
        * 
        * (si-bi-Ki^(-1)ji) is sibi
        * T                 is PiT_this                   (NB PiT_next is PiT_this*T_this)
        * dTi/dx            is PiT_next*dKi/dx*r/2        (from T being scalar/diagonal exponent function)
        * dbi/dx            is dBdx_this/2, dBdx_next/2   (note how outer loop takes care of this)
        * Ki                is K_this/2, K_next/2         (note how outer loop takes care of this)
        * dKi/dx            is dKdx_this, dKdx_next       (note how outer loop takes care of this)
        * ji                is nlte_this/2, nlte_next/2   (note how outer loop takes care of this)
        * dji/dx            is dnltedx_this, dnltedx_next (note how outer loop takes care of this)
        * 
        * Begin solution
        */
      
      /*
        * T * dTi/dx * (si-bi-Ki^(-1)ji)
        */
      const Numeric x = -0.5 * r * PiT_next(i, 0, 0);
      const Numeric y = x * sibi(i, 0);
      // Stokes 1:
      diydx_this(i, 0) += y * dKdx_this.Kjj()[i];
      diydx_next(i, 0) += y * dKdx_next.Kjj()[i];
      // Higher Stokes
      for( Index is = 1; is < stokes_dim; is++ )
      { 
          const Numeric z = x * sibi(i, is); //NB: sibi not iy since ji can be polarized as well?
          diydx_this(i, is) += z * dKdx_this.Kjj()[i];
          diydx_next(i, is) += z * dKdx_next.Kjj()[i];
      }
      
      /*
        * T * (1 - Ti) * dbi/dx
        */
      if(do_Bsource)
      {
          // The terms associated with B-bar:  For extmat_case==1 only works on first term
          const Numeric v = PiT_this(i, 0, 0) * (1.0 - T_this(i, 0, 0));
          diydx_this(i, 0) += 0.5 * v * dBdx_this[i];
          diydx_next(i, 0) += 0.5 * v * dBdx_next[i];
      }
      
      /*
        * T * (1 - Ti) *(K^(-2)dKi/dxji + K^(-1)dji/dx)
        */
      if(do_nonlte)
      {
          const Numeric v = 0.5 * PiT_this(i, 0, 0) * (1.0 - T_this(i, 0, 0));
          
          Numeric invK = 1.0/K_this.Kjj()[i];
          
          diydx_this(i, 0) += -v * invK * (dKdx_this.Kjj()[i] * invK * nlte_this.Kjj()[i] + invK * dnltedx_this.Kjj()[i]);
          
          invK = 1.0/K_next.Kjj()[i];
          diydx_next(i ,0) += -v * invK * (dKdx_next.Kjj()[i] * invK * nlte_next.Kjj()[i] + invK * dnltedx_next.Kjj()[i]);
      }
      
      /*
        * End of solution
        */
      
      // Zero for higher Stokes
      //
      // The terms associated with Delta-s:
      if( do_HSE )
      {
          if(do_nonlte)
              throw std::runtime_error("HSE and non-LTE not yet compatible\n");
              
          // Stokes 1:
          const Numeric kbar = 0.5 * ( K_this.Kjj()[i] + K_next.Kjj()[i] );
          
          diydx_this(i, 0) += y * kbar / temperature_this;
          diydx_next(i, 0) += y * kbar / temperature_next;
          // Higher Stokes
          for( Index is=1; is<stokes_dim; is++ )
          { 
              const Numeric z = x * iy(i, is);
              diydx_this(i, is) += z * kbar / temperature_this;
              diydx_next(i, is) += z * kbar / temperature_next;
          }
      } //hse
    }
    else
    {
      /*
      * Solves ds/dx = T [ dTi/dx * (si-bi-Ki^(-1)ji) + 
      *        (1 - Ti) * ( dbi/dx + K^(-1)dKi/dxK^(-1)ji + K^(-1)dji/dx ) ]
      * 
      * (si-bi-Ki^(-1)ji) is sibi
      * T                 is PiT_this
      * dTi/dx            is dTdx_this, dTdx_next       (note how outer loop takes care of this)
      * Ti                is T_this
      * dbi/dx            is dBdx_this/2, dBdx_next/2   (note how outer loop takes care of this)
      * Ki                is K_this/2, K_next/2         (note how outer loop takes care of this)
      * dKi/dx            is dKdx_this, dKdx_next       (note how outer loop takes care of this)
      * ji                is nlte_this/2, nlte_next/2   (note how outer loop takes care of this)
      * dji/dx            is dnltedx_this, dnltedx_next (note how outer loop takes care of this)
      * 
      * Begin solution
      */
      
      //Helpers
      Vector x(stokes_dim), y(stokes_dim);
      /*
      * This is the T * dTi/dx * (si-bi-Ki^(-1)ji) part
      */
      // Disturb for this
      mult( x, dTdx_this(i, joker, joker), sibi(i, joker) ); // local derivative
      mult( y, PiT_this(i, joker, joker), x ); // propagate derivative to sensor
      diydx_this(i, joker) += y;
      
      // Disturb for next
      mult( x, dTdx_next(i, joker, joker), sibi(i, joker) ); // local derivative
      mult( y, PiT_this(i, joker, joker), x ); // propagate derivative to sensor
      diydx_next(i, joker) += y; 
      
      
      // The terms associated with B-bar:
      /*
      * This is the T * (1-Ti)*dbi/dx part
      */
      if(do_Bsource)
      {
        // Disturb this
        const Numeric v = (1.0 - T_this(i, 0, 0));
        x[0] = 0.5 *v * dBdx_this[i];
        for(Index is = 1; is < stokes_dim; is++) 
        {
          x[is] = -0.5 * T_this(i, is, 0) * dBdx_this[i]; 
        }
        mult(y, PiT_this(i, joker, joker), x);
        diydx_this(i, joker) += y; 
        
        // Disturb next
        x[0] = 0.5 * v * dBdx_next[i];
        for(Index is = 1; is < stokes_dim; is++) 
        { 
          x[is] = -0.5 * T_this(i, is, 0) * dBdx_next[i]; 
        }
        mult(y, PiT_this(i, joker, joker), x);
        diydx_next(i, joker) += y;
      }
      
      /*
      * This is the  T * (1-Ti)* ( K^(-1)dKi/dxK^(-1)ji + K^(-1)dji/dx ) part
      */
      if(do_nonlte)
      {
        Matrix invK(stokes_dim,stokes_dim), mat1(stokes_dim,stokes_dim), OneMinusTransmission(stokes_dim,stokes_dim);
        Vector vec1(stokes_dim), vec2(stokes_dim);
        
        OneMinusTransmission = T_this(i, joker, joker);
        OneMinusTransmission *= -1;
        for(Index is = 0; is < stokes_dim; is++)
          OneMinusTransmission(is, is) += 1;
        
        // For this
        K_this.MatrixInverseAtPosition(invK, i);
        dKdx_this.LeftMultiplyAtPosition(mat1, invK, i);
        
        mult(vec1, mat1, nlte_this.VectorAtPosition(i));
        vec1 -= dnltedx_this.VectorAtPosition(i);
        vec1*=-0.5;
        
        mult(vec2, OneMinusTransmission, vec1);
        mult(vec1, PiT_this(i, joker, joker), vec2);
        
        diydx_this(i, joker) += vec1;
        
        // For next
        K_next.MatrixInverseAtPosition(invK, i);
        dKdx_next.LeftMultiplyAtPosition(mat1, invK, i);
        
        mult(vec1, mat1, nlte_next.VectorAtPosition(i));
        vec1 -= dnltedx_next.VectorAtPosition(i);
        vec1 *= -0.5;
        
        mult(vec2, OneMinusTransmission, vec1);
        mult(vec1, PiT_this(i, joker, joker), vec2);
        
        diydx_next(i, joker) += vec1;
      }
      
      
      /*
      *  End of solution
      * 
      * Below follows hydrostatic equilibrium calculations
      * FIXME:  Add theoretical description
      */
      
      //
      // The terms associated with Delta-s:
      if( do_HSE )
      {
        if(do_nonlte)
          throw std::runtime_error("HSE and non-LTE not yet compatible\n");
        
        Matrix ext_mat(stokes_dim,stokes_dim), dtdx(stokes_dim,stokes_dim);
        for( Index is1=0; is1<stokes_dim; is1++ ){
          for( Index is2=0; is2<stokes_dim; is2++ ){
            ext_mat(is1,is2) = 0.5 * ( K_this(i, is1,is2) + K_next(i, is1,is2) );}}
          // dl for disturbed tbar
          const Numeric tbar = 0.5 * ( temperature_this + temperature_next );
        const Numeric dl = r * ( 1.0 + dt/tbar );
        
        ext_mat *= -dl;
        cayley_hamilton_fitted_method_4x4_propmat_to_transmat__explicit(dtdx, ext_mat);
        
        for( Index is1=0; is1<stokes_dim; is1++ ) 
        {
          for( Index is2=0; is2<stokes_dim; is2++ ) 
          {
            dtdx(is1,is2) = (1/dt) * (dtdx(is1,is2) - T_this(i, is1,is2)); 
          }
        }
            
        mult( x, dtdx, sibi(i, joker) );
        mult( y, PiT_this(i, joker, joker), x );
        
        // Contribution shared between the two
        // points  and is proportional to 1/t
        // See also AUG.
        for( Index is=0; is<stokes_dim; is++ ) 
        {
          diydx_this(i, is) += y[is] * 0.5 * tbar / temperature_this;
          diydx_next(i, is) += y[is] * 0.5 * tbar / temperature_next;
        }
      } // HSE
    }
  }
}

