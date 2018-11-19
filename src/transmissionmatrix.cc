/* Copyright (C) 2018
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
 * \file   transmissionmatrix.c
 * \brief  Stuff related to the transmission matrix.
 * 
 * Using Eigen library to try and speed up computations.
 * 
 * \author Richard Larsson
 * \date   2018-01-30
 */

#include "transmissionmatrix.h"
#include "complex.h"


Eigen::Matrix<double, 1, 1> vector1(const StokesVector& a, const ConstVectorView& B, const StokesVector& da, const ConstVectorView& dB_dT, const StokesVector& dS, bool dT, Index i)
{
  return
  (Eigen::Matrix<double, 1, 1>() <<
    dS.Kjj()[i] + da.Kjj()[i] * B[i]).finished() + 
  (dT ? (Eigen::Matrix<double, 1, 1>() <<
    a.Kjj()[i] * dB_dT[i]).finished() : 
    (Eigen::Matrix<double, 1, 1>() << 0).finished());
}


Eigen::Vector2d vector2(const StokesVector& a, const ConstVectorView& B, const StokesVector& da, const ConstVectorView& dB_dT, const StokesVector& dS, bool dT, Index i)
{
  return
  (Eigen::Vector2d() <<
    dS.Kjj()[i] + da.Kjj()[i] * B[i],
    dS.K12()[i] + da.K12()[i] * B[i]).finished() + 
  (dT ? (Eigen::Vector2d() <<
    a.Kjj()[i] * dB_dT[i],
    a.K12()[i] * dB_dT[i]).finished() : 
    (Eigen::Vector2d() << 0, 0).finished());
}

Eigen::Vector3d vector3(const StokesVector& a, const ConstVectorView& B, const StokesVector& da, const ConstVectorView& dB_dT, const StokesVector& dS, bool dT, Index i)
{
  return
  (Eigen::Vector3d() <<
    dS.Kjj()[i] + da.Kjj()[i] * B[i],
    dS.K12()[i] + da.K12()[i] * B[i],
    dS.K13()[i] + da.K13()[i] * B[i]).finished() + 
  (dT ? (Eigen::Vector3d() <<
    a.Kjj()[i] * dB_dT[i],
    a.K12()[i] * dB_dT[i],
    a.K13()[i] * dB_dT[i]).finished() : 
    (Eigen::Vector3d() << 0, 0, 0).finished());
}

Eigen::Vector4d vector4(const StokesVector& a, const ConstVectorView& B, const StokesVector& da, const ConstVectorView& dB_dT, const StokesVector& dS, bool dT, Index i)
{
  return
  (Eigen::Vector4d() <<
    dS.Kjj()[i] + da.Kjj()[i] * B[i],
    dS.K12()[i] + da.K12()[i] * B[i],
    dS.K13()[i] + da.K13()[i] * B[i],
    dS.K14()[i] + da.K14()[i] * B[i]).finished() + 
  (dT ? (Eigen::Vector4d() <<
    a.Kjj()[i] * dB_dT[i],
    a.K12()[i] * dB_dT[i],
    a.K13()[i] * dB_dT[i],
    a.K14()[i] * dB_dT[i]).finished() : 
    (Eigen::Vector4d() << 0, 0, 0, 0).finished());
}


Eigen::MatrixXd vector_of_dJ(const StokesVector& a, const ConstVectorView& B, const StokesVector& da, const ConstVectorView& dB_dT, const StokesVector& dS, bool dT, Index i) {
  switch(a.StokesDimensions()) {
    case 4: return vector4(a, B, da, dB_dT, dS, dT, i);
    case 3: return vector3(a, B, da, dB_dT, dS, dT, i);
    case 2: return vector2(a, B, da, dB_dT, dS, dT, i);
    case 1: return vector1(a, B, da, dB_dT, dS, dT, i);
  }
  return Eigen::Matrix<double, 0, 0>();
}

Eigen::Matrix<double, 1, 1> matrix1(const Numeric& a)
{
  return (Eigen::Matrix<double, 1, 1>() << a).finished();
}

Eigen::Matrix2d matrix2(const Numeric& a, const Numeric& b)
{
  return (Eigen::Matrix2d() << a, b,
                               b, a).finished();
}

Eigen::Matrix3d matrix3(const Numeric& a, const Numeric& b, const Numeric& c, const Numeric& u)
{
  return (Eigen::Matrix3d() << a, b, c,
                               b, a, u,
                               c,-u, a).finished();
}

Eigen::Matrix4d matrix4(const Numeric& a, const Numeric& b, const Numeric& c, const Numeric& d, const Numeric& u, const Numeric& v, const Numeric& w)
{
  return (Eigen::Matrix4d() << a, b, c, d,
                               b, a, u, v,
                               c,-u, a, w,
                               d,-v,-w, a).finished();
}

Eigen::MatrixXd matrix_of_K(const PropagationMatrix& K, const Index i){
  switch(K.StokesDimensions()) {
    case 4: return matrix4(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K14()[i],
      K.K23()[i], K.K24()[i], K.K34()[i]);
    case 3: return matrix3(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K23()[i]);
    case 2: return matrix2(K.Kjj()[i], K.K12()[i]);
    case 1: return matrix1(K.Kjj()[i]);
  }
  return Eigen::Matrix<double, 0, 0>();
}


Eigen::Matrix<double, 1, 1> inv1(const Numeric& x) 
{
  return Eigen::Matrix<double, 1, 1>(1/x);
}

Eigen::Matrix2d inv2(const Numeric& a, const Numeric& b) 
{
  const Numeric f = a*a - b*b;
  return (Eigen::Matrix2d() << a, -b, a, -b).finished() / f;
}

Eigen::Matrix3d inv3(const Numeric& a, const Numeric& b, const Numeric& c, const Numeric& u)
{
  const Numeric a2 = a*a, b2 = b*b, c2 = c*c, u2 = u*u;
  const Numeric f = a*(a2 - b2 - c2 + u2);
  return  (Eigen::Matrix3d() <<  a2 + u2,   -a*b - c*u, -a*c + b*u, 
           -a*b + c*u,  a2 - c2,   -a*u + b*c,
           -a*c - b*u,  a*u + b*c,  a2 - b2).finished() / f;
}


Eigen::Matrix4d inv4(const Numeric& a, const Numeric& b, const Numeric& c, const Numeric& d, const Numeric& u, const Numeric& v, const Numeric& w)
{
  const Numeric a2 = a*a, b2 = b*b, c2 = c*c, u2 = u*u, d2 = d*d, v2 = v*v, w2 = w*w;
  const Numeric f = a2*a2 - a2*b2 - a2*c2 - a2*d2 + a2*u2 + a2*v2 + a2*w2 - b2*w2 + 2*b*c*v*w - 2*b*d*u*w - c2*v2 + 2*c*d*u*v - d2*u2;
  return (Eigen::Matrix4d() << 
          a*(a2 + u2 + v2 + w2), 
          -a2*b - a*c*u - a*d*v - b*w2 + c*v*w - d*u*w,
          -a2*c + a*b*u - a*d*w + b*v*w - c*v2 + d*u*v,
          -a2*d + a*b*v + a*c*w - b*u*w + c*u*v - d*u2,
          -a2*b + a*c*u + a*d*v - b*w2 + c*v*w - d*u*w,
          a*(a2 - c2 - d2 + w2),
          -a2*u + a*b*c - a*v*w + b*d*w - c*d*v + d2*u,
          -a2*v + a*b*d + a*u*w - b*c*w + c2*v - c*d*u,
          -a2*c - a*b*u + a*d*w + b*v*w - c*v2 + d*u*v,
          a2*u + a*b*c - a*v*w - b*d*w + c*d*v - d2*u,
          a*(a2 - b2 - d2 + v2),
          -a2*w + a*c*d - a*u*v + b2*w - b*c*v + b*d*u,
          -a2*d - a*b*v - a*c*w - b*u*w + c*u*v - d*u2,
          a2*v + a*b*d + a*u*w + b*c*w - c2*v + c*d*u,
          a2*w + a*c*d - a*u*v - b2*w + b*c*v - b*d*u,
          a*(a2 - b2 - c2 + u2)).finished() / f;
}

Eigen::MatrixXd inverse_of_K(const PropagationMatrix& K, const Index i){
  switch(K.StokesDimensions()) {
    case 4: return inv4(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K14()[i],
      K.K23()[i], K.K24()[i], K.K34()[i]);
    case 3: return inv3(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K23()[i]);
    case 2: return inv2(K.Kjj()[i], K.K12()[i]);
    case 1: return inv1(K.Kjj()[i]);
  }
  return Eigen::Matrix<double, 0, 0>();
}


void transmat1(TransmissionMatrix& T,
               const PropagationMatrix& K1,
               const PropagationMatrix& K2,
               const Numeric& r,
               const Index iz=0,
               const Index ia=0)
{
  for(Index i=0; i<K1.NumberOfFrequencies(); i++)
    T.set1(std::exp(-0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])), i);
}

void transmat2(TransmissionMatrix& T,
               const PropagationMatrix& K1,
               const PropagationMatrix& K2,
               const Numeric& r,
               const Index iz=0,
               const Index ia=0)
{
  for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]), 
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    Eigen::Matrix2d F;
    
    if(b == 0.)
      F.setIdentity();
    else {
      const Numeric C0 = (b*std::cosh(b) - a*std::sinh(b))/b;
      const Numeric C1 = std::sinh(b)/b;  
      
      F(0, 0) = F(1, 1) = C0 + C1 * a;
      F(0, 1) = F(1, 0) = C1 * b;
    }
    T.set2(F*exp_a, i);
  }
}


void transmat3(TransmissionMatrix& T,
               const PropagationMatrix& K1,
               const PropagationMatrix& K2,
               const Numeric& r,
               const Index iz=0,
               const Index ia=0)
{
  for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]), 
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]), 
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]), 
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    Eigen::Matrix3d F;
    
    if(b == 0. and c == 0. and u == 0.)
      F.setIdentity();
    else {
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      
      const Numeric x = std::sqrt(b2 + c2 - u2), x2 = x * x, inv_x2 = 1.0/x2;
      const Numeric sinh_x = std::sinh(x), cosh_x = std::cosh(x);
      
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
    }
    T.set3(F*exp_a, i);
  }
}


void transmat4(TransmissionMatrix& T,
               const PropagationMatrix& K1,
               const PropagationMatrix& K2,
               const Numeric& r,
               const Index iz=0,
               const Index ia=0)
{ 
    static const Numeric sqrt_05 = sqrt(0.5);
    for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
      const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]), 
                    b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]), 
                    c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]), 
                    d = -0.5 * r * (K1.K14(iz, ia)[i] + K2.K14(iz, ia)[i]), 
                    u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]), 
                    v = -0.5 * r * (K1.K24(iz, ia)[i] + K2.K24(iz, ia)[i]), 
                    w = -0.5 * r * (K1.K34(iz, ia)[i] + K2.K34(iz, ia)[i]);
      const Numeric exp_a = std::exp(a);
      Eigen::Matrix4d F;
      
      if(b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.)
        F.setIdentity();
      else {
        
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
        
        const Complex sqrt_BpA = std::sqrt(Complex(Const2 + Const1, 0.0));
        const Complex sqrt_BmA = std::sqrt(Complex(Const2 - Const1, 0.0));
        const Numeric x = sqrt_BpA.real() * sqrt_05;
        const Numeric y = sqrt_BmA.imag() * sqrt_05;
        const Numeric x2 = x * x;
        const Numeric y2 = y * y;
        const Numeric cos_y = std::cos(y);
        const Numeric sin_y = std::sin(y);
        const Numeric cosh_x = std::cosh(x);
        const Numeric sinh_x = std::sinh(x);
        const Numeric x2y2 = x2 + y2;
        const Numeric inv_x2y2 = 1.0 / x2y2;
        
        Numeric C0, C1, C2, C3;
        Numeric inv_y = 0.0, inv_x = 0.0;  // Init'd to remove warnings
        
        // X and Y cannot both be zero
        if(x == 0.0) {
          inv_y = 1.0 / y;
          C0 = 1.0;
          C1 = 1.0;
          C2 = (1.0 - cos_y) * inv_x2y2;
          C3 = (1.0 - sin_y*inv_y) * inv_x2y2;
        }
        else if(y == 0.0) {
          inv_x = 1.0 / x;
          C0 = 1.0;
          C1 = 1.0;
          C2 = (cosh_x - 1.0) * inv_x2y2;
          C3 = (sinh_x*inv_x - 1.0) * inv_x2y2;
        }
        else {
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
    }
    T.set4(F*exp_a, i);
  }
}

void dtransmat1(TransmissionMatrix& T,
                ArrayOfTransmissionMatrix& dT1,
                ArrayOfTransmissionMatrix& dT2,
                const PropagationMatrix& K1,
                const PropagationMatrix& K2,
                const ArrayOfPropagationMatrix& dK1,
                const ArrayOfPropagationMatrix& dK2,
                const Numeric& r,
                const Numeric& dr_dT1,
                const Numeric& dr_dT2,
                const Index it,
                const Index iz,
                const Index ia)
{
  for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
    T.set1(std::exp(-0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])), i);
    for(Index j=0; j<dT1.nelem(); j++) {
      if(dK1[j].NumberOfFrequencies())
        dT1[j].set1(T.Mat1(i) * (-0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
      if(dK2[j].NumberOfFrequencies())
        dT2[j].set1(T.Mat1(i) * (-0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
    }
  }
}


void dtransmat2(TransmissionMatrix& T,
                ArrayOfTransmissionMatrix& dT1,
                ArrayOfTransmissionMatrix& dT2,
                const PropagationMatrix& K1,
                const PropagationMatrix& K2,
                const ArrayOfPropagationMatrix& dK1,
                const ArrayOfPropagationMatrix& dK2,
                const Numeric& r,
                const Numeric& dr_dT1,
                const Numeric& dr_dT2,
                const Index it,
                const Index iz,
                const Index ia)
{
  for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]), 
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    Eigen::Matrix2d F;
    
    if(b == 0.) {
      F.setIdentity();
      T.set2(F * exp_a, i);
      for(Index j=0; j<dT1.nelem(); j++) {
        if(dK1[j].NumberOfFrequencies())
          dT1[j].set2(T.Mat2(i) * (-0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
        if(dK2[j].NumberOfFrequencies())
          dT2[j].set2(T.Mat2(i) * (-0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
      }
    }
    else {
      const Numeric C0 = std::cosh(b) - a*std::sinh(b)/b;
      const Numeric C1 = std::sinh(b)/b;  
      
      F(0, 0) = F(1, 1) = C0 + C1 * a;
      F(0, 1) = F(1, 0) = C1 * b;
      T.set2(F * exp_a, i);
      
      for(Index j=0; j<dT2.nelem(); j++) {
        if(not dK2[j].NumberOfFrequencies())
          continue;
        const Numeric da = -0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]) : 0.0)),
                      db = -0.5 * (r * dK2[j].K12(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]) : 0.0));
        const Numeric dC0 = -a*std::cosh(b)*db/b + a*std::sinh(b)*db/b/b + std::sinh(b)*db - std::sinh(b)*da/b;
        const Numeric dC1 = (std::cosh(b) - C1)*db/b;
        Eigen::Matrix2d dF;
        dF(0, 0) = dF(1, 1) = dC0 + C1 * da + dC1 * a + F(0, 0) * da;
        dF(0, 1) = dF(1, 0) = C1 * db + dC1 * b + F(0, 1) * da;
        dT2[j].set2(std::move(dF), i);
      }
      for(Index j=0; j<dT1.nelem(); j++) {
        if(not dK1[j].NumberOfFrequencies())
          continue;
        const Numeric da = -0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]) : 0.0)),
                      db = -0.5 * (r * dK1[j].K12(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]) : 0.0));
        const Numeric dC0 = -a*std::cosh(b)*db/b + a*std::sinh(b)*db/b/b + std::sinh(b)*db - std::sinh(b)*da/b;
        const Numeric dC1 = (std::cosh(b) - C1)*db/b;
        Eigen::Matrix2d dF;
        dF(0, 0) = dF(1, 1) = dC0 + C1 * da + dC1 * a + F(0, 0) * da;
        dF(0, 1) = dF(1, 0) = C1 * db + dC1 * b + F(0, 1) * da;
        dT1[j].set2(std::move(dF), i);
      }
    }
  }
}


void dtransmat3(TransmissionMatrix& T,
                ArrayOfTransmissionMatrix& dT1,
                ArrayOfTransmissionMatrix& dT2,
                const PropagationMatrix& K1,
                const PropagationMatrix& K2,
                const ArrayOfPropagationMatrix& dK1,
                const ArrayOfPropagationMatrix& dK2,
                const Numeric& r,
                const Numeric& dr_dT1,
                const Numeric& dr_dT2,
                const Index it,
                const Index iz,
                const Index ia)
{
  for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]), 
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]), 
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]), 
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    Eigen::Matrix3d F;
    
    if(b == 0. and c == 0. and u == 0.) {
      F.setIdentity();
      T.set3(F*exp_a, i);
      for(Index j=0; j<dT1.nelem(); j++) {
        if(dK1[j].NumberOfFrequencies())
          dT1[j].set3(T.Mat3(i) * (-0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
        if(dK2[j].NumberOfFrequencies())
          dT2[j].set3(T.Mat3(i) * (-0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
      }
    }
    else {
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      
      const Numeric x = std::sqrt(b2 + c2 - u2), x2 = x * x, inv_x2 = 1.0/x2;
      const Numeric sinh_x = std::sinh(x), cosh_x = std::cosh(x);
      
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
      T.set3(F * exp_a, i);
      
      for(Index j=0; j<dK2.nelem(); j++) {
        if(not dK2[j].NumberOfFrequencies())
          continue;
        const Numeric da = -0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]) : 0.0)),
                      db = -0.5 * (r * dK2[j].K12(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]) : 0.0)),
                      dc = -0.5 * (r * dK2[j].K13(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]) : 0.0)),
                      du = -0.5 * (r * dK2[j].K23(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]) : 0.0));
        
        const Numeric da2 = 2 * a * da, db2 = 2 * b * db, dc2 = 2 * c * dc, du2 = 2 * u * du;
        
        const Numeric dx  = (db2 + dc2 -du2)/x/2;
        
        const Numeric dC0 = -2 * (C0 - 1) * dx/x + (da2 * (cosh_x - 1) + a2 * sinh_x*dx - a * b * cosh_x * dx - a * sinh_x * dx - b * sinh_x * da) * inv_x2;
        const Numeric dC1 = -2 * C1 * dx / x + (2 * da * (1 - cosh_x) - 2 * a * sinh_x * dx + x*cosh_x*dx + sinh_x*dx) * inv_x2;
        const Numeric dC2 = (sinh_x/x - 2*C2)*dx/x;
        
        Eigen::Matrix3d dF;
        
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
        dT2[j].set3(std::move(dF), i);
      }
      
      for(Index j=0; j<dK1.nelem(); j++) {
        if(not dK1[j].NumberOfFrequencies())
          continue;
        const Numeric da = -0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]) : 0.0)),
                      db = -0.5 * (r * dK1[j].K12(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]) : 0.0)),
                      dc = -0.5 * (r * dK1[j].K13(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]) : 0.0)),
                      du = -0.5 * (r * dK1[j].K23(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]) : 0.0));
        
        const Numeric da2 = 2 * a * da, db2 = 2 * b * db, dc2 = 2 * c * dc, du2 = 2 * u * du;
        
        const Numeric dx  = (db2 + dc2 -du2)/x/2;
        
        const Numeric dC0 = -2 * (C0 - 1) * dx/x + (da2 * (cosh_x - 1) + a2 * sinh_x*dx - a * b * cosh_x * dx - a * sinh_x * dx - b * sinh_x * da) * inv_x2;
        const Numeric dC1 = -2 * C1 * dx / x + (2 * da * (1 - cosh_x) - 2 * a * sinh_x * dx + x*cosh_x*dx + sinh_x*dx) * inv_x2;
        const Numeric dC2 = (sinh_x/x - 2*C2)*dx/x;
        
        Eigen::Matrix3d dF;
        
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
        dT1[j].set3(std::move(dF), i);
      }
    }
  }
}


void dtransmat4(TransmissionMatrix& T,
                ArrayOfTransmissionMatrix& dT1,
                ArrayOfTransmissionMatrix& dT2,
                const PropagationMatrix& K1,
                const PropagationMatrix& K2,
                const ArrayOfPropagationMatrix& dK1,
                const ArrayOfPropagationMatrix& dK2,
                const Numeric& r,
                const Numeric& dr_dT1,
                const Numeric& dr_dT2,
                const Index it,
                const Index iz,
                const Index ia)
{
  static const Numeric sqrt_05 = sqrt(0.5);
  for(Index i=0; i<K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]), 
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]), 
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]), 
                  d = -0.5 * r * (K1.K14(iz, ia)[i] + K2.K14(iz, ia)[i]), 
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]), 
                  v = -0.5 * r * (K1.K24(iz, ia)[i] + K2.K24(iz, ia)[i]), 
                  w = -0.5 * r * (K1.K34(iz, ia)[i] + K2.K34(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    Eigen::Matrix4d F;

    if(b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.) {
      F.setIdentity();
      T.set4(F * exp_a, i);
      for(Index j=0; j<dK1.nelem(); j++) {
        if(dK1[j].NumberOfFrequencies())
          dT1[j].set4(T.Mat4(i) * (-0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
        if(dK2[j].NumberOfFrequencies())
          dT1[j].set4(T.Mat4(i) * (-0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it)?dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]):0.0))), i);
      }
    }
    else {
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
        Const1 = std::sqrt(Const1);
      else
        Const1 = 0.0;
      
      const Complex sqrt_BpA = std::sqrt(Complex(Const2 + Const1, 0.0));
      const Complex sqrt_BmA = std::sqrt(Complex(Const2 - Const1, 0.0));
      const Numeric x = sqrt_BpA.real() * sqrt_05;
      const Numeric y = sqrt_BmA.imag() * sqrt_05;
      const Numeric x2 = x * x;
      const Numeric y2 = y * y;
      const Numeric cos_y = std::cos(y);
      const Numeric sin_y = std::sin(y);
      const Numeric cosh_x = std::cosh(x);
      const Numeric sinh_x = std::sinh(x);
      const Numeric x2y2 = x2 + y2;
      const Numeric inv_x2y2 = 1.0 / x2y2;
      
      Numeric C0, C1, C2, C3;
      Numeric inv_y=0.0, inv_x=0.0;  // Init'd to remove warnings
      
      // X and Y cannot both be zero
      if(x == 0.0) {
        inv_y = 1.0 / y;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (1.0 - cos_y) * inv_x2y2;
        C3 = (1.0 - sin_y*inv_y) * inv_x2y2;
      }
      else if(y == 0.0) {
        inv_x = 1.0 / x;
        C0 = 1.0;
        C1 = 1.0;
        C2 = (cosh_x - 1.0) * inv_x2y2;
        C3 = (sinh_x*inv_x -1.0) * inv_x2y2;
      }
      else {
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
      
      T.set4(F * exp_a, i);
      
      if(dK1.nelem()) {
        const Numeric inv_x2 = inv_x * inv_x;
        const Numeric inv_y2 = inv_y * inv_y;
        
        for(Index j=0; j<dK1.nelem(); j++) {
          if(not dK1[j].NumberOfFrequencies())
            continue;
          const Numeric da = -0.5 * (r * dK1[j].Kjj(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]) : 0.0)),
                        db = -0.5 * (r * dK1[j].K12(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]) : 0.0)),
                        dc = -0.5 * (r * dK1[j].K13(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]) : 0.0)),
                        dd = -0.5 * (r * dK1[j].K14(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K14(iz, ia)[i] + K2.K14(iz, ia)[i]) : 0.0)),
                        du = -0.5 * (r * dK1[j].K23(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]) : 0.0)),
                        dv = -0.5 * (r * dK1[j].K24(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K24(iz, ia)[i] + K2.K24(iz, ia)[i]) : 0.0)),
                        dw = -0.5 * (r * dK1[j].K34(iz, ia)[i] + ((j==it) ? dr_dT1 * (K1.K34(iz, ia)[i] + K2.K34(iz, ia)[i]) : 0.0));
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c,
                        dd2 = 2 * dd * d, du2 = 2 * du * u,
                        dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          
          Numeric dConst1;
          if(Const1 > 0.) {
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
          if(x == 0.0) {
            const Numeric dy = (0.5 * (dConst2 - dConst1)/sqrt_BmA).imag() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*y*dy * C2 * inv_x2y2 + dy*sin_y * inv_x2y2;
            dC3 = -2*y*dy * C3 * inv_x2y2 + (dy*sin_y*inv_y2 - cos_y*dy*inv_y) * inv_x2y2;;
          }
          else if(y == 0.0) {
            const Numeric dx = (0.5 * (dConst2 + dConst1)/sqrt_BpA).real() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*x*dx * C2 * inv_x2y2 + dx*sinh_x * inv_x2y2;
            dC3 = -2*x*dx * C3 * inv_x2y2 + (cosh_x*dx*inv_x - dx*sinh_x*inv_x2) * inv_x2y2;
          }
          else { 
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
          
          Eigen::Matrix4d dF;
          
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
          
          dT1[j].set4(dF * exp_a + T.Mat4(i) * da, i);
        }
        
        for(Index j=0; j<dK2.nelem(); j++) {
          if(not dK2[j].NumberOfFrequencies())
            continue;
          const Numeric da = -0.5 * (r * dK2[j].Kjj(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]) : 0.0)),
                        db = -0.5 * (r * dK2[j].K12(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]) : 0.0)),
                        dc = -0.5 * (r * dK2[j].K13(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]) : 0.0)),
                        dd = -0.5 * (r * dK2[j].K14(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K14(iz, ia)[i] + K2.K14(iz, ia)[i]) : 0.0)),
                        du = -0.5 * (r * dK2[j].K23(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]) : 0.0)),
                        dv = -0.5 * (r * dK2[j].K24(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K24(iz, ia)[i] + K2.K24(iz, ia)[i]) : 0.0)),
                        dw = -0.5 * (r * dK2[j].K34(iz, ia)[i] + ((j==it) ? dr_dT2 * (K1.K34(iz, ia)[i] + K2.K34(iz, ia)[i]) : 0.0));
          
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c,
                        dd2 = 2 * dd * d, du2 = 2 * du * u,
                        dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          
          Numeric dConst1;
          if(Const1 > 0.) {
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
          if(x == 0.0) {
            const Numeric dy = (0.5 * (dConst2 - dConst1)/sqrt_BmA).imag() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*y*dy * C2 * inv_x2y2 + dy*sin_y * inv_x2y2;
            dC3 = -2*y*dy * C3 * inv_x2y2 + (dy*sin_y*inv_y2 - cos_y*dy*inv_y) * inv_x2y2;;
          }
          else if(y == 0.0) {
            const Numeric dx = (0.5 * (dConst2 + dConst1)/sqrt_BpA).real() * sqrt_05;
            
            dC0 = 0.0;
            dC1 = 0.0;
            dC2 = -2*x*dx * C2 * inv_x2y2 + dx*sinh_x * inv_x2y2;
            dC3 = -2*x*dx * C3 * inv_x2y2 + (cosh_x*dx*inv_x - dx*sinh_x*inv_x2) * inv_x2y2;
          }
          else { 
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
          
          Eigen::Matrix4d dF;
          
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
          
          dT2[j].set4(dF * exp_a + T.Mat4(i) * da, i);
        }
      }
    }
  }
}


void transmat(TransmissionMatrix& T,
              const PropagationMatrix& K1,
              const PropagationMatrix& K2,
              const Numeric& r)
{
  switch(K1.StokesDimensions()) {
    case 4: transmat4(T, K1, K2, r); break;
    case 3: transmat3(T, K1, K2, r); break;
    case 2: transmat2(T, K1, K2, r); break;
    case 1: transmat1(T, K1, K2, r); break;
  }
}


void dtransmat(TransmissionMatrix& T,
               ArrayOfTransmissionMatrix& dT1,
               ArrayOfTransmissionMatrix& dT2,
               const PropagationMatrix& K1,
               const PropagationMatrix& K2,
               const ArrayOfPropagationMatrix& dK1,
               const ArrayOfPropagationMatrix& dK2,
               const Numeric& r,
               const Numeric& dr_dT1=0,
               const Numeric& dr_dT2=0,
               const Index it=-1,
               const Index iz=0,
               const Index ia=0)
{
  switch(K1.StokesDimensions()) {
    case 4: dtransmat4(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia); break;
    case 3: dtransmat3(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia); break;
    case 2: dtransmat2(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia); break;
    case 1: dtransmat1(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia); break;
  }
}



void stepwise_transmission(TransmissionMatrix& PiT,
                           TransmissionMatrix& T,
                           ArrayOfTransmissionMatrix& dT1,
                           ArrayOfTransmissionMatrix& dT2,
                           const TransmissionMatrix& PiT_last,
                           const PropagationMatrix& K1,
                           const PropagationMatrix& K2,
                           const ArrayOfPropagationMatrix& dK1,
                           const ArrayOfPropagationMatrix& dK2,
                           const Numeric& r,
                           const bool& first)
{
  if(first)
    PiT.setIdentity();
  else {
    if(not dT1.nelem())
      transmat(T, K1, K2, r);
    else
      dtransmat(T, dT1, dT2, K1, K2, dK1, dK2, r);  // FIXME:  Add path-derivatives here [Internally, only for temperature for now; expressions valid for other things]
    PiT.mul(PiT_last, T);
  }
}


void stepwise_source(RadiationVector& J,
                     ArrayOfRadiationVector& dJ,
                     const PropagationMatrix& K,
                     const StokesVector& a,
                     const StokesVector& S,
                     const ArrayOfPropagationMatrix& dK,
                     const ArrayOfStokesVector& da,
                     const ArrayOfStokesVector& dS,
                     const ConstVectorView B,
                     const ConstVectorView dB_dT,
                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                     const bool& jacobian_do)
{
  for(Index i=0; i<K.NumberOfFrequencies(); i++) {
    if(K.IsRotational(i)) {
      J.set(0, i);
      for(Index j=0; j<jacobian_quantities.nelem(); j++)
        if(jacobian_quantities[j].Analytical())
          dJ[j].set(0, i);
    }
    else {
      J.setSource(a, B, S, i);
      switch(J.StokesDim()) {
        case 4: {
          const auto invK = inv4(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K14()[i],
                                 K.K23()[i], K.K24()[i], K.K34()[i]);
          J.leftMul4(invK, i);
          if(jacobian_do) {
            for(Index j=0; j<jacobian_quantities.nelem(); j++) {
              if(jacobian_quantities[j].Analytical()) {
                dJ[j].set(J, i);
                dJ[j].leftMul4(-matrix4(dK[j].Kjj()[i], dK[j].K12()[i], dK[j].K13()[i], dK[j].K14()[i],
                                        dK[j].K23()[i], dK[j].K24()[i], dK[j].K34()[i]), i);
                dJ[j].add4(vector4(a, B, da[j], dB_dT, dS[j], jacobian_quantities[j].IsTemperature(), i), i);
                dJ[j].leftMul4(0.5*invK, i);
              }
            }
          }
        } break;
        case 3: {
          const auto invK = inv3(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K23()[i]);
          J.leftMul3(invK, i);
          if(jacobian_do) {
            for(Index j=0; j<jacobian_quantities.nelem(); j++) {
              if(jacobian_quantities[j].Analytical()) {
                dJ[j].set(J, i);
                dJ[j].leftMul3(-matrix3(dK[j].Kjj()[i], dK[j].K12()[i], dK[j].K13()[i],dK[j].K23()[i]), i);
                dJ[j].add3(vector3(a, B, da[j], dB_dT, dS[j], jacobian_quantities[j].IsTemperature(), i), i);
                dJ[j].leftMul3(0.5*invK, i);
              }
            }
          }
        } break;
        case 2: {
          const auto invK = inv2(K.Kjj()[i], K.K12()[i]);
          J.leftMul2(invK, i);
          if(jacobian_do) {
            for(Index j=0; j<jacobian_quantities.nelem(); j++) {
              if(jacobian_quantities[j].Analytical()) {
                dJ[j].set(J, i);
                dJ[j].leftMul2(-matrix2(dK[j].Kjj()[i], dK[j].K12()[i]), i);
                dJ[j].add2(vector2(a, B, da[j], dB_dT, dS[j], jacobian_quantities[j].IsTemperature(), i), i);
                dJ[j].leftMul2(0.5*invK, i);
              }
            }
          }
        } break;
        default: {
          const auto invK = inv1(K.Kjj()[i]);
          J.leftMul1(invK, i);
          if(jacobian_do) {
            for(Index j=0; j<jacobian_quantities.nelem(); j++) {
              if(jacobian_quantities[j].Analytical()) {
                dJ[j].set(J, i);
                dJ[j].leftMul1(-matrix1(dK[j].Kjj()[i]), i);
                dJ[j].add1(vector1(a, B, da[j], dB_dT, dS[j], jacobian_quantities[j].IsTemperature(), i), i);
                dJ[j].leftMul1(0.5*invK, i);
              }
            }
          }
        } break;
      }
    }
  }
}


void update_radiation_vector(RadiationVector& I,
                             ArrayOfRadiationVector& dI1,
                             ArrayOfRadiationVector& dI2,
                             const RadiationVector& J1,
                             const RadiationVector& J2,
                             const ArrayOfRadiationVector& dJ1,
                             const ArrayOfRadiationVector& dJ2,
                             const TransmissionMatrix& T,
                             const TransmissionMatrix& PiT,
                             const ArrayOfTransmissionMatrix& dT1,
                             const ArrayOfTransmissionMatrix& dT2)
{
  I.rem_avg(J1, J2);
  
  for(size_t i=0; i<dI1.size(); i++) {
    dI1[i].setDeriv(PiT, dT1[i], T, I, dJ1[i]);
    dI2[i].addDeriv(PiT, dT2[i], T, I, dJ2[i]);
  }
  
  I.leftMul(T);
  I.add_avg(J1, J2);
}


std::ostream& operator<<(std::ostream& os, const TransmissionMatrix& tm) {
  for(const auto& T: tm.T4) os << T << '\n';
  for(const auto& T: tm.T3) os << T << '\n';
  for(const auto& T: tm.T2) os << T << '\n';
  for(const auto& T: tm.T1) os << T << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfTransmissionMatrix& atm) {
  for(const auto& T: atm) os << T << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfTransmissionMatrix& aatm) {
  for(const auto& T: aatm) os << T << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const RadiationVector& rv) {
  // Write the transpose because it looks better...
  for(const auto& R: rv.R4) os << R.transpose() << '\n';
  for(const auto& R: rv.R3) os << R.transpose() << '\n';
  for(const auto& R: rv.R2) os << R.transpose() << '\n';
  for(const auto& R: rv.R1) os << R.transpose() << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfRadiationVector& arv) {
  for(const auto& R: arv) os << R << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfRadiationVector& aarv) {
  for(const auto& R: aarv) os << R << '\n';
  return os;
}

std::istream& operator>>(std::istream& data, TransmissionMatrix& tm) {
  for(auto& T: tm.T4) data >> T(0, 0) >> T(0, 1) >> T(0, 2) >> T(0, 3) 
                           >> T(1, 0) >> T(1, 1) >> T(1, 2) >> T(1, 3)
                           >> T(2, 0) >> T(2, 1) >> T(2, 2) >> T(2, 3)
                           >> T(3, 0) >> T(3, 1) >> T(3, 2) >> T(3, 3);
  for(auto& T: tm.T3) data >> T(0, 0) >> T(0, 1) >> T(0, 2)
                           >> T(1, 0) >> T(1, 1) >> T(1, 2)
                           >> T(2, 0) >> T(2, 1) >> T(2, 2);
  for(auto& T: tm.T2) data >> T(0, 0) >> T(0, 1)
                           >> T(1, 0) >> T(1, 1);
  for(auto& T: tm.T1) data >> T(0, 0);
  return data;
}

std::istream& operator>>(std::istream& data, RadiationVector& rv) {
  for(auto& R: rv.R4) data >> R[0] >> R[1] >> R[2] >> R[3];
  for(auto& R: rv.R3) data >> R[0] >> R[1] >> R[2];
  for(auto& R: rv.R2) data >> R[0] >> R[1] >> R[2];
  for(auto& R: rv.R1) data >> R[0];
  return data;
}
