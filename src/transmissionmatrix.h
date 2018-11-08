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
 * \file   transmissionmatrix.h
 * \brief  Stuff related to the transmission matrix.
 * 
 * Using Eigen library to speed up computations.
 * 
 * \author Richard Larsson
 * \date   2018-01-30
 */

#ifndef transmissionmatrix_h
#define transmissionmatrix_h

#include "jacobian.h"
#include "propagationmatrix.h"
#include <Eigen/Dense>
#include <Eigen/StdVector>

Eigen::Matrix<double, 1, 1> inv1(const Numeric& x);
Eigen::Matrix2d inv2(const Numeric& a, const Numeric& b);
Eigen::Matrix3d inv3(const Numeric& a, const Numeric& b, const Numeric& c, const Numeric& u);
Eigen::Matrix4d inv4(const Numeric& a, const Numeric& b, const Numeric& c, const Numeric& d, const Numeric& u, const Numeric& v, const Numeric& w);

class TransmissionMatrix {
private:
  Index stokes_dim;
  std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>> T4;
  std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> T3;
  std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> T2;
  std::vector<Eigen::Matrix<double, 1, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 1, 1>>> T1;
public:
  TransmissionMatrix(Index nf=0, Index stokes=1) : stokes_dim(stokes),
  T4(stokes_dim==4?nf:0, Eigen::Matrix4d::Identity()), 
  T3(stokes_dim==3?nf:0, Eigen::Matrix3d::Identity()), 
  T2(stokes_dim==2?nf:0, Eigen::Matrix2d::Identity()), 
  T1(stokes_dim==1?nf:0, Eigen::Matrix<double, 1, 1>::Identity())
  { assert(stokes_dim < 5 and stokes_dim > 0); }
  
  operator Tensor3() const {
    Tensor3 T;
    switch(stokes_dim) {
      case 4:
        T.resize(T4.size(), 4, 4);
        for(size_t i=0; i<T4.size(); i++)
          for(size_t j=0; j<4; j++)
            for(size_t k=0; k<4; k++)
              T(i, j, k) = T4[i](j, k);
        break;
      case 3:
        T.resize(T3.size(), 3, 3);
        for(size_t i=0; i<T3.size(); i++)
          for(size_t j=0; j<3; j++)
            for(size_t k=0; k<3; k++)
              T(i, j, k) = T3[i](j, k);
        break;
      case 2:
        T.resize(T2.size(), 2, 2);
        for(size_t i=0; i<T2.size(); i++)
          for(size_t j=0; j<2; j++)
            for(size_t k=0; k<2; k++)
              T(i, j, k) = T2[i](j, k);
        break;
      case 1:
        T.resize(T1.size(), 1, 1);
        for(size_t i=0; i<T1.size(); i++)
          T(i, 0, 0) = T1[i](0, 0);
        break;
    }
    return T;
  }
  
  const Eigen::Matrix4d& Mat4(size_t i) const {return T4[i];}
  const Eigen::Matrix3d& Mat3(size_t i) const {return T3[i];}
  const Eigen::Matrix2d& Mat2(size_t i) const {return T2[i];}
  const Eigen::Matrix<double, 1, 1>& Mat1(size_t i) const {return T1[i];}
  
  Eigen::Matrix4d& Mat4(size_t i) {return T4[i];}
  Eigen::Matrix3d& Mat3(size_t i) {return T3[i];}
  Eigen::Matrix2d& Mat2(size_t i) {return T2[i];}
  Eigen::Matrix<double, 1, 1>& Mat1(size_t i) {return T1[i];}
  
  void setIdentity() {
    for(size_t i=0; i<T4.size(); i++) T4[i]=Eigen::Matrix4d::Identity();
    for(size_t i=0; i<T3.size(); i++) T3[i]=Eigen::Matrix3d::Identity();
    for(size_t i=0; i<T2.size(); i++) T2[i]=Eigen::Matrix2d::Identity();
    for(size_t i=0; i<T1.size(); i++) T1[i](0, 0)=1;
  }
  
  void mul(const TransmissionMatrix& A, const TransmissionMatrix& B) {
    for(size_t i=0; i<T4.size(); i++) T4[i].noalias() = A.T4[i] * B.T4[i];
    for(size_t i=0; i<T3.size(); i++) T3[i].noalias() = A.T3[i] * B.T3[i];
    for(size_t i=0; i<T2.size(); i++) T2[i].noalias() = A.T2[i] * B.T2[i];
    for(size_t i=0; i<T1.size(); i++) T1[i].noalias() = A.T1[i] * B.T1[i];
  }
  
  void setZero4(Index i) {T4[i].setZero();}
  void setZero3(Index i) {T3[i].setZero();}
  void setZero2(Index i) {T2[i].setZero();}
  void setZero1(Index i) {T1[i].setZero();}
  
  void set4(Eigen::Matrix4d&& O, Index i) {T4[i]=std::move(O);}
  void set4(const Eigen::Matrix4d& O, Index i) {T4[i]=O;}
  
  void set3(Eigen::Matrix3d&& O, Index i) {T3[i]=std::move(O);}
  void set3(const Eigen::Matrix3d& O, Index i) {T3[i]=O;}
  
  void set2(Eigen::Matrix2d&& O, Index i) {T2[i]=std::move(O);}
  void set2(const Eigen::Matrix2d& O, Index i) {T2[i]=O;}
  
  void set1(Eigen::Matrix<double, 1, 1>&& O, Index i) {T1[i]=std::move(O);}
  void set1(const Eigen::Matrix<double, 1, 1>& O, Index i) {T1[i]=O;}
  void set1(const Numeric& O, Index i) {T1[i][0]=O;}
  
  const Numeric& operator()(const Index i, const Index j, const Index k) const {
    switch(stokes_dim) {
      case 4: return T4[i](j, k);
      case 3: return T3[i](j, k);
      case 2: return T2[i](j, k);
      default: return T1[i](j, k);
    }
  }
};


class RadiationVector {
private:
  Index stokes_dim;
  std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d>> R4;
  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> R3;
  std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> R2;
  std::vector<Eigen::Matrix<double, 1, 1>, Eigen::aligned_allocator<Eigen::Matrix<double, 1, 1>>> R1;
public:
  RadiationVector(Index nf=0, Index stokes=1) : stokes_dim(stokes),
  R4(stokes_dim==4?nf:0, Eigen::Vector4d::Zero()), 
  R3(stokes_dim==3?nf:0, Eigen::Vector3d::Zero()), 
  R2(stokes_dim==2?nf:0, Eigen::Vector2d::Zero()), 
  R1(stokes_dim==1?nf:0, Eigen::Matrix<double, 1, 1>::Zero())
  { assert(stokes_dim < 5 and stokes_dim > 0); }
  
  void leftMul(const TransmissionMatrix& T) {
    for(size_t i=0; i<R4.size(); i++) R4[i] = T.Mat4(i) * R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i] = T.Mat3(i) * R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i] = T.Mat2(i) * R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i] = T.Mat1(i) * R1[i];
  }
  
  void leftMul(const Eigen::MatrixXd& X, Index i) {
    switch(stokes_dim) {
      case 4: leftMul4(X, i); break;
      case 3: leftMul3(X, i); break;
      case 2: leftMul2(X, i); break;
      default: leftMul1(X, i); break;
    }
  }
  
  void leftMul4(const Eigen::Matrix4d& X, Index i) {R4[i] = X * R4[i];}
  void leftMul3(const Eigen::Matrix3d& X, Index i) {R3[i] = X * R3[i];}
  void leftMul2(const Eigen::Matrix2d& X, Index i) {R2[i] = X * R2[i];}
  void leftMul1(const Eigen::Matrix<double, 1, 1>& X, Index i) {R1[i] = X * R1[i];}
  
  RadiationVector& operator-=(const RadiationVector& O) {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() -= O.R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() -= O.R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() -= O.R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() -= O.R1[i];
    return *this;
  }
  
  RadiationVector& operator+=(const RadiationVector& O) {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() += O.R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() += O.R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() += O.R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() += O.R1[i];
    return *this;
  }
  
  RadiationVector& operator/=(const Numeric& O) {
    for(size_t i=0; i<R4.size(); i++) R4[i] /= O;
    for(size_t i=0; i<R3.size(); i++) R3[i] /= O;
    for(size_t i=0; i<R2.size(); i++) R2[i] /= O;
    for(size_t i=0; i<R1.size(); i++) R1[i] /= O;
    return *this;
  }
  
  RadiationVector& operator*=(const Numeric& O) {
    for(size_t i=0; i<R4.size(); i++) R4[i] *= O;
    for(size_t i=0; i<R3.size(); i++) R3[i] *= O;
    for(size_t i=0; i<R2.size(); i++) R2[i] *= O;
    for(size_t i=0; i<R1.size(); i++) R1[i] *= O;
    return *this;
  }
  
  RadiationVector& operator=(const StokesVector& s) {
    assert(s.NumberOfAzimuthAngles() == 1);
    assert(s.NumberOfZenithAngles() == 1);
    for(size_t i=0; i<R4.size(); i++) R4[i] = {s.Kjj()[i], s.K12()[i], s.K13()[i], s.K14()[i]};
    for(size_t i=0; i<R3.size(); i++) R3[i] = {s.Kjj()[i], s.K12()[i], s.K13()[i]};
    for(size_t i=0; i<R2.size(); i++) R2[i] = {s.Kjj()[i], s.K12()[i]};
    for(size_t i=0; i<R1.size(); i++) R1[i][0] = s.Kjj()[i];
    return *this;
  }
  
  RadiationVector& operator*=(const ConstVectorView& v) {
    for(size_t i=0; i<R4.size(); i++) R4[i] *= v[i];
    for(size_t i=0; i<R3.size(); i++) R3[i] *= v[i];
    for(size_t i=0; i<R2.size(); i++) R2[i] *= v[i];
    for(size_t i=0; i<R1.size(); i++) R1[i][0] *= v[i];
    return *this;
  }
  
  RadiationVector& operator+=(const StokesVector& s) {
    assert(s.NumberOfAzimuthAngles() == 1);
    assert(s.NumberOfZenithAngles() == 1);
    for(size_t i=0; i<R4.size(); i++) R4[i] += Eigen::Vector4d(s.Kjj()[i], s.K12()[i], s.K13()[i], s.K14()[i]);
    for(size_t i=0; i<R3.size(); i++) R3[i] += Eigen::Vector3d(s.Kjj()[i], s.K12()[i], s.K13()[i]);
    for(size_t i=0; i<R2.size(); i++) R2[i] += Eigen::Vector2d(s.Kjj()[i], s.K12()[i]);
    for(size_t i=0; i<R1.size(); i++) R1[i][0] += s.Kjj()[i];
    return *this;
  }
  
  void set(const RadiationVector& O, Index i) {
    switch(stokes_dim) {
      case 4: R4[i] = O.R4[i]; break;
      case 3: R3[i] = O.R3[i]; break;
      case 2: R2[i] = O.R2[i]; break;
      default: R1[i] = O.R1[i]; break;
    }
  }
  
  void set(const Numeric& O, Index i) {
    switch(stokes_dim) {
      case 4: R4[i] = Eigen::Vector4d::Constant(O); break;
      case 3: R3[i] = Eigen::Vector3d::Constant(O); break;
      case 2: R2[i] = Eigen::Vector2d::Constant(O); break;
      default: R1[i][0] = O; break;
    }
  }
  
  void add(const Eigen::MatrixXd& O, Index i) {
    switch(stokes_dim) {
      case 4: add4(O, i); break;
      case 3: add3(O, i); break;
      case 2: add2(O, i); break;
      default: add1(O, i); break;
    }
  }
  
  void add4(const Eigen::Vector4d& O, Index i) {R4[i] += O;}
  void add3(const Eigen::Vector3d& O, Index i) {R3[i] += O;}
  void add2(const Eigen::Vector2d& O, Index i) {R2[i] += O;}
  void add1(const Eigen::Matrix<double, 1, 1>& O, Index i) {R1[i] += O;}
  
  void rem_avg(const RadiationVector& O1, const RadiationVector& O2) {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() -= 0.5 * (O1.R4[i] + O2.R4[i]);
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() -= 0.5 * (O1.R3[i] + O2.R3[i]);
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() -= 0.5 * (O1.R2[i] + O2.R2[i]);
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() -= 0.5 * (O1.R1[i] + O2.R1[i]);
  }
  
  void add_avg(const RadiationVector& O1, const RadiationVector& O2) {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() += 0.5 * (O1.R4[i] + O2.R4[i]);
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() += 0.5 * (O1.R3[i] + O2.R3[i]);
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() += 0.5 * (O1.R2[i] + O2.R2[i]);
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() += 0.5 * (O1.R1[i] + O2.R1[i]);
  }
  
  void setDeriv(const TransmissionMatrix& PiT,
                const TransmissionMatrix& dT,
                const TransmissionMatrix& T,
                const RadiationVector& ImJ,
                const RadiationVector& dJ)
  {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() = PiT.Mat4(i) * (dT.Mat4(i) * ImJ.R4[i] + dJ.R4[i] - T.Mat4(i) * dJ.R4[i]);
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() = PiT.Mat3(i) * (dT.Mat3(i) * ImJ.R3[i] + dJ.R3[i] - T.Mat3(i) * dJ.R3[i]);
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() = PiT.Mat2(i) * (dT.Mat2(i) * ImJ.R2[i] + dJ.R2[i] - T.Mat2(i) * dJ.R2[i]);
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() = PiT.Mat1(i) * (dT.Mat1(i) * ImJ.R1[i] + dJ.R1[i] - T.Mat1(i) * dJ.R1[i]);
  }
  
  void addDeriv(const TransmissionMatrix& PiT,
                const TransmissionMatrix& dT,
                const TransmissionMatrix& T,
                const RadiationVector& ImJ,
                const RadiationVector& dJ)
  {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() += PiT.Mat4(i) * (dT.Mat4(i) * ImJ.R4[i] + dJ.R4[i] - T.Mat4(i) * dJ.R4[i]);
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() += PiT.Mat3(i) * (dT.Mat3(i) * ImJ.R3[i] + dJ.R3[i] - T.Mat3(i) * dJ.R3[i]);
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() += PiT.Mat2(i) * (dT.Mat2(i) * ImJ.R2[i] + dJ.R2[i] - T.Mat2(i) * dJ.R2[i]);
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() += PiT.Mat1(i) * (dT.Mat1(i) * ImJ.R1[i] + dJ.R1[i] - T.Mat1(i) * dJ.R1[i]);
  }
  
  RadiationVector& operator=(const ConstMatrixView& M) {
    for(size_t i=0; i<R4.size(); i++) {R4[i][0] = M(i, 0); R4[i][1] = M(i, 1); R4[i][2] = M(i, 2); R4[i][3] = M(i, 3);}
    for(size_t i=0; i<R3.size(); i++) {R3[i][0] = M(i, 0); R3[i][1] = M(i, 1); R3[i][2] = M(i, 2);}
    for(size_t i=0; i<R2.size(); i++) {R2[i][0] = M(i, 0); R2[i][1] = M(i, 1);}
    for(size_t i=0; i<R1.size(); i++) {R1[i][0] = M(i, 0);}
    return *this;
  }
  
  const Numeric& operator()(const Index i, const Index j) const {
    switch(stokes_dim) {
      case 4: return R4[i][j];
      case 3: return R3[i][j];
      case 2: return R2[i][j];
      default: return R1[i][j];
    }
  }
  
  RadiationVector& operator+=(Numeric x) {for(size_t i=0; i< R1.size(); i++) R1[i][0]+=x; return *this;}
  
  operator Matrix() const {
    Matrix M;
    switch(stokes_dim) {
      case 4:
        M.resize(R4.size(), 4);
        for(size_t i=0; i<R4.size(); i++)
          for(size_t j=0; j<4; j++)
              M(i, j) = R4[i](j);
        break;
      case 3:
        M.resize(R3.size(), 3);
        for(size_t i=0; i<R3.size(); i++)
          for(size_t j=0; j<3; j++)
              M(i, j) = R3[i](j);
        break;
      case 2:
        M.resize(R2.size(), 2);
        for(size_t i=0; i<R2.size(); i++)
          for(size_t j=0; j<2; j++)
              M(i, j) = R2[i](j);
        break;
      default:
        M.resize(R1.size(), 1);
        for(size_t i=0; i<R1.size(); i++)
          M(i, 0) = R1[i](0);
        break;
    }
    return M;
  }
  
  void setSource(const StokesVector& a, const ConstVectorView& B, const StokesVector& S, Index i)
  {
    switch(stokes_dim) {
      case 4:
        R4[i].noalias() = Eigen::Vector4d(a.Kjj()[i], a.K12()[i], a.K13()[i], a.K14()[i]) * B[i];
        if(not S.IsEmpty())
          R4[i].noalias() += Eigen::Vector4d(S.Kjj()[i], S.K12()[i], S.K13()[i], S.K14()[i]);
        break;
      case 3:
        R3[i].noalias() = Eigen::Vector3d(a.Kjj()[i], a.K12()[i], a.K13()[i]) * B[i];
        if(not S.IsEmpty())
          R3[i].noalias() += Eigen::Vector3d(S.Kjj()[i], S.K12()[i], S.K13()[i]);
        break;
      case 2:
        R2[i].noalias() = Eigen::Vector2d(a.Kjj()[i], a.K12()[i]) * B[i];
        if(not S.IsEmpty())
          R2[i].noalias() += Eigen::Vector2d(S.Kjj()[i], S.K12()[i]);
        break;
      default:
        R1[i][0] = a.Kjj()[i] * B[i];
        if(not S.IsEmpty())
          R1[i][0] += S.Kjj()[i];
    }
  }
};


typedef Array<TransmissionMatrix> ArrayOfTransmissionMatrix;
typedef Array<ArrayOfTransmissionMatrix> ArrayOfArrayOfTransmissionMatrix;
typedef Array<RadiationVector> ArrayOfRadiationVector;
typedef Array<ArrayOfRadiationVector> ArrayOfArrayOfRadiationVector;


inline void update_radiation_vector(RadiationVector& I, const TransmissionMatrix& T, const RadiationVector& J) {
  I -= J;
  I.leftMul(T);
  I += J;
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
                             const ArrayOfTransmissionMatrix& dT2);

void stepwise_source(RadiationVector& J,
                     ArrayOfRadiationVector& dJ,
                     const PropagationMatrix& K,
                     const StokesVector& a,
                     const StokesVector& S,
                     const ArrayOfPropagationMatrix& dK_dx,
                     const ArrayOfStokesVector& da_dx,
                     const ArrayOfStokesVector& dS_dx,
                     const ConstVectorView B,
                     const ConstVectorView dB_dT,
                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                     const bool& jacobian_do);

void stepwise_transmission(TransmissionMatrix& PiT,
                           TransmissionMatrix& T,
                           ArrayOfTransmissionMatrix& dT1,
                           ArrayOfTransmissionMatrix& dT2,
                           const TransmissionMatrix& PiT_last,
                           const PropagationMatrix& K1,
                           const PropagationMatrix& K2,
                           const ArrayOfPropagationMatrix& dK1_dx,
                           const ArrayOfPropagationMatrix& dK2_dx,
                           const Numeric& r,
                           const bool& first);

#endif  // transmissionmatrix_h
