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
  TransmissionMatrix(TransmissionMatrix&& tm)                 = default;
  TransmissionMatrix(const TransmissionMatrix& tm)            = default;
  TransmissionMatrix& operator=(const TransmissionMatrix& tm) = default;
  TransmissionMatrix& operator=(TransmissionMatrix&& tm)      = default;
  
  operator Tensor3() const {
    Tensor3 T(Frequencies(), stokes_dim, stokes_dim);
    for(size_t i=0; i<T4.size(); i++) for(size_t j=0; j<4; j++) for(size_t k=0; k<4; k++) T(i, j, k) = T4[i](j, k);
    for(size_t i=0; i<T3.size(); i++) for(size_t j=0; j<3; j++) for(size_t k=0; k<3; k++) T(i, j, k) = T3[i](j, k);
    for(size_t i=0; i<T2.size(); i++) for(size_t j=0; j<2; j++) for(size_t k=0; k<2; k++) T(i, j, k) = T2[i](j, k);
    for(size_t i=0; i<T1.size(); i++) T(i, 0, 0) = T1[i](0, 0);
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
    for(auto& T: T4) T=Eigen::Matrix4d::Identity();
    for(auto& T: T3) T=Eigen::Matrix3d::Identity();
    for(auto& T: T2) T=Eigen::Matrix2d::Identity();
    for(auto& T: T1) T(0, 0)=1;
  }
  
  void setZero() {
    for(auto& T: T4) T=Eigen::Matrix4d::Zero();
    for(auto& T: T3) T=Eigen::Matrix3d::Zero();
    for(auto& T: T2) T=Eigen::Matrix2d::Zero();
    for(auto& T: T1) T(0, 0)=1;
  }
  
  void mul(const TransmissionMatrix& A, const TransmissionMatrix& B) {
    for(size_t i=0; i<T4.size(); i++) T4[i].noalias() = A.T4[i] * B.T4[i];
    for(size_t i=0; i<T3.size(); i++) T3[i].noalias() = A.T3[i] * B.T3[i];
    for(size_t i=0; i<T2.size(); i++) T2[i].noalias() = A.T2[i] * B.T2[i];
    for(size_t i=0; i<T1.size(); i++) T1[i].noalias() = A.T1[i] * B.T1[i];
  }
  
  void mul_aliased(const TransmissionMatrix& A, const TransmissionMatrix& B) {
    for(size_t i=0; i<T4.size(); i++) T4[i] = A.T4[i] * B.T4[i];
    for(size_t i=0; i<T3.size(); i++) T3[i] = A.T3[i] * B.T3[i];
    for(size_t i=0; i<T2.size(); i++) T2[i] = A.T2[i] * B.T2[i];
    for(size_t i=0; i<T1.size(); i++) T1[i] = A.T1[i] * B.T1[i];
  }
  
  const Numeric& operator()(const Index i, const Index j, const Index k) const {
    switch(stokes_dim) {
      case 4: return T4[i](j, k);
      case 3: return T3[i](j, k);
      case 2: return T2[i](j, k);
      default: return T1[i](j, k);
    }
  }
  
  Index StokesDim() const {return stokes_dim;}
  Index Frequencies() const {
    switch(stokes_dim) {
      case 4: return Index(T4.size());
      case 3: return Index(T3.size());
      case 2: return Index(T2.size());
      default: return Index(T1.size());
    }
  }
  
  friend std::ostream& operator<<(std::ostream& os, const TransmissionMatrix& tm);
  friend std::istream& operator>>(std::istream& data, TransmissionMatrix& tm);
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
  RadiationVector(RadiationVector&& rv)                 = default;
  RadiationVector(const RadiationVector& rv)            = default;
  RadiationVector& operator=(const RadiationVector& rv) = default;
  RadiationVector& operator=(RadiationVector&& rv)      = default;
  
  void leftMul(const TransmissionMatrix& T) {
    for(size_t i=0; i<R4.size(); i++) R4[i] = T.Mat4(i) * R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i] = T.Mat3(i) * R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i] = T.Mat2(i) * R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i] = T.Mat1(i) * R1[i];
  }
  
  void SetZero(size_t i) {
    switch(stokes_dim) {
      case 4: R4[i].noalias()=Eigen::Vector4d::Zero(); break;
      case 3: R3[i].noalias()=Eigen::Vector3d::Zero(); break;
      case 2: R2[i].noalias()=Eigen::Vector2d::Zero(); break;
      case 1: R1[i][0]=0; break;
    }
  }
  
  const Eigen::Vector4d& Vec4(size_t i) const {return R4[i];}
  const Eigen::Vector3d& Vec3(size_t i) const {return R3[i];}
  const Eigen::Vector2d& Vec2(size_t i) const {return R2[i];}
  const Eigen::Matrix<double, 1, 1>& Vec1(size_t i) const {return R1[i];}
  
  Eigen::Vector4d& Vec4(size_t i) {return R4[i];}
  Eigen::Vector3d& Vec3(size_t i) {return R3[i];}
  Eigen::Vector2d& Vec2(size_t i) {return R2[i];}
  Eigen::Matrix<double, 1, 1>& Vec1(size_t i) {return R1[i];}
  
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
  
  void addDerivEmission(const TransmissionMatrix& PiT,
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
  
  void addDerivTransmission(const TransmissionMatrix& PiT,
                            const TransmissionMatrix& dT,
                            const RadiationVector& I)
  {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() += PiT.Mat4(i) * dT.Mat4(i) * I.R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() += PiT.Mat3(i) * dT.Mat3(i) * I.R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() += PiT.Mat2(i) * dT.Mat2(i) * I.R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() += PiT.Mat1(i) * dT.Mat1(i) * I.R1[i];
  }
  
  void setDerivReflection(const RadiationVector& I,
                          const TransmissionMatrix& PiT,
                          const TransmissionMatrix&  Z,
                          const TransmissionMatrix& dZ)
  {
    for(size_t i=0; i<R4.size(); i++) R4[i] = PiT.Mat4(i) * (Z.Mat4(i) * R4[i] + dZ.Mat4(i) * I.R4[i]);
    for(size_t i=0; i<R3.size(); i++) R3[i] = PiT.Mat3(i) * (Z.Mat3(i) * R3[i] + dZ.Mat3(i) * I.R3[i]);
    for(size_t i=0; i<R2.size(); i++) R2[i] = PiT.Mat2(i) * (Z.Mat2(i) * R2[i] + dZ.Mat2(i) * I.R2[i]);
    for(size_t i=0; i<R1.size(); i++) R4[i][0] = PiT.Mat1(i)[0] * (Z.Mat1(i)[0] * R1[i][0] + dZ.Mat1(i)[0] * I.R1[i][0]);
  }
  
  void setBackscatter(const RadiationVector& I0,
                      const TransmissionMatrix& Tr,
                      const TransmissionMatrix& Tf,
                      const TransmissionMatrix& Z)
  {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() = Tr.Mat4(i) * Z.Mat4(i) * Tf.Mat4(i) * I0.R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() = Tr.Mat3(i) * Z.Mat3(i) * Tf.Mat3(i) * I0.R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() = Tr.Mat2(i) * Z.Mat2(i) * Tf.Mat2(i) * I0.R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() = Tr.Mat1(i) * Z.Mat1(i) * Tf.Mat1(i) * I0.R1[i];
  }
  
  void setBackscatterDerivative(const RadiationVector& I0,
                                const TransmissionMatrix& T,
                                const TransmissionMatrix& Tr_past,
                                const TransmissionMatrix& Tf_past,
                                const TransmissionMatrix& Z,
                                const TransmissionMatrix& dT1,
                                const TransmissionMatrix& dT2,
                                const TransmissionMatrix& dZ)
  {
    for(size_t i=0; i<R4.size(); i++) R4[i].noalias() = Tr_past.Mat4(i) * T.Mat4(i) * Z.Mat4(i) * (dT1.Mat4(i) + dT2.Mat4(i)) * Tf_past.Mat4(i) * I0.R4[i] +
                                                        Tr_past.Mat4(i) * T.Mat4(i) * dZ.Mat4(i) * T.Mat4(i) * Tf_past.Mat4(i) * I0.R4[i] + 
                                                        Tr_past.Mat4(i) * (dT1.Mat4(i) + dT2.Mat4(i)) * Z.Mat4(i) * T.Mat4(i) * Tf_past.Mat4(i) * I0.R4[i];
    for(size_t i=0; i<R3.size(); i++) R3[i].noalias() = Tr_past.Mat3(i) * T.Mat3(i) * Z.Mat3(i) * (dT1.Mat3(i) + dT2.Mat3(i)) * Tf_past.Mat3(i) * I0.R3[i] +
                                                        Tr_past.Mat3(i) * T.Mat3(i) * dZ.Mat3(i) * T.Mat3(i) * Tf_past.Mat3(i) * I0.R3[i] + 
                                                        Tr_past.Mat3(i) * (dT1.Mat3(i) + dT2.Mat3(i)) * Z.Mat3(i) * T.Mat3(i) * Tf_past.Mat3(i) * I0.R3[i];
    for(size_t i=0; i<R2.size(); i++) R2[i].noalias() = Tr_past.Mat2(i) * T.Mat2(i) * Z.Mat2(i) * (dT1.Mat2(i) + dT2.Mat2(i)) * Tf_past.Mat2(i) * I0.R2[i] +
                                                        Tr_past.Mat2(i) * T.Mat2(i) * dZ.Mat2(i) * T.Mat2(i) * Tf_past.Mat2(i) * I0.R2[i] + 
                                                        Tr_past.Mat2(i) * (dT1.Mat2(i) + dT2.Mat2(i)) * Z.Mat2(i) * T.Mat2(i) * Tf_past.Mat2(i) * I0.R2[i];
    for(size_t i=0; i<R1.size(); i++) R1[i].noalias() = Tr_past.Mat1(i) * T.Mat1(i) * Z.Mat1(i) * (dT1.Mat1(i) + dT2.Mat1(i)) * Tf_past.Mat1(i) * I0.R1[i] +
                                                        Tr_past.Mat1(i) * T.Mat1(i) * dZ.Mat1(i) * T.Mat1(i) * Tf_past.Mat1(i) * I0.R1[i] + 
                                                        Tr_past.Mat1(i) * (dT1.Mat1(i) + dT2.Mat1(i)) * Z.Mat1(i) * T.Mat1(i) * Tf_past.Mat1(i) * I0.R1[i];
  }
  
  RadiationVector& operator=(const ConstMatrixView& M) {
    assert(M.ncols() == stokes_dim and M.nrows() == Frequencies());
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
  
  operator Matrix() const {
    Matrix M(Frequencies(), stokes_dim);
    for(size_t i=0; i<R4.size(); i++) for(size_t j=0; j<4; j++) M(i, j) = R4[i](j);
    for(size_t i=0; i<R3.size(); i++) for(size_t j=0; j<3; j++) M(i, j) = R3[i](j);
    for(size_t i=0; i<R2.size(); i++) for(size_t j=0; j<2; j++) M(i, j) = R2[i](j);
    for(size_t i=0; i<R1.size(); i++) M(i, 0) = R1[i](0);
    return M;
  }
  
  void setSource(const StokesVector& a, const ConstVectorView& B, const StokesVector& S, Index i) {
    assert(a.NumberOfAzimuthAngles() == 1);
    assert(a.NumberOfZenithAngles() == 1);
    assert(S.NumberOfAzimuthAngles() == 1);
    assert(S.NumberOfZenithAngles() == 1);
    switch(stokes_dim) {
      case 4:
        if(not S.IsEmpty())
          R4[i].noalias() = Eigen::Vector4d(a.Kjj()[i], a.K12()[i], a.K13()[i], a.K14()[i]) * B[i]
                          + Eigen::Vector4d(S.Kjj()[i], S.K12()[i], S.K13()[i], S.K14()[i]);
        else
          R4[i].noalias() = Eigen::Vector4d(a.Kjj()[i], a.K12()[i], a.K13()[i], a.K14()[i]) * B[i];
        break;
      case 3:
        if(not S.IsEmpty())
          R3[i].noalias() = Eigen::Vector3d(a.Kjj()[i], a.K12()[i], a.K13()[i]) * B[i]
                          + Eigen::Vector3d(S.Kjj()[i], S.K12()[i], S.K13()[i]);
        else
          R3[i].noalias() = Eigen::Vector3d(a.Kjj()[i], a.K12()[i], a.K13()[i]) * B[i];
        break;
      case 2:
        if(not S.IsEmpty())
          R2[i].noalias() = Eigen::Vector2d(a.Kjj()[i], a.K12()[i]) * B[i]
                          + Eigen::Vector2d(S.Kjj()[i], S.K12()[i]);
        else
          R2[i].noalias() = Eigen::Vector2d(a.Kjj()[i], a.K12()[i]) * B[i];
        break;
      default:
        if(not S.IsEmpty())
          R1[i][0] = a.Kjj()[i] * B[i] + S.Kjj()[i];
        else
          R1[i][0] = a.Kjj()[i] * B[i];
    }
  }
  
  Index StokesDim() const {return stokes_dim;}
  Index Frequencies() const {
    switch(stokes_dim) {
      case 4: return Index(R4.size());
      case 3: return Index(R3.size());
      case 2: return Index(R2.size());
      default: return Index(R1.size());
    }
  }
  
  friend std::ostream& operator<<(std::ostream& os, const RadiationVector& rv);
  friend std::istream& operator>>(std::istream& data, RadiationVector& rv);
};


typedef Array<TransmissionMatrix> ArrayOfTransmissionMatrix;
typedef Array<ArrayOfTransmissionMatrix> ArrayOfArrayOfTransmissionMatrix;
typedef Array<RadiationVector> ArrayOfRadiationVector;
typedef Array<ArrayOfRadiationVector> ArrayOfArrayOfRadiationVector;


std::ostream& operator<<(std::ostream& os, const TransmissionMatrix& tm);
std::ostream& operator<<(std::ostream& os, const ArrayOfTransmissionMatrix& atm);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfTransmissionMatrix& aatm);
std::ostream& operator<<(std::ostream& os, const RadiationVector& rv);
std::ostream& operator<<(std::ostream& os, const ArrayOfRadiationVector& arv);
std::ostream& operator<<(std::ostream& os, const ArrayOfArrayOfRadiationVector& aarv);

std::istream& operator>>(std::istream& is, TransmissionMatrix& tm);
std::istream& operator>>(std::istream& is, RadiationVector& rv);

enum class BackscatterSolver { Commutative_PureReflectionJacobian, Full, };

enum class CumulativeTransmission { Forward, ForwardReverse, Reflect, };

enum class RadiativeTransferSolver { Emission, Transmission, };

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
                             const ArrayOfTransmissionMatrix& dT2,
                             const RadiativeTransferSolver solver);

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

void stepwise_transmission(TransmissionMatrix& T,
                           ArrayOfTransmissionMatrix& dT1,
                           ArrayOfTransmissionMatrix& dT2,
                           const PropagationMatrix& K1,
                           const PropagationMatrix& K2,
                           const ArrayOfPropagationMatrix& dK1_dx,
                           const ArrayOfPropagationMatrix& dK2_dx,
                           const Numeric& r,
                           const Numeric& dr_dtemp1,
                           const Numeric& dr_dtemp2,
                           const Index temp_deriv_pos);

ArrayOfTransmissionMatrix cumulative_transmission(const ArrayOfTransmissionMatrix& T, const CumulativeTransmission type) /*[[expects: T.nelem()>0]]*/;

void set_backscatter_radiation_vector(ArrayOfRadiationVector& I,
                                      ArrayOfArrayOfRadiationVector& dI,
                                      const ArrayOfTransmissionMatrix& T,
                                      const ArrayOfTransmissionMatrix& PiTf,
                                      const ArrayOfTransmissionMatrix& PiTr,
                                      const ArrayOfTransmissionMatrix& Z,
                                      const ArrayOfArrayOfTransmissionMatrix& dT1,
                                      const ArrayOfArrayOfTransmissionMatrix& dT2,
                                      const ArrayOfArrayOfTransmissionMatrix& dZ,
                                      const BackscatterSolver solver);

ArrayOfTransmissionMatrix cumulative_backscatter(ConstTensor5View t, ConstMatrixView m);

ArrayOfArrayOfTransmissionMatrix cumulative_backscatter_derivative(ConstTensor5View t, const ArrayOfMatrix& aom);

#endif  // transmissionmatrix_h
