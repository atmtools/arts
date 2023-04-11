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
 * @file   transmissionmatrix.c
 * @author Richard Larsson
 * @date   2018-01-30
 * 
 * @brief  Stuff related to the transmission matrix.
 * 
 * Using Eigen library to try and speed up computations.
 */

#include "transmissionmatrix.h"

#include "arts_conversions.h"
#include "debug.h"
#include "double_imanip.h"

TransmissionMatrix::TransmissionMatrix(Index nf, Index stokes)
    : stokes_dim(stokes),
      T4(stokes_dim == 4 ? nf : 0, Eigen::Matrix4d::Identity()),
      T3(stokes_dim == 3 ? nf : 0, Eigen::Matrix3d::Identity()),
      T2(stokes_dim == 2 ? nf : 0, Eigen::Matrix2d::Identity()),
      T1(stokes_dim == 1 ? nf : 0, Eigen::Matrix<double, 1, 1>::Identity()) {
  ARTS_ASSERT(stokes_dim < 5 and stokes_dim > 0);
}

const Eigen::Matrix4d& TransmissionMatrix::Mat4(size_t i) const { return T4[i]; }
const Eigen::Matrix3d& TransmissionMatrix::Mat3(size_t i) const { return T3[i]; }
const Eigen::Matrix2d& TransmissionMatrix::Mat2(size_t i) const { return T2[i]; }
const Eigen::Matrix<double, 1, 1>& TransmissionMatrix::Mat1(size_t i) const { return T1[i]; }

Eigen::Matrix4d& TransmissionMatrix::Mat4(size_t i) { return T4[i]; }
Eigen::Matrix3d& TransmissionMatrix::Mat3(size_t i) { return T3[i]; }
Eigen::Matrix2d& TransmissionMatrix::Mat2(size_t i) { return T2[i]; }
Eigen::Matrix<double, 1, 1>& TransmissionMatrix::Mat1(size_t i) { return T1[i]; }

TransmissionMatrix& TransmissionMatrix::operator=(
    const LazyScale<TransmissionMatrix>& lstm) {
  operator=(lstm.bas);
  operator*=(lstm.scale);
  return *this;
}

TransmissionMatrix::operator Tensor3() const {
  Tensor3 T(Frequencies(), stokes_dim, stokes_dim);
  for (size_t i = 0; i < T4.size(); i++)
    for (size_t j = 0; j < 4; j++)
      for (size_t k = 0; k < 4; k++) T(i, j, k) = T4[i](j, k);
  for (size_t i = 0; i < T3.size(); i++)
    for (size_t j = 0; j < 3; j++)
      for (size_t k = 0; k < 3; k++) T(i, j, k) = T3[i](j, k);
  for (size_t i = 0; i < T2.size(); i++)
    for (size_t j = 0; j < 2; j++)
      for (size_t k = 0; k < 2; k++) T(i, j, k) = T2[i](j, k);
  for (size_t i = 0; i < T1.size(); i++) T(i, 0, 0) = T1[i](0, 0);
  return T;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
[[nodiscard]] Eigen::MatrixXd TransmissionMatrix::Mat(size_t i) const {
  switch (stokes_dim) {
    case 1:
      return Mat1(i);
    case 2:
      return Mat2(i);
    case 3:
      return Mat3(i);
    case 4:
      return Mat4(i);
  }
}
#pragma GCC diagnostic pop

void TransmissionMatrix::setIdentity() {
  std::fill(T4.begin(), T4.end(), Eigen::Matrix4d::Identity());
  std::fill(T3.begin(), T3.end(), Eigen::Matrix3d::Identity());
  std::fill(T2.begin(), T2.end(), Eigen::Matrix2d::Identity());
  std::fill(T1.begin(), T1.end(), Eigen::Matrix<double, 1, 1>::Identity());
}

void TransmissionMatrix::setZero() {
  std::fill(T4.begin(), T4.end(), Eigen::Matrix4d::Zero());
  std::fill(T3.begin(), T3.end(), Eigen::Matrix3d::Zero());
  std::fill(T2.begin(), T2.end(), Eigen::Matrix2d::Zero());
  std::fill(T1.begin(), T1.end(), Eigen::Matrix<double, 1, 1>::Zero());
}

void TransmissionMatrix::mul(const TransmissionMatrix& A,
                             const TransmissionMatrix& B) {
  for (size_t i = 0; i < T4.size(); i++) T4[i].noalias() = A.T4[i] * B.T4[i];
  for (size_t i = 0; i < T3.size(); i++) T3[i].noalias() = A.T3[i] * B.T3[i];
  for (size_t i = 0; i < T2.size(); i++) T2[i].noalias() = A.T2[i] * B.T2[i];
  for (size_t i = 0; i < T1.size(); i++) T1[i].noalias() = A.T1[i] * B.T1[i];
}

void TransmissionMatrix::mul_aliased(const TransmissionMatrix& A,
                                     const TransmissionMatrix& B) {
  for (size_t i = 0; i < T4.size(); i++) T4[i] = A.T4[i] * B.T4[i];
  for (size_t i = 0; i < T3.size(); i++) T3[i] = A.T3[i] * B.T3[i];
  for (size_t i = 0; i < T2.size(); i++) T2[i] = A.T2[i] * B.T2[i];
  for (size_t i = 0; i < T1.size(); i++) T1[i] = A.T1[i] * B.T1[i];
}

Numeric TransmissionMatrix::operator()(const Index i,
                                       const Index j,
                                       const Index k) const {
  switch (stokes_dim) {
    case 4:
      return T4[i](j, k);
    case 3:
      return T3[i](j, k);
    case 2:
      return T2[i](j, k);
    default:
      return T1[i](j, k);
  }
}

Numeric& TransmissionMatrix::operator()(const Index i, const Index j, const Index k) {
  switch (stokes_dim) {
    case 4:
      return T4[i](j, k);
    case 3:
      return T3[i](j, k);
    case 2:
      return T2[i](j, k);
    default:
      return T1[i](j, k);
  }
}

[[nodiscard]] Index TransmissionMatrix::Frequencies() const {
  switch (stokes_dim) {
    case 4:
      return Index(T4.size());
    case 3:
      return Index(T3.size());
    case 2:
      return Index(T2.size());
    default:
      return Index(T1.size());
  }
}

TransmissionMatrix::TransmissionMatrix(const ConstMatrixView& mat): TransmissionMatrix(1, mat.nrows()){
  ARTS_ASSERT(mat.nrows() == mat.ncols());
  for (Index i = 0; i < stokes_dim; i++)
    for (Index j = 0; j < stokes_dim; j++)
      operator()(0, i, j) = mat(i,j);
};

TransmissionMatrix& TransmissionMatrix::operator+=(
    const LazyScale<TransmissionMatrix>& lstm) {
  for (size_t i = 0; i < T4.size(); i++)
    T4[i].noalias() = lstm.scale * lstm.bas.Mat4(i);
  for (size_t i = 0; i < T3.size(); i++)
    T3[i].noalias() = lstm.scale * lstm.bas.Mat3(i);
  for (size_t i = 0; i < T2.size(); i++)
    T2[i].noalias() = lstm.scale * lstm.bas.Mat2(i);
  for (size_t i = 0; i < T1.size(); i++)
    T1[i].noalias() = lstm.scale * lstm.bas.Mat1(i);
  return *this;
}

TransmissionMatrix& TransmissionMatrix::operator*=(const Numeric& scale) {
  std::transform(
      T4.begin(), T4.end(), T4.begin(), [scale](auto& T) { return T * scale; });
  std::transform(
      T3.begin(), T3.end(), T3.begin(), [scale](auto& T) { return T * scale; });
  std::transform(
      T2.begin(), T2.end(), T2.begin(), [scale](auto& T) { return T * scale; });
  std::transform(
      T1.begin(), T1.end(), T1.begin(), [scale](auto& T) { return T * scale; });
  return *this;
}

LazyScale<TransmissionMatrix> operator*(const TransmissionMatrix& tm,
                                        const Numeric& x) {
  return {tm, x};
}

LazyScale<TransmissionMatrix> operator*(const Numeric& x,
                                        const TransmissionMatrix& tm) {
  return {tm, x};
}

RadiationVector::RadiationVector(Index nf, Index stokes)
    : stokes_dim(stokes),
      R4(stokes_dim == 4 ? nf : 0, Eigen::Vector4d::Zero()),
      R3(stokes_dim == 3 ? nf : 0, Eigen::Vector3d::Zero()),
      R2(stokes_dim == 2 ? nf : 0, Eigen::Vector2d::Zero()),
      R1(stokes_dim == 1 ? nf : 0, Eigen::Matrix<double, 1, 1>::Zero()) {
  ARTS_ASSERT(stokes_dim < 5 and stokes_dim > 0);
}

const Eigen::Vector4d& RadiationVector::Vec4(size_t i) const { return R4[i]; }
const Eigen::Vector3d& RadiationVector::Vec3(size_t i) const { return R3[i]; }
const Eigen::Vector2d& RadiationVector::Vec2(size_t i) const { return R2[i]; }
const Eigen::Matrix<double, 1, 1>& RadiationVector::Vec1(size_t i) const {
  return R1[i];
}

Eigen::Vector4d& RadiationVector::Vec4(size_t i) { return R4[i]; }
Eigen::Vector3d& RadiationVector::Vec3(size_t i) { return R3[i]; }
Eigen::Vector2d& RadiationVector::Vec2(size_t i) { return R2[i]; }
Eigen::Matrix<double, 1, 1>& RadiationVector::Vec1(size_t i) { return R1[i]; }

Numeric& RadiationVector::operator()(const Index i, const Index j) {
  switch (stokes_dim) {
    case 4:
      return R4[i][j];
    case 3:
      return R3[i][j];
    case 2:
      return R2[i][j];
    default:
      return R1[i][j];
  }
}

void RadiationVector::leftMul(const TransmissionMatrix& T) {
  for (size_t i = 0; i < R4.size(); i++) R4[i] = T.Mat4(i) * R4[i];
  for (size_t i = 0; i < R3.size(); i++) R3[i] = T.Mat3(i) * R3[i];
  for (size_t i = 0; i < R2.size(); i++) R2[i] = T.Mat2(i) * R2[i];
  for (size_t i = 0; i < R1.size(); i++) R1[i] = T.Mat1(i) * R1[i];
}

void RadiationVector::SetZero(size_t i) {
  switch (stokes_dim) {
    case 4:
      R4[i].noalias() = Eigen::Vector4d::Zero();
      break;
    case 3:
      R3[i].noalias() = Eigen::Vector3d::Zero();
      break;
    case 2:
      R2[i].noalias() = Eigen::Vector2d::Zero();
      break;
    case 1:
      R1[i][0] = 0;
      break;
  }
}

void RadiationVector::SetZero() {
  for (Index i=0; i<Frequencies(); i++) SetZero(i);
}

void RadiationVector::SetUnity() {
  SetZero();
  for (auto& v: R4) v[0] = 1;
  for (auto& v: R3) v[0] = 1;
  for (auto& v: R2) v[0] = 1;
  for (auto& v: R1) v[0] = 1;
}

Eigen::VectorXd RadiationVector::Vec(size_t i) const {
  switch (stokes_dim) {
    case 4:
      return Vec4(i);
    case 3:
      return Vec3(i);
    case 2:
      return Vec2(i);
    default:
      return Vec1(i);
  }
}

void RadiationVector::rem_avg(const RadiationVector& O1,
                              const RadiationVector& O2) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() -= 0.5 * (O1.R4[i] + O2.R4[i]);
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() -= 0.5 * (O1.R3[i] + O2.R3[i]);
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() -= 0.5 * (O1.R2[i] + O2.R2[i]);
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() -= 0.5 * (O1.R1[i] + O2.R1[i]);
}

void RadiationVector::add_avg(const RadiationVector& O1,
                              const RadiationVector& O2) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() += 0.5 * (O1.R4[i] + O2.R4[i]);
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() += 0.5 * (O1.R3[i] + O2.R3[i]);
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() += 0.5 * (O1.R2[i] + O2.R2[i]);
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() += 0.5 * (O1.R1[i] + O2.R1[i]);
}

void RadiationVector::add_weighted(const TransmissionMatrix& T,
                                   const RadiationVector& far,
                                   const RadiationVector& close,
                                   const ConstMatrixView& Kfar,
                                   const ConstMatrixView& Kclose,
                                   const Numeric r) {
  for (size_t i = 0; i < R4.size(); i++) {
    R4[i].noalias() +=
        T.second_order_integration_source<4>(T.TraMat<4>(i),
                                             far.R4[i],
                                             close.R4[i],
                                             prop_matrix<4>(Kfar(i, joker)),
                                             prop_matrix<4>(Kclose(i, joker)),
                                             r);
  }
  for (size_t i = 0; i < R3.size(); i++) {
    R3[i].noalias() +=
        T.second_order_integration_source<3>(T.TraMat<3>(i),
                                             far.R3[i],
                                             close.R3[i],
                                             prop_matrix<3>(Kfar(i, joker)),
                                             prop_matrix<3>(Kclose(i, joker)),
                                             r);
  }
  for (size_t i = 0; i < R2.size(); i++) {
    R2[i].noalias() +=
        T.second_order_integration_source<2>(T.TraMat<2>(i),
                                             far.R2[i],
                                             close.R2[i],
                                             prop_matrix<2>(Kfar(i, joker)),
                                             prop_matrix<2>(Kclose(i, joker)),
                                             r);
  }
  for (size_t i = 0; i < R1.size(); i++) {
    R1[i].noalias() +=
        T.second_order_integration_source<1>(T.TraMat<1>(i),
                                             far.R1[i],
                                             close.R1[i],
                                             prop_matrix<1>(Kfar(i, joker)),
                                             prop_matrix<1>(Kclose(i, joker)),
                                             r);
  }
}

void RadiationVector::addDerivEmission(const TransmissionMatrix& PiT,
                                       const TransmissionMatrix& dT,
                                       const TransmissionMatrix& T,
                                       const RadiationVector& ImJ,
                                       const RadiationVector& dJ) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() += PiT.Mat4(i) * (dT.Mat4(i) * ImJ.R4[i] + dJ.R4[i] -
                                      T.Mat4(i) * dJ.R4[i]);
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() += PiT.Mat3(i) * (dT.Mat3(i) * ImJ.R3[i] + dJ.R3[i] -
                                      T.Mat3(i) * dJ.R3[i]);
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() += PiT.Mat2(i) * (dT.Mat2(i) * ImJ.R2[i] + dJ.R2[i] -
                                      T.Mat2(i) * dJ.R2[i]);
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() += PiT.Mat1(i) * (dT.Mat1(i) * ImJ.R1[i] + dJ.R1[i] -
                                      T.Mat1(i) * dJ.R1[i]);
}

void RadiationVector::addWeightedDerivEmission(const TransmissionMatrix& PiT,
                                               const TransmissionMatrix& dT,
                                               const TransmissionMatrix& T,
                                               const RadiationVector& I,
                                               const RadiationVector& far,
                                               const RadiationVector& close,
                                               const RadiationVector& d,
                                               bool isfar) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() +=
        PiT.Mat4(i) * (dT.Mat4(i) * I.R4[i] +
                       T.second_order_integration_dsource<4>(
                           i, dT, far.R4[i], close.R4[i], d.R4[i], isfar));
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() +=
        PiT.Mat3(i) * (dT.Mat3(i) * I.R3[i] +
                       T.second_order_integration_dsource<3>(
                           i, dT, far.R3[i], close.R3[i], d.R3[i], isfar));
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() +=
        PiT.Mat2(i) * (dT.Mat2(i) * I.R2[i] +
                       T.second_order_integration_dsource<2>(
                           i, dT, far.R2[i], close.R2[i], d.R2[i], isfar));
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() +=
        PiT.Mat1(i) * (dT.Mat1(i) * I.R1[i] +
                       T.second_order_integration_dsource<1>(
                           i, dT, far.R1[i], close.R1[i], d.R1[i], isfar));
}

void RadiationVector::addDerivTransmission(const TransmissionMatrix& PiT,
                                           const TransmissionMatrix& dT,
                                           const RadiationVector& I) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() += PiT.Mat4(i) * dT.Mat4(i) * I.R4[i];
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() += PiT.Mat3(i) * dT.Mat3(i) * I.R3[i];
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() += PiT.Mat2(i) * dT.Mat2(i) * I.R2[i];
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() += PiT.Mat1(i) * dT.Mat1(i) * I.R1[i];
}

void RadiationVector::addMultiplied(const TransmissionMatrix& A,
                                    const RadiationVector& x) {
  for (size_t i = 0; i < R4.size(); i++) R4[i].noalias() += A.Mat4(i) * x.R4[i];
  for (size_t i = 0; i < R3.size(); i++) R3[i].noalias() += A.Mat3(i) * x.R3[i];
  for (size_t i = 0; i < R2.size(); i++) R2[i].noalias() += A.Mat2(i) * x.R2[i];
  for (size_t i = 0; i < R1.size(); i++) R1[i].noalias() += A.Mat1(i) * x.R1[i];
}

void RadiationVector::setDerivReflection(const RadiationVector& I,
                                         const TransmissionMatrix& PiT,
                                         const TransmissionMatrix& Z,
                                         const TransmissionMatrix& dZ) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i] = PiT.Mat4(i) * (Z.Mat4(i) * R4[i] + dZ.Mat4(i) * I.R4[i]);
  for (size_t i = 0; i < R3.size(); i++)
    R3[i] = PiT.Mat3(i) * (Z.Mat3(i) * R3[i] + dZ.Mat3(i) * I.R3[i]);
  for (size_t i = 0; i < R2.size(); i++)
    R2[i] = PiT.Mat2(i) * (Z.Mat2(i) * R2[i] + dZ.Mat2(i) * I.R2[i]);
  for (size_t i = 0; i < R1.size(); i++)
    R4[i][0] =
        PiT.Mat1(i)[0] * (Z.Mat1(i)[0] * R1[i][0] + dZ.Mat1(i)[0] * I.R1[i][0]);
}

void RadiationVector::setBackscatterTransmission(const RadiationVector& I0,
                                                 const TransmissionMatrix& Tr,
                                                 const TransmissionMatrix& Tf,
                                                 const TransmissionMatrix& Z) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() = Tr.Mat4(i) * Z.Mat4(i) * Tf.Mat4(i) * I0.R4[i];
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() = Tr.Mat3(i) * Z.Mat3(i) * Tf.Mat3(i) * I0.R3[i];
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() = Tr.Mat2(i) * Z.Mat2(i) * Tf.Mat2(i) * I0.R2[i];
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() = Tr.Mat1(i) * Z.Mat1(i) * Tf.Mat1(i) * I0.R1[i];
}

void RadiationVector::setBackscatterTransmissionDerivative(
    const RadiationVector& I0,
    const TransmissionMatrix& Tr,
    const TransmissionMatrix& Tf,
    const TransmissionMatrix& dZ) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() += Tr.Mat4(i) * dZ.Mat4(i) * Tf.Mat4(i) * I0.R4[i];
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() += Tr.Mat3(i) * dZ.Mat3(i) * Tf.Mat3(i) * I0.R3[i];
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() += Tr.Mat2(i) * dZ.Mat2(i) * Tf.Mat2(i) * I0.R2[i];
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() += Tr.Mat1(i) * dZ.Mat1(i) * Tf.Mat1(i) * I0.R1[i];
}

RadiationVector::RadiationVector(const ConstMatrixView& M) : RadiationVector(M.nrows(), M.ncols()) {
  for (size_t i = 0; i < R4.size(); i++) {
    R4[i][0] = M(i, 0);
    R4[i][1] = M(i, 1);
    R4[i][2] = M(i, 2);
    R4[i][3] = M(i, 3);
  }
  for (size_t i = 0; i < R3.size(); i++) {
    R3[i][0] = M(i, 0);
    R3[i][1] = M(i, 1);
    R3[i][2] = M(i, 2);
  }
  for (size_t i = 0; i < R2.size(); i++) {
    R2[i][0] = M(i, 0);
    R2[i][1] = M(i, 1);
  }
  for (size_t i = 0; i < R1.size(); i++) {
    R1[i][0] = M(i, 0);
  }
}

RadiationVector& RadiationVector::operator+=(const RadiationVector& rv) {
  for (size_t i = 0; i < R4.size(); i++)
    R4[i].noalias() += rv.R4[i];
  for (size_t i = 0; i < R3.size(); i++)
    R3[i].noalias() += rv.R3[i];
  for (size_t i = 0; i < R2.size(); i++)
    R2[i].noalias() += rv.R2[i];
  for (size_t i = 0; i < R1.size(); i++)
    R1[i].noalias() += rv.R1[i];

  return *this;
}


const Numeric& RadiationVector::operator()(const Index i, const Index j) const {
  switch (stokes_dim) {
    case 4:
      return R4[i][j];
    case 3:
      return R3[i][j];
    case 2:
      return R2[i][j];
    default:
      return R1[i][j];
  }
}

RadiationVector::operator Matrix() const {
  Matrix M(Frequencies(), stokes_dim);
  for (size_t i = 0; i < R4.size(); i++)
    for (size_t j = 0; j < 4; j++) M(i, j) = R4[i](j);
  for (size_t i = 0; i < R3.size(); i++)
    for (size_t j = 0; j < 3; j++) M(i, j) = R3[i](j);
  for (size_t i = 0; i < R2.size(); i++)
    for (size_t j = 0; j < 2; j++) M(i, j) = R2[i](j);
  for (size_t i = 0; i < R1.size(); i++) M(i, 0) = R1[i](0);
  return M;
}

void RadiationVector::setSource(const StokesVector& a,
                                const ConstVectorView& B,
                                const StokesVector& S,
                                Index i) {
  ARTS_ASSERT(a.NumberOfAzimuthAngles() == 1);
  ARTS_ASSERT(a.NumberOfZenithAngles() == 1);
  ARTS_ASSERT(S.NumberOfAzimuthAngles() == 1);
  ARTS_ASSERT(S.NumberOfZenithAngles() == 1);
  switch (stokes_dim) {
    case 4:
      if (not S.IsEmpty())
        R4[i].noalias() =
            Eigen::Vector4d(a.Kjj()[i], a.K12()[i], a.K13()[i], a.K14()[i]) *
                B[i] +
            Eigen::Vector4d(S.Kjj()[i], S.K12()[i], S.K13()[i], S.K14()[i]);
      else
        R4[i].noalias() =
            Eigen::Vector4d(a.Kjj()[i], a.K12()[i], a.K13()[i], a.K14()[i]) *
            B[i];
      break;
    case 3:
      if (not S.IsEmpty())
        R3[i].noalias() =
            Eigen::Vector3d(a.Kjj()[i], a.K12()[i], a.K13()[i]) * B[i] +
            Eigen::Vector3d(S.Kjj()[i], S.K12()[i], S.K13()[i]);
      else
        R3[i].noalias() =
            Eigen::Vector3d(a.Kjj()[i], a.K12()[i], a.K13()[i]) * B[i];
      break;
    case 2:
      if (not S.IsEmpty())
        R2[i].noalias() = Eigen::Vector2d(a.Kjj()[i], a.K12()[i]) * B[i] +
                          Eigen::Vector2d(S.Kjj()[i], S.K12()[i]);
      else
        R2[i].noalias() = Eigen::Vector2d(a.Kjj()[i], a.K12()[i]) * B[i];
      break;
    default:
      if (not S.IsEmpty())
        R1[i][0] = a.Kjj()[i] * B[i] + S.Kjj()[i];
      else
        R1[i][0] = a.Kjj()[i] * B[i];
  }
}

Index RadiationVector::Frequencies() const {
  switch (stokes_dim) {
    case 4:
      return Index(R4.size());
    case 3:
      return Index(R3.size());
    case 2:
      return Index(R2.size());
    default:
      return Index(R1.size());
  }
}

constexpr Numeric lower_is_considered_zero_for_sinc_likes = 1e-4;

template <int N>
constexpr Eigen::Matrix<Numeric, N, 1> source_vector(
    const StokesVector& a,
    const ConstVectorView& B,
    const StokesVector& da,
    const ConstVectorView& dB_dT,
    const StokesVector& dS,
    bool dT,
    Index i) {
  static_assert(N > 0 and N < 5, "Bad stokes dimensions");

  if constexpr (N == 1) {
    if (dT)
      return Eigen::Matrix<Numeric, 1, 1>(dS.Kjj()[i] + da.Kjj()[i] * B[i] +
                                          a.Kjj()[i] * dB_dT[i]);
    return Eigen::Matrix<Numeric, 1, 1>(dS.Kjj()[i] + da.Kjj()[i] * B[i]);
  }

  if constexpr (N == 2) {
    if (dT)
      return Eigen::Vector2d(
          dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
          dS.K12()[i] + da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i]);
    return Eigen::Vector2d(dS.Kjj()[i] + da.Kjj()[i] * B[i],
                           dS.K12()[i] + da.K12()[i] * B[i]);
  }

  if constexpr (N == 3) {
    if (dT)
      return Eigen::Vector3d(
          dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
          dS.K12()[i] + da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i],
          dS.K13()[i] + da.K13()[i] * B[i] + a.K13()[i] * dB_dT[i]);
    return Eigen::Vector3d(dS.Kjj()[i] + da.Kjj()[i] * B[i],
                           dS.K12()[i] + da.K12()[i] * B[i],
                           dS.K13()[i] + da.K13()[i] * B[i]);
  }

  if constexpr (N == 4) {
    if (dT)
      return Eigen::Vector4d(
          dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
          dS.K12()[i] + da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i],
          dS.K13()[i] + da.K13()[i] * B[i] + a.K13()[i] * dB_dT[i],
          dS.K14()[i] + da.K14()[i] * B[i] + a.K14()[i] * dB_dT[i]);
    return Eigen::Vector4d(dS.Kjj()[i] + da.Kjj()[i] * B[i],
                           dS.K12()[i] + da.K12()[i] * B[i],
                           dS.K13()[i] + da.K13()[i] * B[i],
                           dS.K14()[i] + da.K14()[i] * B[i]);
  }
}

template <int N>
constexpr Eigen::Matrix<Numeric, N, 1> source_vector(
    const StokesVector& a,
    const ConstVectorView& B,
    const StokesVector& da,
    const ConstVectorView& dB_dT,
    bool dT,
    Index i) {
  static_assert(N > 0 and N < 5, "Bad stokes dimensions");

  if constexpr (N == 1) {
    if (dT)
      return Eigen::Matrix<Numeric, 1, 1>(da.Kjj()[i] * B[i] +
                                          a.Kjj()[i] * dB_dT[i]);
    return Eigen::Matrix<Numeric, 1, 1>(da.Kjj()[i] * B[i]);
  }

  if constexpr (N == 2) {
    if (dT)
      return Eigen::Vector2d(da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
                             da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i]);
    return Eigen::Vector2d(da.Kjj()[i] * B[i], da.K12()[i] * B[i]);
  }

  if constexpr (N == 3) {
    if (dT)
      return Eigen::Vector3d(da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
                             da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i],
                             da.K13()[i] * B[i] + a.K13()[i] * dB_dT[i]);
    return Eigen::Vector3d(
        da.Kjj()[i] * B[i], da.K12()[i] * B[i], da.K13()[i] * B[i]);
  }

  if constexpr (N == 4) {
    if (dT)
      return Eigen::Vector4d(da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
                             da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i],
                             da.K13()[i] * B[i] + a.K13()[i] * dB_dT[i],
                             da.K14()[i] * B[i] + a.K14()[i] * dB_dT[i]);
    return Eigen::Vector4d(da.Kjj()[i] * B[i],
                           da.K12()[i] * B[i],
                           da.K13()[i] * B[i],
                           da.K14()[i] * B[i]);
  }
}

inline void transmat1(TransmissionMatrix& T,
                      const PropagationMatrix& K1,
                      const PropagationMatrix& K2,
                      const Numeric& r,
                      const Index iz = 0,
                      const Index ia = 0) noexcept {
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++)
    T.Mat1(i)(0, 0) =
        std::exp(-0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]));
}

inline void transmat2(TransmissionMatrix& T,
                      const PropagationMatrix& K1,
                      const PropagationMatrix& K2,
                      const Numeric& r,
                      const Index iz = 0,
                      const Index ia = 0) noexcept {
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]),
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    const Numeric cb = std::cosh(b), sb = std::sinh(b);
    T.Mat2(i).noalias() =
        (Eigen::Matrix2d() << cb, sb, sb, cb).finished() * exp_a;
  }
}

inline void transmat3(TransmissionMatrix& T,
                      const PropagationMatrix& K1,
                      const PropagationMatrix& K2,
                      const Numeric& r,
                      const Index iz = 0,
                      const Index ia = 0) noexcept {
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]),
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]),
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]),
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);

    if (b == 0. and c == 0. and u == 0.) {
      T.Mat3(i).noalias() = Eigen::Matrix3d::Identity() * exp_a;
    } else {
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      const Numeric Const = b2 + c2 - u2;

      const bool real = Const > 0;
      const bool imag = Const < 0;
      const bool either = real or imag;

      const Numeric x =
          std::sqrt(imag ? -Const : Const);  // test to just use real values
      const Numeric x2 =
          (real ? 1 : -1) * x * x;  // test to change sign if imaginary
      const Numeric inv_x2 =
          either ? 1.0 / x2
                 : 1.0;  // test so further calculations are replaced as x→0

      const Numeric sx =
          real ? std::sinh(x) : std::sin(-x);  // -i sin(ix) → sinh(x)
      const Numeric cx =
          real ? std::cosh(x) : std::cos(+x);  //    cos(ix) → cosh(x)

      /* Using:
       *    lim x→0 [(cosh(x) - 1) / x^2] → 1/2
       *    lim x→0 [sinh(x) / x]  → 1
       *    inv_x2 := 1 for x == 0,
       *    C0, C1, C2 ∝ [1/x^2]
       */
      const Numeric C0 =
          either ? a2 * (cx - 1.0) - a * x * sx + x2 : 1.0 + 0.5 * a2 - a;
      const Numeric C1 = either ? 2.0 * a * (1.0 - cx) + x * sx : 1.0 - a;
      const Numeric C2 = either ? cx - 1.0 : 0.5;

      T.Mat3(i).noalias() =
          exp_a * inv_x2 *
          (Eigen::Matrix3d() << C0 + C1 * a + C2 * (a2 + b2 + c2),
           C1 * b + C2 * (2 * a * b - c * u),
           C1 * c + C2 * (2 * a * c + b * u),
           C1 * b + C2 * (2 * a * b + c * u),
           C0 + C1 * a + C2 * (a2 + b2 - u2),
           C1 * u + C2 * (2 * a * u + b * c),
           C1 * c + C2 * (2 * a * c - b * u),
           -C1 * u - C2 * (2 * a * u - b * c),
           C0 + C1 * a + C2 * (a2 + c2 - u2))
              .finished();
    }
  }
}

inline void transmat4(TransmissionMatrix& T,
                      const PropagationMatrix& K1,
                      const PropagationMatrix& K2,
                      const Numeric& r,
                      const Index iz = 0,
                      const Index ia = 0) noexcept {
  static constexpr Numeric sqrt_05 = Constant::inv_sqrt_2;
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]),
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]),
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]),
                  d = -0.5 * r * (K1.K14(iz, ia)[i] + K2.K14(iz, ia)[i]),
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]),
                  v = -0.5 * r * (K1.K24(iz, ia)[i] + K2.K24(iz, ia)[i]),
                  w = -0.5 * r * (K1.K34(iz, ia)[i] + K2.K34(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);

    if (b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.)
      T.Mat4(i).noalias() = Eigen::Matrix4d::Identity() * exp_a;
    else {
      const Numeric b2 = b * b, c2 = c * c, d2 = d * d, u2 = u * u, v2 = v * v,
                    w2 = w * w;

      const Numeric tmp =
          w2 * w2 + 2 * (b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                         c2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                         d2 * (d2 * 0.5 + u2 - v2 - w2) +
                         u2 * (u2 * 0.5 + v2 + w2) + v2 * (v2 * 0.5 + w2) +
                         4 * (b * d * u * w - b * c * v * w - c * d * u * v));
      const Complex Const1 = std::sqrt(Complex(tmp, 0));
      const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;

      const Complex x = std::sqrt(Const2 + Const1) * sqrt_05;
      const Complex y = std::sqrt(Const2 - Const1) * sqrt_05 * Complex(0, 1);
      const Complex x2 = x * x;
      const Complex y2 = y * y;
      const Complex cy = std::cos(y);
      const Complex sy = std::sin(y);
      const Complex cx = std::cosh(x);
      const Complex sx = std::sinh(x);

      const bool x_zero = std::abs(x) < lower_is_considered_zero_for_sinc_likes;
      const bool y_zero = std::abs(y) < lower_is_considered_zero_for_sinc_likes;
      const bool both_zero = y_zero and x_zero;
      const bool either_zero = y_zero or x_zero;

      /* Using:
         *    lim x→0 [({cosh(x),cos(x)} - 1) / x^2] → 1/2
         *    lim x→0 [{sinh(x),sin(x)} / x]  → 1
         *    inv_x2 := 1 for x == 0,
         *    -i sin(ix) → sinh(x)
         *    cos(ix) → cosh(x)
         *    C0, C1, C2 ∝ [1/x^2]
         */
      const Complex ix = x_zero ? 0.0 : 1.0 / x;
      const Complex iy = y_zero ? 0.0 : 1.0 / y;
      const Complex inv_x2y2 =
          both_zero
              ? 1.0
              : 1.0 /
                    (x2 + y2);  // The first "1.0" is the trick for above limits

      const Numeric C0 =
          either_zero ? 1.0 : ((cy * x2 + cx * y2) * inv_x2y2).real();
      const Numeric C1 =
          either_zero ? 1.0 : ((sy * x2 * iy + sx * y2 * ix) * inv_x2y2).real();
      const Numeric C2 = both_zero ? 0.5 : ((cx - cy) * inv_x2y2).real();
      const Numeric C3 = both_zero ? 1.0 / 6.0
                                   : ((x_zero   ? 1.0 - sy * iy
                                       : y_zero ? sx * ix - 1.0
                                                : sx * ix - sy * iy) *
                                      inv_x2y2)
                                         .real();
      T.Mat4(i).noalias() =
          exp_a * (Eigen::Matrix4d() << C0 + C2 * (b2 + c2 + d2),
                   C1 * b + C2 * (-c * u - d * v) +
                       C3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                             v * (b * v + c * w)),
                   C1 * c + C2 * (b * u - d * w) +
                       C3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                             w * (b * v + c * w)),
                   C1 * d + C2 * (b * v + c * w) +
                       C3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                             w * (b * u - d * w)),

                   C1 * b + C2 * (c * u + d * v) +
                       C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                             d * (b * d + u * w)),
                   C0 + C2 * (b2 - u2 - v2),
                   C2 * (b * c - v * w) + C1 * u +
                       C3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                             w * (b * d + u * w)),
                   C2 * (b * d + u * w) + C1 * v +
                       C3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                             w * (b * c - v * w)),

                   C1 * c + C2 * (-b * u + d * w) +
                       C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                             d * (c * d - u * v)),
                   C2 * (b * c - v * w) - C1 * u +
                       C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                             v * (c * d - u * v)),
                   C0 + C2 * (c2 - u2 - w2),
                   C2 * (c * d - u * v) + C1 * w +
                       C3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                             w * (-c2 + u2 + w2)),

                   C1 * d + C2 * (-b * v - c * w) +
                       C3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                             d * (-d2 + v2 + w2)),
                   C2 * (b * d + u * w) - C1 * v +
                       C3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                             v * (-d2 + v2 + w2)),
                   C2 * (c * d - u * v) - C1 * w +
                       C3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                             w * (-d2 + v2 + w2)),
                   C0 + C2 * (d2 - v2 - w2))
                      .finished();
    }
  }
}

inline void dtransmat1(TransmissionMatrix& T,
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
                       const Index ia) noexcept {
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    T.Mat1(i)(0, 0) =
        std::exp(-0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]));
    for (Index j = 0; j < dT1.nelem(); j++) {
      if (dK1[j].NumberOfFrequencies())
        dT1[j].Mat1(i)(0, 0) =
            T.Mat1(i)(0, 0) *
            (-0.5 *
             (r * dK1[j].Kjj(iz, ia)[i] +
              ((j == it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])
                         : 0.0)));
      if (dK2[j].NumberOfFrequencies())
        dT2[j].Mat1(i)(0, 0) =
            T.Mat1(i)(0, 0) *
            (-0.5 *
             (r * dK2[j].Kjj(iz, ia)[i] +
              ((j == it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])
                         : 0.0)));
    }
  }
}

inline void dtransmat2(TransmissionMatrix& T,
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
                       const Index ia) noexcept {
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]),
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);
    const Numeric cb = std::cosh(b), sb = std::sinh(b);
    T.Mat2(i).noalias() =
        (Eigen::Matrix2d() << cb, sb, sb, cb).finished() * exp_a;
    for (Index j = 0; j < dT1.nelem(); j++) {
      if (dK1[j].NumberOfFrequencies()) {
        const Numeric da = -0.5 * (r * dK1[j].Kjj(iz, ia)[i] +
                                   ((j == it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] +
                                                          K2.Kjj(iz, ia)[i])
                                              : 0.0)),
                      db = -0.5 * (r * dK1[j].K12(iz, ia)[i] +
                                   ((j == it) ? dr_dT1 * (K1.K12(iz, ia)[i] +
                                                          K2.K12(iz, ia)[i])
                                              : 0.0));
        dT1[j].Mat2(i).noalias() =
            T.Mat2(i) * da +
            (Eigen::Matrix2d() << sb, cb, cb, sb).finished() * exp_a * db;
      }
      if (dK2[j].NumberOfFrequencies()) {
        const Numeric da = -0.5 * (r * dK2[j].Kjj(iz, ia)[i] +
                                   ((j == it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] +
                                                          K2.Kjj(iz, ia)[i])
                                              : 0.0)),
                      db = -0.5 * (r * dK2[j].K12(iz, ia)[i] +
                                   ((j == it) ? dr_dT2 * (K1.K12(iz, ia)[i] +
                                                          K2.K12(iz, ia)[i])
                                              : 0.0));
        dT2[j].Mat2(i).noalias() =
            T.Mat2(i) * da +
            (Eigen::Matrix2d() << sb, cb, cb, sb).finished() * exp_a * db;
      }
    }
  }
}

inline void dtransmat3(TransmissionMatrix& T,
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
                       const Index ia) noexcept {
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]),
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]),
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]),
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);

    if (b == 0. and c == 0. and u == 0.) {
      T.Mat3(i).noalias() = Eigen::Matrix3d::Identity() * exp_a;
      for (Index j = 0; j < dT1.nelem(); j++) {
        if (dK1[j].NumberOfFrequencies())
          dT1[j].Mat3(i).noalias() =
              T.Mat3(i) *
              (-0.5 *
               (r * dK1[j].Kjj(iz, ia)[i] +
                ((j == it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])
                           : 0.0)));
        if (dK2[j].NumberOfFrequencies())
          dT2[j].Mat3(i).noalias() =
              T.Mat3(i) *
              (-0.5 *
               (r * dK2[j].Kjj(iz, ia)[i] +
                ((j == it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])
                           : 0.0)));
      }
    } else {
      const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u;
      const Numeric Const = b2 + c2 - u2;

      const bool real = Const > 0;
      const bool imag = Const < 0;
      const bool either = real or imag;

      const Numeric x =
          std::sqrt(imag ? -Const : Const);  // test to just use real values
      const Numeric x2 =
          (real ? 1 : -1) * x * x;  // test to change sign if imaginary
      const Numeric inv_x =
          either ? 1.0 / x
                 : 1.0;  // test so further calculations are replaced as x→0
      const Numeric inv_x2 = inv_x * inv_x;

      const Numeric sx =
          real ? std::sinh(x) : std::sin(-x);  // -i sin(ix) → sinh(x)
      const Numeric cx =
          real ? std::cosh(x) : std::cos(+x);  //    cos(ix) → cosh(x)

      /* Using:
       *    lim x→0 [(cosh(x) - 1) / x^2] → 1/2
       *    lim x→0 [sinh(x) / x]  → 1
       *    inv_x2 := 1 for x == 0,
       *    C0, C1, C2 ∝ [1/x^2]
       */
      const Numeric C0 =
          either ? a2 * (cx - 1.0) - a * x * sx + x2 : 1.0 + 0.5 * a2 - a;
      const Numeric C1 = either ? 2.0 * a * (1.0 - cx) + x * sx : 1.0 - a;
      const Numeric C2 = either ? cx - 1.0 : 0.5;

      T.Mat3(i).noalias() =
          exp_a * inv_x2 *
          (Eigen::Matrix3d() << C0 + C1 * a + C2 * (a2 + b2 + c2),
           C1 * b + C2 * (2 * a * b - c * u),
           C1 * c + C2 * (2 * a * c + b * u),
           C1 * b + C2 * (2 * a * b + c * u),
           C0 + C1 * a + C2 * (a2 + b2 - u2),
           C1 * u + C2 * (2 * a * u + b * c),
           C1 * c + C2 * (2 * a * c - b * u),
           -C1 * u - C2 * (2 * a * u - b * c),
           C0 + C1 * a + C2 * (a2 + c2 - u2))
              .finished();

      for (Index j = 0; j < dT1.nelem(); j++) {
        if (dK1[j].NumberOfFrequencies()) {
          const Numeric da = -0.5 * (r * dK1[j].Kjj(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] +
                                                            K2.Kjj(iz, ia)[i])
                                                : 0.0)),
                        db = -0.5 * (r * dK1[j].K12(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K12(iz, ia)[i] +
                                                            K2.K12(iz, ia)[i])
                                                : 0.0)),
                        dc = -0.5 * (r * dK1[j].K13(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K13(iz, ia)[i] +
                                                            K2.K13(iz, ia)[i])
                                                : 0.0)),
                        du = -0.5 * (r * dK1[j].K23(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K23(iz, ia)[i] +
                                                            K2.K23(iz, ia)[i])
                                                : 0.0));
          const Numeric da2 = 2 * a * da, db2 = 2 * b * db, dc2 = 2 * c * dc,
                        du2 = 2 * u * du;
          const Numeric dx = either ? ((db2 + dc2 - du2) * inv_x * 0.5) : 0,
                        dx2 = 2 * x * dx;
          const Numeric dsx = (real ? 1 : -1) * cx * dx;
          const Numeric dcx = sx * dx;

          const Numeric dC0 = either ? da2 * (cx - 1.0) + da2 * (dcx - 1.0) -
                                           da * x * sx - a * dx * sx -
                                           a * x * dsx + dx2
                                     : 0.5 * da2 - da;
          const Numeric dC1 =
              either ? 2.0 * (da * (1.0 - cx) - a * dcx) + dx * sx + x * dsx
                     : -da;
          const Numeric dC2 = either ? dcx : 0;

          dT1[j].Mat3(i).noalias() =
              T.Mat3(i) * (da + dx2 * inv_x2) +
              exp_a * inv_x2 *
                  (Eigen::Matrix3d() << dC0 + dC1 * a + C1 * da +
                                            dC2 * (a2 + b2 + c2) +
                                            C2 * (da2 + db2 + dc2),
                   dC1 * b + C1 * db + dC2 * (2 * a * b - c * u) +
                       C2 * (2 * da * b - dc * u + 2 * a * db - c * du),
                   dC1 * c + C1 * dc + dC2 * (2 * a * c + b * u) +
                       C2 * (2 * da * c + db * u + 2 * a * dc + b * du),
                   dC1 * b + C1 * db + dC2 * (2 * a * b + c * u) +
                       C2 * (2 * da * b + dc * u + 2 * a * db + c * du),
                   dC0 + dC1 * a + C1 * da + dC2 * (a2 + b2 - u2) +
                       C2 * (da2 + db2 - du2),
                   dC1 * u + C1 * du + dC2 * (2 * a * u + b * c) +
                       C2 * (2 * da * u + db * c + 2 * a * du + b * dc),
                   dC1 * c + C1 * dc + dC2 * (2 * a * c - b * u) +
                       C2 * (2 * da * c - db * u + 2 * a * dc - b * du),
                   -dC1 * u - C1 * du - dC2 * (2 * a * u - b * c) -
                       C2 * (2 * da * u - db * c + 2 * a * du - b * dc),
                   dC0 + dC1 * a + C1 * da + dC2 * (a2 + c2 - u2) +
                       C2 * (da2 + dc2 - du2))
                      .finished();
        }
        if (dK2[j].NumberOfFrequencies()) {
          const Numeric da = -0.5 * (r * dK2[j].Kjj(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] +
                                                            K2.Kjj(iz, ia)[i])
                                                : 0.0)),
                        db = -0.5 * (r * dK2[j].K12(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K12(iz, ia)[i] +
                                                            K2.K12(iz, ia)[i])
                                                : 0.0)),
                        dc = -0.5 * (r * dK2[j].K13(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K13(iz, ia)[i] +
                                                            K2.K13(iz, ia)[i])
                                                : 0.0)),
                        du = -0.5 * (r * dK2[j].K23(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K23(iz, ia)[i] +
                                                            K2.K23(iz, ia)[i])
                                                : 0.0));
          const Numeric da2 = 2 * a * da, db2 = 2 * b * db, dc2 = 2 * c * dc,
                        du2 = 2 * u * du;
          const Numeric dx = either ? ((db2 + dc2 - du2) * inv_x * 0.5) : 0,
                        dx2 = 2 * x * dx;
          const Numeric dsx = (real ? 1 : -1) * cx * dx;
          const Numeric dcx = sx * dx;

          const Numeric dC0 = either ? da2 * (cx - 1.0) + da2 * (dcx - 1.0) -
                                           da * x * sx - a * dx * sx -
                                           a * x * dsx + dx2
                                     : 0.5 * da2 - da;
          const Numeric dC1 =
              either ? 2.0 * (da * (1.0 - cx) - a * dcx) + dx * sx + x * dsx
                     : -da;
          const Numeric dC2 = either ? dcx : 0;

          dT2[j].Mat3(i).noalias() =
              T.Mat3(i) * (da + dx2 * inv_x2) +
              exp_a * inv_x2 *
                  (Eigen::Matrix3d() << dC0 + dC1 * a + C1 * da +
                                            dC2 * (a2 + b2 + c2) +
                                            C2 * (da2 + db2 + dc2),
                   dC1 * b + C1 * db + dC2 * (2 * a * b - c * u) +
                       C2 * (2 * da * b - dc * u + 2 * a * db - c * du),
                   dC1 * c + C1 * dc + dC2 * (2 * a * c + b * u) +
                       C2 * (2 * da * c + db * u + 2 * a * dc + b * du),
                   dC1 * b + C1 * db + dC2 * (2 * a * b + c * u) +
                       C2 * (2 * da * b + dc * u + 2 * a * db + c * du),
                   dC0 + dC1 * a + C1 * da + dC2 * (a2 + b2 - u2) +
                       C2 * (da2 + db2 - du2),
                   dC1 * u + C1 * du + dC2 * (2 * a * u + b * c) +
                       C2 * (2 * da * u + db * c + 2 * a * du + b * dc),
                   dC1 * c + C1 * dc + dC2 * (2 * a * c - b * u) +
                       C2 * (2 * da * c - db * u + 2 * a * dc - b * du),
                   -dC1 * u - C1 * du - dC2 * (2 * a * u - b * c) -
                       C2 * (2 * da * u - db * c + 2 * a * du - b * dc),
                   dC0 + dC1 * a + C1 * da + dC2 * (a2 + c2 - u2) +
                       C2 * (da2 + dc2 - du2))
                      .finished();
        }
      }
    }
  }
}

inline void dtransmat4(TransmissionMatrix& T,
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
                       const Index ia) noexcept {
  static constexpr Numeric sqrt_05 = Constant::inv_sqrt_2;
  for (Index i = 0; i < K1.NumberOfFrequencies(); i++) {
    const Numeric a = -0.5 * r * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i]),
                  b = -0.5 * r * (K1.K12(iz, ia)[i] + K2.K12(iz, ia)[i]),
                  c = -0.5 * r * (K1.K13(iz, ia)[i] + K2.K13(iz, ia)[i]),
                  d = -0.5 * r * (K1.K14(iz, ia)[i] + K2.K14(iz, ia)[i]),
                  u = -0.5 * r * (K1.K23(iz, ia)[i] + K2.K23(iz, ia)[i]),
                  v = -0.5 * r * (K1.K24(iz, ia)[i] + K2.K24(iz, ia)[i]),
                  w = -0.5 * r * (K1.K34(iz, ia)[i] + K2.K34(iz, ia)[i]);
    const Numeric exp_a = std::exp(a);

    if (b == 0. and c == 0. and d == 0. and u == 0. and v == 0. and w == 0.) {
      T.Mat4(i).noalias() = Eigen::Matrix4d::Identity() * exp_a;
      for (Index j = 0; j < dK1.nelem(); j++) {
        if (dK1[j].NumberOfFrequencies())
          dT1[j].Mat4(i).noalias() =
              T.Mat4(i) *
              (-0.5 *
               (r * dK1[j].Kjj(iz, ia)[i] +
                ((j == it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])
                           : 0.0)));
        if (dK2[j].NumberOfFrequencies())
          dT2[j].Mat4(i).noalias() =
              T.Mat4(i) *
              (-0.5 *
               (r * dK2[j].Kjj(iz, ia)[i] +
                ((j == it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] + K2.Kjj(iz, ia)[i])
                           : 0.0)));
      }
    } else {
      const Numeric b2 = b * b, c2 = c * c, d2 = d * d, u2 = u * u, v2 = v * v,
                    w2 = w * w;
      const Numeric tmp =
          w2 * w2 + 2 * (b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                         c2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                         d2 * (d2 * 0.5 + u2 - v2 - w2) +
                         u2 * (u2 * 0.5 + v2 + w2) + v2 * (v2 * 0.5 + w2) +
                         4 * (b * d * u * w - b * c * v * w - c * d * u * v));
      const Complex Const1 = std::sqrt(Complex(tmp, 0));
      const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
      const Complex tmp_x_sqrt = std::sqrt(Const2 + Const1);
      const Complex tmp_y_sqrt = std::sqrt(Const2 - Const1);
      const Complex x = tmp_x_sqrt * sqrt_05;
      const Complex y = tmp_y_sqrt * sqrt_05 * Complex(0, 1);
      const Complex x2 = x * x;
      const Complex y2 = y * y;
      const Complex cy = std::cos(y);
      const Complex sy = std::sin(y);
      const Complex cx = std::cosh(x);
      const Complex sx = std::sinh(x);

      const bool x_zero = std::abs(x) < lower_is_considered_zero_for_sinc_likes;
      const bool y_zero = std::abs(y) < lower_is_considered_zero_for_sinc_likes;
      const bool both_zero = y_zero and x_zero;
      const bool either_zero = y_zero or x_zero;

      /* Using:
       *    lim x→0 [({cosh(x),cos(x)} - 1) / x^2] → 1/2
       *    lim x→0 [{sinh(x),sin(x)} / x]  → 1
       *    inv_x2 := 1 for x == 0,
       *    -i sin(ix) → sinh(x)
       *    cos(ix) → cosh(x)
       *    C0, C1, C2 ∝ [1/x^2]
       */
      const Complex ix = x_zero ? 0.0 : 1.0 / x;
      const Complex iy = y_zero ? 0.0 : 1.0 / y;
      const Complex inv_x2y2 =
          both_zero
              ? 1.0
              : 1.0 /
                    (x2 + y2);  // The first "1.0" is the trick for above limits
      const Complex C0c = either_zero ? 1.0 : (cy * x2 + cx * y2) * inv_x2y2;
      const Complex C1c =
          either_zero ? 1.0 : (sy * x2 * iy + sx * y2 * ix) * inv_x2y2;
      const Complex C2c = both_zero ? 0.5 : (cx - cy) * inv_x2y2;
      const Complex C3c = both_zero ? 1.0 / 6.0
                                    : (x_zero   ? 1.0 - sy * iy
                                       : y_zero ? sx * ix - 1.0
                                                : sx * ix - sy * iy) *
                                          inv_x2y2;

      const Numeric& C0 = real_val(C0c);
      const Numeric& C1 = real_val(C1c);
      const Numeric& C2 = real_val(C2c);
      const Numeric& C3 = real_val(C3c);
      T.Mat4(i).noalias() =
          exp_a * (Eigen::Matrix4d() << C0 + C2 * (b2 + c2 + d2),
                   C1 * b + C2 * (-c * u - d * v) +
                       C3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                             v * (b * v + c * w)),
                   C1 * c + C2 * (b * u - d * w) +
                       C3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                             w * (b * v + c * w)),
                   C1 * d + C2 * (b * v + c * w) +
                       C3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                             w * (b * u - d * w)),

                   C1 * b + C2 * (c * u + d * v) +
                       C3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                             d * (b * d + u * w)),
                   C0 + C2 * (b2 - u2 - v2),
                   C2 * (b * c - v * w) + C1 * u +
                       C3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                             w * (b * d + u * w)),
                   C2 * (b * d + u * w) + C1 * v +
                       C3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                             w * (b * c - v * w)),

                   C1 * c + C2 * (-b * u + d * w) +
                       C3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                             d * (c * d - u * v)),
                   C2 * (b * c - v * w) - C1 * u +
                       C3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                             v * (c * d - u * v)),
                   C0 + C2 * (c2 - u2 - w2),
                   C2 * (c * d - u * v) + C1 * w +
                       C3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                             w * (-c2 + u2 + w2)),

                   C1 * d + C2 * (-b * v - c * w) +
                       C3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                             d * (-d2 + v2 + w2)),
                   C2 * (b * d + u * w) - C1 * v +
                       C3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                             v * (-d2 + v2 + w2)),
                   C2 * (c * d - u * v) - C1 * w +
                       C3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                             w * (-d2 + v2 + w2)),
                   C0 + C2 * (d2 - v2 - w2))
                      .finished();

      for (Index j = 0; j < dK1.nelem(); j++) {
        if (dK1[j].NumberOfFrequencies()) {
          const Numeric da = -0.5 * (r * dK1[j].Kjj(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.Kjj(iz, ia)[i] +
                                                            K2.Kjj(iz, ia)[i])
                                                : 0.0)),
                        db = -0.5 * (r * dK1[j].K12(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K12(iz, ia)[i] +
                                                            K2.K12(iz, ia)[i])
                                                : 0.0)),
                        dc = -0.5 * (r * dK1[j].K13(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K13(iz, ia)[i] +
                                                            K2.K13(iz, ia)[i])
                                                : 0.0)),
                        dd = -0.5 * (r * dK1[j].K14(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K14(iz, ia)[i] +
                                                            K2.K14(iz, ia)[i])
                                                : 0.0)),
                        du = -0.5 * (r * dK1[j].K23(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K23(iz, ia)[i] +
                                                            K2.K23(iz, ia)[i])
                                                : 0.0)),
                        dv = -0.5 * (r * dK1[j].K24(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K24(iz, ia)[i] +
                                                            K2.K24(iz, ia)[i])
                                                : 0.0)),
                        dw = -0.5 * (r * dK1[j].K34(iz, ia)[i] +
                                     ((j == it) ? dr_dT1 * (K1.K34(iz, ia)[i] +
                                                            K2.K34(iz, ia)[i])
                                                : 0.0));
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                        du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          const Numeric dtmp =
              2 * w2 * dw2 +
              2 * (db2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                   b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2) +
                   dc2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                   c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2) +
                   dd2 * (d2 * 0.5 + u2 - v2 - w2) +
                   d2 * (dd2 * 0.5 + du2 - dv2 - dw2) +
                   du2 * (u2 * 0.5 + v2 + w2) + u2 * (du2 * 0.5 + dv2 + dw2) +
                   dv2 * (v2 * 0.5 + w2) + v2 * (dv2 * 0.5 + dw2) +
                   4 * (db * d * u * w - db * c * v * w - dc * d * u * v +
                        b * dd * u * w - b * dc * v * w - c * dd * u * v +
                        b * d * du * w - b * c * dv * w - c * d * du * v +
                        b * d * u * dw - b * c * v * dw - c * d * u * dv));
          const Complex dConst1 = 0.5 * dtmp / Const1;
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          const Complex dx =
              x_zero ? 0 : (0.5 * (dConst2 + dConst1) / tmp_x_sqrt) * sqrt_05;
          const Complex dy = y_zero ? 0
                                    : (0.5 * (dConst2 - dConst1) / tmp_y_sqrt) *
                                          sqrt_05 * Complex(0, 1);
          const Complex dx2 = 2 * x * dx;
          const Complex dy2 = 2 * y * dy;
          const Complex dcy = -sy * dy;
          const Complex dsy = cy * dy;
          const Complex dcx = sx * dx;
          const Complex dsx = cx * dx;
          const Complex dix = -dx * ix * ix;
          const Complex diy = -dy * iy * iy;
          const Complex dx2dy2 = dx2 + dy2;
          const Complex dC0c =
              either_zero
                  ? 0.0
                  : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0c * dx2dy2) *
                        inv_x2y2;
          const Complex dC1c =
              either_zero ? 0.0
                          : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy +
                             dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix -
                             C1c * dx2dy2) *
                                inv_x2y2;
          const Complex dC2c =
              both_zero ? 0.0 : (dcx - dcy - C2c * dx2dy2) * inv_x2y2;
          const Complex dC3c =
              both_zero
                  ? 0.0
                  : ((x_zero   ? -dsy * iy - sy * diy
                      : y_zero ? dsx * ix + sx * dix
                               : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                     C3c * dx2dy2) *
                        inv_x2y2;

          const Numeric& dC0 = real_val(dC0c);
          const Numeric& dC1 = real_val(dC1c);
          const Numeric& dC2 = real_val(dC2c);
          const Numeric& dC3 = real_val(dC3c);
          dT1[j].Mat4(i).noalias() =
              T.Mat4(i) * da +
              exp_a *
                  (Eigen::Matrix4d()
                       << dC0 + dC2 * (b2 + c2 + d2) + C2 * (db2 + dc2 + dd2),
                   db * C1 + b * dC1 + dC2 * (-c * u - d * v) +
                       C2 * (-dc * u - dd * v - c * du - d * dv) +
                       dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                              v * (b * v + c * w)) +
                       C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                             dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                             u * (db * u - dd * w) - v * (db * v + dc * w) -
                             u * (b * du - d * dw) - v * (b * dv + c * dw)),
                   dC1 * c + C1 * dc + dC2 * (b * u - d * w) +
                       C2 * (db * u - dd * w + b * du - d * dw) +
                       dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                              w * (b * v + c * w)) +
                       C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                             dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                             u * (dc * u + dd * v) - w * (db * v + dc * w) -
                             u * (c * du + d * dv) - w * (b * dv + c * dw)),
                   dC1 * d + C1 * dd + dC2 * (b * v + c * w) +
                       C2 * (db * v + dc * w + b * dv + c * dw) +
                       dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                              w * (b * u - d * w)) +
                       C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                             dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                             v * (dc * u + dd * v) + w * (db * u - dd * w) -
                             v * (c * du + d * dv) + w * (b * du - d * dw)),

                   db * C1 + b * dC1 + dC2 * (c * u + d * v) +
                       C2 * (dc * u + dd * v + c * du + d * dv) +
                       dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                              d * (b * d + u * w)) +
                       C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                             dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                             c * (db * c - dv * w) + d * (db * d + du * w) +
                             c * (b * dc - v * dw) + d * (b * dd + u * dw)),
                   dC0 + dC2 * (b2 - u2 - v2) + C2 * (db2 - du2 - dv2),
                   dC2 * (b * c - v * w) +
                       C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u +
                       C1 * du +
                       dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                              w * (b * d + u * w)) +
                       C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                             dw * (b * d + u * w) + c * (dc * u + dd * v) -
                             u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                             c * (c * du + d * dv) - w * (b * dd + u * dw)),
                   dC2 * (b * d + u * w) +
                       C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v +
                       C1 * dv +
                       dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                              w * (b * c - v * w)) +
                       C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                             dw * (b * c - v * w) + d * (dc * u + dd * v) -
                             v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                             d * (c * du + d * dv) + w * (b * dc - v * dw)),

                   dC1 * c + C1 * dc + dC2 * (-b * u + d * w) +
                       C2 * (-db * u + dd * w - b * du + d * dw) +
                       dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                              d * (c * d - u * v)) +
                       C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                             dd * (c * d - u * v) + b * (db * c - dv * w) -
                             c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                             b * (b * dc - v * dw) + d * (c * dd - u * dv)),
                   dC2 * (b * c - v * w) +
                       C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u -
                       C1 * du +
                       dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                              v * (c * d - u * v)) +
                       C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                             dv * (c * d - u * v) - b * (db * u - dd * w) +
                             u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                             b * (b * du - d * dw) - v * (c * dd - u * dv)),
                   dC0 + dC2 * (c2 - u2 - w2) + C2 * (dc2 - du2 - dw2),
                   dC2 * (c * d - u * v) +
                       C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w +
                       C1 * dw +
                       dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                              w * (-c2 + u2 + w2)) +
                       C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                             dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                             v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                             d * (b * du - d * dw) + v * (b * dc - v * dw)),

                   dC1 * d + C1 * dd + dC2 * (-b * v - c * w) +
                       C2 * (-db * v - dc * w - b * dv - c * dw) +
                       dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                              d * (-d2 + v2 + w2)) +
                       C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                             dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                             c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                             b * (b * dd + u * dw) + c * (c * dd - u * dv)),
                   dC2 * (b * d + u * w) +
                       C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v -
                       C1 * dv +
                       dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                              v * (-d2 + v2 + w2)) +
                       C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                             dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                             u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                             b * (b * dv + c * dw) - u * (c * dd - u * dv)),
                   dC2 * (c * d - u * v) +
                       C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w -
                       C1 * dw +
                       dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                              w * (-d2 + v2 + w2)) +
                       C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                             dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                             u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                             c * (b * dv + c * dw) + u * (b * dd + u * dw)),
                   dC0 + dC2 * (d2 - v2 - w2) + C2 * (dd2 - dv2 - dw2))
                      .finished();
        }
        if (dK2[j].NumberOfFrequencies()) {
          const Numeric da = -0.5 * (r * dK2[j].Kjj(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.Kjj(iz, ia)[i] +
                                                            K2.Kjj(iz, ia)[i])
                                                : 0.0)),
                        db = -0.5 * (r * dK2[j].K12(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K12(iz, ia)[i] +
                                                            K2.K12(iz, ia)[i])
                                                : 0.0)),
                        dc = -0.5 * (r * dK2[j].K13(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K13(iz, ia)[i] +
                                                            K2.K13(iz, ia)[i])
                                                : 0.0)),
                        dd = -0.5 * (r * dK2[j].K14(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K14(iz, ia)[i] +
                                                            K2.K14(iz, ia)[i])
                                                : 0.0)),
                        du = -0.5 * (r * dK2[j].K23(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K23(iz, ia)[i] +
                                                            K2.K23(iz, ia)[i])
                                                : 0.0)),
                        dv = -0.5 * (r * dK2[j].K24(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K24(iz, ia)[i] +
                                                            K2.K24(iz, ia)[i])
                                                : 0.0)),
                        dw = -0.5 * (r * dK2[j].K34(iz, ia)[i] +
                                     ((j == it) ? dr_dT2 * (K1.K34(iz, ia)[i] +
                                                            K2.K34(iz, ia)[i])
                                                : 0.0));
          const Numeric db2 = 2 * db * b, dc2 = 2 * dc * c, dd2 = 2 * dd * d,
                        du2 = 2 * du * u, dv2 = 2 * dv * v, dw2 = 2 * dw * w;
          const Numeric dtmp =
              2 * w2 * dw2 +
              2 * (db2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2) +
                   b2 * (db2 * 0.5 + dc2 + dd2 - du2 - dv2 + dw2) +
                   dc2 * (c2 * 0.5 + d2 - u2 + v2 - w2) +
                   c2 * (dc2 * 0.5 + dd2 - du2 + dv2 - dw2) +
                   dd2 * (d2 * 0.5 + u2 - v2 - w2) +
                   d2 * (dd2 * 0.5 + du2 - dv2 - dw2) +
                   du2 * (u2 * 0.5 + v2 + w2) + u2 * (du2 * 0.5 + dv2 + dw2) +
                   dv2 * (v2 * 0.5 + w2) + v2 * (dv2 * 0.5 + dw2) +
                   4 * (db * d * u * w - db * c * v * w - dc * d * u * v +
                        b * dd * u * w - b * dc * v * w - c * dd * u * v +
                        b * d * du * w - b * c * dv * w - c * d * du * v +
                        b * d * u * dw - b * c * v * dw - c * d * u * dv));
          const Complex dConst1 = 0.5 * dtmp / Const1;
          const Numeric dConst2 = db2 + dc2 + dd2 - du2 - dv2 - dw2;
          const Complex dx =
              x_zero ? 0 : (0.5 * (dConst2 + dConst1) / tmp_x_sqrt) * sqrt_05;
          const Complex dy = y_zero ? 0
                                    : (0.5 * (dConst2 - dConst1) / tmp_y_sqrt) *
                                          sqrt_05 * Complex(0, 1);
          const Complex dx2 = 2 * x * dx;
          const Complex dy2 = 2 * y * dy;
          const Complex dcy = -sy * dy;
          const Complex dsy = cy * dy;
          const Complex dcx = sx * dx;
          const Complex dsx = cx * dx;
          const Complex dix = -dx * ix * ix;
          const Complex diy = -dy * iy * iy;
          const Complex dx2dy2 = dx2 + dy2;
          const Complex dC0c =
              either_zero
                  ? 0.0
                  : (dcy * x2 + cy * dx2 + dcx * y2 + cx * dy2 - C0c * dx2dy2) *
                        inv_x2y2;
          const Complex dC1c =
              either_zero ? 0.0
                          : (dsy * x2 * iy + sy * dx2 * iy + sy * x2 * diy +
                             dsx * y2 * ix + sx * dy2 * ix + sx * y2 * dix -
                             C1c * dx2dy2) *
                                inv_x2y2;
          const Complex dC2c =
              both_zero ? 0.0 : (dcx - dcy - C2c * dx2dy2) * inv_x2y2;
          const Complex dC3c =
              both_zero
                  ? 0.0
                  : ((x_zero   ? -dsy * iy - sy * diy
                      : y_zero ? dsx * ix + sx * dix
                               : dsx * ix + sx * dix - dsy * iy - sy * diy) -
                     C3c * dx2dy2) *
                        inv_x2y2;

          const Numeric& dC0 = real_val(dC0c);
          const Numeric& dC1 = real_val(dC1c);
          const Numeric& dC2 = real_val(dC2c);
          const Numeric& dC3 = real_val(dC3c);
          dT2[j].Mat4(i).noalias() =
              T.Mat4(i) * da +
              exp_a *
                  (Eigen::Matrix4d()
                       << dC0 + dC2 * (b2 + c2 + d2) + C2 * (db2 + dc2 + dd2),
                   db * C1 + b * dC1 + dC2 * (-c * u - d * v) +
                       C2 * (-dc * u - dd * v - c * du - d * dv) +
                       dC3 * (b * (b2 + c2 + d2) - u * (b * u - d * w) -
                              v * (b * v + c * w)) +
                       C3 * (db * (b2 + c2 + d2) - du * (b * u - d * w) -
                             dv * (b * v + c * w) + b * (db2 + dc2 + dd2) -
                             u * (db * u - dd * w) - v * (db * v + dc * w) -
                             u * (b * du - d * dw) - v * (b * dv + c * dw)),
                   dC1 * c + C1 * dc + dC2 * (b * u - d * w) +
                       C2 * (db * u - dd * w + b * du - d * dw) +
                       dC3 * (c * (b2 + c2 + d2) - u * (c * u + d * v) -
                              w * (b * v + c * w)) +
                       C3 * (dc * (b2 + c2 + d2) - du * (c * u + d * v) -
                             dw * (b * v + c * w) + c * (db2 + dc2 + dd2) -
                             u * (dc * u + dd * v) - w * (db * v + dc * w) -
                             u * (c * du + d * dv) - w * (b * dv + c * dw)),
                   dC1 * d + C1 * dd + dC2 * (b * v + c * w) +
                       C2 * (db * v + dc * w + b * dv + c * dw) +
                       dC3 * (d * (b2 + c2 + d2) - v * (c * u + d * v) +
                              w * (b * u - d * w)) +
                       C3 * (dd * (b2 + c2 + d2) - dv * (c * u + d * v) +
                             dw * (b * u - d * w) + d * (db2 + dc2 + dd2) -
                             v * (dc * u + dd * v) + w * (db * u - dd * w) -
                             v * (c * du + d * dv) + w * (b * du - d * dw)),

                   db * C1 + b * dC1 + dC2 * (c * u + d * v) +
                       C2 * (dc * u + dd * v + c * du + d * dv) +
                       dC3 * (-b * (-b2 + u2 + v2) + c * (b * c - v * w) +
                              d * (b * d + u * w)) +
                       C3 * (-db * (-b2 + u2 + v2) + dc * (b * c - v * w) +
                             dd * (b * d + u * w) - b * (-db2 + du2 + dv2) +
                             c * (db * c - dv * w) + d * (db * d + du * w) +
                             c * (b * dc - v * dw) + d * (b * dd + u * dw)),
                   dC0 + dC2 * (b2 - u2 - v2) + C2 * (db2 - du2 - dv2),
                   dC2 * (b * c - v * w) +
                       C2 * (db * c + b * dc - dv * w - v * dw) + dC1 * u +
                       C1 * du +
                       dC3 * (c * (c * u + d * v) - u * (-b2 + u2 + v2) -
                              w * (b * d + u * w)) +
                       C3 * (dc * (c * u + d * v) - du * (-b2 + u2 + v2) -
                             dw * (b * d + u * w) + c * (dc * u + dd * v) -
                             u * (-db2 + du2 + dv2) - w * (db * d + du * w) +
                             c * (c * du + d * dv) - w * (b * dd + u * dw)),
                   dC2 * (b * d + u * w) +
                       C2 * (db * d + b * dd + du * w + u * dw) + dC1 * v +
                       C1 * dv +
                       dC3 * (d * (c * u + d * v) - v * (-b2 + u2 + v2) +
                              w * (b * c - v * w)) +
                       C3 * (dd * (c * u + d * v) - dv * (-b2 + u2 + v2) +
                             dw * (b * c - v * w) + d * (dc * u + dd * v) -
                             v * (-db2 + du2 + dv2) + w * (db * c - dv * w) +
                             d * (c * du + d * dv) + w * (b * dc - v * dw)),

                   dC1 * c + C1 * dc + dC2 * (-b * u + d * w) +
                       C2 * (-db * u + dd * w - b * du + d * dw) +
                       dC3 * (b * (b * c - v * w) - c * (-c2 + u2 + w2) +
                              d * (c * d - u * v)) +
                       C3 * (db * (b * c - v * w) - dc * (-c2 + u2 + w2) +
                             dd * (c * d - u * v) + b * (db * c - dv * w) -
                             c * (-dc2 + du2 + dw2) + d * (dc * d - du * v) +
                             b * (b * dc - v * dw) + d * (c * dd - u * dv)),
                   dC2 * (b * c - v * w) +
                       C2 * (db * c + b * dc - dv * w - v * dw) - dC1 * u -
                       C1 * du +
                       dC3 * (-b * (b * u - d * w) + u * (-c2 + u2 + w2) -
                              v * (c * d - u * v)) +
                       C3 * (-db * (b * u - d * w) + du * (-c2 + u2 + w2) -
                             dv * (c * d - u * v) - b * (db * u - dd * w) +
                             u * (-dc2 + du2 + dw2) - v * (dc * d - du * v) -
                             b * (b * du - d * dw) - v * (c * dd - u * dv)),
                   dC0 + dC2 * (c2 - u2 - w2) + C2 * (dc2 - du2 - dw2),
                   dC2 * (c * d - u * v) +
                       C2 * (dc * d + c * dd - du * v - u * dv) + dC1 * w +
                       C1 * dw +
                       dC3 * (-d * (b * u - d * w) + v * (b * c - v * w) -
                              w * (-c2 + u2 + w2)) +
                       C3 * (-dd * (b * u - d * w) + dv * (b * c - v * w) -
                             dw * (-c2 + u2 + w2) - d * (db * u - dd * w) +
                             v * (db * c - dv * w) - w * (-dc2 + du2 + dw2) -
                             d * (b * du - d * dw) + v * (b * dc - v * dw)),

                   dC1 * d + C1 * dd + dC2 * (-b * v - c * w) +
                       C2 * (-db * v - dc * w - b * dv - c * dw) +
                       dC3 * (b * (b * d + u * w) + c * (c * d - u * v) -
                              d * (-d2 + v2 + w2)) +
                       C3 * (db * (b * d + u * w) + dc * (c * d - u * v) -
                             dd * (-d2 + v2 + w2) + b * (db * d + du * w) +
                             c * (dc * d - du * v) - d * (-dd2 + dv2 + dw2) +
                             b * (b * dd + u * dw) + c * (c * dd - u * dv)),
                   dC2 * (b * d + u * w) +
                       C2 * (db * d + b * dd + du * w + u * dw) - dC1 * v -
                       C1 * dv +
                       dC3 * (-b * (b * v + c * w) - u * (c * d - u * v) +
                              v * (-d2 + v2 + w2)) +
                       C3 * (-db * (b * v + c * w) - du * (c * d - u * v) +
                             dv * (-d2 + v2 + w2) - b * (db * v + dc * w) -
                             u * (dc * d - du * v) + v * (-dd2 + dv2 + dw2) -
                             b * (b * dv + c * dw) - u * (c * dd - u * dv)),
                   dC2 * (c * d - u * v) +
                       C2 * (dc * d + c * dd - du * v - u * dv) - dC1 * w -
                       C1 * dw +
                       dC3 * (-c * (b * v + c * w) + u * (b * d + u * w) +
                              w * (-d2 + v2 + w2)) +
                       C3 * (-dc * (b * v + c * w) + du * (b * d + u * w) +
                             dw * (-d2 + v2 + w2) - c * (db * v + dc * w) +
                             u * (db * d + du * w) + w * (-dd2 + dv2 + dw2) -
                             c * (b * dv + c * dw) + u * (b * dd + u * dw)),
                   dC0 + dC2 * (d2 - v2 - w2) + C2 * (dd2 - dv2 - dw2))
                      .finished();
        }
      }
    }
  }
}

inline void transmat(TransmissionMatrix& T,
                     const PropagationMatrix& K1,
                     const PropagationMatrix& K2,
                     const Numeric& r) noexcept {
  switch (K1.StokesDimensions()) {
    case 4:
      transmat4(T, K1, K2, r);
      break;
    case 3:
      transmat3(T, K1, K2, r);
      break;
    case 2:
      transmat2(T, K1, K2, r);
      break;
    case 1:
      transmat1(T, K1, K2, r);
      break;
  }
}

inline void dtransmat(TransmissionMatrix& T,
                      ArrayOfTransmissionMatrix& dT1,
                      ArrayOfTransmissionMatrix& dT2,
                      const PropagationMatrix& K1,
                      const PropagationMatrix& K2,
                      const ArrayOfPropagationMatrix& dK1,
                      const ArrayOfPropagationMatrix& dK2,
                      const Numeric& r,
                      const Numeric& dr_dT1 = 0,
                      const Numeric& dr_dT2 = 0,
                      const Index it = -1,
                      const Index iz = 0,
                      const Index ia = 0) noexcept {
  switch (K1.StokesDimensions()) {
    case 4:
      dtransmat4(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia);
      break;
    case 3:
      dtransmat3(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia);
      break;
    case 2:
      dtransmat2(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia);
      break;
    case 1:
      dtransmat1(T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dT1, dr_dT2, it, iz, ia);
      break;
  }
}

void stepwise_transmission(TransmissionMatrix& T,
                           ArrayOfTransmissionMatrix& dT1,
                           ArrayOfTransmissionMatrix& dT2,
                           const PropagationMatrix& K1,
                           const PropagationMatrix& K2,
                           const ArrayOfPropagationMatrix& dK1,
                           const ArrayOfPropagationMatrix& dK2,
                           const Numeric& r,
                           const Numeric& dr_dtemp1,
                           const Numeric& dr_dtemp2,
                           const Index temp_deriv_pos) {
  if (not dT1.nelem())
    transmat(T, K1, K2, r);
  else
    dtransmat(
        T, dT1, dT2, K1, K2, dK1, dK2, r, dr_dtemp1, dr_dtemp2, temp_deriv_pos);
}

void stepwise_source(RadiationVector &J, ArrayOfRadiationVector &dJ,
                     const RadiationVector &J_add, const PropagationMatrix &K,
                     const StokesVector &a, const StokesVector &S,
                     const ArrayOfPropagationMatrix &dK,
                     const ArrayOfStokesVector &da,
                     const ArrayOfStokesVector &dS, const ConstVectorView &B,
                     const ConstVectorView &dB_dT,
                     const ArrayOfRetrievalQuantity &jacobian_quantities,
                     const bool &jacobian_do) {
  for (Index i = 0; i < K.NumberOfFrequencies(); i++) {
    if (K.IsRotational(i)) {
      J.SetZero(i);
      if (jacobian_do) {
        for (Index j = 0; j < jacobian_quantities.nelem(); j++)
          if (dJ[j].Frequencies())
            dJ[j].SetZero(i);
      }
    } else {
      J.setSource(a, B, S, i);

      switch (J.stokes_dim) {
      case 4: {
        const auto invK = inv_prop_matrix<4>(K.Data()(0, 0, i, joker));
        J.Vec4(i) = invK * J.Vec4(i);
        if (J_add.Frequencies())
          J.Vec4(i).noalias() += invK * J_add.Vec4(i);
        // TODO: Add jacobians dJ_add of additional source

        if (jacobian_do) {
          for (Index j = 0; j < jacobian_quantities.nelem(); j++) {
            if (dJ[j].Frequencies() == da[j].NumberOfFrequencies() and
                dJ[j].Frequencies() == dS[j].NumberOfFrequencies()) {
              dJ[j].Vec4(i).noalias() =
                  0.5 * invK *
                  (source_vector<4>(a, B, da[j], dB_dT, dS[j],
                                    jacobian_quantities[j] ==
                                        Jacobian::Atm::Temperature,
                                    i) -
                   prop_matrix<4>(dK[j].Data()(0, 0, i, joker)) * J.Vec4(i));
            }
          }
        }
      } break;
      case 3: {
        const auto invK = inv_prop_matrix<3>(K.Data()(0, 0, i, joker));
        J.Vec3(i) = invK * J.Vec3(i);
        if (J_add.Frequencies())
          J.Vec3(i).noalias() += invK * J_add.Vec3(i);
        // TODO: Add jacobians dJ_add of additional source

        if (jacobian_do) {
          for (Index j = 0; j < jacobian_quantities.nelem(); j++) {
            // Skip others!
            if (dJ[j].Frequencies() == da[j].NumberOfFrequencies() and
                dJ[j].Frequencies() == dS[j].NumberOfFrequencies()) {
              dJ[j].Vec3(i).noalias() =
                  0.5 * invK *
                  (source_vector<3>(a, B, da[j], dB_dT, dS[j],
                                    jacobian_quantities[j] ==
                                        Jacobian::Atm::Temperature,
                                    i) -
                   prop_matrix<3>(dK[j].Data()(0, 0, i, joker)) * J.Vec3(i));
            }
          }
        }
      } break;
      case 2: {
        const auto invK = inv_prop_matrix<2>(K.Data()(0, 0, i, joker));
        J.Vec2(i) = invK * J.Vec2(i);
        if (J_add.Frequencies())
          J.Vec2(i).noalias() += invK * J_add.Vec2(i);
        // TODO: Add jacobians dJ_add of additional source

        if (jacobian_do) {
          for (Index j = 0; j < jacobian_quantities.nelem(); j++) {
            // Skip others!
            if (dJ[j].Frequencies() == da[j].NumberOfFrequencies() and
                dJ[j].Frequencies() == dS[j].NumberOfFrequencies()) {
              dJ[j].Vec2(i).noalias() =
                  0.5 * invK *
                  (source_vector<2>(a, B, da[j], dB_dT, dS[j],
                                    jacobian_quantities[j] ==
                                        Jacobian::Atm::Temperature,
                                    i) -
                   prop_matrix<2>(dK[j].Data()(0, 0, i, joker)) * J.Vec2(i));
            }
          }
        }
      } break;
      default: {
        const auto invK = inv_prop_matrix<1>(K.Data()(0, 0, i, joker));
        J.Vec1(i) = invK * J.Vec1(i);
        if (J_add.Frequencies())
          J.Vec1(i).noalias() += invK * J_add.Vec1(i);
        // TODO: Add jacobians dJ_add of additional source

        if (jacobian_do) {
          for (Index j = 0; j < jacobian_quantities.nelem(); j++) {
            // Skip others!
            if (dJ[j].Frequencies() == da[j].NumberOfFrequencies() and
                dJ[j].Frequencies() == dS[j].NumberOfFrequencies()) {
              dJ[j].Vec1(i).noalias() =
                  0.5 * invK *
                  (source_vector<1>(a, B, da[j], dB_dT, dS[j],
                                    jacobian_quantities[j] ==
                                        Jacobian::Atm::Temperature,
                                    i) -
                   prop_matrix<1>(dK[j].Data()(0, 0, i, joker)) * J.Vec1(i));
            }
          }
        }
      } break;
      }
    }
  }
}

void update_radiation_vector(
    RadiationVector& I,
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
    [[maybe_unused]] const PropagationMatrix& K1,
    [[maybe_unused]] const PropagationMatrix& K2,
    [[maybe_unused]] const ArrayOfPropagationMatrix& dK1,
    [[maybe_unused]] const ArrayOfPropagationMatrix& dK2,
    [[maybe_unused]] const Numeric r,
    [[maybe_unused]] const Vector& dr1,
    [[maybe_unused]] const Vector& dr2,
    [[maybe_unused]] const Index ia,
    [[maybe_unused]] const Index iz,
    const RadiativeTransferSolver solver) {
  switch (solver) {
    case RadiativeTransferSolver::Emission: {
      I.rem_avg(J1, J2);
      for (size_t i = 0; i < dI1.size(); i++) {
        dI1[i].addDerivEmission(PiT, dT1[i], T, I, dJ1[i]);
        dI2[i].addDerivEmission(PiT, dT2[i], T, I, dJ2[i]);
      }
      I.leftMul(T);
      I.add_avg(J1, J2);
    } break;

    case RadiativeTransferSolver::Transmission: {
      for (size_t i = 0; i < dI1.size(); i++) {
        dI1[i].addDerivTransmission(PiT, dT1[i], I);
        dI2[i].addDerivTransmission(PiT, dT2[i], I);
      }
      I.leftMul(T);
    } break;

    case RadiativeTransferSolver::LinearWeightedEmission: {
      ARTS_USER_ERROR_IF(
          dI1.size(),
          "Cannot support derivatives with current integration method\n");

      I.leftMul(T);
      I.add_weighted(T,
                     J1,
                     J2,
                     K1.Data()(ia, iz, joker, joker),
                     K2.Data()(ia, iz, joker, joker),
                     r);

    } break;
    case RadiativeTransferSolver::FINAL: {
      ARTS_ASSERT(false)
    }
  }
}

ArrayOfTransmissionMatrix cumulative_transmission(
    const ArrayOfTransmissionMatrix& T,
    const CumulativeTransmission type) /*[[expects: T.nelem()>0]]*/
{
  const Index n = T.nelem();
  const Index nf = n ? T.front().Frequencies() : 0;
  const Index ns = n ? T.front().stokes_dim : 1;
  ArrayOfTransmissionMatrix PiT(n, TransmissionMatrix(nf, ns));
  switch (type) {
    case CumulativeTransmission::Forward: {
      for (Index i = 1; i < n; i++) PiT[i].mul(PiT[i - 1], T[i]);
    } break;
    case CumulativeTransmission::Reverse: {
      for (Index i = 1; i < n; i++) PiT[i].mul(T[i], PiT[i - 1]);
    } break;
  }
  return PiT;  // Note how the output is such that forward transmission is from -1 to 0
}

// TEST CODE BEGIN

void set_backscatter_radiation_vector(
    ArrayOfRadiationVector& I,
    ArrayOfArrayOfArrayOfRadiationVector& dI,
    const RadiationVector& I_incoming,
    const ArrayOfTransmissionMatrix& T,
    const ArrayOfTransmissionMatrix& PiTf,
    const ArrayOfTransmissionMatrix& PiTr,
    const ArrayOfTransmissionMatrix& Z,
    const ArrayOfArrayOfTransmissionMatrix& dT1,
    const ArrayOfArrayOfTransmissionMatrix& dT2,
    const ArrayOfArrayOfTransmissionMatrix& dZ,
    const BackscatterSolver solver) {
  const Index np = I.nelem();
  const Index nv = np ? I[0].Frequencies() : 0;
  const Index ns = np ? I[0].stokes_dim : 1;
  const Index nq = np ? dI[0][0].nelem() : 0;

  // For all transmission, the I-vector is the same
  for (Index ip = 0; ip < np; ip++)
    I[ip].setBackscatterTransmission(I_incoming, PiTr[ip], PiTf[ip], Z[ip]);

  for (Index ip = 0; ip < np; ip++) {
    for (Index iq = 0; iq < nq; iq++) {
      dI[ip][ip][iq].setBackscatterTransmissionDerivative(
          I_incoming, PiTr[ip], PiTf[ip], dZ[ip][iq]);
    }
  }

  switch (solver) {
    case BackscatterSolver::CommutativeTransmission: {
      // Forward and backwards transmission derivatives
      switch (ns) {
        case 1: {
        BackscatterSolverCommutativeTransmissionStokesDimOne:
          for (Index ip = 0; ip < np; ip++) {
            for (Index j = ip; j < np; j++) {
              for (Index iq = 0; iq < nq; iq++) {
                for (Index iv = 0; iv < nv; iv++) {
                  dI[ip][j][iq].Vec1(iv).noalias() +=
                      T[ip].Mat1(iv).inverse() *
                      (dT1[ip][iq].Mat1(iv) + dT2[ip][iq].Mat1(iv)) *
                      I[j].Vec1(iv);

                  if (j < np - 1 and j > ip)
                    dI[ip][j][iq].Vec1(iv).noalias() +=
                        T[ip + 1].Mat1(iv).inverse() *
                        (dT1[ip][iq].Mat1(iv) + dT2[ip][iq].Mat1(iv)) *
                        I[j].Vec1(iv);
                }
              }
            }
          }
        } break;
        case 2: {
          for (Index ip = 0; ip < np; ip++) {
            for (Index j = ip; j < np; j++) {
              for (Index iq = 0; iq < nq; iq++) {
                for (Index iv = 0; iv < nv; iv++) {
                  dI[ip][j][iq].Vec2(iv).noalias() +=
                      T[ip].Mat2(iv).inverse() *
                      (dT1[ip][iq].Mat2(iv) + dT2[ip][iq].Mat2(iv)) *
                      I[j].Vec2(iv);

                  if (j < np - 1 and j > ip)
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                        T[ip + 1].Mat2(iv).inverse() *
                        (dT1[ip][iq].Mat2(iv) + dT2[ip][iq].Mat2(iv)) *
                        I[j].Vec2(iv);
                }
              }
            }
          }
        } break;
        case 3: {
          for (Index ip = 0; ip < np; ip++) {
            for (Index j = ip; j < np; j++) {
              for (Index iq = 0; iq < nq; iq++) {
                for (Index iv = 0; iv < nv; iv++) {
                  dI[ip][j][iq].Vec3(iv).noalias() +=
                      T[ip].Mat3(iv).inverse() *
                      (dT1[ip][iq].Mat3(iv) + dT2[ip][iq].Mat3(iv)) *
                      I[j].Vec3(iv);

                  if (j < np - 1 and j > ip)
                    dI[ip][j][iq].Vec3(iv).noalias() +=
                        T[ip + 1].Mat3(iv).inverse() *
                        (dT1[ip][iq].Mat3(iv) + dT2[ip][iq].Mat3(iv)) *
                        I[j].Vec3(iv);
                }
              }
            }
          }
        } break;
        case 4: {
          for (Index ip = 0; ip < np; ip++) {
            for (Index j = ip; j < np; j++) {
              for (Index iq = 0; iq < nq; iq++) {
                for (Index iv = 0; iv < nv; iv++) {
                  dI[ip][j][iq].Vec4(iv).noalias() +=
                      T[ip].Mat4(iv).inverse() *
                      (dT1[ip][iq].Mat4(iv) + dT2[ip][iq].Mat4(iv)) *
                      I[j].Vec4(iv);

                  if (j < np - 1 and j > ip)
                    dI[ip][j][iq].Vec4(iv).noalias() +=
                        T[ip + 1].Mat4(iv).inverse() *
                        (dT1[ip][iq].Mat4(iv) + dT2[ip][iq].Mat4(iv)) *
                        I[j].Vec4(iv);
                }
              }
            }
          }
        } break;
      }
    } break;
    case BackscatterSolver::FullTransmission: {
      std::runtime_error("Do not activate this code.  It is not ready yet");
      switch (ns) {
        case 1: {
          // This is the same as CommutativeTransmission, so use that code
          goto BackscatterSolverCommutativeTransmissionStokesDimOne;
        } break;
        case 2: {
          for (Index ip = 0; ip < np; ip++) {
            for (Index j = ip; j < np; j++) {
              for (Index iq = 0; iq < nq; iq++) {
                for (Index iv = 0; iv < nv; iv++) {
                  if (ip > 1) {
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                        PiTr[ip - 2].Mat2(iv) *
                        (dT1[ip - 1][iq].Mat2(iv) + dT2[ip][iq].Mat2(iv)) *
                        PiTr[ip - 1].Mat2(iv).inverse() * I[j].Vec2(iv);

                    if (j < np - 1)
                      dI[ip][j][iq].Vec2(iv).noalias() +=
                          PiTr[ip].Mat2(iv) * Z[ip].Mat2(iv) *
                          PiTf[ip - 2].Mat2(iv) *
                          (dT1[ip][iq].Mat2(iv) + dT2[ip + 1][iq].Mat2(iv)) *
                          PiTf[ip - 1].Mat2(iv).inverse() * PiTf[ip].Mat2(iv) *
                          I[0].Vec2(iv);
                  } else {
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                        (dT1[ip - 1][iq].Mat2(iv) + dT2[ip][iq].Mat2(iv)) *
                        PiTr[ip - 1].Mat2(iv).inverse() * I[j].Vec2(iv);

                    if (j < np - 1)
                      dI[ip][j][iq].Vec2(iv).noalias() +=
                          PiTr[ip].Mat2(iv) * Z[ip].Mat2(iv) *
                          (dT1[ip][iq].Mat2(iv) + dT2[ip + 1][iq].Mat2(iv)) *
                          PiTf[ip - 1].Mat2(iv).inverse() * PiTf[ip].Mat2(iv) *
                          I[0].Vec2(iv);
                  }
                }
              }
            }
          }
        } break;
        case 3: {
        } break;
        case 4: {
        } break;
      }
    } break;
  }
}

ArrayOfTransmissionMatrix bulk_backscatter(const ConstTensor5View& Pe,
                                           const ConstMatrixView& pnd) {
  const Index ns = Pe.ncols();
  const Index nv = Pe.npages();
  const Index np = Pe.nbooks();
  const Index ne = Pe.nshelves();

  ArrayOfTransmissionMatrix aotm(np, TransmissionMatrix(nv, ns));

  for (Index ip = 0; ip < np; ip++) {
    aotm[ip].setZero();

    switch (ns) {
      case 4:
        for (Index iv = 0; iv < nv; iv++)
          for (Index ie = 0; ie < ne; ie++)
            aotm[ip].Mat4(iv).noalias() +=
                pnd(ie, ip) * prop_matrix<4>(Pe(ie, ip, iv, joker, joker));
        break;
      case 3:
        for (Index iv = 0; iv < nv; iv++)
          for (Index ie = 0; ie < ne; ie++)
            aotm[ip].Mat3(iv).noalias() +=
                pnd(ie, ip) * prop_matrix<3>(Pe(ie, ip, iv, joker, joker));
        break;
      case 2:
        for (Index iv = 0; iv < nv; iv++)
          for (Index ie = 0; ie < ne; ie++)
            aotm[ip].Mat2(iv).noalias() +=
                pnd(ie, ip) * prop_matrix<2>(Pe(ie, ip, iv, joker, joker));
        break;
      case 1:
        for (Index iv = 0; iv < nv; iv++)
          for (Index ie = 0; ie < ne; ie++)
            aotm[ip].Mat1(iv).noalias() +=
                pnd(ie, ip) * prop_matrix<1>(Pe(ie, ip, iv, joker, joker));
        break;
    }
  }
  return aotm;
}

ArrayOfArrayOfTransmissionMatrix bulk_backscatter_derivative(
    const ConstTensor5View& Pe, const ArrayOfMatrix& dpnd_dx) {
  const Index ns = Pe.ncols();
  const Index nv = Pe.npages();
  const Index np = Pe.nbooks();
  const Index ne = Pe.nshelves();
  const Index nq = dpnd_dx.nelem();

  ArrayOfArrayOfTransmissionMatrix aoaotm(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nv, ns)));

  for (Index ip = 0; ip < np; ip++) {
    for (Index iq = 0; iq < nq; iq++) {
      aoaotm[ip][iq].setZero();

      if (not dpnd_dx[iq].empty()) {
        switch (ns) {
          case 4:
            for (Index iv = 0; iv < nv; iv++)
              for (Index ie = 0; ie < ne; ie++)
                aoaotm[ip][iq].Mat4(iv).noalias() +=
                    dpnd_dx[iq](ie, ip) *
                    prop_matrix<4>(Pe(ie, ip, iv, joker, joker));
            break;
          case 3:
            for (Index iv = 0; iv < nv; iv++)
              for (Index ie = 0; ie < ne; ie++)
                aoaotm[ip][iq].Mat3(iv).noalias() +=
                    dpnd_dx[iq](ie, ip) *
                    prop_matrix<3>(Pe(ie, ip, iv, joker, joker));
            break;
          case 2:
            for (Index iv = 0; iv < nv; iv++)
              for (Index ie = 0; ie < ne; ie++)
                aoaotm[ip][iq].Mat2(iv).noalias() +=
                    dpnd_dx[iq](ie, ip) *
                    prop_matrix<2>(Pe(ie, ip, iv, joker, joker));
            break;
          case 1:
            for (Index iv = 0; iv < nv; iv++)
              for (Index ie = 0; ie < ne; ie++)
                aoaotm[ip][iq].Mat1(iv).noalias() +=
                    dpnd_dx[iq](ie, ip) *
                    prop_matrix<1>(Pe(ie, ip, iv, joker, joker));
            break;
        }
      }
    }
  }
  return aoaotm;
}

// TEST CODE END

std::ostream& operator<<(std::ostream& os, const TransmissionMatrix& tm) {
  for (const auto& T : tm.T4) os << T << '\n';
  for (const auto& T : tm.T3) os << T << '\n';
  for (const auto& T : tm.T2) os << T << '\n';
  for (const auto& T : tm.T1) os << T << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const RadiationVector& rv) {
  // Write the transpose because it looks better...
  for (const auto& R : rv.R4) os << R.transpose() << '\n';
  for (const auto& R : rv.R3) os << R.transpose() << '\n';
  for (const auto& R : rv.R2) os << R.transpose() << '\n';
  for (const auto& R : rv.R1) os << R.transpose() << '\n';
  return os;
}

std::istream& operator>>(std::istream& is, TransmissionMatrix& tm) {
  for (auto& T : tm.T4)
    is >> double_imanip() >> T(0, 0) >> T(0, 1) >> T(0, 2) >> T(0, 3) >>
        T(1, 0) >> T(1, 1) >> T(1, 2) >> T(1, 3) >> T(2, 0) >> T(2, 1) >>
        T(2, 2) >> T(2, 3) >> T(3, 0) >> T(3, 1) >> T(3, 2) >> T(3, 3);
  for (auto& T : tm.T3)
    is >> double_imanip() >> T(0, 0) >> T(0, 1) >> T(0, 2) >> T(1, 0) >>
        T(1, 1) >> T(1, 2) >> T(2, 0) >> T(2, 1) >> T(2, 2);
  for (auto& T : tm.T2)
    is >> double_imanip() >> T(0, 0) >> T(0, 1) >> T(1, 0) >> T(1, 1);
  for (auto& T : tm.T1) is >> double_imanip() >> T(0, 0);
  return is;
}

std::istream& operator>>(std::istream& is, RadiationVector& rv) {
  for (auto& R : rv.R4) is >> double_imanip() >> R[0] >> R[1] >> R[2] >> R[3];
  for (auto& R : rv.R3) is >> double_imanip() >> R[0] >> R[1] >> R[2];
  for (auto& R : rv.R2) is >> double_imanip() >> R[0] >> R[1] >> R[2];
  for (auto& R : rv.R1) is >> double_imanip() >> R[0];
  return is;
}

TransmissionMatrix::TransmissionMatrix(const PropagationMatrix& pm,
                                       const Numeric& r) {
  *this = TransmissionMatrix(pm.NumberOfFrequencies(), pm.StokesDimensions());
  transmat(*this, pm, pm, r);  // Slower to compute, faster to implement...
}
