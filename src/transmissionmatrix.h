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

/**
 * @file   transmissionmatrix.h
 * @author Richard Larsson
 * @date   2018-01-30
 * 
 * @brief  Stuff related to the transmission matrix.
 * 
 * Using Eigen library to speed up computations.
 */

#ifndef transmissionmatrix_h
#define transmissionmatrix_h

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "jacobian.h"
#include "propagationmatrix.h"

/** Class to keep track of Transmission Matrices for Stokes Dim 1-4 */
class TransmissionMatrix {
 private:
  Index stokes_dim;
  std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>> T4;
  std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> T3;
  std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> T2;
  std::vector<Eigen::Matrix<double, 1, 1>,
              Eigen::aligned_allocator<Eigen::Matrix<double, 1, 1>>>
      T1;

 public:

  /** Construct a new Transmission Matrix object
   * 
   * @param[in] nf Number of frequencies
   * @param[in] stokes Stokes dimension
   */
  TransmissionMatrix(Index nf = 0, Index stokes = 1)
      : stokes_dim(stokes),
        T4(stokes_dim == 4 ? nf : 0, Eigen::Matrix4d::Identity()),
        T3(stokes_dim == 3 ? nf : 0, Eigen::Matrix3d::Identity()),
        T2(stokes_dim == 2 ? nf : 0, Eigen::Matrix2d::Identity()),
        T1(stokes_dim == 1 ? nf : 0, Eigen::Matrix<double, 1, 1>::Identity()) {
    assert(stokes_dim < 5 and stokes_dim > 0);
  }

  /** Construct a new Transmission Matrix object
   * 
   * @param[in] tm Matrix to move from
   */
  TransmissionMatrix(TransmissionMatrix&& tm) noexcept
      : stokes_dim(std::move(tm.stokes_dim)),
        T4(std::move(tm.T4)),
        T3(std::move(tm.T3)),
        T2(std::move(tm.T2)),
        T1(std::move(tm.T1)) {}

  /** Construct a new Transmission Matrix object
   * 
   * @param[in] tm matrix to copy
   */
  TransmissionMatrix(const TransmissionMatrix& tm) = default;

  /** Assignment operator
   * 
   * @param[in] tm matrix to copy
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator=(const TransmissionMatrix& tm) = default;

  /** Move operator
   * 
   * @param[in] tm matrix to move from
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator=(TransmissionMatrix&& tm) noexcept {
    stokes_dim = std::move(tm.stokes_dim);
    T4 = std::move(tm.T4);
    T3 = std::move(tm.T3);
    T2 = std::move(tm.T2);
    T1 = std::move(tm.T1);
    return *this;
  }

  /** Construct a new Transmission Matrix from a Propagation Matrix
   * 
   * @param[in] pm Propagation Matrix
   * @param[in] r Distance
   */
  TransmissionMatrix(const PropagationMatrix& pm, const Numeric& r = 1.0);

  operator Tensor3() const {
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

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix4d& Matrix
   */
  const Eigen::Matrix4d& Mat4(size_t i) const { return T4[i]; }

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix3d& Matrix
   */
  const Eigen::Matrix3d& Mat3(size_t i) const { return T3[i]; }

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix2d& Matrix
   */
  const Eigen::Matrix2d& Mat2(size_t i) const { return T2[i]; }

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix<double, 1, 1>& Matrix
   */
  const Eigen::Matrix<double, 1, 1>& Mat1(size_t i) const { return T1[i]; }

  /** Get Matrix at position by copy
   * 
   * @param[in] i Position
   * @return Right size
   */
  Eigen::MatrixXd Mat(size_t i) const {
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
    std::terminate();
  }

  /** Get Matrix at position
   * 
   * @param [in]i Position
   * @return Eigen::Matrix4d& Matrix
   */
  Eigen::Matrix4d& Mat4(size_t i) { return T4[i]; }

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return Eigen::Matrix3d& Matrix
   */
  Eigen::Matrix3d& Mat3(size_t i) { return T3[i]; }

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return Eigen::Matrix42& Matrix
   */
  Eigen::Matrix2d& Mat2(size_t i) { return T2[i]; }

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return Eigen::Matrix<double, 1, 1>& Matrix
   */
  Eigen::Matrix<double, 1, 1>& Mat1(size_t i) { return T1[i]; }

  /** Set to identity matrix */
  void setIdentity() {
    for (auto& T : T4) T = Eigen::Matrix4d::Identity();
    for (auto& T : T3) T = Eigen::Matrix3d::Identity();
    for (auto& T : T2) T = Eigen::Matrix2d::Identity();
    for (auto& T : T1) T(0, 0) = 1;
  }

  /** Set to zero matrix */
  void setZero() {
    for (auto& T : T4) T = Eigen::Matrix4d::Zero();
    for (auto& T : T3) T = Eigen::Matrix3d::Zero();
    for (auto& T : T2) T = Eigen::Matrix2d::Zero();
    for (auto& T : T1) T(0, 0) = 0;
  }

  /** Set this to a multiple of A by B
   * 
   * *this is not aliased with A or B
   * 
   * @param[in] A Matrix 1
   * @param[in] B Matrix 2
   */
  void mul(const TransmissionMatrix& A, const TransmissionMatrix& B) {
    for (size_t i = 0; i < T4.size(); i++) T4[i].noalias() = A.T4[i] * B.T4[i];
    for (size_t i = 0; i < T3.size(); i++) T3[i].noalias() = A.T3[i] * B.T3[i];
    for (size_t i = 0; i < T2.size(); i++) T2[i].noalias() = A.T2[i] * B.T2[i];
    for (size_t i = 0; i < T1.size(); i++) T1[i].noalias() = A.T1[i] * B.T1[i];
  }


  /** Set this to a multiple of A by B
   * 
   * *this is aliased with A or B
   * 
   * @param[in] A Matrix 1
   * @param[in] B Matrix 2
   */
  void mul_aliased(const TransmissionMatrix& A, const TransmissionMatrix& B) {
    for (size_t i = 0; i < T4.size(); i++) T4[i] = A.T4[i] * B.T4[i];
    for (size_t i = 0; i < T3.size(); i++) T3[i] = A.T3[i] * B.T3[i];
    for (size_t i = 0; i < T2.size(); i++) T2[i] = A.T2[i] * B.T2[i];
    for (size_t i = 0; i < T1.size(); i++) T1[i] = A.T1[i] * B.T1[i];
  }

  /** Access value in matrix
   * 
   * @param[in] i Position in vector
   * @param[in] j Row in matrix
   * @param[in] k Col in matrix
   * @return const Numeric& value
   */
  const Numeric& operator()(const Index i, const Index j, const Index k) const {
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

  /** Stokes dimensionaility */
  Index StokesDim() const { return stokes_dim; }

  /** Number of frequencies */
  Index Frequencies() const {
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

  /** Assign to *this lazily
   * 
   * @param[in] lstm Lazy matrix
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator+=(const LazyScale<TransmissionMatrix>& lstm) {
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

  /** Scale self 
   * 
   * @param[in] scale To scale with 
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator*=(const Numeric& scale) {
    for (auto& T : T4) T *= scale;
    for (auto& T : T3) T *= scale;
    for (auto& T : T2) T *= scale;
    for (auto& T : T1) T *= scale;
    return *this;
  }

  /** Assign lazily
   * 
   * @param[in] lstm Lazy value
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator=(const LazyScale<TransmissionMatrix>& lstm) {
    operator=(lstm.bas);
    operator*=(lstm.scale);
    return *this;
  }
  
  std::pair<Numeric, Numeric> weight(Index i) const { 
    auto e_mtau = this -> operator()(i, 0, 0);
    
    if (e_mtau > 0.9 /* tau approx 0.1 */) {
      return {0.5, 0.5};
    } else {
      auto tau = - std::log(e_mtau);
      return {(1 - (1 + tau) * e_mtau) / tau, (tau - 1 + e_mtau) / tau};
    }
  }

  /** Output operator */
  friend std::ostream& operator<<(std::ostream& os,
                                  const TransmissionMatrix& tm);

  /** Input operator */
  friend std::istream& operator>>(std::istream& data, TransmissionMatrix& tm);
};

/** Lazy scale of Transmission Matrix
 * 
 * @param[in] tm Transmission Matrix
 * @param[in] x Scale
 * @return Lazy Transmission Matrix
 */
inline LazyScale<TransmissionMatrix> operator*(const TransmissionMatrix& tm,
                                               const Numeric& x) {
  return LazyScale<TransmissionMatrix>(tm, x);
}

/** Lazy scale of Transmission Matrix
 * 
 * @param[in] x Scale
 * @param[in] tm Transmission Matrix
 * @return Lazy Transmission Matrix
 */
inline LazyScale<TransmissionMatrix> operator*(const Numeric& x,
                                               const TransmissionMatrix& tm) {
  return LazyScale<TransmissionMatrix>(tm, x);
}

/** Radiation Vector for Stokes dimension 1-4 */
class RadiationVector {
 private:
  Index stokes_dim;
  std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d>> R4;
  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> R3;
  std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> R2;
  std::vector<Eigen::Matrix<double, 1, 1>,
              Eigen::aligned_allocator<Eigen::Matrix<double, 1, 1>>>
      R1;

 public:

  /** Construct a new Radiation Vector object
   * 
   * @param[in] nf Number of frequencies
   * @param[in] stokes Stokes dimension
   */
  RadiationVector(Index nf = 0, Index stokes = 1)
      : stokes_dim(stokes),
        R4(stokes_dim == 4 ? nf : 0, Eigen::Vector4d::Zero()),
        R3(stokes_dim == 3 ? nf : 0, Eigen::Vector3d::Zero()),
        R2(stokes_dim == 2 ? nf : 0, Eigen::Vector2d::Zero()),
        R1(stokes_dim == 1 ? nf : 0, Eigen::Matrix<double, 1, 1>::Zero()) {
    assert(stokes_dim < 5 and stokes_dim > 0);
  }

  /** Construct a new Radiation Vector object
   * 
   * @param[in] rv Old Radiation Vector to move from
   */
  RadiationVector(RadiationVector&& rv) noexcept
      : stokes_dim(std::move(rv.stokes_dim)),
        R4(std::move(rv.R4)),
        R3(std::move(rv.R3)),
        R2(std::move(rv.R2)),
        R1(std::move(rv.R1)) {}

  /** Construct a new Radiation Vector object
   * 
   * @param[in] rv Old Vector to copy from
   */
  RadiationVector(const RadiationVector& rv) = default;

  /** Assign old radiation vector to this
   * 
   * @param[in] rv old vetor to copy
   * @return RadiationVector& *this
   */
  RadiationVector& operator=(const RadiationVector& rv) = default;

  /** Assign old radiation vector to this
   * 
   * @param[in] rv old vetor to move from
   * @return RadiationVector& *this
   */
  RadiationVector& operator=(RadiationVector&& rv) noexcept {
    stokes_dim = std::move(rv.stokes_dim);
    R4 = std::move(rv.R4);
    R3 = std::move(rv.R3);
    R2 = std::move(rv.R2);
    R1 = std::move(rv.R1);
    return *this;
  }

  /** Multiply radiation vector from the left
   * 
   * @param[in] T Tranmission Vector
   */
  void leftMul(const TransmissionMatrix& T) {
    for (size_t i = 0; i < R4.size(); i++) R4[i] = T.Mat4(i) * R4[i];
    for (size_t i = 0; i < R3.size(); i++) R3[i] = T.Mat3(i) * R3[i];
    for (size_t i = 0; i < R2.size(); i++) R2[i] = T.Mat2(i) * R2[i];
    for (size_t i = 0; i < R1.size(); i++) R1[i] = T.Mat1(i) * R1[i];
  }

  /** Set Radiation Vector to Zero at position
   * 
   * @param[in] i position
   */
  void SetZero(size_t i) {
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

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Vector4d& Vector
   */
  const Eigen::Vector4d& Vec4(size_t i) const { return R4[i]; }

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Vector3d& Vector
   */
  const Eigen::Vector3d& Vec3(size_t i) const { return R3[i]; }

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Vector2d& Vector
   */
  const Eigen::Vector2d& Vec2(size_t i) const { return R2[i]; }

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Matrix<double, 1, 1>& Vector
   */
  const Eigen::Matrix<double, 1, 1>& Vec1(size_t i) const { return R1[i]; }
  
  /** Return Vector at position by copy
   * 
   * @param[in] i position
   * @return const Eigen::Matrix<double, 1, 1>& Vector
   */
  Eigen::VectorXd Vec(size_t i) const { 
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


  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Vector4d& Vector
   */
  Eigen::Vector4d& Vec4(size_t i) { return R4[i]; }

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Vector3d& Vector
   */
  Eigen::Vector3d& Vec3(size_t i) { return R3[i]; }

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Vector2d& Vector
   */
  Eigen::Vector2d& Vec2(size_t i) { return R2[i]; }

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Matrix<double, 1, 1>& Vector
   */
  Eigen::Matrix<double, 1, 1>& Vec1(size_t i) { return R1[i]; }

  /** Remove the average of two other RadiationVector from *this
   * 
   * @param[in] O1 Input 1
   * @param[in] O2 Input 2
   */
  void rem_avg(const RadiationVector& O1, const RadiationVector& O2) {
    for (size_t i = 0; i < R4.size(); i++)
      R4[i].noalias() -= 0.5 * (O1.R4[i] + O2.R4[i]);
    for (size_t i = 0; i < R3.size(); i++)
      R3[i].noalias() -= 0.5 * (O1.R3[i] + O2.R3[i]);
    for (size_t i = 0; i < R2.size(); i++)
      R2[i].noalias() -= 0.5 * (O1.R2[i] + O2.R2[i]);
    for (size_t i = 0; i < R1.size(); i++)
      R1[i].noalias() -= 0.5 * (O1.R1[i] + O2.R1[i]);
  }


  /** Add the average of two other RadiationVector to *this
   * 
   * @param[in] O1 Input 1
   * @param[in] O2 Input 2
   */
  void add_avg(const RadiationVector& O1, const RadiationVector& O2) {
    for (size_t i = 0; i < R4.size(); i++)
      R4[i].noalias() += 0.5 * (O1.R4[i] + O2.R4[i]);
    for (size_t i = 0; i < R3.size(); i++)
      R3[i].noalias() += 0.5 * (O1.R3[i] + O2.R3[i]);
    for (size_t i = 0; i < R2.size(); i++)
      R2[i].noalias() += 0.5 * (O1.R2[i] + O2.R2[i]);
    for (size_t i = 0; i < R1.size(); i++)
      R1[i].noalias() += 0.5 * (O1.R1[i] + O2.R1[i]);
  }
  
  
  /** Remove two weighted RadiationVector to *this
   * 
   * @param[in] O1 Input 1
   * @param[in] O2 Input 2
   * @param[in] T Transmission matrix from O1 to O2
   */
  void rem_weighted(const RadiationVector& O1, const RadiationVector& O2, const TransmissionMatrix& T) {
    for (size_t i = 0; i < R4.size(); i++) {
      auto w = T.weight(i);
      R4[i].noalias() -= w.first * O1.R4[i] + w.second * O2.R4[i];
    } for (size_t i = 0; i < R3.size(); i++) {
      auto w = T.weight(i);
      R3[i].noalias() -= w.first * O1.R3[i] + w.second * O2.R3[i];
    } for (size_t i = 0; i < R2.size(); i++) {
      auto w = T.weight(i);
      R2[i].noalias() -= w.first * O1.R2[i] + w.second * O2.R2[i];
    } for (size_t i = 0; i < R1.size(); i++) {
      auto w = T.weight(i);
      R1[i].noalias() -= w.first * O1.R1[i] + w.second * O2.R1[i];
    }
  }
  
  
  /** Add two weighted RadiationVector to *this
   * 
   * @param[in] O1 Input 1
   * @param[in] O2 Input 2
   * @param[in] T Transmission matrix from O1 to O2
   */
  void add_weighted(const RadiationVector& O1, const RadiationVector& O2, const TransmissionMatrix& T) {
    for (size_t i = 0; i < R4.size(); i++) {
      auto w = T.weight(i);
      R4[i].noalias() += w.first * O1.R4[i] + w.second * O2.R4[i];
    } for (size_t i = 0; i < R3.size(); i++) {
      auto w = T.weight(i);
      R3[i].noalias() += w.first * O1.R3[i] + w.second * O2.R3[i];
    } for (size_t i = 0; i < R2.size(); i++) {
      auto w = T.weight(i);
      R2[i].noalias() += w.first * O1.R2[i] + w.second * O2.R2[i];
    } for (size_t i = 0; i < R1.size(); i++) {
      auto w = T.weight(i);
      R1[i].noalias() += w.first * O1.R1[i] + w.second * O2.R1[i];
    }
  }

  /** Add the emission derivative to this
   * 
   * @param[in] PiT Accumulated transmission to space
   * @param[in] dT Derivative of transmission matrix
   * @param[in] T Transmission matrix
   * @param[in] ImJ Intensity minus the emission vector
   * @param[in] dJ Derivative of the emission vector
   */
  void addDerivEmission(const TransmissionMatrix& PiT,
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


  /** Add the transmission derivative to this
   * 
   * @param[in] PiT Accumulated transmission to space
   * @param[in] dT Derivative of transmission matrix
   * @param[in] I Intensity vector
   */
  void addDerivTransmission(const TransmissionMatrix& PiT,
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


  /** Add multiply
   * 
   * Performs essentially this += A * x
   * 
   * Assumes this and x are not aliases
   * 
   * @param[in] A Derivative of transmission matrix
   * @param[in] x Intensity vector
   */
  void addMultiplied(const TransmissionMatrix& A,
                     const RadiationVector& x) {
    for (size_t i = 0; i < R4.size(); i++)
      R4[i].noalias() += A.Mat4(i) * x.R4[i];
    for (size_t i = 0; i < R3.size(); i++)
      R3[i].noalias() += A.Mat3(i) * x.R3[i];
    for (size_t i = 0; i < R2.size(); i++)
      R2[i].noalias() += A.Mat2(i) * x.R2[i];
    for (size_t i = 0; i < R1.size(); i++)
      R1[i].noalias() += A.Mat1(i) * x.R1[i];
  }

  /** Sets *this to the reflection derivative
   * 
   * @param[in] I Intensity vector
   * @param[in] PiT Accumulated transmission to space
   * @param[in] Z Reflection matrix
   * @param[in] dZ Derivative of reflection matrix
   */
  void setDerivReflection(const RadiationVector& I,
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
      R4[i][0] = PiT.Mat1(i)[0] *
                 (Z.Mat1(i)[0] * R1[i][0] + dZ.Mat1(i)[0] * I.R1[i][0]);
  }

  /** Set this to backscatter transmission
   * 
   * @param[in] I0 Incoming intensity vector
   * @param[in] Tr Reflection transmission
   * @param[in] Tf Forward transmission
   * @param[in] Z Reflection matrix
   */
  void setBackscatterTransmission(const RadiationVector& I0,
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

  /** Set this to backscatter transmission scatter derivative
   * 
   * @param[in] I0 Incoming intensity vector
   * @param[in] Tr Reflection transmission
   * @param[in] Tf Forward transmission
   * @param[in] dZ Reflection matrix derivative
   */
  void setBackscatterTransmissionDerivative(const RadiationVector& I0,
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

  /** Set *this from matrix
   * 
   * @param[in] M Matrix
   * @return RadiationVector& *this
   */
  RadiationVector& operator=(const ConstMatrixView& M) {
    assert(M.ncols() == stokes_dim and M.nrows() == Frequencies());
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
    return *this;
  }

  /** Access operator
   * 
   * @param[in] i Position in outer vector
   * @param[in] j Position in inner vector
   * @return const Numeric& 
   */
  const Numeric& operator()(const Index i, const Index j) const {
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

  /** Convert *this to Matrix class
   * 
   * @return Matrix
   */
  operator Matrix() const {
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

  /** Set this to source vector at position
   * 
   * @param[in] a Absorption vector
   * @param[in] B Planck vector
   * @param[in] S Scattering source vector
   * @param[in] i Position
   */
  void setSource(const StokesVector& a,
                 const ConstVectorView& B,
                 const StokesVector& S,
                 Index i) {
    assert(a.NumberOfAzimuthAngles() == 1);
    assert(a.NumberOfZenithAngles() == 1);
    assert(S.NumberOfAzimuthAngles() == 1);
    assert(S.NumberOfZenithAngles() == 1);
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

  /** Get Stokes dimension */
  Index StokesDim() const { return stokes_dim; }

  /** Get frequency count */
  Index Frequencies() const {
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

  /** Output operator */
  friend std::ostream& operator<<(std::ostream& os, const RadiationVector& rv);

  /** Input operator */
  friend std::istream& operator>>(std::istream& data, RadiationVector& rv);
};

typedef Array<TransmissionMatrix> ArrayOfTransmissionMatrix;
typedef Array<ArrayOfTransmissionMatrix> ArrayOfArrayOfTransmissionMatrix;
typedef Array<ArrayOfArrayOfTransmissionMatrix> ArrayOfArrayOfArrayOfTransmissionMatrix;
typedef Array<RadiationVector> ArrayOfRadiationVector;
typedef Array<ArrayOfRadiationVector> ArrayOfArrayOfRadiationVector;
typedef Array<ArrayOfArrayOfRadiationVector> ArrayOfArrayOfArrayOfRadiationVector;

  /** Output operator */
std::ostream& operator<<(std::ostream& os, const TransmissionMatrix& tm);

  /** Output operator */
std::ostream& operator<<(std::ostream& os,
                         const ArrayOfTransmissionMatrix& atm);

  /** Output operator */
std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfTransmissionMatrix& aatm);

  /** Output operator */
std::ostream& operator<<(std::ostream& os, const RadiationVector& rv);

  /** Output operator */
std::ostream& operator<<(std::ostream& os, const ArrayOfRadiationVector& arv);

  /** Output operator */
std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfRadiationVector& aarv);

  /** Input operator */
std::istream& operator>>(std::istream& is, TransmissionMatrix& tm);

  /** Input operator */
std::istream& operator>>(std::istream& is, RadiationVector& rv);

/** Intended to hold various backscatter solvers */
enum class BackscatterSolver {
  CommutativeTransmission,
  FullTransmission,
};

/** Intended to hold various ways to accumulate the transmission matrix */
enum class CumulativeTransmission {
  Forward,
  Reverse,
};

/** Intended to hold various forward solvers */
enum class RadiativeTransferSolver {
  Emission,
  Transmission,
  WeightedEmission,
};

/** Update the Radiation Vector
 * 
 * @param[in,out] I Radiation vector
 * @param[in,out] dI1 Radiation vector derivatives to level 1
 * @param[in,out] dI2 Radiation vector derivatives to level 2
 * @param[in] J1 Source vector from level 1
 * @param[in] J2 Source vector from level 2
 * @param[in] dJ1 Source vector derivative from level 1
 * @param[in] dJ2 Source vector derivative from level 2
 * @param[in] T Transmission matrix through layer
 * @param[in] PiT Accumulated transmission matrix to space
 * @param[in] dT1 Transmission matrix derivatives through layer from level 1
 * @param[in] dT2 Transmission matrix derivatives through layer from level 2
 * @param[in] solver Type of solver to use
 */
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

/** Set the stepwise source
 * 
 * @param[in,out] J Source vector
 * @param[in,out] dJ Source vector derivatives
 * @param[in] K Propagation matrix
 * @param[in] a Absorption vector
 * @param[in] S Scattering source vector
 * @param[in] dK_dx Propagation matrix derivatives
 * @param[in] da_dx Absorption vector derivatives
 * @param[in] dS_dx Scattering source vector derivatives
 * @param[in] B Planck vector
 * @param[in] dB_dT Planck vector derivative wrt temperature
 * @param[in] jacobian_quantities As WSV
 * @param[in] jacobian_do Do Jacobian?
 */
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

/** Set the stepwise transmission matrix
 * 
 * @param[in,out] T Transmission matrix
 * @param[in,out] dT1 Transmission matrix derivative wrt level 1
 * @param[in,out] dT2 Transmission matrix derivative wrt level 2
 * @param[in] K1 Propagation matrix wrt level 1
 * @param[in] K2 Propagation matrix wrt level 2
 * @param[in] dK1_dx Propagation matrix derivative wrt level 1
 * @param[in] dK2_dx Propagation matrix derivative wrt level 2
 * @param[in] r Distance through layer
 * @param[in] dr_dtemp1 Distance through layer derivative wrt temperature of level 1
 * @param[in] dr_dtemp2 Distance through layer derivative wrt temperature of level 2
 * @param[in] temp_deriv_pos Position of derivative of temperature (-1 if not present)
 */
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

/** Accumulate the transmission matrix over all layers
 * 
 * @param[in] T Transmission matrix through all layers
 * @param[in] type Type of accumulation to target
 * @return ArrayOfTransmissionMatrix Transmission to target
 */
ArrayOfTransmissionMatrix cumulative_transmission(
    const ArrayOfTransmissionMatrix& T,
    const CumulativeTransmission type) /*[[expects: T.nelem()>0]]*/;

/** Set the backscatter radiation vector
 * 
 * @param[in,out] I Radiation vector of all layers
 * @param[in,out] dI Radiation vector derivative of all layers
 * @param[in] I_incoming Incoming radiation vector
 * @param[in] T Transmission matrix of all layers
 * @param[in] PiTf Forwards accumulated transmission of all layers
 * @param[in] PiTr Backwards accumulated transmission of all layers
 * @param[in] Z  Reflection matrix of all layers
 * @param[in] dT1 Transmission matrix derivative for level 1 of all layers
 * @param[in] dT2 Transmission matrix derivative for level 2 of all layers
 * @param[in] dZ  erivative of reflection matrix of all layers
 * @param[in] solver Type of backscattering of all layers
 */
void set_backscatter_radiation_vector(
    ArrayOfRadiationVector& I,
    ArrayOfArrayOfArrayOfRadiationVector& dI,
    const RadiationVector & I_incoming,
    const ArrayOfTransmissionMatrix& T,
    const ArrayOfTransmissionMatrix& PiTf,
    const ArrayOfTransmissionMatrix& PiTr,
    const ArrayOfTransmissionMatrix& Z,
    const ArrayOfArrayOfTransmissionMatrix& dT1,
    const ArrayOfArrayOfTransmissionMatrix& dT2,
    const ArrayOfArrayOfTransmissionMatrix& dZ,
    const BackscatterSolver solver);

/** Accumulated backscatter (???)
 * 
 * FIXMEDOC Patrick, these are translated from other functions that accumulate
 * the back-scattering.  I do not remember how this worked.
 * 
 * @param t Tensor5 of backscattering
 * @param m (???)
 * @return ArrayOfTransmissionMatrix cumulative backscattering
 */
ArrayOfTransmissionMatrix cumulative_backscatter(ConstTensor5View t,
                                                 ConstMatrixView m);

/** Accumulated backscatter derivative (???)
 * 
 * FIXMEDOC Patrick, these are translated from other functions that accumulate
 * the back-scattering.  I do not remember how this worked.
 * 
 * @param t Tensor5 of backscattering
 * @param m (???)
 * @return ArrayOfTArrayOfTransmissionMatrix cumulative backscattering
 */
ArrayOfArrayOfTransmissionMatrix cumulative_backscatter_derivative(
    ConstTensor5View t, const ArrayOfMatrix& aom);

#endif  // transmissionmatrix_h
