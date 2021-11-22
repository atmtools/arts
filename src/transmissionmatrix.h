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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#if !defined(__clang__)
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#include <Eigen/Dense>
#include <Eigen/StdVector>

#pragma GCC diagnostic pop

#include "jacobian.h"
#include "propagationmatrix.h"

/** Class to keep track of Transmission Matrices for Stokes Dim 1-4 */
struct TransmissionMatrix {
  Index stokes_dim;
  std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>> T4;
  std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>> T3;
  std::vector<Eigen::Matrix2d, Eigen::aligned_allocator<Eigen::Matrix2d>> T2;
  std::vector<Eigen::Matrix<double, 1, 1>,
              Eigen::aligned_allocator<Eigen::Matrix<double, 1, 1>>>
      T1;

  /** Construct a new Transmission Matrix object
   * 
   * @param[in] nf Number of frequencies
   * @param[in] stokes Stokes dimension
   */
  TransmissionMatrix(Index nf = 0, Index stokes = 1);

  /** Construct a new Transmission Matrix object
   * 
   * @param[in] tm Matrix to move from
   */
  TransmissionMatrix(TransmissionMatrix&& tm) noexcept = default;

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
  TransmissionMatrix& operator=(TransmissionMatrix&& tm) noexcept = default;

  /** Construct a new Transmission Matrix from a Propagation Matrix
   * 
   * @param[in] pm Propagation Matrix
   * @param[in] r Distance
   */
  explicit TransmissionMatrix(const PropagationMatrix& pm,
                              const Numeric& r = 1.0);

  operator Tensor3() const;

  explicit TransmissionMatrix(const ConstMatrixView& mat): TransmissionMatrix(1, mat.nrows()){
   ARTS_ASSERT(mat.nrows() == mat.ncols());
   for (Index i = 0; i < stokes_dim; i++)
    for (Index j = 0; j < stokes_dim; j++)
      operator()(0, i, j) = mat(i,j);
  };

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix4d& Matrix
   */
  [[nodiscard]] const Eigen::Matrix4d& Mat4(size_t i) const;

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix3d& Matrix
   */
  [[nodiscard]] const Eigen::Matrix3d& Mat3(size_t i) const;

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix2d& Matrix
   */
  [[nodiscard]] const Eigen::Matrix2d& Mat2(size_t i) const;

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return const Eigen::Matrix<double, 1, 1>& Matrix
   */
  [[nodiscard]] const Eigen::Matrix<double, 1, 1>& Mat1(size_t i) const;

  /** Get Matrix at position by copy
   * 
   * @param[in] i Position
   * @return Right size
   */
  [[nodiscard]] Eigen::MatrixXd Mat(size_t i) const;

  /** Get Matrix at position
   * 
   * @param [in]i Position
   * @return Eigen::Matrix4d& Matrix
   */
  Eigen::Matrix4d& Mat4(size_t i);

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return Eigen::Matrix3d& Matrix
   */
  Eigen::Matrix3d& Mat3(size_t i);

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return Eigen::Matrix42& Matrix
   */
  Eigen::Matrix2d& Mat2(size_t i);

  /** Get Matrix at position
   * 
   * @param[in] i Position
   * @return Eigen::Matrix<double, 1, 1>& Matrix
   */
  Eigen::Matrix<double, 1, 1>& Mat1(size_t i);

  /** Set to identity matrix */
  void setIdentity();

  /** Set to zero matrix */
  void setZero();

  /** Set this to a multiple of A by B
   * 
   * *this is not aliased with A or B
   * 
   * @param[in] A Matrix 1
   * @param[in] B Matrix 2
   */
  void mul(const TransmissionMatrix& A, const TransmissionMatrix& B);

  /** Set this to a multiple of A by B
   * 
   * *this is aliased with A or B
   * 
   * @param[in] A Matrix 1
   * @param[in] B Matrix 2
   */
  void mul_aliased(const TransmissionMatrix& A, const TransmissionMatrix& B);

  /** Access value in matrix
   * 
   * @param[in] i Position in vector
   * @param[in] j Row in matrix
   * @param[in] k Col in matrix
   * @return const Numeric& value
   */
  [[nodiscard]] Numeric operator()(const Index i,
                                   const Index j,
                                   const Index k) const;

  /** Access value in matrix
   *
   * @param[in] i Position in vector
   * @param[in] j Row in matrix
   * @param[in] k Col in matrix
   * @return const Numeric& value
   */
  Numeric& operator()(const Index i, const Index j, const Index k) {
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



  /** Number of frequencies */
  [[nodiscard]] Index Frequencies() const;

  /** Assign to *this lazily
   * 
   * @param[in] lstm Lazy matrix
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator+=(const LazyScale<TransmissionMatrix>& lstm);

  /** Scale self 
   * 
   * @param[in] scale To scale with 
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator*=(const Numeric& scale);

  /** Assign lazily
   * 
   * @param[in] lstm Lazy value
   * @return TransmissionMatrix& *this
   */
  TransmissionMatrix& operator=(const LazyScale<TransmissionMatrix>& lstm);

  /** Output operator */
  friend std::ostream& operator<<(std::ostream& os,
                                  const TransmissionMatrix& tm);

  /** Input operator */
  friend std::istream& operator>>(std::istream& data, TransmissionMatrix& tm);

  /** Simple template access for the transmission */
  template <int N>
  auto& TraMat(size_t i) noexcept {
    static_assert(N > 0 and N < 5, "Bad size N");
    if constexpr (N == 1)
      return T1[i];
    else if constexpr (N == 2)
      return T2[i];
    else if constexpr (N == 3)
      return T3[i];
    else if constexpr (N == 4)
      return T4[i];
  }

  /** Simple template access for the transmission */
  template <int N>
  [[nodiscard]] auto& TraMat(size_t i) const noexcept {
    static_assert(N > 0 and N < 5, "Bad size N");
    if constexpr (N == 1)
      return T1[i];
    else if constexpr (N == 2)
      return T2[i];
    else if constexpr (N == 3)
      return T3[i];
    else if constexpr (N == 4)
      return T4[i];
  }

  /** Simple template access for the optical depth */
  template <int N>
  [[nodiscard]] Eigen::Matrix<Numeric, N, N> OptDepth(size_t i) const noexcept {
    return (-TraMat<N>(i).diagonal().array().log().matrix()).asDiagonal();
  }

  /*! Return the weighted source term using second order integration
   * 
   \f[ far = \frac{1-\left(1+\log{T_{00}}\right) T}{\log{T_{00}}} \f]
   \f[ close = \frac{\log{T_{00}} - 1 + T}{\log{T_{00}}} \f]
   * 
   * This follows definition of equation 3.34 of http://www.ita.uni-heidelberg.de/~dullemond/lectures/radtrans_2013/Chapter_3.pdf.
   
   One key change is that we consider polarization but only based on unpolarized radiation
   *
   * FIXME: This function is not done properly for Stokes Dim > 1.  The results might be correct but the derivation is not understood.
   *
   * @param[in] T The transmission matrix
   * @param[in] far The source at the destination of the RT step
   * @param[in] close The source at the start of the RT step
   * @param[in] Kfar The propagation matrix at the destination of the RT step
   * @param[in] Kclose The propagation matrix at the start of the RT step
   * @param[in] r The distance of the RT step
   * @return Linear Weights
   */
  template <int N>
  [[nodiscard]] Eigen::Matrix<Numeric, N, 1> second_order_integration_source(
      const Eigen::Matrix<Numeric, N, N> T,
      const Eigen::Matrix<Numeric, N, 1> far,
      const Eigen::Matrix<Numeric, N, 1> close,
      const Eigen::Matrix<Numeric, N, N> Kfar,
      const Eigen::Matrix<Numeric, N, N> Kclose,
      const Numeric r) const noexcept {
    static_assert(N > 0 and N < 5, "Bad size N");

    const auto I = Eigen::Matrix<Numeric, N, N>::Identity();
    if (T(0, 0) < 0.99) {
      const auto od = 0.5 * r * (Kfar + Kclose);
      Eigen::Matrix<Numeric, N, 1> second =
          od.inverse() * ((I - (I + od) * T) * far + (od - I + T) * close);
      if ((far[0] < close[0] and Kfar(0, 0) > Kclose(0, 0)) or
          (far[0] > close[0] and Kfar(0, 0) < Kclose(0, 0))) {
        Eigen::Matrix<Numeric, N, 1> second_limit =
            0.5 * (I - T) * (far + close);

        // FIXME: This is the equation given in the source material...
        // const Eigen::Matrix<Numeric, N, 1> second_limit = 0.25 * r * (Kfar * far + Kclose * close);

        if (second_limit[0] > second[0]) return second;
        return second_limit;
      }

      return second;
    }

    return 0.5 * (I - T) * (far + close);
  }

  template <int N>
  [[nodiscard]] Eigen::Matrix<Numeric, N, 1> second_order_integration_dsource(
      size_t i,
      const TransmissionMatrix& dx,
      const Eigen::Matrix<Numeric, N, 1> far,
      const Eigen::Matrix<Numeric, N, 1> close,
      const Eigen::Matrix<Numeric, N, 1> d,
      bool isfar) const noexcept {
    static_assert(N > 0 and N < 5, "Bad size N");

    const auto I = Eigen::Matrix<Numeric, N, N>::Identity();
    const auto T = TraMat<N>(i);
    const auto& dTdx = dx.TraMat<N>(i);
    const Eigen::Matrix<Numeric, N, 1> first = 0.5 * (I - T) * (far + close);
    if (T(0, 0) < 0.99) {
      const auto od = OptDepth<N>(i);
      const auto doddx = dx.OptDepth<N>(i);
      const Eigen::Matrix<Numeric, N, 1> second =
          od.inverse() * ((I - (I + od) * T) * far + (od - I + T) * close);
      if (first[0] > second[0]) {
        if (isfar) {
          return od.inverse() *
                 ((I - (I + od) * T) * d - (I + od) * dTdx * far +
                  (doddx + T) * close - doddx * second);
        }

        return od.inverse() * (-(I + od) * dTdx * far + (od - I + T) * d +
                               (doddx + T) * close - doddx * second);
      }

      return 0.5 * ((I - T) * d - dTdx * (far + close));
    }

    return 0.5 * ((I - T) * d - dTdx * (far + close));
  }
};

/** Lazy scale of Transmission Matrix
 * 
 * @param[in] tm Transmission Matrix
 * @param[in] x Scale
 * @return Lazy Transmission Matrix
 */
[[nodiscard]] LazyScale<TransmissionMatrix> operator*(
    const TransmissionMatrix& tm, const Numeric& x);

/** Lazy scale of Transmission Matrix
 * 
 * @param[in] x Scale
 * @param[in] tm Transmission Matrix
 * @return Lazy Transmission Matrix
 */
[[nodiscard]] LazyScale<TransmissionMatrix> operator*(
    const Numeric& x, const TransmissionMatrix& tm);

/** Radiation Vector for Stokes dimension 1-4 */
struct RadiationVector {
  Index stokes_dim;
  std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d>> R4;
  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> R3;
  std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d>> R2;
  std::vector<Eigen::Matrix<double, 1, 1>,
              Eigen::aligned_allocator<Eigen::Matrix<double, 1, 1>>>
      R1;

  /** Construct a new Radiation Vector object
   * 
   * @param[in] nf Number of frequencies
   * @param[in] stokes Stokes dimension
   */
  RadiationVector(Index nf = 0, Index stokes = 1);

  /** Construct a new Radiation Vector object
   * 
   * @param[in] rv Old Radiation Vector to move from
   */
  RadiationVector(RadiationVector&& rv) noexcept = default;

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
  RadiationVector& operator=(RadiationVector&& rv) noexcept = default;

  /** Multiply radiation vector from the left
   * 
   * @param[in] T Tranmission Vector
   */
  void leftMul(const TransmissionMatrix& T);

  /** Set Radiation Vector to Zero at position
   * 
   * @param[in] i position
   */
  void SetZero(size_t i);

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Vector4d& Vector
   */
  [[nodiscard]] const Eigen::Vector4d& Vec4(size_t i) const;

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Vector3d& Vector
   */
  [[nodiscard]] const Eigen::Vector3d& Vec3(size_t i) const;

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Vector2d& Vector
   */
  [[nodiscard]] const Eigen::Vector2d& Vec2(size_t i) const;

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return const Eigen::Matrix<double, 1, 1>& Vector
   */
  [[nodiscard]] const Eigen::Matrix<double, 1, 1>& Vec1(size_t i) const;

  /** Return Vector at position by copy
   * 
   * @param[in] i position
   * @return const Eigen::Matrix<double, 1, 1>& Vector
   */
  [[nodiscard]] Eigen::VectorXd Vec(size_t i) const;

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Vector4d& Vector
   */
  Eigen::Vector4d& Vec4(size_t i);

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Vector3d& Vector
   */
  Eigen::Vector3d& Vec3(size_t i);

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Vector2d& Vector
   */
  Eigen::Vector2d& Vec2(size_t i);

  /** Return Vector at position
   * 
   * @param[in] i position
   * @return Eigen::Matrix<double, 1, 1>& Vector
   */
  Eigen::Matrix<double, 1, 1>& Vec1(size_t i);

  /** Remove the average of two other RadiationVector from *this
   * 
   * @param[in] O1 Input 1
   * @param[in] O2 Input 2
   */
  void rem_avg(const RadiationVector& O1, const RadiationVector& O2);

  /** Add the average of two other RadiationVector to *this
   * 
   * @param[in] O1 Input 1
   * @param[in] O2 Input 2
   */
  void add_avg(const RadiationVector& O1, const RadiationVector& O2);

  /** Add the weighted source of two RadiationVector to *this
   * 
   * @param[in] T The transmission matrix
   * @param[in] far   Input 1
   * @param[in] close Input 2
   */
  void add_weighted(const TransmissionMatrix& T,
                    const RadiationVector& far,
                    const RadiationVector& close,
                    const ConstMatrixView& Kfar,
                    const ConstMatrixView& Kclose,
                    const Numeric r);

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
                        const RadiationVector& dJ);

  /** Add the emission derivative to this
   * 
   * @param[in] PiT Accumulated transmission to space
   * @param[in] dT Derivative of transmission matrix
   * @param[in] T Transmission matrix
   * @param[in] ImJ Intensity minus the emission vector
   * @param[in] dJ Derivative of the emission vector
   */
  void addWeightedDerivEmission(const TransmissionMatrix& PiT,
                                const TransmissionMatrix& dT,
                                const TransmissionMatrix& T,
                                const RadiationVector& I,
                                const RadiationVector& far,
                                const RadiationVector& close,
                                const RadiationVector& d,
                                bool isfar);

  /** Add the transmission derivative to this
   * 
   * @param[in] PiT Accumulated transmission to space
   * @param[in] dT Derivative of transmission matrix
   * @param[in] I Intensity vector
   */
  void addDerivTransmission(const TransmissionMatrix& PiT,
                            const TransmissionMatrix& dT,
                            const RadiationVector& I);

  /** Add multiply
   * 
   * Performs essentially this += A * x
   * 
   * Assumes this and x are not aliases
   * 
   * @param[in] A Derivative of transmission matrix
   * @param[in] x Intensity vector
   */
  void addMultiplied(const TransmissionMatrix& A, const RadiationVector& x);

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
                          const TransmissionMatrix& dZ);

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
                                  const TransmissionMatrix& Z);

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
                                            const TransmissionMatrix& dZ);

  /** Set *this from matrix
   * 
   * @param[in] M Matrix
   * @return RadiationVector& *this
   */
  RadiationVector& operator=(const ConstMatrixView& M);

  /** Access operator
   * 
   * @param[in] i Position in outer vector
   * @param[in] j Position in inner vector
   * @return const Numeric& 
   */
  const Numeric& operator()(const Index i, const Index j) const;

  /** Convert *this to Matrix class
   * 
   * @return Matrix
   */
  operator Matrix() const;

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
                 Index i);

  /** Get frequency count */
  [[nodiscard]] Index Frequencies() const;

  /** Output operator */
  friend std::ostream& operator<<(std::ostream& os, const RadiationVector& rv);

  /** Input operator */
  friend std::istream& operator>>(std::istream& data, RadiationVector& rv);
};

using ArrayOfTransmissionMatrix = Array<TransmissionMatrix>;
using ArrayOfArrayOfTransmissionMatrix = Array<ArrayOfTransmissionMatrix>;
using ArrayOfArrayOfArrayOfTransmissionMatrix =
    Array<ArrayOfArrayOfTransmissionMatrix>;
using ArrayOfRadiationVector = Array<RadiationVector>;
using ArrayOfArrayOfRadiationVector = Array<ArrayOfRadiationVector>;
using ArrayOfArrayOfArrayOfRadiationVector =
    Array<ArrayOfArrayOfRadiationVector>;

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
  LinearWeightedEmission,
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
                             const PropagationMatrix& K1,
                             const PropagationMatrix& K2,
                             const ArrayOfPropagationMatrix& dK1,
                             const ArrayOfPropagationMatrix& dK2,
                             const Numeric r,
                             const Vector& dr1,
                             const Vector& dr2,
                             const Index ia,
                             const Index iz,
                             const RadiativeTransferSolver solver);

/** Set the stepwise source
 * 
 * @param[in,out] J Source vector
 * @param[in,out] dJ Source vector derivatives
 * @param[in] K Propagation matrix
 * @param[in] a Absorption vector
 * @param[in] S Scattering source vector
 * @param[in] dK Propagation matrix derivatives
 * @param[in] da Absorption vector derivatives
 * @param[in] dS Scattering source vector derivatives
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
                     const ArrayOfPropagationMatrix& dK,
                     const ArrayOfStokesVector& da,
                     const ArrayOfStokesVector& dS,
                     const ConstVectorView& B,
                     const ConstVectorView& dB_dT,
                     const ArrayOfRetrievalQuantity& jacobian_quantities,
                     const bool& jacobian_do);

/** Set the stepwise transmission matrix
 * 
 * @param[in,out] T Transmission matrix
 * @param[in,out] dT1 Transmission matrix derivative wrt level 1
 * @param[in,out] dT2 Transmission matrix derivative wrt level 2
 * @param[in] K1 Propagation matrix wrt level 1
 * @param[in] K2 Propagation matrix wrt level 2
 * @param[in] dK1 Propagation matrix derivative wrt level 1
 * @param[in] dK2 Propagation matrix derivative wrt level 2
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
                           const ArrayOfPropagationMatrix& dK1,
                           const ArrayOfPropagationMatrix& dK2,
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
    const RadiationVector& I_incoming,
    const ArrayOfTransmissionMatrix& T,
    const ArrayOfTransmissionMatrix& PiTf,
    const ArrayOfTransmissionMatrix& PiTr,
    const ArrayOfTransmissionMatrix& Z,
    const ArrayOfArrayOfTransmissionMatrix& dT1,
    const ArrayOfArrayOfTransmissionMatrix& dT2,
    const ArrayOfArrayOfTransmissionMatrix& dZ,
    const BackscatterSolver solver);

/** Bulk back-scattering 
 * 
 * Sums up the back-scattering per element with particle number densities.
 *
 * Below ns is Stokes dim, nf the number of frequencies, np number of
 * atmospheric positions, and ne the number of scattering elements.
 * 
 * @param Pe Back-scattering on scattering element basis (ne,np,nf,ns,ns)
 * @param pnd Particle number densities (ne,np)
 * @return ArrayOfTransmissionMatrix Bulk back-scattering matrices
 */
ArrayOfTransmissionMatrix bulk_backscatter(const ConstTensor5View& Pe,
                                           const ConstMatrixView& pnd);

/** Derivatives of bulk back-scattering  
 * 
 * Below ns is Stokes dim, nf the number of frequencies, np number of
 * atmospheric positions, and ne the number of scattering elements.
 * 
 * @param Pe Back-scattering on scattering element basis (ne,np,nf,ns,ns)
 * @param dpnd_dx Derivatives of pnd with respect to Jacobian quantities 
 * @return ArrayOfTArrayOfTransmissionMatrix Derivatives of bulk back-scattering
 */
ArrayOfArrayOfTransmissionMatrix bulk_backscatter_derivative(
    const ConstTensor5View& Pe, const ArrayOfMatrix& dpnd_dx);

#endif  // transmissionmatrix_h
