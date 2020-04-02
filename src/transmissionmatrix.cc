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
#include "complex.h"
#include "constants.h"

constexpr Numeric lower_is_considered_zero_for_sinc_likes = 1e-4;

inline Numeric vector1(const StokesVector& a,
                       const ConstVectorView& B,
                       const StokesVector& da,
                       const ConstVectorView& dB_dT,
                       const StokesVector& dS,
                       bool dT,
                       Index i) noexcept {
  if (dT)
    return dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i];
  else
    return dS.Kjj()[i] + da.Kjj()[i] * B[i];
}

inline Eigen::Vector2d vector2(const StokesVector& a,
                               const ConstVectorView& B,
                               const StokesVector& da,
                               const ConstVectorView& dB_dT,
                               const StokesVector& dS,
                               bool dT,
                               Index i) noexcept {
  if (dT)
    return Eigen::Vector2d(
        dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
        dS.K12()[i] + da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i]);
  else
    return Eigen::Vector2d(dS.Kjj()[i] + da.Kjj()[i] * B[i],
                           dS.K12()[i] + da.K12()[i] * B[i]);
}

inline Eigen::Vector3d vector3(const StokesVector& a,
                               const ConstVectorView& B,
                               const StokesVector& da,
                               const ConstVectorView& dB_dT,
                               const StokesVector& dS,
                               bool dT,
                               Index i) noexcept {
  if (dT)
    return Eigen::Vector3d(
        dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
        dS.K12()[i] + da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i],
        dS.K13()[i] + da.K13()[i] * B[i] + a.K13()[i] * dB_dT[i]);
  else
    return Eigen::Vector3d(dS.Kjj()[i] + da.Kjj()[i] * B[i],
                           dS.K12()[i] + da.K12()[i] * B[i],
                           dS.K13()[i] + da.K13()[i] * B[i]);
}

inline Eigen::Vector4d vector4(const StokesVector& a,
                               const ConstVectorView& B,
                               const StokesVector& da,
                               const ConstVectorView& dB_dT,
                               const StokesVector& dS,
                               bool dT,
                               Index i) noexcept {
  if (dT)
    return Eigen::Vector4d(
        dS.Kjj()[i] + da.Kjj()[i] * B[i] + a.Kjj()[i] * dB_dT[i],
        dS.K12()[i] + da.K12()[i] * B[i] + a.K12()[i] * dB_dT[i],
        dS.K13()[i] + da.K13()[i] * B[i] + a.K13()[i] * dB_dT[i],
        dS.K14()[i] + da.K14()[i] * B[i] + a.K14()[i] * dB_dT[i]);
  else
    return Eigen::Vector4d(dS.Kjj()[i] + da.Kjj()[i] * B[i],
                           dS.K12()[i] + da.K12()[i] * B[i],
                           dS.K13()[i] + da.K13()[i] * B[i],
                           dS.K14()[i] + da.K14()[i] * B[i]);
}

inline Eigen::Matrix<double, 1, 1> matrix1(const Numeric& a) noexcept {
  return Eigen::Matrix<double, 1, 1>(a);
}

inline Eigen::Matrix2d matrix2(const Numeric& a, const Numeric& b) noexcept {
  return (Eigen::Matrix2d() << a, b, b, a).finished();
}

inline Eigen::Matrix3d matrix3(const Numeric& a,
                               const Numeric& b,
                               const Numeric& c,
                               const Numeric& u) noexcept {
  return (Eigen::Matrix3d() << a, b, c, b, a, u, c, -u, a).finished();
}

inline Eigen::Matrix4d matrix4(const Numeric& a,
                               const Numeric& b,
                               const Numeric& c,
                               const Numeric& d,
                               const Numeric& u,
                               const Numeric& v,
                               const Numeric& w) noexcept {
  return (Eigen::Matrix4d() << a,
          b,
          c,
          d,
          b,
          a,
          u,
          v,
          c,
          -u,
          a,
          w,
          d,
          -v,
          -w,
          a)
      .finished();
}

inline Eigen::Matrix<double, 1, 1> matrix1(const ConstMatrixView m) noexcept {
  return Eigen::Matrix<double, 1, 1>(m(0, 0));
}

inline Eigen::Matrix2d matrix2(const ConstMatrixView m) noexcept {
  return (Eigen::Matrix2d() << m(0, 0), m(0, 1), m(1, 0), m(1, 1)).finished();
}

inline Eigen::Matrix3d matrix3(const ConstMatrixView m) noexcept {
  return (Eigen::Matrix3d() << m(0, 0),
          m(0, 1),
          m(0, 2),
          m(1, 0),
          m(1, 1),
          m(1, 2),
          m(2, 0),
          m(2, 1),
          m(2, 2))
      .finished();
}

inline Eigen::Matrix4d matrix4(const ConstMatrixView m) noexcept {
  return (Eigen::Matrix4d() << m(0, 0),
          m(0, 1),
          m(0, 2),
          m(0, 3),
          m(1, 0),
          m(1, 1),
          m(1, 2),
          m(1, 3),
          m(2, 0),
          m(2, 1),
          m(2, 2),
          m(2, 3),
          m(3, 0),
          m(3, 1),
          m(3, 2),
          m(3, 3))
      .finished();
}

inline Eigen::Matrix<double, 1, 1> inv1(const Numeric& a) noexcept {
  return Eigen::Matrix<double, 1, 1>(1 / a);
}

inline Eigen::Matrix2d inv2(const Numeric& a, const Numeric& b) noexcept {
  return (Eigen::Matrix2d() << a, -b, -b, a).finished() / (a * a - b * b);
}

inline Eigen::Matrix3d inv3(const Numeric& a,
                            const Numeric& b,
                            const Numeric& c,
                            const Numeric& u) noexcept {
  return (Eigen::Matrix3d() << a * a + u * u,
          -a * b - c * u,
          -a * c + b * u,
          -a * b + c * u,
          a * a - c * c,
          -a * u + b * c,
          -a * c - b * u,
          a * u + b * c,
          a * a - b * b)
             .finished() /
         (a * a * a - a * b * b - a * c * c + a * u * u);
}

inline Eigen::Matrix4d inv4(const Numeric& a,
                            const Numeric& b,
                            const Numeric& c,
                            const Numeric& d,
                            const Numeric& u,
                            const Numeric& v,
                            const Numeric& w) noexcept {
  return (Eigen::Matrix4d() << a * a * a + a * u * u + a * v * v + a * w * w,
          -a * a * b - a * c * u - a * d * v - b * w * w + c * v * w -
              d * u * w,
          -a * a * c + a * b * u - a * d * w + b * v * w - c * v * v +
              d * u * v,
          -a * a * d + a * b * v + a * c * w - b * u * w + c * u * v -
              d * u * u,
          -a * a * b + a * c * u + a * d * v - b * w * w + c * v * w -
              d * u * w,
          a * a * a - a * c * c - a * d * d + a * w * w,
          -a * a * u + a * b * c - a * v * w + b * d * w - c * d * v +
              d * d * u,
          -a * a * v + a * b * d + a * u * w - b * c * w + c * c * v -
              c * d * u,
          -a * a * c - a * b * u + a * d * w + b * v * w - c * v * v +
              d * u * v,
          a * a * u + a * b * c - a * v * w - b * d * w + c * d * v - d * d * u,
          a * a * a - a * b * b - a * d * d + a * v * v,
          -a * a * w + a * c * d - a * u * v + b * b * w - b * c * v +
              b * d * u,
          -a * a * d - a * b * v - a * c * w - b * u * w + c * u * v -
              d * u * u,
          a * a * v + a * b * d + a * u * w + b * c * w - c * c * v + c * d * u,
          a * a * w + a * c * d - a * u * v - b * b * w + b * c * v - b * d * u,
          a * a * a - a * b * b - a * c * c + a * u * u)
             .finished() /
         (a * a * a * a - a * a * b * b - a * a * c * c - a * a * d * d +
          a * a * u * u + a * a * v * v + a * a * w * w - b * b * w * w +
          2 * b * c * v * w - 2 * b * d * u * w - c * c * v * v +
          2 * c * d * u * v - d * d * u * u);
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

    if (b == 0. and c == 0. and u == 0.)
      T.Mat3(i).noalias() = Eigen::Matrix3d::Identity() * exp_a;
    else {
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
      const Numeric C3 =
          both_zero ? 1.0 / 6.0
                    : ((x_zero ? 1.0 - sy * iy
                               : y_zero ? sx * ix - 1.0 : sx * ix - sy * iy) *
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
          const Numeric dsx = real ? cx * dx : -cx * dx;
          const Numeric dcx = real ? sx * dx : sx * dx;

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
          const Numeric dsx = real ? cx * dx : -cx * dx;
          const Numeric dcx = real ? sx * dx : sx * dx;

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
      const Complex C3c =
          both_zero ? 1.0 / 6.0
                    : (x_zero ? 1.0 - sy * iy
                              : y_zero ? sx * ix - 1.0 : sx * ix - sy * iy) *
                          inv_x2y2;

      const Numeric& C0 = reinterpret_cast<const Numeric(&)[2]>(C0c)[0];
      const Numeric& C1 = reinterpret_cast<const Numeric(&)[2]>(C1c)[0];
      const Numeric& C2 = reinterpret_cast<const Numeric(&)[2]>(C2c)[0];
      const Numeric& C3 = reinterpret_cast<const Numeric(&)[2]>(C3c)[0];
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
              both_zero ? 0.0
                        : ((x_zero ? -dsy * iy - sy * diy
                                   : y_zero ? dsx * ix + sx * dix
                                            : dsx * ix + sx * dix - dsy * iy -
                                                  sy * diy) -
                           C3c * dx2dy2) *
                              inv_x2y2;

          const Numeric& dC0 = reinterpret_cast<const Numeric(&)[2]>(dC0c)[0];
          const Numeric& dC1 = reinterpret_cast<const Numeric(&)[2]>(dC1c)[0];
          const Numeric& dC2 = reinterpret_cast<const Numeric(&)[2]>(dC2c)[0];
          const Numeric& dC3 = reinterpret_cast<const Numeric(&)[2]>(dC3c)[0];
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
              both_zero ? 0.0
                        : ((x_zero ? -dsy * iy - sy * diy
                                   : y_zero ? dsx * ix + sx * dix
                                            : dsx * ix + sx * dix - dsy * iy -
                                                  sy * diy) -
                           C3c * dx2dy2) *
                              inv_x2y2;

          const Numeric& dC0 = reinterpret_cast<const Numeric(&)[2]>(dC0c)[0];
          const Numeric& dC1 = reinterpret_cast<const Numeric(&)[2]>(dC1c)[0];
          const Numeric& dC2 = reinterpret_cast<const Numeric(&)[2]>(dC2c)[0];
          const Numeric& dC3 = reinterpret_cast<const Numeric(&)[2]>(dC3c)[0];
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
                     const bool& jacobian_do) {
  for (Index i = 0; i < K.NumberOfFrequencies(); i++) {
    if (K.IsRotational(i)) {
      J.SetZero(i);
      if (jacobian_do) {
        for (Index j = 0; j < jacobian_quantities.nelem(); j++)
          if (jacobian_quantities[j].Analytical()) dJ[j].SetZero(i);
      }
    } else {
      J.setSource(a, B, S, i);
      switch (J.StokesDim()) {
        case 4: {
          const auto invK = inv4(K.Kjj()[i],
                                 K.K12()[i],
                                 K.K13()[i],
                                 K.K14()[i],
                                 K.K23()[i],
                                 K.K24()[i],
                                 K.K34()[i]);
          J.Vec4(i) = invK * J.Vec4(i);
          if (jacobian_do)
            for (Index j = 0; j < jacobian_quantities.nelem(); j++)
              if (jacobian_quantities[j].Analytical())
                dJ[j].Vec4(i).noalias() =
                    0.5 * invK *
                    (vector4(a,
                             B,
                             da[j],
                             dB_dT,
                             dS[j],
                             jacobian_quantities[j] == JacPropMatType::Temperature,
                             i) -
                     matrix4(dK[j].Kjj()[i],
                             dK[j].K12()[i],
                             dK[j].K13()[i],
                             dK[j].K14()[i],
                             dK[j].K23()[i],
                             dK[j].K24()[i],
                             dK[j].K34()[i]) *
                         J.Vec4(i));
        } break;
        case 3: {
          const auto invK =
              inv3(K.Kjj()[i], K.K12()[i], K.K13()[i], K.K23()[i]);
          J.Vec3(i) = invK * J.Vec3(i);
          if (jacobian_do)
            for (Index j = 0; j < jacobian_quantities.nelem(); j++)
              if (jacobian_quantities[j].Analytical())
                dJ[j].Vec3(i).noalias() =
                    0.5 * invK *
                    (vector3(a,
                             B,
                             da[j],
                             dB_dT,
                             dS[j],
                             jacobian_quantities[j] == JacPropMatType::Temperature,
                             i) -
                     matrix3(dK[j].Kjj()[i],
                             dK[j].K12()[i],
                             dK[j].K13()[i],
                             dK[j].K23()[i]) *
                         J.Vec3(i));
        } break;
        case 2: {
          const auto invK = inv2(K.Kjj()[i], K.K12()[i]);
          J.Vec2(i) = invK * J.Vec2(i);
          if (jacobian_do)
            for (Index j = 0; j < jacobian_quantities.nelem(); j++)
              if (jacobian_quantities[j].Analytical())
                dJ[j].Vec2(i).noalias() =
                    0.5 * invK *
                    (vector2(a,
                             B,
                             da[j],
                             dB_dT,
                             dS[j],
                             jacobian_quantities[j] == JacPropMatType::Temperature,
                             i) -
                     matrix2(dK[j].Kjj()[i], dK[j].K12()[i]) * J.Vec2(i));
        } break;
        default: {
          const auto invK = 1 / K.Kjj()[i];
          J.Vec1(i)[0] *= invK;
          if (jacobian_do)
            for (Index j = 0; j < jacobian_quantities.nelem(); j++)
              if (jacobian_quantities[j].Analytical())
                dJ[j].Vec1(i)[0] =
                    0.5 * invK *
                    (vector1(a,
                             B,
                             da[j],
                             dB_dT,
                             dS[j],
                             jacobian_quantities[j] == JacPropMatType::Temperature,
                             i) -
                     dK[j].Kjj()[i] * J.Vec1(i)[0]);
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
                             const ArrayOfTransmissionMatrix& dT2,
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
    
    case RadiativeTransferSolver::WeightedEmission: {
      I.rem_weighted(J1, J2, T);
      /* ADD DERIVATIVES */
      I.leftMul(T);
      I.add_weighted(J1, J2, T);
    } break;
  }
}

ArrayOfTransmissionMatrix cumulative_transmission(
    const ArrayOfTransmissionMatrix& T,
    const CumulativeTransmission type) /*[[expects: T.nelem()>0]]*/
{
  const Index n = T.nelem();
  const Index nf = n ? T[0].Frequencies() : 1;
  const Index ns = n ? T[0].StokesDim() : 1;
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
    const RadiationVector &I_incoming,
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
  const Index ns = np ? I[0].StokesDim() : 1;
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
  
  switch(solver) {
    case BackscatterSolver::CommutativeTransmission: {
      // Forward and backwards transmission derivatives
      switch(ns) { 
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
                      T[ip+1].Mat1(iv).inverse() *
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
                      T[ip+1].Mat2(iv).inverse() *
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
                      T[ip+1].Mat3(iv).inverse() *
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
                      T[ip+1].Mat4(iv).inverse() *
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
      switch(ns) {
        case 1: {
          // This is the same as CommutativeTransmission, so use that code
          goto BackscatterSolverCommutativeTransmissionStokesDimOne;
        } break;
        case 2: {
          for (Index ip = 0; ip < np; ip++) {
            for (Index j = ip; j < np; j++) {
              for (Index iq = 0; iq < nq; iq++) {
                for (Index iv = 0; iv < nv; iv++) {
                  if(ip > 1) {
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                      PiTr[ip-2].Mat2(iv) * 
                      (dT1[ip-1][iq].Mat2(iv) + dT2[ip][iq].Mat2(iv)) * 
                      PiTr[ip-1].Mat2(iv).inverse() * I[j].Vec2(iv);
                      
                    if(j < np - 1)
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                      PiTr[ip].Mat2(iv) * Z[ip].Mat2(iv) * PiTf[ip-2].Mat2(iv) *
                      (dT1[ip][iq].Mat2(iv) + dT2[ip+1][iq].Mat2(iv)) * 
                      PiTf[ip-1].Mat2(iv).inverse() * PiTf[ip].Mat2(iv) *
                      I[0].Vec2(iv);
                  }
                  else {
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                      (dT1[ip-1][iq].Mat2(iv) + dT2[ip][iq].Mat2(iv)) * 
                      PiTr[ip-1].Mat2(iv).inverse() * I[j].Vec2(iv);
                    
                    if(j < np - 1)
                    dI[ip][j][iq].Vec2(iv).noalias() +=
                      PiTr[ip].Mat2(iv) * Z[ip].Mat2(iv) *
                      (dT1[ip][iq].Mat2(iv) + dT2[ip+1][iq].Mat2(iv)) * 
                      PiTf[ip-1].Mat2(iv).inverse() * PiTf[ip].Mat2(iv) *
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

ArrayOfTransmissionMatrix cumulative_backscatter(ConstTensor5View t,
                                                 ConstMatrixView m) {
  const Index ns = t.ncols();
  const Index nv = t.npages();
  const Index np = t.nbooks();
  const Index nd = t.nshelves();

  ArrayOfTransmissionMatrix aotm(np, TransmissionMatrix(nv, ns));
  for (Index ip = 0; ip < np; ip++) {
    aotm[ip].setZero();

    switch (ns) {
      case 4:
        for (Index iv = 0; iv < nv; iv++)
          for (Index id = 0; id < nd; id++)
            aotm[ip].Mat4(iv).noalias() +=
                m(id, ip) * matrix4(t(id, ip, iv, joker, joker));
        break;
      case 3:
        for (Index iv = 0; iv < nv; iv++)
          for (Index id = 0; id < nd; id++)
            aotm[ip].Mat3(iv).noalias() +=
                m(id, ip) * matrix3(t(id, ip, iv, joker, joker));
        break;
      case 2:
        for (Index iv = 0; iv < nv; iv++)
          for (Index id = 0; id < nd; id++)
            aotm[ip].Mat2(iv).noalias() +=
                m(id, ip) * matrix2(t(id, ip, iv, joker, joker));
        break;
      case 1:
        for (Index iv = 0; iv < nv; iv++)
          for (Index id = 0; id < nd; id++)
            aotm[ip].Mat1(iv).noalias() +=
                m(id, ip) * matrix1(t(id, ip, iv, joker, joker));
        break;
    }
  }
  return aotm;
}

ArrayOfArrayOfTransmissionMatrix cumulative_backscatter_derivative(
    ConstTensor5View t, const ArrayOfMatrix& aom) {
  const Index ns = t.ncols();
  const Index nv = t.npages();
  const Index np = t.nbooks();
  const Index nd = t.nshelves();
  const Index nq = aom.nelem();

  ArrayOfArrayOfTransmissionMatrix aoaotm(
      np, ArrayOfTransmissionMatrix(nq, TransmissionMatrix(nv, ns)));
  for (Index ip = 0; ip < np; ip++) {
    for (Index iq = 0; iq < nq; iq++) {
      aoaotm[ip][iq].setZero();

      if(not aom[iq].empty()) {
        switch (ns) {
          case 4:
            for (Index iv = 0; iv < nv; iv++)
              for (Index id = 0; id < nd; id++)
                aoaotm[ip][iq].Mat4(iv).noalias() +=
                    aom[iq](id, ip) * matrix4(t(id, ip, iv, joker, joker));
            break;
          case 3:
            for (Index iv = 0; iv < nv; iv++)
              for (Index id = 0; id < nd; id++)
                aoaotm[ip][iq].Mat3(iv).noalias() +=
                    aom[iq](id, ip) * matrix3(t(id, ip, iv, joker, joker));
            break;
          case 2:
            for (Index iv = 0; iv < nv; iv++)
              for (Index id = 0; id < nd; id++)
                aoaotm[ip][iq].Mat2(iv).noalias() +=
                    aom[iq](id, ip) * matrix2(t(id, ip, iv, joker, joker));
            break;
          case 1:
            for (Index iv = 0; iv < nv; iv++)
              for (Index id = 0; id < nd; id++)
                aoaotm[ip][iq].Mat1(iv).noalias() +=
                    aom[iq](id, ip) * matrix1(t(id, ip, iv, joker, joker));
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

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfTransmissionMatrix& atm) {
  for (const auto& T : atm) os << T << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfTransmissionMatrix& aatm) {
  for (const auto& T : aatm) os << T << '\n';
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

std::ostream& operator<<(std::ostream& os, const ArrayOfRadiationVector& arv) {
  for (const auto& R : arv) os << R << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os,
                         const ArrayOfArrayOfRadiationVector& aarv) {
  for (const auto& R : aarv) os << R << '\n';
  return os;
}

std::istream& operator>>(std::istream& data, TransmissionMatrix& tm) {
  for (auto& T : tm.T4)
    data >> T(0, 0) >> T(0, 1) >> T(0, 2) >> T(0, 3) >> T(1, 0) >> T(1, 1) >>
        T(1, 2) >> T(1, 3) >> T(2, 0) >> T(2, 1) >> T(2, 2) >> T(2, 3) >>
        T(3, 0) >> T(3, 1) >> T(3, 2) >> T(3, 3);
  for (auto& T : tm.T3)
    data >> T(0, 0) >> T(0, 1) >> T(0, 2) >> T(1, 0) >> T(1, 1) >> T(1, 2) >>
        T(2, 0) >> T(2, 1) >> T(2, 2);
  for (auto& T : tm.T2) data >> T(0, 0) >> T(0, 1) >> T(1, 0) >> T(1, 1);
  for (auto& T : tm.T1) data >> T(0, 0);
  return data;
}

std::istream& operator>>(std::istream& data, RadiationVector& rv) {
  for (auto& R : rv.R4) data >> R[0] >> R[1] >> R[2] >> R[3];
  for (auto& R : rv.R3) data >> R[0] >> R[1] >> R[2];
  for (auto& R : rv.R2) data >> R[0] >> R[1] >> R[2];
  for (auto& R : rv.R1) data >> R[0];
  return data;
}

TransmissionMatrix::TransmissionMatrix(const PropagationMatrix& pm,
                                       const Numeric& r) {
  *this = TransmissionMatrix(pm.NumberOfFrequencies(), pm.StokesDimensions());
  transmat(*this, pm, pm, r);  // Slower to compute, faster to implement...
}
