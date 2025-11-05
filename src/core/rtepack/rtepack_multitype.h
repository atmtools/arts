#pragma once

#include "rtepack_concepts.h"
#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
//! Treat the propagation matrix as a Mueller matrix
constexpr muelmat to_muelmat(const propmat &k) {
  const auto [a, b, c, d, u, v, w] = k.data;
  return muelmat{a, b, c, d, b, a, u, v, c, -u, a, w, d, -v, -w, a};
}

//! Treat the propagation matrix as a Stokes (absorption) vector
constexpr stokvec absvec(const propmat &k) {
  return stokvec{k.A(), k.B(), k.C(), k.D()};
}

//! Treat a list of propagation matrices as a list of Stokes (absorption) vectors
stokvec_vector absvec(const propmat_vector_const_view &k);

//! Get the inverse of a propagation matrix as a Mueller matrix
constexpr muelmat inv(const propmat &k) {
  const auto &[a, b, c, d, u, v, w] = k.data;

  const Numeric div =
      1.0 / (a * a * (a * a - b * b - c * c - d * d + u * u + v * v + w * w) -
             b * w * (b * w - 2 * c * v + 2 * d * u) +
             c * v * (-c * v + 2 * d * u) - d * d * u * u);

  return {a * (a * a + u * u + v * v + w * w) * div,
          -(a * (a * b + c * u + d * v) + w * (b * w - c * v + d * u)) * div,
          (a * (-a * c + b * u - d * w) + v * (b * w - c * v + d * u)) * div,
          (a * (-a * d + b * v + c * w) - u * (b * w - c * v + d * u)) * div,
          (a * (-a * b + c * u + d * v) - w * (b * w - c * v + d * u)) * div,
          a * (a * a - c * c - d * d + w * w) * div,
          (d * (b * w - c * v + d * u) - a * (a * u - b * c + v * w)) * div,
          (a * (-a * v + b * d + u * w) - c * (b * w - c * v + d * u)) * div,
          (v * (b * w - c * v + d * u) - a * (a * c + b * u - d * w)) * div,
          (a * (a * u + b * c - v * w) - d * (b * w - c * v + d * u)) * div,
          a * (a * a - b * b - d * d + v * v) * div,
          (b * (b * w - c * v + d * u) - a * (a * w - c * d + u * v)) * div,
          -(a * (a * d + b * v + c * w) + u * (b * w - c * v + d * u)) * div,
          (a * (a * v + b * d + u * w) + c * (b * w - c * v + d * u)) * div,
          -(a * (-a * w - c * d + u * v) + b * (b * w - c * v + d * u)) * div,
          a * (a * a - b * b - c * c + u * u) * div};
}

//! muelmat matrix multiplied by a stokvec vector
constexpr stokvec operator*(const muelmat &a, const stokvec &b) {
  const auto
      &[m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16] =
          a;
  const auto &[s1, s2, s3, s4] = b;
  return {m1 * s1 + m2 * s2 + m3 * s3 + m4 * s4,
          m5 * s1 + m6 * s2 + m7 * s3 + m8 * s4,
          m10 * s2 + m11 * s3 + m12 * s4 + m9 * s1,
          m13 * s1 + m14 * s2 + m15 * s3 + m16 * s4};
}

//! Mutliply a propmat with a muelmat matrix
constexpr muelmat operator*(const propmat &k, const muelmat &m) {
  const auto [a, b, c, d, u, v, w] = k.data;
  const auto
      [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16] =
          m.data;

  return muelmat{a * m1 + b * m5 + c * m9 + d * m13,
                 a * m2 + b * m6 + c * m10 + d * m14,
                 a * m3 + b * m7 + c * m11 + d * m15,
                 a * m4 + b * m8 + c * m12 + d * m16,
                 a * m5 + b * m1 + m13 * v + m9 * u,
                 a * m6 + b * m2 + m10 * u + m14 * v,
                 a * m7 + b * m3 + m11 * u + m15 * v,
                 a * m8 + b * m4 + m12 * u + m16 * v,
                 a * m9 + c * m1 + m13 * w - m5 * u,
                 a * m10 + c * m2 + m14 * w - m6 * u,
                 a * m11 + c * m3 + m15 * w - m7 * u,
                 a * m12 + c * m4 + m16 * w - m8 * u,
                 a * m13 + d * m1 - m5 * v - m9 * w,
                 a * m14 + d * m2 - m10 * w - m6 * v,
                 a * m15 + d * m3 - m11 * w - m7 * v,
                 a * m16 + d * m4 - m12 * w - m8 * v};
}

//! Mutliply a muelmat matrix with a propmat
constexpr muelmat operator*(const muelmat &m, const propmat &k) {
  const auto [a, b, c, d, u, v, w] = k.data;
  const auto
      [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16] =
          m.data;

  return muelmat{a * m1 + b * m2 + c * m3 + d * m4,
                 a * m2 + b * m1 - m3 * u - m4 * v,
                 a * m3 + c * m1 + m2 * u - m4 * w,
                 a * m4 + d * m1 + m2 * v + m3 * w,
                 a * m5 + b * m6 + c * m7 + d * m8,
                 a * m6 + b * m5 - m7 * u - m8 * v,
                 a * m7 + c * m5 + m6 * u - m8 * w,
                 a * m8 + d * m5 + m6 * v + m7 * w,
                 a * m9 + b * m10 + c * m11 + d * m12,
                 a * m10 + b * m9 - m11 * u - m12 * v,
                 a * m11 + c * m9 + m10 * u - m12 * w,
                 a * m12 + d * m9 + m10 * v + m11 * w,
                 a * m13 + b * m14 + c * m15 + d * m16,
                 a * m14 + b * m13 - m15 * u - m16 * v,
                 a * m15 + c * m13 + m14 * u - m16 * w,
                 a * m16 + d * m13 + m14 * v + m15 * w};
}

//! Multiply a propagation matrix with a stokvec vector
constexpr stokvec operator*(const propmat &k, const stokvec s) {
  const auto [a, b, c, d, u, v, w] = k.data;
  const auto [s1, s2, s3, s4]      = s.data;

  return stokvec{a * s1 + b * s2 + c * s3 + d * s4,
                 a * s2 + b * s1 + s3 * u + s4 * v,
                 a * s3 + c * s1 - s2 * u + s4 * w,
                 a * s4 + d * s1 - s2 * v - s3 * w};
}

//! Transform a matrix of shape (N, 4) to a list of Stokes (absorption) vectors
stokvec_vector to_stokvec_vector(const ConstMatrixView &v);

//! Transform a list Mueller matrices to a Tensor3 of shape (N, 4, 4)
Tensor3 to_tensor3(const muelmat_vector_const_view &m);

//! Transform a list of Stokes (absorption) vectors to a matrix of shape (N, 4)
Matrix to_matrix(const stokvec_vector_const_view &v);

//! Transform a propagation matrix to a matrix of shape (4, 4)
Matrix to_matrix(const propmat &v);

//! Transform a Stokes (absorption) vector to a vector of shape (4)
Vector to_vector(const stokvec &v);

//! Transform a propmat to a stokvec
stokvec to_stokvec(const propmat &a);

//! Transform a matpack type to a stokvec
stokvec to_stokvec(const ConstVectorView &a);

//! Transform a matpack type to a propmat
propmat to_propmat(const ConstMatrixView &a);

//! Transform a matpack type to a muelmat
muelmat to_muelmat(const ConstMatrixView &a);
}  // namespace rtepack
