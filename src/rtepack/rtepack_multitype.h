#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {
constexpr muelmat to_muelmat(const propmat &k) {
  const auto [a, b, c, d, u, v, w] = k.data;
  return {a, b, c, d, b, a, u, v, c, -u, a, w, d, -v, -w, a};
}
constexpr stokvec absvec(const propmat &k) {
  return {k.A(), k.B(), k.C(), k.D()};
}

constexpr muelmat inv(const propmat &k) {
  const auto [a, b, c, d, u, v, w] = k.data;
  const Numeric a2 = a * a, b2 = b * b, c2 = c * c, u2 = u * u, d2 = d * d,
                v2 = v * v, w2 = w * w;

  const Numeric f = a2 * a2 - a2 * b2 - a2 * c2 - a2 * d2 + a2 * u2 + a2 * v2 +
                    a2 * w2 - b2 * w2 + 2 * b * c * v * w - 2 * b * d * u * w -
                    c2 * v2 + 2 * c * d * u * v - d2 * u2;

  const Numeric div = 1.0 / f;

  return {
      a * (a2 + u2 + v2 + w2) * div,
      (-a2 * b - a * c * u - a * d * v - b * w2 + c * v * w - d * u * w) * div,
      (-a2 * c + a * b * u - a * d * w + b * v * w - c * v2 + d * u * v) * div,
      (-a2 * d + a * b * v + a * c * w - b * u * w + c * u * v - d * u2) * div,
      (-a2 * b + a * c * u + a * d * v - b * w2 + c * v * w - d * u * w) * div,
      a * (a2 - c2 - d2 + w2) * div,
      (-a2 * u + a * b * c - a * v * w + b * d * w - c * d * v + d2 * u) * div,
      (-a2 * v + a * b * d + a * u * w - b * c * w + c2 * v - c * d * u) * div,
      (-a2 * c - a * b * u + a * d * w + b * v * w - c * v2 + d * u * v) * div,
      (a2 * u + a * b * c - a * v * w - b * d * w + c * d * v - d2 * u) * div,
      a * (a2 - b2 - d2 + v2) * div,
      (-a2 * w + a * c * d - a * u * v + b2 * w - b * c * v + b * d * u) * div,
      (-a2 * d - a * b * v - a * c * w - b * u * w + c * u * v - d * u2) * div,
      (a2 * v + a * b * d + a * u * w + b * c * w - c2 * v + c * d * u) * div,
      (a2 * w + a * c * d - a * u * v - b2 * w + b * c * v - b * d * u) * div,
      a * (a2 - b2 - c2 + u2) * div};
}

//! muelmat matrix multiplied by a lazy stokvec vector
constexpr auto operator*(const muelmat &a, const lazy_stokvec auto &b) {
  return smvmul{constexpr_smat_data<Numeric, 4>(a), b};
}

//! muelmat matrix multiplied by a stokvec vector
constexpr auto operator*(const muelmat &a, const stokvec &b) {
  return smvmul{constexpr_smat_data{a}, constexpr_vec_data{b}};
}

//! Lazy muelmat matrix multiplied by a stokvec vector
constexpr auto operator*(const lazy_muelmat auto &a, const stokvec &b) {
  return smvmul{a, constexpr_vec_data{b}};
}

//! Mutliply a propmat with a muelmat matrix
constexpr muelmat operator*(const propmat &k, const muelmat &m) {
  const auto [a, b, c, d, u, v, w] = k.data;
  const auto [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15,
              m16] = m.data;

  return {
      a * m1 + b * m5 + c * m9 + d * m13,  a * m2 + b * m6 + c * m10 + d * m14,
      a * m3 + b * m7 + c * m11 + d * m15, a * m4 + b * m8 + c * m12 + d * m16,
      a * m5 + b * m1 + m13 * v + m9 * u,  a * m6 + b * m2 + m10 * u + m14 * v,
      a * m7 + b * m3 + m11 * u + m15 * v, a * m8 + b * m4 + m12 * u + m16 * v,
      a * m9 + c * m1 + m13 * w - m5 * u,  a * m10 + c * m2 + m14 * w - m6 * u,
      a * m11 + c * m3 + m15 * w - m7 * u, a * m12 + c * m4 + m16 * w - m8 * u,
      a * m13 + d * m1 - m5 * v - m9 * w,  a * m14 + d * m2 - m10 * w - m6 * v,
      a * m15 + d * m3 - m11 * w - m7 * v, a * m16 + d * m4 - m12 * w - m8 * v};
}

//! Mutliply a muelmat matrix with a propmat
constexpr muelmat operator*(const muelmat &m, const propmat &k) {
  const auto [a, b, c, d, u, v, w] = k.data;
  const auto [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15,
              m16] = m.data;

  return {a * m1 + b * m2 + c * m3 + d * m4,
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
  const auto [s1, s2, s3, s4] = s.data;

  return {a * s1 + b * s2 + c * s3 + d * s4, a * s2 + b * s1 + s3 * u + s4 * v,
          a * s3 + c * s1 - s2 * u + s4 * w, a * s4 + d * s1 - s2 * v - s3 * w};
}
} // namespace rtepack
