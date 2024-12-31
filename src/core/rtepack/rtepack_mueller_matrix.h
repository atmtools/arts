#pragma once

#include "array.h"
#include "rtepack_concepts.h"

namespace rtepack {

//! A 4x4 matrix of Numeric values to be used as a Mueller Matrix
struct muelmat final : mat44 {
  constexpr muelmat(Numeric tau = 1.0) noexcept
      : mat44{tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau} {}

  constexpr muelmat(Numeric a,
                    Numeric b,
                    Numeric c,
                    Numeric d,
                    Numeric e,
                    Numeric f,
                    Numeric g,
                    Numeric h,
                    Numeric i,
                    Numeric j,
                    Numeric k,
                    Numeric l,
                    Numeric m,
                    Numeric n,
                    Numeric o,
                    Numeric p) noexcept
      : mat44{a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p} {}

  constexpr muelmat(std::array<Numeric, 16> data) noexcept : mat44{data} {}

  //! The identity matrix
  static constexpr muelmat id() { return muelmat{1.0}; }

  constexpr muelmat &operator+=(const muelmat &b) {
    data[0]  += b.data[0];
    data[1]  += b.data[1];
    data[2]  += b.data[2];
    data[3]  += b.data[3];
    data[4]  += b.data[4];
    data[5]  += b.data[5];
    data[6]  += b.data[6];
    data[7]  += b.data[7];
    data[8]  += b.data[8];
    data[9]  += b.data[9];
    data[10] += b.data[10];
    data[11] += b.data[11];
    data[12] += b.data[12];
    data[13] += b.data[13];
    data[14] += b.data[14];
    data[15] += b.data[15];
    return *this;
  }

  constexpr muelmat &operator-=(const muelmat &b) {
    data[0]  -= b.data[0];
    data[1]  -= b.data[1];
    data[2]  -= b.data[2];
    data[3]  -= b.data[3];
    data[4]  -= b.data[4];
    data[5]  -= b.data[5];
    data[6]  -= b.data[6];
    data[7]  -= b.data[7];
    data[8]  -= b.data[8];
    data[9]  -= b.data[9];
    data[10] -= b.data[10];
    data[11] -= b.data[11];
    data[12] -= b.data[12];
    data[13] -= b.data[13];
    data[14] -= b.data[14];
    data[15] -= b.data[15];
    return *this;
  }

  constexpr muelmat &operator*=(const muelmat &b) {
    const auto [a00,
                a01,
                a02,
                a03,
                a10,
                a11,
                a12,
                a13,
                a20,
                a21,
                a22,
                a23,
                a30,
                a31,
                a32,
                a33] = data;
    const auto [b00,
                b01,
                b02,
                b03,
                b10,
                b11,
                b12,
                b13,
                b20,
                b21,
                b22,
                b23,
                b30,
                b31,
                b32,
                b33] = b;

    return *this = {a00 * b00 + a01 * b10 + a02 * b20 + a03 * b30,
                    a00 * b01 + a01 * b11 + a02 * b21 + a03 * b31,
                    a00 * b02 + a01 * b12 + a02 * b22 + a03 * b32,
                    a00 * b03 + a01 * b13 + a02 * b23 + a03 * b33,
                    a10 * b00 + a11 * b10 + a12 * b20 + a13 * b30,
                    a10 * b01 + a11 * b11 + a12 * b21 + a13 * b31,
                    a10 * b02 + a11 * b12 + a12 * b22 + a13 * b32,
                    a10 * b03 + a11 * b13 + a12 * b23 + a13 * b33,
                    a20 * b00 + a21 * b10 + a22 * b20 + a23 * b30,
                    a20 * b01 + a21 * b11 + a22 * b21 + a23 * b31,
                    a20 * b02 + a21 * b12 + a22 * b22 + a23 * b32,
                    a20 * b03 + a21 * b13 + a22 * b23 + a23 * b33,
                    a30 * b00 + a31 * b10 + a32 * b20 + a33 * b30,
                    a30 * b01 + a31 * b11 + a32 * b21 + a33 * b31,
                    a30 * b02 + a31 * b12 + a32 * b22 + a33 * b32,
                    a30 * b03 + a31 * b13 + a32 * b23 + a33 * b33};
  }
};

//! Addition between muelmat matrices
constexpr muelmat operator+(muelmat a, const muelmat &b) { return a += b; }

constexpr muelmat operator+(Numeric a, muelmat b) { return muelmat{a} + b; }

//! Subtraction between muelmat matrices
constexpr muelmat operator-(muelmat a, const muelmat &b) { return a -= b; }

constexpr muelmat operator-(Numeric a, muelmat b) { return muelmat{a} - b; }

//! Scaling a muelmat matrix
constexpr muelmat operator*(muelmat a, const Numeric &b) { return a *= b; }

//! Scaling a muelmat matrix
constexpr muelmat operator*(const Numeric &a, muelmat b) { return b * a; }

//! Scaling a muelmat matrix
constexpr muelmat operator/(muelmat a, const Numeric &b) {
  a /= b;
  return a;
}

//! Scaling a muelmat matrix
constexpr muelmat operator*(muelmat a, const muelmat &b) { return a *= b; }

//! Take the average of two muelmat matrices
constexpr muelmat avg(muelmat a, const muelmat &b) {
  a += b;
  a *= 0.5;
  return a;
}

constexpr muelmat inv(const muelmat &A) {
  const auto [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p] = A;

  const Numeric div =
      1.0 / (a * f * k * p - a * f * l * o - a * g * j * p + a * g * l * n +
             a * h * j * o - a * h * k * n - b * e * k * p + b * e * l * o +
             b * g * i * p - b * g * l * m - b * h * i * o + b * h * k * m +
             c * e * j * p - c * e * l * n - c * f * i * p + c * f * l * m +
             c * h * i * n - c * h * j * m - d * e * j * o + d * e * k * n +
             d * f * i * o - d * f * k * m - d * g * i * n + d * g * j * m);

  return div * muelmat{(f * k * p - f * l * o - g * j * p + g * l * n +
                        h * j * o - h * k * n),
                       (-b * k * p + b * l * o + c * j * p - c * l * n -
                        d * j * o + d * k * n),
                       (b * g * p - b * h * o - c * f * p + c * h * n +
                        d * f * o - d * g * n),
                       (-b * g * l + b * h * k + c * f * l - c * h * j -
                        d * f * k + d * g * j),
                       (-e * k * p + e * l * o + g * i * p - g * l * m -
                        h * i * o + h * k * m),
                       (a * k * p - a * l * o - c * i * p + c * l * m +
                        d * i * o - d * k * m),
                       (-a * g * p + a * h * o + c * e * p - c * h * m -
                        d * e * o + d * g * m),
                       (a * g * l - a * h * k - c * e * l + c * h * i +
                        d * e * k - d * g * i),
                       (e * j * p - e * l * n - f * i * p + f * l * m +
                        h * i * n - h * j * m),
                       (-a * j * p + a * l * n + b * i * p - b * l * m -
                        d * i * n + d * j * m),
                       (a * f * p - a * h * n - b * e * p + b * h * m +
                        d * e * n - d * f * m),
                       (-a * f * l + a * h * j + b * e * l - b * h * i -
                        d * e * j + d * f * i),
                       (-e * j * o + e * k * n + f * i * o - f * k * m -
                        g * i * n + g * j * m),
                       (a * j * o - a * k * n - b * i * o + b * k * m +
                        c * i * n - c * j * m),
                       (-a * f * o + a * g * n + b * e * o - b * g * m -
                        c * e * n + c * f * m),
                       (a * f * k - a * g * j - b * e * k + b * g * i +
                        c * e * j - c * f * i)};
}

using muelmat_vector            = matpack::data_t<muelmat, 1>;
using muelmat_vector_view       = matpack::view_t<muelmat, 1>;
using muelmat_vector_const_view = matpack::view_t<const muelmat, 1>;

using muelmat_matrix            = matpack::data_t<muelmat, 2>;
using muelmat_matrix_view       = matpack::view_t<muelmat, 2>;
using muelmat_matrix_const_view = matpack::view_t<const muelmat, 2>;

using muelmat_tensor3            = matpack::data_t<muelmat, 3>;
using muelmat_tensor3_view       = matpack::view_t<muelmat, 3>;
using muelmat_tensor3_const_view = matpack::view_t<const muelmat, 3>;

void forward_cumulative_transmission(Array<muelmat_vector> &Pi,
                                     const Array<muelmat_vector> &T);

Array<muelmat_vector> forward_cumulative_transmission(
    const Array<muelmat_vector> &T);
}  // namespace rtepack

template <>

struct std::formatter<rtepack::muelmat> {
  std::formatter<rtepack::mat44> fmt;

  [[nodiscard]] constexpr auto &inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto &inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const rtepack::muelmat &v,
                              FmtContext &ctx) const {
    return fmt.format(v, ctx);
  }
};
