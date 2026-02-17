#pragma once

#include <array.h>
#include <matpack.h>
#include <xml.h>

namespace rtepack {

//! A 4x4 matrix of Complex values to be used as a Mueller Matrix
struct specmat final : ComplexMatrix44 {
  constexpr specmat(Complex tau = 1.0) noexcept
      : ComplexMatrix44{
            tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau} {}

  constexpr specmat(Complex m00,
                    Complex m01,
                    Complex m02,
                    Complex m03,
                    Complex m10,
                    Complex m11,
                    Complex m12,
                    Complex m13,
                    Complex m20,
                    Complex m21,
                    Complex m22,
                    Complex m23,
                    Complex m30,
                    Complex m31,
                    Complex m32,
                    Complex m33) noexcept
      : ComplexMatrix44{m00,
                        m01,
                        m02,
                        m03,
                        m10,
                        m11,
                        m12,
                        m13,
                        m20,
                        m21,
                        m22,
                        m23,
                        m30,
                        m31,
                        m32,
                        m33} {}

  constexpr specmat(std::array<Complex, 16> data) noexcept
      : ComplexMatrix44{data} {}

  //! The identity matrix
  static constexpr specmat id() { return specmat{1.0}; }

  constexpr specmat &operator+=(const specmat &b) {
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

  constexpr specmat &operator-=(const specmat &b) {
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

  constexpr specmat &operator*=(const specmat &b) {
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

  constexpr specmat operator-() const {
    return specmat{-data[0],
                   -data[1],
                   -data[2],
                   -data[3],
                   -data[4],
                   -data[5],
                   -data[6],
                   -data[7],
                   -data[8],
                   -data[9],
                   -data[10],
                   -data[11],
                   -data[12],
                   -data[13],
                   -data[14],
                   -data[15]};
  }
};

//! Addition between specmat matrices
constexpr specmat operator+(specmat a, const specmat &b) { return a += b; }

constexpr specmat operator+(Complex a, specmat b) { return specmat{a} + b; }

//! Subtraction between specmat matrices
constexpr specmat operator-(specmat a, const specmat &b) { return a -= b; }

constexpr specmat operator-(Complex a, specmat b) { return specmat{a} - b; }

//! Scaling a specmat matrix
constexpr specmat operator*(specmat a, const Complex &b) { return a *= b; }

//! Scaling a specmat matrix
constexpr specmat operator*(const Complex &a, specmat b) { return b * a; }

//! Scaling a specmat matrix
constexpr specmat operator/(specmat a, const Complex &b) {
  a *= (1.0 / b);
  return a;
}

//! Scaling a specmat matrix
constexpr specmat operator*(specmat a, const specmat &b) { return a *= b; }

//! Take the average of two specmat matrices
constexpr specmat avg(specmat a, const specmat &b) {
  a += b;
  a *= Complex{0.5};
  return a;
}

constexpr Complex det(const specmat &A) {
  const auto [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p] = A;

  return a * (f * (k * p - l * o) + g * (l * n - j * p) + h * (j * o - k * n)) +
         b * (e * (l * o - k * p) + g * (i * p - l * m) + h * (k * m - i * o)) +
         c * (e * (j * p - l * n) + f * (l * m - i * p) + h * (i * n - j * m)) +
         d * (e * (k * n - j * o) + f * (i * o - k * m) + g * (j * m - i * n));
}

constexpr Complex tr(const specmat &A) {
  return A[0, 0] + A[1, 1] + A[2, 2] + A[3, 3];
}

constexpr Complex midtr(const specmat &A) {
  return {std::midpoint(std::midpoint(A[0, 0].real(), A[1, 1].real()),
                        std::midpoint(A[2, 2].real(), A[3, 3].real())),
          std::midpoint(std::midpoint(A[0, 0].imag(), A[1, 1].imag()),
                        std::midpoint(A[2, 2].imag(), A[3, 3].imag()))};
}

constexpr specmat adj(const specmat &A) {
  const auto [a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p] = A;

  return specmat{
      f * (k * p - l * o) + g * (l * n - j * p) + h * (j * o - k * n),
      b * (l * o - k * p) + c * (j * p - l * n) + d * (k * n - j * o),
      b * (g * p - h * o) + c * (h * n - f * p) + d * (f * o - g * n),
      b * (h * k - g * l) + c * (f * l - h * j) + d * (g * j - f * k),
      e * (l * o - k * p) + g * (i * p - l * m) + h * (k * m - i * o),
      a * (k * p - l * o) + c * (l * m - i * p) + d * (i * o - k * m),
      a * (h * o - g * p) + c * (e * p - h * m) + d * (g * m - e * o),
      a * (g * l - h * k) + c * (h * i - e * l) + d * (e * k - g * i),
      e * (j * p - l * n) + f * (l * m - i * p) + h * (i * n - j * m),
      a * (l * n - j * p) + b * (i * p - l * m) + d * (j * m - i * n),
      a * (f * p - h * n) + b * (h * m - e * p) + d * (e * n - f * m),
      a * (h * j - f * l) + b * (e * l - h * i) + d * (f * i - e * j),
      e * (k * n - j * o) + f * (i * o - k * m) + g * (j * m - i * n),
      a * (j * o - k * n) + b * (k * m - i * o) + c * (i * n - j * m),
      a * (g * n - f * o) + b * (e * o - g * m) + c * (f * m - e * n),
      a * (f * k - g * j) + b * (g * i - e * k) + c * (e * j - f * i)};
}

constexpr specmat inv(const specmat &A) { return adj(A) / det(A); }

//! Element-wise Dawson function of a specmat matrix (FIXME: make matrix-wise?)
specmat dawson(const specmat &A);

//! Element-wise erfcx function of a specmat matrix (FIXME: make matrix-wise?)
specmat erfcx(const specmat &m);

constexpr specmat elem_prod(const specmat &a, const specmat &b) {
  specmat c;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      c[i, j] = a[i, j] * b[i, j];
    }
  }
  return c;
}

using specmat_vector            = matpack::data_t<specmat, 1>;
using specmat_vector_view       = matpack::view_t<specmat, 1>;
using specmat_vector_const_view = matpack::view_t<const specmat, 1>;

using specmat_matrix            = matpack::data_t<specmat, 2>;
using specmat_matrix_view       = matpack::view_t<specmat, 2>;
using specmat_matrix_const_view = matpack::view_t<const specmat, 2>;

using specmat_tensor3            = matpack::data_t<specmat, 3>;
using specmat_tensor3_view       = matpack::view_t<specmat, 3>;
using specmat_tensor3_const_view = matpack::view_t<const specmat, 3>;
}  // namespace rtepack

template <>

struct std::formatter<rtepack::specmat> : std::formatter<ComplexMatrix44> {};

template <>
struct xml_io_stream<rtepack::specmat> : xml_io_stream<ComplexMatrix44> {};
