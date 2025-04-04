#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"

namespace rtepack {
  struct tran {
  Numeric a{}, b{}, c{}, d{}, u{}, v{}, w{};   // To not repeat input
  Numeric exp_a{};                             // To not repeat exp(a)
  Numeric b2{}, c2{}, d2{}, u2{}, v2{}, w2{};  // To shorten expressions
  Numeric B, C, S;  // From L^4 + BL^2 + C = 0; S = sqrt(B^2 - 4C)
  Numeric x2{}, y2{}, x{}, y{}, cy{}, sy{}, cx{},
      sx{};  // Eigenvalues and their used trigonometric functions
  Numeric ix{}, iy{}, inv_x2y2{};  // Computational helpers
  Numeric C0{}, C1{}, C2{}, C3{};  // The Cayley-Hamilton coefficients
  bool unpolarized{}, x_zero{}, y_zero{}, both_zero{}, either_zero{};

  constexpr tran() = default;

  tran(const propmat &k1, const propmat &k2, const Numeric r);

  muelmat operator()() const noexcept;

  [[nodiscard]] muelmat deriv(const muelmat &t,
                              const propmat &k1,
                              const propmat &k2,
                              const propmat &dk,
                              const Numeric r,
                              const Numeric dr) const;
};

void two_level_exp(muelmat &t,
                   muelmat_vector_view dt1,
                   muelmat_vector_view dt2,
                   const propmat &k1,
                   const propmat &k2,
                   const propmat_vector_const_view &dk1,
                   const propmat_vector_const_view &dk2,
                   const Numeric r,
                   const ConstVectorView &dr1,
                   const ConstVectorView &dr2);

void two_level_exp(muelmat_vector_view t,
                   muelmat_matrix_view dt1,
                   muelmat_matrix_view dt2,
                   const propmat_vector_const_view &k1,
                   const propmat_vector_const_view &k2,
                   const propmat_matrix_const_view &dk1,
                   const propmat_matrix_const_view &dk2,
                   const Numeric r,
                   const ConstVectorView &dr1,
                   const ConstVectorView &dr2);

muelmat exp (propmat k, Numeric r=1.0);

void two_level_exp(muelmat_vector_view t,
                   const propmat_vector_const_view &k1,
                   const propmat_vector_const_view &k2,
                   const Numeric r);

void two_level_exp(std::vector<muelmat_vector> &T,
                   std::vector<muelmat_tensor3> &dT,
                   const std::vector<propmat_vector> &K,
                   const std::vector<propmat_matrix> &dK,
                   const Vector &r,
                   const Tensor3 &dr);
}  // namespace rtepack
