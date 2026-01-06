#pragma once

#include "rtepack_mueller_matrix.h"
#include "rtepack_propagation_matrix.h"
#include "rtepack_spectral_matrix.h"

namespace rtepack {
struct tran {
  Numeric a, exp_a;
  Numeric b, c, d, u, v, w;
  Numeric b2, c2, d2, u2, v2, w2;
  Numeric B, C, S;
  Numeric x2, y2, x, y, cy, sy, cx, sx;
  Numeric ix, iy, inv_x2y2;
  Numeric C0, C1, C2, C3;
  bool polarized, x_zero, y_zero, both_zero, either_zero;

  constexpr tran() = default;

  tran(const propmat &k1, const propmat &k2, const Numeric r);

  [[nodiscard]] muelmat operator()() const noexcept;
  [[nodiscard]] muelmat expm1() const noexcept;
  [[nodiscard]] muelmat evolve_operator() const noexcept;

  [[nodiscard]] muelmat deriv(const muelmat &t,
                              const propmat &k1,
                              const propmat &k2,
                              const propmat &dk,
                              const Numeric r,
                              const Numeric dr) const;
  [[nodiscard]] muelmat evolve_operator_deriv(const muelmat &l,
                                              const propmat &dk,
                                              const muelmat &dt,
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

muelmat exp(propmat k, Numeric r = 1.0);

propmat logK(const muelmat &m);

specmat sqrt(const propmat &pm);

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

void two_level_exp(std::vector<muelmat_vector> &T,
                   std::vector<muelmat_vector> &L,
                   std::vector<muelmat_tensor3> &dT,
                   std::vector<muelmat_tensor3> &dL,
                   const std::vector<propmat_vector> &K,
                   const std::vector<propmat_matrix> &dK,
                   const Vector &r,
                   const Tensor3 &dr);

void two_level_exp(std::vector<muelmat_vector> &T,
                   std::vector<muelmat_vector> &L0,
                   std::vector<muelmat_vector> &L1,
                   std::vector<muelmat_tensor3> &dT,
                   std::vector<muelmat_tensor3> &dL0,
                   std::vector<muelmat_tensor3> &dL1,
                   const std::vector<propmat_vector> &K,
                   const std::vector<propmat_matrix> &dK,
                   const Vector &r,
                   const Tensor3 &dr);
}  // namespace rtepack
