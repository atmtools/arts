/* The half-size eigenproblem of a mirror-symmetric medium, and the full one
   otherwise.

   VDISORT diagonalizes a layer's operator through its half-size product
   problem when the phase matrix of the Fourier mode is mirror symmetric about
   the horizontal plane (except in the cosine system of order zero, which
   holds the energy-conservation pair).  That problem reads only the upward
   rows of the operator.  Here a problem driven only by boundary radiance,
   whose particular solutions are zero and so read nothing else, is solved
   with symmetric phase matrices and again with the downward rows changed by a
   relative 1e-9 in every mode but the cosine system of order zero.  Taking
   the half-size problem for the changed medium would ignore the change and
   reproduce the symmetric result bit for bit; the full problem it must take
   instead differs by about the change.  The agreement to that level is also
   the agreement of the half-size and full eigensolutions. */

#include <legendre.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <algorithm>
#include <format>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include "lab-frame.h"

namespace {
using vdisort_test::fourier_modes;
using vdisort_test::legendre_matrix;

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

constexpr Index NQ = 8, NF = 4, NL = 3;

Vector streams() {
  Vector mu(NQ), w(NQ / 2);
  Legendre::PositiveDoubleGaussLegendre(mu[Range{0, NQ / 2}], w);
  for (Index i = 0; i < NQ / 2; i++) mu[NQ / 2 + i] = -mu[i];
  return mu;
}

//! Rayleigh with F34 = 0.1 x (1 - x^2) = 0.04 (P1 - P3), regular at forward and backward scattering
CompactPlanarMuelmatVector series() {
  CompactPlanarMuelmatVector c(4);
  c[0] = {1.0, -0.5, 1.0, 0.0, 0.0, 0.0};
  c[1] = {0.0, 0.0, 0.0, 1.5, 0.04, 1.5};
  c[2] = {0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
  c[3] = {0.0, 0.0, 0.0, 0.0, -0.04, 0.0};
  return c;
}

struct phases {
  vdisort::phase_matrix_data      diffuse;
  vdisort::beam_phase_matrix_data beam;
};

phases phase_matrices(Numeric mu0) {
  const Vector             mu = streams();
  const auto               F  = legendre_matrix(series());
  const auto               d  = fourier_modes(F, mu, mu, NF, 64);
  const auto               b  = fourier_modes(F, mu, Vector{-mu0}, NF, 64);
  rtepack::muelmat_tensor4 Cd(NF, NL, NQ, NQ), Sd(NF, NL, NQ, NQ);
  rtepack::muelmat_tensor3 Cb(NF, NL, NQ), Sb(NF, NL, NQ);
  for (Index m = 0; m < NF; m++)
    for (Index l = 0; l < NL; l++)
      for (Index i = 0; i < NQ; i++) {
        for (Index j = 0; j < NQ; j++) {
          Cd[m, l, i, j] = d.C[m, i, j];
          Sd[m, l, i, j] = d.S[m, i, j];
        }
        Cb[m, l, i] = b.C[m, i, 0];
        Sb[m, l, i] = b.S[m, i, 0];
      }
  return {.diffuse = vdisort::combine_phase_matrices(Cd, Sd), .beam = vdisort::combine_beam_phase_matrices(Cb, Sb)};
}

rtepack::stokvec_tensor3 boundary(Numeric sign) {
  rtepack::stokvec_tensor3 b(2, NF, NQ / 2);
  for (Index i = 0; i < NQ / 2; i++) {
    b[0, 0, i] = rtepack::stokvec{0.5 * (1 + sign * 0.1 * static_cast<Numeric>(i)), 0.01, 0.0, 0.0};
    b[1, 1, i] = rtepack::stokvec{0.0, 0.0, 0.003 * sign, 0.001};
  }
  return b;
}

Tensor4 radiance(vdisort::phase_matrix_data diffuse) {
  const AscendingGrid tau{Vector{0.3, 0.8, 1.6}};
  Vector              streams(NQ / 2), weights(NQ / 2);
  Legendre::PositiveDoubleGaussLegendre(streams, weights);
  const auto          v = vdisort::main_data(NQ,
                                             NF,
                                             tau,
                                             Vector{0.5, 0.95, 0.7},
                                             std::move(diffuse),
                                             boundary(1),
                                             boundary(-1),
                                             rtepack::stokvec_matrix(NL, 0),
                                             vdisort::brdf::lambertian_fourier_modes(
                                                 0.3, 1, streams, streams, std::numeric_limits<Numeric>::quiet_NaN()),
                                             0.6,
                                             rtepack::stokvec{0.0, 0.0, 0.0, 0.0},
                                             0.3);
  const AscendingGrid levels{Vector{0.0, 0.15, 0.3, 0.8, 1.2, 1.6}};
  const Vector        phi{0.0, 0.7, 2.5};
  Tensor4             out(levels.size(), phi.size(), NQ, 4);
  v.ungridded_u(out, levels, phi);
  return out;
}
}  // namespace

int main() try {
  auto       symmetric = phase_matrices(0.6).diffuse;
  auto       changed   = symmetric;
  const auto N         = NQ / 2;
  for (Index alpha = 0; alpha < 2; alpha++)
    for (Index m = 0; m < NF; m++) {
      if (alpha == vdisort::cosine_mode and m == 0) continue;
      for (Index l = 0; l < NL; l++)
        for (Index i = 0; i < N; i++)
          for (Index j = 0; j < NQ; j++) changed[alpha, m, l, N + i, j] *= 1.0 + 1e-9;
    }

  const Tensor4 a = radiance(symmetric), b = radiance(changed);
  Numeric       scale = 0.0, difference = 0.0;
  for (Size n = 0; n < a.size(); n++) {
    scale      = std::max(scale, std::abs(a.data_handle()[n]));
    difference = std::max(difference, std::abs(a.data_handle()[n] - b.data_handle()[n]));
  }
  std::cout << std::format("downward rows changed by 1e-9: radiance changed by {:.3e} of max\n", difference / scale);
  require(difference > 0.0, "An asymmetric phase matrix must not take the half-size eigenproblem");
  require(difference < 1e-8 * scale, "The half-size and full eigensolutions must agree");
  std::cout << "vdisort mirror test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
