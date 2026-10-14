/* One VDISORT solver reused for several problems, as a thread does for many
   frequencies, against a fresh solver for each problem.

   update_all() allocates nothing and resets only what is accumulated or
   written conditionally, so every reused result must equal the fresh one bit
   for bit.  The problems change every input through the solver's accessors
   and exercise each conditional path: a layer without scattering (no
   scattered beam, emitting), an exactly conservative layer (the Jordan pair,
   no emission), a polarized beam, no beam, and a surface in the first
   Fourier mode only.  Without the reset of the beam terms of a layer without
   scattering, the reused radiances miss by 0.14. */

#include <legendre.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <bit>
#include <cstdint>
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

constexpr Index NQ = 8, NF = 4, NL = 3, NS = 2;

struct problem {
  Vector           tau, omega;
  rtepack::stokvec beam;
  Numeric          mu0, phi0, albedo, source_scale, boundary;
};

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

rtepack::stokvec_matrix sources(const problem& p) {
  rtepack::stokvec_matrix s(NL, NS);
  for (Index l = 0; l < NL; l++) {
    s[l, 0] = rtepack::stokvec{p.source_scale * static_cast<Numeric>(1 + l), 0.01 * static_cast<Numeric>(l), 0.0, 0.0};
    s[l, 1] = rtepack::stokvec{0.2, 0.0, 0.0, 0.0};
  }
  return s;
}

rtepack::stokvec_tensor3 boundary(const problem& p, Numeric sign) {
  rtepack::stokvec_tensor3 b(2, NF, NQ / 2);
  for (Index i = 0; i < NQ / 2; i++) {
    b[0, 0, i] = rtepack::stokvec{p.boundary * (1 + sign * 0.1 * static_cast<Numeric>(i)), 0.01, 0.0, 0.0};
    b[1, 1, i] = rtepack::stokvec{0.0, 0.0, 0.003 * sign, 0.001};
  }
  return b;
}

//! The Lambertian surface on the streams, for the problem's beam
vdisort::surface_reflection surface(const problem& p) {
  Vector streams(NQ / 2), weights(NQ / 2);
  Legendre::PositiveDoubleGaussLegendre(streams, weights);
  const Numeric mu0 = p.beam.I() > 0.0 ? p.mu0 : std::numeric_limits<Numeric>::quiet_NaN();
  return vdisort::brdf::lambertian_fourier_modes(p.albedo, 1, streams, streams, mu0);
}

vdisort::main_data fresh(const problem& p) {
  auto ph = phase_matrices(p.mu0);
  return vdisort::main_data(NQ,
                            NF,
                            AscendingGrid{p.tau},
                            p.omega,
                            std::move(ph.diffuse),
                            boundary(p, 1),
                            boundary(p, -1),
                            sources(p),
                            surface(p),
                            p.mu0,
                            p.beam,
                            p.phi0,
                            std::move(ph.beam));
}

void reuse(vdisort::main_data& v, const problem& p) {
  const auto ph = phase_matrices(p.mu0);
  v.tau(p.tau);
  v.omega()                   = p.omega;
  v.all_phase_matrices()      = ph.diffuse;
  v.all_beam_phase_matrices() = ph.beam;
  v.upward_boundary()         = boundary(p, 1);
  v.downward_boundary()       = boundary(p, -1);
  v.source_poly()             = sources(p);
  v.surface(surface(p));
  v.solar_zenith()            = p.mu0;
  v.beam_azimuth()            = p.phi0;
  v.set_beam_source(p.beam);
  v.update_all();
}

bool identical(const Numeric* a, const Numeric* b, Size n) {
  for (Size i = 0; i < n; i++)
    if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
  return true;
}
}  // namespace

int main() try {
  const problem A{.tau          = Vector{0.3, 0.8, 1.6},
                  .omega        = Vector{0.5, 0.95, 0.7},
                  .beam         = {1.0, 0.1, 0.0, 0.0},
                  .mu0          = 0.6,
                  .phi0         = 0.3,
                  .albedo       = 0.3,
                  .source_scale = 1.0,
                  .boundary     = 0.5};
  const problem B{.tau          = Vector{0.2, 0.5, 2.0},
                  .omega        = Vector{0.0, 1.0, 0.6},
                  .beam         = {1.0, -0.2, 0.1, 0.05},
                  .mu0          = 0.45,
                  .phi0         = 1.1,
                  .albedo       = 0.1,
                  .source_scale = 0.7,
                  .boundary     = 0.2};
  const problem C = [&B] {
    problem c = B;
    c.beam    = rtepack::stokvec{0.0, 0.0, 0.0, 0.0};
    return c;
  }();

  vdisort::main_data v = fresh(A);
  for (const problem* p : {&B, &C, &A, &B}) {
    reuse(v, *p);
    const auto          w      = fresh(*p);
    const Numeric       bottom = p->tau[NL - 1];
    const AscendingGrid tau{Vector{0.0, 0.05 * bottom, 0.1 * bottom, 0.3 * bottom, 0.6 * bottom, 0.9 * bottom, bottom}};
    const Vector        phi{0.0, 0.7, 2.5};
    Tensor4             a(tau.size(), phi.size(), NQ, 4), b(tau.size(), phi.size(), NQ, 4);
    v.ungridded_u(a, tau, phi);
    w.ungridded_u(b, tau, phi);
    Matrix fa(4, tau.size()), fb(4, tau.size());
    v.ungridded_flux(fa[0], fa[1], fa[2], fa[3], tau);
    w.ungridded_flux(fb[0], fb[1], fb[2], fb[3], tau);
    require(identical(a.data_handle(), b.data_handle(), a.size()) and
                identical(fa.data_handle(), fb.data_handle(), fa.size()),
            "A reused solver must give the radiances and fluxes of a fresh one");
  }
  std::cout << "vdisort reuse test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
