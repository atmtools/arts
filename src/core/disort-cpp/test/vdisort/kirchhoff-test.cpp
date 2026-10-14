/* Kirchhoff's law in DISORT and VDISORT.

   An isothermal atmosphere over a surface at the same temperature, under a sky
   at it too, holds the blackbody radiance B e_I everywhere, e_I = [1, 0, 0, 0].
   Each solver is given the Planck radiance beneath the surface as its upward
   boundary and transmits it with 1 - R, R its reflection on the streams
   (kernel and specular part): for B e_I that is the emission B (e_I - R e_I).
   The discrete equations then balance exactly, so on every stream and at
   every user angle the radiance must be B e_I to rounding, whatever the
   surface reflects, polarizing or not.  The layers scatter (Henyey-Greenstein
   in DISORT, Rayleigh with polarization in VDISORT) and emit (1 - omega) B.

   Negative control: given the emission (e_I - R e_I) B as its boundary instead,
   the convention before the solvers applied Kirchhoff's law, the surface emits
   it once more reduced, and the radiance misses B by about the surface's
   albedo times its emissivity. */

#include <disort-brdf.h>
#include <disort.h>
#include <legendre.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

#include "lab-frame.h"

namespace {
void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

constexpr Index   NQ = 8, N = NQ / 2, NF = 4, NL = 3;
constexpr Numeric B  = 3.7;  // any radiance: the cavity is linear in it

const Vector        omega{0.5, 0.95, 0.3};
const AscendingGrid tau{0.2, 0.7, 1.5};
const AscendingGrid output_tau{0.0, 0.1, 0.2, 0.45, 0.7, 1.1, 1.5};
const Vector        phi{0.0, 0.8, 2.6, 4.4};
const Vector        user_mu{0.9, 0.33, 0.05, -0.05, -0.33, -0.9};  // upward with their downward partners
const Vector        user_abs_mu{0.9, 0.33, 0.05, 0.05, 0.33, 0.9};
constexpr Numeric   no_beam = std::numeric_limits<Numeric>::quiet_NaN();

//! The cosines of the upward streams
Vector upward_streams() {
  Vector mu(N), w(N);
  Legendre::PositiveDoubleGaussLegendre(mu, w);
  return mu;
}

//! A surface's reflection from the streams into the given upward directions
template <typename surface_reflection>
using surface_model = std::function<surface_reflection(const ConstVectorView& mu_out)>;

//! The largest |radiance / B - 1| of the I component and |Q|, |U|, |V| / B
struct deviation {
  Numeric I{}, polarized{};
};

void report(std::string_view what, deviation d, Numeric tolerance) {
  std::cout << std::format("{:<52} I {:9.2e}  Q, U, V {:9.2e}\n", what, d.I, d.polarized);
  require(d.I <= tolerance and d.polarized <= tolerance,
          std::format("{}: an isothermal cavity must hold B e_I to {:.0e}, got {:.2e} in I and {:.2e} in Q, U, V",
                      what,
                      tolerance,
                      d.I,
                      d.polarized));
}

void report_miss(std::string_view what, deviation d, Numeric least) {
  std::cout << std::format("{:<52} I {:9.2e}  (negative control)\n", what, d.I);
  require(d.I >= least,
          std::format("{}: the emission as the boundary must miss B by at least {:.0e}, got {:.2e}", what, least, d.I));
}

/////////////////////////////////////////////////////////////////////////////
// DISORT
/////////////////////////////////////////////////////////////////////////////

disort::main_data disort_cavity(const disort::surface_reflection& surface, bool emission_as_boundary) {
  Matrix legendre(NL, NQ);
  for (Index l = 0; l < NL; l++)
    for (Index k = 0; k < NQ; k++) legendre[l, k] = std::pow(0.5, static_cast<Numeric>(k));
  Matrix source(NL, 1, B), up(NF, N, 0.0), down(NF, N, 0.0);
  down[0] = B;
  up[0]   = B;
  if (emission_as_boundary) {
    // The emission B (1 - sum_j R_ij) on the streams, R the solver's reflection of mode 0
    Vector mu(N), W(N);
    Legendre::PositiveDoubleGaussLegendre(mu, W);
    const auto D = surface.brdf()[0];
    for (Index i = 0; i < N; i++)
      for (Index j = 0; j < N; j++) up[0, i] -= 2.0 * D[i, j] * mu[j] * W[j] * B;
  }
  return disort::main_data(
      NQ, NQ, NF, tau, omega, legendre, up, down, Vector(NL, 0.0), source, surface, 0.5, 0.0, 0.0);
}

deviation disort_deviation(const disort::main_data& dis, const disort::surface_reflection& user_surface) {
  deviation d;
  Tensor3   streams(output_tau.size(), phi.size(), NQ);
  dis.ungridded_u(streams, output_tau, phi);
  for (const Numeric x : streams | by_elem) d.I = std::max(d.I, std::abs(x / B - 1.0));

  disort::user_u_data user;
  for (const Numeric t : output_tau)
    for (const Numeric p : phi) {
      dis.u_user(user, t, p, user_mu, user_surface);
      for (const Numeric x : user.intensities) d.I = std::max(d.I, std::abs(x / B - 1.0));
    }
  return d;
}

void test_disort() {
  using namespace disort::brdf;
  const Vector mu = upward_streams();
  const auto   raw = [&mu](const RawFunction& f) -> surface_model<disort::surface_reflection> {
    return [&mu, f](const ConstVectorView& out) { return fourier_modes(f, NF, out, mu, no_beam); };
  };
  const std::vector<std::pair<std::string_view, surface_model<disort::surface_reflection>>> surfaces{
      {"Lambertian A = 0.3",
       [&mu](const ConstVectorView& out) { return lambertian_fourier_modes(0.3, NF, out, mu, no_beam); }},
      {"Hapke (opposition hotspot)", raw(Hapke{})},
      {"Cox-Munk, shadowed", raw(CoxMunk{.shadowing = true})},
      {"RPV", raw(RPV{})},
      {"Ross-Li", raw(RossLi{})},
      {"Cox-Munk 0.7 / Lambertian A = 0.4",
       [&mu](const ConstVectorView& out) {
         return cox_munk_lambertian_fourier_modes(0.7, 0.4, 7.0, Complex{1.34, 0.0}, true, NF, out, mu, no_beam);
       }},
  };
  for (const auto& [name, surface] : surfaces)
    report(std::format("DISORT, {}", name),
           disort_deviation(disort_cavity(surface(mu), false), surface(user_abs_mu)),
           1e-12);

  const auto& lambertian = surfaces[0].second;
  report_miss("DISORT, Lambertian A = 0.3, emission as the boundary",
              disort_deviation(disort_cavity(lambertian(mu), true), lambertian(user_abs_mu)),
              0.1 * 0.3 * 0.7);
}

/////////////////////////////////////////////////////////////////////////////
// VDISORT
/////////////////////////////////////////////////////////////////////////////

//! Rayleigh scattering: F11 = 3/4 (1 + x^2), F12 = -3/4 (1 - x^2), F33 = F44 = 3/2 x
CompactPlanarMuelmatVector rayleigh() {
  CompactPlanarMuelmatVector c(3);
  c[0] = {1.0, -0.5, 1.0, 0.0, 0.0, 0.0};
  c[1] = {0.0, 0.0, 0.0, 1.5, 0.0, 1.5};
  c[2] = {0.5, 0.5, 0.5, 0.0, 0.0, 0.0};
  return c;
}

Vector streams() {
  Vector mu(NQ), w(N);
  Legendre::PositiveDoubleGaussLegendre(mu[Range{0, N}], w);
  for (Index i = 0; i < N; i++) mu[N + i] = -mu[i];
  return mu;
}

//! The combined phase matrices into the directions mu_out from the streams, the same in every layer
vdisort::phase_matrix_data phase_matrices(const Vector& mu_out) {
  const auto modes = vdisort_test::fourier_modes(vdisort_test::legendre_matrix(rayleigh()), mu_out, streams(), NF, 64);
  const Index              no    = static_cast<Index>(mu_out.size());
  rtepack::muelmat_tensor4 C(NF, NL, no, NQ), S(NF, NL, no, NQ);
  for (Index m = 0; m < NF; m++)
    for (Index l = 0; l < NL; l++)
      for (Index o = 0; o < no; o++)
        for (Index j = 0; j < NQ; j++) {
          C[m, l, o, j] = modes.C[m, o, j];
          S[m, l, o, j] = modes.S[m, o, j];
        }
  return vdisort::combine_phase_matrices(C, S);
}

vdisort::main_data vdisort_cavity(const vdisort::surface_reflection& surface, bool emission_as_boundary) {
  rtepack::stokvec_tensor3 up(2, NF, N), down(2, NF, N);
  up   = rtepack::stokvec{};
  down = rtepack::stokvec{};
  for (Index i = 0; i < N; i++) {
    up[vdisort::cosine_mode, 0, i]   = {B, 0.0, 0.0, 0.0};
    down[vdisort::cosine_mode, 0, i] = {B, 0.0, 0.0, 0.0};
  }
  if (emission_as_boundary) {
    // The emission B (e_I - sum_j R_ij e_I) on the streams, R the solver's reflection of mode 0
    Vector mu(N), W(N);
    Legendre::PositiveDoubleGaussLegendre(mu, W);
    const auto R = surface.brdf()[vdisort::cosine_mode, 0];
    for (Index i = 0; i < N; i++) {
      rtepack::stokvec reflected{};
      for (Index j = 0; j < N; j++) reflected += Constant::pi * W[j] * mu[j] * (R[i, j] * rtepack::stokvec{B, 0, 0, 0});
      if (surface.has_specular()) reflected += surface.specular()[i] * rtepack::stokvec{B, 0, 0, 0};
      up[vdisort::cosine_mode, 0, i] -= reflected;
    }
  }
  rtepack::stokvec_matrix source(NL, 1);
  source = rtepack::stokvec{B, 0.0, 0.0, 0.0};
  return vdisort::main_data(
      NQ, NF, tau, omega, phase_matrices(streams()), up, down, source, surface, 0.5, rtepack::stokvec{}, 0.0);
}

deviation vdisort_deviation(const vdisort::main_data& v, const vdisort::surface_reflection& user_surface) {
  deviation d;
  const auto add = [&d](const rtepack::stokvec& s) {
    d.I         = std::max(d.I, std::abs(s.I() / B - 1.0));
    d.polarized = std::max({d.polarized, std::abs(s.Q() / B), std::abs(s.U() / B), std::abs(s.V() / B)});
  };

  Tensor4 field(output_tau.size(), phi.size(), NQ, 4);
  v.ungridded_u(field, output_tau, phi);
  for (Index t = 0; t < field.extent(0); t++)
    for (Index p = 0; p < field.extent(1); p++)
      for (Index i = 0; i < NQ; i++)
        add(rtepack::stokvec{field[t, p, i, 0], field[t, p, i, 1], field[t, p, i, 2], field[t, p, i, 3]});

  rtepack::stokvec_tensor3 user(output_tau.size(), phi.size(), user_mu.size());
  v.ungridded_u_user(user, output_tau, phi, user_mu, phase_matrices(user_mu), {}, {}, user_surface);
  for (Index t = 0; t < user.extent(0); t++)
    for (Index p = 0; p < user.extent(1); p++)
      for (Index u = 0; u < user.extent(2); u++) add(user[t, p, u]);
  return d;
}

void test_vdisort() {
  using namespace vdisort::brdf;
  const Vector mu = upward_streams();
  const std::vector<std::pair<std::string_view, surface_model<vdisort::surface_reflection>>> surfaces{
      {"Lambertian A = 0.3",
       [&mu](const ConstVectorView& out) { return lambertian_fourier_modes(0.3, NF, out, mu, no_beam); }},
      {"Fresnel n = 3 + 0.2i (specular)",
       [&mu](const ConstVectorView& out) { return fresnel_fourier_modes(Complex{3.0, 0.2}, NF, out, mu, no_beam); }},
      {"Cox-Munk n = 1.34 + 0.01i, shadowed",
       [&mu](const ConstVectorView& out) {
         return cox_munk_fourier_modes(5.0, Complex{1.34, 0.01}, true, NF, out, mu, no_beam, 100);
       }},
      {"Fresnel 0.6 / Lambertian A = 0.3",
       [&mu](const ConstVectorView& out) {
         return fresnel_lambertian_fourier_modes(0.6, 0.3, Complex{2.0, 0.5}, NF, out, mu, no_beam);
       }},
      {"Hapke, depolarizing",
       [&mu](const ConstVectorView& out) { return hapke_fourier_modes(1.0, 0.06, 0.6, NF, out, mu, no_beam, 100); }},
  };
  for (const auto& [name, surface] : surfaces)
    report(std::format("VDISORT, {}", name),
           vdisort_deviation(vdisort_cavity(surface(mu), false), surface(user_abs_mu)),
           1e-12);

  const auto& fresnel = surfaces[1].second;
  report_miss("VDISORT, Fresnel n = 3 + 0.2i, emission as the boundary",
              vdisort_deviation(vdisort_cavity(fresnel(mu), true), fresnel(user_abs_mu)),
              1e-2);
}
}  // namespace

int main() try {
  test_disort();
  test_vdisort();
  std::cout << "Kirchhoff test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
