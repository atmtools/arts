/* Tests of rt4::scattering_optics and rt4::problem_from_path (rt4_arts.h).

   Where the expected answers come from:
   - The m = 0 azimuthal mean of Rayleigh scattering in RT4's [I, Q] basis
     is the closed form of Chandrasekhar (1960), written out below.  ARTS's
     GasScatterer gives the phase matrix through its own laboratory-frame
     rotation formulas (phase_matrix.h, Mishchenko-style spherical
     trigonometry), so the comparison tests the direction mapping (RT4's
     hemispheres to ARTS's propagation zenith angles), the in/out and Stokes
     order, the hemisphere quadrants and the normalisation.  GasScatterer
     gives the azimuthal mean exactly (its m = 0 Fourier mode at the
     streams), so the tolerance is round-off.
   - A3: strongly forward-peaked Henyey-Greenstein scattering (g up to
     0.99): the azimuthal mean of Z11 is, by the addition theorem,
     sigma sum_l (2 l + 1) / (4 pi) g^l P_l(mu_o) P_l(mu_i) with signed stream
     cosines, summed until g^l < 1e-17.
   - The path builder must reproduce the path data it is given (heights,
     temperatures, midpoint gas extinction) and the level optics of
     scattering_optics() averaged per layer. */
#include <arts_constants.h>
#include <legendre.h>
#include <physics_funcs.h>
#include <rt4_arts.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace rt4 = polradtran::rt4;

namespace {
constexpr Numeric pi = Constant::pi;

using rt4::down;
using rt4::up;

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

void require_error(const std::function<void()>& f, const std::string& what) {
  bool threw = false;
  try {
    f();
  } catch (const std::exception&) { threw = true; }
  require(threw, std::format("Expected an error: {}", what));
}

constexpr Numeric cross_section = 1e-30;  // m^2

ArrayOfScatteringSpecies rayleigh_species() {
  ArrayOfScatteringSpecies s;
  s.add(GasScatterer{scattering::ConstantGasScattering{cross_section}, scattering::RayleighGasScattering{0.0}});
  return s;
}

AtmPoint air(Numeric p, Numeric t) {
  AtmPoint a;
  a.pressure    = p;
  a.temperature = t;
  return a;
}

/* Rayleigh, m = 0, normalised to 1 over 4 pi, mo outgoing, mi incident (signed cosines), Q = I_v - I_h:
     P_II = 3/8 (3 - a - b + 3 a b),  P_IQ = 3/8 (1 - 3 a)(1 - b),
     P_QI = 3/8 (1 - a)(1 - 3 b),     P_QQ = 9/8 (1 - a)(1 - b),   a = mo^2, b = mi^2. */
Numeric rayleigh_m0(Index so, Index si, Numeric mo, Numeric mi) {
  const Numeric a = mo * mo, b = mi * mi;
  if (so == 0 and si == 0) return 3.0 / 8.0 * (3 - a - b + 3 * a * b);
  if (so == 0) return 3.0 / 8.0 * (1 - 3 * a) * (1 - b);
  if (si == 0) return 3.0 / 8.0 * (1 - a) * (1 - 3 * b);
  return 9.0 / 8.0 * (1 - a) * (1 - b);
}

Numeric signed_mu(Index h, Numeric mu) { return h == up ? mu : -mu; }

//! A1: GasScatterer Rayleigh through scattering_optics against the closed form
void test_rayleigh() {
  const auto species = rayleigh_species();
  const auto atm     = air(8e4, 260.0);
  const auto sigma   = cross_section * number_density(atm.pressure, atm.temperature);

  const auto q = polradtran::get_quadrature(8, polradtran::quadrature_type::double_gauss);
  Vector     mu(q.mu.size() + 2);
  std::ranges::copy(q.mu, mu.begin());
  mu[q.mu.size()]     = 0.35;
  mu[q.mu.size() + 1] = 1.0;
  const Index n       = static_cast<Index>(mu.size());

  const Numeric tol = 1e-13 * sigma / (4 * pi);

  {
    const auto o = rt4::scattering_optics(species, atm, 89e9, mu, 2);
    require(o.phase.shape() == (std::array<Index, 6>{2, 2, n, n, 2, 2}) and
                o.extinction.shape() == (std::array<Index, 4>{2, n, 2, 2}) and
                o.absorption.shape() == (std::array<Index, 3>{2, n, 2}),
            "A1: layer_optics shapes");

    Numeric dz = 0.0, dk = 0.0;
    for (Index ho = 0; ho < 2; ho++) {
      for (Index hi = 0; hi < 2; hi++) {
        for (Index io = 0; io < n; io++) {
          for (Index ii = 0; ii < n; ii++) {
            for (Index so = 0; so < 2; so++) {
              for (Index si = 0; si < 2; si++) {
                const Numeric ref =
                    sigma / (4 * pi) * rayleigh_m0(so, si, signed_mu(ho, mu[io]), signed_mu(hi, mu[ii]));
                dz = std::max(dz, std::abs(o.phase[ho, hi, io, ii, so, si] - ref));
              }
            }
          }
        }
      }
    }
    for (Index h0 = 0; h0 < 2; h0++) {
      for (Index i = 0; i < n; i++) {
        dk = std::max({dk,
                       std::abs(o.extinction[h0, i, 0, 0] - sigma),
                       std::abs(o.extinction[h0, i, 1, 1] - sigma),
                       std::abs(o.extinction[h0, i, 0, 1]),
                       std::abs(o.extinction[h0, i, 1, 0]),
                       std::abs(o.absorption[h0, i, 0]),
                       std::abs(o.absorption[h0, i, 1])});
      }
    }
    std::cout << std::format(
        "A1 Rayleigh GasScatterer: max |phase - closed form| / (sigma / 4 pi) {:.2e} (tolerance {:.2e}); K, a "
        "{:.1e}, mu = 1 included\n",
        dz / (sigma / (4 * pi)),
        tol / (sigma / (4 * pi)),
        dk / sigma);
    require(dz <= tol,
            std::format("A1: the RT4 phase quadrants of ARTS's Rayleigh GasScatterer must equal sigma / (4 pi) times "
                        "the m = 0 closed form to {:.2e}, got {:.2e}",
                        tol,
                        dz));
    require(dk <= 1e-14 * sigma, "A1: the extinction must be sigma on the diagonal and the absorption zero");
  }

  // nstokes 1 is the I element alone
  const auto o1 = rt4::scattering_optics(species, atm, 89e9, mu, 1);
  const auto o2 = rt4::scattering_optics(species, atm, 89e9, mu, 2);
  Numeric    d1 = 0.0;
  for (Index ho = 0; ho < 2; ho++)
    for (Index hi = 0; hi < 2; hi++)
      for (Index io = 0; io < n; io++)
        for (Index ii = 0; ii < n; ii++)
          d1 = std::max(d1, std::abs(o1.phase[ho, hi, io, ii, 0, 0] - o2.phase[ho, hi, io, ii, 0, 0]));
  require(d1 == 0.0, "A1: nstokes 1 must give the I element of nstokes 2");

  // No species: zero optics of the right shape
  const auto z = rt4::scattering_optics(ArrayOfScatteringSpecies{}, atm, 89e9, mu, 2);
  require(z.phase.shape() == (std::array<Index, 6>{2, 2, n, n, 2, 2}) and
              stdr::all_of(z.phase | by_elem, [](Numeric x) { return x == 0.0; }),
          "A1: an empty species array must give zero optics");

  require_error([&] { (void)rt4::scattering_optics(species, atm, 89e9, mu, 3); }, "nstokes 3");
  require_error([&] { (void)rt4::scattering_optics(species, atm, 89e9, Vector{0.5, 0.0}, 2); }, "mu = 0");
  require_error([&] { (void)rt4::scattering_optics(species, atm, -1.0, mu, 2); }, "negative frequency");
}

//! A3: forward-peaked Henyey-Greenstein against the addition theorem
void test_forward_peaked_hg() {
  const auto atm = air(8e4, 260.0);
  const auto q   = polradtran::get_quadrature(8, polradtran::quadrature_type::double_gauss);
  Vector     mu(q.mu.size() + 1);
  std::ranges::copy(q.mu, mu.begin());
  mu[q.mu.size()]     = 1.0;
  const Index   n     = static_cast<Index>(mu.size());
  const Numeric sigma = 0.9e-4;

  for (const Numeric g : {0.9, 0.95, 0.99}) {
    ArrayOfScatteringSpecies species;
    species.add(HenyeyGreensteinScatterer{
        ExtSSACallback{[](Numeric, const AtmPoint&) { return std::pair<Numeric, Numeric>{1e-4, 0.9}; }}, g});
    const auto o = rt4::scattering_optics(species, atm, 89e9, mu, 1);

    const auto L = static_cast<Index>(std::ceil(std::log(1e-17) / std::log(g)));
    Vector     po(L + 1), pi_(L + 1);
    Numeric    d = 0.0, scale = 0.0;
    for (Index ho = 0; ho < 2; ho++) {
      for (Index hi = 0; hi < 2; hi++) {
        for (Index io = 0; io < n; io++) {
          for (Index ii = 0; ii < n; ii++) {
            Legendre::legendre_polynomials(po, signed_mu(ho, mu[io]));
            Legendre::legendre_polynomials(pi_, signed_mu(hi, mu[ii]));
            Numeric ref = 0.0, gl = 1.0;
            for (Index l = 0; l <= L; l++, gl *= g)
              ref += (2.0 * static_cast<Numeric>(l) + 1.0) / (4 * pi) * gl * po[l] * pi_[l];
            ref   *= sigma;
            scale  = std::max(scale, std::abs(ref));
            d      = std::max(d, std::abs(o.phase[ho, hi, io, ii, 0, 0] - ref));
          }
        }
      }
    }
    std::cout << std::format(
        "A3 Henyey-Greenstein g = {}, 8 double-Gauss streams and mu = 1: max |phase_II - addition theorem| / max "
        "{:.1e} (tolerance 1e-10)\n",
        g,
        d / scale);
    require(d <= 1e-10 * scale,
            std::format("A3: RT4's azimuthal mean of Henyey-Greenstein scattering with g = {} must be that of the "
                        "addition theorem to 1e-10, got {:.2e}",
                        g,
                        d / scale));
  }
}

//! A path of nlev levels, top first, with gas extinction 1e-4 * (1 + level) per metre
struct path_data {
  ArrayOfPropagationPathPoint ray_path;
  ArrayOfAtmPoint             atm_path;
  ArrayOfPropmatVector        propmat;
  AscendingGrid               freq_grid{Vector{50e9, 89e9}};
};

path_data make_path(const Vector& altitude, const Vector& temperature) {
  path_data d;
  for (Index l = 0; l < static_cast<Index>(altitude.size()); l++) {
    PropagationPathPoint pp;
    pp.pos = {altitude[l], 0.0, 0.0};
    pp.los = {180.0, 0.0};
    d.ray_path.push_back(pp);
    d.atm_path.push_back(air(1e5 * std::exp(-altitude[l] / 8e3), temperature[l]));
    PropmatVector pm(2);
    pm[0].A() = 5e-5;
    pm[1].A() = 1e-4 * static_cast<Numeric>(1 + l);
    d.propmat.push_back(pm);
  }
  return d;
}

//! A2: the path builder against its inputs and scattering_optics
void test_path() {
  const auto species = rayleigh_species();
  const auto d       = make_path(Vector{3000.0, 2000.0, 800.0, 0.0}, Vector{230.0, 245.0, 270.0, 288.0});

  const rt4::path_settings s{.nstokes = 2, .nmu = 6, .extra_mu = Vector{1.0}};
  const auto               p = rt4::problem_from_path(d.ray_path,
                                                      d.atm_path,
                                                      d.propmat,
                                                      d.freq_grid,
                                                      1,
                                                      species,
                                                      s,
                                                      polradtran::lambertian_surface{.albedo = 0.2},
                                                      290.0,
                                                      2.7);

  require(p.frequency == 89e9 and p.nmu == 6 and p.extra_mu.size() == 1 and p.surface_temperature == 290.0 and
              p.sky_temperature == 2.7,
          "A2: settings and boundary values must be copied");
  Numeric dev = 0.0;
  for (Index l = 0; l < 4; l++)
    dev = std::max({dev,
                    std::abs(p.height[l] - d.ray_path[l].altitude()),
                    std::abs(p.temperature[l] - d.atm_path[l].temperature)});
  for (Index l = 0; l < 3; l++)
    dev = std::max(dev, std::abs(p.gas_extinction[l] - 1e-4 * (1.5 + static_cast<Numeric>(l))) / 1e-4);
  require(dev < 1e-14, std::format("A2: heights, temperatures and midpoint gas extinction, deviation {:.1e}", dev));
  require(p.optics.size() == 3 and p.layer_optics_index == ArrayOfIndex({0, 1, 2}),
          "A2: every layer with Rayleigh scattering must have its own optics set");

  // Layer 1 is the mean of levels 1 and 2
  const auto q = polradtran::get_quadrature(6, polradtran::quadrature_type::double_gauss);
  Vector     mu(7);
  std::ranges::copy(q.mu, mu.begin());
  mu[6]        = 1.0;
  const auto a = rt4::scattering_optics(species, d.atm_path[1], 89e9, mu, 2);
  const auto b = rt4::scattering_optics(species, d.atm_path[2], 89e9, mu, 2);
  Numeric    m = 0.0, scale = 0.0;
  for (Index x = 0; x < static_cast<Index>(a.phase.size()); x++) {
    const Numeric ref = 0.5 * (a.phase.data_handle()[x] + b.phase.data_handle()[x]);
    m                 = std::max(m, std::abs(p.optics[1].phase.data_handle()[x] - ref));
    scale             = std::max(scale, std::abs(ref));
  }
  Numeric mk = 0.0, sk = 0.0;
  for (Index x = 0; x < static_cast<Index>(a.extinction.size()); x++) {
    const Numeric ref = 0.5 * (a.extinction.data_handle()[x] + b.extinction.data_handle()[x]);
    mk                = std::max(mk, std::abs(p.optics[1].extinction.data_handle()[x] - ref));
    sk                = std::max(sk, std::abs(ref));
  }
  require(m <= 1e-15 * scale and mk <= 1e-15 * sk,
          std::format("A2: layer optics must be the mean of the level optics, {:.1e} and {:.1e}", m / scale, mk / sk));

  // Without species every layer is gas-only
  const auto g = rt4::problem_from_path(d.ray_path,
                                        d.atm_path,
                                        d.propmat,
                                        d.freq_grid,
                                        0,
                                        ArrayOfScatteringSpecies{},
                                        s,
                                        polradtran::lambertian_surface{},
                                        290.0,
                                        2.7);
  require(g.optics.empty() and g.layer_optics_index == ArrayOfIndex({-1, -1, -1}) and g.frequency == 50e9 and
              std::abs(g.gas_extinction[2] - 5e-5) < 1e-20,
          "A2: without scattering species every layer must be gas-only");

  const auto build = [&](const path_data& x, Index iv) {
    (void)rt4::problem_from_path(
        x.ray_path, x.atm_path, x.propmat, x.freq_grid, iv, species, s, polradtran::lambertian_surface{}, 290.0, 2.7);
  };
  auto e = d;
  require_error([&] { build(e, 2); }, "freq_index out of range");
  e.propmat[2][1].U() = 1e-6;
  require_error([&] { build(e, 1); }, "a polarized gas propagation matrix");
  e = d;
  std::swap(e.ray_path[1], e.ray_path[2]);
  require_error([&] { build(e, 1); }, "altitudes not decreasing");
  e = d;
  e.atm_path.pop_back();
  require_error([&] { build(e, 1); }, "atm_path of the wrong size");
  e            = d;
  e.propmat[0] = PropmatVector(1);
  require_error([&] { build(e, 0); }, "a propagation-matrix level of the wrong size");
  std::cout << "A2 path builder: heights, temperatures, midpoint gas extinction, layer means and 5 error paths\n";
}
}  // namespace

int main() try {
  test_rayleigh();
  test_forward_peaked_hg();
  test_path();
  std::cout << "rt4-arts test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}
