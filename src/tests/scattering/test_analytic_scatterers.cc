/** Bulk scattering properties of the analytic species against closed forms.
 *
 * The Henyey-Greenstein scattering matrix (F11 = F22 = p, F33 = F44 =
 * p (3 cos(Theta) - cos^3(Theta)) / 2) and Chandrasekhar's Rayleigh scattering matrix are written
 * out here, independently of the species code; the spherical-harmonics
 * coefficients are checked by summing them with Legendre polynomials from
 * their recurrence.  The lab-frame (ARO) phase matrix is checked through quantities that
 * do not depend on the rotation into the lab frame: with Theta the angle
 * between the incident and scattered directions, Z11 = F11(Theta),
 * Z44 = F44(Theta), Z12^2 + Z13^2 = F12^2 and Z21^2 + Z31^2 = F12^2.  The
 * directions are chosen away from the tabulation points of any grid.
 */
#include <arts_conversions.h>
#include <legendre.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include "gas_scattering.h"
#include "henyey_greenstein.h"

namespace {
using namespace scattering;

bool close(Numeric a, Numeric b, Numeric relative) {
  return std::abs(a - b) <= relative * std::max({Numeric{1e-300}, std::abs(a), std::abs(b)});
}

/** Henyey-Greenstein phase function per steradian, normalised to 1 over the sphere */
Numeric henyey_greenstein(Numeric g, Numeric cos_theta) {
  return (1.0 - g * g) / (4.0 * Constant::pi * std::pow(1.0 + g * g - 2.0 * g * cos_theta, 1.5));
}

/** F33 / F11 = F44 / F11 of the Henyey-Greenstein scattering matrix */
Numeric hg_ratio(Numeric c) { return (3.0 * c - c * c * c) / 2.0; }

/** Rayleigh [F11, F12, F44] (Chandrasekhar, 1950), F11 normalised to 4 pi over the sphere */
std::array<Numeric, 3> rayleigh(Numeric cos_theta) {
  return {0.75 * (1.0 + cos_theta * cos_theta), -0.75 * (1.0 - cos_theta * cos_theta), 1.5 * cos_theta};
}

Numeric cos_scattering_angle(Numeric za_inc, Numeric delta_aa, Numeric za_scat) {
  using Conversion::cosd, Conversion::sind;
  return cosd(za_inc) * cosd(za_scat) + sind(za_inc) * sind(za_scat) * cosd(delta_aa);
}

const Vector za_inc{12.3, 77.7, 140.1};
const Vector delta_aa{21.4, 47.3, 133.9, 251.9};
const Vector za_scat{33.3, 61.1, 101.9, 166.2};

std::shared_ptr<ZenithAngleGrid> za_scat_grid() {
  return std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat));
}

/** Every frequency of a multi-frequency TRO request has its own scattering coefficient and asymmetry parameter */
bool test_henyey_greenstein_frequencies() {
  const Vector                    f_grid{1e9, 2e9, 3e9};
  const auto                      ext = [](Numeric f) { return 1e-3 * f / 1e9; };
  const auto                      ssa = [](Numeric f) { return 0.5 + 0.1 * f / 1e9; };
  const auto                      g   = [](Numeric f) { return 0.8 - 0.2 * f / 1e9; };
  const HenyeyGreensteinScatterer hg{"extinction", "ssa", "g"};
  Vector                          ext_grid(f_grid.size()), ssa_grid(f_grid.size()), g_grid(f_grid.size());
  for (Size iv = 0; iv < f_grid.size(); iv++) {
    ext_grid[iv] = ext(f_grid[iv]);
    ssa_grid[iv] = ssa(f_grid[iv]);
    g_grid[iv]   = g(f_grid[iv]);
  }

  // Gauss-Legendre scattering angles, ascending
  const Index n = 100;
  Vector      x(n), w(n), angles(n);
  Legendre::GaussLegendre(x, w);
  for (Index i = 0; i < n; i++) angles[i] = Conversion::acosd(x[n - 1 - i]);

  const auto bulk = hg.tro_gridded(
      ext_grid, ssa_grid, g_grid, f_grid, std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(angles)));
  const auto& pm = *bulk.phase_matrix;
  for (Index iv = 0; iv < 3; iv++) {
    const Numeric scattering = ext(f_grid[iv]) * ssa(f_grid[iv]);
    if (not close(bulk.extinction_matrix[0, iv].A(), ext(f_grid[iv]), 1e-14)) return false;
    if (not close(bulk.absorption_vector[0, iv].I(), ext(f_grid[iv]) - scattering, 1e-14)) return false;

    Numeric integral = 0.0;
    for (Index i = 0; i < n; i++) {
      const Numeric c = x[n - 1 - i];
      const Numeric p = scattering * henyey_greenstein(g(f_grid[iv]), c);
      for (auto [is, ref] :
           {std::pair{0, p}, std::pair{2, p}, std::pair{3, p * hg_ratio(c)}, std::pair{5, p * hg_ratio(c)}}) {
        if (not close(pm[0, iv, i].data[is], ref, 1e-12)) {
          std::cout << "f = " << f_grid[iv] << " Hz, F" << is << ": " << pm[0, iv, i].data[is] << " != " << ref << '\n';
          return false;
        }
      }
      if (pm[0, iv, i].F12() != 0.0 or pm[0, iv, i].F34() != 0.0) return false;
      integral += 2.0 * Constant::pi * w[n - 1 - i] * pm[0, iv, i].F11();
    }
    if (not close(integral, scattering, 1e-12)) {
      std::cout << "f = " << f_grid[iv] << " Hz: the phase function integrates to " << integral << ", not "
                << scattering << '\n';
      return false;
    }
  }
  return true;
}

/** Checks Z11 = s * F11(Theta), Z44 = s * F44(Theta) and |F12| invariants for every direction pair */
template <typename ScatteringMatrix>
bool check_lab_frame(const BulkScatteringProperties<Format::ARO, Representation::Gridded>& bulk,
                     Numeric                                                               scale,
                     ScatteringMatrix&&                                                    f,
                     const char*                                                           name) {
  const auto& pm = *bulk.phase_matrix;
  for (Size i = 0; i < za_inc.size(); i++) {
    for (Size k = 0; k < delta_aa.size(); k++) {
      for (Size j = 0; j < za_scat.size(); j++) {
        const auto [f11, f12, f44] = f(cos_scattering_angle(za_inc[i], delta_aa[k], za_scat[j]));
        const auto z               = [&](Index row, Index col) { return pm[0, 0, i, k, j][row, col]; };
        const bool ok              = close(z(0, 0), scale * f11, 1e-12) and close(z(3, 3), scale * f44, 1e-12) and
                                     close(std::hypot(z(0, 1), z(0, 2)), scale * std::abs(f12), 1e-12) and
                                     close(std::hypot(z(1, 0), z(2, 0)), scale * std::abs(f12), 1e-12);
        if (not ok) {
          std::cout << name << " at za_inc = " << za_inc[i] << ", delta_aa = " << delta_aa[k]
                    << ", za_scat = " << za_scat[j] << ": Z11 = " << z(0, 0) << " != " << scale * f11 << '\n';
          return false;
        }
      }
    }
  }
  return true;
}

bool test_henyey_greenstein_lab_frame() {
  // Strongly forward peaked, where tabulating Theta is least accurate
  const Numeric                   g = 0.9;
  const HenyeyGreensteinScatterer hg{"extinction", "ssa", "g"};
  AtmSpectralPoint                spectral;
  spectral[AtmSpectralPropertyTag{"extinction"}] = Vector{2e-3};
  spectral[AtmSpectralPropertyTag{"ssa"}]        = Vector{0.7};
  spectral[AtmSpectralPropertyTag{"g"}]          = Vector{g};
  const auto                      f = [g](Numeric c) {
    return std::array{henyey_greenstein(g, c), 0.0, hg_ratio(c) * henyey_greenstein(g, c)};
  };

  const auto bulk = hg.aro_gridded(AtmPoint{}, spectral, Vector{1e9}, za_inc, delta_aa, za_scat_grid());
  if (not check_lab_frame(bulk, 2e-3 * 0.7, f, "HG")) return false;
  for (Size i = 0; i < za_inc.size(); i++) {
    if (not close(bulk.extinction_matrix[0, 0, i].A(), 2e-3, 1e-14)) return false;
    if (not close(bulk.absorption_vector[0, 0, i].I(), 2e-3 * 0.3, 1e-14)) return false;
  }

  const Vector extinction{2e-3}, ssa{0.7}, gs{g}, one{1.0}, zero{0.0};
  const auto   d_ext =
      hg.aro_gridded_derivative(extinction, ssa, gs, one, zero, zero, Vector{1e9}, za_inc, delta_aa, za_scat_grid());
  if (not check_lab_frame(d_ext, 0.7, f, "dHG/dextinction")) return false;

  const auto d_ssa =
      hg.aro_gridded_derivative(extinction, ssa, gs, zero, one, zero, Vector{1e9}, za_inc, delta_aa, za_scat_grid());
  return check_lab_frame(d_ssa, 2e-3, f, "dHG/dssa");
}

/** The derivative by g against central finite differences, and by the spectral properties through the wrappers
 *
 * Two frequencies have different g, so the derivative at each frequency must
 * be by the g of that frequency alone.
 */
bool test_henyey_greenstein_g_derivative() {
  const HenyeyGreensteinScatterer hg{"extinction", "ssa", "g"};
  const Vector                    f_grid{1e9, 2e9}, extinction{2e-3, 1e-3}, ssa{0.7, 0.4}, g{0.6, -0.3};
  const Vector                    zero(2, 0.0), one(2, 1.0);
  constexpr Numeric               h = 1e-6;

  const auto with_g = [&](Numeric dg) {
    Vector x  = g;
    x        += dg;
    return x;
  };

  // TRO: every element of the phase matrix
  const auto grid  = [] { return std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat)); };
  const auto d_tro = hg.tro_gridded_derivative(extinction, ssa, g, zero, zero, one, f_grid, grid());
  const auto plus  = hg.tro_gridded(extinction, ssa, with_g(h), f_grid, grid());
  const auto minus = hg.tro_gridded(extinction, ssa, with_g(-h), f_grid, grid());
  for (Size iv = 0; iv < f_grid.size(); iv++) {
    for (Size ia = 0; ia < za_scat.size(); ia++) {
      for (Index e = 0; e < 6; e++) {
        const Numeric fd =
            ((*plus.phase_matrix)[0, iv, ia].data[e] - (*minus.phase_matrix)[0, iv, ia].data[e]) / (2 * h);
        const Numeric an = (*d_tro.phase_matrix)[0, iv, ia].data[e];
        if (std::abs(an - fd) > 1e-7 * std::max(1e-6, std::abs(fd))) {
          std::cout << "TRO dHG/dg at f = " << f_grid[iv] << " Hz, za = " << za_scat[ia] << ", element " << e << ": "
                    << an << " != " << fd << '\n';
          return false;
        }
      }
    }
    if (d_tro.extinction_matrix[0, iv].A() != 0.0 or d_tro.absorption_vector[0, iv].I() != 0.0) return false;
  }

  // ARO: every element of the laboratory-frame phase matrix
  const auto d_aro =
      hg.aro_gridded_derivative(extinction, ssa, g, zero, zero, one, f_grid, za_inc, delta_aa, za_scat_grid());
  const auto aro_plus  = hg.aro_gridded(extinction, ssa, with_g(h), f_grid, za_inc, delta_aa, za_scat_grid());
  const auto aro_minus = hg.aro_gridded(extinction, ssa, with_g(-h), f_grid, za_inc, delta_aa, za_scat_grid());
  Numeric    worst     = 0.0;
  for (Size iv = 0; iv < f_grid.size(); iv++) {
    for (Size i = 0; i < za_inc.size(); i++) {
      for (Size k = 0; k < delta_aa.size(); k++) {
        for (Size j = 0; j < za_scat.size(); j++) {
          for (Index e = 0; e < 16; e++) {
            const Numeric fd =
                ((*aro_plus.phase_matrix)[0, iv, i, k, j].data[e] - (*aro_minus.phase_matrix)[0, iv, i, k, j].data[e]) /
                (2 * h);
            const Numeric an = (*d_aro.phase_matrix)[0, iv, i, k, j].data[e];
            worst            = std::max(worst, std::abs(an - fd) / std::max(1e-6, std::abs(fd)));
          }
        }
      }
    }
  }
  if (worst > 1e-7) {
    std::cout << "ARO dHG/dg deviates from finite differences by " << worst << '\n';
    return false;
  }

  // The wrappers select the property by its name, and a property the scatterer does not read gives zero
  AtmSpectralPoint spectral;
  spectral[AtmSpectralPropertyTag{"extinction"}] = extinction;
  spectral[AtmSpectralPropertyTag{"ssa"}]        = ssa;
  spectral[AtmSpectralPropertyTag{"g"}]          = g;
  const auto by_g                                = hg.aro_gridded_derivative(
      AtmPoint{}, spectral, f_grid, za_inc, delta_aa, za_scat_grid(), AtmSpectralPropertyTag{"g"});
  const auto by_other = hg.aro_gridded_derivative(
      AtmPoint{}, spectral, f_grid, za_inc, delta_aa, za_scat_grid(), AtmSpectralPropertyTag{"refractive index"});
  for (Size iv = 0; iv < f_grid.size(); iv++) {
    for (Index e = 0; e < 16; e++) {
      if ((*by_g.phase_matrix)[0, iv, 1, 1, 1].data[e] != (*d_aro.phase_matrix)[0, iv, 1, 1, 1].data[e]) return false;
      if ((*by_other.phase_matrix)[0, iv, 1, 1, 1].data[e] != 0.0) return false;
    }
  }
  return true;
}

/** The spherical-harmonics coefficients sum to the closed form: F11 = F22 = p and F33 = F44 = p f */
bool test_henyey_greenstein_spectral() {
  const Numeric                   g = 0.5, scattering = 3e-4;
  const Index                     degree = 80;  // g^degree is far below round-off
  const HenyeyGreensteinScatterer hg{"extinction", "ssa", "g"};
  const auto  spectral = hg.tro_spectral(Vector{2.0 * scattering}, Vector{0.5}, Vector{g}, Vector{1e9}, degree);
  const auto& pm       = *spectral.phase_matrix;
  for (Numeric c : {-1.0, -0.73, -0.2, 0.0, 0.31, 0.88, 1.0}) {
    std::array<Numeric, 4> sum{};
    Numeric                p_prev = 0.0, p_k = 1.0;  // P_{k-1}(c), P_k(c)
    for (Index k = 0; k <= degree; k++) {
      const Numeric y = std::sqrt((2.0 * static_cast<Numeric>(k) + 1.0) / (4.0 * Constant::pi)) * p_k;
      for (Index i = 0; i < 4; i++) sum[i] += pm[0, k].expand()[i, i] * y;
      const Numeric p_next = ((2.0 * static_cast<Numeric>(k) + 1.0) * c * p_k - static_cast<Numeric>(k) * p_prev) /
                             static_cast<Numeric>(k + 1);
      p_prev               = p_k;
      p_k                  = p_next;
    }
    const Numeric p = scattering * henyey_greenstein(g, c), pf = p * hg_ratio(c);
    // F33 and F44 vanish at 90 deg, so they are compared relative to p
    if (not(close(sum[0], p, 1e-12) and close(sum[1], p, 1e-12) and std::abs(sum[2] - pf) <= 1e-12 * p and
            std::abs(sum[3] - pf) <= 1e-12 * p)) {
      std::cout << "cos(Theta) = " << c << ": " << sum[0] << ", " << sum[1] << ", " << sum[2] << ", " << sum[3]
                << " != " << p << ", " << p << ", " << pf << ", " << pf << '\n';
      return false;
    }
  }
  return true;
}

/** The laboratory-frame matrix is continuous at backscattering
 *
 * Approaching the backscattering direction from different azimuths must
 * give the same limit, the matrix at exact backscattering, so the deviation
 * must shrink in proportion to the distance.  This holds only for a matrix
 * with F33 = -F22 at 180 deg; F33 = F22 everywhere gives an O(1) deviation at
 * any distance.
 */
bool test_henyey_greenstein_backscatter() {
  const Numeric                   g = 0.5;
  const HenyeyGreensteinScatterer hg{"extinction", "ssa", "g"};
  const Numeric                   za_in = 40.0, za_back = 180.0 - za_in;

  const auto lab = [&](Numeric delta_aa, Numeric za) {
    const auto bulk = hg.aro_gridded(Vector{1.0},
                                     Vector{1.0},
                                     Vector{g},
                                     Vector{1e9},
                                     Vector{za_in},
                                     Vector{delta_aa},
                                     std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(Vector{za})));
    Vector     z(16);
    for (Index e = 0; e < 16; e++) z[e] = (*bulk.phase_matrix)[0, 0, 0, 0, 0].data[e];
    return z;
  };

  const Vector  back  = lab(180.0, za_back);
  const Numeric scale = henyey_greenstein(g, -1.0);
  for (Numeric eps : {0.8, 0.2}) {
    Numeric worst = 0.0;
    for (Numeric psi : {30.0, 120.0, 210.0, 300.0}) {
      const Vector z =
          lab(180.0 + eps * Conversion::sind(psi) / Conversion::sind(za_back), za_back + eps * Conversion::cosd(psi));
      for (Index e = 0; e < 16; e++) worst = std::max(worst, std::abs(z[e] - back[e]) / scale);
    }
    std::cout << "Henyey-Greenstein, " << eps
              << " deg from backscattering: max |Z - Z(180 deg)| / p(180 deg) = " << worst << '\n';
    if (worst > 0.05 * eps) return false;
  }
  return close(back[0], scale, 1e-12) and close(back[5], scale, 1e-12) and close(back[10], -scale, 1e-12) and
         close(back[15], -scale, 1e-12);
}

bool test_rayleigh_lab_frame() {
  const AtmPoint     point{101325.0, 288.15};
  const GasScatterer gas{ConstantGasScattering{4.65e-31}, RayleighGasScattering{0.0}};

  const auto    bulk = gas.aro_gridded(point, AtmSpectralPoint{}, Vector{6e14}, za_inc, delta_aa, za_scat_grid());
  const Numeric k    = bulk.extinction_matrix[0, 0, 0].A();
  if (not check_lab_frame(bulk, k / (4.0 * Constant::pi), rayleigh, "Rayleigh")) return false;

  const auto d_p = gas.aro_gridded_derivative(
      point, AtmSpectralPoint{}, Vector{6e14}, za_inc, delta_aa, za_scat_grid(), AtmKeyVal{AtmKey::p});
  return check_lab_frame(d_p, k / (4.0 * Constant::pi * point.pressure), rayleigh, "dRayleigh/dp");
}
}  // namespace

int main() {
  if (not test_henyey_greenstein_frequencies()) {
    std::cerr << "Henyey-Greenstein TRO data at several frequencies failed\n";
    return 1;
  }
  if (not test_henyey_greenstein_lab_frame()) {
    std::cerr << "Henyey-Greenstein lab-frame data failed\n";
    return 1;
  }
  if (not test_henyey_greenstein_g_derivative()) {
    std::cerr << "Henyey-Greenstein derivative by g failed\n";
    return 1;
  }
  if (not test_henyey_greenstein_spectral()) {
    std::cerr << "Henyey-Greenstein spectral data failed\n";
    return 1;
  }
  if (not test_henyey_greenstein_backscatter()) {
    std::cerr << "Henyey-Greenstein lab-frame data at backscattering failed\n";
    return 1;
  }
  if (not test_rayleigh_lab_frame()) {
    std::cerr << "Rayleigh lab-frame data failed\n";
    return 1;
  }
  return 0;
}
