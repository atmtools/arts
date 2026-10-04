/** Bulk scattering properties of the analytic species against closed forms.
 *
 * The Henyey-Greenstein phase function and Chandrasekhar's Rayleigh
 * scattering matrix are written out here, independently of the species
 * code.  The lab-frame (ARO) phase matrix is checked through quantities that
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

/** Every frequency of a multi-frequency TRO request is scaled by its own scattering coefficient */
bool test_henyey_greenstein_frequencies() {
  const Numeric                   g = 0.6;
  const Vector                    f_grid{1e9, 2e9, 3e9};
  const auto                      ext = [](Numeric f) { return 1e-3 * f / 1e9; };
  const auto                      ssa = [](Numeric f) { return 0.5 + 0.1 * f / 1e9; };
  const HenyeyGreensteinScatterer hg{
      ExtSSACallback{[&](Numeric f, const AtmPoint&) { return std::pair{ext(f), ssa(f)}; }}, g};

  // Gauss-Legendre scattering angles, ascending
  const Index n = 100;
  Vector      x(n), w(n), angles(n);
  Legendre::GaussLegendre(x, w);
  for (Index i = 0; i < n; i++) angles[i] = Conversion::acosd(x[n - 1 - i]);

  const auto bulk = hg.get_bulk_scattering_properties_tro_gridded(
      AtmPoint{}, f_grid, std::make_shared<ZenithAngleGrid>(IrregularZenithAngleGrid(angles)));
  const auto& pm = *bulk.phase_matrix;
  for (Index iv = 0; iv < 3; iv++) {
    const Numeric scattering = ext(f_grid[iv]) * ssa(f_grid[iv]);
    if (not close(bulk.extinction_matrix[0, iv, 0], ext(f_grid[iv]), 1e-14)) return false;
    if (not close(bulk.absorption_vector[0, iv, 0], ext(f_grid[iv]) - scattering, 1e-14)) return false;

    Numeric integral = 0.0;
    for (Index i = 0; i < n; i++) {
      const Numeric p = scattering * henyey_greenstein(g, x[n - 1 - i]);
      for (Index is : {0, 2, 3, 5}) {
        if (not close(pm[0, iv, i, is], p, 1e-12)) {
          std::cout << "f = " << f_grid[iv] << " Hz, F" << is << ": " << pm[0, iv, i, is] << " != " << p << '\n';
          return false;
        }
      }
      if (pm[0, iv, i, 1] != 0.0 or pm[0, iv, i, 4] != 0.0) return false;
      integral += 2.0 * Constant::pi * w[n - 1 - i] * pm[0, iv, i, 0];
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
        const auto z               = [&](Index row, Index col) { return pm[0, 0, i, k, j, 4 * row + col]; };
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
  const auto hg_ext = ScatteringSpeciesProperty{"hg", ParticulateProperty::Extinction};
  const auto hg_ssa = ScatteringSpeciesProperty{"hg", ParticulateProperty::SingleScatteringAlbedo};
  AtmPoint   point;
  point[hg_ext] = 2e-3;
  point[hg_ssa] = 0.7;

  // Strongly forward peaked, where tabulating Theta is least accurate
  const Numeric                   g = 0.9;
  const HenyeyGreensteinScatterer hg{hg_ext, hg_ssa, g};
  const auto f = [g](Numeric c) { return std::array{henyey_greenstein(g, c), 0.0, henyey_greenstein(g, c)}; };

  const auto bulk = hg.get_bulk_scattering_properties_aro_gridded(point, Vector{1e9}, za_inc, delta_aa, za_scat_grid());
  if (not check_lab_frame(bulk, 2e-3 * 0.7, f, "HG")) return false;
  for (Size i = 0; i < za_inc.size(); i++) {
    if (not close(bulk.extinction_matrix[0, 0, i, 0], 2e-3, 1e-14)) return false;
    if (not close(bulk.absorption_vector[0, 0, i, 0], 2e-3 * 0.3, 1e-14)) return false;
  }

  const auto d_ext = hg.get_bulk_scattering_properties_aro_gridded_derivative(
      point, Vector{1e9}, za_inc, delta_aa, za_scat_grid(), AtmKeyVal{hg_ext});
  if (not check_lab_frame(d_ext, 0.7, f, "dHG/dextinction")) return false;

  const auto d_ssa = hg.get_bulk_scattering_properties_aro_gridded_derivative(
      point, Vector{1e9}, za_inc, delta_aa, za_scat_grid(), AtmKeyVal{hg_ssa});
  return check_lab_frame(d_ssa, 2e-3, f, "dHG/dssa");
}

bool test_rayleigh_lab_frame() {
  const AtmPoint     point{101325.0, 288.15};
  const GasScatterer gas{ConstantGasScattering{4.65e-31}, RayleighGasScattering{0.0}};

  const auto bulk =
      gas.get_bulk_scattering_properties_aro_gridded(point, Vector{6e14}, za_inc, delta_aa, za_scat_grid());
  const Numeric k = bulk.extinction_matrix[0, 0, 0, 0];
  if (not check_lab_frame(bulk, k / (4.0 * Constant::pi), rayleigh, "Rayleigh")) return false;

  const auto d_p = gas.get_bulk_scattering_properties_aro_gridded_derivative(
      point, Vector{6e14}, za_inc, delta_aa, za_scat_grid(), AtmKeyVal{AtmKey::p});
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
  if (not test_rayleigh_lab_frame()) {
    std::cerr << "Rayleigh lab-frame data failed\n";
    return 1;
  }
  return 0;
}
