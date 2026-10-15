#include <arts_conversions.h>
#include <physics_funcs.h>

#include <cmath>
#include <iostream>

#include "gas_scattering.h"

namespace {

bool close(Numeric a, Numeric b, Numeric relative = 1e-10) {
  return std::abs(a - b) <= relative * std::max({Numeric{1e-300}, std::abs(a), std::abs(b)});
}

bool test_constant_isotropic() {
  const AtmPoint point{101325.0, 288.15};
  const Numeric  cross_section = 4.65e-31;
  const Numeric  expected      = cross_section * number_density(point.pressure, point.temperature);
  const Vector   frequencies{6e14};
  auto           angles =
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{0.0, 90.0, 180.0}));

  const scattering::GasScatterer gas{scattering::ConstantGasScattering{cross_section},
                                     scattering::IsotropicGasScattering{}};
  auto bulk = gas.get_bulk_scattering_properties_tro_gridded(point, frequencies, std::move(angles));

  if (not close(bulk.extinction_matrix[0, 0, 0], expected)) return false;
  if (bulk.absorption_vector[0, 0, 0] != 0.0) return false;
  for (Index ia = 0; ia < 3; ++ia) {
    for (Index is = 0; is < 6; ++is) {
      const Numeric normalized = (is == 0 or is == 2 or is == 3 or is == 5) ? 1.0 : 0.0;
      if (not close((*bulk.phase_matrix)[0, 0, ia, is], expected * normalized / (4.0 * Constant::pi))) { return false; }
    }
  }
  return true;
}

bool test_air_simple_rayleigh() {
  const AtmPoint point{101325.0, 288.15};
  const Vector   frequencies{6e14};
  auto           angles =
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{0.0, 90.0, 180.0}));

  const scattering::GasScatterer gas{scattering::AirSimpleGasScattering{}, scattering::RayleighGasScattering{0.0}};
  auto bulk = gas.get_bulk_scattering_properties_tro_gridded(point, frequencies, std::move(angles));

  const Numeric wavelength_um     = Conversion::freq2wavelen(frequencies[0]) * 1e6;
  const Numeric inv_wavelength_sq = 1.0 / Math::pow2(wavelength_um);
  const Numeric polynomial =
      3.9729066 +
      inv_wavelength_sq * (4.6547659e-2 + inv_wavelength_sq * (4.5055995e-4 + inv_wavelength_sq * 2.3229848e-5));
  const Numeric expected =
      1e-32 * polynomial / Math::pow4(wavelength_um) * number_density(point.pressure, point.temperature);
  if (not close(bulk.extinction_matrix[0, 0, 0], expected)) return false;

  const Numeric scale = expected / (4.0 * Constant::pi);
  const Vector  at_90 = calc_rayleighPhaMat(Constant::pi / 2.0);
  for (Index is = 0; is < 6; ++is) {
    if (not close((*bulk.phase_matrix)[0, 0, 1, is], scale * at_90[is], 1e-9)) return false;
  }

  // The species must also provide the spectral TRO representation consumed by
  // the scattering-species agenda and DISORT.
  const auto spectral = gas.get_bulk_scattering_properties_tro_spectral(point, frequencies, 2);
  if (not close(spectral.extinction_matrix[0].A(), expected)) return false;
  if (spectral.absorption_vector[0].I() != 0.0) return false;
  if (not spectral.phase_matrix or spectral.phase_matrix->ncols() != 3) return false;
  const Numeric amp = 0.5 * Constant::inv_sqrt_pi * expected;
  if (not close((*spectral.phase_matrix)[0, 0][0, 0].real(), amp)) return false;
  if (not close((*spectral.phase_matrix)[0, 1][2, 2].real(), 0.5 * std::sqrt(3.0) * amp)) return false;
  if (not close((*spectral.phase_matrix)[0, 2][0, 0].real(), 0.5 * amp / std::sqrt(5.0))) return false;

  return true;
}

/* The SHT of Rayleigh scattering (aro_spectral) and the Fourier modes taken from it, against the exact modes
   (aro_fourier).  The laboratory-frame Z11 is a spherical-harmonic series of degree 2 over the scattering
   directions, so an SHT of degree and order 8 holds it exactly, and its modes at any scattering zenith angle and
   its phase integral are exact.  The polarized elements are not compared: referenced to the meridian, which
   turns with the azimuth at the poles, they are not finite series in scalar spherical harmonics. */
#ifndef ARTS_NO_SHTNS
bool test_rayleigh_sht_to_fourier() {
  const AtmPoint point{8e4, 260.0};
  const Vector   frequencies{89e9};
  const Vector   za_inc{0.0, 35.0, 120.0};
  const Vector   za_scat{0.0, 20.0, 35.0, 90.0, 157.0, 180.0};
  const Index    M = 3;

  const scattering::GasScatterer gas{scattering::ConstantGasScattering{1e-30}, scattering::RayleighGasScattering{0.0}};
  const auto sht   = gas.get_bulk_scattering_properties_aro_spectral(point, frequencies, za_inc, 8, 8);
  const auto exact = gas.get_bulk_scattering_properties_aro_fourier(point, frequencies, za_inc, za_scat, M);
  const auto modes = sht.phase_matrix->to_fourier(
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(za_scat)), M);

  Numeric scale = 0.0;
  for (Numeric x : *exact.phase_matrix | by_elem) scale = std::max(scale, std::abs(x));
  if (not(scale > 0.0) or modes[0, 0, 1, 2, 2, 0, 0] == 0.0) return false;  // a non-trivial comparison, m = 2 included
  for (Size i = 0; i < za_inc.size(); ++i) {
    if (not close(modes.get_phase_integral()[0, 0, i], exact.phase_matrix->get_phase_integral()[0, 0, i], 1e-12))
      return false;
    for (Size j = 0; j < za_scat.size(); ++j)
      for (Index m = 0; m <= M; ++m)
        for (Index cs = 0; cs < 2; ++cs)
          if (std::abs(modes[0, 0, i, j, m, cs, 0] - (*exact.phase_matrix)[0, 0, i, j, m, cs, 0]) > 1e-12 * scale)
            return false;
  }
  return true;
}
#endif

}  // namespace

int main() {
#ifndef ARTS_NO_SHTNS
  if (not test_rayleigh_sht_to_fourier()) {
    std::cerr << "Rayleigh SHT to Fourier modes failed\n";
    return 1;
  }
#endif
  if (not test_constant_isotropic()) {
    std::cerr << "Constant/isotropic gas scattering failed\n";
    return 1;
  }
  if (not test_air_simple_rayleigh()) {
    std::cerr << "AirSimple/Rayleigh gas scattering failed\n";
    return 1;
  }
  return 0;
}
