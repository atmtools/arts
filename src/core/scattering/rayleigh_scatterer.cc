#include "rayleigh_scatterer.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <physics_funcs.h>
#include <species.h>

#include <array>
#include <cmath>

#include "mie.h"

namespace scattering {

RayleighScatterer::RayleighScatterer(RayleighType t, Numeric d)
    : type(t), diameter(d) {
  if (type == RayleighType::WaterDrop) {
    ARTS_USER_ERROR_IF(
        diameter <= 0, "WaterDrop requires diameter > 0, got {}", diameter);
  }
}

// -----------------------------------------------------------------------
//  AirSimple cross-section model
// -----------------------------------------------------------------------
namespace {
/// Empirical air Rayleigh coefficients (same as spectral_propmat_scatAirSimple).

std::pair<Numeric, Numeric> earth_air_simple(Numeric freq_hz,
                                             const AtmPoint& atm_point) {
  constexpr std::array air_coefficients{
      3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};

  const Numeric nd = number_density(atm_point.pressure, atm_point.temperature);
  const Numeric wavelen = Conversion::freq2wavelen(freq_hz) * 1e6;
  Numeric sum           = 0;
  Numeric pows          = 1;
  for (auto& c : air_coefficients) {
    sum  += c * pows;
    pows /= Math::pow2(wavelen);
  }
  return {1e-32 * nd * sum / Math::pow4(wavelen), 0.0};
}

// -----------------------------------------------------------------------
//  WaterDrop cross-section model (Rayleigh limit, Liebe 1993 dielectric)
// -----------------------------------------------------------------------
std::pair<Numeric, Numeric> water_drop(Numeric freq_hz,
                                       const AtmPoint& atm_point,
                                       Numeric diameter) {
  constexpr Numeric rho_water = 1.0e3;  // [kg/m³]
  const Numeric lwc           = atm_point["liquidcloud"_spec];
  if (lwc <= 0) return {0.0, 0.0};

  const Numeric vol = Constant::pi / 6.0 * diameter * diameter * diameter;
  const Numeric nd  = lwc / (rho_water * vol);

  const Numeric T = atm_point.temperature;
  RayleighSphere<Numeric> sphere(
      freq_hz, T, diameter, refr_index_water_liebe93(freq_hz, T));

  return {nd * sphere.get_scattering_coeff(),
          nd * sphere.get_absorption_coeff()};
}
}  // namespace

std::pair<Numeric, Numeric> RayleighScatterer::cross_sections(
    Numeric freq_hz, const AtmPoint& atm_point) const {
  switch (type) {
    case RayleighType::EarthAir: return earth_air_simple(freq_hz, atm_point);
    case RayleighType::WaterDrop:
      return water_drop(freq_hz, atm_point, diameter);
  }
  std::unreachable();
}

// -----------------------------------------------------------------------
//  Interface methods — all dispatch through cross_sections()
// -----------------------------------------------------------------------

ScatteringTroSpectralVector
RayleighScatterer::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index l) const {
  const Index nf = f_grid.size();
  SpecmatMatrix pm(nf, l + 1);
  PropmatVector emd(nf);
  StokvecVector av(nf);

  // Rayleigh Legendre coefficients: β_0 = 1, β_2 = 1/10
  // ARTS SHT normalization: f_l = 1/(2√π) × √(2l+1) × β_l
  constexpr Numeric inv_sphere = 0.5 * Constant::inv_sqrt_pi;
  Vector coeffs(l + 1, 0.0);
  coeffs[0] = inv_sphere;
  if (l >= 2) coeffs[2] = inv_sphere * std::sqrt(Numeric{5}) * 0.1;

  for (Index iv = 0; iv < nf; iv++) {
    const auto [sca, abs] = cross_sections(f_grid[iv], atm_point);

    for (Index ind = 0; ind <= l; ind++) {
      for (Index i = 0; i < 4; i++) {
        pm[iv, ind][i, i] = coeffs[ind] * sca;
      }
    }

    emd[iv].A() = sca + abs;
    av[iv].I()  = abs;
  }

  return {.phase_matrix      = std::move(pm),
          .extinction_matrix = std::move(emd),
          .absorption_vector = std::move(av)};
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
RayleighScatterer::get_bulk_scattering_properties_tro_gridded(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  auto t_grid     = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);

  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> pm{
      t_grid, f_grid_ptr, za_scat_grid};
  ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> emd{
      t_grid, f_grid_ptr};
  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> av{
      t_grid, f_grid_ptr};

  const auto za  = grid_vector(*za_scat_grid);
  const auto nza = za.size();

  // Precompute angle-dependent Rayleigh phase function elements.
  // P_11 = (3/4)(1 + cos²θ),  P_12 = -(3/4)(1 - cos²θ),
  // P_33 = P_44 = (3/2)cosθ,  P_34 = 0.
  // To obtain phase matrix per unit solid angle: divide by 4π.
  constexpr Numeric inv_4pi = 1.0 / (4.0 * Constant::pi);
  std::vector<std::array<Numeric, 6>> pf(nza);
  for (Size ia = 0; ia < nza; ia++) {
    const Numeric ct  = std::cos(Conversion::deg2rad(za[ia]));
    const Numeric ct2 = ct * ct;
    pf[ia]            = {0.75 * (1.0 + ct2) * inv_4pi,   // Z11
                         -0.75 * (1.0 - ct2) * inv_4pi,  // Z12
                         0.75 * (1.0 + ct2) * inv_4pi,   // Z22
                         1.5 * ct * inv_4pi,             // Z33 = Z44
                         0.0,                            // Z34
                         1.5 * ct * inv_4pi};            // Z44
  }

  for (Size iv = 0; iv < f_grid.size(); iv++) {
    const auto [sca, abs] = cross_sections(f_grid[iv], atm_point);
    emd[0, iv, 0]         = sca + abs;
    av[0, iv, 0]          = abs;
    for (Size ia = 0; ia < nza; ia++) {
      for (Index k = 0; k < 6; k++) {
        pm[0, iv, ia, k] = pf[ia][k] * sca;
      }
    }
  }

  return {
      .phase_matrix = pm, .extinction_matrix = emd, .absorption_vector = av};
}

void RayleighScatterer::get_radar_single_scat(MuelmatVector& Z_back,
                                              PropmatVector& K_ext,
                                              const AtmPoint& atm_point,
                                              const Vector& f_grid) const {
  const Index nf = f_grid.size();
  Z_back         = MuelmatVector(nf, rtepack::muelmat{0.0});
  K_ext          = PropmatVector(nf, rtepack::propmat{});

  for (Index iv = 0; iv < nf; iv++) {
    const auto [sca, abs] = cross_sections(f_grid[iv], atm_point);
    if (sca <= 0 && abs <= 0) continue;

    // Rayleigh backscatter phase matrix: P(180°) = 3/2 for (1,1) and (2,2);
    // z11 = sca * P(180°) / (4π) = sca * 3 / (8π)
    const Numeric z11 = sca * 3.0 / (8.0 * Constant::pi);
    Z_back[iv]        = rtepack::muelmat{
        z11, 0, 0, 0, 0, z11, 0, 0, 0, 0, -z11, 0, 0, 0, 0, -z11};
    K_ext[iv].A() = sca + abs;
  }
}

std::ostream& operator<<(std::ostream& os, const RayleighScatterer& s) {
  return os << "RayleighScatterer(type=" << s.type
            << ", diameter=" << s.diameter << ")";
}
}  // namespace scattering

void xml_io_stream<scattering::RayleighScatterer>::write(
    std::ostream& os,
    const scattering::RayleighScatterer& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.diameter, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<scattering::RayleighScatterer>::read(
    std::istream& is, scattering::RayleighScatterer& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.diameter, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
