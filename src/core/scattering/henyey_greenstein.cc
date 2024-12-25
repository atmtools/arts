#include "henyey_greenstein.h"

#include <cmath>

#include "math_funcs.h"

namespace scattering {
HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(
    ExtSSACallback ext_ssa_callback_, const Numeric& g_)
    : ext_ssa_callback(std::move(ext_ssa_callback_)), g(g_) {
  if (std::abs(g) > 1)
    throw std::runtime_error(
        "The Henyey-Greenstein asymmetry parameter g must be in the range [-1, 1].");
}

HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(
    ScatteringSpeciesProperty extinction_field,
    ScatteringSpeciesProperty ssa_field,
    const Numeric& g_)
    : ext_ssa_callback(ExtinctionSSALookup(std::move(extinction_field),
                                           std::move(ssa_field))),
      g(g_) {
  if (std::abs(g) > 1)
    throw std::runtime_error(
        "The Henyey-Greenstein asymmetry parameter g must be in the range [-1, 1].");
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const {
  auto t_grid     = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);

  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> pm{
      t_grid, f_grid_ptr, zenith_angle_grid};
  ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> emd{
      t_grid, f_grid_ptr};

  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> av{
      t_grid, f_grid_ptr};

  auto zenith_angles = grid_vector(*zenith_angle_grid);
  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    float extinction, ssa;
    std::tie(extinction, ssa) = ext_ssa_callback(f_grid[f_ind], atm_point);
    float scattering_xsec     = extinction * ssa;

    emd[0, f_ind, 0] = extinction;
    av[0, f_ind, 0]  = extinction - scattering_xsec;
    Numeric g2       = g * g;
    for (Size ind = 0; ind < zenith_angles.size(); ++ind) {
      pm[0, f_ind, ind, 0]  = (1.0 - g2);
      pm[0, f_ind, ind, 0] /= std::pow(
          1.0 + g2 - 2.0 * g * cos(Conversion::deg2rad(zenith_angles[ind])),
          3.0 / 2.0);
    }
    pm *= scattering_xsec / (4.0 * Constant::pi);
  }

  return BulkScatteringProperties<Format::TRO, Representation::Gridded>{
      pm, emd, av};
}

ScatteringTroSpectralVector
HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index l) const {
  auto t_grid     = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);

  SpecmatMatrix pm(f_grid.size(), l + 1);
  PropmatVector emd(f_grid.size());
  StokvecVector av(f_grid.size());

  constexpr Numeric inv_sphere = 0.5 * Constant::inv_sqrt_pi;  // sqrt(1/4pi)
  Vector f(l + 1, inv_sphere);
  for (Index ind = 1; ind <= l; ind++) {
    f[ind] *= std::sqrt(2 * ind + 1) * std::pow(g, ind);
  }

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    const auto [extinction, ssa] = ext_ssa_callback(f_grid[f_ind], atm_point);
    const auto scattering_xsec   = extinction * ssa;

    for (Index ind = 0; ind <= l; ++ind) {
      for (Index i = 0; i < 4; i++)
        pm[f_ind, ind][i, i] = f[ind] * scattering_xsec;
    }

    emd[f_ind].A() = extinction;
    av[f_ind].I()  = extinction - scattering_xsec;
  }

  return {.phase_matrix      = std::move(pm),
          .extinction_matrix = std::move(emd),
          .absorption_vector = std::move(av)};
}

BulkScatteringProperties<scattering::Format::ARO,
                         scattering::Representation::Gridded>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& delta_aa_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  auto bsp_tro = get_bulk_scattering_properties_tro_gridded(
      atm_point, f_grid, za_scat_grid);
  return bsp_tro.to_lab_frame(std::make_shared<Vector>(za_inc_grid),
                              std::make_shared<Vector>(delta_aa_grid),
                              za_scat_grid);
}

BulkScatteringProperties<scattering::Format::ARO,
                         scattering::Representation::Spectral>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(
    const AtmPoint& atm_point,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    Index degree,
    Index order) const {
  auto sht_ptr          = sht::provider.get_instance(degree, order);
  auto aa_scat_grid_ptr = sht_ptr->get_aa_grid_ptr();
  auto za_scat_grid_ptr =
      std::make_shared<ZenithAngleGrid>(sht_ptr->get_zenith_angle_grid());
  auto bsp_tro = get_bulk_scattering_properties_tro_gridded(
      atm_point, f_grid, za_scat_grid_ptr);
  auto bsp_aro = bsp_tro.to_lab_frame(std::make_shared<Vector>(za_inc_grid),
                                      aa_scat_grid_ptr,
                                      za_scat_grid_ptr);
  return bsp_aro.to_spectral(degree, order);
}

std::ostream& operator<<(std::ostream& os,
                         const HenyeyGreensteinScatterer& scatterer) {
  return os << "HenyeyGreensteinScatterer(g = " << scatterer.g << ")";
}
}  // namespace scattering
