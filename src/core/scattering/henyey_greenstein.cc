#include "henyey_greenstein.h"

#include <cmath>

#include "math_funcs.h"

namespace scattering {
ExtinctionSSALookup::ExtinctionSSALookup(ScatteringSpeciesProperty extinction_field_,
                                         ScatteringSpeciesProperty ssa_field_)
    : extinction_field(std::move(extinction_field_)), ssa_field(std::move(ssa_field_)) {}

std::pair<Numeric, Numeric> ExtinctionSSALookup::operator()(Numeric, const AtmPoint& atm_point) {
  return {atm_point[extinction_field], atm_point[ssa_field]};
}

HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(ExtSSACallback ext_ssa_callback_, const Numeric& g_)
    : ext_ssa_callback(std::move(ext_ssa_callback_)), g(g_) {
  if (std::abs(g) > 1)
    throw std::runtime_error("The Henyey-Greenstein asymmetry parameter g must be in the range [-1, 1].");
}

HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(ScatteringSpeciesProperty extinction_field,
                                                     ScatteringSpeciesProperty ssa_field,
                                                     const Numeric&            g_)
    : ext_ssa_callback(ExtinctionSSALookup(std::move(extinction_field), std::move(ssa_field))), g(g_) {
  if (std::abs(g) > 1)
    throw std::runtime_error("The Henyey-Greenstein asymmetry parameter g must be in the range [-1, 1].");
}

namespace {
/** The Henyey-Greenstein phase function at scattering angle theta [rad], normalised to 1 over the sphere */
Numeric phase_function(Numeric g, Numeric theta) {
  const Numeric g2 = g * g;
  return (1.0 - g2) / (4.0 * Constant::pi * std::pow(1.0 + g2 - 2.0 * g * std::cos(theta), 1.5));
}

/** Extinction, absorption and scattering coefficient per frequency, or their derivatives */
struct Coefficients {
  std::shared_ptr<const Vector>                                       t_grid;
  std::shared_ptr<const Vector>                                       f_grid;
  ExtinctionMatrixData<Numeric, Format::TRO, Representation::Gridded> extinction;
  AbsorptionVectorData<Numeric, Format::TRO, Representation::Gridded> absorption;
  Vector                                                              scattering;

  explicit Coefficients(const Vector& f)
      : t_grid(std::make_shared<const Vector>(Vector{0.0})),
        f_grid(std::make_shared<const Vector>(f)),
        extinction(t_grid, f_grid),
        absorption(t_grid, f_grid),
        scattering(f.size()) {}
};

Coefficients coefficients(const ExtSSACallback& ext_ssa_callback, const AtmPoint& atm_point, const Vector& f_grid) {
  Coefficients out(f_grid);
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    const auto [extinction, ssa] = ext_ssa_callback(f_grid[iv], atm_point);
    out.extinction[0, iv, 0]     = extinction;
    out.scattering[iv]           = extinction * ssa;
    out.absorption[0, iv, 0]     = extinction - out.scattering[iv];
  }
  return out;
}

Coefficients coefficient_derivatives(const ExtSSACallback& ext_ssa_callback,
                                     const AtmPoint&       atm_point,
                                     const Vector&         f_grid,
                                     const AtmKeyVal&      target) {
  const auto* lookup = ext_ssa_callback.f.target<ExtinctionSSALookup>();
  ARTS_USER_ERROR_IF(not lookup,
                     "Analytical derivatives of a HenyeyGreensteinScatterer require its "
                     "extinction/SSA atmospheric-field constructor, not an arbitrary callback")

  const bool   d_ext = target == AtmKeyVal{lookup->extinction_field};
  const bool   d_ssa = target == AtmKeyVal{lookup->ssa_field};
  Coefficients out(f_grid);
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    const auto [extinction, ssa] = ext_ssa_callback(f_grid[iv], atm_point);
    out.extinction[0, iv, 0]     = d_ext ? 1.0 : 0.0;
    out.absorption[0, iv, 0]     = d_ext ? 1.0 - ssa : (d_ssa ? -extinction : 0.0);
    out.scattering[iv]           = d_ext ? ssa : (d_ssa ? extinction : 0.0);
  }
  return out;
}

BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(Coefficients                     c,
                                                                           Numeric                          g,
                                                                           std::shared_ptr<ZenithAngleGrid> za_grid) {
  PhaseMatrixData<Numeric, Format::TRO, Representation::Gridded> phase{c.t_grid, c.f_grid, za_grid};
  const auto                                                     angles = grid_vector(*za_grid);
  for (Size iv = 0; iv < c.f_grid->size(); ++iv) {
    for (Size ia = 0; ia < angles.size(); ++ia) {
      const Numeric p     = c.scattering[iv] * phase_function(g, Conversion::deg2rad(angles[ia]));
      phase[0, iv, ia, 0] = p;
      phase[0, iv, ia, 2] = p;
      phase[0, iv, ia, 3] = p;
      phase[0, iv, ia, 5] = p;
    }
  }
  return {std::move(phase), std::move(c.extinction), std::move(c.absorption)};
}

/** The lab-frame phase matrix at the exact scattering angle of every direction pair */
BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
    const Coefficients&              c,
    Numeric                          g,
    const Vector&                    za_inc_grid,
    const Vector&                    delta_aa_grid,
    std::shared_ptr<ZenithAngleGrid> za_scat_grid) {
  auto       za_inc            = std::make_shared<const Vector>(za_inc_grid);
  const auto henyey_greenstein = [&](Numeric theta, matpack::data_t<Numeric, 3>& scattering_matrix) {
    const Numeric p = phase_function(g, theta);
    for (Size iv = 0; iv < c.f_grid->size(); ++iv) {
      const Numeric z             = c.scattering[iv] * p;
      scattering_matrix[0, iv, 0] = z;
      scattering_matrix[0, iv, 1] = 0.0;
      scattering_matrix[0, iv, 2] = z;
      scattering_matrix[0, iv, 3] = z;
      scattering_matrix[0, iv, 4] = 0.0;
      scattering_matrix[0, iv, 5] = z;
    }
  };
  auto phase = tro_lab_frame<Numeric>(c.t_grid,
                                      c.f_grid,
                                      za_inc,
                                      std::make_shared<const Vector>(delta_aa_grid),
                                      std::move(za_scat_grid),
                                      henyey_greenstein);
  return {std::move(phase), c.extinction.to_lab_frame(za_inc), c.absorption.to_lab_frame(za_inc)};
}
}  // namespace

BulkScatteringProperties<Format::TRO, Representation::Gridded>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded(
    const AtmPoint& atm_point, const Vector& f_grid, std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const {
  return tro_gridded(coefficients(ext_ssa_callback, atm_point, f_grid), g, std::move(zenith_angle_grid));
}

BulkScatteringProperties<Format::TRO, Representation::Gridded>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_gridded_derivative(
    const AtmPoint&                  atm_point,
    const Vector&                    f_grid,
    std::shared_ptr<ZenithAngleGrid> zenith_angle_grid,
    const AtmKeyVal&                 target) const {
  return tro_gridded(
      coefficient_derivatives(ext_ssa_callback, atm_point, f_grid, target), g, std::move(zenith_angle_grid));
}

ScatteringTroSpectralVector HenyeyGreensteinScatterer::get_bulk_scattering_properties_tro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, Index l) const {
  auto t_grid     = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);

  SpecmatMatrix pm(f_grid.size(), l + 1, Specmat(0.0));
  PropmatVector emd(f_grid.size());
  StokvecVector av(f_grid.size());

  constexpr Numeric inv_sphere = 0.5 * Constant::inv_sqrt_pi;  // sqrt(1/4pi)
  Vector            f(l + 1, inv_sphere);
  for (Index ind = 1; ind <= l; ind++) { f[ind] *= std::sqrt(2 * ind + 1) * std::pow(g, ind); }

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    const auto [extinction, ssa] = ext_ssa_callback(f_grid[f_ind], atm_point);
    const auto scattering_xsec   = extinction * ssa;

    for (Index ind = 0; ind <= l; ++ind) {
      for (Index i = 0; i < 4; i++) pm[f_ind, ind][i, i] = f[ind] * scattering_xsec;
    }

    emd[f_ind].A() = extinction;
    av[f_ind].I()  = extinction - scattering_xsec;
  }

  return {.phase_matrix = std::move(pm), .extinction_matrix = std::move(emd), .absorption_vector = std::move(av)};
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded(
    const AtmPoint&                              atm_point,
    const Vector&                                f_grid,
    const Vector&                                za_inc_grid,
    const Vector&                                delta_aa_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid) const {
  return aro_gridded(
      coefficients(ext_ssa_callback, atm_point, f_grid), g, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_gridded_derivative(
    const AtmPoint&                              atm_point,
    const Vector&                                f_grid,
    const Vector&                                za_inc_grid,
    const Vector&                                delta_aa_grid,
    std::shared_ptr<scattering::ZenithAngleGrid> za_scat_grid,
    const AtmKeyVal&                             target) const {
  return aro_gridded(coefficient_derivatives(ext_ssa_callback, atm_point, f_grid, target),
                     g,
                     za_inc_grid,
                     delta_aa_grid,
                     std::move(za_scat_grid));
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_spectral(
    const AtmPoint& atm_point, const Vector& f_grid, const Vector& za_inc_grid, Index degree, Index order) const {
  auto sht_ptr          = sht::provider.get_instance(degree, order);
  auto aa_scat_grid_ptr = sht_ptr->get_aa_grid_ptr();
  auto za_scat_grid_ptr = std::make_shared<ZenithAngleGrid>(sht_ptr->get_zenith_angle_grid());
  auto bsp_tro          = get_bulk_scattering_properties_tro_gridded(atm_point, f_grid, za_scat_grid_ptr);
  auto bsp_aro = bsp_tro.to_lab_frame(std::make_shared<Vector>(za_inc_grid), aa_scat_grid_ptr, za_scat_grid_ptr);
  return bsp_aro.to_spectral(degree, order);
}

std::ostream& operator<<(std::ostream& os, const HenyeyGreensteinScatterer& scatterer) {
  return os << "HenyeyGreensteinScatterer(g = " << scatterer.g << ")";
}
}  // namespace scattering

void xml_io_stream<scattering::HenyeyGreensteinScatterer>::write(std::ostream&                                os,
                                                                 const scattering::HenyeyGreensteinScatterer& x,
                                                                 bofstream*                                   pbofs,
                                                                 std::string_view                             name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.ext_ssa_callback, pbofs);
  xml_write_to_stream(os, x.g, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<scattering::HenyeyGreensteinScatterer>::read(std::istream&                          is,
                                                                scattering::HenyeyGreensteinScatterer& x,
                                                                bifstream*                             pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.ext_ssa_callback, pbifs);
  xml_read_from_stream(is, x.g, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
