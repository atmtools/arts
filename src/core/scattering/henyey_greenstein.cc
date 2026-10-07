#include "henyey_greenstein.h"

#include <cmath>

#include "math_funcs.h"
#include "sht.h"

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

/** F33 / F11 = F44 / F11, odd in cos(Theta), 1 forward and -1 backward
 *
 * 1 + f = (1 + c)^2 (2 - c) / 2 and 1 - f = (1 - c)^2 (2 + c) / 2 have the
 * double zeros that make the matrix regular at forward and backward scattering.
 */
Numeric polarization_ratio(Numeric cos_theta) { return 0.5 * cos_theta * (3.0 - cos_theta * cos_theta); }

/** The coefficients on Y_k0 of cos(Theta) times the function with coefficients a on Y_k0, one degree fewer
 *
 * cos(Theta) Y_k0 = A(k) Y_{k+1,0} + A(k - 1) Y_{k-1,0} with A(k) = (k + 1) / sqrt((2 k + 1) (2 k + 3)).
 */
Vector times_cos(const Vector& a) {
  const auto A = [](Index k) {
    const auto x = static_cast<Numeric>(k);
    return (x + 1.0) / std::sqrt((2.0 * x + 1.0) * (2.0 * x + 3.0));
  };
  const Index n = static_cast<Index>(a.size()) - 1;
  Vector      b(n);
  for (Index k = 0; k < n; k++) b[k] = (k == 0 ? 0.0 : a[k - 1] * A(k - 1)) + a[k + 1] * A(k);
  return b;
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
      const Numeric theta = Conversion::deg2rad(angles[ia]);
      const Numeric p     = c.scattering[iv] * phase_function(g, theta);
      phase[0, iv, ia, 0] = p;
      phase[0, iv, ia, 2] = p;
      phase[0, iv, ia, 3] = p * polarization_ratio(std::cos(theta));
      phase[0, iv, ia, 5] = p * polarization_ratio(std::cos(theta));
    }
  }
  return {std::move(phase), std::move(c.extinction), std::move(c.absorption)};
}

/** The scattering matrix at scattering angle theta [rad] as tro_lab_frame takes it, from the closed form */
auto scattering_matrix(const Coefficients& c, Numeric g) {
  return [&c, g](Numeric theta, matpack::data_t<Numeric, 3>& scattering_matrix) {
    const Numeric p = phase_function(g, theta);
    const Numeric f = polarization_ratio(std::cos(theta));
    for (Size iv = 0; iv < c.f_grid->size(); ++iv) {
      const Numeric z             = c.scattering[iv] * p;
      scattering_matrix[0, iv, 0] = z;
      scattering_matrix[0, iv, 1] = 0.0;
      scattering_matrix[0, iv, 2] = z;
      scattering_matrix[0, iv, 3] = z * f;
      scattering_matrix[0, iv, 4] = 0.0;
      scattering_matrix[0, iv, 5] = z * f;
    }
  };
}

/** The lab-frame phase matrix at the exact scattering angle of every direction pair */
BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
    const Coefficients&              c,
    Numeric                          g,
    const Vector&                    za_inc_grid,
    const Vector&                    delta_aa_grid,
    std::shared_ptr<ZenithAngleGrid> za_scat_grid) {
  auto za_inc = std::make_shared<const Vector>(za_inc_grid);
  auto phase  = tro_lab_frame<Numeric>(c.t_grid,
                                       c.f_grid,
                                       za_inc,
                                       std::make_shared<const Vector>(delta_aa_grid),
                                       std::move(za_scat_grid),
                                       scattering_matrix(c, g));
  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = c.extinction.to_lab_frame(za_inc),
          .absorption_vector = c.absorption.to_lab_frame(za_inc)};
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

  // The coefficients of p and of p f on Y_k0 = sqrt((2 k + 1) / 4 pi) P_k(cos(Theta)), from
  // p = sum_k sqrt((2 k + 1) / 4 pi) g^k Y_k0 to degree l + 3, which the cos^3 term of f needs
  constexpr Numeric inv_sphere = 0.5 * Constant::inv_sqrt_pi;  // sqrt(1/4pi)
  Vector            f(l + 4);
  for (Index k = 0; k <= l + 3; k++) f[k] = inv_sphere * std::sqrt(static_cast<Numeric>(2 * k + 1)) * std::pow(g, k);
  const Vector c1 = times_cos(f), c3 = times_cos(times_cos(c1));
  Vector       f_pol(l + 1);
  for (Index k = 0; k <= l; k++) f_pol[k] = 0.5 * (3.0 * c1[k] - c3[k]);

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    const auto [extinction, ssa] = ext_ssa_callback(f_grid[f_ind], atm_point);
    const auto scattering_xsec   = extinction * ssa;

    for (Index ind = 0; ind <= l; ++ind) {
      pm[f_ind, ind][0, 0] = f[ind] * scattering_xsec;
      pm[f_ind, ind][1, 1] = f[ind] * scattering_xsec;
      pm[f_ind, ind][2, 2] = f_pol[ind] * scattering_xsec;
      pm[f_ind, ind][3, 3] = f_pol[ind] * scattering_xsec;
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
  auto sht = sht::provider.get_instance_lm(degree, order);
  return get_bulk_scattering_properties_aro_gridded(atm_point,
                                                    f_grid,
                                                    za_inc_grid,
                                                    *sht->get_aa_grid_ptr(),
                                                    std::make_shared<ZenithAngleGrid>(sht->get_zenith_angle_grid()))
      .to_spectral(degree, order);
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>
HenyeyGreensteinScatterer::get_bulk_scattering_properties_aro_fourier(const AtmPoint& atm_point,
                                                                      const Vector&   f_grid,
                                                                      const Vector&   za_inc_grid,
                                                                      const Vector&   za_scat_grid,
                                                                      Index           max_mode) const {
  const auto c      = coefficients(ext_ssa_callback, atm_point, f_grid);
  auto       za_inc = std::make_shared<const Vector>(za_inc_grid);
  Matrix     integral(1, f_grid.size());  // p integrates to 1, so the phase integral is the scattering coefficient
  integral[0, joker] = c.scattering;
  auto phase         = tro_lab_frame_fourier_modes<Numeric>(
      c.t_grid,
      c.f_grid,
      za_inc,
      std::make_shared<const ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat_grid)),
      max_mode,
      integral,
      scattering_matrix(c, g));
  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = c.extinction.to_lab_frame(za_inc).to_fourier(),
          .absorption_vector = c.absorption.to_lab_frame(za_inc).to_fourier()};
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
