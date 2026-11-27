#include "henyey_greenstein.h"

#include <cmath>

#include "math_funcs.h"
#include "sht.h"

namespace scattering {
HenyeyGreensteinScatterer::HenyeyGreensteinScatterer(String extinction_name_, String ssa_name_, String g_name_)
    : extinction_name(std::move(extinction_name_)), ssa_name(std::move(ssa_name_)), g_name(std::move(g_name_)) {}

namespace {
/** The Henyey-Greenstein phase function at scattering angle theta [rad], normalised to 1 over the sphere */
Numeric phase_function(Numeric g, Numeric theta) {
  const Numeric g2 = g * g;
  return (1.0 - g2) / (4.0 * Constant::pi * std::pow(1.0 + g2 - 2.0 * g * std::cos(theta), 1.5));
}

/** The derivative of the Henyey-Greenstein phase function by g at scattering angle theta [rad]
 *
 * With D = 1 + g^2 - 2 g cos(Theta), p = (1 - g^2) / (4 pi D^(3/2)) and
 * dp/dg = (-2 g D - 3 (1 - g^2) (g - cos(Theta))) / (4 pi D^(5/2)).
 */
Numeric phase_function_derivative(Numeric g, Numeric theta) {
  const Numeric c = std::cos(theta), D = 1.0 + g * g - 2.0 * g * c;
  return (-2.0 * g * D - 3.0 * (1.0 - g * g) * (g - c)) / (4.0 * Constant::pi * std::pow(D, 2.5));
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
  std::shared_ptr<const Vector>                              t_grid;
  std::shared_ptr<const Vector>                              f_grid;
  ExtinctionMatrixData<Format::TRO, Representation::Gridded> extinction;
  AbsorptionVectorData<Format::TRO, Representation::Gridded> absorption;
  Vector                                                     scattering;

  explicit Coefficients(const Vector& f)
      : t_grid(std::make_shared<const Vector>(Vector{0.0})),
        f_grid(std::make_shared<const Vector>(f)),
        extinction(t_grid, f_grid),
        absorption(t_grid, f_grid),
        scattering(f.size()) {}
};

void check_sizes(const Vector& extinction, const Vector& ssa, const Vector& g, const Vector& f_grid) {
  if (extinction.size() != f_grid.size() or ssa.size() != f_grid.size() or g.size() != f_grid.size())
    arts::user_error(
        "There are {} extinctions, {} single scattering albedos and {} asymmetry parameters, but {} frequencies",
        extinction.size(),
        ssa.size(),
        g.size(),
        f_grid.size());
}

Coefficients coefficients(const Vector& extinction, const Vector& ssa, const Vector& f_grid) {
  Coefficients out(f_grid);
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    out.extinction[0, iv].A() = extinction[iv];
    out.scattering[iv]        = extinction[iv] * ssa[iv];
    out.absorption[0, iv].I() = extinction[iv] - out.scattering[iv];
  }
  return out;
}

//! The derivatives of the coefficients, from those of the extinction and the single scattering albedo
Coefficients coefficient_derivatives(const Vector& extinction,
                                     const Vector& ssa,
                                     const Vector& d_extinction,
                                     const Vector& d_ssa,
                                     const Vector& f_grid) {
  Coefficients out(f_grid);
  for (Size iv = 0; iv < f_grid.size(); ++iv) {
    out.extinction[0, iv].A() = d_extinction[iv];
    out.scattering[iv]        = d_extinction[iv] * ssa[iv] + extinction[iv] * d_ssa[iv];
    out.absorption[0, iv].I() = d_extinction[iv] - out.scattering[iv];
  }
  return out;
}

//! The scattering coefficient times the phase function at frequency iv and scattering angle theta [rad]
struct Scattered {
  const Vector& scattering;
  const Vector& g;

  Numeric operator()(Size iv, Numeric theta) const { return scattering[iv] * phase_function(g[iv], theta); }
};

//! The derivative of Scattered, from those of the scattering coefficient and of g
struct ScatteredDerivative {
  const Vector& scattering;
  const Vector& d_scattering;
  const Vector& g;
  const Vector& d_g;

  Numeric operator()(Size iv, Numeric theta) const {
    return d_scattering[iv] * phase_function(g[iv], theta) +
           scattering[iv] * phase_function_derivative(g[iv], theta) * d_g[iv];
  }
};

/** The scattering-plane properties: the extinction and absorption of c, and the phase matrix of scattered */
template <typename ScatteredFunction>
BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_properties(
    const Coefficients& c, const ScatteredFunction& scattered, std::shared_ptr<ZenithAngleGrid> za_grid) {
  PhaseMatrixData<Format::TRO, Representation::Gridded> phase{c.t_grid, c.f_grid, za_grid};
  const auto                                            angles = grid_vector(*za_grid);
  for (Size iv = 0; iv < c.f_grid->size(); ++iv) {
    for (Size ia = 0; ia < angles.size(); ++ia) {
      const Numeric theta = Conversion::deg2rad(angles[ia]);
      const Numeric p     = scattered(iv, theta);
      const Numeric pf    = p * polarization_ratio(std::cos(theta));
      phase[0, iv, ia]    = rtepack::compact_planar_muelmat{p, 0.0, p, pf, 0.0, pf};
    }
  }
  return {.phase_matrix = std::move(phase), .extinction_matrix = c.extinction, .absorption_vector = c.absorption};
}

/** The scattering matrix at scattering angle theta [rad] as tro_lab_frame takes it, from the closed form */
template <typename ScatteredFunction> auto scattering_matrix(Size nf, const ScatteredFunction& scattered) {
  return [nf, &scattered](Numeric theta, rtepack::compact_planar_muelmat_matrix& scattering_matrix) {
    const Numeric f = polarization_ratio(std::cos(theta));
    for (Size iv = 0; iv < nf; ++iv) {
      const Numeric z          = scattered(iv, theta);
      scattering_matrix[0, iv] = rtepack::compact_planar_muelmat{z, 0.0, z, z * f, 0.0, z * f};
    }
  };
}

/** The lab-frame phase matrix of scattered at the exact scattering angle of every direction pair */
template <typename ScatteredFunction>
BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_properties(
    const Coefficients&              c,
    const ScatteredFunction&         scattered,
    const Vector&                    za_inc_grid,
    const Vector&                    delta_aa_grid,
    std::shared_ptr<ZenithAngleGrid> za_scat_grid) {
  auto za_inc = std::make_shared<const Vector>(za_inc_grid);
  auto phase  = tro_lab_frame(c.t_grid,
                              c.f_grid,
                              za_inc,
                              std::make_shared<const Vector>(delta_aa_grid),
                              std::move(za_scat_grid),
                              scattering_matrix(c.f_grid->size(), scattered));
  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = c.extinction.to_lab_frame(za_inc),
          .absorption_vector = c.absorption.to_lab_frame(za_inc)};
}
}  // namespace

BulkScatteringProperties<Format::TRO, Representation::Gridded> HenyeyGreensteinScatterer::tro_gridded(
    const Vector&                    extinction,
    const Vector&                    ssa,
    const Vector&                    g,
    const Vector&                    f_grid,
    std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const {
  check_sizes(extinction, ssa, g, f_grid);
  const auto c = coefficients(extinction, ssa, f_grid);
  return tro_gridded_properties(c, Scattered{c.scattering, g}, std::move(zenith_angle_grid));
}

BulkScatteringProperties<Format::TRO, Representation::Gridded> HenyeyGreensteinScatterer::tro_gridded_derivative(
    const Vector&                    extinction,
    const Vector&                    ssa,
    const Vector&                    g,
    const Vector&                    d_extinction,
    const Vector&                    d_ssa,
    const Vector&                    d_g,
    const Vector&                    f_grid,
    std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const {
  check_sizes(extinction, ssa, g, f_grid);
  check_sizes(d_extinction, d_ssa, d_g, f_grid);
  const auto c = coefficients(extinction, ssa, f_grid);
  const auto d = coefficient_derivatives(extinction, ssa, d_extinction, d_ssa, f_grid);
  return tro_gridded_properties(
      d, ScatteredDerivative{c.scattering, d.scattering, g, d_g}, std::move(zenith_angle_grid));
}

ScatteringTroSpectralVector HenyeyGreensteinScatterer::tro_spectral(
    const Vector& extinction_, const Vector& ssa_, const Vector& g, const Vector& f_grid, Index l) const {
  check_sizes(extinction_, ssa_, g, f_grid);

  auto t_grid     = std::make_shared<Vector>(Vector{0.0});
  auto f_grid_ptr = std::make_shared<Vector>(f_grid);

  CompactPlanarMuelmatMatrix pm(f_grid.size(), l + 1, CompactPlanarMuelmat{});
  CompactPlanarPropmatVector emd(f_grid.size());
  CompactPlanarStokvecVector av(f_grid.size());

  constexpr Numeric inv_sphere = 0.5 * Constant::inv_sqrt_pi;  // sqrt(1/4pi)
  Vector            f(l + 4), f_pol(l + 1);

  for (Size f_ind = 0; f_ind < f_grid.size(); ++f_ind) {
    // The coefficients of p and of p f on Y_k0 = sqrt((2 k + 1) / 4 pi) P_k(cos(Theta)), from
    // p = sum_k sqrt((2 k + 1) / 4 pi) g^k Y_k0 to degree l + 3, which the cos^3 term of f needs
    for (Index k = 0; k <= l + 3; k++)
      f[k] = inv_sphere * std::sqrt(static_cast<Numeric>(2 * k + 1)) * std::pow(g[f_ind], k);
    const Vector c1 = times_cos(f), c3 = times_cos(times_cos(c1));
    for (Index k = 0; k <= l; k++) f_pol[k] = 0.5 * (3.0 * c1[k] - c3[k]);

    const Numeric extinction      = extinction_[f_ind];
    const auto    scattering_xsec = extinction * ssa_[f_ind];

    for (Index ind = 0; ind <= l; ++ind) {
      pm[f_ind, ind].F11() = f[ind] * scattering_xsec;
      pm[f_ind, ind].F22() = f[ind] * scattering_xsec;
      pm[f_ind, ind].F33() = f_pol[ind] * scattering_xsec;
      pm[f_ind, ind].F44() = f_pol[ind] * scattering_xsec;
    }

    emd[f_ind].A() = extinction;
    av[f_ind].I()  = extinction - scattering_xsec;
  }

  return {.phase_matrix = std::move(pm), .extinction_matrix = std::move(emd), .absorption_vector = std::move(av)};
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
HenyeyGreensteinScatterer::aro_gridded(const Vector&                    extinction,
                                       const Vector&                    ssa,
                                       const Vector&                    g,
                                       const Vector&                    f_grid,
                                       const Vector&                    za_inc_grid,
                                       const Vector&                    delta_aa_grid,
                                       std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  check_sizes(extinction, ssa, g, f_grid);
  const auto c = coefficients(extinction, ssa, f_grid);
  return aro_gridded_properties(c, Scattered{c.scattering, g}, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Gridded>
HenyeyGreensteinScatterer::aro_gridded_derivative(const Vector&                    extinction,
                                                  const Vector&                    ssa,
                                                  const Vector&                    g,
                                                  const Vector&                    d_extinction,
                                                  const Vector&                    d_ssa,
                                                  const Vector&                    d_g,
                                                  const Vector&                    f_grid,
                                                  const Vector&                    za_inc_grid,
                                                  const Vector&                    delta_aa_grid,
                                                  std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
  check_sizes(extinction, ssa, g, f_grid);
  check_sizes(d_extinction, d_ssa, d_g, f_grid);
  const auto c = coefficients(extinction, ssa, f_grid);
  const auto d = coefficient_derivatives(extinction, ssa, d_extinction, d_ssa, f_grid);
  return aro_gridded_properties(
      d, ScatteredDerivative{c.scattering, d.scattering, g, d_g}, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Spectral>
HenyeyGreensteinScatterer::aro_spectral(const Vector& extinction,
                                        const Vector& ssa,
                                        const Vector& g,
                                        const Vector& f_grid,
                                        const Vector& za_inc_grid,
                                        Index         degree,
                                        Index         order) const {
  auto sht = sht::provider.get_instance_lm(degree, order);
  return aro_gridded(extinction,
                     ssa,
                     g,
                     f_grid,
                     za_inc_grid,
                     *sht->get_aa_grid_ptr(),
                     std::make_shared<ZenithAngleGrid>(sht->get_zenith_angle_grid()))
      .to_spectral(degree, order);
}

BulkScatteringProperties<scattering::Format::ARO, scattering::Representation::Fourier>
HenyeyGreensteinScatterer::aro_fourier(const Vector& extinction,
                                       const Vector& ssa,
                                       const Vector& g,
                                       const Vector& f_grid,
                                       const Vector& za_inc_grid,
                                       const Vector& za_scat_grid,
                                       Index         max_mode) const {
  check_sizes(extinction, ssa, g, f_grid);
  const auto c      = coefficients(extinction, ssa, f_grid);
  auto       za_inc = std::make_shared<const Vector>(za_inc_grid);
  Matrix     integral(1, f_grid.size());  // p integrates to 1, so the phase integral is the scattering coefficient
  integral[0, joker] = c.scattering;
  auto phase =
      tro_lab_frame_fourier_modes(c.t_grid,
                                  c.f_grid,
                                  za_inc,
                                  std::make_shared<const ZenithAngleGrid>(IrregularZenithAngleGrid(za_scat_grid)),
                                  max_mode,
                                  integral,
                                  scattering_matrix(f_grid.size(), Scattered{c.scattering, g}));
  return {.phase_matrix      = std::move(phase),
          .extinction_matrix = c.extinction.to_lab_frame(za_inc).to_fourier(),
          .absorption_vector = c.absorption.to_lab_frame(za_inc).to_fourier()};
}

std::ostream& operator<<(std::ostream& os, const HenyeyGreensteinScatterer& scatterer) {
  return os << "HenyeyGreensteinScatterer(" << scatterer.extinction_name << ", " << scatterer.ssa_name << ", "
            << scatterer.g_name << ")";
}
}  // namespace scattering

void xml_io_stream<scattering::HenyeyGreensteinScatterer>::write(std::ostream&                                os,
                                                                 const scattering::HenyeyGreensteinScatterer& x,
                                                                 bofstream*                                   pbofs,
                                                                 std::string_view                             name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.extinction_name, pbofs);
  xml_write_to_stream(os, x.ssa_name, pbofs);
  xml_write_to_stream(os, x.g_name, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<scattering::HenyeyGreensteinScatterer>::read(std::istream&                          is,
                                                                scattering::HenyeyGreensteinScatterer& x,
                                                                bifstream*                             pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.extinction_name, pbifs);
  xml_read_from_stream(is, x.ssa_name, pbifs);
  xml_read_from_stream(is, x.g_name, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
