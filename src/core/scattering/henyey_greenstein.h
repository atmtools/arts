#ifndef HENYEY_GREENSTEIN_H_
#define HENYEY_GREENSTEIN_H_

#include <configtypes.h>

#include "bulk_scattering_properties.h"
#include "general_tro_spectral.h"
#include "properties.h"
#include "scattering_medium.h"

namespace scattering {
/** A Henyey-Greenstein scatterer with asymmetry parameter g
 *
 * The scattering matrix is F11 = F22 = p, F33 = F44 = p f and F12 = F34 = 0,
 * with p the Henyey-Greenstein phase function and
 * f = (3 cos(Theta) - cos^3(Theta)) / 2.  It produces no polarization, and it
 * is regular at forward and backward scattering, as a physical scattering
 * matrix is: F22 + F33 = p (1 + cos(Theta))^2 (2 - cos(Theta)) / 2 and
 * F22 - F33 = p (1 - cos(Theta))^2 (2 + cos(Theta)) / 2, so F33 = F22 at
 * 0 deg, F33 = -F22 and F44 = F11 - 2 F22 at 180 deg, and its expansion in
 * generalized spherical functions converges as fast as the Legendre series
 * of p.
 *
 * The explicit methods take the extinction [m^-1], the single scattering
 * albedo and the asymmetry parameter g, each at every frequency, with g in
 * [-1, 1].  The template wrappers take the point of a medium and its spectral
 * point, and read them from the spectral point as the spectral properties
 * named extinction_name, ssa_name and g_name.
 */
struct HenyeyGreensteinScatterer {
  String extinction_name{};
  String ssa_name{};
  String g_name{};

  HenyeyGreensteinScatterer() = default;
  HenyeyGreensteinScatterer(String extinction_name, String ssa_name, String g_name);

  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(
      const Vector&                    extinction,
      const Vector&                    ssa,
      const Vector&                    g,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  //! The derivative, from the derivatives of the extinction and the single scattering albedo at each frequency
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_derivative(
      const Vector&                    extinction,
      const Vector&                    ssa,
      const Vector&                    g,
      const Vector&                    d_extinction,
      const Vector&                    d_ssa,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const;

  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(
      const Vector& extinction, const Vector& ssa, const Vector& g, const Vector& f_grid, Index l) const;

  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
      const Vector&                    extinction,
      const Vector&                    ssa,
      const Vector&                    g,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  //! The derivative, from the derivatives of the extinction and the single scattering albedo at each frequency
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_derivative(
      const Vector&                    extinction,
      const Vector&                    ssa,
      const Vector&                    g,
      const Vector&                    d_extinction,
      const Vector&                    d_ssa,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  /** The SHT of degree and order of the laboratory-frame phase matrix over the scattering directions, per
   *  incidence zenith angle [deg]: the closed form sampled exactly on the SHT's own grid and transformed */
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Spectral> aro_spectral(const Vector& extinction,
                                                                                             const Vector& ssa,
                                                                                             const Vector& g,
                                                                                             const Vector& f_grid,
                                                                                             const Vector& za_inc_grid,
                                                                                             Index         degree,
                                                                                             Index order) const;

  /** The azimuthal Fourier modes m = 0..max_mode of the laboratory-frame phase matrix, exactly
   *
   * See tro_lab_frame_fourier_modes; the phase matrix is the closed form at
   * every scattering angle.  The extinction matrix and absorption vector are
   * at the incidence zenith angles.
   *
   * @param za_inc_grid The incidence zenith angles [deg]
   * @param za_scat_grid The scattering zenith angles [deg], ascending
   * @param max_mode The highest mode, >= 0
   */
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(const Vector& extinction,
                                                                                           const Vector& ssa,
                                                                                           const Vector& g,
                                                                                           const Vector& f_grid,
                                                                                           const Vector& za_inc_grid,
                                                                                           const Vector& za_scat_grid,
                                                                                           Index max_mode) const;

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(
      const Point&,
      const SpectralPoint&             spectral_point,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> zenith_angle_grid) const {
    return tro_gridded(spectral_property(spectral_point, extinction_name, f_grid.size()),
                       spectral_property(spectral_point, ssa_name, f_grid.size()),
                       spectral_property(spectral_point, g_name, f_grid.size()),
                       f_grid,
                       std::move(zenith_angle_grid));
  }

  //! The properties depend on the spectral point only, so their derivative by a key of the medium is zero
  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_derivative(
      const Point&,
      const SpectralPoint&             spectral_point,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> zenith_angle_grid,
      const medium_key<Point>&) const {
    const Vector zero(f_grid.size(), 0.0);
    return tro_gridded_derivative(spectral_property(spectral_point, extinction_name, f_grid.size()),
                                  spectral_property(spectral_point, ssa_name, f_grid.size()),
                                  spectral_property(spectral_point, g_name, f_grid.size()),
                                  zero,
                                  zero,
                                  f_grid,
                                  std::move(zenith_angle_grid));
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(const Point&,
                                                         const SpectralPoint& spectral_point,
                                                         const Vector&        f_grid,
                                                         Index                l) const {
    return tro_spectral(spectral_property(spectral_point, extinction_name, f_grid.size()),
                        spectral_property(spectral_point, ssa_name, f_grid.size()),
                        spectral_property(spectral_point, g_name, f_grid.size()),
                        f_grid,
                        l);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
      const Point&,
      const SpectralPoint&             spectral_point,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
    return aro_gridded(spectral_property(spectral_point, extinction_name, f_grid.size()),
                       spectral_property(spectral_point, ssa_name, f_grid.size()),
                       spectral_property(spectral_point, g_name, f_grid.size()),
                       f_grid,
                       za_inc_grid,
                       delta_aa_grid,
                       std::move(za_scat_grid));
  }

  //! The properties depend on the spectral point only, so their derivative by a key of the medium is zero
  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_derivative(
      const Point&,
      const SpectralPoint&             spectral_point,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const medium_key<Point>&) const {
    const Vector zero(f_grid.size(), 0.0);
    return aro_gridded_derivative(spectral_property(spectral_point, extinction_name, f_grid.size()),
                                  spectral_property(spectral_point, ssa_name, f_grid.size()),
                                  spectral_property(spectral_point, g_name, f_grid.size()),
                                  zero,
                                  zero,
                                  f_grid,
                                  za_inc_grid,
                                  delta_aa_grid,
                                  std::move(za_scat_grid));
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Spectral> aro_spectral(
      const Point&,
      const SpectralPoint& spectral_point,
      const Vector&        f_grid,
      const Vector&        za_inc_grid,
      Index                degree,
      Index                order) const {
    return aro_spectral(spectral_property(spectral_point, extinction_name, f_grid.size()),
                        spectral_property(spectral_point, ssa_name, f_grid.size()),
                        spectral_property(spectral_point, g_name, f_grid.size()),
                        f_grid,
                        za_inc_grid,
                        degree,
                        order);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(
      const Point&,
      const SpectralPoint& spectral_point,
      const Vector&        f_grid,
      const Vector&        za_inc_grid,
      const Vector&        za_scat_grid,
      Index                max_mode) const {
    return aro_fourier(spectral_property(spectral_point, extinction_name, f_grid.size()),
                       spectral_property(spectral_point, ssa_name, f_grid.size()),
                       spectral_property(spectral_point, g_name, f_grid.size()),
                       f_grid,
                       za_inc_grid,
                       za_scat_grid,
                       max_mode);
  }

  friend std::ostream& operator<<(std::ostream& os, const HenyeyGreensteinScatterer& scatterer);
};
}  // namespace scattering

template <> struct std::formatter<scattering::HenyeyGreensteinScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const scattering::HenyeyGreensteinScatterer& v, FmtContext& ctx) const {
    if (tags.names) { return tags.format(ctx, "HenyeyGreensteinScatterer"sv); }

    return tags.format(ctx, v.extinction_name, tags.sep(), v.ssa_name, tags.sep(), v.g_name);
  }
};

template <> struct xml_io_stream<scattering::HenyeyGreensteinScatterer> {
  static constexpr std::string_view type_name = "HenyeyGreensteinScatterer";

  static void write(std::ostream&,
                    const scattering::HenyeyGreensteinScatterer&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv);

  static void read(std::istream&, scattering::HenyeyGreensteinScatterer&, bifstream* = nullptr);
};

#endif  // HENYEY_GREENSTEIN_H_
