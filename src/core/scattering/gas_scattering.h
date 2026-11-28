/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <atm.h>
#include <format_tags.h>
#include <matpack.h>

#include <variant>

#include "bulk_scattering_properties.h"
#include "general_tro_spectral.h"
#include "scattering_medium.h"

/*!
  \file   gas_scattering.h
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>,
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-03-15

  \brief  Header file for functions related to gas scattering.
*/

#ifndef gas_scattering_h
#define gas_scattering_h

CompactPlanarMuelmat calc_rayleighPhaMat(const Numeric& theta_rad, const Numeric& depolarization_factor = 0.0);

namespace scattering {

struct AirSimpleGasScattering {
  //! The scattering cross-section of air at a frequency [m^2]
  [[nodiscard]] Numeric operator()(Numeric frequency) const;
};

struct ConstantGasScattering {
  Numeric cross_section{};

  ConstantGasScattering() = default;
  explicit ConstantGasScattering(Numeric cross_section);

  //! The scattering cross-section [m^2]
  [[nodiscard]] Numeric operator()(Numeric frequency) const;
};

using GasScatteringCoefficient = std::variant<AirSimpleGasScattering, ConstantGasScattering>;

struct IsotropicGasScattering {};

struct RayleighGasScattering {
  Numeric depolarization_factor{};

  RayleighGasScattering() = default;
  explicit RayleighGasScattering(Numeric depolarization_factor);
};

using GasScatteringPhaseMatrix = std::variant<IsotropicGasScattering, RayleighGasScattering>;

/** A molecular gas-scattering species.
 *
 * The extinction coefficient and normalized phase matrix are deliberately
 * independent, matching the ARTS2 gas-scattering agenda contract.  Gas
 * scattering is conservative: its absorption vector is zero.
 */
struct GasScatterer {
  GasScatteringCoefficient coefficient{AirSimpleGasScattering{}};
  GasScatteringPhaseMatrix phase_matrix{RayleighGasScattering{}};

  GasScatterer() = default;
  GasScatterer(GasScatteringCoefficient coefficient, GasScatteringPhaseMatrix phase_matrix);

  /* The explicit methods take the number density of the gas [m^-3]: every property is the cross-section times it.
     The template wrappers take the point of a medium and its number density from gas_number_density.  Gas
     scattering has no spectral properties, so they do not read the spectral point. */

  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(
      Numeric number_density, const Vector& f_grid, std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(Numeric       number_density,
                                                         const Vector& f_grid,
                                                         Index         degree) const;

  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
      Numeric                          number_density,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const;

  /** The SHT of degree and order of the laboratory-frame phase matrix over the scattering directions, per
   *  incidence zenith angle [deg]: the closed form sampled exactly on the SHT's own grid and transformed */
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Spectral> aro_spectral(
      Numeric number_density, const Vector& f_grid, const Vector& za_inc_grid, Index degree, Index order) const;

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
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(Numeric       number_density,
                                                                                           const Vector& f_grid,
                                                                                           const Vector& za_inc_grid,
                                                                                           const Vector& za_scat_grid,
                                                                                           Index max_mode) const;

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
    return tro_gridded(gas_number_density(point), f_grid, std::move(za_scat_grid));
  }

  //! The properties are proportional to the number density, so their derivative is them times its relative derivative
  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_derivative(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const medium_key<Point>&         target) const {
    auto out  = tro_gridded(gas_number_density(point), f_grid, std::move(za_scat_grid));
    out      *= gas_number_density_relative_derivative(point, target);
    return out;
  }

  //! Gas scattering reads no spectral property, so its derivative by one is zero
  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::TRO, Representation::Gridded> tro_gridded_derivative(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const spectral_medium_key<Point>&) const {
    auto out  = tro_gridded(gas_number_density(point), f_grid, std::move(za_scat_grid));
    out      *= 0.0;
    return out;
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] ScatteringTroSpectralVector tro_spectral(const Point& point,
                                                         const SpectralPoint&,
                                                         const Vector& f_grid,
                                                         Index         degree) const {
    return tro_spectral(gas_number_density(point), f_grid, degree);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid) const {
    return aro_gridded(gas_number_density(point), f_grid, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
  }

  //! The properties are proportional to the number density, so their derivative is them times its relative derivative
  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_derivative(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const medium_key<Point>&         target) const {
    auto out  = aro_gridded(gas_number_density(point), f_grid, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
    out      *= gas_number_density_relative_derivative(point, target);
    return out;
  }

  //! Gas scattering reads no spectral property, so its derivative by one is zero
  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Gridded> aro_gridded_derivative(
      const Point& point,
      const SpectralPoint&,
      const Vector&                    f_grid,
      const Vector&                    za_inc_grid,
      const Vector&                    delta_aa_grid,
      std::shared_ptr<ZenithAngleGrid> za_scat_grid,
      const spectral_medium_key<Point>&) const {
    auto out  = aro_gridded(gas_number_density(point), f_grid, za_inc_grid, delta_aa_grid, std::move(za_scat_grid));
    out      *= 0.0;
    return out;
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Spectral> aro_spectral(const Point& point,
                                                                                             const SpectralPoint&,
                                                                                             const Vector& f_grid,
                                                                                             const Vector& za_inc_grid,
                                                                                             Index         degree,
                                                                                             Index order) const {
    return aro_spectral(gas_number_density(point), f_grid, za_inc_grid, degree, order);
  }

  template <medium_point Point, spectral_medium_point<Point> SpectralPoint>
  [[nodiscard]] BulkScatteringProperties<Format::ARO, Representation::Fourier> aro_fourier(const Point& point,
                                                                                           const SpectralPoint&,
                                                                                           const Vector& f_grid,
                                                                                           const Vector& za_inc_grid,
                                                                                           const Vector& za_scat_grid,
                                                                                           Index max_mode) const {
    return aro_fourier(gas_number_density(point), f_grid, za_inc_grid, za_scat_grid, max_mode);
  }
};

}  // namespace scattering

template <> struct std::formatter<scattering::AirSimpleGasScattering> : std::formatter<std::string_view> {
  template <class FmtContext>
  FmtContext::iterator format(const scattering::AirSimpleGasScattering&, FmtContext& ctx) const {
    return std::formatter<std::string_view>::format("AirSimpleGasScattering", ctx);
  }
};

template <> struct std::formatter<scattering::ConstantGasScattering> {
  format_tags    tags;
  constexpr auto parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class FmtContext>
  FmtContext::iterator format(const scattering::ConstantGasScattering& x, FmtContext& ctx) const {
    return tags.format(ctx, x.cross_section);
  }
};

template <> struct std::formatter<scattering::IsotropicGasScattering> : std::formatter<std::string_view> {
  template <class FmtContext>
  FmtContext::iterator format(const scattering::IsotropicGasScattering&, FmtContext& ctx) const {
    return std::formatter<std::string_view>::format("IsotropicGasScattering", ctx);
  }
};

template <> struct std::formatter<scattering::RayleighGasScattering> {
  format_tags    tags;
  constexpr auto parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class FmtContext>
  FmtContext::iterator format(const scattering::RayleighGasScattering& x, FmtContext& ctx) const {
    return tags.format(ctx, x.depolarization_factor);
  }
};

template <> struct std::formatter<scattering::GasScatterer> {
  format_tags    tags;
  constexpr auto parse(std::format_parse_context& ctx) { return parse_format_tags(tags, ctx); }
  template <class FmtContext> FmtContext::iterator format(const scattering::GasScatterer& x, FmtContext& ctx) const {
    return tags.format(ctx, x.coefficient, tags.sep(), x.phase_matrix);
  }
};

#define GAS_SCATTERING_XML_STREAM(Type)                                                                       \
  template <> struct xml_io_stream<scattering::Type> {                                                        \
    static constexpr std::string_view type_name = #Type;                                                      \
    static void write(std::ostream&, const scattering::Type&, bofstream* = nullptr, std::string_view = ""sv); \
    static void read(std::istream&, scattering::Type&, bifstream* = nullptr);                                 \
  }

GAS_SCATTERING_XML_STREAM(AirSimpleGasScattering);
GAS_SCATTERING_XML_STREAM(ConstantGasScattering);
GAS_SCATTERING_XML_STREAM(IsotropicGasScattering);
GAS_SCATTERING_XML_STREAM(RayleighGasScattering);
GAS_SCATTERING_XML_STREAM(GasScatterer);

#undef GAS_SCATTERING_XML_STREAM

#endif /* gas_scattering_h */
