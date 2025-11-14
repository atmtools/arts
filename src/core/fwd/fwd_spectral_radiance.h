#pragma once

#include <atm.h>
#include <fwd_path.h>
#include <fwd_propmat.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <surf.h>
#include <xml.h>

#include <memory>

namespace fwd {
struct spectral_rad {
  AscendingGrid alt;
  LatGrid lat;
  LonGrid lon;
  matpack::data_t<std::shared_ptr<AtmPoint>, 3> atm;
  matpack::data_t<propmat, 3> pm;

  matpack::data_t<std::function<Stokvec(Numeric, Vector2)>, 2>
      spectral_rad_surface;
  matpack::data_t<std::function<Stokvec(Numeric, Vector2)>, 2>
      spectral_rad_space;

  Vector2 ellipsoid;

  struct as_vector {};
  struct weighted_position {
    Numeric w{0.};
    Index i{0}, j{0}, k{0};
  };

  spectral_rad();
  spectral_rad(const spectral_rad&);
  spectral_rad(spectral_rad&&) noexcept;
  spectral_rad& operator=(const spectral_rad&);
  spectral_rad& operator=(spectral_rad&&) noexcept;

  spectral_rad(AscendingGrid alt,
               LatGrid lat,
               LonGrid lon,
               const AtmField& atm,
               const SurfaceField& surf,
               const std::shared_ptr<AbsorptionBands>& lines,
               const std::shared_ptr<ArrayOfCIARecord>& cia,
               const std::shared_ptr<ArrayOfXsecRecord>& xsec,
               const std::shared_ptr<PredefinedModelData>& predef,
               Numeric ciaextrap = {},
               Index ciarobust   = {});

  Stokvec operator()(const Numeric f,
                     const std::vector<path>& path_points,
                     const Numeric cutoff_transmission = 1e-6) const;

  StokvecVector operator()(const Numeric f,
                           const std::vector<path>& path_points,
                           spectral_rad::as_vector) const;

  [[nodiscard]] const AscendingGrid& altitude() const { return alt; }
  [[nodiscard]] const LatGrid& latitude() const { return lat; }
  [[nodiscard]] const LonGrid& longitude() const { return lon; }

  [[nodiscard]] std::vector<path> geometric_planar(const Vector3 pos,
                                                   const Vector2 los) const;
  [[nodiscard]] std::vector<path> from_path(
      const ArrayOfPropagationPathPoint& propagation_path) const;
  void from_path(std::vector<path>& out,
                 const ArrayOfPropagationPathPoint& propagation_path) const;

  [[nodiscard]] std::array<weighted_position, 8> pos_weights(
      const path& pp) const;

  [[nodiscard]] Stokvec B(const Numeric f,
                          const std::array<weighted_position, 8>& pos) const;

  [[nodiscard]] Stokvec Iback(const Numeric f,
                              const std::array<weighted_position, 8>& pos,
                              const path& pp) const;

  [[nodiscard]] std::pair<Propmat, Stokvec> PM(
      const Numeric f,
      const std::array<weighted_position, 8>& pos,
      const path& pp) const;
};
}  // namespace fwd

using SpectralRadianceOperator = fwd::spectral_rad;

template <>
struct std::formatter<SpectralRadianceOperator> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SpectralRadianceOperator& v,
                              FmtContext& ctx) const {
    return tags.format(ctx,
                       "operator over\n  alt: "sv,
                       v.alt,
                       "\n  lat: "sv,
                       v.lat,
                       "\n  lon: "sv,
                       v.lon);
  }
};

template <>
struct xml_io_stream<SpectralRadianceOperator> {
  static constexpr std::string_view type_name = "SpectralRadianceOperator"sv;

  static void write(std::ostream& os,
                    const SpectralRadianceOperator& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SpectralRadianceOperator& x,
                   bifstream* pbifs = nullptr);
};
