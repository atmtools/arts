#pragma once

#include <artstime.h>
#include <enumsFieldComponent.h>
#include <matpack.h>
#include <xml.h>

namespace Atm {
struct IGRF {
  Time           time{};
  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;
};

struct SchmidthLegendre {
  /** The g and h coefficients that affect the field component, in the order of the IGRF model
   *
   * The order is g_1^0, g_1^1, h_1^1, g_2^0, g_2^1, h_2^1, g_2^2, h_2^2, and so on to
   * degree N - 1.  The g_0^0, h_0^0, and h_n^0 coefficients do not affect the field and
   * are not kept, so there are N * N - 1 coefficients.  The east component, u, does not
   * depend on the g_n^0 either, so it has N * (N - 1) coefficients.
   */
  Vector gh;

  //! The number of rows of the triangular g and h matrices, the largest degree plus one
  Size N;

  Numeric        r0{};
  Vector2        ell{};
  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;

  //! The g and h coefficients as N x N lower triangular matrices, NaN where not in gh, throws if gh has the wrong size
  [[nodiscard]] std::pair<TriangularMatrix, TriangularMatrix> coefficients() const;

  [[nodiscard]] ConstVectorView                        x() const;
  [[nodiscard]] VectorView                             x();
  [[nodiscard]] std::vector<std::pair<Index, Numeric>> w(Numeric alt, Numeric lat, Numeric lon) const;
};

SchmidthLegendre from(const IGRF&);
}  // namespace Atm

template <> struct xml_io_stream_name<Atm::IGRF> {
  static constexpr std::string_view name = "IGRF";
};

template <> struct xml_io_stream_aggregate<Atm::IGRF> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_name<Atm::SchmidthLegendre> {
  static constexpr std::string_view name = "SchmidthLegendre";
};

template <> struct xml_io_stream_aggregate<Atm::SchmidthLegendre> {
  static constexpr bool value = true;
};
