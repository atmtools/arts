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
  //! The g (first row) and h (second row) coefficients of N x N lower triangular matrices, see matpack::triangular_index
  Matrix         gh;
  Size           N;
  Numeric        r0{};
  Vector2        ell{};
  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;

  //! The g and h coefficients as triangular matrices, throws if gh does not hold two N x N lower triangular matrices
  [[nodiscard]] std::pair<ConstTriangularMatrixView, ConstTriangularMatrixView> coefficients() const;

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
