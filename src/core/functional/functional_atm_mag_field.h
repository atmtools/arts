#pragma once

#include <artstime.h>
#include <enumsFieldComponent.h>
#include <matpack.h>
#include <xml.h>

namespace Atm {
struct IGRF13 {
  Time time{};
  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;
};

struct SchmidthLegendre {
  Matrix gh;
  Size N;
  Numeric r0{};
  Vector2 ell{};
  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;

  [[nodiscard]] ConstVectorView x() const;
  [[nodiscard]] VectorView x();
  [[nodiscard]] std::vector<std::pair<Index, Numeric>> w(Numeric alt,
                                                         Numeric lat,
                                                         Numeric lon) const;
};

SchmidthLegendre from(const IGRF13&);
}  // namespace Atm

template <>
struct xml_io_stream_name<Atm::IGRF13> {
  static constexpr std::string_view name = "IGRF13";
};

template <>
struct xml_io_stream_aggregate<Atm::IGRF13> {
  static constexpr bool value = true;
};

template <>
struct xml_io_stream_name<Atm::SchmidthLegendre> {
  static constexpr std::string_view name = "SchmidthLegendre";
};

template <>
struct xml_io_stream_aggregate<Atm::SchmidthLegendre> {
  static constexpr bool value = true;
};
