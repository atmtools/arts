#pragma once

#include <artstime.h>
#include <enumsFieldComponent.h>
#include <matpack.h>
#include <xml.h>

namespace Atm {
struct IGRF13 {
  static constexpr Vector2 ell{6378137., 6356752.314245};

  Time time{};
  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;
};
}  // namespace Atm

template <>
struct xml_io_stream_name<Atm::IGRF13> {
  static constexpr std::string_view name = "IGRF13";
};

template <>
struct xml_io_stream_aggregate<Atm::IGRF13> {
  static constexpr bool value = true;
};
