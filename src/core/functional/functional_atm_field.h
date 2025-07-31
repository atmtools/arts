#pragma once

#include <artstime.h>
#include <enumsFieldComponent.h>
#include <matpack.h>
#include <xml.h>

namespace Atm {
struct MagnitudeField {
  SortedGriddedField3 magnitude{};
  SortedGriddedField3 theta{};
  SortedGriddedField3 phi{};

  FieldComponent component{};

  Numeric operator()(Numeric, Numeric, Numeric) const;

  [[nodiscard]] ConstVectorView x() const;
  [[nodiscard]] VectorView x();
  [[nodiscard]] std::vector<std::pair<Index, Numeric>> w(Numeric alt,
                                                         Numeric lat,
                                                         Numeric lon) const;
};
}  // namespace Atm

template <>
struct xml_io_stream_name<Atm::MagnitudeField> {
  static constexpr std::string_view name = "MagnitudeField";
};

template <>
struct xml_io_stream_aggregate<Atm::MagnitudeField> {
  static constexpr bool value = true;
};
