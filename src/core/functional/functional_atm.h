#pragma once

#include "functional_atm_field.h"
#include "functional_atm_mag_field.h"
#include "functional_atm_pressure_field.h"

namespace Atm {
struct External {
  bool is_field_function = false;
  std::function<Numeric(Numeric, Numeric, Numeric)> f;
  std::function<std::span<Numeric>()> mx;
  std::function<std::span<const Numeric>()> cx;
  std::function<std::vector<std::pair<Index, Numeric>>(Numeric, Numeric, Numeric)> cw;

  [[nodiscard]] Numeric operator()(Numeric, Numeric, Numeric) const;
  [[nodiscard]] VectorView x();
  [[nodiscard]] ConstVectorView x() const;
  [[nodiscard]] std::vector<std::pair<Index, Numeric>> w(Numeric,
                                                         Numeric,
                                                         Numeric) const;
};
};  // namespace Atm

template <>
struct xml_io_stream_name<Atm::External> {
  static constexpr std::string_view name = "External";
};

template <>
struct xml_io_stream_aggregate<Atm::External> {
  static constexpr bool value = true;
};
