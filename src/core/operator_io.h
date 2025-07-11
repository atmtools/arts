#pragma once

#include <atm.h>
#include <gravity.h>

template <>
struct xml_io_stream_name<std::function<Numeric(Numeric, Numeric, Numeric)>> {
  constexpr static std::string_view name = "NumericTernary"sv;
};

template <>
struct xml_io_stream_functional<
    std::function<Numeric(Numeric, Numeric, Numeric)>> {
  using func_t = Numeric (*)(Numeric, Numeric, Numeric);
  using structs_t =
      std::variant<Atm::IGRF13, Atm::HydrostaticPressure, Gravity>;
  static constexpr std::array<func_t*, 0> funcs{};
};
