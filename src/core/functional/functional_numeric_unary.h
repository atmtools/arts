#pragma once

#include <matpack.h>
#include <xml.h>
#include <xml_io_stream_functional.h>

#include <functional>
#include <variant>

//! Water equivalent saturation pressure [Pa] from temperature [K], by Murphy and Koop (2005)
struct WaterEquivalentPressureMK05 {
  //! Use the pressure over liquid water also below 0 C, otherwise over ice
  bool only_liquid{false};

  Numeric operator()(Numeric t) const;
};

template <> struct xml_io_stream_name<WaterEquivalentPressureMK05> {
  static constexpr std::string_view name = "WaterEquivalentPressureMK05"sv;
};

template <> struct xml_io_stream_aggregate<WaterEquivalentPressureMK05> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_functional<std::function<Numeric(Numeric)>> {
  using func_t    = Numeric (*)(Numeric);
  using structs_t = std::variant<WaterEquivalentPressureMK05>;
  static constexpr std::array<func_t*, 0> funcs{};
};
