#pragma once

#include <artstime.h>
#include <enumsFieldComponent.h>
#include <matpack.h>
#include <xml.h>

namespace Atm {
struct HydrostaticPressure {
  SortedGriddedField3 grad_p{};
  SortedGriddedField3 pre{};
  HydrostaticPressureOption option{};

  [[nodiscard]] Numeric step(Numeric p, Numeric h, Numeric d) const;

  HydrostaticPressure()                                          = default;
  HydrostaticPressure(const HydrostaticPressure&)                = default;
  HydrostaticPressure(HydrostaticPressure&&) noexcept            = default;
  HydrostaticPressure& operator=(const HydrostaticPressure&)     = default;
  HydrostaticPressure& operator=(HydrostaticPressure&&) noexcept = default;

  HydrostaticPressure(Tensor3 in_grad_p,
                      const SortedGriddedField2& pre0,
                      AscendingGrid in_alt,
                      HydrostaticPressureOption option);

  [[nodiscard]] std::pair<Size, Numeric> find_alt(Numeric al) const;

  [[nodiscard]] std::pair<Numeric, Numeric> level(Index alt_ind,
                                                  Numeric la,
                                                  Numeric lo) const;

  Numeric operator()(Numeric al, Numeric la, Numeric lo) const;
};
}  // namespace Atm

template <>
struct xml_io_stream_name<Atm::HydrostaticPressure> {
  static constexpr std::string_view name = "HydrostaticPressure";
};

template <>
struct xml_io_stream<Atm::HydrostaticPressure> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<Atm::HydrostaticPressure>;

  static void write(std::ostream& os,
                    const Atm::HydrostaticPressure& a,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   Atm::HydrostaticPressure& a,
                   bifstream* pbifs = nullptr);
};
