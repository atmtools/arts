#include "functional_atm_pressure_field.h"

#include "functional_atm_field_interp.h"

namespace Atm {
Numeric HydrostaticPressure::step(Numeric p, Numeric h, Numeric d) const {
  switch (option) {
    using enum HydrostaticPressureOption;
    case HypsometricEquation: return p * std::exp(-h * d);
    case HydrostaticEquation: return std::max<Numeric>(0, p * (1.0 - h * d));
  }
  return NAN;
}

HydrostaticPressure::HydrostaticPressure(Tensor3 in_grad_p,
                                         const SortedGriddedField2& pre0,
                                         AscendingGrid in_alt,
                                         HydrostaticPressureOption option_)
    : grad_p{.data_name  = "Pressure Gradient",
             .data       = std::move(in_grad_p),
             .grid_names = {"Altitude", "Latitude", "Longitude"},
             .grids      = {std::move(in_alt), pre0.grid<0>(), pre0.grid<1>()}},
      pre(grad_p),
      option{option_} {
  pre.data_name         = "Pressure";
  pre[0]                = pre0.data;
  auto& [alt, lat, lon] = grad_p.grids;

  for (Size i = 1; i < alt.size(); i++) {
    for (Size j = 0; j < lat.size(); j++) {
      for (Size k = 0; k < lon.size(); k++) {
        const Numeric h  = alt[i] - alt[i - 1];
        const Numeric p0 = pre[i - 1, j, k];
        const Numeric d0 = grad_p[i - 1, j, k];

        pre[i, j, k] = step(p0, h, d0);
      }
    }
  }
}

std::pair<Size, Numeric> HydrostaticPressure::find_alt(Numeric al) const {
  auto& alt  = std::get<0>(grad_p.grids);
  Size i     = std::distance(alt.begin(), std::ranges::upper_bound(alt, al));
  i         -= (i == alt.size());
  while (i > 0 and alt[i] > al) i--;
  return {i, al - alt[i]};
}

std::pair<Numeric, Numeric> HydrostaticPressure::level(Index alt_ind,
                                                       Numeric la,
                                                       Numeric lo) const {
  const auto latlag = interp::latlag(pre.grid<1>(), la);
  const auto lonlag = interp::lonlag(pre.grid<2>(), lo);
  const auto iw     = interp::interpweights(latlag, lonlag);
  const Numeric p   = interp::get(pre, alt_ind, iw, latlag, lonlag);
  const Numeric d   = interp::get(grad_p, alt_ind, iw, latlag, lonlag);
  return {p, d};
}

Numeric HydrostaticPressure::operator()(Numeric al,
                                        Numeric la,
                                        Numeric lo) const {
  const auto [i, h] = find_alt(al);
  const auto [p, d] = level(i, la, lo);
  return step(p, h, d);
}
}  // namespace Atm

void xml_io_stream<Atm::HydrostaticPressure>::write(
    std::ostream& os,
    const Atm::HydrostaticPressure& a,
    bofstream* pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, a.grad_p, pbofs);
  xml_write_to_stream(os, a.pre, pbofs);
  xml_write_to_stream(os, a.option, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Atm::HydrostaticPressure>::read(std::istream& is,
                                                   Atm::HydrostaticPressure& a,
                                                   bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, a.grad_p, pbifs);
  xml_read_from_stream(is, a.pre, pbifs);
  xml_read_from_stream(is, a.option, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
