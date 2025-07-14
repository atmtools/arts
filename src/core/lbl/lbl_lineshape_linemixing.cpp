#include "lbl_lineshape_linemixing.h"

namespace lbl::linemixing {
Numeric species_data::Q(const Rational J,
                        const Numeric T,
                        const Numeric T0,
                        const Numeric energy) const {
  return std::exp(-beta(T0, T) * energy / (Constant::k * T)) * scaling(T0, T) /
         pow(J * (J + 1), lambda(T0, T));
}

Numeric species_data::Omega(const Numeric T,
                            const Numeric T0,
                            const Numeric mass,
                            const Numeric other_mass,
                            const Numeric energy_x,
                            const Numeric energy_xm2) const {
  using Constant::h;
  using Constant::h_bar;
  using Constant::k;
  using Constant::m_u;
  using Constant::pi;
  using Math::pow2;

  // Constants for the expression
  constexpr Numeric fac = 8 * k / (m_u * pi);

  const Numeric wnnm2 = (energy_x - energy_xm2) / h_bar;

  const Numeric inv_eff_mass = 1 / mass + 1 / other_mass;
  const Numeric v_bar_pow2   = fac * T * inv_eff_mass;
  const Numeric tauc_pow2    = pow2(collisional_distance(T0, T)) / v_bar_pow2;

  return 1.0 / pow2(1 + pow2(wnnm2) * tauc_pow2 / 24.0);
}
}  // namespace lbl::linemixing

void xml_io_stream<LinemixingSingleEcsData>::write(
    std::ostream &os,
    const LinemixingSingleEcsData &x,
    bofstream *pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.scaling, pbofs);
  xml_write_to_stream(os, x.beta, pbofs);
  xml_write_to_stream(os, x.lambda, pbofs);
  xml_write_to_stream(os, x.collisional_distance, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<LinemixingSingleEcsData>::read(std::istream &is,
                                                  LinemixingSingleEcsData &x,
                                                  bifstream *pbifs) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.scaling, pbifs);
  xml_read_from_stream(is, x.beta, pbifs);
  xml_read_from_stream(is, x.lambda, pbifs);
  xml_read_from_stream(is, x.collisional_distance, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
} catch (const std::exception &e) {
  throw std::runtime_error(
      std::format("Error reading {}:\n{}", type_name, e.what()));
}
