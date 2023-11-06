#include "lbl_lineshape_voigt.h"

#include "arts_constants.h"
#include "atm.h"
#include "isotopologues.h"
#include "lbl_data.h"
#include "partfun.h"
#include "quantum_numbers.h"

namespace lbl::voigt {
Complex line_strength_calc(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const Complex lm{1 + line.ls.G(atm), -line.ls.Y(atm)};

  const auto n = atm[spec] * atm[spec.spec];

  return n * lm * s;
}

Complex line_strength_calc(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm,
                           const Size ispec) {
  const auto s =
      line.s(atm.temperature, PartitionFunctions::Q(atm.temperature, spec));

  const auto& ls = line.ls.single_models[ispec];

  const Complex lm{1 + ls.G(line.ls.T0, atm.temperature, atm.pressure),
                   -ls.Y(line.ls.T0, atm.temperature, atm.pressure)};

  const auto n = atm[spec] * atm[spec.spec] * atm[ls.species];

  return n * lm * s;
}

Numeric line_center_calc(const line& line, const AtmPoint& atm) {
  return line.f0 + line.ls.D0(atm) + line.ls.DV(atm);
}

Numeric line_center_calc(const line& line, const AtmPoint& atm, Size ispec) {
  return line.f0 +
         line.ls.single_models[ispec].D0(
             line.ls.T0, atm.temperature, atm.pressure) +
         line.ls.single_models[ispec].DV(
             line.ls.T0, atm.temperature, atm.pressure);
}

namespace lte {
single_shape::single_shape(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm)
    : s(line_strength_calc(spec, line, atm)),
      f0(line_center_calc(line, atm)),
      inv_gd(1.0 / (std::sqrt(Constant::doppler_broadening_const_squared *
                              atm.temperature / spec.mass) *
                    f0)),
      z_imag(line.ls.G0(atm) * inv_gd) {}

single_shape::single_shape(const SpeciesIsotopeRecord& spec,
                           const line& line,
                           const AtmPoint& atm,
                          const  Size ispec)
    : s(line_strength_calc(spec, line, atm, ispec)),
      f0(line_center_calc(line, atm, ispec)),
      inv_gd(1.0 / (std::sqrt(Constant::doppler_broadening_const_squared *
                              atm.temperature / spec.mass) *
                    f0)),
      z_imag(line.ls.single_models[ispec].G0(
                 line.ls.T0, atm.temperature, atm.pressure) *
             inv_gd) {}

void single_shape::as_zeeman(const line& line, zeeman::pol pol, Index iz) {
  s *= line.z.Strength(line.qn.val, pol, iz);
  f0 += line.z.Splitting(line.qn.val, pol, iz);
}
}  // namespace lte
}  // namespace lbl::voigt
