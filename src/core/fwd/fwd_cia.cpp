#include "fwd_cia.h"

#include <arts_constexpr_math.h>
#include <physics_funcs.h>

#include <algorithm>
#include <functional>
#include <numeric>

#include "debug.h"

namespace fwd::cia {
full::single::single(Numeric p,
                     Numeric t,
                     Numeric VMR1,
                     Numeric VMR2,
                     CIARecord* cia,
                     Numeric extrap,
                     Index robust)
    : scl(VMR1 * VMR2 * Math::pow2(number_density(p, t))),
      T(t),
      extrapol(extrap),
      ignore_errors(robust),
      ciarecords(cia) {}

Complex full::single::at(Numeric f) const {
  return scl * ciarecords->Extract(f, T, extrapol, ignore_errors);
}

void full::single::at(ExhaustiveComplexVectorView abs, const Vector& fs) const {
  std::transform(
      fs.begin(),
      fs.end(),
      abs.begin(),
      abs.begin(),
      [this](const Numeric& f, const Complex& s) { return s + at(f); });
}

ComplexVector full::single::at(const Vector& fs) const {
  ComplexVector abs(fs.size());
  at(abs, fs);
  return abs;
}

full::full(const AtmPoint& atm_point,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const std::shared_ptr<std::vector<CIARecord>>& cia,
           Numeric extrap,
           Index robust)
    : ciarecords(cia) {
  for (auto& specs : allspecs) {
    for (auto& spec : specs) {
      if (spec.type == Species::TagType::Cia) {
        auto* data = cia_get_data(cia, spec.Spec(), spec.cia_2nd_species);
        ARTS_USER_ERROR_IF(not data, "Cannot find CIA data for tag: ", spec)

        const Numeric VMR1 = atm_point[data->Species(0)];
        const Numeric VMR2 = atm_point[data->Species(1)];

        models.emplace_back(atm_point.pressure,
                            atm_point.temperature,
                            VMR1,
                            VMR2,
                            data,
                            extrap,
                            robust);
      }
    }
  }
}

Complex full::operator()(Numeric f) const {
  return std::transform_reduce(
      models.begin(), models.end(), Complex{}, std::plus<>{}, [f](auto& mod) {
        return mod.at(f);
      });
}

void full::operator()(ExhaustiveComplexVectorView abs, const Vector& fs) const {
  for (auto& mod : models) {
    mod.at(abs, fs);
  }
}

ComplexVector full::operator()(const Vector& fs) const {
  ComplexVector abs(fs.size());
  operator()(abs, fs);
  return abs;
}
}  // namespace fwd::cia
