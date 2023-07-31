#include "fwd_cia.h"

#include <arts_constexpr_math.h>
#include <physics_funcs.h>

#include <algorithm>
#include <functional>
#include <numeric>

#include "debug.h"

namespace fwd::cia {
single::single(Numeric p,
               Numeric t,
               Numeric VMR1,
               Numeric VMR2,
               const std::shared_ptr<CIARecord>& cia,
               Numeric extrap,
               Index robust)
    : scl(VMR1 * VMR2 * Math::pow2(number_density(p, t))),
      T(t),
      extrapol(extrap),
      ignore_errors(robust),
      ciarecords(cia) {}

Complex single::at(Numeric f) const {
  return scl * ciarecords->Extract(f, T, extrapol, ignore_errors);
}

void single::at(ExhaustiveComplexVectorView abs, const Vector& fs) const {
  std::transform(fs.begin(), fs.end(), abs.begin(), [this](const auto& f) {
    return at(f);
  });
}

ComplexVector single::at(const Vector& fs) const {
  ComplexVector abs(fs.size());
  at(abs, fs);
  return abs;
}

full::full(Numeric p,
           Numeric t,
           const Vector& vmrs,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const std::vector<std::shared_ptr<CIARecord>>& cia,
           Numeric extrap,
           Index robust) {
  for (auto& specs : allspecs) {
    for (auto& spec : specs) {
      if (spec.type == Species::TagType::Cia) {
        const auto data = cia_get_data(cia, spec.Spec(), spec.cia_2nd_species);
        ARTS_USER_ERROR_IF(not data, "Cannot find CIA data for tag: ", spec)

        const Numeric VMR1 =
            vmrs[find_first_species(allspecs, data->Species(0))];
        const Numeric VMR2 =
            vmrs[find_first_species(allspecs, data->Species(1))];

        models.emplace_back(p, t, VMR1, VMR2, data, extrap, robust);
      }
    }
  }
}

Complex full::at(Numeric f) const {
  return std::transform_reduce(
      models.begin(), models.end(), Complex{}, std::plus<>{}, [f](auto& mod) {
        return mod.at(f);
      });
}

void full::at(ExhaustiveComplexVectorView abs, const Vector& fs) const {
  for (auto& mod : models) {
    mod.at(abs, fs);
  }
}

ComplexVector full::at(const Vector& fs) const {
  ComplexVector abs(fs.size());
  at(abs, fs);
  return abs;
}
}  // namespace fwd::cia
