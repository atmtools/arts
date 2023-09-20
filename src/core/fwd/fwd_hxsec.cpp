#include "fwd_hxsec.h"

#include "physics_funcs.h"

namespace fwd::hxsec {
single::single(Numeric p,
               Numeric t,
               Numeric VMR,
               const std::shared_ptr<XsecRecord>& cia)
    : scl{number_density(p, t) * VMR},
      P{p},
      T{t},
      xsecrec{cia} {}

Complex single::at(Numeric f) const {
  Numeric out{};
  xsecrec->Extract(ExhaustiveVectorView{out}, Vector{f}, P, T);
  return scl * out;
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

full::full(
       const AtmPoint& atm_point,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::vector<std::shared_ptr<XsecRecord>>& xsec) {
  for (auto& specs : allspecs) {
    for (auto& spec : specs) {
      if (spec.type == Species::TagType::XsecFit) {
        const auto data = hitran_xsec_get_data(xsec, spec.Spec());
        ARTS_USER_ERROR_IF(not data, "Cannot find XSEC data for tag: ", spec)

        const Numeric VMR = atm_point[spec.Spec()];

        models.emplace_back(atm_point.pressure, atm_point.temperature, VMR, data);
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
}  // namespace fwd::hxsec
