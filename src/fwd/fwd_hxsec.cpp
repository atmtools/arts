#include "fwd_hxsec.h"

#include "physics_funcs.h"

namespace fwd::hxsec {
single::single(Numeric p,
               Numeric t,
               Numeric VMR,
               const std::shared_ptr<XsecRecord>& cia,
               Verbosity verb)
    : scl{number_density(p, t) * VMR},
      P{p},
      T{t},
      verbosity{verb},
      xsecrec{cia} {}

Complex single::at(Numeric f) const {
  Numeric out{};
  xsecrec->Extract(ExhaustiveVectorView{out}, Vector{f}, P, T, verbosity);
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

full::full(Numeric p,
       Numeric t,
       const Vector& vmrs,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::vector<std::shared_ptr<XsecRecord>>& xsec,
       Verbosity verb) {
  for (auto& specs : allspecs) {
    for (auto& spec : specs) {
      if (spec.type == Species::TagType::Cia) {
        const auto data = hitran_xsec_get_data(xsec, spec.Spec());
        ARTS_USER_ERROR_IF(not data, "Cannot find XSEC data for tag: ", spec)

        const Numeric VMR =
            vmrs[find_first_species(allspecs, spec.Spec())];

        models.emplace_back(p, t, VMR, data, verb);
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
