#include "fwd_hxsec.h"

#include "physics_funcs.h"

namespace fwd::hxsec {
full::single::single(Numeric p, Numeric t, Numeric VMR, XsecRecord* xsec)
    : scl{number_density(p, t) * VMR}, P{p}, T{t}, xsecrec{xsec} {}

Complex full::single::at(const Numeric frequency) const {
  Numeric out{};
  xsecrec->Extract(ExhaustiveVectorView{out}, Vector{frequency}, P, T);
  return scl * out;
}

void full::adapt() {
  models.resize(0);

  if (not xsecrec) {
    return;
  }

  ARTS_USER_ERROR_IF(not atm, "Must have an atmosphere")

  models.reserve(xsecrec->size());
  for (auto& model : *xsecrec) {
    models.emplace_back(atm->pressure,
                        atm->temperature,
                        atm->operator[](model.Species()),
                        &model);
  }
}

full::full(std::shared_ptr<AtmPoint> atm_,
           std::shared_ptr<ArrayOfXsecRecord> xsecrec_)
    : atm(std::move(atm_)), xsecrec(std::move(xsecrec_)) {
  adapt();
}

Complex full::operator()(const Numeric frequency) const {
  return std::transform_reduce(
      models.begin(),
      models.end(),
      Complex{},
      std::plus<>{},
      [f = frequency](auto& mod) { return mod.at(f); });
}

void full::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}

void full::set_model(std::shared_ptr<ArrayOfXsecRecord> xsecrec_) {
  xsecrec = std::move(xsecrec_);
  adapt();
}
}  // namespace fwd::hxsec
