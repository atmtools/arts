#include "fwd_hxsec.h"

#include "physics_funcs.h"

namespace fwd::hxsec {
full::single::single(Numeric p, Numeric t, Numeric VMR, XsecRecord* xsec)
    : scl{number_density(p, t) * VMR}, P{p}, T{t}, xsecrec{xsec} {}

Complex full::single::at(Numeric f) const {
  Numeric out{};
  xsecrec->Extract(ExhaustiveVectorView{out}, Vector{f}, P, T);
  return scl * out;
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

void full::adapt() {
  models.resize(0);
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

void full::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}

void full::set_model(std::shared_ptr<ArrayOfXsecRecord> xsecrec_) {
  xsecrec = std::move(xsecrec_);
  adapt();
}
}  // namespace fwd::hxsec
