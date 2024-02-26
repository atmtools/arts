#include "fwd_cia.h"

#include <arts_constexpr_math.h>
#include <physics_funcs.h>

#include <algorithm>
#include <functional>
#include <numeric>

#include "cia.h"
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

Complex full::single::at(const Numeric frequency) const {
  return scl * ciarecords->Extract(frequency, T, extrapol, ignore_errors);
}

void full::adapt() try {
  models.resize(0);

  if (not ciarecords) {
    return;
  }

  if (ciarecords->empty()) {
    return;
  }

  ARTS_USER_ERROR_IF(not atm, "Must have an atmosphere")

  models.reserve(ciarecords->size());
  for (CIARecord& data : *ciarecords) {
    const Numeric VMR1 = atm->operator[](data.Species(0));
    const Numeric VMR2 = atm->operator[](data.Species(1));

    models.emplace_back(
        atm->pressure, atm->temperature, VMR1, VMR2, &data, extrap, robust);
  }
}
ARTS_METHOD_ERROR_CATCH

full::full(std::shared_ptr<AtmPoint> atm_,
           std::shared_ptr<ArrayOfCIARecord> cia,
           Numeric extrap_,
           Index robust_)
    : atm(std::move(atm_)),
      ciarecords(std::move(cia)),
      extrap(extrap_),
      robust(robust_) {
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

void full::set_extrap(Numeric extrap_) {
  extrap = extrap_;
  adapt();
}

void full::set_robust(Index robust_) {
  robust = robust_;
  adapt();
}

void full::set_model(std::shared_ptr<ArrayOfCIARecord> cia) {
  ciarecords = std::move(cia);
  adapt();
}

void full::set_atm(std::shared_ptr<AtmPoint> atm_) {
  atm = std::move(atm_);
  adapt();
}
}  // namespace fwd::cia
