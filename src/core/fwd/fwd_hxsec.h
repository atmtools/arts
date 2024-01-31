#pragma once

#include <atm.h>
#include <species_tags.h>
#include <xsec_fit.h>

#include <memory>

namespace fwd::hxsec {

struct full {
  struct single {
    Numeric scl{};
    Numeric P{};
    Numeric T{};
    XsecRecord* xsecrec{};

    single() = default;

    single(Numeric p, Numeric t, Numeric VMR, XsecRecord* xsec);

    [[nodiscard]] Complex at(Numeric f) const;
    void at(ExhaustiveComplexVectorView abs, const Vector& fs) const;
    [[nodiscard]] ComplexVector at(const Vector& fs) const;
  };

  std::shared_ptr<ArrayOfXsecRecord> xsecrec;
  std::vector<single> models{};

  full() = default;

  full(const AtmPoint& atm_point,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::shared_ptr<ArrayOfXsecRecord>& xsec);

  [[nodiscard]] Complex operator()(Numeric f) const;
  void operator()(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector operator()(const Vector& fs) const;
};
}  // namespace fwd::hxsec
