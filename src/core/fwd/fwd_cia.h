#pragma once

#include <atm.h>
#include <cia.h>

#include <memory>
#include <vector>

#include "matpack_concepts.h"
#include "species_tags.h"

namespace fwd::cia {

struct full {
  struct single {
    Numeric scl{};
    Numeric T;
    Numeric extrapol;
    Index ignore_errors;
    CIARecord* ciarecords;

    single(Numeric p,
           Numeric t,
           Numeric VMR1,
           Numeric VMR2,
           CIARecord* cia,
           Numeric extrap = {},
           Index robust = {});

    [[nodiscard]] Complex at(Numeric f) const;
    void at(ExhaustiveComplexVectorView abs, const Vector& fs) const;
    [[nodiscard]] ComplexVector at(const Vector& fs) const;
  };

  std::shared_ptr<ArrayOfCIARecord> ciarecords;
  std::vector<single> models{};

  full() = default;

  full(const AtmPoint& atm_point,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::shared_ptr<ArrayOfCIARecord>& cia,
       Numeric extrap = {},
       Index robust = {});

  [[nodiscard]] Complex operator()(Numeric f) const;
  void operator()(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector operator()(const Vector& fs) const;
};
}  // namespace fwd::cia
