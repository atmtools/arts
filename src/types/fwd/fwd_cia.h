#pragma once

#include <memory>
#include <vector>

#include <atm.h>
#include <cia.h>
#include "matpack_concepts.h"
#include "species_tags.h"

namespace fwd::cia {
struct single {
  Numeric scl{};
  Numeric T;
  Numeric extrapol;
  Index ignore_errors;
  std::shared_ptr<CIARecord> ciarecords{};

  single() = default;

  single(Numeric p,
         Numeric t,
         Numeric VMR1,
         Numeric VMR2,
         const std::shared_ptr<CIARecord>& cia,
         Numeric extrap = {},
         Index robust = {});

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};

struct full {
  std::vector<single> models{};

  full() = default;

  full(const AtmPoint& atm_point,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::vector<std::shared_ptr<CIARecord>>& cia,
       Numeric extrap = {},
       Index robust = {});

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};
}  // namespace fwd::cia
