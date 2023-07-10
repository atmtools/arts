#pragma once

#include <memory>
#include <vector>
#include "../cia.h"
#include "matpack_concepts.h"
#include "messages.h"
#include "species_tags.h"


namespace cia {
struct single {
  Numeric scl{};
  Numeric T;
  Numeric extrapol;
  Index ignore_errors;
  Verbosity verbosity;
  std::shared_ptr<CIARecord> ciarecords{};

  single() = default;

  single(Numeric p,
         Numeric t,
         Numeric VMR1,
         Numeric VMR2,
         const std::shared_ptr<CIARecord>& cia,
         Numeric extrap={},
         Index robust={},
         Verbosity verb={});

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ComplexVector& abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};

struct full {
  std::vector<single> models{};

  full() = default;

  full(Numeric p,
       Numeric t,
       const Vector& vmrs,
       const ArrayOfArrayOfSpeciesTag& allspecs,
       const std::vector<std::shared_ptr<CIARecord>>& cia,
       Numeric extrap={},
       Index robust={},
       Verbosity verb={});

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ComplexVector& abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};
}  // namespace cia