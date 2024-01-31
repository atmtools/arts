#pragma once

#include <atm.h>
#include <cia.h>

#include <memory>
#include <vector>

namespace fwd::cia {
class full {
  struct single {
    Numeric scl{};
    Numeric T;
    Numeric extrapol;
    Index ignore_errors;
    CIARecord* ciarecords;

    single() = default;

    single(Numeric p,
           Numeric t,
           Numeric VMR1,
           Numeric VMR2,
           CIARecord* cia,
           Numeric extrap,
           Index robust);

    [[nodiscard]] Complex at(Numeric f) const;
    void at(ExhaustiveComplexVectorView abs, const Vector& fs) const;
    [[nodiscard]] ComplexVector at(const Vector& fs) const;
  };

  std::shared_ptr<AtmPoint> atm;
  std::shared_ptr<ArrayOfCIARecord> ciarecords;
  Numeric extrap;
  Index robust;

  std::vector<single> models{};

  void adapt();

 public:
  full() = default;

  full(std::shared_ptr<AtmPoint> atm,
       std::shared_ptr<ArrayOfCIARecord> cia,
       Numeric extrap = {},
       Index robust = {});

  [[nodiscard]] Complex operator()(Numeric f) const;
  void operator()(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector operator()(const Vector& fs) const;

  void set_extrap(Numeric extrap);
  void set_robust(Index robust);
  void set_lines(std::shared_ptr<ArrayOfCIARecord> cia);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};
}  // namespace fwd::cia
