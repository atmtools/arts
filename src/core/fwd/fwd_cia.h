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
    single(const single&) = default;
    single(single&&) = default;
    single& operator=(const single&) = default;
    single& operator=(single&&) = default;

    single(Numeric p,
           Numeric t,
           Numeric VMR1,
           Numeric VMR2,
           CIARecord* cia,
           Numeric extrap,
           Index robust);

    [[nodiscard]] Complex at(const Numeric frequency) const;
  };

  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<ArrayOfCIARecord> ciarecords{};
  Numeric extrap{};
  Index robust{};

  std::vector<single> models{};

  void adapt();

 public:
  full() = default;
  full(const full&) = default;
  full(full&&) = default;
  full& operator=(const full&) = default;
  full& operator=(full&&) = default;

  full(std::shared_ptr<AtmPoint> atm,
       std::shared_ptr<ArrayOfCIARecord> cia,
       Numeric extrap = {},
       Index robust = {});

  [[nodiscard]] Complex operator()(const Numeric frequency) const;

  void set_extrap(Numeric extrap);
  void set_robust(Index robust);
  void set_model(std::shared_ptr<ArrayOfCIARecord> cia);
  void set_atm(std::shared_ptr<AtmPoint> atm);
};
}  // namespace fwd::cia
