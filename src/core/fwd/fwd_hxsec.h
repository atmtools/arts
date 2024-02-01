#pragma once

#include <atm.h>
#include <species_tags.h>
#include <xsec_fit.h>

#include <memory>

namespace fwd::hxsec {

class full {
  struct single {
    Numeric scl{};
    Numeric P{};
    Numeric T{};
    XsecRecord* xsecrec{};

    single() = default;

    single(Numeric p, Numeric t, Numeric VMR, XsecRecord* xsec);

    [[nodiscard]] Complex at(const Numeric frequency) const;
  };

  std::shared_ptr<AtmPoint> atm{};
  std::shared_ptr<ArrayOfXsecRecord> xsecrec{};
  std::vector<single> models{};

  void adapt();

public:
  full() = default;

  full(std::shared_ptr<AtmPoint> atm,
       std::shared_ptr<ArrayOfXsecRecord> xsecrec);

  [[nodiscard]] Complex operator()(const Numeric frequency) const;

  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_model(std::shared_ptr<ArrayOfXsecRecord> xsecrec);
};
}  // namespace fwd::hxsec
