#pragma once

#include <lbl.h>

#include <memory>

#include "atm.h"
#include "fwd_cia.h"
#include "fwd_hxsec.h"
#include "fwd_predef.h"
#include "lbl_data.h"
#include "rtepack.h"

namespace fwd {
class propmat {
  std::shared_ptr<AtmPoint> atm{};

  lbl::fwd::line_storage lines{};

  cia::full cia{};
  predef::full predef{};
  hxsec::full xsec{};

 public:
  propmat() = default;
  propmat(const propmat&) = default;
  propmat(propmat&&) = default;
  propmat& operator=(const propmat&) = default;
  propmat& operator=(propmat&&) = default;

  propmat(std::shared_ptr<AtmPoint> atm,
          std::shared_ptr<ArrayOfAbsorptionBand> lines,
          std::shared_ptr<ArrayOfCIARecord> cia,
          std::shared_ptr<ArrayOfXsecRecord> xsec,
          std::shared_ptr<PredefinedModelData> predef,
          Numeric ciaextrap = {},
          Index ciarobust = {});

  std::pair<Propmat, Stokvec> operator()(const Numeric frequency,
                                         const Vector2 los) const;

  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_ciaextrap(Numeric extrap);
  void set_ciarobust(Index robust);
  void set_bands(std::shared_ptr<ArrayOfAbsorptionBand> lines);
  void set_cia(std::shared_ptr<ArrayOfCIARecord> cia);
  void set_predef(std::shared_ptr<PredefinedModelData> predef);
  void set_model(std::shared_ptr<ArrayOfXsecRecord> xsec);
};  // struct propmat
}  // namespace fwd
