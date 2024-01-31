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
class propmat_operator {
  std::shared_ptr<AtmPoint> atm{};
  lbl::zeeman::pol pol{lbl::zeeman::pol::no};

  lbl::fwd::line_storage lines{};
  cia::full cia{};
  predef::full predef{};
  hxsec::full xsec{};

 public:
  propmat_operator() = default;
  propmat_operator(const propmat_operator&) = default;
  propmat_operator(propmat_operator&&) = default;
  propmat_operator& operator=(const propmat_operator&) = default;
  propmat_operator& operator=(propmat_operator&&) = default;

  propmat_operator(std::shared_ptr<AtmPoint> atm,
                   std::shared_ptr<AbsorptionBands> lines,
                   std::shared_ptr<ArrayOfCIARecord> cia,
                   std::shared_ptr<ArrayOfXsecRecord> xsec,
                   std::shared_ptr<PredefinedModelData> predef,
                   Numeric ciaextrap = {},
                   Index ciarobust = {},
                   lbl::zeeman::pol pol = lbl::zeeman::pol::no);

  std::pair<Propmat, Stokvec> operator()(const Numeric frequency,
                                         const Vector2 los) const;

  void set_atm(std::shared_ptr<AtmPoint> atm);
  void set_pol(lbl::zeeman::pol pol);
  void set_ciaextrap(Numeric extrap);
  void set_ciarobust(Index robust);
  void set_bands(std::shared_ptr<AbsorptionBands> lines);
  void set_cia(std::shared_ptr<ArrayOfCIARecord> cia);
  void set_predef(std::shared_ptr<PredefinedModelData> predef);
  void set_model(std::shared_ptr<ArrayOfXsecRecord> xsec);

};  // struct propmat_operator
}  // namespace fwd