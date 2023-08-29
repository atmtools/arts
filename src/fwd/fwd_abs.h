#pragma once

#include "fwd_cia.h"
#include "fwd_hxsec.h"
#include "fwd_lbl.h"
#include "fwd_predef.h"

namespace fwd {
struct full_absorption {
  cia::full cia;
  predef::full predef;
  lbl::full lbl;
  hxsec::full hxsec;

  full_absorption() = default;

  full_absorption(const AtmPoint& atm_point,
                  const ArrayOfArrayOfSpeciesTag& allspecs,
                  const std::shared_ptr<PredefinedModelData>& predef_data,
                  const std::vector<std::shared_ptr<CIARecord>>& cia_data,
                  const std::vector<std::shared_ptr<XsecRecord>>& hxsec_data,
                  const SpeciesIsotopologueRatios& isotopologue_ratios,
                  const ArrayOfArrayOfAbsorptionLines& lbl_data,
                  Numeric cia_extrap = {},
                  Index cia_robust = {});

  [[nodiscard]] Complex at(Numeric f) const;
  void at(ExhaustiveComplexVectorView abs, const Vector& fs) const;
  [[nodiscard]] ComplexVector at(const Vector& fs) const;
};
}  // namespace fwd
