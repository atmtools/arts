#include "fwd_abs.h"

namespace fwd {
full_absorption::full_absorption(
    const AtmPoint& atm_point,
    const ArrayOfArrayOfSpeciesTag& allspecs,
    const std::shared_ptr<PredefinedModelData>& predef_data,
    const std::vector<std::shared_ptr<CIARecord>>& cia_data,
    const std::vector<std::shared_ptr<XsecRecord>>& hxsec_data,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ArrayOfArrayOfAbsorptionLines& lbl_data,
    Numeric cia_extrap,
    Index cia_robust)
    : cia(atm_point, allspecs, cia_data, cia_extrap, cia_robust),
      predef(atm_point, allspecs, predef_data),
      lbl(atm_point, isotopologue_ratios, lbl_data),
      hxsec(atm_point, allspecs, hxsec_data) {}

Complex full_absorption::at(Numeric f) const {
  return cia.at(f) + predef.at(f) + lbl.at(f) + hxsec.at(f);
}

void full_absorption::at(ExhaustiveComplexVectorView abs,
                         const Vector& fs) const {
  cia.at(abs, fs);
  predef.at(abs, fs);
  lbl.at(abs, fs);
  hxsec.at(abs, fs);
}

ComplexVector full_absorption::at(const Vector& fs) const {
  ComplexVector abs(fs.size());
  at(abs, fs);
  return abs;
}
}  // namespace fwd
