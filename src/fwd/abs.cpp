#include "abs.h"

namespace fwd {
full::full(Numeric p,
           Numeric t,
           const Vector& allvmrs,
           const ArrayOfArrayOfSpeciesTag& allspecs,
           const std::shared_ptr<PredefinedModelData>& predef_data,
           const std::vector<std::shared_ptr<CIARecord>>& cia_data,
           const SpeciesIsotopologueRatios& isotopologue_ratios,
           const ArrayOfArrayOfAbsorptionLines& lbl_data,
           Numeric cia_extrap,
           Index cia_robust,
           Verbosity cia_verb)
    : cia(p, t, allvmrs, allspecs, cia_data, cia_extrap, cia_robust, cia_verb),
      predef(p, t, allvmrs, allspecs, predef_data),
      lbl(p, t, isotopologue_ratios, allspecs, allvmrs, lbl_data) {}

Complex full::at(Numeric f) const {
  return cia.at(f) + predef.at(f) + lbl.at(f);
}

void full::at(ExhaustiveComplexVectorView abs, const Vector& fs) const {
  cia.at(abs, fs);
  predef.at(abs, fs);
  lbl.at(abs, fs);
}

ComplexVector full::at(const Vector& fs) const {
  ComplexVector abs(fs.size());
  at(abs, fs);
  return abs;
}
}  // namespace fwd
