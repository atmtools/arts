#ifndef partfun_h
#define partfun_h

#include "auto_partfun.h"

#include "enums.h"
#include "matpackI.h"
#include "isotopologues.h"

namespace PartitionFunctions {

namespace detail {
template <Derivatives d>
constexpr Numeric partfun_impl(Numeric T, const Species::IsotopeRecord& ir) {
  using Species::Species;
  
#define deal_with_spec(SPEC) case Species::SPEC: return compute##SPEC<d>(T, ir.isotname);
  
  switch (ir.spec) {
    case Species::Bath: break;
    deal_with_spec(Water)
    deal_with_spec(CarbonDioxide)
    deal_with_spec(Ozone)
    deal_with_spec(NitrogenOxide)
    deal_with_spec(CarbonMonoxide)
    deal_with_spec(Methane)
    deal_with_spec(Oxygen)
    deal_with_spec(NitricOxide)
    deal_with_spec(SulfurDioxide)
    deal_with_spec(NitrogenDioxide)
    deal_with_spec(Ammonia)
    deal_with_spec(NitricAcid)
    deal_with_spec(Hydroxyl)
    deal_with_spec(HydrogenFluoride)
    deal_with_spec(HydrogenChloride)
    deal_with_spec(HydrogenBromide)
    deal_with_spec(HydrogenIodide)
    deal_with_spec(ChlorineMonoxide)
    deal_with_spec(CarbonylSulfide)
    deal_with_spec(Formaldehyde)
    deal_with_spec(HeavyFormaldehyde)
    deal_with_spec(VeryHeavyFormaldehyde)
    deal_with_spec(HypochlorousAcid)
    deal_with_spec(Nitrogen)
    deal_with_spec(HydrogenCyanide)
    deal_with_spec(MethylChloride)
    deal_with_spec(HydrogenPeroxide)
    deal_with_spec(Acetylene)
    deal_with_spec(Ethane)
    deal_with_spec(Phosphine)
    deal_with_spec(CarbonylFluoride)
    deal_with_spec(SulfurHexafluoride)
    deal_with_spec(HydrogenSulfide)
    deal_with_spec(FormicAcid)
    deal_with_spec(LeftHeavyFormicAcid)
    deal_with_spec(RightHeavyFormicAcid)
    deal_with_spec(Hydroperoxyl)
    deal_with_spec(OxygenAtom)
    deal_with_spec(ChlorineNitrate)
    deal_with_spec(NitricOxideCation)
    deal_with_spec(HypobromousAcid)
    deal_with_spec(Ethylene)
    deal_with_spec(Methanol)
    deal_with_spec(MethylBromide)
    deal_with_spec(Acetonitrile)
    deal_with_spec(HeavyAcetonitrile)
    deal_with_spec(CarbonTetrafluoride)
    deal_with_spec(Diacetylene)
    deal_with_spec(Cyanoacetylene)
    deal_with_spec(Hydrogen)
    deal_with_spec(CarbonMonosulfide)
    deal_with_spec(SulfurTrioxide)
    deal_with_spec(Cyanogen)
    deal_with_spec(Phosgene)
    deal_with_spec(SulfurMonoxide)
    deal_with_spec(CarbonDisulfide)
    deal_with_spec(Methyl)
    deal_with_spec(Cyclopropene)
    deal_with_spec(SulfuricAcid)
    deal_with_spec(HydrogenIsocyanide)
    deal_with_spec(BromineMonoxide)
    deal_with_spec(ChlorineDioxide)
    deal_with_spec(Propane)
    deal_with_spec(Helium)
    deal_with_spec(ChlorineMonoxideDimer)
    deal_with_spec(HydrogenAtom)
    deal_with_spec(Argon)
    deal_with_spec(Hexafluoroethane)
    deal_with_spec(Perfluoropropane)
    deal_with_spec(Perfluorobutane)
    deal_with_spec(Perfluoropentane)
    deal_with_spec(Perfluorohexane)
    deal_with_spec(Perfluorooctane)
    deal_with_spec(Perfluorocyclobutane)
    deal_with_spec(CarbonTetrachloride)
    deal_with_spec(CFC11)
    deal_with_spec(CFC113)
    deal_with_spec(CFC114)
    deal_with_spec(CFC115)
    deal_with_spec(CFC12)
    deal_with_spec(Dichloromethane)
    deal_with_spec(Trichloroethane)
    deal_with_spec(Trichloromethane)
    deal_with_spec(Bromochlorodifluoromethane)
    deal_with_spec(Bromotrifluoromethane)
    deal_with_spec(Dibromotetrafluoroethane)
    deal_with_spec(HCFC141b)
    deal_with_spec(HCFC142b)
    deal_with_spec(HCFC22)
    deal_with_spec(HFC125)
    deal_with_spec(HFC134a)
    deal_with_spec(HFC143a)
    deal_with_spec(HFC152a)
    deal_with_spec(HFC227ea)
    deal_with_spec(HFC23)
    deal_with_spec(HFC245fa)
    deal_with_spec(HFC32)
    deal_with_spec(NitrogenTrifluoride)
    deal_with_spec(SulfurylFluoride)
    deal_with_spec(HFC4310mee)
    deal_with_spec(liquidcloud)
    deal_with_spec(icecloud)
    deal_with_spec(rain)
    deal_with_spec(free_electrons)
    deal_with_spec(particles)
    case Species::FINAL: { /* leave last */
    }
  }
  
#undef deal_with_spec
  
  ARTS_USER_ERROR("This is not a valid IsotopeRecord:\n", ir)
}
}

constexpr Numeric Q(Numeric T, const Species::IsotopeRecord& ir) {
  return detail::partfun_impl<Derivatives::No>(T, ir);
}

constexpr Numeric dQdT(Numeric T, const Species::IsotopeRecord& ir) {
  return detail::partfun_impl<Derivatives::Yes>(T, ir);
}

constexpr bool has_partfun(const Species::IsotopeRecord& ir) noexcept {
  using Species::Species;
  
  #define deal_with_spec(SPEC) case Species::SPEC: for (auto& x: has##SPEC) if (x == ir.isotname) return true; break;
  
  switch (ir.spec) {
    case Species::Bath: break;
    deal_with_spec(Water)
    deal_with_spec(CarbonDioxide)
    deal_with_spec(Ozone)
    deal_with_spec(NitrogenOxide)
    deal_with_spec(CarbonMonoxide)
    deal_with_spec(Methane)
    deal_with_spec(Oxygen)
    deal_with_spec(NitricOxide)
    deal_with_spec(SulfurDioxide)
    deal_with_spec(NitrogenDioxide)
    deal_with_spec(Ammonia)
    deal_with_spec(NitricAcid)
    deal_with_spec(Hydroxyl)
    deal_with_spec(HydrogenFluoride)
    deal_with_spec(HydrogenChloride)
    deal_with_spec(HydrogenBromide)
    deal_with_spec(HydrogenIodide)
    deal_with_spec(ChlorineMonoxide)
    deal_with_spec(CarbonylSulfide)
    deal_with_spec(Formaldehyde)
    deal_with_spec(HeavyFormaldehyde)
    deal_with_spec(VeryHeavyFormaldehyde)
    deal_with_spec(HypochlorousAcid)
    deal_with_spec(Nitrogen)
    deal_with_spec(HydrogenCyanide)
    deal_with_spec(MethylChloride)
    deal_with_spec(HydrogenPeroxide)
    deal_with_spec(Acetylene)
    deal_with_spec(Ethane)
    deal_with_spec(Phosphine)
    deal_with_spec(CarbonylFluoride)
    deal_with_spec(SulfurHexafluoride)
    deal_with_spec(HydrogenSulfide)
    deal_with_spec(FormicAcid)
    deal_with_spec(LeftHeavyFormicAcid)
    deal_with_spec(RightHeavyFormicAcid)
    deal_with_spec(Hydroperoxyl)
    deal_with_spec(OxygenAtom)
    deal_with_spec(ChlorineNitrate)
    deal_with_spec(NitricOxideCation)
    deal_with_spec(HypobromousAcid)
    deal_with_spec(Ethylene)
    deal_with_spec(Methanol)
    deal_with_spec(MethylBromide)
    deal_with_spec(Acetonitrile)
    deal_with_spec(HeavyAcetonitrile)
    deal_with_spec(CarbonTetrafluoride)
    deal_with_spec(Diacetylene)
    deal_with_spec(Cyanoacetylene)
    deal_with_spec(Hydrogen)
    deal_with_spec(CarbonMonosulfide)
    deal_with_spec(SulfurTrioxide)
    deal_with_spec(Cyanogen)
    deal_with_spec(Phosgene)
    deal_with_spec(SulfurMonoxide)
    deal_with_spec(CarbonDisulfide)
    deal_with_spec(Methyl)
    deal_with_spec(Cyclopropene)
    deal_with_spec(SulfuricAcid)
    deal_with_spec(HydrogenIsocyanide)
    deal_with_spec(BromineMonoxide)
    deal_with_spec(ChlorineDioxide)
    deal_with_spec(Propane)
    deal_with_spec(Helium)
    deal_with_spec(ChlorineMonoxideDimer)
    deal_with_spec(HydrogenAtom)
    deal_with_spec(Argon)
    deal_with_spec(Hexafluoroethane)
    deal_with_spec(Perfluoropropane)
    deal_with_spec(Perfluorobutane)
    deal_with_spec(Perfluoropentane)
    deal_with_spec(Perfluorohexane)
    deal_with_spec(Perfluorooctane)
    deal_with_spec(Perfluorocyclobutane)
    deal_with_spec(CarbonTetrachloride)
    deal_with_spec(CFC11)
    deal_with_spec(CFC113)
    deal_with_spec(CFC114)
    deal_with_spec(CFC115)
    deal_with_spec(CFC12)
    deal_with_spec(Dichloromethane)
    deal_with_spec(Trichloroethane)
    deal_with_spec(Trichloromethane)
    deal_with_spec(Bromochlorodifluoromethane)
    deal_with_spec(Bromotrifluoromethane)
    deal_with_spec(Dibromotetrafluoroethane)
    deal_with_spec(HCFC141b)
    deal_with_spec(HCFC142b)
    deal_with_spec(HCFC22)
    deal_with_spec(HFC125)
    deal_with_spec(HFC134a)
    deal_with_spec(HFC143a)
    deal_with_spec(HFC152a)
    deal_with_spec(HFC227ea)
    deal_with_spec(HFC23)
    deal_with_spec(HFC245fa)
    deal_with_spec(HFC32)
    deal_with_spec(NitrogenTrifluoride)
    deal_with_spec(SulfurylFluoride)
    deal_with_spec(HFC4310mee)
    deal_with_spec(liquidcloud)
    deal_with_spec(icecloud)
    deal_with_spec(rain)
    deal_with_spec(free_electrons)
    deal_with_spec(particles)
    case Species::FINAL: { /* leave last */
    }
  }
  
  #undef deal_with_spec
  
  return false;
}
}

#endif  // partfun_h
