#include "hitran_species.h"

#include "abs_species_tags.h"
#include "enums.h"

namespace Hitran {
QuantumIdentifier from_mol_iso(Index molnum, Index isonum) {
  switch (molnum) {
    case 1:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("H2O-161");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("H2O-181");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("H2O-171");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("H2O-162");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 5: {
          static const SpeciesTag spec("H2O-182");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 6: {
          static const SpeciesTag spec("H2O-172");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 7: {
          static const SpeciesTag spec("H2O-262");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 2:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CO2-626");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CO2-636");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("CO2-628");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("CO2-627");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 5: {
          static const SpeciesTag spec("CO2-638");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 6: {
          static const SpeciesTag spec("CO2-637");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 7: {
          static const SpeciesTag spec("CO2-828");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 8: {
          static const SpeciesTag spec("CO2-728");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 9: {
          static const SpeciesTag spec("CO2-727");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 10: {
          static const SpeciesTag spec("CO2-838");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 11: {
          static const SpeciesTag spec("CO2-837");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 12: {
          static const SpeciesTag spec("CO2-737");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 3:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("O3-666");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("O3-668");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("O3-686");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("O3-667");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 5: {
          static const SpeciesTag spec("O3-676");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 4:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("N2O-446");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("N2O-456");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("N2O-546");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("N2O-448");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 5: {
          static const SpeciesTag spec("N2O-447");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 5:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CO-26");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CO-36");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("CO-28");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("CO-27");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 5: {
          static const SpeciesTag spec("CO-38");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 6: {
          static const SpeciesTag spec("CO-37");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 6:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CH4-211");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CH4-311");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("CH4-212");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("CH4-312");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 7:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("O2-66");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("O2-68");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("O2-67");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 8:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("NO-46");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("NO-56");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("NO-48");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 9:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("SO2-626");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("SO2-646");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 10:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("NO2-646");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("NO2-656");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 11:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("NH3-4111");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("NH3-5111");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 12:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HNO3-146");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HNO3-156");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 13:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("OH-61");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("OH-81");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("OH-62");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 14:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HF-19");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HF-29");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 15:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HCl-15");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HCl-17");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("HCl-25");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("HCl-27");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 16:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HBr-19");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HBr-11");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("HBr-29");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("HBr-21");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 17:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HI-17");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HI-27");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 18:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("ClO-56");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("ClO-76");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 19:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("OCS-622");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("OCS-624");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("OCS-632");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("OCS-623");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 5: {
          static const SpeciesTag spec("OCS-822");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 6: {
          static const SpeciesTag spec("OCS-634");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 20:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("H2CO-1126");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("H2CO-1136");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("H2CO-1128");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 21:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HOCl-165");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HOCl-167");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 22:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("N2-44");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("N2-45");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 23:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HCN-124");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HCN-134");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("HCN-125");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 24:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CH3Cl-215");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CH3Cl-217");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 25:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("H2O2-1661");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 26:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("C2H2-1221");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("C2H2-1231");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("C2H2-1222");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 27:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("C2H6-1221");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("C2H6-1231");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 28:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("PH3-1111");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 29:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("COF2-269");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("COF2-369");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 30:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("SF6-29");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 31:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("H2S-121");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("H2S-141");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("H2S-131");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 32:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HCOOH-1261");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 33:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HO2-166");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 34:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("O-6");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 35:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("ClONO2-5646");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("ClONO2-7646");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 36:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("NO+-46");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 37:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HOBr-169");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("HOBr-161");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 38:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("C2H4-221");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("C2H4-231");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 39:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CH3OH-2161");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 40:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CH3Br-219");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CH3Br-211");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 41:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CH3CN-211124");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 42:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CF4-29");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 43:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("C4H2-2211");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 44:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("HC3N-12224");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 45:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("H2-11");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("H2-12");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 46:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CS-22");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CS-24");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("CS-32");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("CS-23");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 47:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("SO3-26");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 48:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("C2N2-4224");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 49:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("COCl2-2655");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("COCl2-2657");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      ARTS_USER_ERROR("Cannot understand isotopologue for: ", isonum, " for species with first entry: ", from_mol_iso(molnum, 1).SpeciesName())
      break;
    case 53:
      switch (isonum) {
        case 1: {
          static const SpeciesTag spec("CS2-222");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 2: {
          static const SpeciesTag spec("CS2-222");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 3: {
          static const SpeciesTag spec("CS2-223");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
        case 4: {
          static const SpeciesTag spec("CS2-232");
          return QuantumIdentifier(QuantumIdentifier::TRANSITION,
                                   spec.Species(), spec.Isotopologue());
        } break;
      }
      break;
  }

  ARTS_USER_ERROR("[HITRAN SPECIES] Does not recognize species ", molnum, " iso ", isonum);
}

QuantumIdentifier from_lookup(Index mol, char isochar) {
  Index iso;
  if (isochar == '1') iso = 1;
  if (isochar == '2') iso = 2;
  if (isochar == '3') iso = 3;
  if (isochar == '4') iso = 4;
  if (isochar == '5') iso = 5;
  if (isochar == '6') iso = 6;
  if (isochar == '7') iso = 7;
  if (isochar == '8') iso = 8;
  if (isochar == '9') iso = 9;
  if (isochar == '0') iso = 10;
  if (isochar == 'A') iso = 11;
  if (isochar == 'B') iso = 12;
  return from_mol_iso(mol, iso);
}

}  // namespace Hitran
