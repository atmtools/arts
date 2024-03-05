/** Contains the line shape namespace
 * @file   lineshapemodel.cc
 * @author Richard Larsson
 * @date   2018-09-19
 * 
 * @brief Implementations of lineshapemodel.h
 * 
 * This mostly contains functions that either did not
 * compile while placed in the header or for other 
 * reasons are not there.  This means most of the 
 * real work happens in the header file and not here
 **/

#include "lineshapemodel.h"

#include <limits>

#include "atm.h"
#include "debug.h"
#include "double_imanip.h"
#include "enums.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "species.h"

std::istream& LineShape::from_artscat4(std::istream& is,
                                       LineShapeTypeOld& mtype,
                                       bool& self,
                                       bool& bath,
                                       Model& m,
                                       ArrayOfSpeciesEnum& species,
                                       const QuantumIdentifier& qid) {
  // Set or reset variables
  mtype = LineShapeTypeOld::VP;
  self = true;
  bath = false;
  m.mdata = std::vector<SingleSpeciesModel>(7);
  species = ArrayOfSpeciesEnum(7);

  // Set species
  species[1] = to<SpeciesEnum>("N2");
  species[2] = to<SpeciesEnum>("O2");
  species[3] = to<SpeciesEnum>("H2O");
  species[4] = to<SpeciesEnum>("CO2");
  species[5] = to<SpeciesEnum>("H2");
  species[6] = to<SpeciesEnum>("He");

  // Temperature types
  for (auto& v : m.mdata) {
    v.G0().type = LineShapeTemperatureModelOld::T1;
    v.D0().type = LineShapeTemperatureModelOld::T5;
  }

  // G0 main coefficient
  for (auto& v : m.mdata) is >> double_imanip() >> v.G0().X0;

  // G0 exponent is same as D0 exponent
  for (auto& v : m.mdata) {
    is >> double_imanip() >> v.G0().X1;
    v.D0().X1 = v.G0().X1;
  }

  // D0 coefficient
  m.mdata.front().D0().X0 = 0;
  for (int k = 1; k < 7; k++) is >> double_imanip() >> m.mdata[k].D0().X0;

  // Special case when self is part of this list, it needs to be removed
  for (int k = 1; k < 7; k++) {
    if (qid.Species() == species[k]) {
      ARTS_USER_ERROR_IF(
          m.mdata.front().G0().X0 not_eq m.mdata[k].G0().X0 or
              m.mdata.front().G0().X1 not_eq m.mdata[k].G0().X1 or
              m.mdata.front().D0().X1 not_eq m.mdata[k].D0().X1,
          "Species is ",
          qid.Isotopologue(),
          " and this is a broadening species in ARTSCAT-4.\n"
          "Despite this, values representing self and ",
          qid.Isotopologue(),
          " does not match "
          "in input string\n")
      m.mdata.front().D0().X0 = m.mdata[k].D0().X0;
      m.Remove(k, species);
      break;
    }
  }

  return is;
}

std::istream& LineShape::from_linefunctiondata(std::istream& data,
                                               LineShapeTypeOld& mtype,
                                               bool& self,
                                               bool& bath,
                                               Model& m,
                                               ArrayOfSpeciesEnum& species) {
  self = bath = false;
  Index specs;
  String s;

  // The first tag should give the line shape scheme
  data >> mtype;

  // Order of elements for line shape
  const auto shapeparams = LegacyLineFunctionData::lineshapetag2variablesvector(
      String{toString(mtype)});

  // The second tag should give the line mixing scheme
  data >> s;

  // Order of elements for line mixing
  const auto mixingparams =
      LegacyLineFunctionData::linemixingtag2variablesvector(s);

  // The third tag should contain the number of species
  data >> specs;
  species.resize(specs);
  m.mdata.resize(specs);

  ARTS_USER_ERROR_IF(not specs and mtype not_eq LineShapeTypeOld::DP,
                     "Need at least one species for non-Doppler line shapes");

  // For all species, we need to set the methods to compute them
  for (Index i = 0; i < specs; i++) {
    // This should be a species tag or one of the specials, SELF or BATH
    data >> s;
    if (s == self_broadening) {
      // If the species is self, then  we need to flag this
      self = true;
      ARTS_USER_ERROR_IF(i not_eq 0,
                         "Self broadening must be first, it is not\n");
    } else if (s == bath_broadening) {
      // If the species is air, then we need to flag this
      bath = true;
      species[i] = SpeciesEnum::Bath;
      ARTS_USER_ERROR_IF(i not_eq specs - 1,
                         "Air/bath broadening must be last, it is not\n");
    } else {
      // Otherwise, we hope we find a species
      species[i] = to<SpeciesEnum>(s);
      ARTS_USER_ERROR_IF(
          not good_enum(species[i]),
          "Encountered ",
          s,
          " in a position where a species should have been "
          "defined.\nPlease check your pressure broadening data structure and ensure "
          "that it follows the correct conventions.\n")
    }

    // For all parameters
    for (auto& params : {shapeparams, mixingparams}) {
      for (auto& param : params) {
        data >> s;  // Should contain a temperature tag

        const auto type = to<LineShapeTemperatureModelOld>(s);
        const Index ntemp =
            LegacyLineFunctionData::temperaturemodel2legacysize(type);

        m.mdata[i].Data()[Index(param)].type = type;
        if (ntemp <= ModelParameters::N) {
          switch (ntemp) {
            case 1:
              data >> double_imanip() >> m.mdata[i].Data()[Index(param)].X0;
              break;
            case 2:
              data >> double_imanip() >> m.mdata[i].Data()[Index(param)].X0 >>
                  m.mdata[i].Data()[Index(param)].X1;
              break;
            case 3:
              data >> double_imanip() >> m.mdata[i].Data()[Index(param)].X0 >>
                  m.mdata[i].Data()[Index(param)].X1 >>
                  m.mdata[i].Data()[Index(param)].X2;
              break;
            case 0:
              break;
            default:
              ARTS_USER_ERROR(
                  "Unknown number of input parameters in Legacy mode.");
          }
        } else {  // Has to be the only allowed interpolation case
          ARTS_USER_ERROR_IF(
              ntemp > 12,
              "Too many input parameters in interpolation results Legacy mode.");
          Numeric temp;
          data >> double_imanip() >> temp;  // should be 200
          data >> double_imanip() >> temp;  // should be 250
          data >> double_imanip() >> temp;  // should be 296
          data >> double_imanip() >> temp;  // should be 340
          data >> double_imanip() >> m.mdata[i].Y().X0 >> m.mdata[i].Y().X1 >>
              m.mdata[i].Y().X2 >> m.mdata[i].Y().X3 >> m.mdata[i].G().X0 >>
              m.mdata[i].G().X1 >> m.mdata[i].G().X2 >> m.mdata[i].G().X3;
        }
      }
    }
  }

  return data;
}

std::istream& LineShape::from_pressurebroadeningdata(
    std::istream& data,
    LineShapeTypeOld& mtype,
    bool& self,
    bool& bath,
    Model& m,
    ArrayOfSpeciesEnum& species,
    const QuantumIdentifier& qid) {
  String s;
  data >> s;

  const auto type = LegacyPressureBroadeningData::string2typepb(s);
  const auto n = LegacyPressureBroadeningData::typepb2size(type);
  const auto self_in_list =
      LegacyPressureBroadeningData::self_listed(qid, type);

  Vector x(n);
  for (auto& num : x) data >> double_imanip() >> num;

  LegacyPressureBroadeningData::vector2modelpb(
      mtype, self, bath, m, species, x, type, self_in_list, qid.Species());

  return data;
}

std::istream& LineShape::from_linemixingdata(std::istream& data,
                                             LineShape::Model& lsc) {
  String s;
  data >> s;

  const auto type = LegacyLineMixingData::string2typelm(s);
  const auto n = LegacyLineMixingData::typelm2size(type);

  Vector x(n);
  for (auto& num : x) data >> double_imanip() >> num;

  lsc = LegacyLineMixingData::vector2modellm(x, type);

  return data;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
void LineShape::LegacyPressureBroadeningData::vector2modelpb(
    LineShapeTypeOld& mtype,
    bool& self,
    bool& bath,
    Model& m,
    ArrayOfSpeciesEnum& species,
    Vector x,
    LineShape::LegacyPressureBroadeningData::TypePB type,
    bool self_in_list,
    SpeciesEnum self_spec) {
  switch (type) {
    case TypePB::PB_NONE:
      mtype = LineShapeTypeOld::DP;
      self = bath = false;
      m = Model();
      species.resize(0);
      return;
    case TypePB::PB_AIR_BROADENING:
      mtype = LineShapeTypeOld::VP;
      self = bath = true;
      m = Model(x[0], x[1], x[2], x[3], x[4]);
      species.resize(2);
      species[0] = self_spec;
      species[1] = SpeciesEnum::Bath;
      return;
    case TypePB::PB_AIR_AND_WATER_BROADENING:
      if (self_in_list) {
        mtype = LineShapeTypeOld::VP;
        self = false;
        bath = true;
        m.Data().resize(2);
        m.Data()[0].G0() = {LineShapeTemperatureModelOld::T1, x[0], x[1], 0, 0};
        m.Data()[0].D0() = {LineShapeTemperatureModelOld::T5, x[2], x[1], 0, 0};
        m.Data()[1].G0() = {LineShapeTemperatureModelOld::T1, x[3], x[4], 0, 0};
        m.Data()[1].D0() = {LineShapeTemperatureModelOld::T5, x[5], x[4], 0, 0};
        species.resize(2);
        species[0] = to<SpeciesEnum>("H2O");
        species[1] = SpeciesEnum::Bath;
        return;
      } else {
        mtype = LineShapeTypeOld::VP;
        self = bath = true;
        m.Data().resize(2);
        m.Data()[0].G0() = {LineShapeTemperatureModelOld::T1, x[0], x[1], 0, 0};
        m.Data()[0].D0() = {LineShapeTemperatureModelOld::T5, x[2], x[1], 0, 0};
        m.Data()[2].G0() = {LineShapeTemperatureModelOld::T1, x[3], x[4], 0, 0};
        m.Data()[2].D0() = {LineShapeTemperatureModelOld::T5, x[5], x[4], 0, 0};
        m.Data()[1].G0() = {LineShapeTemperatureModelOld::T1, x[6], x[7], 0, 0};
        m.Data()[1].D0() = {LineShapeTemperatureModelOld::T5, x[8], x[7], 0, 0};
        species.resize(3);
        species[0] = self_spec;
        species[1] = to<SpeciesEnum>("H2O");
        species[2] = SpeciesEnum::Bath;
        return;
      }
    case TypePB::PB_PLANETARY_BROADENING:
      if (self_in_list) {
        mtype = LineShapeTypeOld::VP;
        self = bath = false;
        m.Data().resize(6);
        m.Data()[0].G0() = {LineShapeTemperatureModelOld::T1, x[1], x[8], 0, 0};
        m.Data()[0].D0() = {LineShapeTemperatureModelOld::T5, x[14], x[8], 0, 0};
        m.Data()[1].G0() = {LineShapeTemperatureModelOld::T1, x[2], x[9], 0, 0};
        m.Data()[1].D0() = {LineShapeTemperatureModelOld::T5, x[15], x[9], 0, 0};
        m.Data()[2].G0() = {LineShapeTemperatureModelOld::T1, x[3], x[10], 0, 0};
        m.Data()[2].D0() = {LineShapeTemperatureModelOld::T5, x[16], x[10], 0, 0};
        m.Data()[3].G0() = {LineShapeTemperatureModelOld::T1, x[4], x[11], 0, 0};
        m.Data()[3].D0() = {LineShapeTemperatureModelOld::T5, x[17], x[11], 0, 0};
        m.Data()[4].G0() = {LineShapeTemperatureModelOld::T1, x[5], x[12], 0, 0};
        m.Data()[4].D0() = {LineShapeTemperatureModelOld::T5, x[18], x[12], 0, 0};
        m.Data()[5].G0() = {LineShapeTemperatureModelOld::T1, x[6], x[13], 0, 0};
        m.Data()[5].D0() = {LineShapeTemperatureModelOld::T5, x[19], x[13], 0, 0};
        species = {to<SpeciesEnum>("N2"),
                   to<SpeciesEnum>("O2"),
                   to<SpeciesEnum>("H2O"),
                   to<SpeciesEnum>("CO2"),
                   to<SpeciesEnum>("H2"),
                   to<SpeciesEnum>("He")};
        return;
      } else {
        mtype = LineShapeTypeOld::VP;
        self = true;
        bath = false;
        m.Data().resize(7);
        m.Data()[0].G0() = {LineShapeTemperatureModelOld::T1, x[0], x[7], 0, 0};
        //          ssm[0].D0() = ...
        m.Data()[1].G0() = {LineShapeTemperatureModelOld::T1, x[1], x[8], 0, 0};
        m.Data()[1].D0() = {LineShapeTemperatureModelOld::T5, x[14], x[8], 0, 0};
        m.Data()[2].G0() = {LineShapeTemperatureModelOld::T1, x[2], x[9], 0, 0};
        m.Data()[2].D0() = {LineShapeTemperatureModelOld::T5, x[15], x[9], 0, 0};
        m.Data()[3].G0() = {LineShapeTemperatureModelOld::T1, x[3], x[10], 0, 0};
        m.Data()[3].D0() = {LineShapeTemperatureModelOld::T5, x[16], x[10], 0, 0};
        m.Data()[4].G0() = {LineShapeTemperatureModelOld::T1, x[4], x[11], 0, 0};
        m.Data()[4].D0() = {LineShapeTemperatureModelOld::T5, x[17], x[11], 0, 0};
        m.Data()[5].G0() = {LineShapeTemperatureModelOld::T1, x[5], x[12], 0, 0};
        m.Data()[5].D0() = {LineShapeTemperatureModelOld::T5, x[18], x[12], 0, 0};
        m.Data()[6].G0() = {LineShapeTemperatureModelOld::T1, x[6], x[13], 0, 0};
        m.Data()[6].D0() = {LineShapeTemperatureModelOld::T5, x[19], x[13], 0, 0};
        species.resize(7);
        species[0] = self_spec;
        species[1] = to<SpeciesEnum>("N2");
        species[2] = to<SpeciesEnum>("O2");
        species[3] = to<SpeciesEnum>("H2O");
        species[4] = to<SpeciesEnum>("CO2");
        species[5] = to<SpeciesEnum>("H2");
        species[6] = to<SpeciesEnum>("He");
        return;
      }
  }
}
#pragma GCC diagnostic pop

LineShape::Model LineShape::LegacyLineMixingData::vector2modellm(
    Vector x, LineShape::LegacyLineMixingData::TypeLM type) {
  Model y(1);
  switch (type) {
    case TypeLM::LM_NONE:
      break;
    case TypeLM::LM_LBLRTM:
      y.Data().front().Y().type = LineShapeTemperatureModelOld::LM_AER;
      y.Data().front().G().type = LineShapeTemperatureModelOld::LM_AER;
      y.Data().front().Y().X0 = x[4];
      y.Data().front().Y().X1 = x[5];
      y.Data().front().Y().X2 = x[6];
      y.Data().front().Y().X3 = x[7];
      y.Data().front().G().X0 = x[8];
      y.Data().front().G().X1 = x[9];
      y.Data().front().G().X2 = x[10];
      y.Data().front().G().X3 = x[11];
      break;
    case TypeLM::LM_LBLRTM_O2NonResonant:
      y.Data().front().G().type = LineShapeTemperatureModelOld::T0;
      y.Data().front().G().X0 = x[0];
      break;
    case TypeLM::LM_2NDORDER:
      y.Data().front().Y().type = LineShapeTemperatureModelOld::T4;
      y.Data().front().Y().X0 = x[0];
      y.Data().front().Y().X1 = x[1];
      y.Data().front().Y().X2 = x[7];
      y.Data().front().G().type = LineShapeTemperatureModelOld::T4;
      y.Data().front().G().X0 = x[2];
      y.Data().front().G().X1 = x[3];
      y.Data().front().G().X2 = x[8];
      y.Data().front().DV().type = LineShapeTemperatureModelOld::T4;
      y.Data().front().DV().X0 = x[4];
      y.Data().front().DV().X1 = x[5];
      y.Data().front().DV().X2 = x[9];
      break;
    case TypeLM::LM_1STORDER:
      y.Data().front().Y().type = LineShapeTemperatureModelOld::T1;
      y.Data().front().Y().X0 = x[1];
      y.Data().front().Y().X1 = x[2];
      break;
    case TypeLM::LM_BYBAND:
      break;
  }
  return y;
}

Vector LineShape::vmrs(const ConstVectorView& atmospheric_vmrs,
                       const ArrayOfArrayOfSpeciesTag& atmospheric_species,
                       const ArrayOfSpeciesEnum& lineshape_species) ARTS_NOEXCEPT {
  ARTS_ASSERT(
      atmospheric_species.size() == static_cast<Size>(atmospheric_vmrs.size()),
      "Bad atmospheric inputs");

  const Size n = lineshape_species.size();

  // Initialize list of VMRS to 0
  Vector line_vmrs(n, 0);

  // We need to know if bath is an actual species
  const bool bath = n and lineshape_species.back() == SpeciesEnum::Bath;

  // Loop species
  for (Size i = 0; i < n - bath; i++) {
    const SpeciesEnum target = lineshape_species[i];

    // Find species in list or do nothing at all
    Index this_species_index = -1;
    for (Size j = 0; j < atmospheric_species.size(); j++) {
      if (atmospheric_species[j].Species() == target) {
        this_species_index = j;
      }
    }

    // Set to non-zero in-case species exists
    if (this_species_index not_eq -1) {
      line_vmrs[i] = atmospheric_vmrs[this_species_index];
    }
  }

  // Renormalize, if bath-species exist this is automatic.
  if (auto sl = sum(line_vmrs); bath) {
    line_vmrs[n - 1] = 1.0 - sl;
  } else if (sl == 0) {  // Special case
  } else {
    line_vmrs /= sl;
  }

  return line_vmrs;
}

Vector LineShape::vmrs(const AtmPoint& atm_point,
                       const ArrayOfSpeciesEnum& lineshape_species) ARTS_NOEXCEPT {
  const Index n = lineshape_species.size();

  // We need to know if bath is an actual species
  const bool bath = n and lineshape_species.back() == SpeciesEnum::Bath;

  // Initialize list of VMRS to 0
  Vector line_vmrs(n, 0);

  // Extract the VMR of the first species
  std::transform(
      lineshape_species.begin(),
      lineshape_species.end() - bath,
      line_vmrs.begin(),
      [&atm_point](const SpeciesEnum spec) { return atm_point[spec]; });

  // Renormalize, if bath-species exist this is automatic.
  if (auto sl = sum(line_vmrs); bath) {
    line_vmrs[n - 1] = 1.0 - sl;
  } else if (sl == 0) {  // Special case
  } else {
    line_vmrs /= sl;
  }

  return line_vmrs;
}

Vector LineShape::mass(const ConstVectorView& atmospheric_vmrs,
                       const ArrayOfArrayOfSpeciesTag& atmospheric_species,
                       const ArrayOfSpeciesEnum& lineshape_species,
                       const SpeciesIsotopologueRatios& ir) ARTS_NOEXCEPT {
  ARTS_ASSERT(
      atmospheric_species.size() == static_cast<Size>(atmospheric_vmrs.size()),
      "Bad atmospheric inputs");

  const Size n = lineshape_species.size();

  // Initialize list of VMRS to 0
  Vector line_vmrs(n, 0);
  Vector line_mass(n, 0);

  // We need to know if bath is an actual species
  const bool bath = n and lineshape_species.back() == SpeciesEnum::Bath;

  // Loop species ignoring self and bath
  for (Size i = 0; i < n - bath; i++) {
    // Select target in-case this is self-broadening
    const SpeciesEnum target = lineshape_species[i];

    // Find species in list or do nothing at all
    Index this_species_index = -1;
    for (Size j = 0; j < atmospheric_species.size(); j++) {
      if (atmospheric_species[j].Species() == target) {
        this_species_index = j;
      }
    }

    // Set to non-zero in-case species exists
    if (this_species_index not_eq -1) {
      line_vmrs[i] = atmospheric_vmrs[this_species_index];
      line_mass[i] = Species::mean_mass(target, ir);
    }
  }

  // Renormalize, if bath-species exist this is automatic.
  if (auto sl = sum(line_vmrs); sl == 0) {  // Special case
  } else if (bath) {
    line_mass[n - 1] = (line_vmrs * line_mass) / sl;
  }

  return line_mass;
}

namespace LineShape {
std::ostream& operator<<(std::ostream& os, const Model& m) {
  for (auto& data : m.Data()) os << data;
  return os;
}

std::istream& operator>>(std::istream& is, Model& m) {
  for (auto& data : m.Data()) is >> data;
  return is;
}

String ModelShape2MetaData(const Model& m) {
  String out;

  for (auto& var : enumtyps::LineShapeVariableOldTypes) {
    if (std::any_of(m.Data().cbegin(), m.Data().cend(), [var](auto& x) {
          return x.Get(var).type not_eq LineShapeTemperatureModelOld::None;
        })) {
      out += String(toString(var)) + ' ';
      for (auto& ssm : m.Data())
        out += String(toString(ssm.Get(var).type)) + ' ';
    }
  }

  if (out.size()) out.pop_back();

  return out;
}

Model MetaData2ModelShape(const String& s) {
  if (s.size() == 0) return {};

  const auto& names = enumstrs::LineShapeVariableOldNames<>;

  std::istringstream str(s);
  String part;
  LineShapeVariableOld var = LineShapeVariableOld::ETA;
  LineShapeTemperatureModelOld tm = LineShapeTemperatureModelOld::None;
  Index i = -100000;

  std::vector<SingleSpeciesModel> ssms(0);
  while (not str.eof()) {
    str >> part;
    if (std::any_of(names.cbegin(), names.cend(), [part](auto x) {
          return part == x;
        })) {
      i = -1;
      var = to<LineShapeVariableOld>(part);
    } else {
      i++;
      tm = to<LineShapeTemperatureModelOld>(part);
    }

    if (i < 0) continue;
    if (i < Index(ssms.size()))
      goto add_var;
    else {
      ssms.emplace_back();
    add_var:
      auto mp = ssms[i].Get(var);
      mp.type = tm;
      ssms[i].Set(var, mp);
    }
  }

  return Model(ssms);
}

String modelparameters2metadata(const ModelParameters mp, const Numeric T0) {
  std::ostringstream os;
  switch (mp.type) {
    case LineShapeTemperatureModelOld::None:
      os << 0;
      break;
    case LineShapeTemperatureModelOld::T0:
      os << mp.X0;
      break;
    case LineShapeTemperatureModelOld::T1:
      os << mp.X0 << " * (" << T0 << "/T)^" << mp.X1;
      break;
    case LineShapeTemperatureModelOld::T2:
      os << mp.X0 << " * (" << T0 << "/T)^" << mp.X1 << " / (1 + " << mp.X2
         << " * log(T/" << T0 << "))";
      break;
    case LineShapeTemperatureModelOld::T3:
      os << mp.X0 << " + " << mp.X1 << " * (" << T0 << " - T)";
      break;
    case LineShapeTemperatureModelOld::T4:
      os << "(" << mp.X0 << " + " << mp.X1 << " * (" << T0 << "/T - 1)) * ("
         << T0 << "/T)^" << mp.X2;
      break;
    case LineShapeTemperatureModelOld::T5:
      os << mp.X0 << " * (" << T0 << "/T)^(0.25 + 1.5 * " << mp.X1 << ")";
      break;
    case LineShapeTemperatureModelOld::LM_AER:
      os << '('
         << "Linear interpolation to y(x) from x-ref = [200, 250, 296, 340] and y-ref = ["
         << mp.X0 << ", " << mp.X1 << ", " << mp.X2 << ", " << mp.X3 << ']'
         << ')';
      break;
    case LineShapeTemperatureModelOld::DPL:
      os << '(' << mp.X0 << " * (" << T0 << "/T)^" << mp.X1 << " + " << mp.X2
         << " * (" << T0 << "/T)^" << mp.X3 << ')';
      break;
    case LineShapeTemperatureModelOld::POLY:
      os << '(' << mp.X0 << " + " << mp.X1 << " * T  + " << mp.X2
         << " * T * T + " << mp.X3 << " * T * T * T)";
  }

  return os.str();
}

ArrayOfString ModelMetaDataArray(const LineShape::Model&,
                                 const bool,
                                 const ArrayOfSpeciesEnum&,
                                 const Numeric) {
  // const auto& vars = enumtyps::VariableTypes;
  ArrayOfString as(0);

  // for (Index i = 0; i < Index(LineShapeVariableOld::FINAL); i++) {
  //   Variable var = vars[i];

  //   if (std::any_of(m.Data().cbegin(), m.Data().cend(), [var](auto& x) {
  //         return x.Get(var).type not_eq LineShapeTemperatureModelOld::None;
  //       })) {
  //     std::ostringstream os;
  //     os << var << " ~ ";
  //     for (Size j = 0; j < sts.size(); j++) {
  //       if (j == 0 and self)
  //         os << "VMR(" << self_broadening << ") * "
  //            << modelparameters2metadata(m.Data().front().Get(var), T0);
  //       else
  //         os << "VMR(" << toShortName(sts[j]) << ") * "
  //            << modelparameters2metadata(m.Data()[j].Get(var), T0);

  //       if (sts[j] not_eq sts.back()) os << " + ";
  //     }
  //     as.push_back(os.str());
  //   }
  // }

  return as;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
Numeric& SingleModelParameter(ModelParameters& mp, const String& type) {
  if (type == "X0") return mp.X0;
  if (type == "X1") return mp.X1;
  if (type == "X2") return mp.X2;
  if (type == "X3") return mp.X3;
  ARTS_USER_ERROR("Type: ",
                  type,
                  ", is not accepted.  "
                  "See documentation for accepted types\n")
}
#pragma GCC diagnostic pop

std::ostream& operator<<(std::ostream& os, const ModelParameters& mp) {
  os << mp.type << ' ' << mp.X0 << ' ' << mp.X1 << ' ' << mp.X2 << ' ' << mp.X3
     << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, ModelParameters& mp) {
  is >> mp.type >> double_imanip() >> mp.X0 >> mp.X1 >> mp.X2 >> mp.X3;
  return is;
}

Numeric ModelParameters::at(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return X0;
    case LineShapeTemperatureModelOld::T1:
      return X0 * pow(T0 / T, X1);
    case LineShapeTemperatureModelOld::T2:
      return X0 * pow(T0 / T, X1) * (1 + X2 * log(T / T0));
    case LineShapeTemperatureModelOld::T3:
      return X0 + X1 * (T - T0);
    case LineShapeTemperatureModelOld::T4:
      return (X0 + X1 * (T0 / T - 1.)) * pow(T0 / T, X2);
    case LineShapeTemperatureModelOld::T5:
      return X0 * pow(T0 / T, 0.25 + 1.5 * X1);
    case LineShapeTemperatureModelOld::LM_AER:
      return special_linemixing_aer(T);
    case LineShapeTemperatureModelOld::DPL:
      return X0 * pow(T0 / T, X1) + X2 * pow(T0 / T, X3);
    case LineShapeTemperatureModelOld::POLY:
      return X0 + X1 * T + X2 * T * T + X3 * T * T * T;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric ModelParameters::dX0(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return 1;
    case LineShapeTemperatureModelOld::T1:
      return pow(T0 / T, X1);
    case LineShapeTemperatureModelOld::T2:
      return pow(T0 / T, X1) * (1 + X2 * log(T / T0));
    case LineShapeTemperatureModelOld::T3:
      return 1;
    case LineShapeTemperatureModelOld::T4:
      return pow(T0 / T, X2);
    case LineShapeTemperatureModelOld::T5:
      return pow(T0 / T, 1.5 * X1 + 0.25);
    case LineShapeTemperatureModelOld::LM_AER:
      return special_linemixing_aer_dX0(T);
    case LineShapeTemperatureModelOld::DPL:
      return pow(T0 / T, X1);
    case LineShapeTemperatureModelOld::POLY:
      return 1;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric ModelParameters::dX1(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return 0;
    case LineShapeTemperatureModelOld::T1:
      return X0 * pow(T0 / T, X1) * log(T0 / T);
    case LineShapeTemperatureModelOld::T2:
      return X0 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) * log(T0 / T);
    case LineShapeTemperatureModelOld::T3:
      return (T - T0);
    case LineShapeTemperatureModelOld::T4:
      return pow(T0 / T, X2) * (T0 / T - 1.);
    case LineShapeTemperatureModelOld::T5:
      return 1.5 * X0 * pow(T0 / T, 1.5 * X1 + 0.25) * log(T0 / T);
    case LineShapeTemperatureModelOld::LM_AER:
      return special_linemixing_aer_dX1(T);
    case LineShapeTemperatureModelOld::DPL:
      return X0 * pow(T0 / T, X1) * log(T0 / T);
    case LineShapeTemperatureModelOld::POLY:
      return T;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric ModelParameters::dX2(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return 0;
    case LineShapeTemperatureModelOld::T1:
      return 0;
    case LineShapeTemperatureModelOld::T2:
      return X0 * pow(T0 / T, X1) * log(T / T0);
    case LineShapeTemperatureModelOld::T3:
      return 0;
    case LineShapeTemperatureModelOld::T4:
      return pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1)) * log(T0 / T);
    case LineShapeTemperatureModelOld::T5:
      return 0;
    case LineShapeTemperatureModelOld::LM_AER:
      return special_linemixing_aer_dX2(T);
    case LineShapeTemperatureModelOld::DPL:
      return pow(T0 / T, X3);
    case LineShapeTemperatureModelOld::POLY:
      return T * T;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric ModelParameters::dX3(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return 0;
    case LineShapeTemperatureModelOld::T1:
      return 0;
    case LineShapeTemperatureModelOld::T2:
      return 0;
    case LineShapeTemperatureModelOld::T3:
      return 0;
    case LineShapeTemperatureModelOld::T4:
      return 0;
    case LineShapeTemperatureModelOld::T5:
      return 0;
    case LineShapeTemperatureModelOld::LM_AER:
      return special_linemixing_aer_dX3(T);
    case LineShapeTemperatureModelOld::DPL:
      return X2 * pow(T0 / T, X3) * log(T0 / T);
    case LineShapeTemperatureModelOld::POLY:
      return T * T * T;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric ModelParameters::dT(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return 0;
    case LineShapeTemperatureModelOld::T1:
      return -X0 * X1 * pow(T0 / T, X1) / T;
    case LineShapeTemperatureModelOld::T2:
      return -X0 * X1 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) / T +
             X0 * X2 * pow(T0 / T, X1) / T;
    case LineShapeTemperatureModelOld::T3:
      return X1;
    case LineShapeTemperatureModelOld::T4:
      return -X2 * pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1.)) / T -
             T0 * X1 * pow(T0 / T, X2) / pow(T, 2);
    case LineShapeTemperatureModelOld::T5:
      return -X0 * pow(T0 / T, 1.5 * X1 + 0.25) * (1.5 * X1 + 0.25) / T;
    case LineShapeTemperatureModelOld::LM_AER:
      return special_linemixing_aer_dT(T);
    case LineShapeTemperatureModelOld::DPL:
      return -X0 * X1 * pow(T0 / T, X1) / T + -X2 * X3 * pow(T0 / T, X3) / T;
    case LineShapeTemperatureModelOld::POLY:
      return X1 + 2 * X2 * T + 3 * X3 * T * T;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

Numeric ModelParameters::dT0(Numeric T, Numeric T0) const noexcept {
  using std::log;
  using std::pow;

  switch (type) {
    case LineShapeTemperatureModelOld::None:
      return 0;
    case LineShapeTemperatureModelOld::T0:
      return 0;
    case LineShapeTemperatureModelOld::T1:
      return X0 * X1 * pow(T0 / T, X1) / T0;
    case LineShapeTemperatureModelOld::T2:
      return X0 * X1 * pow(T0 / T, X1) * (X2 * log(T / T0) + 1.) / T0 -
             X0 * X2 * pow(T0 / T, X1) / T0;
    case LineShapeTemperatureModelOld::T3:
      return -X1;
    case LineShapeTemperatureModelOld::T4:
      return X2 * pow(T0 / T, X2) * (X0 + X1 * (T0 / T - 1.)) / T0 +
             X1 * pow(T0 / T, X2) / T;
    case LineShapeTemperatureModelOld::T5:
      return X0 * pow(T0 / T, 1.5 * X1 + 0.25) * (1.5 * X1 + 0.25) / T0;
    case LineShapeTemperatureModelOld::LM_AER:
      return 0;
    case LineShapeTemperatureModelOld::DPL:
      return X0 * X1 * pow(T0 / T, X1) / T0 + X2 * X3 * pow(T0 / T, X3) / T0;
    case LineShapeTemperatureModelOld::POLY:
      return 0;
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

bifstream& SingleSpeciesModel::read(bifstream& bif) {
  for (auto& data : X) {
    Index x;
    bif >> x >> data.X0 >> data.X1 >> data.X2 >> data.X3;
    data.type = LineShapeTemperatureModelOld(x);
  }
  return bif;
}

bofstream& SingleSpeciesModel::write(bofstream& bof) const {
  for (auto& data : X) {
    bof << Index(data.type) << data.X0 << data.X1 << data.X2 << data.X3;
  }
  return bof;
}

Output SingleSpeciesModel::at(Numeric T, Numeric T0, Numeric P) const noexcept {
  static_assert(nVars == 9);
  return {P * G0().at(T, T0),
          P * D0().at(T, T0),
          P * G2().at(T, T0),
          P * D2().at(T, T0),
          P * FVC().at(T, T0),
          ETA().at(T, T0),
          P * Y().at(T, T0),
          P * P * G().at(T, T0),
          P * P * DV().at(T, T0)};
}

Output SingleSpeciesModel::dT(Numeric T, Numeric T0, Numeric P) const noexcept {
  static_assert(nVars == 9);
  return {P * G0().dT(T, T0),
          P * D0().dT(T, T0),
          P * G2().dT(T, T0),
          P * D2().dT(T, T0),
          P * FVC().dT(T, T0),
          ETA().dT(T, T0),
          P * Y().dT(T, T0),
          P * P * G().dT(T, T0),
          P * P * DV().dT(T, T0)};
}

Output SingleSpeciesModel::dT0(Numeric T,
                               Numeric T0,
                               Numeric P) const noexcept {
  static_assert(nVars == 9);
  return {P * G0().dT0(T, T0),
          P * D0().dT0(T, T0),
          P * G2().dT0(T, T0),
          P * D2().dT0(T, T0),
          P * FVC().dT0(T, T0),
          ETA().dT0(T, T0),
          P * Y().dT0(T, T0),
          P * P * G().dT0(T, T0),
          P * P * DV().dT0(T, T0)};
}

#define FUNC(X, PVAR)                                                        \
  Numeric Model::X(                                                          \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Vector& vmrs) \
      const ARTS_NOEXCEPT {                                                  \
    ARTS_ASSERT(size() == static_cast<Size>(vmrs.size()))                    \
                                                                             \
    return PVAR * std::transform_reduce(begin(),                             \
                                        end(),                               \
                                        vmrs.data_handle(),                  \
                                        0.0,                                 \
                                        std::plus<>(),                       \
                                        [T, T0](auto& ls, auto& xi) {        \
                                          return xi * ls.X().at(T, T0);      \
                                        });                                  \
  }
FUNC(G0, P)
FUNC(D0, P)
FUNC(G2, P)
FUNC(D2, P)
FUNC(ETA, 1)
FUNC(FVC, P)
FUNC(Y, P)
FUNC(G, P* P)
FUNC(DV, P* P)
#undef FUNC

#define FUNC(X, PVAR)                                                        \
  Numeric Model::d##X##dT(                                                   \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Vector& vmrs) \
      const ARTS_NOEXCEPT {                                                  \
    ARTS_ASSERT(size() == static_cast<Size>(vmrs.size()))                    \
                                                                             \
    return PVAR * std::transform_reduce(begin(),                             \
                                        end(),                               \
                                        vmrs.data_handle(),                  \
                                        0.0,                                 \
                                        std::plus<>(),                       \
                                        [T, T0](auto& ls, auto& xi) {        \
                                          return xi * ls.X().dT(T, T0);      \
                                        });                                  \
  }
FUNC(G0, P)
FUNC(D0, P)
FUNC(G2, P)
FUNC(D2, P)
FUNC(ETA, 1)
FUNC(FVC, P)
FUNC(Y, P)
FUNC(G, P* P)
FUNC(DV, P* P)
#undef FUNC

std::pair<bool, bool> SingleSpeciesModel::MatchTypes(
    const SingleSpeciesModel& other) const noexcept {
  bool match = true, nullable = true;
  for (Index i = 0; i < nVars; i++) {
    match = match and other.X[i].type == X[i].type;
    nullable = nullable and
               (modelparameterEmpty(other.X[i]) or modelparameterEmpty(X[i]) or
                other.X[i].type == X[i].type);
  }
  return {match, nullable};
}

std::ostream& operator<<(std::ostream& os, const SingleSpeciesModel& ssm) {
  for (const auto& mp : ssm.Data())
    if (mp.type not_eq LineShapeTemperatureModelOld::None)
      os << mp.X0 << ' ' << mp.X1 << ' ' << mp.X2 << ' ' << mp.X3 << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, SingleSpeciesModel& ssm) {
  for (auto& mp : ssm.Data())
    if (mp.type not_eq LineShapeTemperatureModelOld::None)
      is >> double_imanip() >> mp.X0 >> mp.X1 >> mp.X2 >> mp.X3;
  return is;
}

std::ostream& operator<<(std::ostream& os, Output x) {
  return os << "G0: " << x.G0 << " D0: " << x.D0 << " G2: " << x.G2
            << " D2: " << x.D2 << " FVC: " << x.FVC << " ETA: " << x.ETA
            << " Y: " << x.Y << " G: " << x.G << " DV: " << x.DV;
}

Model::Model(Numeric sgam,
             Numeric nself,
             Numeric agam,
             Numeric nair,
             Numeric psf,
             std::array<Numeric, 12> aer_interp) noexcept
    : mdata(2) {
  mdata.front().G0() = {LineShapeTemperatureModelOld::T1, sgam, nself, 0, 0};
  mdata.front().D0() = {LineShapeTemperatureModelOld::T5, psf, nair, 0, 0};

  mdata.back().G0() = {LineShapeTemperatureModelOld::T1, agam, nair, 0, 0};
  mdata.back().D0() = {LineShapeTemperatureModelOld::T5, psf, nair, 0, 0};

  if (std::any_of(aer_interp.cbegin(), aer_interp.cend(), [](auto x) {
        return x not_eq 0;
      })) {
    mdata.front().Y().type = LineShapeTemperatureModelOld::LM_AER;
    mdata.front().Y().X0 = aer_interp[4];
    mdata.front().Y().X1 = aer_interp[5];
    mdata.front().Y().X2 = aer_interp[6];
    mdata.front().Y().X3 = aer_interp[7];
    mdata.front().G().type = LineShapeTemperatureModelOld::LM_AER;
    mdata.front().G().X0 = aer_interp[8];
    mdata.front().G().X1 = aer_interp[9];
    mdata.front().G().X2 = aer_interp[10];
    mdata.front().G().X3 = aer_interp[11];

    mdata.back().Y().type = LineShapeTemperatureModelOld::LM_AER;
    mdata.back().Y().X0 = aer_interp[4];
    mdata.back().Y().X1 = aer_interp[5];
    mdata.back().Y().X2 = aer_interp[6];
    mdata.back().Y().X3 = aer_interp[7];
    mdata.back().G().type = LineShapeTemperatureModelOld::LM_AER;
    mdata.back().G().X0 = aer_interp[8];
    mdata.back().G().X1 = aer_interp[9];
    mdata.back().G().X2 = aer_interp[10];
    mdata.back().G().X3 = aer_interp[11];
  }
}

Model hitran_model(
    Numeric sgam, Numeric nself, Numeric agam, Numeric nair, Numeric psf) {
  Model m(2);

  m.Data().front().G0() = {LineShapeTemperatureModelOld::T1, sgam, nself, 0, 0};
  m.Data().front().D0() = {LineShapeTemperatureModelOld::T0, psf, 0, 0, 0};

  m.Data().back().G0() = {LineShapeTemperatureModelOld::T1, agam, nair, 0, 0};
  m.Data().back().D0() = {LineShapeTemperatureModelOld::T0, psf, 0, 0, 0};

  return m;
}

Model lblrtm_model(Numeric sgam,
                   Numeric nself,
                   Numeric agam,
                   Numeric nair,
                   Numeric psf,
                   std::array<Numeric, 12> aer_interp) {
  Model m(2);

  m.Data().front().G0() = {LineShapeTemperatureModelOld::T1, sgam, nself, 0, 0};
  m.Data().front().D0() = {LineShapeTemperatureModelOld::T0, psf, 0, 0, 0};

  m.Data().back().G0() = {LineShapeTemperatureModelOld::T1, agam, nair, 0, 0};
  m.Data().back().D0() = {LineShapeTemperatureModelOld::T0, psf, 0, 0, 0};

  if (std::any_of(aer_interp.cbegin(), aer_interp.cend(), [](auto x) {
        return x not_eq 0;
      })) {
    m.Data().front().Y().type = LineShapeTemperatureModelOld::LM_AER;
    m.Data().front().Y().X0 = aer_interp[4];
    m.Data().front().Y().X1 = aer_interp[5];
    m.Data().front().Y().X2 = aer_interp[6];
    m.Data().front().Y().X3 = aer_interp[7];
    m.Data().front().G().type = LineShapeTemperatureModelOld::LM_AER;
    m.Data().front().G().X0 = aer_interp[8];
    m.Data().front().G().X1 = aer_interp[9];
    m.Data().front().G().X2 = aer_interp[10];
    m.Data().front().G().X3 = aer_interp[11];

    m.Data().back().Y().type = LineShapeTemperatureModelOld::LM_AER;
    m.Data().back().Y().X0 = aer_interp[4];
    m.Data().back().Y().X1 = aer_interp[5];
    m.Data().back().Y().X2 = aer_interp[6];
    m.Data().back().Y().X3 = aer_interp[7];
    m.Data().back().G().type = LineShapeTemperatureModelOld::LM_AER;
    m.Data().back().G().X0 = aer_interp[8];
    m.Data().back().G().X1 = aer_interp[9];
    m.Data().back().G().X2 = aer_interp[10];
    m.Data().back().G().X3 = aer_interp[11];
  }

  return m;
}

bool Model::OK(LineShapeTypeOld type,
               bool self,
               bool bath,
               const std::size_t nspecies) const noexcept {
  Index n = mdata.size();
  Index m = Index(self) + Index(bath);
  bool needs_any = type not_eq LineShapeTypeOld::DP;
  return not(n not_eq Index(nspecies) or m > n or (needs_any and not n));
}

void Model::Remove(Index i, ArrayOfSpeciesTag& specs) {
  mdata.erase(mdata.begin() + i);
  specs.erase(specs.begin() + i);
}

void Model::Remove(Index i, ArrayOfSpeciesEnum& specs) {
  mdata.erase(mdata.begin() + i);
  specs.erase(specs.begin() + i);
}

void Model::SetLineMixingModel(SingleSpeciesModel x) {
  for (auto& ssm : mdata) {
    ssm.Y() = x.Y();
    ssm.G() = x.G();
    ssm.DV() = x.DV();
  }
}

std::pair<bool, bool> Model::Match(const Model& other) const noexcept {
  const Size n = size();
  if (other.size() not_eq n) return {false, false};

  bool match = true, nullable = true;
  for (Size i = 0; i < n; i++) {
    const auto x = mdata[i].MatchTypes(other[i]);
    match = match and x.first;
    nullable = nullable and x.second;
  }

  return {match, nullable};
}

bifstream& Model::read(bifstream& bif) {
  for (auto& data : mdata) data.read(bif);
  return bif;
}

bofstream& Model::write(bofstream& bof) const {
  for (auto& data : mdata) data.write(bof);
  return bof;
}

namespace LegacyLineFunctionData {
std::vector<LineShapeVariableOld> lineshapetag2variablesvector(String type) {
  if (type == String("DP")) return {};
  if (type == String("LP")) return {LineShapeVariableOld::G0, LineShapeVariableOld::D0};
  if (type == String("VP")) return {LineShapeVariableOld::G0, LineShapeVariableOld::D0};
  if (type == String("SDVP"))
    return {LineShapeVariableOld::G0, LineShapeVariableOld::D0, LineShapeVariableOld::G2, LineShapeVariableOld::D2};
  if (type == String("HTP"))
    return {LineShapeVariableOld::G0,
            LineShapeVariableOld::D0,
            LineShapeVariableOld::G2,
            LineShapeVariableOld::D2,
            LineShapeVariableOld::FVC,
            LineShapeVariableOld::ETA};
  ARTS_USER_ERROR("Type: ",
                  type,
                  ", is not accepted.  "
                  "See documentation for accepted types\n")
}

std::vector<LineShapeVariableOld> linemixingtag2variablesvector(String type) {
  if (type == "#") return {};
  if (type == "LM1") return {LineShapeVariableOld::Y};
  if (type == "LM2") return {LineShapeVariableOld::Y, LineShapeVariableOld::G, LineShapeVariableOld::DV};
  if (type == "INT") return {};
  if (type == "ConstG") return {LineShapeVariableOld::G};
  ARTS_USER_ERROR("Type: ",
                  type,
                  ", is not accepted.  "
                  "See documentation for accepted types\n")
}
}  // namespace LegacyLineFunctionData

namespace LegacyLineMixingData {
LegacyLineMixingData::TypeLM string2typelm(String type) {
  if (type == "NA")  // The standard case
    return TypeLM::LM_NONE;
  if (type == "LL")  // The LBLRTM case
    return TypeLM::LM_LBLRTM;
  if (type == "NR")  // The LBLRTM O2 non-resonant case
    return TypeLM::LM_LBLRTM_O2NonResonant;
  if (type == "L2")  // The 2nd order case
    return TypeLM::LM_2NDORDER;
  if (type == "L1")  // The 2nd order case
    return TypeLM::LM_1STORDER;
  if (type == "BB")  // The band class
    return TypeLM::LM_BYBAND;
  ARTS_USER_ERROR("Type: ",
                  type,
                  ", is not accepted.  "
                  "See documentation for accepted types\n")
}
}  // namespace LegacyLineMixingData

namespace LegacyPressureBroadeningData {
LegacyPressureBroadeningData::TypePB string2typepb(String type) {
  if (type == "NA")  // The none case
    return TypePB::PB_NONE;
  if (type == "N2")  // Air Broadening is N2 broadening mostly...
    return TypePB::PB_AIR_BROADENING;
  if (type == "WA")  // Water and Air Broadening
    return TypePB::PB_AIR_AND_WATER_BROADENING;
  if (type == "AP")  // Planetary broadening
    return TypePB::PB_PLANETARY_BROADENING;
  ARTS_USER_ERROR("Type: ",
                  type,
                  ", is not accepted.  "
                  "See documentation for accepted types\n")
}

Index self_listed(const QuantumIdentifier& qid,
                  LegacyPressureBroadeningData::TypePB t) {
  if (t == TypePB::PB_PLANETARY_BROADENING and
      (qid.Species() == to<SpeciesEnum>("N2") or
       qid.Species() == to<SpeciesEnum>("O2") or
       qid.Species() == to<SpeciesEnum>("H2O") or
       qid.Species() == to<SpeciesEnum>("CO2") or
       qid.Species() == to<SpeciesEnum>("H2") or
       qid.Species() == to<SpeciesEnum>("He")))
    return true;
  if (t == TypePB::PB_AIR_AND_WATER_BROADENING and
      qid.Species() == to<SpeciesEnum>("H2O"))
    return true;
  return false;
}
}  // namespace LegacyPressureBroadeningData
}  // namespace LineShape
