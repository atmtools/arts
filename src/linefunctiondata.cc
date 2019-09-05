/* Copyright (C) 2018
 Richard Larsson <larsson@mps.mpg.de>

 
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any
 later version.

 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 USA. */

/** Contains the line function data class
 * @file   linefunctiondata.cc
 * @author Richard Larsson
 * @date   2018-09-19
 * 
 * @brief Implementations of linefunctiondata.h
 * 
 * This mostly contains functions that either did not
 * compile while placed in the header or for other 
 * reasons are not there.  This means most of the 
 * real work happens in the header file and not here
 * 
 * FIXME:  Rename this and its header to match newer content, lineshapemodel.cc/h
 **/

#include "linefunctiondata.h"

//! {"X0", "X1", "X2"}
ArrayOfString all_coefficientsLineFunctionData() { return {"X0", "X1", "X2"}; }

//! {"G0", "D0", "G2", "D2", "ETA", "FVC", "Y", "G", "DV"}
ArrayOfString all_variablesLineFunctionData() {
  return {"G0", "D0", "G2", "D2", "FVC", "ETA", "Y", "G", "DV"};
}

JacPropMatType select_derivativeLineShape(const String& var,
                                          const String& coeff) {
  // Test viability of model variables
  static const ArrayOfString vars = all_variablesLineFunctionData();
  bool var_OK = false;
  for (auto& v : vars)
    if (var == v) var_OK = true;

  // Test viability of model coefficients
  static const ArrayOfString coeffs = all_coefficientsLineFunctionData();
  bool coeff_OK = false;
  for (auto& c : coeffs)
    if (coeff == c) coeff_OK = true;

  // Fails either when the user has bad input or when the developer fails to update all_variablesLineFunctionData or all_coefficientsLineFunctionData
  if (not var_OK or not coeff_OK) {
    std::ostringstream os;
    os << "At least one of your variable and/or your coefficient is not OK\n";
    os << "Your variable: \"" << var << "\".  OK variables include: " << vars
       << "\n";
    os << "Your coefficient: \"" << coeff
       << "\".  OK coefficients include: " << coeffs << "\n";
    throw std::runtime_error(os.str());
  }

// Define a repetitive pattern.  Update if/when there are more coefficients in the future
#define ReturnJacPropMatType(ID)                \
  (var == #ID) {                                \
    if (coeff == "X0")                          \
      return JacPropMatType::LineShape##ID##X0; \
    else if (coeff == "X1")                     \
      return JacPropMatType::LineShape##ID##X1; \
    else if (coeff == "X2")                     \
      return JacPropMatType::LineShape##ID##X2; \
  }

  if
    ReturnJacPropMatType(G0) else if ReturnJacPropMatType(D0) else if ReturnJacPropMatType(G2) else if ReturnJacPropMatType(D2) else if ReturnJacPropMatType(
        FVC) else if ReturnJacPropMatType(ETA) else if ReturnJacPropMatType(Y) else if ReturnJacPropMatType(G) else if ReturnJacPropMatType(DV)
#undef ReturnJacPropMatType

        std::terminate();
}

std::istream& LineShape::from_artscat4(std::istream& is,
                                       Model& m,
                                       const QuantumIdentifier& qid) {
  // Set or reset variables
  m.mtype = Type::VP;
  m.mself = true;
  m.mbath = false;
  m.mdata = std::vector<SingleSpeciesModel>(7);
  m.mspecies = ArrayOfSpeciesTag(7);

  // Set species
  m.mspecies[1] = SpeciesTag("N2");
  m.mspecies[2] = SpeciesTag("O2");
  m.mspecies[3] = SpeciesTag("H2O");
  m.mspecies[4] = SpeciesTag("CO2");
  m.mspecies[5] = SpeciesTag("H2");
  m.mspecies[6] = SpeciesTag("He");

  // Temperature types
  for (auto& v : m.mdata) {
    v.G0().type = TemperatureModel::T1;
    v.D0().type = TemperatureModel::T5;
  }

  // G0 main coefficient
  for (auto& v : m.mdata) is >> v.G0().X0;

  // G0 exponent is same as D0 exponent
  for (auto& v : m.mdata) {
    is >> v.G0().X1;
    v.D0().X1 = v.G0().X1;
  }

  // D0 coefficient
  m.mdata[0].D0().X0 = 0;
  for (int k = 1; k < 7; k++)
    is >> m.mdata[k].D0().X0;

  // Special case when self is part of this list, it needs to be removed
  for (int k = 1; k < 7; k++) {
    if (qid.Species() == m.mspecies[k].Species()) {
      if(m.mdata[0].G0().X0 not_eq m.mdata[k].G0().X0 or
         m.mdata[0].G0().X1 not_eq m.mdata[k].G0().X1 or
         m.mdata[0].D0().X1 not_eq m.mdata[k].D0().X1) {
        std::ostringstream os;
        os << "Species is " << qid.SpeciesName() << " and this is a broadening species in ARTSCAT-4.\n"
           << "Despite this, values representing self and " << qid.SpeciesName() << " does not match "
           << "in input string\n";
        throw std::runtime_error(os.str());
      }
      m.mdata[0].D0().X0 = m.mdata[k].D0().X0;
      m.Remove(k);
      break;
    }
  }
  
  return is;
}

std::istream& LineShape::from_linefunctiondata(std::istream& data, Model& m) {
  m.mself = m.mbath = false;
  Index specs;
  String s;

  // The first tag should give the line shape scheme
  data >> s;
  m.mtype = LineShape::string2shapetype(s);

  // Order of elements for line shape
  const auto shapeparams =
      LegacyLineFunctionData::lineshapetag2variablesvector(s);

  // The second tag should give the line mixing scheme
  data >> s;

  // Order of elements for line mixing
  const auto mixingparams =
      LegacyLineFunctionData::linemixingtag2variablesvector(s);

  // The third tag should contain the number of species
  data >> specs;
  m.mspecies.resize(specs);
  m.mdata.resize(specs);

  if (not specs and m.mtype not_eq Type::DP)
    throw std::runtime_error(
        "Need at least one species for non-Doppler line shapes");

  // For all species, we need to set the methods to compute them
  for (Index i = 0; i < specs; i++) {
    // This should be a species tag or one of the specials, SELF or BATH
    data >> s;
    if (s == self_broadening) {
      // If the species is self, then  we need to flag this
      m.mself = true;
      if (i not_eq 0)  // but self has to be first for consistent behavior
        throw std::runtime_error("Self broadening must be first, it is not\n");
    }

    else if (s == bath_broadening) {
      // If the species is air, then we need to flag this
      m.mbath = true;
      if (i not_eq
          specs - 1)  // but air has to be last because it needs the rest's VMR
        throw std::runtime_error(
            "Air/bath broadening must be last, it is not\n");
    } else {
      // Otherwise, we hope we find a species
      try {
        m.mspecies[i] = SpeciesTag(s);
      } catch (const std::runtime_error& e) {
        ostringstream os;
        os << "Encountered " << s
           << " in a position where a species should have been ";
        os << "defined.\nPlease check your pressure broadening data structure and ensure ";
        os << "that it follows the correct conventions.\n";
        os << "SpeciesTag error reads:  " << e.what();
        throw std::runtime_error(os.str());
      }
    }

    // For all parameters
    for (auto& params : {shapeparams, mixingparams}) {
      for (auto& param : params) {
        data >> s;  // Should contain a temperature tag

        const auto type = string2temperaturemodel(s);
        const Index ntemp =
            LegacyLineFunctionData::temperaturemodel2legacynelem(type);

        m.mdata[i].Data()[Index(param)].type = type;
        if (ntemp <= nmaxTempModelParams) {
          switch (ntemp) {
            case 1:
              data >> m.mdata[i].Data()[Index(param)].X0;
              break;
            case 2:
              data >> m.mdata[i].Data()[Index(param)].X0 >>
                  m.mdata[i].Data()[Index(param)].X1;
              break;
            case 3:
              data >> m.mdata[i].Data()[Index(param)].X0 >>
                  m.mdata[i].Data()[Index(param)].X1 >>
                  m.mdata[i].Data()[Index(param)].X2;
              break;
            case 0:
              break;
            default:
              throw std::runtime_error(
                  "Unknown number of input parameters in Legacy mode.");
          }
        } else {  // Has to be the only allowed interpolation case
          if (ntemp > nmaxInterpModels)
            throw std::runtime_error(
                "Too many input parameters in interpolation results Legacy mode.");
          for (Index k = 0; k < ntemp; k++) data >> m.mdata[i].Interp()[k];
        }
      }
    }
  }

  return data;
}

std::istream& LineShape::from_pressurebroadeningdata(
    std::istream& data, LineShape::Model& lsc, const QuantumIdentifier& qid) {
  String s;
  data >> s;

  const auto type = LegacyPressureBroadeningData::string2typepb(s);
  const auto n = LegacyPressureBroadeningData::typepb2nelem(type);
  const auto self_in_list =
      LegacyPressureBroadeningData::self_listed(qid, type);

  Vector x(n);
  for (auto& num : x) data >> num;

  lsc = LegacyPressureBroadeningData::vector2modelpb(x, type, self_in_list);

  return data;
}

std::istream& LineShape::from_linemixingdata(std::istream& data,
                                             LineShape::Model& lsc) {
  String s;
  data >> s;

  const auto type = LegacyLineMixingData::string2typelm(s);
  const auto n = LegacyLineMixingData::typelm2nelem(type);

  Vector x(n);
  for (auto& num : x) data >> num;

  lsc = LegacyLineMixingData::vector2modellm(x, type);

  return data;
}

LineShape::Model LineShape::LegacyPressureBroadeningData::vector2modelpb(
    Vector x,
    LineShape::LegacyPressureBroadeningData::TypePB type,
    bool self_in_list) {
  switch (type) {
    case TypePB::PB_NONE:
      return Model();
    case TypePB::PB_AIR_BROADENING:
      return Model(x[0], x[1], x[2], x[3], x[4]);
    case TypePB::PB_AIR_AND_WATER_BROADENING:
      if (self_in_list) {
        ArrayOfSpeciesTag spec(2);
        spec[0] = SpeciesTag("H2O");
        std::vector<SingleSpeciesModel> ssm(2);
        ssm[0].G0() = {TemperatureModel::T1, x[0], x[1], 0};
        ssm[0].D0() = {TemperatureModel::T5, x[2], x[1], 0};
        ssm[1].G0() = {TemperatureModel::T1, x[3], x[4], 0};
        ssm[1].D0() = {TemperatureModel::T5, x[5], x[4], 0};
        return Model(LineShape::Type::VP, false, true, spec, ssm);
      } else {
        ArrayOfSpeciesTag spec(3);
        spec[1] = SpeciesTag("H2O");
        std::vector<SingleSpeciesModel> ssm(3);
        ssm[0].G0() = {TemperatureModel::T1, x[0], x[1], 0};
        ssm[0].D0() = {TemperatureModel::T5, x[2], x[1], 0};
        ssm[2].G0() = {TemperatureModel::T1, x[3], x[4], 0};
        ssm[2].D0() = {TemperatureModel::T5, x[5], x[4], 0};
        ssm[1].G0() = {TemperatureModel::T1, x[6], x[7], 0};
        ssm[1].D0() = {TemperatureModel::T5, x[8], x[7], 0};
        return Model(LineShape::Type::VP, true, true, spec, ssm);
      }
    case TypePB::PB_PLANETARY_BROADENING:
      if (self_in_list) {
        const ArrayOfSpeciesTag spec = {SpeciesTag(String("N2")),
                                        SpeciesTag(String("O2")),
                                        SpeciesTag(String("H2O")),
                                        SpeciesTag(String("CO2")),
                                        SpeciesTag(String("H2")),
                                        SpeciesTag(String("He"))};
        std::vector<SingleSpeciesModel> ssm(6);
        ssm[0].G0() = {TemperatureModel::T1, x[1], x[8], 0};
        ssm[0].D0() = {TemperatureModel::T5, x[14], x[8], 0};
        ssm[1].G0() = {TemperatureModel::T1, x[2], x[9], 0};
        ssm[1].D0() = {TemperatureModel::T5, x[15], x[9], 0};
        ssm[2].G0() = {TemperatureModel::T1, x[3], x[10], 0};
        ssm[2].D0() = {TemperatureModel::T5, x[16], x[10], 0};
        ssm[3].G0() = {TemperatureModel::T1, x[4], x[11], 0};
        ssm[3].D0() = {TemperatureModel::T5, x[17], x[11], 0};
        ssm[4].G0() = {TemperatureModel::T1, x[5], x[12], 0};
        ssm[4].D0() = {TemperatureModel::T5, x[18], x[12], 0};
        ssm[5].G0() = {TemperatureModel::T1, x[6], x[13], 0};
        ssm[5].D0() = {TemperatureModel::T5, x[19], x[13], 0};
        return Model(LineShape::Type::VP, false, false, spec, ssm);
      } else {
        ArrayOfSpeciesTag spec(7);
        spec[1] = SpeciesTag(String("N2"));
        spec[2] = SpeciesTag(String("O2"));
        spec[3] = SpeciesTag(String("H2O"));
        spec[4] = SpeciesTag(String("CO2"));
        spec[5] = SpeciesTag(String("H2"));
        spec[6] = SpeciesTag(String("He"));
        std::vector<SingleSpeciesModel> ssm(7);
        ssm[0].G0() = {TemperatureModel::T1, x[0], x[7], 0};
        //          ssm[0].D0() = ...
        ssm[1].G0() = {TemperatureModel::T1, x[1], x[8], 0};
        ssm[1].D0() = {TemperatureModel::T5, x[14], x[8], 0};
        ssm[2].G0() = {TemperatureModel::T1, x[2], x[9], 0};
        ssm[2].D0() = {TemperatureModel::T5, x[15], x[9], 0};
        ssm[3].G0() = {TemperatureModel::T1, x[3], x[10], 0};
        ssm[3].D0() = {TemperatureModel::T5, x[16], x[10], 0};
        ssm[4].G0() = {TemperatureModel::T1, x[4], x[11], 0};
        ssm[4].D0() = {TemperatureModel::T5, x[17], x[11], 0};
        ssm[5].G0() = {TemperatureModel::T1, x[5], x[12], 0};
        ssm[5].D0() = {TemperatureModel::T5, x[18], x[12], 0};
        ssm[6].G0() = {TemperatureModel::T1, x[6], x[13], 0};
        ssm[6].D0() = {TemperatureModel::T5, x[19], x[13], 0};
        return Model(LineShape::Type::VP, true, false, spec, ssm);
      }
  }
  std::terminate();
}

LineShape::Model LineShape::LegacyLineMixingData::vector2modellm(
    Vector x, LineShape::LegacyLineMixingData::TypeLM type) {
  auto y = Model();
  y.resize(1);
  switch (type) {
    case TypeLM::LM_NONE:
      break;
    case TypeLM::LM_LBLRTM:
      y.Data()[0].Y().type = LineShape::TemperatureModel::LM_AER;
      y.Data()[0].G().type = LineShape::TemperatureModel::LM_AER;
      std::copy(x.begin(), x.end(), y.Data()[0].Interp().begin());
      break;
    case TypeLM::LM_LBLRTM_O2NonResonant:
      y.Data()[0].G().type = LineShape::TemperatureModel::T0;
      y.Data()[0].G().X0 = x[0];
      break;
    case TypeLM::LM_2NDORDER:
      y.Data()[0].Y().type = LineShape::TemperatureModel::T4;
      y.Data()[0].Y().X0 = x[0];
      y.Data()[0].Y().X1 = x[1];
      y.Data()[0].Y().X2 = x[7];
      y.Data()[0].G().type = LineShape::TemperatureModel::T4;
      y.Data()[0].G().X0 = x[2];
      y.Data()[0].G().X1 = x[3];
      y.Data()[0].G().X2 = x[8];
      y.Data()[0].DV().type = LineShape::TemperatureModel::T4;
      y.Data()[0].DV().X0 = x[4];
      y.Data()[0].DV().X1 = x[5];
      y.Data()[0].DV().X2 = x[9];
      break;
    case TypeLM::LM_1STORDER:
      y.Data()[0].Y().type = LineShape::TemperatureModel::T1;
      y.Data()[0].Y().X0 = x[1];
      y.Data()[0].Y().X1 = x[2];
      break;
    case TypeLM::LM_BYBAND:
      break;
  }
  return y;
}

void LineShape::Model::Set(const LineShape::ModelParameters& param,
                           const String& spec,
                           const LineShape::Variable var) {
  bool self = spec == self_broadening;
  bool bath = spec == bath_broadening;
  if (mself and self)
    mdata.front().Set(var, param);
  else if (self)
    throw std::runtime_error(
        "No self species but trying to set self in line shape model");
  else if (mbath and bath) {
    mdata.back().Set(var, param);
  } else if (bath)
    throw std::runtime_error(
        "No bath species but trying to set bath in line shape model");
  else {
    const SpeciesTag sp(spec);
    bool found = false;
    for (Index i = Index(mself); i < nelem() - Index(mbath); i++) {
      if (sp.Species() == mspecies[i].Species()) {
        found = true;
        mdata[i].Set(var, param);
      }
    }
    if (not found) {
      std::ostringstream os;
      os << "No species of type " << spec << " found in line shape model\n";
      os << "Available species are: " << mspecies << "\n";
      throw std::runtime_error(os.str());
    }
  }
}

Vector LineShape::Model::vmrs(const ConstVectorView& atmospheric_vmrs,
                              const ArrayOfArrayOfSpeciesTag& atmospheric_species,
                              const QuantumIdentifier& self) const {
  if (atmospheric_species.nelem() != atmospheric_vmrs.nelem())
    throw std::runtime_error("Bad atmospheric inputs");
  
  // Initialize list of VMRS to 0
  Vector line_vmrs(mspecies.nelem(), 0);
  const Index back = mspecies.nelem() - 1;  // Last index
  
  if (mtype == Type::DP) return line_vmrs;
  
  // Loop species ignoring self and bath
  for (Index i = 0; i < mspecies.nelem(); i++) {
    if (mbath and i == back) {
    } else {
      // Select target in-case this is self-broadening
      const auto target =
      (mself and i == 0) ? self.Species() : mspecies[i].Species();
      
      // Find species in list or do nothing at all
      Index this_species_index = -1;
      for (Index j = 0; j < atmospheric_species.nelem(); j++)
        if (atmospheric_species[j][0].Species() == target)
          this_species_index = j;
        
      // Set to non-zero in-case species exists
      if (this_species_index not_eq -1)
        line_vmrs[i] = atmospheric_vmrs[this_species_index];
    }
  }
  
  // Renormalize, if bath-species exist this is automatic.
  if (mbath)
    line_vmrs[back] = 1.0 - line_vmrs.sum();
  else if(line_vmrs.sum() == 0)  // Special case, there should be no atmosphere if this happens???
    return line_vmrs;
  else
    line_vmrs /= line_vmrs.sum();
  
  // The result must be non-zero, a real number, and finite
  if (not std::isnormal(line_vmrs.sum()))
    throw std::runtime_error(
      "Bad VMRs, your atmosphere does not support the line of interest");
    
  return line_vmrs;
}
