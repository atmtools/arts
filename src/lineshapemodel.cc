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

ArrayOfString AllLineShapeCoeffs() { return {"X0", "X1", "X2", "X3"}; }

ArrayOfString AllLineShapeVars() {
  return {"G0", "D0", "G2", "D2", "FVC", "ETA", "Y", "G", "DV"};
}

Jacobian::Line select_derivativeLineShape(const String& var,
                                          const String& coeff) {
  // Test viability of model variables
  static const ArrayOfString vars = AllLineShapeVars();
  bool var_OK = false;
  for (auto& v : vars)
    if (var == v) var_OK = true;

  // Test viability of model coefficients
  static const ArrayOfString coeffs = AllLineShapeCoeffs();
  bool coeff_OK = false;
  for (auto& c : coeffs)
    if (coeff == c) coeff_OK = true;

  // Fails either when the user has bad input or when the developer fails to update AllLineShapeVars or AllLineShapeCoeffs
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
      return Jacobian::Line::Shape##ID##X0; \
    else if (coeff == "X1")                     \
      return Jacobian::Line::Shape##ID##X1; \
    else if (coeff == "X2")                     \
      return Jacobian::Line::Shape##ID##X2; \
    else if (coeff == "X3")                     \
      return Jacobian::Line::Shape##ID##X3; \
  }

  if ReturnJacPropMatType(G0)
  else if ReturnJacPropMatType(D0)
  else if ReturnJacPropMatType(G2)
  else if ReturnJacPropMatType(D2)
  else if ReturnJacPropMatType(FVC)
  else if ReturnJacPropMatType(ETA)
  else if ReturnJacPropMatType(Y)
  else if ReturnJacPropMatType(G)
  else if ReturnJacPropMatType(DV)
#undef ReturnJacPropMatType
  
  return Jacobian::Line::FINAL;
}

std::istream& LineShape::from_artscat4(std::istream& is,
                                       Type& mtype,
                                       bool& self,
                                       bool& bath,
                                       Model& m,
                                       ArrayOfSpeciesTag& species,
                                       const QuantumIdentifier& qid) {
  // Set or reset variables
  mtype = Type::VP;
  self = true;
  bath = false;
  m.mdata = std::vector<SingleSpeciesModel>(7);
  species = ArrayOfSpeciesTag(7);

  // Set species
  species[1] = SpeciesTag("N2");
  species[2] = SpeciesTag("O2");
  species[3] = SpeciesTag("H2O");
  species[4] = SpeciesTag("CO2");
  species[5] = SpeciesTag("H2");
  species[6] = SpeciesTag("He");

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
  m.mdata.front().D0().X0 = 0;
  for (int k = 1; k < 7; k++)
    is >> m.mdata[k].D0().X0;

  // Special case when self is part of this list, it needs to be removed
  for (int k = 1; k < 7; k++) {
    if (qid.Species() == species[k].Species()) {
      if(m.mdata.front().G0().X0 not_eq m.mdata[k].G0().X0 or
        m.mdata.front().G0().X1 not_eq m.mdata[k].G0().X1 or
        m.mdata.front().D0().X1 not_eq m.mdata[k].D0().X1) {
        std::ostringstream os;
        os << "Species is " << qid.SpeciesName() << " and this is a broadening species in ARTSCAT-4.\n"
           << "Despite this, values representing self and " << qid.SpeciesName() << " does not match "
           << "in input string\n";
        throw std::runtime_error(os.str());
      }
      m.mdata.front().D0().X0 = m.mdata[k].D0().X0;
      m.Remove(k, species);
      break;
    }
  }
  
  return is;
}

std::istream& LineShape::from_linefunctiondata(std::istream& data,
                                               Type& mtype,
                                               bool& self,
                                               bool& bath,
                                               Model& m,
                                               ArrayOfSpeciesTag& species) {
  self = bath = false;
  Index specs;
  String s;

  // The first tag should give the line shape scheme
  data >> mtype;
  check_enum_error(mtype, "Bad Data");

  // Order of elements for line shape
  const auto shapeparams =
      LegacyLineFunctionData::lineshapetag2variablesvector(toString(mtype));

  // The second tag should give the line mixing scheme
  data >> s;

  // Order of elements for line mixing
  const auto mixingparams =
      LegacyLineFunctionData::linemixingtag2variablesvector(s);

  // The third tag should contain the number of species
  data >> specs;
  species.resize(specs);
  m.mdata.resize(specs);

  if (not specs and mtype not_eq Type::DP)
    throw std::runtime_error(
        "Need at least one species for non-Doppler line shapes");

  // For all species, we need to set the methods to compute them
  for (Index i = 0; i < specs; i++) {
    // This should be a species tag or one of the specials, SELF or BATH
    data >> s;
    if (s == self_broadening) {
      // If the species is self, then  we need to flag this
      self = true;
      if (i not_eq 0)  // but self has to be first for consistent behavior
        throw std::runtime_error("Self broadening must be first, it is not\n");
    }

    else if (s == bath_broadening) {
      // If the species is air, then we need to flag this
      bath = true;
      if (i not_eq
          specs - 1)  // but air has to be last because it needs the rest's VMR
        throw std::runtime_error(
            "Air/bath broadening must be last, it is not\n");
    } else {
      // Otherwise, we hope we find a species
      try {
        species[i] = SpeciesTag(s);
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

        const auto type = toTemperatureModel(s);
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
          if (ntemp > 12) {
            throw std::runtime_error(
                "Too many input parameters in interpolation results Legacy mode.");
          }
          Numeric temp;
          data >> temp;  // should be 200
          data >> temp;  // should be 250
          data >> temp;  // should be 296
          data >> temp;  // should be 340
          data >> m.mdata[i].Y().X0 >> m.mdata[i].Y().X1 >> m.mdata[i].Y().X2 >> m.mdata[i].Y().X3 
               >> m.mdata[i].G().X0 >> m.mdata[i].G().X1 >> m.mdata[i].G().X2 >> m.mdata[i].G().X3;
        }
      }
    }
  }

  return data;
}

std::istream& LineShape::from_pressurebroadeningdata(
  std::istream& data,
  LineShape::Type& mtype,
  bool& self,
  bool& bath,
  Model& m,
  ArrayOfSpeciesTag& species,
  const QuantumIdentifier& qid) {
  String s;
  data >> s;

  const auto type = LegacyPressureBroadeningData::string2typepb(s);
  const auto n = LegacyPressureBroadeningData::typepb2nelem(type);
  const auto self_in_list =
      LegacyPressureBroadeningData::self_listed(qid, type);

  Vector x(n);
  for (auto& num : x) data >> num;

  LegacyPressureBroadeningData::vector2modelpb(mtype,
                                               self,
                                               bath,
                                               m,
                                               species,
                                               x,
                                               type,
                                               self_in_list);

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
void LineShape::LegacyPressureBroadeningData::vector2modelpb(
  LineShape::Type& mtype,
  bool& self,
  bool& bath,
  Model& m,
  ArrayOfSpeciesTag& species,
  Vector x,
  LineShape::LegacyPressureBroadeningData::TypePB type,
  bool self_in_list) {
  switch (type) {
    case TypePB::PB_NONE:
      mtype = LineShape::Type::DP;
      self = bath = false;
      m = Model();
      species.resize(0);
      return;
    case TypePB::PB_AIR_BROADENING:
      mtype = LineShape::Type::VP;
      self = bath = true;
      m = Model(x[0], x[1], x[2], x[3], x[4]);
      species.resize(2);
      return;
    case TypePB::PB_AIR_AND_WATER_BROADENING:
      if (self_in_list) {
        mtype = LineShape::Type::VP;
        self = false;
        bath = true;
        m.Data().resize(2);
        m.Data()[0].G0() = {TemperatureModel::T1, x[0], x[1], 0, 0};
        m.Data()[0].D0() = {TemperatureModel::T5, x[2], x[1], 0, 0};
        m.Data()[1].G0() = {TemperatureModel::T1, x[3], x[4], 0, 0};
        m.Data()[1].D0() = {TemperatureModel::T5, x[5], x[4], 0, 0};
        species.resize(2);
        species[0] = SpeciesTag("H2O");
        return;
      } else {
        mtype = LineShape::Type::VP;
        self = bath = true;
        m.Data().resize(2);
        m.Data()[0].G0() = {TemperatureModel::T1, x[0], x[1], 0, 0};
        m.Data()[0].D0() = {TemperatureModel::T5, x[2], x[1], 0, 0};
        m.Data()[2].G0() = {TemperatureModel::T1, x[3], x[4], 0, 0};
        m.Data()[2].D0() = {TemperatureModel::T5, x[5], x[4], 0, 0};
        m.Data()[1].G0() = {TemperatureModel::T1, x[6], x[7], 0, 0};
        m.Data()[1].D0() = {TemperatureModel::T5, x[8], x[7], 0, 0};
        species.resize(3);
        species[1] = SpeciesTag("H2O");
        return;
      }
    case TypePB::PB_PLANETARY_BROADENING:
      if (self_in_list) {
        mtype = LineShape::Type::VP;
        self = bath = false;
        m.Data().resize(6);
        m.Data()[0].G0() = {TemperatureModel::T1, x[1], x[8], 0, 0};
        m.Data()[0].D0() = {TemperatureModel::T5, x[14], x[8], 0, 0};
        m.Data()[1].G0() = {TemperatureModel::T1, x[2], x[9], 0, 0};
        m.Data()[1].D0() = {TemperatureModel::T5, x[15], x[9], 0, 0};
        m.Data()[2].G0() = {TemperatureModel::T1, x[3], x[10], 0, 0};
        m.Data()[2].D0() = {TemperatureModel::T5, x[16], x[10], 0, 0};
        m.Data()[3].G0() = {TemperatureModel::T1, x[4], x[11], 0, 0};
        m.Data()[3].D0() = {TemperatureModel::T5, x[17], x[11], 0, 0};
        m.Data()[4].G0() = {TemperatureModel::T1, x[5], x[12], 0, 0};
        m.Data()[4].D0() = {TemperatureModel::T5, x[18], x[12], 0, 0};
        m.Data()[5].G0() = {TemperatureModel::T1, x[6], x[13], 0, 0};
        m.Data()[5].D0() = {TemperatureModel::T5, x[19], x[13], 0, 0};
        species = {SpeciesTag(String("N2")),
                   SpeciesTag(String("O2")),
                   SpeciesTag(String("H2O")),
                   SpeciesTag(String("CO2")),
                   SpeciesTag(String("H2")),
                   SpeciesTag(String("He"))};
        return;
      } else {
        mtype = LineShape::Type::VP;
        self = true;
        bath = false;
        m.Data().resize(7);
        m.Data()[0].G0() = {TemperatureModel::T1, x[0], x[7], 0, 0};
        //          ssm[0].D0() = ...
        m.Data()[1].G0() = {TemperatureModel::T1, x[1], x[8], 0, 0};
        m.Data()[1].D0() = {TemperatureModel::T5, x[14], x[8], 0, 0};
        m.Data()[2].G0() = {TemperatureModel::T1, x[2], x[9], 0, 0};
        m.Data()[2].D0() = {TemperatureModel::T5, x[15], x[9], 0, 0};
        m.Data()[3].G0() = {TemperatureModel::T1, x[3], x[10], 0, 0};
        m.Data()[3].D0() = {TemperatureModel::T5, x[16], x[10], 0, 0};
        m.Data()[4].G0() = {TemperatureModel::T1, x[4], x[11], 0, 0};
        m.Data()[4].D0() = {TemperatureModel::T5, x[17], x[11], 0, 0};
        m.Data()[5].G0() = {TemperatureModel::T1, x[5], x[12], 0, 0};
        m.Data()[5].D0() = {TemperatureModel::T5, x[18], x[12], 0, 0};
        m.Data()[6].G0() = {TemperatureModel::T1, x[6], x[13], 0, 0};
        m.Data()[6].D0() = {TemperatureModel::T5, x[19], x[13], 0, 0};
        species.resize(7);
        species[1] = SpeciesTag(String("N2"));
        species[2] = SpeciesTag(String("O2"));
        species[3] = SpeciesTag(String("H2O"));
        species[4] = SpeciesTag(String("CO2"));
        species[5] = SpeciesTag(String("H2"));
        species[6] = SpeciesTag(String("He"));
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
      y.Data().front().Y().type = LineShape::TemperatureModel::LM_AER;
      y.Data().front().G().type = LineShape::TemperatureModel::LM_AER;
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
      y.Data().front().G().type = LineShape::TemperatureModel::T0;
      y.Data().front().G().X0 = x[0];
      break;
    case TypeLM::LM_2NDORDER:
      y.Data().front().Y().type = LineShape::TemperatureModel::T4;
      y.Data().front().Y().X0 = x[0];
      y.Data().front().Y().X1 = x[1];
      y.Data().front().Y().X2 = x[7];
      y.Data().front().G().type = LineShape::TemperatureModel::T4;
      y.Data().front().G().X0 = x[2];
      y.Data().front().G().X1 = x[3];
      y.Data().front().G().X2 = x[8];
      y.Data().front().DV().type = LineShape::TemperatureModel::T4;
      y.Data().front().DV().X0 = x[4];
      y.Data().front().DV().X1 = x[5];
      y.Data().front().DV().X2 = x[9];
      break;
    case TypeLM::LM_1STORDER:
      y.Data().front().Y().type = LineShape::TemperatureModel::T1;
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
                       const QuantumIdentifier& self,
                       const ArrayOfSpeciesTag& lineshape_species,
                       bool self_in_list,
                       bool bath_in_list,
                       Type type) {
  if (atmospheric_species.nelem() != atmospheric_vmrs.nelem())
    throw std::runtime_error("Bad atmospheric inputs");
  
  // Initialize list of VMRS to 0
  Vector line_vmrs(lineshape_species.nelem(), 0);
  const Index back = lineshape_species.nelem() - 1;  // Last index
  
  if (type == Type::DP) return line_vmrs;
  
  // Loop species ignoring self and bath
  for (Index i = 0; i < lineshape_species.nelem()-bath_in_list; i++) {
    // Select target in-case this is self-broadening
    const auto target =
    (self_in_list and  &lineshape_species[i] == &lineshape_species.front()) ? self.Species() : lineshape_species[i].Species();
    
    // Find species in list or do nothing at all
    Index this_species_index = -1;
    for (Index j = 0; j < atmospheric_species.nelem(); j++)
      if (atmospheric_species[j][0].Species() == target)
        this_species_index = j;
      
    // Set to non-zero in-case species exists
    if (this_species_index not_eq -1)
      line_vmrs[i] = atmospheric_vmrs[this_species_index];
  }
  
  // Renormalize, if bath-species exist this is automatic.
  if (bath_in_list)
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

Vector LineShape::mass(const ConstVectorView& atmospheric_vmrs,
                       const ArrayOfArrayOfSpeciesTag& atmospheric_species,
                       const QuantumIdentifier& self,
                       const ArrayOfSpeciesTag& lineshape_species,
                       bool self_in_list,
                       bool bath_in_list) {
  if (atmospheric_species.nelem() != atmospheric_vmrs.nelem())
    throw std::runtime_error("Bad atmospheric inputs");
  
  // Initialize list of VMRS to 0
  Vector line_vmrs(lineshape_species.nelem(), 0);
  Vector line_mass(lineshape_species.nelem(), 0);
  const Index back = lineshape_species.nelem() - 1;  // Last index
  
  // Loop species ignoring self and bath
  for (Index i = 0; i < lineshape_species.nelem()-bath_in_list; i++) {
    // Select target in-case this is self-broadening
    const auto target =
    (self_in_list and  &lineshape_species[i] == &lineshape_species.front()) ? self.Species() : lineshape_species[i].Species();
    
    // Find species in list or do nothing at all
    Index this_species_index = -1;
    for (Index j = 0; j < atmospheric_species.nelem(); j++)
      if (atmospheric_species[j][0].Species() == target)
        this_species_index = j;
      
    // Set to non-zero in-case species exists
    if (this_species_index not_eq -1) {
      line_vmrs[i] = atmospheric_vmrs[this_species_index];
      line_mass[i] = atmospheric_species[this_species_index][0].SpeciesMass();
    }
  }
  
  // Renormalize, if bath-species exist this is automatic.
  if(line_vmrs.sum() == 0) {  // Special case, there should be no atmosphere if this happens???
    return line_mass;
  } else if (bath_in_list) {
    line_mass[back] = (line_vmrs * line_mass) / line_vmrs.sum();
  }
  
  if (self_in_list) {
    line_mass[0] = self.SpeciesMass();
  }
  
  // The result must be non-zero, a real number, and finite
  if (not std::isnormal(line_mass.sum())) {
    throw std::runtime_error(
      "Bad VMRs, your atmosphere does not support the line of interest");
  }
    
  return line_vmrs;
}

std::ostream& LineShape::operator<<(std::ostream& os, const LineShape::Model& m)
{
  for(auto& data: m.Data())
    os << data;
  return os;
}

std::istream& LineShape::operator>>(std::istream& is, Model& m)
{
  for(auto& data: m.Data())
    is >> data;
  return is;
}


String LineShape::ModelShape2MetaData(const Model& m)
{
  String out = "";
  
  const auto names = AllLineShapeVars();
  std::vector<Variable> vars(0);
  for (auto& n: names)
    vars.push_back(toVariable(n));
  
  for (auto& var: vars) {
    if (std::any_of(m.Data().cbegin(), m.Data().cend(),
      [var](auto& x){return x.Get(var).type not_eq TemperatureModel::None;})) {
      out += String(toString(var)) + ' ';
      for (auto& ssm: m.Data())
        out += String(toString(ssm.Get(var).type)) + ' ';
    }
  }
  
  if(out.size())
    out.pop_back();
  
  return out;
}


LineShape::Model LineShape::MetaData2ModelShape(const String& s)
{
  if (s.nelem() == 0)
    return LineShape::Model();
  
  const auto names = AllLineShapeVars();
  
  std::istringstream str(s);
  String part;
  Variable var=Variable::ETA;
  TemperatureModel tm=TemperatureModel::None;
  Index i=-100000;
  
  std::vector<SingleSpeciesModel> ssms(0);
  while (not str.eof()) {
    str >> part;
    if(std::any_of(names.cbegin(), names.cend(),
      [part](auto x){return part == x;})) {
      i=-1;
      var = toVariable(part);
    }
    else {
      i++;
      tm = toTemperatureModel(part);
    }
    
    if (i < 0)
      continue;
    else if (i < Index(ssms.size()))
      goto add_var;
    else {
      ssms.push_back(SingleSpeciesModel());
      add_var:
      auto mp = ssms[i].Get(var);
      mp.type = tm;
      ssms[i].Set(var, mp);
    }
  }
  
  return Model(std::move(ssms));
}

String LineShape::modelparameters2metadata(const LineShape::ModelParameters mp, const Numeric T0)
{
  std::ostringstream os;
  switch (mp.type) {
    case TemperatureModel::None:
      os << 0;
      break;
    case TemperatureModel::T0:
      os << mp.X0;
      break;
    case TemperatureModel::T1:
      os << mp.X0 << " * (" << T0 << "/T)^" << mp.X1;
      break;
    case TemperatureModel::T2:
      os << mp.X0 << " * (" << T0 << "/T)^" << mp.X1 << " / (1 + " << mp.X2 << " * log(T/" << T0 << "))";
      break;
    case TemperatureModel::T3:
      os << mp.X0 << " + " << mp.X1 << " * (" << T0 << " - T)";
      break;
    case TemperatureModel::T4:
      os << "(" << mp.X0 << " + " << mp.X1 << " * (" << T0 << "/T - 1)) * (" << T0 << "/T)^" << mp.X2;
      break;
    case TemperatureModel::T5:
      os << mp.X0 << " * (" << T0 << "/T)^(0.25 + 1.5 * " << mp.X1 << ")";
      break;
    case TemperatureModel::LM_AER:
      os << '(' << "Linear interpolation to y(x) from x-ref = [200, 250, 296, 340] and y-ref = [" << mp.X0 << ", " << mp.X1 << ", " << mp.X2 << ", " << mp.X3 << ']' << ')';
      break;
    case TemperatureModel::DPL:
      os << '(' << mp.X0 << " * (" << T0 << "/T)^" << mp.X1 << " + "  << mp.X2 << " * (" << T0 << "/T)^" << mp.X3 << ')';
      break;
    case TemperatureModel::FINAL: break;
  }
  
  return os.str();
}

ArrayOfString LineShape::ModelMetaDataArray(const LineShape::Model& m,
                                const bool self,
                                const bool bath,
                                const ArrayOfSpeciesTag& sts,
                                const Numeric T0)
{
  const auto names = AllLineShapeVars();
  std::vector<Variable> vars(0);
  for (auto& n: names)
    vars.push_back(toVariable(n));
  
  ArrayOfString as(0);
  
  for (Index i=0; i<names.nelem(); i++) {
    Variable var = vars[i];
    
    if (std::any_of(m.Data().cbegin(), m.Data().cend(),
      [var](auto& x){return x.Get(var).type not_eq TemperatureModel::None;})) {
      
      std::ostringstream os;
      os << names[i] << " ~ ";
      for (Index j=0; j<sts.nelem(); j++) {
        if (&sts[j] == &sts.front() and self)
          os << "VMR(" << self_broadening << ") * "
             << modelparameters2metadata(m.Data().front().Get(var), T0);
        else if (&sts[j] == &sts.back() and bath)
          os << "VMR(" << bath_broadening << ") * "
             << modelparameters2metadata(m.Data().back().Get(var), T0);
        else 
          os << "VMR(" << sts[j].SpeciesNameMain() << ") * "
             << modelparameters2metadata(m.Data()[j].Get(var), T0);
             
        if (&sts[j] not_eq &sts.back())
          os << " + ";
      }
      as.push_back(os.str());
    }
  }
  
  return as;
}


namespace LineShape {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
Numeric& SingleModelParameter(ModelParameters& mp, const String& type) {
  if (type == "X0")
    return mp.X0;
  else if (type == "X1")
    return mp.X1;
  else if (type == "X2")
    return mp.X2;
  else if (type == "X3")
    return mp.X3;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
    << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}
#pragma GCC diagnostic pop

std::ostream& operator<<(std::ostream& os, const ModelParameters& mp) {
  os << mp.type << ' ' << mp.X0 << ' ' << mp.X1 << ' '
  << mp.X2 << ' ' << mp.X3 << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, ModelParameters& mp) {
  is >> mp.type >> mp.X0 >> mp.X1 >> mp.X2 >> mp.X3;
  return is;
}

#define x0 X[Index(var)].X0
#define x1 X[Index(var)].X1
#define x2 X[Index(var)].X2
#define x3 X[Index(var)].X3

Numeric SingleSpeciesModel::compute(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = x0; break;
    case TemperatureModel::T1:
      out = x0 * pow(T0 / T, x1); break;
    case TemperatureModel::T2:
      out = x0 * pow(T0 / T, x1) * (1 + x2 * log(T / T0)); break;
    case TemperatureModel::T3:
      out = x0 + x1 * (T - T0); break;
    case TemperatureModel::T4:
      out = (x0 + x1 * (T0 / T - 1.)) * pow(T0 / T, x2); break;
    case TemperatureModel::T5:
      out = x0 * pow(T0 / T, 0.25 + 1.5 * x1); break;
    case TemperatureModel::LM_AER:
      out = special_linemixing_aer(T, X[Index(var)]); break;
    case TemperatureModel::DPL:
      out = x0 * pow(T0 / T, x1) + x2 * pow(T0 / T, x3); break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

Numeric SingleSpeciesModel::compute_dX0(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = 1; break;
    case TemperatureModel::T1:
      out = pow(T0 / T, x1); break;
    case TemperatureModel::T2:
      out = pow(T0 / T, x1) * (1 + x2 * log(T / T0)); break;
    case TemperatureModel::T3:
      out = 1; break;
    case TemperatureModel::T4:
      out = pow(T0 / T, x2); break;
    case TemperatureModel::T5:
      out = pow(T0 / T, 1.5 * x1 + 0.25); break;
    case TemperatureModel::LM_AER:
      out = 0; break;
    case TemperatureModel::DPL:
      out = pow(T0 / T, x1); break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

Numeric SingleSpeciesModel::compute_dX1(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = 0; break;
    case TemperatureModel::T1:
      out = x0 * pow(T0 / T, x1) * log(T0 / T); break;
    case TemperatureModel::T2:
      out = x0 * pow(T0 / T, x1) * (x2 * log(T / T0) + 1.) * log(T0 / T); break;
    case TemperatureModel::T3:
      out = (T - T0); break;
    case TemperatureModel::T4:
      out = pow(T0 / T, x2) * (T0 / T - 1.); break;
    case TemperatureModel::T5:
      out = 1.5 * x0 * pow(T0 / T, 1.5 * x1 + 0.25) * log(T0 / T); break;
    case TemperatureModel::LM_AER:
      out = 0; break;
    case TemperatureModel::DPL:
      out = x0 * pow(T0 / T, x1) * log(T0 / T); break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

Numeric SingleSpeciesModel::compute_dX2(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = 0; break;
    case TemperatureModel::T1:
      out = 0; break;
    case TemperatureModel::T2:
      out = x0 * pow(T0 / T, x1) * log(T / T0); break;
    case TemperatureModel::T3:
      out = 0; break;
    case TemperatureModel::T4:
      out = pow(T0 / T, x2) * (x0 + x1 * (T0 / T - 1)) * log(T0 / T); break;
    case TemperatureModel::T5:
      out = 0; break;
    case TemperatureModel::LM_AER:
      out = 0; break;
    case TemperatureModel::DPL:
      out = pow(T0 / T, x3); break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

Numeric SingleSpeciesModel::compute_dX3(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = 0; break;
    case TemperatureModel::T1:
      out = 0; break;
    case TemperatureModel::T2:
      out = 0; break;
    case TemperatureModel::T3:
      out = 0; break;
    case TemperatureModel::T4:
      out = 0; break;
    case TemperatureModel::T5:
      out = 0; break;
    case TemperatureModel::LM_AER:
      out = 0; break;
    case TemperatureModel::DPL:
      out = x2 * pow(T0 / T, x3) * log(T0 / T); break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

Numeric SingleSpeciesModel::compute_dT(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = 0; break;
    case TemperatureModel::T1:
      out = -x0 * x1 * pow(T0 / T, x1) / T; break;
    case TemperatureModel::T2:
      out = -x0 * x1 * pow(T0 / T, x1) * (x2 * log(T / T0) + 1.) / T +
      x0 * x2 * pow(T0 / T, x1) / T; break;
    case TemperatureModel::T3:
      out = x1; break;
    case TemperatureModel::T4:
      out = -x2 * pow(T0 / T, x2) * (x0 + x1 * (T0 / T - 1.)) / T -
      T0 * x1 * pow(T0 / T, x2) / pow(T, 2); break;
    case TemperatureModel::T5:
      out = -x0 * pow(T0 / T, 1.5 * x1 + 0.25) * (1.5 * x1 + 0.25) / T; break;
    case TemperatureModel::LM_AER:
      out = special_linemixing_aer_dT(T, X[Index(var)]); break;
    case TemperatureModel::DPL:
      out = -x0 * x1 * pow(T0 / T, x1) / T + -x2 * x3 * pow(T0 / T, x3) / T; break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

Numeric SingleSpeciesModel::compute_dT0(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  Numeric out=std::numeric_limits<Numeric>::quiet_NaN();
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      out = 0; break;
    case TemperatureModel::T0:
      out = 0; break;
    case TemperatureModel::T1:
      out = x0 * x1 * pow(T0 / T, x1) / T0; break;
    case TemperatureModel::T2:
      out = x0 * x1 * pow(T0 / T, x1) * (x2 * log(T / T0) + 1.) / T0 -
      x0 * x2 * pow(T0 / T, x1) / T0; break;
    case TemperatureModel::T3:
      out = -x1; break;
    case TemperatureModel::T4:
      out = x2 * pow(T0 / T, x2) * (x0 + x1 * (T0 / T - 1.)) / T0 +
      x1 * pow(T0 / T, x2) / T; break;
    case TemperatureModel::T5:
      out = x0 * pow(T0 / T, 1.5 * x1 + 0.25) * (1.5 * x1 + 0.25) / T0; break;
    case TemperatureModel::LM_AER:
      out = 0; break;
    case TemperatureModel::DPL:
      out = x0 * x1 * pow(T0 / T, x1) / T0 + x2 * x3 * pow(T0 / T, x3) / T0; break;
    case TemperatureModel::FINAL: break;
  }
  return out;
}

#undef x0
#undef x1
#undef x2
#undef x3

bifstream & SingleSpeciesModel::read(bifstream& bif) {
  for (auto& data: X) {
    Index x;
    bif >> x >> data.X0 >> data.X1 >> data.X2 >> data.X3;
    data.type = TemperatureModel(x);
  }
  return bif;
}

bofstream & SingleSpeciesModel::write(bofstream& bof) const {
  for (auto& data: X) {
    bof << Index(data.type) << data.X0 << data.X1 << data.X2 << data.X3;
  }
  return bof;
}

bool SingleSpeciesModel::MatchTypes(const SingleSpeciesModel& other) const noexcept {
  return std::equal (X.cbegin(), X.cend(), other.X.cbegin(), other.X.cend(), [](const auto& a, const auto& b){return a.type == b.type;});
}

std::ostream& operator<<(std::ostream& os, const SingleSpeciesModel& ssm) {
  for (const auto& mp : ssm.Data())
    if (mp.type not_eq TemperatureModel::None)
      os << mp.X0 << ' ' << mp.X1 << ' ' << mp.X2 << ' ' << mp.X3 << ' ';
  return os;
}

std::istream& operator>>(std::istream& is, SingleSpeciesModel& ssm) {
  for (auto& mp : ssm.Data())
    if(mp.type not_eq TemperatureModel::None)
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
             std::array<Numeric, 12> aer_interp) noexcept : mdata(2) {
  mdata.front().G0() = {TemperatureModel::T1, sgam, nself, 0, 0};
  mdata.front().D0() = {TemperatureModel::T5, psf, nair, 0, 0};

  mdata.back().G0() = {TemperatureModel::T1, agam, nair, 0, 0};
  mdata.back().D0() = {TemperatureModel::T5, psf, nair, 0, 0};
  
  if (std::any_of(aer_interp.cbegin(), aer_interp.cend(), [](auto x){return x not_eq 0;})) {
    mdata.front().Y().type = TemperatureModel::LM_AER;
    mdata.front().Y().X0 = aer_interp[4];
    mdata.front().Y().X1 = aer_interp[5];
    mdata.front().Y().X2 = aer_interp[6];
    mdata.front().Y().X3 = aer_interp[7];
    mdata.front().G().type = TemperatureModel::LM_AER;
    mdata.front().G().X0 = aer_interp[8];
    mdata.front().G().X1 = aer_interp[9];
    mdata.front().G().X2 = aer_interp[10];
    mdata.front().G().X3 = aer_interp[11];
    
    mdata.back().Y().type = TemperatureModel::LM_AER;
    mdata.back().Y().X0 = aer_interp[4];
    mdata.back().Y().X1 = aer_interp[5];
    mdata.back().Y().X2 = aer_interp[6];
    mdata.back().Y().X3 = aer_interp[7];
    mdata.back().G().type = TemperatureModel::LM_AER;
    mdata.back().G().X0 = aer_interp[8];
    mdata.back().G().X1 = aer_interp[9];
    mdata.back().G().X2 = aer_interp[10];
    mdata.back().G().X3 = aer_interp[11];
  }
}

bool Model::OK(Type type, bool self, bool bath,
               const std::vector<SpeciesTag>& species) const noexcept {
  Index n = mdata.size();
  Index k = species.size();
  Index m = Index(self) + Index(bath);
  bool needs_any = type not_eq Type::DP;
  if (n not_eq k or m > n or (needs_any and not n))
    return false;
  else
    return true;
}

#define LSPC(XVAR, PVAR)                                                       \
  Numeric Model::XVAR(                                                         \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], ConstVectorView vmrs) \
      const noexcept {                                                         \
    return PVAR *                                                              \
           std::inner_product(                                                 \
               mdata.cbegin(),                                                 \
               mdata.cend(),                                                   \
               vmrs.begin(),                                                   \
               0.0,                                                            \
               std::plus<Numeric>(),                                           \
               [=](auto& x, auto vmr) -> Numeric {                             \
                 return vmr * x.compute(T, T0, Variable::XVAR);                \
               });                                                             \
  }
  LSPC(G0, P)
  LSPC(D0, P)
  LSPC(G2, P)
  LSPC(D2, P)
  LSPC(FVC, P)
  LSPC(ETA, 1)
  LSPC(Y, P)
  LSPC(G, P* P)
  LSPC(DV, P* P)
#undef LSPC

#define LSPCV(XVAR, PVAR)                                            \
  Numeric Model::d##XVAR##_dVMR(Numeric T,                           \
                         Numeric T0,                                 \
                         Numeric P [[maybe_unused]],                 \
                         const Index deriv_pos) const noexcept {     \
    if (deriv_pos not_eq -1)                                         \
      return PVAR * mdata[deriv_pos].compute(T, T0, Variable::XVAR); \
    else                                                             \
      return 0;                                                      \
  }
  LSPCV(G0, P)
  LSPCV(D0, P)
  LSPCV(G2, P)
  LSPCV(D2, P)
  LSPCV(FVC, P)
  LSPCV(ETA, 1)
  LSPCV(Y, P)
  LSPCV(G, P* P)
  LSPCV(DV, P* P)
#undef LSPCV

#define LSPCT(XVAR, PVAR)                                                       \
  Numeric Model::d##XVAR##_dT(                                                  \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], ConstVectorView vmrs)  \
      const noexcept {                                                          \
    return PVAR *                                                               \
           std::inner_product(                                                  \
               mdata.cbegin(),                                                  \
               mdata.cend(),                                                    \
               vmrs.begin(),                                                    \
               0.0,                                                             \
               std::plus<Numeric>(),                                            \
               [=](auto& x, auto vmr) -> Numeric {                              \
                 return vmr * x.compute_dT(T, T0, Variable::XVAR);              \
               });                                                              \
  }
  LSPCT(G0, P)
  LSPCT(D0, P)
  LSPCT(G2, P)
  LSPCT(D2, P)
  LSPCT(FVC, P)
  LSPCT(ETA, 1)
  LSPCT(Y, P)
  LSPCT(G, P* P)
  LSPCT(DV, P* P)
#undef LSPCT

#define LSPDC(XVAR, DERIV, PVAR)                                     \
  Numeric Model::d##XVAR##DERIV(Numeric T,                           \
                         Numeric T0,                                 \
                         Numeric P [[maybe_unused]],                 \
                         Index deriv_pos,                            \
                         ConstVectorView vmrs) const noexcept {      \
    if (deriv_pos not_eq -1)                                         \
      return vmrs[deriv_pos] * PVAR *                                \
             mdata[deriv_pos].compute##DERIV(T, T0, Variable::XVAR); \
    else                                                             \
      return 0;                                                      \
  }
  LSPDC(G0, _dT0, P)
  LSPDC(G0, _dX0, P)
  LSPDC(G0, _dX1, P)
  LSPDC(G0, _dX2, P)
  LSPDC(G0, _dX3, P)
  LSPDC(D0, _dT0, P)
  LSPDC(D0, _dX0, P)
  LSPDC(D0, _dX1, P)
  LSPDC(D0, _dX2, P)
  LSPDC(D0, _dX3, P)
  LSPDC(G2, _dT0, P)
  LSPDC(G2, _dX0, P)
  LSPDC(G2, _dX1, P)
  LSPDC(G2, _dX2, P)
  LSPDC(G2, _dX3, P)
  LSPDC(D2, _dT0, P)
  LSPDC(D2, _dX0, P)
  LSPDC(D2, _dX1, P)
  LSPDC(D2, _dX2, P)
  LSPDC(D2, _dX3, P)
  LSPDC(FVC, _dT0, P)
  LSPDC(FVC, _dX0, P)
  LSPDC(FVC, _dX1, P)
  LSPDC(FVC, _dX2, P)
  LSPDC(FVC, _dX3, P)
  LSPDC(ETA, _dT0, 1)
  LSPDC(ETA, _dX0, 1)
  LSPDC(ETA, _dX1, 1)
  LSPDC(ETA, _dX2, 1)
  LSPDC(ETA, _dX3, 1)
  LSPDC(Y, _dT0, P)
  LSPDC(Y, _dX0, P)
  LSPDC(Y, _dX1, P)
  LSPDC(Y, _dX2, P)
  LSPDC(Y, _dX3, P)
  LSPDC(G, _dT0, P* P)
  LSPDC(G, _dX0, P* P)
  LSPDC(G, _dX1, P* P)
  LSPDC(G, _dX2, P* P)
  LSPDC(G, _dX3, P* P)
  LSPDC(DV, _dT0, P* P)
  LSPDC(DV, _dX0, P* P)
  LSPDC(DV, _dX1, P* P)
  LSPDC(DV, _dX2, P* P)
  LSPDC(DV, _dX3, P* P)
#undef LSPDC

Output Model::GetParams(Numeric T,
                        Numeric T0,
                        Numeric P,
                        ConstVectorView vmrs) const noexcept {
  return {G0(T, T0, P, vmrs),
          D0(T, T0, P, vmrs),
          G2(T, T0, P, vmrs),
          D2(T, T0, P, vmrs),
          FVC(T, T0, P, vmrs),
          ETA(T, T0, P, vmrs),
          Y(T, T0, P, vmrs),
          G(T, T0, P, vmrs),
          DV(T, T0, P, vmrs)};
}

Output Model::GetParams(Numeric T,
                        Numeric T0,
                        Numeric P,
                        size_t k) const noexcept {
  return {P * mdata[k].compute(T, T0, Variable::G0),
          P * mdata[k].compute(T, T0, Variable::D0),
          P * mdata[k].compute(T, T0, Variable::G2),
          P * mdata[k].compute(T, T0, Variable::D2),
          P * mdata[k].compute(T, T0, Variable::FVC),
          mdata[k].compute(T, T0, Variable::ETA),
          P * mdata[k].compute(T, T0, Variable::Y),
          P * P * mdata[k].compute(T, T0, Variable::G),
          P * P * mdata[k].compute(T, T0, Variable::DV)};
}

Output Model::GetTemperatureDerivs(Numeric T,
                                   Numeric T0,
                                   Numeric P,
                                   ConstVectorView vmrs) const noexcept {
  return {dG0_dT(T, T0, P, vmrs),
          dD0_dT(T, T0, P, vmrs),
          dG2_dT(T, T0, P, vmrs),
          dD2_dT(T, T0, P, vmrs),
          dFVC_dT(T, T0, P, vmrs),
          dETA_dT(T, T0, P, vmrs),
          dY_dT(T, T0, P, vmrs),
          dG_dT(T, T0, P, vmrs),
          dDV_dT(T, T0, P, vmrs)};
}

Output Model::GetVMRDerivs(Numeric T, Numeric T0, Numeric P, const Index pos) const noexcept {
  return {dG0_dVMR(T, T0, P, pos),
    dD0_dVMR(T, T0, P, pos),
    dG2_dVMR(T, T0, P, pos),
    dD2_dVMR(T, T0, P, pos),
    dFVC_dVMR(T, T0, P, pos),
    dETA_dVMR(T, T0, P, pos),
    dY_dVMR(T, T0, P, pos),
    dG_dVMR(T, T0, P, pos),
    dDV_dVMR(T, T0, P, pos)};
}

Numeric Model::GetInternalDeriv(Numeric T,
                                Numeric T0,
                                Numeric P,
                                Index pos,
                                const Vector& vmrs,
                                Jacobian::Line deriv) const noexcept {
  if (pos < 0) return 0;

#define RETURNINTERNALDERIVATIVE(TYPE)         \
  case Jacobian::Line::Shape##TYPE##X0:        \
    return d##TYPE##_dX0(T, T0, P, pos, vmrs); \
  case Jacobian::Line::Shape##TYPE##X1:        \
    return d##TYPE##_dX1(T, T0, P, pos, vmrs); \
  case Jacobian::Line::Shape##TYPE##X2:        \
    return d##TYPE##_dX2(T, T0, P, pos, vmrs); \
  case Jacobian::Line::Shape##TYPE##X3:        \
    return d##TYPE##_dX3(T, T0, P, pos, vmrs)
  switch (deriv) {
    RETURNINTERNALDERIVATIVE(G0);
    RETURNINTERNALDERIVATIVE(D0);
    RETURNINTERNALDERIVATIVE(G2);
    RETURNINTERNALDERIVATIVE(D2);
    RETURNINTERNALDERIVATIVE(FVC);
    RETURNINTERNALDERIVATIVE(ETA);
    RETURNINTERNALDERIVATIVE(Y);
    RETURNINTERNALDERIVATIVE(G);
    RETURNINTERNALDERIVATIVE(DV);
    default:
      return std::numeric_limits<Numeric>::quiet_NaN();
  }
#undef RETURNINTERNALDERIVATIVE
}

void Model::Remove(Index i, ArrayOfSpeciesTag& specs) {
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

bool Model::Match(const Model& other) const noexcept {
  return std::equal (mdata.cbegin(), mdata.cend(), other.mdata.cbegin(), other.mdata.cend(), [](auto& a, auto& b) {return a.MatchTypes(b);});
}

bifstream& Model::read(bifstream& bif) {
  for (auto& data: mdata)
    data.read(bif);
  return bif;
}

bofstream& Model::write(bofstream& bof) const {
  for (auto& data: mdata)
    data.write(bof);
  return bof;
}

namespace LegacyLineFunctionData {
std::vector<Variable> lineshapetag2variablesvector(String type) {
  if (type == String("DP"))
    return {};
  else if (type == String("LP"))
    return {Variable::G0, Variable::D0};
  else if (type == String("VP"))
    return {Variable::G0, Variable::D0};
  else if (type == String("SDVP"))
    return {Variable::G0, Variable::D0, Variable::G2, Variable::D2};
  else if (type == String("HTP"))
    return {Variable::G0,
      Variable::D0,
      Variable::G2,
      Variable::D2,
      Variable::FVC,
      Variable::ETA};
      else {
        std::ostringstream os;
        os << "Type: " << type << ", is not accepted.  "
        << "See documentation for accepted types\n";
        throw std::runtime_error(os.str());
      }
}

std::vector<Variable> linemixingtag2variablesvector(String type) {
  if (type == "#")
    return {};
  else if (type == "LM1")
    return {Variable::Y};
  else if (type == "LM2")
    return {Variable::Y, Variable::G, Variable::DV};
  else if (type == "INT")
    return {};
  else if (type == "ConstG")
    return {Variable::G};
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
    << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}
}  // LegacyLineFunctionData

namespace LegacyLineMixingData {
LegacyLineMixingData::TypeLM string2typelm(String type) {
  if (type == "NA")  // The standard case
    return TypeLM::LM_NONE;
  else if (type == "LL")  // The LBLRTM case
    return TypeLM::LM_LBLRTM;
  else if (type == "NR")  // The LBLRTM O2 non-resonant case
    return TypeLM::LM_LBLRTM_O2NonResonant;
  else if (type == "L2")  // The 2nd order case
    return TypeLM::LM_2NDORDER;
  else if (type == "L1")  // The 2nd order case
    return TypeLM::LM_1STORDER;
  else if (type == "BB")  // The band class
    return TypeLM::LM_BYBAND;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
    << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}
}  // LegacyLineMixingData


namespace LegacyPressureBroadeningData {
LegacyPressureBroadeningData::TypePB string2typepb(String type) {
  if (type == "NA")  // The none case
    return TypePB::PB_NONE;
  else if (type == "N2")  // Air Broadening is N2 broadening mostly...
    return TypePB::PB_AIR_BROADENING;
  else if (type == "WA")  // Water and Air Broadening
    return TypePB::PB_AIR_AND_WATER_BROADENING;
  else if (type == "AP")  // Planetary broadening
    return TypePB::PB_PLANETARY_BROADENING;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
    << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}

Index self_listed(const QuantumIdentifier& qid,
                  LegacyPressureBroadeningData::TypePB t) {
  if (t == TypePB::PB_PLANETARY_BROADENING and
      (qid.Species() == SpeciesTag(String("N2")).Species() or
       qid.Species() == SpeciesTag(String("O2")).Species() or
       qid.Species() == SpeciesTag(String("H2O")).Species() or
       qid.Species() == SpeciesTag(String("CO2")).Species() or
       qid.Species() == SpeciesTag(String("H2")).Species() or
       qid.Species() == SpeciesTag(String("He")).Species()))
    return true;
  else if (t == TypePB::PB_AIR_AND_WATER_BROADENING and
           qid.Species() == SpeciesTag(String("H2O")).Species())
    return true;
  else
    return false;
}
}  // LegacyPressureBroadeningData
}  // LineShape
