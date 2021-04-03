#include "absorption.h"
#include "artstime.h"
#include "lineshape.h"
#include "gui/plot.h"

namespace {
  constexpr Index N = 10001;
  constexpr Index M = 25;
  constexpr Numeric dG0 = 10;
  constexpr Numeric dD0 = 1;
  constexpr Numeric dG2 = 1;
  constexpr Numeric dD2 = 1;
  constexpr Numeric dETA = 1e-4;
  constexpr Numeric dFVC = 1e-3;
  constexpr Numeric dY = 1e-10;
  constexpr Numeric dG = 1e-14;
  constexpr Numeric dDV = 1e-1;
  constexpr Numeric dVMR = 1e-4;
  constexpr Numeric dT = 0.1;
  constexpr Numeric dH = 1e-8;
  constexpr Numeric df = 1e1;
  constexpr Numeric dF0 = 1e1;
  constexpr Numeric dI0 = 1e1;
  static_assert (N > 1);
}

ENUMCLASS (ModificationInternal, char,
           None,
           T,
           f,
           H,
           F0,
           I0,
           SelfG0,
           SelfD0,
           SelfG2,
           SelfD2,
           SelfETA,
           SelfFVC,
           SelfY,
           SelfG,
           SelfDV,
           AirG0,
           AirD0,
           AirG2,
           AirD2,
           AirETA,
           AirFVC,
           AirY,
           AirG,
           AirDV,
           SelfVMR,
           AirVMR
          )

template <ModificationInternal mod>
LineShape::Model lineshape_model() {
  LineShape::Model model(2);
  model[0].G0() = LineShape::ModelParameters(LineShape::TemperatureModel::T1, 20e3, 0.8);
  model[0].D0() = LineShape::ModelParameters(LineShape::TemperatureModel::T5, 200, 0.8);
  model[0].G2() = LineShape::ModelParameters(LineShape::TemperatureModel::T3, 200, 20);
  model[0].D2() = LineShape::ModelParameters(LineShape::TemperatureModel::T0, 20);
  model[0].G2() = LineShape::ModelParameters(LineShape::TemperatureModel::T3, 200, 20);
  model[0].ETA() = LineShape::ModelParameters(LineShape::TemperatureModel::T0, 0.1);
  model[0].FVC() = LineShape::ModelParameters(LineShape::TemperatureModel::T2, 1e2, .2, 1);
  model[0].Y() = LineShape::ModelParameters(LineShape::TemperatureModel::LM_AER, 1e-6, 1e-7, 1e-8, 1e-9);
  model[0].G() = LineShape::ModelParameters(LineShape::TemperatureModel::LM_AER, 1e-12, 1e-11, 1e-10, 1e-9);
  model[0].DV() = LineShape::ModelParameters(LineShape::TemperatureModel::DPL, 1e2, 0.8, -1e1, 0.3);
  
  model[1].G0() = LineShape::ModelParameters(LineShape::TemperatureModel::T1, 1.2*20e3, 0.8);
  model[1].D0() = LineShape::ModelParameters(LineShape::TemperatureModel::T5, 1.2*200, 0.8);
  model[1].G2() = LineShape::ModelParameters(LineShape::TemperatureModel::T3, 1.2*200, 1.2*20);
  model[1].D2() = LineShape::ModelParameters(LineShape::TemperatureModel::T0, 1.2*20);
  model[1].G2() = LineShape::ModelParameters(LineShape::TemperatureModel::T3, 1.2*200, 1.2*20);
  model[1].ETA() = LineShape::ModelParameters(LineShape::TemperatureModel::T0, 1.2*0.1);
  model[1].FVC() = LineShape::ModelParameters(LineShape::TemperatureModel::T2, 1.2*1e2, 1.2*.2, 1.2*1);
  model[1].Y() = LineShape::ModelParameters(LineShape::TemperatureModel::LM_AER, 1.2*1e-6, 1.2*1e-7, 1.2*1e-8, 1.2*1e-9);
  model[1].G() = LineShape::ModelParameters(LineShape::TemperatureModel::LM_AER, 1.2*1e-12, 1.2*1e-11, 1.2*1e-10, 1.2*1e-9);
  model[1].DV() = LineShape::ModelParameters(LineShape::TemperatureModel::DPL, 1.2*1e2, 0.8, 1.2*-1e1, 0.3);
  
  // Modify
  if constexpr (mod == ModificationInternal::SelfG0)  model[0].G0().X0 += ::dG0;
  if constexpr (mod == ModificationInternal::SelfD0)  model[0].D0().X0 += ::dD0;
  if constexpr (mod == ModificationInternal::SelfG2)  model[0].G2().X0 += ::dG2;
  if constexpr (mod == ModificationInternal::SelfD2)  model[0].D2().X0 += ::dD2;
  if constexpr (mod == ModificationInternal::SelfETA) model[0].ETA().X0 += ::dETA;
  if constexpr (mod == ModificationInternal::SelfFVC) model[0].FVC().X0 += ::dFVC;
  if constexpr (mod == ModificationInternal::SelfY)   model[0].Y().X1 += ::dY;
  if constexpr (mod == ModificationInternal::SelfG)   model[0].G().X1 += ::dG;
  if constexpr (mod == ModificationInternal::SelfDV)  model[0].DV().X0 += ::dDV;
  if constexpr (mod == ModificationInternal::AirG0)   model[1].G0().X0 += ::dG0;
  if constexpr (mod == ModificationInternal::AirD0)   model[1].D0().X0 += ::dD0;
  if constexpr (mod == ModificationInternal::AirG2)   model[1].G2().X0 += ::dG2;
  if constexpr (mod == ModificationInternal::AirD2)   model[1].D2().X0 += ::dD2;
  if constexpr (mod == ModificationInternal::AirETA)  model[1].ETA().X0 += ::dETA;
  if constexpr (mod == ModificationInternal::AirFVC)  model[1].FVC().X0 += ::dFVC;
  if constexpr (mod == ModificationInternal::AirY)    model[1].Y().X1 += ::dY;
  if constexpr (mod == ModificationInternal::AirG)    model[1].G().X1 += ::dG;
  if constexpr (mod == ModificationInternal::AirDV)   model[1].DV().X0 += ::dDV;
    
  return model;
}

struct InternalData {
  QuantumIdentifier qid;
  AbsorptionLines band;
  Vector f_grid;
  Vector vmr;
  Numeric P;
  Numeric T;
  Numeric H;
};

template <ModificationInternal mod, LineShape::Type ls_type>
InternalData line_model() {
  QuantumNumbers upp, low;
  upp[QuantumNumberType::J] = 3;
  low[QuantumNumberType::J] = 2;
  QuantumIdentifier qid = QuantumIdentifier(0, 0, upp, low);
  
  Absorption::SingleLine line(1e9, 1e10, 1e-20, 1., 3., 1e-14,
                              Zeeman::Model(2, 1.5), lineshape_model<mod>(),
                              {low[QuantumNumberType::J]}, {upp[QuantumNumberType::J]});
  
  if constexpr (mod == ModificationInternal::F0) line.F0() += ::dF0;
  if constexpr (mod == ModificationInternal::I0) line.I0() += ::dI0;
  
  AbsorptionLines band(true, true,
                       Absorption::CutoffType::None, Absorption::MirroringType::None,
                       Absorption::PopulationType::LTE, Absorption::NormalizationType::None,
                       ls_type, 296, -1, -1, qid, {QuantumNumberType::J},
                       {SpeciesTag(), SpeciesTag()}, {line});
  
  const Numeric f0 = ((mod == ModificationInternal::f) ? 500e6 + ::df : 500e6) + (ls_type == LineShape::Type::DP ? 498e6 : 0);
  const Numeric df_ = ls_type == LineShape::Type::DP ? 4e6/Numeric(::N - 1) : 1e9/Numeric(::N - 1);
  
  const Vector f_grid(f0, ::N, df_);
  Vector vmr({0.3, 0.7});
  if constexpr (mod == ModificationInternal::SelfVMR) vmr[0] += ::dVMR;
  else if constexpr (mod == ModificationInternal::AirVMR) vmr[1] += ::dVMR;
  
  const Numeric P = 1e3;
  const Numeric T = (mod == ModificationInternal::T) ? 255 + ::dT : 255;
  const Numeric H = (mod == ModificationInternal::H) ? 50e-6 + ::dH : 50e-6;
  
  return {qid, band, f_grid, vmr, P, T, H};
}

template <LineShape::Type ls_type>
void test_ls() {
  std::cout << "Line Shape: " << ls_type << "\n";
  Index j=0;
  
  Vector f_grid(::N);
  ComplexVector F(::N), N(::N);
  ComplexMatrix dF(::N, ::M, 0), dN(::N, ::M, 0);
  ArrayOfComplexVector dF_mod(::M, ComplexVector(::N, 0)), dN_mod(::M, ComplexVector(::N, 0));
  ComplexMatrix dFnull(0, 0), dNnull(0, 0);
  ArrayOfRetrievalQuantity jacobian_quantities(::M);
  ArrayOfRetrievalQuantity jacobian_quantities_null(0);
  EnergyLevelMap nlte;
  SpeciesAuxData partition_functions;
  fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(partition_functions);
  
  // Ignore sparsity here
  const Vector f_grid_sparse(0);
  LineShape::ComputeData sparse_com(f_grid_sparse, {}, false);
  
  constexpr Index TN = 1;
  
  // Compute
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::None, ls_type>();
    f_grid = f_g; f_grid -= 1e9; f_grid /= 1e6;
    jacobian_quantities[0].Target() = JacobianTarget(Jacobian::Atm::Temperature);
    jacobian_quantities[1].Target() = JacobianTarget(Jacobian::Atm::WindMagnitude);
    jacobian_quantities[2].Target() = JacobianTarget(Jacobian::Atm::MagneticMagnitude);
    jacobian_quantities[3].Target() = JacobianTarget(Jacobian::Line::Center, qid);
    jacobian_quantities[4].Target() = JacobianTarget(Jacobian::Line::Strength, qid);
    jacobian_quantities[5].Target() = JacobianTarget(Jacobian::Line::ShapeG0X0, qid, -1);
    jacobian_quantities[6].Target() = JacobianTarget(Jacobian::Line::ShapeD0X0, qid, -1);
    jacobian_quantities[7].Target() = JacobianTarget(Jacobian::Line::ShapeG2X0, qid, -1);
    jacobian_quantities[8].Target() = JacobianTarget(Jacobian::Line::ShapeD2X0, qid, -1);
    jacobian_quantities[9].Target() = JacobianTarget(Jacobian::Line::ShapeETAX0, qid, -1);
    jacobian_quantities[10].Target() = JacobianTarget(Jacobian::Line::ShapeFVCX0, qid, -1);
    jacobian_quantities[11].Target() = JacobianTarget(Jacobian::Line::ShapeYX1, qid, -1);
    jacobian_quantities[12].Target() = JacobianTarget(Jacobian::Line::ShapeGX1, qid, -1);
    jacobian_quantities[13].Target() = JacobianTarget(Jacobian::Line::ShapeDVX0, qid, -1);
    jacobian_quantities[14].Target() = JacobianTarget(Jacobian::Line::ShapeG0X0, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[15].Target() = JacobianTarget(Jacobian::Line::ShapeD0X0, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[16].Target() = JacobianTarget(Jacobian::Line::ShapeG2X0, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[17].Target() = JacobianTarget(Jacobian::Line::ShapeD2X0, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[18].Target() = JacobianTarget(Jacobian::Line::ShapeETAX0, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[19].Target() = JacobianTarget(Jacobian::Line::ShapeFVCX0, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[20].Target() = JacobianTarget(Jacobian::Line::ShapeYX1, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[21].Target() = JacobianTarget(Jacobian::Line::ShapeGX1, qid, std::numeric_limits<Index>::max());
    jacobian_quantities[22].Target() = JacobianTarget(Jacobian::Line::ShapeDVX0, qid, std::numeric_limits<Index>::max());
    qid.SetAll();
    jacobian_quantities[23].Target() = JacobianTarget(Jacobian::Line::VMR, qid);  // Self is increasing
    jacobian_quantities[24].Target() = JacobianTarget(Jacobian::Line::VMR, QuantumIdentifier(QuantumIdentifier::ALL, 1, 0));  // Another species... main point: it's air
    LineShape::ComputeData com(f_g, jacobian_quantities, false);
    ArrayOfTimeStep dt(TN);
    std::cout << "all derivs...\n";
    for (Index i=0; i<TN; i++) {
      Time start;
      LineShape::compute(com, sparse_com, band, jacobian_quantities, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
      dt[i] = Time() - start;
    }
    std::sort(dt.begin(), dt.end());
    std::cerr << "Slow: " << dt[TN-1] << '\n' << "Mean: " << mean(dt) << '\n' << "Med.: " << median(dt) << '\n' << "Fast: " << dt[0] << '\n';
    
    LineShape::ComputeData com2(f_g, jacobian_quantities, false);
    LineShape::compute(com2, sparse_com, band, jacobian_quantities, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    F = com2.F;
    N = com2.N;
    dF = com2.dF;
    dN = com2.dN;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::T, ls_type>();
    
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    ArrayOfTimeStep dt(TN);
    std::cout << "\njust forward calcs...\n";
    for (Index i=0; i<TN; i++) {
      Time start;
      LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
      dt[i] = Time() - start;
    }
    std::sort(dt.begin(), dt.end());
    std::cerr << "Slow: " << dt[TN-1] << '\n' << "Mean: " << mean(dt) << '\n' << "Med.: " << median(dt) << '\n' << "Fast: " << dt[0] << '\n';
    
    LineShape::ComputeData com2(f_g, jacobian_quantities_null, false);
    LineShape::compute(com2, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com2.F;
    dN_mod[j] = com2.N;
    
    dF_mod[j] -= F;
    dF_mod[j] /= ::dT;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::f, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::df;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::H, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dH;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::F0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dF0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::I0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dI0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfD0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfD2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfETA, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dETA;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfFVC, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dFVC;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfY, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dY;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfDV, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dDV;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirD0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirD2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirETA, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dETA;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirFVC, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dFVC;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirY, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dY;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirDV, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dDV;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfVMR, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dVMR;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirVMR, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte,  partition_functions.getParamType(band.QuantumIdentity()), partition_functions.getParam(band.QuantumIdentity()), vmr, 1, 1, P, T, H, 0, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dVMR;
    j++;
  }
  
  ARTSGUI::plot(f_grid, F.real(), f_grid, F.imag());
  for (Index i=0; i<::M; i++) {
    bool all_zero = true;
    for (Index iv=0; iv<::N; iv++) all_zero = all_zero and dF(iv, i) == Complex(0, 0) and dF_mod[i][iv] == Complex(0, 0);
    if (all_zero) continue;  // plot only if some are non-zero
    
    std::cout << jacobian_quantities[i].Target() << '\n';
    ARTSGUI::plot(f_grid, dF(joker, i).real(), f_grid, dF(joker, i).imag(), f_grid, dF_mod[i].real(), f_grid, dF_mod[i].imag());
  }
  std::cout << '\n';
}

void test_norm() {
  
  [[maybe_unused]]
  auto [qid, band, f_grid, vmr, P, T, H] = line_model<ModificationInternal::None, LineShape::Type::VP>();
  Vector N_nonorm(::N);
  Vector N_vvh(::N);
  Vector N_vvw(::N);
  Vector N_rq(::N);
  
  ArrayOfVector dN_nonorm(3, Vector(::N));
  ArrayOfVector dN_vvh(3, Vector(::N));
  ArrayOfVector dN_vvw(3, Vector(::N));
  ArrayOfVector dN_rq(3, Vector(::N));
  
  ArrayOfVector dN_nonorm_mod(3, Vector(::N));
  ArrayOfVector dN_vvh_mod(3, Vector(::N));
  ArrayOfVector dN_vvw_mod(3, Vector(::N));
  ArrayOfVector dN_rq_mod(3, Vector(::N));
  
  LineShape::Normalizer lsn_nonorm = LineShape::Nonorm{};
  LineShape::Normalizer lsn_vvh = LineShape::VanVleckHuber(band.Line(0).F0(), T);
  LineShape::Normalizer lsn_vvw = LineShape::VanVleckWeisskopf(band.Line(0).F0());
  LineShape::Normalizer lsn_rq = LineShape::RosenkranzQuadratic(band.Line(0).F0(), T);
  for (Index i=0; i<::N; i++) {
    N_nonorm[i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_nonorm);
    N_vvh[i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_vvh);
    N_vvw[i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_vvw);
    N_rq[i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_rq);
    dN_nonorm[0][i] = std::visit([T=T,f=f_grid[i]](auto&& Norm){ return Norm.dNdT(T, f);}, lsn_nonorm);
    dN_vvh[0][i] = std::visit([T=T,f=f_grid[i]](auto&& Norm){ return Norm.dNdT(T, f);}, lsn_vvh);
    dN_vvw[0][i] = std::visit([T=T,f=f_grid[i]](auto&& Norm){ return Norm.dNdT(T, f);}, lsn_vvw);
    dN_rq[0][i] = std::visit([T=T,f=f_grid[i]](auto&& Norm){ return Norm.dNdT(T, f);}, lsn_rq);
    dN_nonorm[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm.dNdf(f);}, lsn_nonorm);
    dN_vvh[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm.dNdf(f);}, lsn_vvh);
    dN_vvw[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm.dNdf(f);}, lsn_vvw);
    dN_rq[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm.dNdf(f);}, lsn_rq);
    dN_nonorm[2][i] = std::visit([](auto&& Norm){ return Norm.dNdF0();}, lsn_nonorm);
    dN_vvh[2][i] = std::visit([](auto&& Norm){ return Norm.dNdF0();}, lsn_vvh);
    dN_vvw[2][i] = std::visit([](auto&& Norm){ return Norm.dNdF0();}, lsn_vvw);
    dN_rq[2][i] = std::visit([](auto&& Norm){ return Norm.dNdF0();}, lsn_rq);
    
    // Do derivatives last...
    dN_nonorm_mod[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f+::df);}, lsn_nonorm);
    dN_vvh_mod[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f+::df);}, lsn_vvh);
    dN_vvw_mod[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f+::df);}, lsn_vvw);
    dN_rq_mod[1][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f+::df);}, lsn_rq);
  }
  
  LineShape::Normalizer lsn_nonormdT = LineShape::Nonorm{};
  LineShape::Normalizer lsn_vvhdT = LineShape::VanVleckHuber(band.Line(0).F0(), T+::dT);
  LineShape::Normalizer lsn_vvwdT = LineShape::VanVleckWeisskopf(band.Line(0).F0());
  LineShape::Normalizer lsn_rqdT = LineShape::RosenkranzQuadratic(band.Line(0).F0(), T+::dT);
  for (Index i=0; i<::N; i++) {
    dN_nonorm_mod[0][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_nonormdT);
    dN_vvh_mod[0][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_vvhdT);
    dN_vvw_mod[0][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_vvwdT);
    dN_rq_mod[0][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_rqdT);
  }
  
  LineShape::Normalizer lsn_nonormdF0 = LineShape::Nonorm{};
  LineShape::Normalizer lsn_vvhdF0 = LineShape::VanVleckHuber(band.Line(0).F0()+::dF0, T);
  LineShape::Normalizer lsn_vvwdF0 = LineShape::VanVleckWeisskopf(band.Line(0).F0()+::dF0);
  LineShape::Normalizer lsn_rqdF0 = LineShape::RosenkranzQuadratic(band.Line(0).F0()+::dF0, T);
  for (Index i=0; i<::N; i++) {
    dN_nonorm_mod[2][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_nonormdF0);
    dN_vvh_mod[2][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_vvhdF0);
    dN_vvw_mod[2][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_vvwdF0);
    dN_rq_mod[2][i] = std::visit([f=f_grid[i]](auto&& Norm){ return Norm(f);}, lsn_rqdF0);
  }
  
  dN_nonorm_mod[0] -= N_nonorm; dN_nonorm_mod[0] /= ::dT;
  dN_nonorm_mod[1] -= N_nonorm; dN_nonorm_mod[1] /= ::df;
  dN_nonorm_mod[2] -= N_nonorm; dN_nonorm_mod[2] /= ::dF0;
  dN_vvh_mod[0] -= N_vvh; dN_vvh_mod[0] /= ::dT;
  dN_vvh_mod[1] -= N_vvh; dN_vvh_mod[1] /= ::df;
  dN_vvh_mod[2] -= N_vvh; dN_vvh_mod[2] /= ::dF0;
  dN_vvw_mod[0] -= N_vvw; dN_vvw_mod[0] /= ::dT;
  dN_vvw_mod[1] -= N_vvw; dN_vvw_mod[1] /= ::df;
  dN_vvw_mod[2] -= N_vvw; dN_vvw_mod[2] /= ::dF0;
  dN_rq_mod[0] -= N_rq; dN_rq_mod[0] /= ::dT;
  dN_rq_mod[1] -= N_rq; dN_rq_mod[1] /= ::df;
  dN_rq_mod[2] -= N_rq; dN_rq_mod[2] /= ::dF0;
  
  ARTSGUI::plot(f_grid, N_nonorm, f_grid, N_vvh, f_grid, N_vvw, f_grid, N_rq);
  for (Index i=0; i<3; i++) {
    std::cout << i << '\n';
    if (max(dN_nonorm[i]) not_eq 0 or min(dN_nonorm[i]) not_eq 0 or max(dN_nonorm_mod[i]) not_eq 0 or min(dN_nonorm_mod[i]) not_eq 0) {
      std::cout << "NONORM" << '\n';
      ARTSGUI::plot(f_grid, dN_nonorm[i], f_grid, dN_nonorm_mod[i]);
    }
    if (max(dN_vvh[i]) not_eq 0 or min(dN_vvh[i]) not_eq 0 or max(dN_vvh_mod[i]) not_eq 0 or min(dN_vvh_mod[i]) not_eq 0) {
      std::cout << "VVH" << '\n';
      ARTSGUI::plot(f_grid, dN_vvh[i], f_grid, dN_vvh_mod[i]);
    }
    if (max(dN_vvw[i]) not_eq 0 or min(dN_vvw[i]) not_eq 0 or max(dN_vvw_mod[i]) not_eq 0 or min(dN_vvw_mod[i]) not_eq 0) {
      std::cout << "VVW" << '\n';
      ARTSGUI::plot(f_grid, dN_vvw[i], f_grid, dN_vvw_mod[i]);
    }
    if (max(dN_rq[i]) not_eq 0 or min(dN_rq[i]) not_eq 0 or max(dN_rq_mod[i]) not_eq 0 or min(dN_rq_mod[i]) not_eq 0) {
      std::cout << "RQ" << '\n';
      ARTSGUI::plot(f_grid, dN_rq[i], f_grid, dN_rq_mod[i]);
    }
  }
}

void test_lte_strength() {
  SpeciesAuxData partition_functions;
  fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(partition_functions);
  
  Numeric I0=1e10, T0=296, T=255, F0=1e9, E0=1e-20, r=1.0;
  LineShape::LocalThermodynamicEquilibrium x(I0, T0, T, F0, E0,
                                             single_partition_function(T, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                             single_partition_function(T0, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                             dsingle_partition_function_dT(T, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)), r, 0, 0);
  LineShape::LocalThermodynamicEquilibrium dxdI0(I0+1e4, T0, T, F0, E0,
                                                 single_partition_function(T, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                                 single_partition_function(T0, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                                 dsingle_partition_function_dT(T, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)), r, 0, 0);
  LineShape::LocalThermodynamicEquilibrium dxdT(I0, T0, T+0.01, F0, E0,
                                                single_partition_function(T+0.01, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                                single_partition_function(T0, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                                dsingle_partition_function_dT(T+0.01, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)), r, 0, 0);
  LineShape::LocalThermodynamicEquilibrium dxdF0(I0, T0, T, F0+100e3, E0,
                                                 single_partition_function(T, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                                 single_partition_function(T0, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)),
                                                 dsingle_partition_function_dT(T, partition_functions.getParamType(0, 0), partition_functions.getParam(0, 0)), r, 0, 0);
  
  std::cout << x.S << '\n';
  std::cout << x.dSdI0() << ' ' << 'v' << 's' << ' ' << (dxdI0.S - x.S) / 1e4 << '\n';
  std::cout << x.dSdT() << ' ' << 'v' << 's' << ' ' << (dxdT.S - x.S) / 0.01 << '\n';
  std::cout << x.dSdF0() << ' ' << 'v' << 's' << ' ' << (dxdF0.S - x.S) / 100e3 << '\n';
}

void test_vp_sparse() {
  auto [qid, band, f_gp, vmr, P, T, H] = line_model<ModificationInternal::None, LineShape::Type::HTP>();
  P *= 1000;
  band.F0(0) = 1500e9;
  const Vector f_g(700e9, 20001, 10e7);
  auto band2 = band;
  band2.Cutoff(Absorption::CutoffType::ByLine);
  band2.CutoffFreqValue(750e9);
  auto band3 = band2;
  band3.CutoffFreqValue(333.33e9);
  
  EnergyLevelMap nlte;
  SpeciesAuxData partition_functions;
  fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(partition_functions);
  
  const Numeric dfreq = 45e9;
  
  LineShape::ComputeData com_full(f_g, {}, false);
  LineShape::ComputeData sparse_com_full(Vector(0), {}, false);
  LineShape::compute(com_full, sparse_com_full, band, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
  LineShape::compute(com_full, sparse_com_full, band2, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
  LineShape::compute(com_full, sparse_com_full, band3, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::None);
  
  Vector sparse3_f_grid = LineShape::triple_sparse_f_grid(f_g, 30e9);
  LineShape::ComputeData com3(f_g, {}, false);
  LineShape::ComputeData sparse_com3(sparse3_f_grid, {}, false);
  LineShape::compute(com3, sparse_com3, band, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::QuadraticIndependent);
  LineShape::compute(com3, sparse_com3, band2, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::QuadraticIndependent);
  LineShape::compute(com3, sparse_com3, band3, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::QuadraticIndependent);
  
  Vector sparsel_f_grid = LineShape::linear_sparse_f_grid(f_g, 10e9);
  LineShape::ComputeData coml(f_g, {}, false);
  LineShape::ComputeData sparse_coml(sparsel_f_grid, {}, false);
  LineShape::compute(coml, sparse_coml, band, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::LinearEven);
  LineShape::compute(coml, sparse_coml, band2, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::LinearEven);
  LineShape::compute(coml, sparse_coml, band3, {}, nlte,
                     partition_functions.getParamType(band.QuantumIdentity()),
                     partition_functions.getParam(band.QuantumIdentity()),
                     vmr, 1, 1, P, T, H, dfreq, true, Zeeman::Polarization::Pi, Options::LblSpeedup::LinearEven);
  
  std::cout << f_g[0]/1e9 << ' ' << f_g[10000]/1e9 << '\n';
  std::cout << com_full.f_grid.nelem()<< ' ' << sparse_com3.f_grid.nelem()<< ' ' << sparse_coml.f_grid.nelem() << '\n';
  
//   ARTSGUI::plot(com_full.f_grid, com_full.F.real(), com3.f_grid, com3.F.real(), coml.f_grid, coml.F.real());
  
//   ARTSGUI::plot(sparse_com3.f_grid, sparse_com3.F.real());
  
  com3.interp_add_triplequad(sparse_com3);
  coml.interp_add_even(sparse_coml);
  
  ARTSGUI::plot(com_full.f_grid, com_full.F.real(),
                com3.f_grid, com3.F.real(),
                coml.f_grid, coml.F.real());
}

int main() try {
  define_species_data();
  define_species_map();
  
  make_wigner_ready(20, 30, 6);
  std::cout << std::setprecision(15);
  
//   test_ls<LineShape::Type::DP>();
//   test_ls<LineShape::Type::LP>();
//   test_ls<LineShape::Type::VP>();
//   test_ls<LineShape::Type::SDVP>();
//   test_ls<LineShape::Type::HTP>();
// 
//   test_norm();
//   
//   test_lte_strength();
  
  
  test_vp_sparse();
  
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
