#include "absorption.h"
#include "artstime.h"
#include "gui/plot.h"
#include "lineshape.h"
#include "species.h"
#include "wigner_functions.h"

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
} // namespace

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
  model[0].G2() = LineShape::ModelParameters(LineShape::TemperatureModel::T3, 200, 20.);
  model[0].D2() = LineShape::ModelParameters(LineShape::TemperatureModel::T0, 20.);
  model[0].G2() = LineShape::ModelParameters(LineShape::TemperatureModel::T3, 200., 20.);
  model[0].ETA() = LineShape::ModelParameters(LineShape::TemperatureModel::T0, 0.1);
  model[0].FVC() = LineShape::ModelParameters(LineShape::TemperatureModel::T2, 1e2, .2, 1.);
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
  QuantumIdentifier qid = QuantumIdentifier("H2O-161 J 3 2");
  
  Absorption::SingleLine line(1e9, 1e10, 1e-20, 1., 3., 1e-14,
                              Zeeman::Model(2, 1.5), lineshape_model<mod>(),
                              "J 3 2");
  
  if constexpr (mod == ModificationInternal::F0) line.F0 += ::dF0;
  if constexpr (mod == ModificationInternal::I0) line.I0 += ::dI0;
  
  AbsorptionLines band(true, true,
                       Absorption::CutoffType::None, Absorption::MirroringType::None,
                       Absorption::PopulationType::LTE, Absorption::NormalizationType::None,
                       ls_type, 296, -1, -1, qid,
                       {Species::Species::Water, Species::Species::Bath}, {line});
  
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
    jacobian_quantities[3].Target() = JacobianTarget(Jacobian::Line::Center, qid, qid.Species());
    jacobian_quantities[4].Target() = JacobianTarget(Jacobian::Line::Strength, qid, qid.Species());
    jacobian_quantities[5].Target() = JacobianTarget(Jacobian::Line::ShapeG0X0, qid, qid.Species());
    jacobian_quantities[6].Target() = JacobianTarget(Jacobian::Line::ShapeD0X0, qid, qid.Species());
    jacobian_quantities[7].Target() = JacobianTarget(Jacobian::Line::ShapeG2X0, qid, qid.Species());
    jacobian_quantities[8].Target() = JacobianTarget(Jacobian::Line::ShapeD2X0, qid, qid.Species());
    jacobian_quantities[9].Target() = JacobianTarget(Jacobian::Line::ShapeETAX0, qid, qid.Species());
    jacobian_quantities[10].Target() = JacobianTarget(Jacobian::Line::ShapeFVCX0, qid, qid.Species());
    jacobian_quantities[11].Target() = JacobianTarget(Jacobian::Line::ShapeYX1, qid, qid.Species());
    jacobian_quantities[12].Target() = JacobianTarget(Jacobian::Line::ShapeGX1, qid, qid.Species());
    jacobian_quantities[13].Target() = JacobianTarget(Jacobian::Line::ShapeDVX0, qid, qid.Species());
    jacobian_quantities[14].Target() = JacobianTarget(Jacobian::Line::ShapeG0X0, qid, Species::Species::Bath);
    jacobian_quantities[15].Target() = JacobianTarget(Jacobian::Line::ShapeD0X0, qid, Species::Species::Bath);
    jacobian_quantities[16].Target() = JacobianTarget(Jacobian::Line::ShapeG2X0, qid, Species::Species::Bath);
    jacobian_quantities[17].Target() = JacobianTarget(Jacobian::Line::ShapeD2X0, qid, Species::Species::Bath);
    jacobian_quantities[18].Target() = JacobianTarget(Jacobian::Line::ShapeETAX0, qid, Species::Species::Bath);
    jacobian_quantities[19].Target() = JacobianTarget(Jacobian::Line::ShapeFVCX0, qid, Species::Species::Bath);
    jacobian_quantities[20].Target() = JacobianTarget(Jacobian::Line::ShapeYX1, qid, Species::Species::Bath);
    jacobian_quantities[21].Target() = JacobianTarget(Jacobian::Line::ShapeGX1, qid, Species::Species::Bath);
    jacobian_quantities[22].Target() = JacobianTarget(Jacobian::Line::ShapeDVX0, qid, Species::Species::Bath);
    qid = QuantumIdentifier(qid.Isotopologue());
    jacobian_quantities[23].Target() = JacobianTarget(Jacobian::Line::VMR, qid, qid.Species());  // Self is increasing
    jacobian_quantities[24].Target() = JacobianTarget(Jacobian::Line::VMR, QuantumIdentifier(Species::Isotopologues[0]), Species::Species::Bath);  // Another species... main point: it's air
    LineShape::ComputeData com(f_g, jacobian_quantities, false);
    ArrayOfTimeStep dt(TN);
    std::cout << "all derivs...\n";
    for (Index i=0; i<TN; i++) {
      Time start;
      LineShape::compute(com, sparse_com, band, jacobian_quantities, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
      dt[i] = Time() - start;
    }
    std::sort(dt.begin(), dt.end());
    std::cerr << "Slow: " << dt[TN-1] << '\n' << "Mean: " << mean(dt) << '\n' << "Med.: " << median(dt) << '\n' << "Fast: " << dt[0] << '\n';
    
    LineShape::ComputeData com2(f_g, jacobian_quantities, false);
    LineShape::compute(com2, sparse_com, band, jacobian_quantities, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
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
      LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
      dt[i] = Time() - start;
    }
    std::sort(dt.begin(), dt.end());
    std::cerr << "Slow: " << dt[TN-1] << '\n' << "Mean: " << mean(dt) << '\n' << "Med.: " << median(dt) << '\n' << "Fast: " << dt[0] << '\n';
    
    LineShape::ComputeData com2(f_g, jacobian_quantities_null, false);
    LineShape::compute(com2, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com2.F;
    dN_mod[j] = com2.N;
    
    dF_mod[j] -= F;
    dF_mod[j] /= ::dT;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::f, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::df;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::H, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dH;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::F0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dF0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::I0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dI0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfD0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfD2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfETA, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dETA;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfFVC, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dFVC;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfY, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dY;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfDV, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dDV;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirD0, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirD2, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirETA, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dETA;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirFVC, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dFVC;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirY, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dY;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirDV, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dDV;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfVMR, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
    dF_mod[j] = com.F;
    dN_mod[j] = com.N;
    dF_mod[j] -= F;
    dF_mod[j] /= ::dVMR;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirVMR, ls_type>();
    LineShape::ComputeData com(f_g, jacobian_quantities_null, false);
    LineShape::compute(com, sparse_com, band, jacobian_quantities_null, nlte, vmr, {}, 1, 1, P, T, H, 0, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
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
  auto [qid, band, f_grid2, vmr, P, T, H] = line_model<ModificationInternal::None, LineShape::Type::VP>();
  Vector f_grid(::N);
  Numeric f0 = f_grid2[0];
  for (Index i=0; i<::N; i++) {
    f_grid[i] = f0;
    f0 += 3 * (f_grid2[1] - f_grid2[0]);
  }
  Vector N_nonorm(::N);
  Vector N_vvh(::N);
  Vector N_vvw(::N);
  Vector N_rq(::N);
  Vector N_sfs(::N);
  
  ArrayOfVector dN_nonorm(3, Vector(::N));
  ArrayOfVector dN_vvh(3, Vector(::N));
  ArrayOfVector dN_vvw(3, Vector(::N));
  ArrayOfVector dN_rq(3, Vector(::N));
  ArrayOfVector dN_sfs(3, Vector(::N));
  
  ArrayOfVector dN_nonorm_mod(3, Vector(::N));
  ArrayOfVector dN_vvh_mod(3, Vector(::N));
  ArrayOfVector dN_vvw_mod(3, Vector(::N));
  ArrayOfVector dN_rq_mod(3, Vector(::N));
  ArrayOfVector dN_sfs_mod(3, Vector(::N));
  
  LineShape::Normalizer lsn_nonorm(Absorption::NormalizationType::None, band.lines.front().F0, T);
  LineShape::Normalizer lsn_vvh(Absorption::NormalizationType::VVH, band.lines.front().F0, T);
  LineShape::Normalizer lsn_vvw(Absorption::NormalizationType::VVW, band.lines.front().F0, T);
  LineShape::Normalizer lsn_rq(Absorption::NormalizationType::RQ, band.lines.front().F0, T);
  LineShape::Normalizer lsn_sfs(Absorption::NormalizationType::SFS, band.lines.front().F0, T);
  for (Index i=0; i<::N; i++) {
    const Numeric& f = f_grid[i];
    
    N_nonorm[i] = lsn_nonorm(f);
    N_vvh[i] = lsn_vvh(f);
    N_vvw[i] = lsn_vvw(f);
    N_rq[i] = lsn_rq(f);
    N_sfs[i] = lsn_sfs(f);
    
    dN_nonorm[0][i] = lsn_nonorm.dNdT(T, f);
    dN_vvh[0][i] = lsn_vvh.dNdT(T, f);
    dN_vvw[0][i] = lsn_vvw.dNdT(T, f);
    dN_rq[0][i] = lsn_rq.dNdT(T, f);
    dN_sfs[0][i] = lsn_sfs.dNdT(T, f);
    
    dN_nonorm[1][i] = lsn_nonorm.dNdf(f);
    dN_vvh[1][i] = lsn_vvh.dNdf(f);
    dN_vvw[1][i] = lsn_vvw.dNdf(f);
    dN_rq[1][i] = lsn_rq.dNdf(f);
    dN_sfs[1][i] = lsn_sfs.dNdf(f);
    
    dN_nonorm[2][i] = lsn_nonorm.dNdF0();
    dN_vvh[2][i] = lsn_vvh.dNdF0();
    dN_vvw[2][i] = lsn_vvw.dNdF0();
    dN_rq[2][i] = lsn_rq.dNdF0();
    dN_sfs[2][i] = lsn_sfs.dNdF0();
    
    // Do derivatives last...
    dN_nonorm_mod[1][i] = lsn_nonorm(f+::df);
    dN_vvh_mod[1][i] = lsn_vvh(f+::df);
    dN_vvw_mod[1][i] = lsn_vvw(f+::df);
    dN_rq_mod[1][i] = lsn_rq(f+::df);
    dN_sfs_mod[1][i] = lsn_sfs(f+::df);
  }
  
  LineShape::Normalizer lsn_nonormdT(Absorption::NormalizationType::None, band.lines.front().F0, T+::dT);
  LineShape::Normalizer lsn_vvhdT(Absorption::NormalizationType::VVH, band.lines.front().F0, T+::dT);
  LineShape::Normalizer lsn_vvwdT(Absorption::NormalizationType::VVW, band.lines.front().F0, T+::dT);
  LineShape::Normalizer lsn_rqdT(Absorption::NormalizationType::RQ, band.lines.front().F0, T+::dT);
  LineShape::Normalizer lsn_sfsdT(Absorption::NormalizationType::SFS, band.lines.front().F0, T+::dT);
  for (Index i=0; i<::N; i++) {
    const Numeric& f = f_grid[i];

    dN_nonorm_mod[0][i] = lsn_nonormdT(f);
    dN_vvh_mod[0][i] = lsn_vvhdT(f);
    dN_vvw_mod[0][i] = lsn_vvwdT(f);
    dN_rq_mod[0][i] = lsn_rqdT(f);
    dN_sfs_mod[0][i] = lsn_sfsdT(f);
  }
  
  LineShape::Normalizer lsn_nonormdF0(Absorption::NormalizationType::None, band.lines.front().F0+::dF0, T);
  LineShape::Normalizer lsn_vvhdF0(Absorption::NormalizationType::VVH, band.lines.front().F0+::dF0, T);
  LineShape::Normalizer lsn_vvwdF0(Absorption::NormalizationType::VVW, band.lines.front().F0+::dF0, T);
  LineShape::Normalizer lsn_rqdF0(Absorption::NormalizationType::RQ, band.lines.front().F0+::dF0, T);
  LineShape::Normalizer lsn_sfsdF0(Absorption::NormalizationType::SFS, band.lines.front().F0+::dF0, T);
  for (Index i=0; i<::N; i++) {
    const Numeric& f = f_grid[i];

    dN_nonorm_mod[2][i] = lsn_nonormdF0(f);
    dN_vvh_mod[2][i] = lsn_vvhdF0(f);
    dN_vvw_mod[2][i] = lsn_vvwdF0(f);
    dN_rq_mod[2][i] = lsn_rqdF0(f);
    dN_sfs_mod[2][i] = lsn_sfsdF0(f);
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
  dN_sfs_mod[0] -= N_sfs; dN_sfs_mod[0] /= ::dT;
  dN_sfs_mod[1] -= N_sfs; dN_sfs_mod[1] /= ::df;
  dN_sfs_mod[2] -= N_sfs; dN_sfs_mod[2] /= ::dF0;
  
  //! Plot results, so modify frequency to be readable
  f_grid /= band.lines.front().F0;
  
  ARTSGUI::PlotConfig::X = "Frequency [line center]";
  ARTSGUI::PlotConfig::Y = "Norm";
  ARTSGUI::PlotConfig::Legend = {"None", "VVH", "VVW", "RQ", "SFS"};
  ARTSGUI::PlotConfig::Title = "Normalizers";
  ARTSGUI::plot(f_grid, N_nonorm, f_grid, N_vvh, f_grid, N_vvw, f_grid, N_rq, f_grid, N_sfs);
  
  ARTSGUI::PlotConfig::Legend = {"Analytical", "Perturbed"};
  for (Index i=0; i<3; i++) {
    if (i == 0) {
      ARTSGUI::PlotConfig::Y = "Temperature Derivative";
    } else if (i == 1) {
      ARTSGUI::PlotConfig::Y = "Frequency Derivative";
    } else if (i == 2) {
      ARTSGUI::PlotConfig::Y = "Line Center Derivative";
    } else {
      ARTSGUI::PlotConfig::Y = "Unknown Derivative";
    }
    
    if (max(dN_nonorm[i]) not_eq 0 or min(dN_nonorm[i]) not_eq 0 or max(dN_nonorm_mod[i]) not_eq 0 or min(dN_nonorm_mod[i]) not_eq 0) {
      ARTSGUI::PlotConfig::Title = "No Normalizer";
      ARTSGUI::plot(f_grid, dN_nonorm[i], f_grid, dN_nonorm_mod[i]);
    }
    
    if (max(dN_vvh[i]) not_eq 0 or min(dN_vvh[i]) not_eq 0 or max(dN_vvh_mod[i]) not_eq 0 or min(dN_vvh_mod[i]) not_eq 0) {
      ARTSGUI::PlotConfig::Title = "VVH Normalizer";
      ARTSGUI::plot(f_grid, dN_vvh[i], f_grid, dN_vvh_mod[i]);
    }
    
    if (max(dN_vvw[i]) not_eq 0 or min(dN_vvw[i]) not_eq 0 or max(dN_vvw_mod[i]) not_eq 0 or min(dN_vvw_mod[i]) not_eq 0) {
      ARTSGUI::PlotConfig::Title = "VVW Normalizer";
      ARTSGUI::plot(f_grid, dN_vvw[i], f_grid, dN_vvw_mod[i]);
    }
    
    if (max(dN_rq[i]) not_eq 0 or min(dN_rq[i]) not_eq 0 or max(dN_rq_mod[i]) not_eq 0 or min(dN_rq_mod[i]) not_eq 0) {
      ARTSGUI::PlotConfig::Title = "RQ Normalizer";
      ARTSGUI::plot(f_grid, dN_rq[i], f_grid, dN_rq_mod[i]);
    }
    
    if (max(dN_sfs[i]) not_eq 0 or min(dN_sfs[i]) not_eq 0 or max(dN_sfs_mod[i]) not_eq 0 or min(dN_sfs_mod[i]) not_eq 0) {
      ARTSGUI::PlotConfig::Title = "SFS Normalizer";
      ARTSGUI::plot(f_grid, dN_sfs[i], f_grid, dN_sfs_mod[i]);
    }
  }
}

void test_lte_strength() {
  Numeric I0=1e10, T0=296, T=255, F0=1e9, E0=1e-20, r=1.0;
  LineShape::LocalThermodynamicEquilibrium x(I0, T0, T, F0, E0,
                                             single_partition_function(T, Species::Tag("H2O-161").Isotopologue()),
                                             single_partition_function(T0, Species::Tag("H2O-161").Isotopologue()),
                                             dsingle_partition_function_dT(T, Species::Tag("H2O-161").Isotopologue()), r, 0, 0);
  LineShape::LocalThermodynamicEquilibrium dxdI0(I0+1e4, T0, T, F0, E0,
                                                 single_partition_function(T, Species::Tag("H2O-161").Isotopologue()),
                                                 single_partition_function(T0, Species::Tag("H2O-161").Isotopologue()),
                                                 dsingle_partition_function_dT(T, Species::Tag("H2O-161").Isotopologue()), r, 0, 0);
  LineShape::LocalThermodynamicEquilibrium dxdT(I0, T0, T+0.01, F0, E0,
                                                single_partition_function(T+0.01, Species::Tag("H2O-161").Isotopologue()),
                                                single_partition_function(T0, Species::Tag("H2O-161").Isotopologue()),
                                                dsingle_partition_function_dT(T+0.01, Species::Tag("H2O-161").Isotopologue()), r, 0, 0);
  LineShape::LocalThermodynamicEquilibrium dxdF0(I0, T0, T, F0+100e3, E0,
                                                 single_partition_function(T, Species::Tag("H2O-161").Isotopologue()),
                                                 single_partition_function(T0, Species::Tag("H2O-161").Isotopologue()),
                                                 dsingle_partition_function_dT(T, Species::Tag("H2O-161").Isotopologue()), r, 0, 0);
  
  std::cout << x.S << '\n';
  std::cout << x.dSdI0() << ' ' << 'v' << 's' << ' ' << (dxdI0.S - x.S) / 1e4 << '\n';
  std::cout << x.dSdT() << ' ' << 'v' << 's' << ' ' << (dxdT.S - x.S) / 0.01 << '\n';
  std::cout << x.dSdF0() << ' ' << 'v' << 's' << ' ' << (dxdF0.S - x.S) / 100e3 << '\n';
}

void test_sparse() {
  auto [qid, band, f_gp, vmr, P, T, H] = line_model<ModificationInternal::None, LineShape::Type::HTP>();
  P *= 1000;
  band.lines.front().F0 = 1500e9;
  Vector f_g(500e9, 20001, 10e7);
  auto band2 = band;
  band2.cutoff = Absorption::CutoffType::ByLine;
  band2.cutofffreq = 750e9;
  auto band3 = band2;
  band3.cutofffreq = 333.33e9;
  
  EnergyLevelMap nlte;
  
  const Numeric dfreq = 45e9;
  
  LineShape::ComputeData com_full(f_g, {}, false);
  LineShape::ComputeData sparse_com_full(Vector(0), {}, false);
  LineShape::compute(com_full, sparse_com_full, band, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
  LineShape::compute(com_full, sparse_com_full, band2, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
  LineShape::compute(com_full, sparse_com_full, band3, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::None, false);
  
  Vector sparse3_f_grid = LineShape::triple_sparse_f_grid(f_g, 30e9);
  LineShape::ComputeData com3(f_g, {}, false);
  LineShape::ComputeData sparse_com3(sparse3_f_grid, {}, false);
  LineShape::compute(com3, sparse_com3, band, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::QuadraticIndependent, false);
  LineShape::compute(com3, sparse_com3, band2, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::QuadraticIndependent, false);
  LineShape::compute(com3, sparse_com3, band3, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::QuadraticIndependent, false);
  
  Vector sparsel_f_grid = LineShape::linear_sparse_f_grid(f_g, 10e9);
  std::cout << sparsel_f_grid <<'\n';
  LineShape::ComputeData coml(f_g, {}, false);
  LineShape::ComputeData sparse_coml(sparsel_f_grid, {}, false);
  LineShape::compute(coml, sparse_coml, band, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::LinearIndependent, false);
  LineShape::compute(coml, sparse_coml, band2, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::LinearIndependent, false);
  LineShape::compute(coml, sparse_coml, band3, {}, nlte,
                     vmr, {}, 1, 1, P, T, H, dfreq, Zeeman::Polarization::None, Options::LblSpeedup::LinearIndependent, false);
  Vector f = f_g;
  f /= 1e9;
  
  std::cout << f[0] << ' ' << f[f_g.nelem() - 1] << '\n';
  std::cout << com_full.f_grid.nelem()<< ' ' << sparse_com3.f_grid.nelem()<< ' ' << sparse_coml.f_grid.nelem() << '\n';
  
  ARTSGUI::plot(f, com_full.F.real(), f, com3.F.real(), f, coml.F.real());
  
  com3.interp_add_triplequad(sparse_com3);
  coml.interp_add_even(sparse_coml);
  
  ARTSGUI::plot(f, com_full.F.real(),
                f, com3.F.real(),
                f, coml.F.real());
}

int main() try {
  make_wigner_ready(20, 30, 6);
  std::cout << std::setprecision(15);
  
//   test_ls<LineShape::Type::DP>();
//   test_ls<LineShape::Type::LP>();
   test_ls<LineShape::Type::VP>();
//   test_ls<LineShape::Type::SDVP>();
//   test_ls<LineShape::Type::HTP>();
// 
  test_norm();
//   
   test_lte_strength();
//   
   test_sparse();
//   
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
