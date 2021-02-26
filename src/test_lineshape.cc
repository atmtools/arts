#include "lineshape.h"
#include "gui/plot.h"

namespace {
  constexpr Index N = 10001;
  constexpr Index M = 24;
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
  constexpr Numeric df = 1e2;
  constexpr Numeric dF0 = 1e2;
  static_assert (N > 1);
}

ENUMCLASS (ModificationInternal, char,
           None,
           T,
           f,
           H,
           F0,
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
  
  Absorption::SingleLine line(1e9, 1e-10, 1e-12, 1., 3., 1e-14,
                              Zeeman::Model(2, 1.5), lineshape_model<mod>(),
                              {low[QuantumNumberType::J]}, {upp[QuantumNumberType::J]});
  
  if constexpr (mod == ModificationInternal::F0) line.F0() += ::dF0;
  
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
  Index j=0;
  
  Vector f_grid(::N);
  ComplexVector F(::N);
  ArrayOfComplexVector dF(::M, ComplexVector(::N));
  ArrayOfComplexVector dF_mod(::M, ComplexVector(::N));
  ArrayOfComplexVector dFnull(0, ComplexVector(::N));
  ArrayOfRetrievalQuantity jacobian_quantities(::M);
  ArrayOfRetrievalQuantity jacobian_quantities_null(0);
  
  // Compute
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::None, ls_type>();
    f_grid = f_g; f_grid -= 1e9; f_grid /= 1e6;
    jacobian_quantities[0].Target() = JacobianTarget(Jacobian::Atm::Temperature);
    jacobian_quantities[1].Target() = JacobianTarget(Jacobian::Atm::WindMagnitude);
    jacobian_quantities[2].Target() = JacobianTarget(Jacobian::Atm::MagneticMagnitude);
    jacobian_quantities[3].Target() = JacobianTarget(Jacobian::Line::Center, qid);
    jacobian_quantities[4].Target() = JacobianTarget(Jacobian::Line::ShapeG0X0, qid, 0);
    jacobian_quantities[5].Target() = JacobianTarget(Jacobian::Line::ShapeD0X0, qid, 0);
    jacobian_quantities[6].Target() = JacobianTarget(Jacobian::Line::ShapeG2X0, qid, 0);
    jacobian_quantities[7].Target() = JacobianTarget(Jacobian::Line::ShapeD2X0, qid, 0);
    jacobian_quantities[8].Target() = JacobianTarget(Jacobian::Line::ShapeETAX0, qid, 0);
    jacobian_quantities[9].Target() = JacobianTarget(Jacobian::Line::ShapeFVCX0, qid, 0);
    jacobian_quantities[10].Target() = JacobianTarget(Jacobian::Line::ShapeYX1, qid, 0);
    jacobian_quantities[11].Target() = JacobianTarget(Jacobian::Line::ShapeGX1, qid, 0);
    jacobian_quantities[12].Target() = JacobianTarget(Jacobian::Line::ShapeDVX0, qid, 0);
    jacobian_quantities[13].Target() = JacobianTarget(Jacobian::Line::ShapeG0X0, qid, 1);
    jacobian_quantities[14].Target() = JacobianTarget(Jacobian::Line::ShapeD0X0, qid, 1);
    jacobian_quantities[15].Target() = JacobianTarget(Jacobian::Line::ShapeG2X0, qid, 1);
    jacobian_quantities[16].Target() = JacobianTarget(Jacobian::Line::ShapeD2X0, qid, 1);
    jacobian_quantities[17].Target() = JacobianTarget(Jacobian::Line::ShapeETAX0, qid, 1);
    jacobian_quantities[18].Target() = JacobianTarget(Jacobian::Line::ShapeFVCX0, qid, 1);
    jacobian_quantities[19].Target() = JacobianTarget(Jacobian::Line::ShapeYX1, qid, 1);
    jacobian_quantities[20].Target() = JacobianTarget(Jacobian::Line::ShapeGX1, qid, 1);
    jacobian_quantities[21].Target() = JacobianTarget(Jacobian::Line::ShapeDVX0, qid, 1);
    qid.SetAll();
    jacobian_quantities[22].Target() = JacobianTarget(Jacobian::Line::VMR, qid);  // Self is increasing
    jacobian_quantities[23].Target() = JacobianTarget(Jacobian::Line::VMR, QuantumIdentifier(QuantumIdentifier::ALL, 1, 0));  // Another species... main point: it's air
    LineShape::compute(F, dF, f_g, band, jacobian_quantities, vmr, P, T, H, true, Zeeman::Polarization::Pi);
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::T, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dT;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::f, ls_type>();
    LineShape::compute(dF_mod[1], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::df;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::H, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dH;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::F0, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dF0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG0, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfD0, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG2, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfD2, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfETA, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dETA;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfFVC, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dFVC;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfY, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dY;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfG, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfDV, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dDV;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG0, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirD0, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD0;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG2, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirD2, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dD2;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirETA, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dETA;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirFVC, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dFVC;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirY, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dY;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirG, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dG;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirDV, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dDV;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::SelfVMR, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dVMR;
    j++;
  }
  {
    auto [qid, band, f_g, vmr, P, T, H] = line_model<ModificationInternal::AirVMR, ls_type>();
    LineShape::compute(dF_mod[j], dFnull, f_g, band, jacobian_quantities_null, vmr, P, T, H, true, Zeeman::Polarization::Pi);
    dF_mod[j] -= F;
    dF_mod[j] /= ::dVMR;
    j++;
  }
  
  std::cout << "Line Shape: " << ls_type << "\n";
  ARTSGUI::plot(f_grid, F.real(), f_grid, F.imag());
  for (Index i=0; i<::M; i++) {
    bool all_zero = true;
    for (Index iv=0; iv<::N; iv++) all_zero = all_zero and dF[i][iv] == Complex(0, 0) and dF_mod[i][iv] == Complex(0, 0);
    if (all_zero) continue;  // plot only if some are non-zero
    
    std::cout << jacobian_quantities[i].Target() << '\n';
    ARTSGUI::plot(f_grid, dF[i].real(), f_grid, dF[i].imag(), f_grid, dF_mod[i].real(), f_grid, dF_mod[i].imag());
  }
  std::cout << '\n';
}

int main() try {
  define_species_data();
  define_species_map();
  make_wigner_ready(20, 30, 6);
  std::cout << std::setprecision(15);
//   test_ls<LineShape::Type::DP>();
//   test_ls<LineShape::Type::LP>();
//   test_ls<LineShape::Type::VP>();
  test_ls<LineShape::Type::SDVP>();
  test_ls<LineShape::Type::HTP>();
  
  return EXIT_SUCCESS;
} catch (std::exception& e) {
  std::cerr << e.what() << '\n';
  return EXIT_FAILURE;
}
