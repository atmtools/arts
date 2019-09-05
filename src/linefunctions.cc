/* Copyright (C) 2017
 * Richard Larsson <ric.larsson@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * @file   linefunctions.h
 * @author Richard Larsson
 * @date   2017-05-16
 * 
 * @brief  Stuff related to lineshape functions.
 * 
 * This file should contain complete handling of individual lines.
 * The reason is that the old methods are cumbersome to adapt and need redesigning
 * 
 * Example usage for simple Lorentz line shape with line
 * strength scaled to the correct integration
 * set_lorentz(...)
 * apply_linestrength_scaling(...)
 */

#include "linefunctions.h"
#include <Eigen/Core>
#include <Faddeeva/Faddeeva.hh>
#include "constants.h"
#include "linescaling.h"

/** The Faddeeva function */
inline Complex w(Complex z) noexcept { return Faddeeva::w(z); }

/** The Faddeeva function partial derivative */
constexpr Complex dw(Complex z, Complex w) noexcept {
  return Complex(0, 2) * (Constant::inv_sqrt_pi - z * w);
}

/** Conversion from CGS-style lineshape to ARTS */
constexpr Complex pCqSDHC_to_arts(Complex x) noexcept {
  using Constant::c;
  using Constant::pow2;
  using Conversion::hitran2arts_linestrength;
  
  return conj(hitran2arts_linestrength(x) / pow2(c));
}

/** Conversion from CGS-style lineshape to ARTS for frequncy derivatives */
constexpr Complex pCqSDHC_to_arts_freq_deriv(Complex x) noexcept {
  using Constant::c;
  using Constant::pow4;
  using Conversion::hitran2arts_linestrength;
  
  return hitran2arts_linestrength(hitran2arts_linestrength(x)) / pow4(c);
}

/** Conversion from CGS-style lineshape to ARTS G2 derivative */
constexpr Complex pCqSDHC_to_arts_G2_deriv(Complex x) noexcept {
  using Constant::c;
  using Constant::pow4;
  using Conversion::hitran2arts_linestrength;
  
  return Complex(0, -1) *
  hitran2arts_linestrength(hitran2arts_linestrength(x)) / pow4(c);
}

/** Conversion from CGS-style lineshape to ARTS D2 derivative */
constexpr Complex pCqSDHC_to_arts_D2_deriv(Complex x) noexcept {
  using Constant::c;
  using Constant::pow4;
  using Conversion::hitran2arts_linestrength;
  
  return Complex(0, 1) * hitran2arts_linestrength(hitran2arts_linestrength(x)) /
  pow4(c);
}

void Linefunctions::set_lineshape(
    Eigen::Ref<Eigen::VectorXcd> F,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const LineRecord& line,
    const ConstVectorView vmrs,
    const Numeric& temperature,
    const Numeric& pressure,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const ArrayOfArrayOfSpeciesTag& abs_species) {
  // Pressure broadening and line mixing terms
  const auto X = line.GetShapeParams(temperature, pressure, vmrs, abs_species);

  Eigen::MatrixXcd dF(0, 0), data(F.size(), Linefunctions::ExpectedDataSize());
  const Numeric doppler_constant =
      DopplerConstant(temperature, line.IsotopologueData().Mass());

  switch (line.GetLineShapeType()) {
    case LineShape::Type::HTP:
    case LineShape::Type::SDVP:
      set_htp(F,
              dF,
              f_grid,
              zeeman_df,
              magnetic_magnitude,
              line.F(),
              doppler_constant,
              X);
      break;
    case LineShape::Type::VP:
      set_voigt(F,
                dF,
                data,
                f_grid,
                zeeman_df,
                magnetic_magnitude,
                line.F(),
                doppler_constant,
                X);
      break;
    case LineShape::Type::DP:
      set_doppler(F,
                  dF,
                  data,
                  f_grid,
                  zeeman_df,
                  magnetic_magnitude,
                  line.F(),
                  doppler_constant);
      break;
    case LineShape::Type::LP:
      set_lorentz(
          F, dF, data, f_grid, zeeman_df, magnetic_magnitude, line.F(), X);
      break;
  }

  switch (line.GetMirroringType()) {
    case MirroringType::None:
    case MirroringType::Manual:
      break;
    case MirroringType::Lorentz: {
      // Set the mirroring computational vectors and size them as needed
      Eigen::VectorXcd Fm(F.size());

      set_lorentz(Fm,
                  dF,
                  data,
                  f_grid,
                  -zeeman_df,
                  magnetic_magnitude,
                  -line.F(),
                  LineShape::mirroredOutput(X));

      // Apply mirroring;  FIXME: Add conjugate?
      F.noalias() += Fm;
    } break;
    case MirroringType::SameAsLineShape: {
      // Set the mirroring computational vectors and size them as needed
      Eigen::VectorXcd Fm(F.size());

      switch (line.GetLineShapeType()) {
        case LineShape::Type::DP:
          set_doppler(Fm,
                      dF,
                      data,
                      f_grid,
                      -zeeman_df,
                      magnetic_magnitude,
                      -line.F(),
                      -doppler_constant);
          break;
        case LineShape::Type::LP:
          set_lorentz(Fm,
                      dF,
                      data,
                      f_grid,
                      -zeeman_df,
                      magnetic_magnitude,
                      -line.F(),
                      LineShape::mirroredOutput(X));
          break;
        case LineShape::Type::VP:
          set_voigt(Fm,
                    dF,
                    data,
                    f_grid,
                    -zeeman_df,
                    magnetic_magnitude,
                    -line.F(),
                    -doppler_constant,
                    LineShape::mirroredOutput(X));
          break;
        case LineShape::Type::HTP:
        case LineShape::Type::SDVP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm,
                  dF,
                  f_grid,
                  -zeeman_df,
                  magnetic_magnitude,
                  -line.F(),
                  -doppler_constant,
                  LineShape::mirroredOutput(X));
          break;
      }

      F.noalias() += Fm;
      break;
    }
  }

  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by a type
  // that is internally interpreted to mean some kind of lineshape normalization
  switch (line.GetLineNormalizationType()) {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::VVH:
      apply_VVH_scaling(F, dF, data, f_grid, line.F(), temperature);
      break;
    case LineNormalizationType::VVW:
      apply_VVW_scaling(F, dF, f_grid, line.F());
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F, dF, f_grid, line.F(), temperature);
      break;
  }
}

void Linefunctions::set_lorentz(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& F0_noshift,
    const LineShape::Output& x,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const LineShape::Output& dT,
    const LineShape::Output& dVMR) {
  constexpr Complex cpi(0, Constant::pi);
  constexpr Complex iz(0.0, 1.0);

  // Size of the problem
  auto nppd = derivatives_data_position.nelem();

  // The central frequency
  const Numeric F0 = F0_noshift + x.D0 + zeeman_df * magnetic_magnitude + x.DV;

  // Naming req data blocks
  auto z = data.col(0);
  auto dw = data.col(1);

  z.noalias() =
      (Constant::pi * Complex(x.G0, F0) - cpi * f_grid.array()).matrix();
  F.noalias() = z.cwiseInverse();

  if (nppd) {
    dw.noalias() = -Constant::pi * F.array().square().matrix();

    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (deriv == JacPropMatType::Temperature)
        dF.col(iq).noalias() = Complex(dT.G0, dT.D0 + dT.DV) * dw;
      else if (is_frequency_parameter(deriv))
        dF.col(iq).noalias() = -iz * dw;
      else if ((deriv == JacPropMatType::LineCenter or
                is_pressure_broadening_DV(deriv)) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = iz * dw;
      else if (is_pressure_broadening_G0(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = dw;
      else if (is_pressure_broadening_D0(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = iz * dw;
      else if (is_magnetic_parameter(deriv))
        dF.col(iq).noalias() = iz * zeeman_df * dw;
      else if (deriv == JacPropMatType::VMR) {
        if (deriv.QuantumIdentity().In(quantum_identity))
          dF.col(iq).noalias() = Complex(dVMR.G0, dVMR.D0 + dVMR.DV) * dw;
        else
          dF.col(iq).setZero();
      }
    }
  }
}

void Linefunctions::set_voigt(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& F0_noshift,
    const Numeric& GD_div_F0,
    const LineShape::Output& x,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const Numeric& dGD_div_F0_dT,
    const LineShape::Output& dxdT,
    const LineShape::Output& dxdVMR) {
  constexpr Complex iz(0.0, 1.0);

  // Size of problem
  auto nppd = derivatives_data_position.nelem();

  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + x.D0 + x.DV;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0 - GD_div_F0 * (dxdT.D0 + dxdT.DV);

  // constant normalization factor for Voigt
  const Numeric fac = Constant::inv_sqrt_pi * invGD;

  // Naming req data blocks
  auto z = data.col(0);
  auto dw = data.col(1);

  // Frequency grid
  z.noalias() = invGD * (Complex(-F0, x.G0) + f_grid.array()).matrix();

  // Line shape
  F.noalias() = fac * z.unaryExpr(&w);

  if (nppd) {
    dw.noalias() = 2 * (Complex(0, fac * Constant::inv_sqrt_pi) -
                        z.cwiseProduct(F).array())
                           .matrix();

    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (is_frequency_parameter(deriv))
        dF.col(iq).noalias() = invGD * dw;
      else if (deriv == JacPropMatType::Temperature)
        dF.col(iq).noalias() =
            dw * Complex(-dxdT.D0 - dxdT.DV, dxdT.G0) * invGD -
            F * dGD_dT * invGD - dw.cwiseProduct(z) * dGD_dT * invGD;
      else if ((deriv == JacPropMatType::LineCenter or
                is_pressure_broadening_DV(deriv)) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = -F / F0 - dw * invGD - dw.cwiseProduct(z) / F0;
      else if (is_pressure_broadening_G0(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = iz * dw * invGD;
      else if (is_pressure_broadening_D0(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = -dw * invGD;
      else if (is_magnetic_parameter(deriv))
        dF.col(iq).noalias() = dw * (-zeeman_df * invGD);
      else if (deriv == JacPropMatType::VMR and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() =
            dw * Complex(-dxdVMR.D0 - dxdVMR.DV, dxdVMR.G0) * invGD;
      else
        dF.col(iq).setZero();
    }
  }
}

void Linefunctions::set_doppler(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& F0_noshift,
    const Numeric& GD_div_F0,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const Numeric& dGD_div_F0_dT) {
  auto nppd = derivatives_data_position.nelem();

  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude;
  const Numeric invGD = 1.0 / (GD_div_F0 * F0);

  // Naming data blocks
  auto x = data.col(0);
  auto mx2 = data.col(1);

  x.noalias() = (f_grid.array() - F0).matrix() * invGD;
  mx2.noalias() = -data.col(0).cwiseAbs2();
  F.noalias() = invGD * Constant::inv_sqrt_pi * mx2.array().exp().matrix();

  for (auto iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (is_frequency_parameter(deriv))
      dF.col(iq).noalias() = -2 * invGD * F.cwiseProduct(x);
    else if (deriv == JacPropMatType::Temperature)
      dF.col(iq).noalias() =
          -dGD_div_F0_dT * F0 * invGD * (2.0 * F.cwiseProduct(mx2) + F);
    else if (deriv == JacPropMatType::LineCenter and
             deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() = -F.cwiseProduct(mx2) * (2 / F0) +
                             F.cwiseProduct(x) * 2 * (invGD - 1 / F0);
    else if (deriv == JacPropMatType::VMR)
      dF.col(iq).setZero();  // Must reset incase other lineshapes are mixed in
  }
}

void Linefunctions::apply_linemixing_scaling_and_mirroring(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<Eigen::VectorXcd> Fm,
    const Eigen::Ref<Eigen::MatrixXcd> dFm,
    const LineShape::Output& X,
    const bool with_mirroring,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const LineShape::Output& dT,
    const LineShape::Output& dVMR) {
  auto nppd = derivatives_data_position.nelem();

  const Complex LM = Complex(1.0 + X.G, -X.Y);

  dF *= LM;
  if (with_mirroring) dF.noalias() += dFm * conj(LM);

  if (with_mirroring) {
    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (deriv == JacPropMatType::Temperature) {
        const auto c = Complex(dT.G, -dT.Y);
        dF.col(iq).noalias() += F * c + Fm * conj(c);
      } else if (is_pressure_broadening_G(deriv) and
                 deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = F + Fm;
      else if (is_pressure_broadening_Y(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = Complex(0, -1) * (F - Fm);
      else if (deriv == JacPropMatType::VMR and
               deriv.QuantumIdentity().In(quantum_identity)) {
        const auto c = Complex(dVMR.G, -dVMR.Y);
        dF.col(iq).noalias() += F * c + Fm * conj(c);
      }
    }
  } else {
    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (deriv == JacPropMatType::Temperature)
        dF.col(iq).noalias() += F * Complex(dT.G, -dT.Y);
      else if (is_pressure_broadening_G(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = F;
      else if (is_pressure_broadening_Y(deriv) and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = Complex(0, -1) * F;
      else if (deriv == JacPropMatType::VMR and
               deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() += F * Complex(dVMR.G, -dVMR.Y);
    }
  }

  F *= LM;
  if (with_mirroring) F.noalias() += Fm * std::conj(LM);
}

void Linefunctions::apply_rosenkranz_quadratic_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const Numeric& T,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity) {
  auto nf = f_grid.size(), nppd = derivatives_data_position.nelem();

  const Numeric invF0 = 1.0 / F0;
  const Numeric mafac =
      (Constant::h) / (2.0 * Constant::k * T) /
      std::sinh((Constant::h * F0) / (2.0 * Constant::k * T)) * invF0;

  Numeric dmafac_dT_div_fun = 0;
  if (do_temperature_jacobian(derivatives_data))
    dmafac_dT_div_fun =
        -(Constant::k * T -
          F0 * Constant::h /
              (2.0 * std::tanh(F0 * Constant::h / (2.0 * Constant::k * T)))) /
        (Constant::k * T * T);

  Numeric fun;

  for (auto iv = 0; iv < nf; iv++) {
    fun = mafac * (f_grid[iv] * f_grid[iv]);
    F[iv] *= fun;

    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      dF(iv, iq) *= fun;
      if (deriv == JacPropMatType::Temperature)
        dF(iv, iq) += dmafac_dT_div_fun * F[iv];
      else if (deriv == JacPropMatType::LineCenter and
               deriv.QuantumIdentity().In(quantum_identity))
        dF(iv, iq) +=
            (-invF0 - Constant::h / (2.0 * Constant::k * T *
                                     std::tanh(F0 * Constant::h /
                                               (2.0 * Constant::k * T)))) *
            F[iv];
      else if (is_frequency_parameter(deriv))
        dF(iv, iq) += (2.0 / f_grid[iv]) * F[iv];
    }
  }
}

void Linefunctions::apply_VVH_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const Numeric& T,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity) {
  auto nppd = derivatives_data_position.nelem();

  // 2kT is constant for the loop
  const Numeric kT = 2.0 * Constant::k * T;
  const Numeric c1 = Constant::h / kT;

  // denominator is constant for the loop
  const Numeric tanh_f0part = std::tanh(c1 * F0);
  const Numeric denom = F0 * tanh_f0part;

  // Name the data
  auto tanh_f = data.col(0);
  auto ftanh_f = data.col(1);

  tanh_f.noalias() = (c1 * f_grid.array()).tanh().matrix();
  ftanh_f.noalias() = f_grid.cwiseProduct(tanh_f) / denom;
  F.array() *= ftanh_f.array();

  dF.array().colwise() *= ftanh_f.array();
  for (auto iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (deriv == JacPropMatType::Temperature)
      dF.col(iq).noalias() +=
          -c1 / T *
          ((denom - F0 / tanh_f0part) * F - denom * ftanh_f.cwiseProduct(F) +
           f_grid.cwiseProduct(F).cwiseQuotient(tanh_f));
    else if (is_frequency_parameter(deriv))
      dF.col(iq).noalias() +=
          F.cwiseQuotient(f_grid) +
          c1 * (F.cwiseQuotient(tanh_f) - F.cwiseProduct(tanh_f));
    else if (deriv == JacPropMatType::LineCenter and
             deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() +=
          (-1.0 / F0 + c1 * tanh_f0part - c1 / tanh_f0part) * F;
  }
}

void Linefunctions::apply_VVW_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity) {
  auto nf = f_grid.size(), nppd = derivatives_data_position.nelem();

  // denominator is constant for the loop
  const Numeric invF0 = 1.0 / F0;
  const Numeric invF02 = invF0 * invF0;

  for (auto iv = 0; iv < nf; iv++) {
    // Set the factor
    const Numeric fac = f_grid[iv] * f_grid[iv] * invF02;

    // Set the line shape
    F[iv] *= fac;

    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      // The factor is applied to all partial derivatives
      dF(iv, iq) *= fac;

      // These partial derivatives are special
      if (deriv == JacPropMatType::LineCenter and
          deriv.QuantumIdentity().In(quantum_identity))
        dF(iv, iq) -= 2.0 * invF0 * F[iv];
      else if (is_frequency_parameter(deriv))
        dF(iv, iq) += 2.0 / f_grid[iv] * F[iv];
    }
  }
}

Numeric Linefunctions::lte_linestrength(Numeric S0,
                                        Numeric E0,
                                        Numeric F0,
                                        Numeric QT0,
                                        Numeric T0,
                                        Numeric QT,
                                        Numeric T) {
  const Numeric gamma = stimulated_emission(T, F0);
  const Numeric gamma_ref = stimulated_emission(T0, F0);
  const Numeric K1 = boltzman_ratio(T, T0, E0);
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
  return S0 * K1 * K2 * QT0 / QT;
}

void Linefunctions::apply_linestrength_scaling_by_lte(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const LineRecord& line,
    const Numeric& T,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const Numeric& dQT_dT) {
  auto nppd = derivatives_data_position.nelem();

  const Numeric gamma = stimulated_emission(T, line.F());
  const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
  const Numeric K1 = boltzman_ratio(T, line.Ti0(), line.Elow());
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);

  const Numeric invQT = 1.0 / QT;
  const Numeric S = line.I0() * isotopic_ratio * QT0 * invQT * K1 * K2;

  F *= S;
  dF *= S;
  for (auto iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (deriv == JacPropMatType::Temperature)
      dF.col(iq).noalias() +=
          F * (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T) /
                   K2 +
               dboltzman_ratio_dT(K1, T, line.Elow()) / K1 - invQT * dQT_dT);
    else if (deriv == JacPropMatType::LineStrength and
             deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() = F / line.I0();  //nb. overwrite
    else if (deriv == JacPropMatType::LineCenter and
             deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() +=
          F *
          dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0()) /
          K2;
  }

  // No NLTE variables
  N.setZero();
  dN.setZero();
}

void Linefunctions::apply_linestrength_scaling_by_vibrational_nlte(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const LineRecord& line,
    const Numeric& T,
    const Numeric& Tu,
    const Numeric& Tl,
    const Numeric& Evu,
    const Numeric& Evl,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const Numeric& dQT_dT) {
  auto nppd = derivatives_data_position.nelem();

  const Numeric gamma = stimulated_emission(T, line.F());
  const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
  const Numeric r_low = boltzman_ratio(Tl, T, Evl);
  const Numeric r_upp = boltzman_ratio(Tu, T, Evu);

  const Numeric K1 = boltzman_ratio(T, line.Ti0(), line.Elow());
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
  const Numeric K3 = absorption_nlte_ratio(gamma, r_upp, r_low);
  const Numeric K4 = r_upp;

  const Numeric invQT = 1.0 / QT;
  const Numeric QT_ratio = QT0 * invQT;

  const Numeric dS_dS0_abs = isotopic_ratio * QT_ratio * K1 * K2 * K3;
  const Numeric S_abs = line.I0() * dS_dS0_abs;
  const Numeric dS_dS0_src = isotopic_ratio * QT_ratio * K1 * K2 * K4;
  const Numeric S_src = line.I0() * dS_dS0_src;

  dN.noalias() = dF * (S_src - S_abs);
  dF *= S_abs;
  for (auto iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (deriv == JacPropMatType::Temperature) {
      const Numeric dS_dT_abs =
          S_abs *
          (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T) /
               K2 +
           dboltzman_ratio_dT(K1, T, line.Elow()) / K1 +
           dabsorption_nlte_rate_dT(gamma, T, line.F(), Evl, Evu, K4, r_low) /
               K3 -
           invQT * dQT_dT);
      const Numeric dS_dT_src =
          S_src *
          (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T) /
               K2 +
           dboltzman_ratio_dT(K1, T, line.Elow()) / K1 -
           dboltzman_ratio_dT(K4, T, Evu) / K4 - invQT * dQT_dT);

      dN.col(iq).noalias() += F * (dS_dT_src - dS_dT_abs);
      dF.col(iq).noalias() += F * dS_dT_abs;
    } else if (deriv == JacPropMatType::LineStrength and
               deriv.QuantumIdentity().In(quantum_identity)) {
      dF.col(iq).noalias() = F * dS_dS0_abs;
      dN.col(iq).noalias() = F * (dS_dS0_src - dS_dS0_abs);
    } else if (deriv == JacPropMatType::LineCenter and
               deriv.QuantumIdentity().In(quantum_identity)) {
      const Numeric dS_dF0_abs =
          S_abs *
          (dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0()) /
               K2 +
           dabsorption_nlte_rate_dF0(gamma, T, K4, r_low) / K3);
      const Numeric dS_dF0_src =
          S_src *
          dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0()) /
          K2;

      dN.col(iq).noalias() += F * (dS_dF0_src - dS_dF0_abs);
      dF.col(iq).noalias() += F * dS_dF0_abs;
    } else if (deriv == JacPropMatType::NLTE) {
      if (deriv.QuantumIdentity().InLower(quantum_identity)) {
        const Numeric dS_dTl_abs =
            S_abs * dabsorption_nlte_rate_dTl(gamma, T, Tl, Evl, r_low) / K3;

        dN.col(iq).noalias() = -F * dS_dTl_abs;
        dF.col(iq).noalias() = -dN.col(iq);
      } else if (deriv.QuantumIdentity().InUpper(quantum_identity)) {
        const Numeric dS_dTu_abs =
            S_abs * dabsorption_nlte_rate_dTu(gamma, T, Tu, Evu, K4) / K3;
        const Numeric dS_dTu_src = S_src * dboltzman_ratio_dT(K4, Tu, Evu) / K4;

        dN.col(iq).noalias() = F * (dS_dTu_src - dS_dTu_abs);
        dF.col(iq).noalias() = F * dS_dTu_abs;
      }
    }
  }

  N.noalias() = F * (S_src - S_abs);
  F *= S_abs;
}

void Linefunctions::apply_linestrength_from_full_linemixing(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Numeric& F0,
    const Numeric& T,
    const Complex& S_LM,
    const Numeric& isotopic_ratio,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const Complex& dS_LM_dT) {
  const Index nppd = derivatives_data_position.nelem();

  constexpr Numeric C1 = -Constant::h / Constant::k;

  const Numeric invT = 1.0 / T, F0_invT = F0 * invT,
                exp_factor = std::exp(C1 * F0_invT),
                f0_factor = F0 * (1.0 - exp_factor);

  const Complex S = S_LM * f0_factor * isotopic_ratio;

  for (Index iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (deriv == JacPropMatType::Temperature) {
      dF.col(iq) *= S_LM;
      dF.col(iq).noalias() +=
          F *
          (dS_LM_dT * f0_factor + S_LM * C1 * F0_invT * F0_invT * exp_factor) *
          isotopic_ratio;
    } else if (deriv == JacPropMatType::LineStrength) {
      if (deriv.QuantumIdentity().In(quantum_identity))
        throw std::runtime_error("Not working yet");
    } else if (deriv == JacPropMatType::LineCenter or
               is_pressure_broadening_DV(deriv)) {
      if (deriv.QuantumIdentity().In(quantum_identity))
        throw std::runtime_error("Not working yet");
    } else {
      dF.col(iq) *= S;
    }
  }

  F *= S;
}

void Linefunctions::apply_dipole(Eigen::Ref<Eigen::VectorXcd> F,
                                 Eigen::Ref<Eigen::MatrixXcd> /*dF*/,
                                 const Numeric& F0,
                                 const Numeric& T,
                                 const Numeric& d0,
                                 const Numeric& rho,
                                 const Numeric& isotopic_ratio,
                                 const ArrayOfRetrievalQuantity&,
                                 const ArrayOfIndex& derivatives_data_position,
                                 const QuantumIdentifier&,
                                 const Numeric&) {
  // Output is d0^2 * rho * F * isotopic_ratio * F0 * (1-e^(hF0/kT))

  constexpr Numeric C1 = -Constant::h / Constant::k;

  const Index nppd = derivatives_data_position.nelem();

  const Numeric S = d0 * d0 * rho * isotopic_ratio, invT = 1.0 / T,
                F0_invT = F0 * invT, exp_factor = exp(C1 * F0_invT),
                f0_factor = F0 * (1.0 - exp_factor);

  if (nppd)
    throw std::runtime_error(
        "Cannot support Jacobian from dipole calculations yet");

  F *= S * f0_factor;
}

void Linefunctions::apply_linefunctiondata_jacobian_scaling(
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const LineRecord& line,
    const Numeric& T,
    const Numeric& P,
    const Vector& vmrs) {
  auto nppd = derivatives_data_position.nelem();

  for (auto iq = 0; iq < nppd; iq++) {
    const RetrievalQuantity& rt =
        derivatives_data[derivatives_data_position[iq]];
    if (is_lineshape_parameter(rt) and
        rt.QuantumIdentity().In(quantum_identity))
      dF.col(iq) *= line.GetPrepInternalDerivative(T, P, vmrs, rt);
  }
}

Numeric Linefunctions::DopplerConstant(Numeric T, Numeric mass) {
  return std::sqrt(Constant::doppler_broadening_const_squared * T / mass);
}

Numeric Linefunctions::dDopplerConstant_dT(const Numeric& T,
                                           const Numeric& dc) {
  return dc / (2 * T);
}

void Linefunctions::set_cross_section_for_single_line(
    Eigen::Ref<Eigen::VectorXcd> F_full,
    Eigen::Ref<Eigen::MatrixXcd> dF_full,
    Eigen::Ref<Eigen::VectorXcd> N_full,
    Eigen::Ref<Eigen::MatrixXcd> dN_full,
    Eigen::Ref<Eigen::MatrixXcd> data_block_full,
    Index& start_cutoff,
    Index& nelem_cutoff,
    const Eigen::Ref<const Eigen::VectorXd> f_grid_full,
    const LineRecord& line,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const Vector& volume_mixing_ratio_of_lineshape,
    const ConstVectorView nlte_distribution,
    const Numeric& pressure,
    const Numeric& temperature,
    const Numeric& doppler_constant,
    const Numeric& isotopologue_ratio,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& ddoppler_constant_dT,
    const Numeric& partition_function_at_temperature,
    const Numeric& dpartition_function_at_temperature_dT,
    const Numeric& partition_function_at_line_temperature,
    const bool cutoff_call) {
  /* Single line shape solver
     
     The equation being solved here:
     F = LM NORM K1 K2 K3 QT/QT0 S0 f(...) - CUT(F),
     N = LM NORM K1 K2 (K4/K3 - 1) QT/QT0 S0 f(...) - CUT(N),
     dF = d(F)/dx
     dN = d(N)/dx
     
     or if using external population distribution calculations
     
     F = A f(...) - CUT(F)
     N = r f(...) = CUT(N)
     
     Where 
     LM:  Line mixing
     NORM: Normalization
     K1: Boltzmann ratio
     K2: Stimulated emission ratio
     K3: NLTE absorption ratio
     K4: NLTE emission ratio
     QT: Partition function
     QT0: Partition function at reference temperature
     S0: Line strength at reference temperature
     f(...): Line shape, including mirroring
     CUT(X): Computation of F and N at some provided cutoff frequency
     A is the absorption from external level distributions
     r is the emission ratio addition to Planck from external level distributions
  */

  /* Get the cutoff range if applicable
   * The line will have had a number set either by default or by finding the LSM CUT flag
   * followed by the numeric that is used to set the cutoff frequency.  By default, the 
   * cutoff frequency is set to a negative number, and it will only be used if there is
   * a non-negative number in the frequency range.
   */

  // Size and type of problem
  const Numeric& cutoff = line.CutOff();
  if (not cutoff_call)
    find_cutoff_ranges(
        start_cutoff, nelem_cutoff, f_grid_full, line.F(), cutoff);
  else {
    start_cutoff = 0;
    nelem_cutoff = 1;
  }

  const bool need_cutoff = cutoff > 0 and not cutoff_call;
  const bool do_temperature = do_temperature_jacobian(derivatives_data);

  // Leave this function if there is nothing to compute
  if (nelem_cutoff == 0) return;

  // Extract the quantum identify of the line to be used in-case there are derivatives
  const QuantumIdentifier& QI = line.QuantumIdentity();

  // Pressure broadening and line mixing terms
  const auto X = line.GetPrepShapeParams(
      temperature, pressure, volume_mixing_ratio_of_lineshape);

  constexpr LineShape::Output empty_output = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Partial derivatives for temperature
  const auto dXdT =
      do_temperature
          ? line.GetPrepShapeParams_dT(
                temperature, pressure, volume_mixing_ratio_of_lineshape)
          : empty_output;

  // Partial derivatives for VMR... the first function
  auto do_vmr = do_vmr_jacobian(derivatives_data,
                                line.QuantumIdentity());  // At all, Species
  const auto dXdVMR =
      do_vmr.test ? line.GetShapeParams_dVMR(temperature, pressure, do_vmr.qid)
                  : empty_output;

  // Arrays on which all computations happen are segments of the full input, and the segmenting is part of the output
  auto F = F_full.segment(start_cutoff, nelem_cutoff);
  auto N = N_full.segment(start_cutoff, nelem_cutoff);
  auto dF = dF_full.middleRows(start_cutoff, nelem_cutoff);
  auto dN = dN_full.middleRows(start_cutoff, nelem_cutoff);
  auto data = data_block_full.middleRows(start_cutoff, nelem_cutoff);
  auto f_grid = f_grid_full.middleRows(start_cutoff, nelem_cutoff);

  switch (line.GetLineShapeType()) {
    case LineShape::Type::DP:
      set_doppler(F,
                  dF,
                  data,
                  f_grid,
                  zeeman_df,
                  magnetic_magnitude,
                  line.F(),
                  doppler_constant,
                  derivatives_data,
                  derivatives_data_position,
                  QI,
                  ddoppler_constant_dT);
      break;
    case LineShape::Type::HTP:
    case LineShape::Type::SDVP:
      set_htp(F,
              dF,
              f_grid,
              zeeman_df,
              magnetic_magnitude,
              line.F(),
              doppler_constant,
              X,
              derivatives_data,
              derivatives_data_position,
              QI,
              ddoppler_constant_dT,
              dXdT,
              dXdVMR);
      break;
    case LineShape::Type::LP:
      set_lorentz(F,
                  dF,
                  data,
                  f_grid,
                  zeeman_df,
                  magnetic_magnitude,
                  line.F(),
                  X,
                  derivatives_data,
                  derivatives_data_position,
                  QI,
                  dXdT,
                  dXdVMR);
      break;
    case LineShape::Type::VP:
      set_voigt(F,
                dF,
                data,
                f_grid,
                zeeman_df,
                magnetic_magnitude,
                line.F(),
                doppler_constant,
                X,
                derivatives_data,
                derivatives_data_position,
                QI,
                ddoppler_constant_dT,
                dXdT,
                dXdVMR);
      break;
  }

  // Set the mirroring by repeating computations above using
  // negative numbers for frequency of line related terms
  // The user sets if they want mirroring by LSM MTM followed by an index
  // that is interpreted as either mirroring by the same line shape or as
  // mirroring by Lorentz lineshape
  const bool with_mirroring =
      line.GetMirroringType() not_eq MirroringType::None and
      line.GetMirroringType() not_eq MirroringType::Manual;
  switch (line.GetMirroringType()) {
    case MirroringType::None:
    case MirroringType::Manual:
      break;
    case MirroringType::Lorentz:
      set_lorentz(N,
                  dN,
                  data,
                  f_grid,
                  -zeeman_df,
                  magnetic_magnitude,
                  -line.F(),
                  LineShape::mirroredOutput(X),
                  derivatives_data,
                  derivatives_data_position,
                  QI,
                  LineShape::mirroredOutput(dXdT),
                  LineShape::mirroredOutput(dXdVMR));
      break;
    case MirroringType::SameAsLineShape:
      switch (line.GetLineShapeType()) {
        case LineShape::Type::DP:
          set_doppler(N,
                      dN,
                      data,
                      f_grid,
                      -zeeman_df,
                      magnetic_magnitude,
                      -line.F(),
                      -doppler_constant,
                      derivatives_data,
                      derivatives_data_position,
                      QI,
                      -ddoppler_constant_dT);
          break;
        case LineShape::Type::LP:
          set_lorentz(N,
                      dN,
                      data,
                      f_grid,
                      -zeeman_df,
                      magnetic_magnitude,
                      -line.F(),
                      LineShape::mirroredOutput(X),
                      derivatives_data,
                      derivatives_data_position,
                      QI,
                      LineShape::mirroredOutput(dXdT),
                      LineShape::mirroredOutput(dXdVMR));
          break;
        case LineShape::Type::VP:
          set_voigt(N,
                    dN,
                    data,
                    f_grid,
                    -zeeman_df,
                    magnetic_magnitude,
                    -line.F(),
                    -doppler_constant,
                    LineShape::mirroredOutput(X),
                    derivatives_data,
                    derivatives_data_position,
                    QI,
                    -ddoppler_constant_dT,
                    LineShape::mirroredOutput(dXdT),
                    LineShape::mirroredOutput(dXdVMR));
          break;
        case LineShape::Type::HTP:
        case LineShape::Type::SDVP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(N,
                  dN,
                  f_grid,
                  -zeeman_df,
                  magnetic_magnitude,
                  -line.F(),
                  -doppler_constant,
                  LineShape::mirroredOutput(X),
                  derivatives_data,
                  derivatives_data_position,
                  QI,
                  -ddoppler_constant_dT,
                  LineShape::mirroredOutput(dXdT),
                  LineShape::mirroredOutput(dXdVMR));
          break;
      }
      break;
  }

  // Mixing and mirroring can only apply to non-Doppler shapes
  if (line.GetLineShapeType() not_eq LineShape::Type::DP) {
    apply_linemixing_scaling_and_mirroring(F,
                                           dF,
                                           N,
                                           dN,
                                           X,
                                           with_mirroring,
                                           derivatives_data,
                                           derivatives_data_position,
                                           QI,
                                           dXdT,
                                           dXdVMR);

    // Apply line mixing and pressure broadening partial derivatives
    apply_linefunctiondata_jacobian_scaling(dF,
                                            derivatives_data,
                                            derivatives_data_position,
                                            QI,
                                            line,
                                            temperature,
                                            pressure,
                                            volume_mixing_ratio_of_lineshape);
  }

  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by an index
  // that is internally interpreted to mean some kind of lineshape normalization
  switch (line.GetLineNormalizationType()) {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::VVH:
      apply_VVH_scaling(F,
                        dF,
                        data,
                        f_grid,
                        line.F(),
                        temperature,
                        derivatives_data,
                        derivatives_data_position,
                        QI);
      break;
    case LineNormalizationType::VVW:
      apply_VVW_scaling(F,
                        dF,
                        f_grid,
                        line.F(),
                        derivatives_data,
                        derivatives_data_position,
                        QI);
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F,
                                         dF,
                                         f_grid,
                                         line.F(),
                                         temperature,
                                         derivatives_data,
                                         derivatives_data_position,
                                         QI);
      break;
  }

  // Apply line strength by whatever method is necessary
  switch (line.GetLinePopulationType()) {
    case LinePopulationType::ByLTE:
      apply_linestrength_scaling_by_lte(F,
                                        dF,
                                        N,
                                        dN,
                                        line,
                                        temperature,
                                        isotopologue_ratio,
                                        partition_function_at_temperature,
                                        partition_function_at_line_temperature,
                                        derivatives_data,
                                        derivatives_data_position,
                                        QI,
                                        dpartition_function_at_temperature_dT);
      break;
    case LinePopulationType::ByVibrationalTemperatures:
      if (line.NLTELowerIndex() >= 0 and
          nlte_distribution.nelem() <= line.NLTELowerIndex())
        throw std::runtime_error(
            "Bad lower vibrational temperature record.  It contains fewer temperatures than LineRecord require");

      if (line.NLTEUpperIndex() >= 0 and
          nlte_distribution.nelem() <= line.NLTEUpperIndex())
        throw std::runtime_error(
            "Bad upper vibrational temperature record.  It contains fewer temperatures than LineRecord require");

      apply_linestrength_scaling_by_vibrational_nlte(
          F,
          dF,
          N,
          dN,
          line,
          temperature,
          line.NLTEUpperIndex() < 0 ? temperature
                                    : nlte_distribution[line.NLTEUpperIndex()],
          line.NLTELowerIndex() < 0 ? temperature
                                    : nlte_distribution[line.NLTELowerIndex()],
          line.Evupp() < 0 ? 0 : line.Evupp(),
          line.Evlow() < 0 ? 0 : line.Evlow(),
          isotopologue_ratio,
          partition_function_at_temperature,
          partition_function_at_line_temperature,
          derivatives_data,
          derivatives_data_position,
          QI,
          dpartition_function_at_temperature_dT);
      break;
    case LinePopulationType::ByPopulationDistribution:
      if (line.NLTELowerIndex() < 0)
        throw std::runtime_error(
            "No lower NLTE distribution data for line marked to require it");
      if (nlte_distribution.nelem() <= line.NLTELowerIndex())
        throw std::runtime_error(
            "Bad lower NLTE distribution record.  It contains fewer values than LineRecord require");
      if (nlte_distribution[line.NLTELowerIndex()] <= 0)
        throw std::runtime_error(
            "Bad lower NLTE distribution number.  It should be strictly above 0.");

      if (line.NLTEUpperIndex() < 0)
        throw std::runtime_error(
            "No upper NLTE distribution data for line marked to require it");
      if (nlte_distribution.nelem() <= line.NLTEUpperIndex())
        throw std::runtime_error(
            "Bad upper NLTE distribution record.  It contains fewer values than LineRecord require");
      if (nlte_distribution[line.NLTEUpperIndex()] <= 0)
        throw std::runtime_error(
            "Bad upper NLTE distribution number.  It should be strictly above 0.");

      apply_linestrength_from_nlte_level_distributions(
          F,
          dF,
          N,
          dN,
          nlte_distribution[line.NLTELowerIndex()],
          nlte_distribution[line.NLTEUpperIndex()],
          line.G_lower(),
          line.G_upper(),
          line.A(),
          line.F(),
          temperature,
          derivatives_data,
          derivatives_data_position,
          QI);
      break;
  }

  // Cutoff frequency is applied at the end because
  // the entire process above is complicated and applying
  // cutoff last means that the code is kept cleaner

  if (need_cutoff)
    apply_cutoff(F,
                 dF,
                 N,
                 dN,
                 derivatives_data,
                 derivatives_data_position,
                 line,
                 volume_mixing_ratio_of_lineshape,
                 nlte_distribution,
                 pressure,
                 temperature,
                 doppler_constant,
                 isotopologue_ratio,
                 zeeman_df,
                 magnetic_magnitude,
                 ddoppler_constant_dT,
                 partition_function_at_temperature,
                 dpartition_function_at_temperature_dT,
                 partition_function_at_line_temperature);
}

void Linefunctions::apply_cutoff(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const LineRecord& line,
    const Vector& volume_mixing_ratio_of_lineshape,
    const ConstVectorView nlte_distribution,
    const Numeric& pressure,
    const Numeric& temperature,
    const Numeric& doppler_constant,
    const Numeric& isotopologue_ratio,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& ddoppler_constant_dT,
    const Numeric& partition_function_at_temperature,
    const Numeric& dpartition_function_at_temperature_dT,
    const Numeric& partition_function_at_line_temperature) {
  // Size of derivatives
  auto nj = dF.cols();

  // Setup compute variables

  const auto v = Vector(
      1, (line.F() > 0) ? line.F() + line.CutOff() : line.F() - line.CutOff());
  const auto f_grid_cutoff = MapToEigen(v);
  Eigen::VectorXcd Fc(1), Nc(1);
  Eigen::RowVectorXcd dFc(nj), dNc(nj), data(Linefunctions::ExpectedDataSize());
  Index _tmp1, _tmp2;

  // Recompute the line for a single frequency
  set_cross_section_for_single_line(Fc,
                                    dFc,
                                    Nc,
                                    dNc,
                                    data,
                                    _tmp1,
                                    _tmp2,
                                    f_grid_cutoff,
                                    line,
                                    derivatives_data,
                                    derivatives_data_position,
                                    volume_mixing_ratio_of_lineshape,
                                    nlte_distribution,
                                    pressure,
                                    temperature,
                                    doppler_constant,
                                    isotopologue_ratio,
                                    zeeman_df,
                                    magnetic_magnitude,
                                    ddoppler_constant_dT,
                                    partition_function_at_temperature,
                                    dpartition_function_at_temperature_dT,
                                    partition_function_at_line_temperature,
                                    true);

  // Apply cutoff values
  F.array() -= Fc[0];
  N.array() -= Nc[0];
  for (Index i = 0; i < nj; i++) dF.col(i).array() -= dFc[i];
  for (Index i = 0; i < nj; i++) dN.col(i).array() -= dNc[i];
}

void Linefunctions::find_cutoff_ranges(
    Index& start_cutoff,
    Index& nelem_cutoff,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const Numeric& cutoff) {
  auto nf = f_grid.size();

  const bool need_cutoff = (cutoff > 0);
  if (need_cutoff) {
    // Find range of simulations
    start_cutoff = 0;
    Index i_f_max = nf - 1;

    // Loop over positions to compute the line shape cutoff point
    while (start_cutoff < nf and (F0 - cutoff) > f_grid[start_cutoff])
      ++start_cutoff;
    while (i_f_max >= start_cutoff and (F0 + cutoff) < f_grid[i_f_max])
      --i_f_max;

    //  The extent is one more than the difference between the indices of interest
    nelem_cutoff = i_f_max - start_cutoff + 1;  // min is 0, max is nf
  } else {
    start_cutoff = 0;
    nelem_cutoff = nf;
  }
}

void Linefunctions::apply_linestrength_from_nlte_level_distributions(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const Numeric& r1,
    const Numeric& r2,
    const Numeric& g1,
    const Numeric& g2,
    const Numeric& A21,
    const Numeric& F0,
    const Numeric& T,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity) {
  // Size of the problem
  auto nppd = derivatives_data_position.nelem();

  // Physical constants
  constexpr Numeric c0 = 2.0 * Constant::h / Constant::pow2(Constant::c);
  constexpr Numeric c1 = Constant::h / (4 * Constant::pi);

  // Constants based on input
  const Numeric c2 = c0 * F0 * F0 * F0;
  const Numeric c3 = c1 * F0;
  const Numeric x = g2 / g1;

  /*
    Einstein 'active' coefficients are as follows:
    const Numeric B21 =  A21 / c2;
    const Numeric B12 = x * B21;
  */

  // Absorption strength
  const Numeric k = c3 * (r1 * x - r2) * (A21 / c2);

  // Emission strength
  const Numeric e = c3 * r2 * A21;

  // Planck function of this line
  const Numeric exp_T = std::exp(Constant::h / Constant::k * F0 / T),
                b = c2 / (exp_T - 1);

  // Ratio between emission and absorption constant
  const Numeric ratio = e / b - k;

  // Constants ALMOST everywhere inside these loops
  dN.noalias() = dF * ratio;
  dF *= k;

  // Partial derivatives
  for (auto iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (deriv == JacPropMatType::Temperature)
      dN.col(iq).noalias() +=
          F * e * Constant::h * F0 * exp_T / (c2 * Constant::k * T * T);
    else if (deriv == JacPropMatType::LineCenter and
             deriv.QuantumIdentity().In(quantum_identity)) {
      Numeric done_over_b_df0 =
          Constant::h * exp_T / (c2 * Constant::k * T) - 3.0 * b / F0;
      Numeric de_df0 = c1 * r2 * A21;
      Numeric dk_df0 = c1 * (r1 * x - r2) * (A21 / c2) - 3.0 * k / F0;

      dN.col(iq).noalias() += F * (e * done_over_b_df0 + de_df0 / b - dk_df0);
      dF.col(iq).noalias() += F * dk_df0;
    } else if (deriv == JacPropMatType::NLTE) {
      if (deriv.QuantumIdentity().InLower(quantum_identity)) {
        const Numeric dk_dr2 = -c3 * A21 / c2, de_dr2 = c3 * A21,
                      dratio_dr2 = de_dr2 / b - dk_dr2;
        dN.col(iq).noalias() = F * dratio_dr2;
        dF.col(iq).noalias() += F * dk_dr2;
      } else if (deriv.QuantumIdentity().InUpper(quantum_identity)) {
        const Numeric dk_dr1 = c3 * x * A21 / c2;
        dF.col(iq).noalias() += F * dk_dr1;
      }
    }
  }

  // Set source function to be the relative amount emitted by the line divided by the Planck function
  N.noalias() = F * ratio;

  // Set absorption
  F *= k;
}

void Linefunctions::set_htp(Eigen::Ref<Eigen::VectorXcd> F,
                            Eigen::Ref<Eigen::MatrixXcd> dF,
                            const Eigen::Ref<const Eigen::VectorXd> f_grid,
                            const Numeric& zeeman_df_si,
                            const Numeric& magnetic_magnitude_si,
                            const Numeric& F0_noshift_si,
                            const Numeric& GD_div_F0_si,
                            const LineShape::Output& x_si,
                            const ArrayOfRetrievalQuantity& derivatives_data,
                            const ArrayOfIndex& derivatives_data_position,
                            const QuantumIdentifier& quantum_identity,
                            const Numeric& dGD_div_F0_dT_si,
                            const LineShape::Output& dxdT_si,
                            const LineShape::Output& dxdVMR_si) {
  using Constant::inv_sqrt_pi;
  using Constant::pi;
  using Constant::pow2;
  using Constant::pow3;
  using Constant::pow4;
  using Constant::sqrt_ln_2;
  using Constant::sqrt_pi;
  using Conversion::freq2kaycm;
  using std::abs;
  using std::imag;
  using std::real;
  using std::sqrt;

  // Convert to CGS for using original function
  const Numeric sg0 =
      freq2kaycm(F0_noshift_si + zeeman_df_si * magnetic_magnitude_si);
  const Numeric GamD = GD_div_F0_si * sg0 / sqrt_ln_2;
  const auto x = si2cgs(x_si);
  const auto dT = si2cgs(dxdT_si);
  const auto dV = si2cgs(dxdVMR_si);

  // General normalization
  const Numeric cte = sqrt_ln_2 / GamD;
  constexpr Complex iz(0, 1);

  // Calculating the different parameters
  const Complex c0(x.G0, -x.D0);
  const Complex c2(x.G2, -x.D2);
  const Complex c0t = (1 - x.ETA) * (c0 - 1.5 * c2) + x.FVC;
  const Complex c2t = (1 - x.ETA) * c2;
  const Complex Y = pow2(1 / (2 * cte * c2t));
  const Complex sqrtY = sqrt(Y);

  // For all frequencies
  for (auto iv = 0; iv < f_grid.size(); iv++) {
    const Complex X = (iz * (freq2kaycm(f_grid[iv]) - sg0) + c0t) / c2t;
    const Complex sqrtXY = sqrt(X + Y);
    const Complex sqrtX = sqrt(X);

    // Declare terms required for the derivatives as external to the if-statement
    Complex Z1, Z2, Zb, W1, W2, Wb;
    Complex Aterm, Bterm;
    if (abs(c2t) ==
        0) {  // If this method does not require G2, D2, or ETA==1, then FVC matters.
      Z1 = (iz * (freq2kaycm(f_grid[iv]) - sg0) + c0t) * cte;
      W1 = w(iz * Z1);

      Aterm = sqrt_pi * cte * W1;
      if (abs(Z1) <=
          4e3)  // For very large Z1 (i.e., very large broadening or very large distance from line-center
        Bterm = sqrt_pi * cte * ((1 - pow2(Z1)) * W1 + Z1 * inv_sqrt_pi);
      else
        Bterm = cte * (sqrt_pi * W1 + 0.5 / Z1 - 0.75 / pow3(Z1));
    } else if (
        abs(X) <=
        3e-8 *
            abs(Y)) {  // If this method is executed very close to the line center
      Z1 = (iz * (freq2kaycm(f_grid[iv]) - sg0) + c0t) * cte;
      Z2 = sqrtXY + sqrtY;

      W1 = w(iz * Z1);
      W2 = w(iz * Z2);

      Aterm = sqrt_pi * cte * (W1 - W2);
      Bterm = (-1 + sqrt_pi / (2 * sqrtY) * (1 - pow2(Z1)) * W1 -
               sqrt_pi / (2 * sqrtY) * (1 - pow2(Z2)) * W2) /
              c2t;
    } else if (
        abs(Y) <=
        1e-15 *
            abs(X)) {  // If this method is executed very far from the line center
      Z1 = sqrtXY;

      W1 = w(iz * Z1);

      if (abs(sqrtX) <= 4e3) {  // If X is small still
        Zb = sqrtX;
        Wb = w(iz * Zb);

        Aterm = (2 * sqrt_pi / c2t) * (inv_sqrt_pi - sqrtX * Wb);
        Bterm =
            (1 / c2t) *
            (-1 + 2 * sqrt_pi * (1 - X - 2 * Y) * (inv_sqrt_pi - sqrtX * Wb) +
             2 * sqrt_pi * sqrtXY * W1);
      } else {
        Aterm = (1 / c2t) * (1 / X - 1.5 / pow2(X));
        Bterm = (1 / c2t) * (-1 + (1 - X - 2 * Y) * (1 / X - 1.5 / pow2(X)) +
                             2 * sqrt_pi * sqrtXY * W1);
      }
    } else {  // General calculations
      Z1 = sqrtXY - sqrtY;
      Z2 = Z1 + 2 * sqrtY;

      // NOTE: the region of w might matter according to original code!  So this might need changing...
      W1 = w(iz * Z1);
      W2 = w(iz * Z2);

      Aterm = sqrt_pi * cte * (W1 - W2);
      Bterm = (-1 + sqrt_pi / (2 * sqrtY) * (1 - pow2(Z1)) * W1 -
               sqrt_pi / (2 * sqrtY) * (1 - pow2(Z2)) * W2) /
              c2t;
    }

    F(iv) = Aterm / (pi * (((c0 - 1.5 * c2) * x.ETA - x.FVC) * Aterm +
                           Bterm * c2 * x.ETA + 1));

    for (auto iq = 0; iq < derivatives_data_position.nelem(); iq++) {
      const RetrievalQuantity& rt =
          derivatives_data[derivatives_data_position[iq]];

      Numeric dcte = 0;
      Complex dc0 = 0, dc2 = 0, dc0t = 0, dc2t = 0;
      if (rt == JacPropMatType::Temperature) {
        dcte = (-(dGD_div_F0_dT_si * sg0 / Constant::sqrt_ln_2) / GamD) *
               cte;                    // for T
        dc0 = Complex(dT.G0, -dT.D0);  // for T
        dc2 = Complex(dT.G2, -dT.D2);  // for T
        dc0t = (-(c0 - 1.5 * c2) * dT.ETA + dT.FVC) +
               ((1 - x.ETA) * (dc0 - 1.5 * dc2));     // for T
        dc2t = (-c2 * dT.ETA) + ((1 - x.ETA) * dc2);  // for T
      } else if (rt == JacPropMatType::VMR and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc0 = Complex(dV.G0, -dV.D0);  // for VMR
        dc2 = Complex(dV.G2, -dV.D2);  // for VMR
        dc0t = (-(c0 - 1.5 * c2) * dV.ETA + dV.FVC) +
               ((1 - x.ETA) * (dc0 - 1.5 * dc2));     // for VMR
        dc2t = (-c2 * dV.ETA) + ((1 - x.ETA) * dc2);  // for VMR
      } else if (is_pressure_broadening_G2(rt) and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc2 =
            1;  // for Gam2(T, VMR) and for Shift2(T, VMR), the derivative is wrt c2 for later computations
        dc0t = ((1 - x.ETA) * (dc0 - 1.5 * dc2));
        dc2t = ((1 - x.ETA) * dc2);
      } else if (is_pressure_broadening_D2(rt) and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc2 =
            -iz;  // for Gam2(T, VMR) and for Shift2(T, VMR), the derivative is wrt c2 for later computations
        dc0t = ((1 - x.ETA) * (dc0 - 1.5 * dc2));
        dc2t = ((1 - x.ETA) * dc2);
      } else if (is_pressure_broadening_G0(rt) and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc0 =
            1;  // for Gam0(T, VMR) and for Shift0(T, VMR), the derivative is wrt c0 for later computations
        dc0t = ((1 - x.ETA) * (dc0 - 1.5 * dc2));
      } else if (is_pressure_broadening_D0(rt) and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc0 =
            -iz;  // for Gam0(T, VMR) and for Shift0(T, VMR), the derivative is wrt c0 for later computations
        dc0t = ((1 - x.ETA) * (dc0 - 1.5 * dc2));
      } else if (is_pressure_broadening_FVC(rt) and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc0t = (1);  // for FVC(T, VMR)
      } else if (is_pressure_broadening_ETA(rt) and
                 rt.QuantumIdentity().In(quantum_identity)) {
        dc0t = (-c0 + 1.5 * c2);  // for eta(T, VMR)
        dc2t = (-c2);             // for eta(T, VMR)
      }

      const Complex dY = (-2 * dcte / cte - 2 * dc2t / c2t) * Y;  // for all

      Complex dX;
      if (is_magnetic_parameter(rt))
        dX = -iz * freq2kaycm(zeeman_df_si) / c2t;  // for H
      else if (rt == JacPropMatType::LineCenter and
               rt.QuantumIdentity().In(quantum_identity))
        dX = -iz / c2t;  // for sg0
      else if (is_frequency_parameter(rt))
        dX = iz / c2t;  // for sg
      else
        dX =
            (-(iz * (freq2kaycm(f_grid[iv]) - sg0) + c0t) * dc2t + c2t * dc0t) /
            pow2(c2t);  // for c0t and c2t

      Complex dAterm, dBterm;
      if (abs(c2t) == 0) {
        Complex dZ1;
        if (is_magnetic_parameter(rt))
          dZ1 = -iz * cte * freq2kaycm(zeeman_df_si);  // for H
        else if (rt == JacPropMatType::LineCenter and
                 rt.QuantumIdentity().In(quantum_identity))
          dZ1 = -iz * cte;  // for sg0
        else if (is_frequency_parameter(rt))
          dZ1 = iz * cte;  // for sg
        else
          dZ1 = (iz * (freq2kaycm(f_grid[iv]) - sg0) + c0t) * dcte +
                cte * dc0t;  // for c0t

        const Complex dW1 = iz * dZ1 * dw(Z1, W1);  // NEED TO CHECK DW!

        dAterm = sqrt_pi * (W1 * dcte + cte * dW1);  // for all

        if (abs(Z1) <= 4e3)
          dBterm =
              -(sqrt_pi * ((pow2(Z1) - 1) * dW1 + 2 * W1 * Z1 * dZ1) - dZ1) *
                  cte -
              (sqrt_pi * (pow2(Z1) - 1) * W1 - Z1) * dcte;  // for all
        else
          dBterm =
              ((sqrt_pi * W1 * pow3(Z1) + 0.5 * pow2(Z1) - 0.75) * Z1 * dcte +
               (sqrt_pi * pow4(Z1) * dW1 - 0.5 * pow2(Z1) * dZ1 + 2.25 * dZ1) *
                   cte) /
              pow4(Z1);  // for all
      } else if (abs(X) <= 3e-8 * abs(Y)) {
        Complex dZ1;
        if (is_magnetic_parameter(rt))
          dZ1 = -iz * cte * freq2kaycm(zeeman_df_si);  // for H
        else if (rt == JacPropMatType::LineCenter and
                 rt.QuantumIdentity().In(quantum_identity))
          dZ1 = -iz * cte;  // for sg0
        else if (is_frequency_parameter(rt))
          dZ1 = iz * cte;  // for sg
        else
          dZ1 = (iz * (freq2kaycm(f_grid[iv]) - sg0) + c0t) * dcte +
                cte * dc0t;  // for c0t

        const Complex dZ2 = dY / (2 * sqrtY) + dX / (2 * sqrtXY) +
                            dY / (2 * sqrtXY);  // for all

        const Complex dW1 = iz * dZ1 * dw(Z1, W1);  // NEED TO CHECK DW!
        const Complex dW2 = iz * dZ2 * dw(Z2, W2);  // NEED TO CHECK DW!

        dAterm = sqrt_pi * ((W1 - W2) * dcte + (dW1 - dW2) * cte);  // for all

        dBterm = (sqrt_pi *
                      (((pow2(Z1) - 1) * W1 - (pow2(Z2) - 1) * W2) * dY +
                       2 *
                           (-(pow2(Z1) - 1) * dW1 + (pow2(Z2) - 1) * dW2 -
                            2 * W1 * Z1 * dZ1 + 2 * W2 * Z2 * dZ2) *
                           Y) *
                      c2t +
                  2 *
                      (sqrt_pi * (pow2(Z1) - 1) * W1 -
                       sqrt_pi * (pow2(Z2) - 1) * W2 + 2 * sqrtY) *
                      Y * dc2t) /
                 (4 * Y * sqrtY * pow2(c2t));  // for all
      } else if (abs(Y) <= 1e-15 * abs(X)) {
        const Complex dZ1 = (dX + dY) / (2 * sqrtXY);  // for all
        const Complex dW1 = iz * dZ1 * dw(Z1, W1);     // NEED TO CHECK DW!
        if (abs(sqrtX) <= 4e3) {
          const Complex dZb = dX / (2 * sqrtX);       // for all
          const Complex dWb = iz * dZb * dw(Zb, Wb);  // NEED TO CHECK DW!

          dAterm = (-sqrt_pi * (Wb * dX + 2 * X * dWb) * c2t +
                    2 * (sqrt_pi * Wb * sqrtX - 1) * sqrtX * dc2t) /
                   (sqrtX * pow2(c2t));  // for all

          dBterm =
              (-sqrtXY *
                   (2 * (sqrt_pi * Wb * sqrtX - 1) * (X + 2 * Y - 1) +
                    2 * sqrt_pi * sqrtXY * W1 - 1) *
                   sqrtX * dc2t +
               (2 *
                    ((sqrt_pi * Wb * sqrtX - 1) * (dX + 2 * dY) +
                     sqrt_pi * sqrtXY * dW1) *
                    sqrtXY * sqrtX +
                sqrt_pi * (Wb * dX + 2 * X * dWb) * sqrtXY * (X + 2 * Y - 1) +
                sqrt_pi * (dX + dY) * W1 * sqrtX) *
                   c2t) /
              (sqrtXY * sqrtX * pow2(c2t));  // for all
        } else {
          dAterm = ((-X + 3.0) * c2t * dX - (X - 1.5) * X * dc2t) /
                   (pow3(X) * pow2(c2t));  // for all

          dBterm = (((-2 * sqrt_pi * sqrtXY * W1 + 1) * pow2(X) +
                     (X - 1.5) * (X + 2 * Y - 1)) *
                        sqrtXY * X * dc2t +
                    ((X - 3.0) * sqrtXY * (X + 2 * Y - 1) * dX -
                     (X - 1.5) * sqrtXY * (dX + 2 * dY) * X +
                     2 * sqrt_pi * (X + Y) * pow3(X) * dW1 +
                     sqrt_pi * (dX + dY) * W1 * pow3(X)) *
                        c2t) /
                   (sqrtXY * pow3(X) * pow2(c2t));  // for all
        }
      } else {
        const Complex dZ1 = -dY / (2 * sqrtY) + dX / (2 * sqrtXY) +
                            dY / (2 * sqrtXY);  // for all
        const Complex dZ2 = dY / (2 * sqrtY) + dX / (2 * sqrtXY) +
                            dY / (2 * sqrtXY);  // for all

        const Complex dW1 = iz * dZ1 * dw(Z1, W1);  // NEED TO CHECK DW!
        const Complex dW2 = iz * dZ2 * dw(Z2, W2);  // NEED TO CHECK DW!

        dAterm = sqrt_pi * ((W1 - W2) * dcte + (dW1 - dW2) * cte);  // for all

        dBterm = (sqrt_pi *
                      (((pow2(Z1) - 1) * W1 - (pow2(Z2) - 1) * W2) * dY +
                       2 *
                           (-(pow2(Z1) - 1) * dW1 + (pow2(Z2) - 1) * dW2 -
                            2 * W1 * Z1 * dZ1 + 2 * W2 * Z2 * dZ2) *
                           Y) *
                      c2t +
                  2 *
                      (sqrt_pi * (pow2(Z1) - 1) * W1 -
                       sqrt_pi * (pow2(Z2) - 1) * W2 + 2 * sqrtY) *
                      Y * dc2t) /
                 (4 * Y * sqrtY * pow2(c2t));  // for all
      }

      Numeric dFVC = 0, dETA = 0;
      if (rt == JacPropMatType::Temperature) {
        dETA = dT.ETA;
        dFVC = dT.FVC;
      } else if (rt == JacPropMatType::VMR) {
        dETA = dV.ETA;
        dFVC = dV.FVC;
      } else if (is_pressure_broadening_ETA(rt) and
                 rt.QuantumIdentity().In(quantum_identity))
        dETA = 1;
      else if (is_pressure_broadening_FVC(rt) and
               rt.QuantumIdentity().In(quantum_identity))
        dFVC = 1;

      dF.col(iq)(iv) =
          ((((c0 - 1.5 * c2) * x.ETA - x.FVC) * Aterm + Bterm * c2 * x.ETA +
            1) *
               dAterm -
           (((c0 - 1.5 * c2) * x.ETA - x.FVC) * dAterm +
            ((c0 - 1.5 * c2) * dETA + (dc0 - 1.5 * dc2) * x.ETA - dFVC) *
                Aterm +
            Bterm * c2 * dETA + Bterm * x.ETA * dc2 + c2 * x.ETA * dBterm) *
               Aterm) /
          (pi * pow2(((c0 - 1.5 * c2) * x.ETA - x.FVC) * Aterm +
                     Bterm * c2 * x.ETA + 1));  // for all
    }
  }

  // Convert back to ARTS
  F = F.unaryExpr(&pCqSDHC_to_arts);
  for (auto iq = 0; iq < derivatives_data_position.nelem(); iq++) {
    const RetrievalQuantity& rt =
        derivatives_data[derivatives_data_position[iq]];

    if (rt == JacPropMatType::LineCenter or is_frequency_parameter(rt) or
        is_pressure_broadening_G0(rt) or is_pressure_broadening_D0(rt) or
        is_pressure_broadening_FVC(rt))
      dF.col(iq) = dF.col(iq).unaryExpr(&pCqSDHC_to_arts_freq_deriv);
    else if (is_pressure_broadening_G2(rt))
      dF.col(iq) = dF.col(iq).unaryExpr(&pCqSDHC_to_arts_G2_deriv);
    else if (is_pressure_broadening_D2(rt))
      dF.col(iq) = dF.col(iq).unaryExpr(&pCqSDHC_to_arts_D2_deriv);
    else
      dF.col(iq) = dF.col(iq).unaryExpr(&pCqSDHC_to_arts);
  }
}
