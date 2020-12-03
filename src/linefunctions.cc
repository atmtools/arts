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
    const Absorption::SingleLine& line,
    const Numeric& temperature,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& doppler_constant,
    const LineShape::Output& X,
    const LineShape::Type lineshape_type,
    const Absorption::MirroringType mirroring_type,
    const Absorption::NormalizationType norm_type)
{
  Eigen::MatrixXcd dF(0, 0), data(F.size(), Linefunctions::ExpectedDataSize());

  switch (lineshape_type) {
    case LineShape::Type::HTP:
    case LineShape::Type::SDVP:
      set_htp(F,
              dF,
              f_grid,
              zeeman_df,
              magnetic_magnitude,
              line.F0(),
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
                line.F0(),
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
                  line.F0(),
                  doppler_constant);
      break;
    case LineShape::Type::LP:
      set_lorentz(
          F, dF, data, f_grid, zeeman_df, magnetic_magnitude, line.F0(), X);
      break;
  }

  switch (mirroring_type) {
    case Absorption::MirroringType::None:
    case Absorption::MirroringType::Manual:
      break;
    case Absorption::MirroringType::Lorentz: {
      // Set the mirroring computational vectors and size them as needed
      Eigen::VectorXcd Fm(F.size());

      set_lorentz(Fm,
                  dF,
                  data,
                  f_grid,
                  -zeeman_df,
                  magnetic_magnitude,
                  -line.F0(),
                  LineShape::mirroredOutput(X));

      // Apply mirroring;  FIXME: Add conjugate?
      F.noalias() += Fm;
    } break;
    case Absorption::MirroringType::SameAsLineShape: {
      // Set the mirroring computational vectors and size them as needed
      Eigen::VectorXcd Fm(F.size());

      switch (lineshape_type) {
        case LineShape::Type::DP:
          set_doppler(Fm,
                      dF,
                      data,
                      f_grid,
                      -zeeman_df,
                      magnetic_magnitude,
                      -line.F0(),
                      -doppler_constant);
          break;
        case LineShape::Type::LP:
          set_lorentz(Fm,
                      dF,
                      data,
                      f_grid,
                      -zeeman_df,
                      magnetic_magnitude,
                      -line.F0(),
                      LineShape::mirroredOutput(X));
          break;
        case LineShape::Type::VP:
          set_voigt(Fm,
                    dF,
                    data,
                    f_grid,
                    -zeeman_df,
                    magnetic_magnitude,
                    -line.F0(),
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
                  -line.F0(),
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
  switch (norm_type) {
    case Absorption::NormalizationType::None:
      break;
    case Absorption::NormalizationType::VVH:
      apply_VVH_scaling(F, dF, data, f_grid, line.F0(), temperature);
      break;
    case Absorption::NormalizationType::VVW:
      apply_VVW_scaling(F, dF, f_grid, line.F0());
      break;
    case Absorption::NormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F, dF, f_grid, line.F0(), temperature);
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
    const LineShape::Output& lso,
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const LineShape::Output& dT,
    const LineShape::Output& dVMR) {
  constexpr Complex cpi(0, Constant::pi);
  constexpr Complex iz(0.0, 1.0);

  // Size of the problem
  auto nppd = derivatives_data_position.nelem();

  // The central frequency
  const Numeric F0 = F0_noshift + lso.D0 + zeeman_df * magnetic_magnitude + lso.DV;

  // Naming req data blocks
  auto z = data.col(0);
  auto dw = data.col(1);

  z.noalias() =
      (Constant::pi * Complex(lso.G0, F0) - cpi * f_grid.array()).matrix();
  F.noalias() = z.cwiseInverse();

  if (nppd) {
    dw.noalias() = -Constant::pi * F.array().square().matrix();

    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (deriv == Jacobian::Atm::Temperature)
        dF.col(iq).noalias() = Complex(dT.G0, dT.D0 + dT.DV) * dw;
      else if (is_frequency_parameter(deriv))
        dF.col(iq).noalias() = -iz * dw;
      else if ((deriv == Jacobian::Line::Center or
                is_pressure_broadening_DV(deriv)) and
                Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = iz * dw;
      else if (is_pressure_broadening_G0(deriv) and
        Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = dw;
      else if (is_pressure_broadening_D0(deriv) and
        Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = iz * dw;
      else if (is_magnetic_parameter(deriv))
        dF.col(iq).noalias() = iz * zeeman_df * dw;
      else if (deriv == Jacobian::Line::VMR) {
        if (Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
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
    const LineShape::Output& lso,
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const Numeric& dGD_div_F0_dT,
    const LineShape::Output& dT,
    const LineShape::Output& dVMR) {
  constexpr Complex iz(0.0, 1.0);

  // Size of problem
  auto nppd = derivatives_data_position.nelem();

  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + lso.D0 + lso.DV;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0 - GD_div_F0 * (dT.D0 + dT.DV);

  // constant normalization factor for Voigt
  const Numeric fac = Constant::inv_sqrt_pi * invGD;

  // Naming req data blocks
  auto z = data.col(0);
  auto dw = data.col(1);

  // Frequency grid
  z.noalias() = invGD * (Complex(-F0, lso.G0) + f_grid.array()).matrix();

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
      else if (deriv == Jacobian::Atm::Temperature)
        dF.col(iq).noalias() =
            dw * Complex(-dT.D0 - dT.DV, dT.G0) * invGD -
            F * dGD_dT * invGD - dw.cwiseProduct(z) * dGD_dT * invGD;
      else if ((deriv == Jacobian::Line::Center or
                is_pressure_broadening_DV(deriv)) and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = -F / F0 - dw * invGD - dw.cwiseProduct(z) / F0;
      else if (is_pressure_broadening_G0(deriv) and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = iz * dw * invGD;
      else if (is_pressure_broadening_D0(deriv) and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = -dw * invGD;
      else if (is_magnetic_parameter(deriv))
        dF.col(iq).noalias() = dw * (-zeeman_df * invGD);
      else if (deriv == Jacobian::Line::VMR and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() =
            dw * Complex(-dVMR.D0 - dVMR.DV, dVMR.G0) * invGD;
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
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
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
    else if (deriv == Jacobian::Atm::Temperature)
      dF.col(iq).noalias() =
          -dGD_div_F0_dT * F0 * invGD * (2.0 * F.cwiseProduct(mx2) + F);
    else if (deriv == Jacobian::Line::Center and
             Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
      dF.col(iq).noalias() = -F.cwiseProduct(mx2) * (2 / F0) +
                             F.cwiseProduct(x) * 2 * (invGD - 1 / F0);
    else if (deriv == Jacobian::Line::VMR)
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
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const LineShape::Output& dT,
    const LineShape::Output& dVMR) {
  auto nppd = derivatives_data_position.nelem();

  const Complex LM = Complex(1.0 + X.G, -X.Y);

  dF *= LM;
  if (with_mirroring) {
    dF.noalias() += dFm * conj(LM);
    
    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (deriv == Jacobian::Atm::Temperature) {
        const auto c = Complex(dT.G, -dT.Y);
        dF.col(iq).noalias() += F * c + Fm * conj(c);
      } else if (is_pressure_broadening_G(deriv) and
                 Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = F + Fm;
      else if (is_pressure_broadening_Y(deriv) and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = Complex(0, -1) * (F - Fm);
      else if (deriv == Jacobian::Line::VMR and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind)) {
        const auto c = Complex(dVMR.G, -dVMR.Y);
        dF.col(iq).noalias() += F * c + Fm * conj(c);
      }
    }
  } else {
    for (auto iq = 0; iq < nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];

      if (deriv == Jacobian::Atm::Temperature)
        dF.col(iq).noalias() += F * Complex(dT.G, -dT.Y);
      else if (is_pressure_broadening_G(deriv) and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = F;
      else if (is_pressure_broadening_Y(deriv) and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
        dF.col(iq).noalias() = Complex(0, -1) * F;
      else if (deriv == Jacobian::Line::VMR and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
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
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position) {
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
      if (deriv == Jacobian::Atm::Temperature)
        dF(iv, iq) += dmafac_dT_div_fun * F[iv];
      else if (deriv == Jacobian::Line::Center and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
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
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position) {
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

    if (deriv == Jacobian::Atm::Temperature)
      dF.col(iq).noalias() +=
          -c1 / T *
          ((denom - F0 / tanh_f0part) * F - denom * ftanh_f.cwiseProduct(F) +
           f_grid.cwiseProduct(F).cwiseQuotient(tanh_f));
    else if (is_frequency_parameter(deriv))
      dF.col(iq).noalias() +=
          F.cwiseQuotient(f_grid) +
          c1 * (F.cwiseQuotient(tanh_f) - F.cwiseProduct(tanh_f));
    else if (deriv == Jacobian::Line::Center and
             Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
      dF.col(iq).noalias() +=
          (-1.0 / F0 + c1 * tanh_f0part - c1 / tanh_f0part) * F;
  }
}

void Linefunctions::apply_VVW_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position) {
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
      if (deriv == Jacobian::Line::Center and
          Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
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
    const Absorption::SingleLine& line,
    const Numeric& T,
    const Numeric& T0,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const Numeric& dQT_dT) {
  auto nppd = derivatives_data_position.nelem();

  const Numeric gamma = stimulated_emission(T, line.F0());
  const Numeric gamma_ref = stimulated_emission(T0, line.F0());
  const Numeric K1 = boltzman_ratio(T, T0, line.E0());
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);

  const Numeric invQT = 1.0 / QT;
  const Numeric S = line.I0() * isotopic_ratio * QT0 * invQT * K1 * K2;

  F *= S;
  dF *= S;
  for (auto iq = 0; iq < nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];

    if (deriv == Jacobian::Atm::Temperature)
      dF.col(iq).noalias() +=
      F * (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F0(), T) /
                   K2 +
                   dboltzman_ratio_dT_div_boltzmann_ratio(T, line.E0()) - invQT * dQT_dT);
    else if (deriv == Jacobian::Line::Strength and
             Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
      dF.col(iq).noalias() = F / line.I0();  //nb. overwrite
    else if (deriv == Jacobian::Line::Center and
             Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind))
      dF.col(iq).noalias() +=
          F *
          dstimulated_relative_emission_dF0(gamma, gamma_ref, T, T0) /
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
    const Absorption::SingleLine& line,
    const Numeric& T,
    const Numeric& T0,
    const Numeric& Tu,
    const Numeric& Tl,
    const Numeric& Evu,
    const Numeric& Evl,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const Numeric& dQT_dT) {
  auto nppd = derivatives_data_position.nelem();

  const Numeric gamma = stimulated_emission(T,line.F0());
  const Numeric gamma_ref = stimulated_emission(T0,line.F0());
  const Numeric r_low = boltzman_ratio(Tl,T,Evl);
  const Numeric r_upp = boltzman_ratio(Tu,T,Evu);

  const Numeric K1 = boltzman_ratio(T,T0,line.E0());
  const Numeric K2 = stimulated_relative_emission(gamma,gamma_ref);
  const Numeric K3 = absorption_nlte_ratio(gamma,r_upp,r_low);
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

    if (deriv == Jacobian::Atm::Temperature) {
      const Numeric dS_dT_abs =
          S_abs *
          (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F0(), T) /
               K2 +
               dboltzman_ratio_dT(K1, T, line.E0()) / K1 +
               dabsorption_nlte_rate_dT(gamma, T, line.F0(), Evl, Evu, K4, r_low) /
               K3 -
           invQT * dQT_dT);
      const Numeric dS_dT_src =
          S_src *
          (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F0(), T) /
               K2 +
               dboltzman_ratio_dT(K1, T, line.E0()) / K1 -
           dboltzman_ratio_dT(K4, T, Evu) / K4 - invQT * dQT_dT);

      dN.col(iq).noalias() += F * (dS_dT_src - dS_dT_abs);
      dF.col(iq).noalias() += F * dS_dT_abs;
    } else if (deriv == Jacobian::Line::Strength and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind)) {
      dF.col(iq).noalias() = F * dS_dS0_abs;
      dN.col(iq).noalias() = F * (dS_dS0_src - dS_dS0_abs);
    } else if (deriv == Jacobian::Line::Center and
               Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind)) {
      const Numeric dS_dF0_abs =
          S_abs *
          (dstimulated_relative_emission_dF0(gamma, gamma_ref, T, T0) /
               K2 +
           dabsorption_nlte_rate_dF0(gamma, T, K4, r_low) / K3);
      const Numeric dS_dF0_src =
          S_src *
          dstimulated_relative_emission_dF0(gamma, gamma_ref, T, T0) /
          K2;

      dN.col(iq).noalias() += F * (dS_dF0_src - dS_dF0_abs);
      dF.col(iq).noalias() += F * dS_dF0_abs;
    } else if (deriv == Jacobian::Line::NLTE) {
      if (Absorption::id_in_line_lower(band, deriv.QuantumIdentity(), line_ind)) {
        const Numeric dS_dTl_abs =
            S_abs * dabsorption_nlte_rate_dTl(gamma, T, Tl, Evl, r_low) / K3;

        dN.col(iq).noalias() = -F * dS_dTl_abs;
        dF.col(iq).noalias() = -dN.col(iq);
      } else if (Absorption::id_in_line_upper(band, deriv.QuantumIdentity(), line_ind)) {
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

void Linefunctions::apply_lineshapemodel_jacobian_scaling(
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const Numeric& T,
    const Numeric& P,
    const Vector& vmrs) {
  auto nppd = derivatives_data_position.nelem();

  for (auto iq = 0; iq < nppd; iq++) {
    const RetrievalQuantity& rt =
        derivatives_data[derivatives_data_position[iq]];
    if (is_lineshape_parameter(rt) and
        Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind))
      dF.col(iq) *= band.ShapeParameter_dInternal(line_ind, T, P, vmrs, rt);
  }
}

Numeric Linefunctions::DopplerConstant(Numeric T, Numeric mass) {
  return std::sqrt(Constant::doppler_broadening_const_squared * T / mass);
}

Numeric Linefunctions::dDopplerConstant_dT(const Numeric& T,
                                           const Numeric& dc) {
  return dc / (2 * T);
}

void Linefunctions::find_cutoff_ranges(
    Index& start_cutoff,
    Index& nelem_cutoff,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& fmin,
    const Numeric& fmax) {
  auto nf = f_grid.size();

  const bool need_cutoff = (fmax > fmin);
  if (need_cutoff) {
    // Find range of simulations
    start_cutoff = 0;
    Index i_f_max = nf - 1;

    // Loop over positions to compute the line shape cutoff point
    while (start_cutoff < nf and fmin > f_grid[start_cutoff])
      ++start_cutoff;
    while (i_f_max >= start_cutoff and fmax < f_grid[i_f_max])
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
    const AbsorptionLines& band,
    const Index& line_ind,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position) {
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

    if (deriv == Jacobian::Atm::Temperature)
      dN.col(iq).noalias() +=
          F * e * Constant::h * F0 * exp_T / (c2 * Constant::k * T * T);
    else if (deriv == Jacobian::Line::Center and
             Absorption::id_in_line(band, deriv.QuantumIdentity(), line_ind)) {
      Numeric done_over_b_df0 =
          Constant::h * exp_T / (c2 * Constant::k * T) - 3.0 * b / F0;
      Numeric de_df0 = c1 * r2 * A21;
      Numeric dk_df0 = c1 * (r1 * x - r2) * (A21 / c2) - 3.0 * k / F0;

      dN.col(iq).noalias() += F * (e * done_over_b_df0 + de_df0 / b - dk_df0);
      dF.col(iq).noalias() += F * dk_df0;
    } else if (deriv == Jacobian::Line::NLTE) {
      if (Absorption::id_in_line_lower(band, deriv.QuantumIdentity(), line_ind)) {
        const Numeric dk_dr2 = -c3 * A21 / c2, de_dr2 = c3 * A21,
                      dratio_dr2 = de_dr2 / b - dk_dr2;
        dN.col(iq).noalias() = F * dratio_dr2;
        dF.col(iq).noalias() = F * dk_dr2;
      } else if (Absorption::id_in_line_upper(band, deriv.QuantumIdentity(), line_ind)) {
        const Numeric dk_dr1 = c3 * x * A21 / c2;
        dF.col(iq).noalias() = F * dk_dr1;
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
                            const LineShape::Output& lso_si,
                            const AbsorptionLines& band,
                            const Index& line_ind,
                            const ArrayOfRetrievalQuantity& derivatives_data,
                            const ArrayOfIndex& derivatives_data_position,
                            const Numeric& dGD_div_F0_dT_si,
                            const LineShape::Output& dT_si,
                            const LineShape::Output& dVMR_si) {
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
  const auto lso = si2cgs(lso_si);
  const auto dT = si2cgs(dT_si);
  const auto dV = si2cgs(dVMR_si);

  // General normalization
  const Numeric cte = sqrt_ln_2 / GamD;
  constexpr Complex iz(0, 1);

  // Calculating the different parameters
  const Complex c0(lso.G0, -lso.D0);
  const Complex c2(lso.G2, -lso.D2);
  const Complex c0t = (1 - lso.ETA) * (c0 - 1.5 * c2) + lso.FVC;
  const Complex c2t = (1 - lso.ETA) * c2;
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

    F(iv) = Aterm / (pi * (((c0 - 1.5 * c2) * lso.ETA - lso.FVC) * Aterm +
                           Bterm * c2 * lso.ETA + 1));

    for (auto iq = 0; iq < derivatives_data_position.nelem(); iq++) {
      const RetrievalQuantity& rt =
          derivatives_data[derivatives_data_position[iq]];

      Numeric dcte = 0;
      Complex dc0 = 0, dc2 = 0, dc0t = 0, dc2t = 0;
      if (rt == Jacobian::Atm::Temperature) {
        dcte = (-(dGD_div_F0_dT_si * sg0 / Constant::sqrt_ln_2) / GamD) *
               cte;                    // for T
        dc0 = Complex(dT.G0, -dT.D0);  // for T
        dc2 = Complex(dT.G2, -dT.D2);  // for T
        dc0t = (-(c0 - 1.5 * c2) * dT.ETA + dT.FVC) +
               ((1 - lso.ETA) * (dc0 - 1.5 * dc2));     // for T
        dc2t = (-c2 * dT.ETA) + ((1 - lso.ETA) * dc2);  // for T
      } else if (rt == Jacobian::Line::VMR and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc0 = Complex(dV.G0, -dV.D0);  // for VMR
        dc2 = Complex(dV.G2, -dV.D2);  // for VMR
        dc0t = (-(c0 - 1.5 * c2) * dV.ETA + dV.FVC) +
               ((1 - lso.ETA) * (dc0 - 1.5 * dc2));     // for VMR
        dc2t = (-c2 * dV.ETA) + ((1 - lso.ETA) * dc2);  // for VMR
      } else if (is_pressure_broadening_G2(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc2 =
            1;  // for Gam2(T, VMR) and for Shift2(T, VMR), the derivative is wrt c2 for later computations
        dc0t = ((1 - lso.ETA) * (dc0 - 1.5 * dc2));
        dc2t = ((1 - lso.ETA) * dc2);
      } else if (is_pressure_broadening_D2(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc2 =
            -iz;  // for Gam2(T, VMR) and for Shift2(T, VMR), the derivative is wrt c2 for later computations
        dc0t = ((1 - lso.ETA) * (dc0 - 1.5 * dc2));
        dc2t = ((1 - lso.ETA) * dc2);
      } else if (is_pressure_broadening_G0(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc0 =
            1;  // for Gam0(T, VMR) and for Shift0(T, VMR), the derivative is wrt c0 for later computations
        dc0t = ((1 - lso.ETA) * (dc0 - 1.5 * dc2));
      } else if (is_pressure_broadening_D0(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc0 =
            -iz;  // for Gam0(T, VMR) and for Shift0(T, VMR), the derivative is wrt c0 for later computations
        dc0t = ((1 - lso.ETA) * (dc0 - 1.5 * dc2));
      } else if (is_pressure_broadening_FVC(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc0t = (1);  // for FVC(T, VMR)
      } else if (is_pressure_broadening_ETA(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind)) {
        dc0t = (-c0 + 1.5 * c2);  // for eta(T, VMR)
        dc2t = (-c2);             // for eta(T, VMR)
      }

      const Complex dY = (-2 * dcte / cte - 2 * dc2t / c2t) * Y;  // for all

      Complex dX;
      if (is_magnetic_parameter(rt))
        dX = -iz * freq2kaycm(zeeman_df_si) / c2t;  // for H
      else if (rt == Jacobian::Line::Center and
               Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind))
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
        else if (rt == Jacobian::Line::Center and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind))
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
        else if (rt == Jacobian::Line::Center and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind))
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
      if (rt == Jacobian::Atm::Temperature) {
        dETA = dT.ETA;
        dFVC = dT.FVC;
      } else if (rt == Jacobian::Line::VMR) {
        dETA = dV.ETA;
        dFVC = dV.FVC;
      } else if (is_pressure_broadening_ETA(rt) and
                 Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind))
        dETA = 1;
      else if (is_pressure_broadening_FVC(rt) and
               Absorption::id_in_line(band, rt.QuantumIdentity(), line_ind))
        dFVC = 1;

      dF.col(iq)(iv) =
          ((((c0 - 1.5 * c2) * lso.ETA - lso.FVC) * Aterm + Bterm * c2 * lso.ETA +
            1) *
               dAterm -
           (((c0 - 1.5 * c2) * lso.ETA - lso.FVC) * dAterm +
            ((c0 - 1.5 * c2) * dETA + (dc0 - 1.5 * dc2) * lso.ETA - dFVC) *
                Aterm +
            Bterm * c2 * dETA + Bterm * lso.ETA * dc2 + c2 * lso.ETA * dBterm) *
               Aterm) /
          (pi * pow2(((c0 - 1.5 * c2) * lso.ETA - lso.FVC) * Aterm +
                     Bterm * c2 * lso.ETA + 1));  // for all
    }
  }

  // Convert back to ARTS
  F = F.unaryExpr(&pCqSDHC_to_arts);
  for (auto iq = 0; iq < derivatives_data_position.nelem(); iq++) {
    const RetrievalQuantity& rt =
        derivatives_data[derivatives_data_position[iq]];

    if (rt == Jacobian::Line::Center or is_frequency_parameter(rt) or
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

void Linefunctions::set_cross_section_of_band(
    InternalData& scratch,
    InternalData& sum,
    const ConstVectorView f_grid,
    const AbsorptionLines& band,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_active,
    const Vector& vmrs,
    const EnergyLevelMap& nlte,
    const Numeric& P,
    const Numeric& T,
    const Numeric& isot_ratio,
    const Numeric& H,
    const Numeric& DC,
    const Numeric& dDCdT,
    const Numeric& QT,
    const Numeric& dQTdT,
    const Numeric& QT0,
    const bool no_negatives,
    const bool zeeman,
    const Zeeman::Polarization zeeman_polarization)
{
  const Index nj = derivatives_data_active.nelem();
  const bool do_temperature = do_temperature_jacobian(derivatives_data);
  
  // Sum up variable reset
  sum.SetZero();
  
  if (band.NumLines() == 0 or (Absorption::relaxationtype_relmat(band.Population()) and band.LinemixingLimit() > P)) {
    return;  // No line-by-line computations required/wanted
  }
  
  // Cutoff for Eigen-library types
  Eigen::Matrix<Numeric, 1, 1> fc;
  auto& Fc = scratch.Fc;
  auto& Nc = scratch.Nc;
  auto& dFc = scratch.dFc;
  auto& dNc = scratch.dNc;
  auto& datac = scratch.datac;
  
  // Frequency grid as Eigen type
  const auto f_full = MapToEigen(f_grid);
  
  // Cut off range
  const Numeric fmean = (band.Cutoff() == Absorption::CutoffType::BandFixedFrequency) ? band.F_mean() : 0;
  Numeric fcut_upp, fcut_low;
  Index start, nelem;
  fcut_upp = band.CutoffFreq(0);
  fcut_low = band.CutoffFreqMinus(0, fmean);
  find_cutoff_ranges(start, nelem, f_full, fcut_low, fcut_upp);
  fc[0] = fcut_upp;
  
  // VMR Jacobian check
  auto do_vmr = do_vmr_jacobian(derivatives_data, band.QuantumIdentity());
  
  // Placeholder nothingness
  constexpr LineShape::Output empty_output = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  
  for (Index i=0; i<band.NumLines(); i++) {
    
    // Select the range of cutoff if different for each line
    if (band.Cutoff() == Absorption::CutoffType::LineByLineOffset and i>0) {
      fcut_upp = band.CutoffFreq(i);
      fcut_low = band.CutoffFreqMinus(i, fmean);
      find_cutoff_ranges(start, nelem, f_full, fcut_low, fcut_upp);
      fc[0] = fcut_upp;
    }
    
    // Relevant range FIXME: By Band and no-cutoff does not need this...
    auto F = scratch.F.segment(start, nelem);
    auto N = scratch.N.segment(start, nelem);
    auto dF = scratch.dF.middleRows(start, nelem);
    auto dN = scratch.dN.middleRows(start, nelem);
    auto data = scratch.data.middleRows(start, nelem);
    const auto f = f_full.middleRows(start, nelem);
    
    // Pressure broadening and line mixing terms
    const auto X = band.ShapeParameters(i, T, P, vmrs);
    
    // Partial derivatives for temperature
    const auto dXdT = do_temperature ?
      band.ShapeParameters_dT(i, T, P, vmrs) : empty_output;
    
    // Partial derivatives for VMR of self (function works for any species but only do self for now)
    const auto dXdVMR = do_vmr.test ?
      band.ShapeParameters_dVMR(i, T, P, do_vmr.qid) : empty_output;
    
    // Zeeman lines if necessary
    const Index nz = zeeman ?
      band.ZeemanCount(i, zeeman_polarization) : 1;
    
    for (Index iz=0; iz<nz; iz++) {
      
      // Zeeman values for this sub-line
      const Numeric Sz = zeeman ?
        band.ZeemanStrength(i, zeeman_polarization, iz) : 1;
      const Numeric dfdH = zeeman ?
        band.ZeemanSplitting(i, zeeman_polarization, iz) : 0;
      
      // Set the line shape and its derivatives
      switch (band.LineShapeType()) {
        case LineShape::Type::DP:
          set_doppler(F, dF, data, f, dfdH, H, band.F0(i), DC, band, i, derivatives_data, derivatives_data_active, dDCdT);
          if (band.Cutoff() not_eq Absorption::CutoffType::None)
            set_doppler(Fc, dFc, datac, fc, dfdH, H, band.F0(i), DC, band, i, derivatives_data, derivatives_data_active, dDCdT);
          break;
        case LineShape::Type::HTP:
        case LineShape::Type::SDVP:
          set_htp(F, dF, f, dfdH, H, band.F0(i), DC, X, band, i, derivatives_data, derivatives_data_active, dDCdT, dXdT, dXdVMR);
          if (band.Cutoff() not_eq Absorption::CutoffType::None)
            set_htp(Fc, dFc, fc, dfdH, H, band.F0(i), DC, X, band, i, derivatives_data, derivatives_data_active, dDCdT, dXdT, dXdVMR);
          break;
        case LineShape::Type::LP:
          set_lorentz(F, dF, data, f, dfdH, H, band.F0(i), X, band, i, derivatives_data, derivatives_data_active, dXdT, dXdVMR);
          if (band.Cutoff() not_eq Absorption::CutoffType::None)
            set_lorentz(Fc, dFc, datac, fc, dfdH, H, band.F0(i), X, band, i, derivatives_data, derivatives_data_active, dXdT, dXdVMR);
          break;
        case LineShape::Type::VP:
          set_voigt(F, dF, data, f, dfdH, H, band.F0(i), DC, X, band, i, derivatives_data, derivatives_data_active, dDCdT, dXdT, dXdVMR);
          if (band.Cutoff() not_eq Absorption::CutoffType::None)
            set_voigt(Fc, dFc, datac, fc, dfdH, H, band.F0(i), DC, X, band, i, derivatives_data, derivatives_data_active, dDCdT, dXdT, dXdVMR);
          break;
      }
      
      // Remove the cutoff values
      if (band.Cutoff() not_eq Absorption::CutoffType::None) {
        F.array() -= Fc[0];
        for (Index ij = 0; ij < nj; ij++) {
          dF.col(ij).array() -= dFc[ij];
        }
      }

      // Set the mirrored line shape
      const bool with_mirroring =
      band.Mirroring() not_eq Absorption::MirroringType::None and
      band.Mirroring() not_eq Absorption::MirroringType::Manual;
      switch (band.Mirroring()) {
        case Absorption::MirroringType::None:
        case Absorption::MirroringType::Manual:
          break;
        case Absorption::MirroringType::Lorentz:
          set_lorentz(N, dN, data, f, -dfdH, H, -band.F0(i), LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
          if (band.Cutoff() not_eq Absorption::CutoffType::None)
            set_lorentz(Nc, dNc, datac, fc, -dfdH, H, -band.F0(i), LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
          break;
        case Absorption::MirroringType::SameAsLineShape:
          switch (band.LineShapeType()) {
            case LineShape::Type::DP:
              set_doppler(N, dN, data, f, -dfdH, H, -band.F0(i), -DC, band, i, derivatives_data, derivatives_data_active, -dDCdT);
              if (band.Cutoff() not_eq Absorption::CutoffType::None)
                set_doppler(Nc, dNc, datac, fc, -dfdH, H, -band.F0(i), -DC, band, i, derivatives_data, derivatives_data_active, -dDCdT);
              break;
            case LineShape::Type::LP:
              set_lorentz(N, dN, data, f, -dfdH, H, -band.F0(i), LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
              if (band.Cutoff() not_eq Absorption::CutoffType::None)
                set_lorentz(Nc, dNc, datac, fc, -dfdH, H, -band.F0(i), LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
              break;
            case LineShape::Type::VP:
              set_voigt(N, dN, data, f, -dfdH, H, -band.F0(i), -DC, LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, -dDCdT, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
              if (band.Cutoff() not_eq Absorption::CutoffType::None)
                set_voigt(Nc, dNc, datac, fc, -dfdH, H, -band.F0(i), -DC, LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, -dDCdT, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
              break;
            case LineShape::Type::HTP:
            case LineShape::Type::SDVP:
              // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
              set_htp(N, dN, f, -dfdH, H, -band.F0(i), -DC, LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, -dDCdT, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
              if (band.Cutoff() not_eq Absorption::CutoffType::None)
                set_htp(Nc, dNc, fc, -dfdH, H, -band.F0(i), -DC, LineShape::mirroredOutput(X), band, i, derivatives_data, derivatives_data_active, -dDCdT, do_temperature ? LineShape::mirroredOutput(dXdT) : empty_output, do_vmr.test ? LineShape::mirroredOutput(dXdVMR) : empty_output);
              break;
          }
          break;
      }
      
      // Remove the mirrored cutoff values
      if (band.Cutoff() not_eq Absorption::CutoffType::None and with_mirroring) {
        N.array() -= Nc[0];
        for (Index ij = 0; ij < nj; ij++) {
          dN.col(ij).array() -= dNc[ij];
        }
      }

      // Mirror and and line mixing is added together (because of conjugate)
      if (band.LineShapeType() not_eq LineShape::Type::DP) {
        apply_linemixing_scaling_and_mirroring(F, dF, N, dN, X, with_mirroring, band, i, derivatives_data, derivatives_data_active, dXdT, dXdVMR);

        // Apply line mixing and pressure broadening partial derivatives
        apply_lineshapemodel_jacobian_scaling(dF, band, i, derivatives_data, derivatives_data_active, T, P, vmrs);
      }

      // Normalize the lines
      switch (band.Normalization()) {
        case Absorption::NormalizationType::None:
          break;
        case Absorption::NormalizationType::VVH:
          apply_VVH_scaling(F, dF, data, f, band.F0(i), T, band, i, derivatives_data, derivatives_data_active);
          break;
        case Absorption::NormalizationType::VVW:
          apply_VVW_scaling(F, dF, f, band.F0(i), band, i, derivatives_data, derivatives_data_active);
          break;
        case Absorption::NormalizationType::RosenkranzQuadratic:
          apply_rosenkranz_quadratic_scaling(F, dF, f, band.F0(i), T, band, i, derivatives_data, derivatives_data_active);
          break;
      }

      // Apply line strength by whatever method is necessary
      switch (band.Population()) {
        case Absorption::PopulationType::ByHITRANFullRelmat:
        case Absorption::PopulationType::ByMakarovFullRelmat:
        case Absorption::PopulationType::ByHITRANRosenkranzRelmat:
        case Absorption::PopulationType::ByLTE:
          apply_linestrength_scaling_by_lte(F, dF, N, dN, band.Line(i), T, band.T0(), isot_ratio, QT, QT0, band, i, derivatives_data, derivatives_data_active, dQTdT);
          break;
        case Absorption::PopulationType::ByNLTEVibrationalTemperatures: {
          auto nlte_data = nlte.get_vibtemp_params(band, i, T);
          apply_linestrength_scaling_by_vibrational_nlte(F, dF, N, dN, band.Line(i), T, band.T0(), nlte_data.T_upp, nlte_data.T_low, nlte_data.E_upp, nlte_data.E_low, isot_ratio, QT, QT0, band, i, derivatives_data, derivatives_data_active, dQTdT);
        } break;
        case Absorption::PopulationType::ByNLTEPopulationDistribution: {
          auto nlte_data = nlte.get_ratio_params(band, i);
          apply_linestrength_from_nlte_level_distributions(F, dF, N, dN, nlte_data.r_low, nlte_data.r_upp, band.g_low(i), band.g_upp(i), band.A(i), band.F0(i), T, band, i, derivatives_data, derivatives_data_active);
        } break;
      }
      
      // Zeeman-adjusted strength
      if (zeeman) {
        F *= Sz;
        N *= Sz;
        dF *= Sz;
        dN *= Sz;
      }
      
      // Sum up the contributions
      sum.F.segment(start, nelem).noalias() += F;
      sum.N.segment(start, nelem).noalias() += N;
      sum.dF.middleRows(start, nelem).noalias() += dF;
      sum.dN.middleRows(start, nelem).noalias() += dN;
    }
  }
  
  // Set negative values to zero incase this is requested
  if (no_negatives) {
    auto reset_zeroes = (sum.F.array().real() < 0);
    
    sum.N = reset_zeroes.select(Complex(0, 0), sum.N);
    for (Index ij=0; ij<nj; ij++)
      sum.dF.col(ij) = reset_zeroes.select(Complex(0, 0), sum.dF.col(ij));
    for (Index ij=0; ij<nj; ij++)
      sum.dN.col(ij) = reset_zeroes.select(Complex(0, 0), sum.dN.col(ij));
    sum.F = reset_zeroes.select(Complex(0, 0), sum.F);
  }
}
