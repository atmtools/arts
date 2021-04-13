/* Copyright (C) 2020
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
 * @file   predefined_absorption_models.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */

#include "predefined_absorption_models.h"
#include "lin_alg.h"
#include "linescaling.h"
#include "wigner_functions.h"

constexpr std::size_t nlines_mpm2020 = 44;


constexpr LineShape::SingleSpeciesModel init_mpm2020_slsm(Numeric g00,
                                                          Numeric y0,
                                                          Numeric y1,
                                                          Numeric g0,
                                                          Numeric g1,
                                                          Numeric dv0,
                                                          Numeric dv1,
                                                          Numeric x) noexcept
{
  LineShape::SingleSpeciesModel ssm;
  ssm.G0() = {LineShape::TemperatureModel::T1, g00, x,   NAN, NAN};
  ssm.Y()  = {LineShape::TemperatureModel::T4, y0,  y1,    x, NAN};
  ssm.G()  = {LineShape::TemperatureModel::T4, g0,  g1,  2*x, NAN};
  ssm.DV() = {LineShape::TemperatureModel::T4, dv0, dv1, 2*x, NAN};
  return ssm;
}


constexpr std::array<LineShape::SingleSpeciesModel, nlines_mpm2020> init_mpm2020_lsm() noexcept
{
  // Pressure broadening [1/Pa] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> g00 =
    {1.685E+4, 1.703E+4, 1.513E+4, 1.495E+4, 1.433E+4, 
      1.408E+4, 1.353E+4, 1.353E+4, 1.303E+4, 1.319E+4, 
      1.262E+4, 1.265E+4, 1.238E+4, 1.217E+4, 1.207E+4, 
      1.207E+4, 1.137E+4, 1.137E+4, 1.101E+4, 1.101E+4, 
      1.037E+4, 1.038E+4, 9.96E+3, 9.96E+3, 9.55E+3, 
      9.55E+3, 9.06E+3, 9.06E+3, 8.58E+3, 8.58E+3, 
      8.11E+3, 8.11E+3, 7.64E+3, 7.64E+3, 7.17E+3, 
      7.17E+3, 6.69E+3, 6.69E+3, 1.64E+4, 1.64E+4, 
      1.60E+4, 1.60E+4, 1.62E+4, 1.47E+4,};
  
  // First order line mixing first coefficient [1/Pa] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> y0 = 
    {-4.1E-7, 0.00000277, -0.00000372, 0.00000559, -0.00000573, 
      0.00000618, -0.00000366, 0.00000278, -8.9E-7, -2.1E-7, 
      6.0E-7, -0.00000152, 0.00000216, -0.00000293, 0.00000373, 
      -0.00000436, 0.00000491, -0.00000542, 0.00000571, -0.00000613, 
      0.00000636, -0.00000670, 0.00000690, -0.00000718, 0.00000740, 
      -0.00000763, 0.00000788, -0.00000807, 0.00000834, -0.00000849, 
      0.00000876, -0.00000887, 0.00000915, -0.00000922, 0.00000950, 
      -0.00000955, 0.00000987, -0.00000988, 0.00000, 0.00000, 
      0.00000, 0.00000, 0.00000, 0.00000,};
  
  // First order line mixing second coefficient [1/Pa] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> y1 =
    {0.00000, 0.00000124, -2E-8, 8E-8, 4.5E-7, 
      -9.3E-7, 0.00000264, -0.00000351, 0.00000359, -0.00000416, 
      0.00000326, -0.00000353, 0.00000484, -0.00000503, 0.00000579, 
      -0.00000590, 0.00000616, -0.00000619, 0.00000611, -0.00000609, 
      0.00000574, -0.00000568, 0.00000574, -0.00000566, 0.0000060, 
      -0.0000059, 0.0000063, -0.0000062, 0.0000064, -0.0000063, 
      0.0000065, -0.0000064, 0.0000065, -0.0000064, 0.0000065, 
      -0.0000064, 0.0000064, -0.0000062, 0.00000, 0.00000, 
      0.00000, 0.00000, 0.00000, 0.00000, };
  
  // Second order line mixing strength adjustment first coefficient [1/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> g0 =
    {-6.95E-14, -9.0E-12, -1.03E-11, -2.39E-11, -1.72E-11, 
      -1.71E-11, 2.8E-12, 1.50E-11, 1.32E-11, 1.70E-11, 
      8.7E-12, 6.9E-12, 8.3E-12, 6.7E-12, 7E-13, 
      1.6E-12, -2.1E-12, -6.6E-12, -9.5E-12, -1.15E-11, 
      -1.18E-11, -1.40E-11, -1.73E-11, -1.86E-11, -2.17E-11, 
      -2.27E-11, -2.34E-11, -2.42E-11, -2.66E-11, -2.72E-11, 
      -3.01E-11, -3.04E-11, -3.34E-11, -3.33E-11, -3.61E-11, 
      -3.58E-11, -3.48E-11, -3.44E-11, 0E-10, 0E-10, 
      0E-10, 0E-10, 0E-10, 0E-10,};
  
  // Second order line mixing strength adjustment second coefficient [1/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> g1 =
    {0E-10, -4.5E-12, 7E-13, 3.3E-12, 8.1E-12, 
      1.62E-11, 1.79E-11, 2.25E-11, 5.4E-12, 3E-13, 
      4E-14, -4.7E-12, -3.4E-12, -7.1E-12, -1.80E-11, 
      -2.10E-11, -2.85E-11, -3.23E-11, -3.63E-11, -3.80E-11, 
      -3.78E-11, -3.87E-11, -3.92E-11, -3.94E-11, -4.24E-11, 
      -4.22E-11, -4.65E-11, -4.6E-11, -5.1E-11, -5.0E-11, 
      -5.5E-11, -5.4E-11, -5.8E-11, -5.6E-11, -6.2E-11, 
      -5.9E-11, -6.8E-11, -6.5E-11, 0E-10, 0E-10, 
      0E-10, 0E-10, 0E-10, 0E-10, };
  
  // Second order line mixing frequency adjustment first coefficient [Hz/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> dv0 =
    {-0.000028, 0.000597, -0.00195, 0.0032, -0.00475, 
      0.00541, -0.00232, 0.00154, 0.00007, -0.00084, 
      -0.00025, -0.00014, -0.00004, -0.00020, 0.0005, 
      -0.00066, 0.00072, -0.0008, 0.00064, -0.00070, 
      0.00056, -0.00060, 0.00047, -0.00049, 0.00040, 
      -0.00041, 0.00036, -0.00037, 0.00033, -0.00034, 
      0.00032, -0.00032, 0.00030, -0.00030, 0.00028, 
      -0.00029, 0.00029, -0.00029, 0.0, 0.0, 
      0.0, 0.0, 0.0, 0.0, };
  
  // Second order line mixing frequency adjustment second coefficient [Hz/Pa^2] at reference temperature
  constexpr std::array<Numeric, nlines_mpm2020> dv1 =
    {-0.000039, 0.0009, -0.0012, 0.0016, -0.0027, 
      0.0029, 0.0006, -0.0015, 0.0010, -0.0014, 
      -0.0013, 0.0013, 0.0004, -0.0005, 0.0010, 
      -0.0010, 0.0010, -0.0011, 0.0008, -0.0009, 
      0.0003, -0.0003, 0.00009, -0.00009, 0.00017, 
      -0.00016, 0.00024, -0.00023, 0.00024, -0.00024, 
      0.00024, -0.00020, 0.00017, -0.00016, 0.00013, 
      -0.00012, 0.00005, -0.00004, 0.0, 0.0, 
      0.0, 0.0, 0.0, 0.0, };
  
  // Temperature scaling exponent [-]
  constexpr Numeric x = 0.754;
  
  // Init all the values
  std::array<LineShape::SingleSpeciesModel, nlines_mpm2020> out;
  for (std::size_t i=0; i<nlines_mpm2020; i++) {
    out[i] = init_mpm2020_slsm(g00[i], y0[i], y1[i], g0[i], g1[i], dv0[i], dv1[i], x);
  }
  return out;
}


constexpr QuantumIdentifier init_mpm2020_qid(Index species, Index isot, Rational Jup, Rational Jlow, Rational Nup, Rational Nlow) noexcept
{
  QuantumNumbers upp;
  QuantumNumbers low;
  upp[QuantumNumberType::J] = Jup;
  upp[QuantumNumberType::N] = Nup;
  upp[QuantumNumberType::v1] = 0;
  low[QuantumNumberType::J] = Jlow;
  low[QuantumNumberType::N] = Nlow;
  low[QuantumNumberType::v1] = 0;
  return QuantumIdentifier(species, isot, upp, low);
}


constexpr std::array<QuantumIdentifier, nlines_mpm2020> init_mpm2020_qids(const Index& species, const Index& isot) noexcept
{
  // N of upper level
  constexpr std::array<Index, nlines_mpm2020> Np = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 1, 1, 1, 3, 3, 3};
  
  // N of lower level
  constexpr std::array<Index, nlines_mpm2020> Npp = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 3, 3, 3, 5, 5, 5};
  
  // J of upper level
  constexpr std::array<Index, nlines_mpm2020> Jp = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 1, 2, 2, 3, 4, 4};
  
  // J of lower level
  constexpr std::array<Index, nlines_mpm2020> Jpp = {
    0, 2, 2, 4, 4, 6, 6, 8, 8, 10,
    10, 12, 12, 14, 14, 16, 16, 18,
    18, 20,  20, 22, 22, 24, 24, 26,
    26, 28, 28, 30, 30, 32, 32, 34,
    34, 36, 36, 38, 2, 2, 3, 4, 4, 5};
  
  // Init all the values
  std::array<QuantumIdentifier, nlines_mpm2020> out;
  for (std::size_t i=0; i<nlines_mpm2020; i++) {
    out[i] = init_mpm2020_qid(species, isot, Rational(Jp[i]), Rational(Np[i]), Rational(Jpp[i]), Rational(Npp[i]));
  }
  return out;
}


void Absorption::PredefinedModel::makarov2020_o2_lines_mpm(Matrix& xsec,
                                                           ArrayOfMatrix& dxsec,
                                                           const Vector& f,
                                                           const Vector& p,
                                                           const Vector& t,
                                                           const Vector& water_vmr,
                                                           const ArrayOfRetrievalQuantity& jacs)
{
  using Constant::pi;
  using Constant::sqrt_pi;
  using Constant::inv_sqrt_pi;
  using Constant::pow2;
  using Constant::pow3;
  
  // Central frequency [Hz]
  constexpr std::array<Numeric, nlines_mpm2020> f0 = {
    1.18750334E+11, 5.6264774E+10, 6.2486253E+10, 5.8446588E+10, 6.0306056E+10, 
    5.9590983E+10, 5.9164204E+10, 6.0434778E+10, 5.8323877E+10, 6.1150562E+10, 
    5.7612486E+10, 6.1800158E+10, 5.6968211E+10, 6.2411220E+10, 5.6363399E+10, 
    6.2997984E+10, 5.5783815E+10, 6.3568526E+10, 5.5221384E+10, 6.4127775E+10, 
    5.4671180E+10, 6.4678910E+10, 5.4130025E+10, 6.5224078E+10, 5.3595775E+10, 
    6.5764779E+10, 5.3066934E+10, 6.6302096E+10, 5.2542418E+10, 6.6836834E+10, 
    5.2021429E+10, 6.7369601E+10, 5.1503360E+10, 6.7900868E+10, 5.0987745E+10, 
    6.8431006E+10, 5.0474214E+10, 6.8960312E+10, 3.68498246E+11, 4.24763020E+11, 
    4.87249273E+11, 7.15392902E+11, 7.73839490E+11, 8.34145546E+11, };
  
  // Intensity [1 / Pa] (rounded to 10 digits because at most 9 digits exist in f0)
  constexpr std::array<Numeric, nlines_mpm2020> intens = {
    1.591521878E-21, 1.941172240E-21, 4.834543970E-21, 4.959264029E-21, 7.010386457E-21, 
    7.051673348E-21, 8.085012578E-21, 8.108262250E-21, 8.145673278E-21, 8.149757320E-21, 
    7.396406085E-21, 7.401923754E-21, 6.162286575E-21, 6.168475265E-21, 4.749226167E-21, 
    4.754435107E-21, 3.405982896E-21, 3.408455562E-21, 2.282498656E-21, 2.283934341E-21, 
    1.432692459E-21, 1.433513473E-21, 8.439995690E-22, 8.443521837E-22, 4.672706507E-22, 
    4.676049313E-22, 2.435008301E-22, 2.437304596E-22, 1.195038747E-22, 1.196873412E-22, 
    5.532759045E-23, 5.537261239E-23, 2.416832398E-23, 2.418989865E-23, 9.969285671E-24, 
    9.977543709E-24, 3.882541154E-24, 3.888101811E-24, 3.676253816E-23, 3.017524005E-22, 
    9.792882227E-23, 2.756166168E-23, 1.486462215E-22, 4.411918954E-23, };
  
  // Temperature intensity modifier
  constexpr std::array<Numeric, nlines_mpm2020> a2 = {
    0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.386,
    0.621, 0.621, 0.910, 0.910, 1.255, 1.255, 1.654, 1.654,
    2.109, 2.108, 2.618, 2.617, 3.182, 3.181, 3.800, 3.800,
    4.474, 4.473, 5.201, 5.200, 5.983, 5.982, 6.819, 6.818,
    7.709, 7.708, 8.653, 8.652, 9.651, 9.650, 0.048, 0.044,
    0.049, 0.145, 0.141, 0.145};
  
  // Line shape model in SI units
  constexpr auto lsm = init_mpm2020_lsm();
  
  // Reference temperature [K]
  constexpr Numeric t0 = 300;
  
  // QuantumIdentifier if we need it
  auto species = SpeciesTag("O2-66");
  const std::array<QuantumIdentifier, nlines_mpm2020> qids = init_mpm2020_qids(species.Species(), species.Isotopologue());
  
  // Model setting
  const bool do_temp_deriv = do_temperature_jacobian(jacs);
  
  // Per pressure level
  for (Index ip=0; ip<p.nelem(); ip++) {
    const Numeric theta = t0 / t[ip];
    const Numeric theta_m1 = theta - 1;
    const Numeric theta_3 = pow3(theta);
    const Numeric GD_div_F0 = Linefunctions::DopplerConstant(t[ip], species.SpeciesMass());
    
    for (std::size_t i=0; i<nlines_mpm2020; i++) {
      const Numeric invGD = 1 / (GD_div_F0 * f0[i]);
      const Numeric fac = sqrt_pi * invGD;
      const Numeric ST = theta_3 * p[ip] * intens[i] * std::exp(-a2[i] * theta_m1);
      const Numeric G0 = (1 + 0.1*water_vmr[ip]) * p[ip] * lsm[i].compute(t[ip], t0, LineShape::Variable::G0);
      const Numeric Y = p[ip] * lsm[i].compute(t[ip], t0, LineShape::Variable::Y);
      const Numeric G = pow2( p[ip]) * lsm[i].compute(t[ip], t0, LineShape::Variable::G);
      const Numeric DV = pow2(p[ip]) * lsm[i].compute(t[ip], t0, LineShape::Variable::DV);
      
      const Numeric dinvGD_dT = do_temp_deriv ? - invGD * Linefunctions::dDopplerConstant_dT(t[ip], GD_div_F0) : 0;
      const Numeric dST_dT = do_temp_deriv ? (a2[i]*t0 - 3*t[ip]) / pow2(t[ip]) * ST : 0;
      const Numeric dG0_dT = do_temp_deriv ? (1 + 0.1*water_vmr[ip]) * p[ip] * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::G0) : 0;
      const Numeric dY_dT = do_temp_deriv ? p[ip] * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::Y) : 0;
      const Numeric dG_dT = do_temp_deriv ? pow2(p[ip]) * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::G) : 0;
      const Numeric dDV_dT = do_temp_deriv ? pow2(p[ip]) * lsm[i].compute_dT(t[ip], t0, LineShape::Variable::DV) : 0;
      
      for (Index j=0; j<f.nelem(); j++) {
        const Complex z = Complex(f0[i] + DV - f[j], G0) * invGD;
        const Complex Fv = fac * Faddeeva::w(z);
        const Complex Flm = 1 / Complex(G0, f[j] + f0[i] + DV);
        
        const Complex abs = std::real(
          Complex(1 + G, Y) * Fv +
          Complex(1 + G, -Y) * Flm);
        
        xsec(j, ip) += ST * pow2(f[j]) * abs.real();
        
        if (jacs.nelem()) {
          const Complex dw = 2 * (Complex(0, fac * inv_sqrt_pi) - z * Fv);
          const Complex dm = - pi * pow2(Flm);
          
          for (Index iq=0; iq<jacs.nelem(); iq++) {
            const auto& deriv = jacs[iq];
            
            if (not propmattype(deriv)) continue;
            
            
            if (deriv == Jacobian::Atm::Temperature) {
              const Complex dFv = dw * (invGD * Complex(dDV_dT, dG0_dT) - dinvGD_dT) + Fv * dinvGD_dT;
              const Complex dFlm = dm * Complex(dG0_dT, dDV_dT);
              dxsec[iq](j, ip) += pow2(f[j]) * (ST * std::real(
                Complex(1 + G, Y) * dFv + Complex(dG_dT, dY_dT) * Fv +
                Complex(1 + G, -Y) * dFlm + Complex(G, -dY_dT) * Flm) + abs.real() * dST_dT);
            } else if (is_frequency_parameter(deriv)) {
              const Complex dFv = - dw * invGD;
              const Complex dFlm = Complex(0, 1) * dm;
              dxsec[iq](j, ip) += ST * (pow2(f[j]) * std::real(
                Complex(1 + G, Y) * dFv +
                Complex(1 + G, -Y) * dFlm) + 2 * abs.real() * f[j]);
            } else if (deriv.Target().needQuantumIdentity()) {
              const Absorption::QuantumIdentifierLineTarget lt(deriv.QuantumIdentity(), qids[i]);
              
              //NOTE: This is a special case where each line must be seen as a "Band" by themselves.
              //NOTE: (cont) This is because we never check for "Line" unless a full Absorption::Lines
              //NOTE: (cont) is used in the QuantumIdentifierLineTarget struct.
              if (lt not_eq Absorption::QuantumIdentifierLineTargetType::Band) continue;
              
              if (deriv == Jacobian::Line::ShapeG0X0) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  Complex(1 + G, Y) * Complex(0, 1) * dw * invGD +
                  Complex(1 + G, -Y) * dm) * 
                  lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::G0);
              } else if (deriv == Jacobian::Line::ShapeG0X1) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  Complex(1 + G, Y) * Complex(0, 1) * dw * invGD +
                  Complex(1 + G, -Y) * dm) * 
                  lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::DV);
              } else if (deriv == Jacobian::Line::ShapeDVX0) {
                const Complex dFv = dw * invGD;
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  Complex(1 + G, Y) * dFv +
                  Complex(1 + G, -Y) * dFlm) * 
                  lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::DV);
              } else if (deriv == Jacobian::Line::ShapeDVX1) {
                const Complex dFv = dw * invGD;
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  Complex(1 + G, Y) * dFv +
                  Complex(1 + G, -Y) * dFlm) * 
                  lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::DV);
              } else if (deriv == Jacobian::Line::ShapeDVX2) {
                const Complex dFv = dw * invGD;
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  Complex(1 + G, Y) * dFv +
                  Complex(1 + G, -Y) * dFlm) * 
                  lsm[i].compute_dX2(t[ip], t0, LineShape::Variable::DV);
              } else if (deriv == Jacobian::Line::ShapeGX0) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(Fv + Flm) * 
                  lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::G);
              } else if (deriv == Jacobian::Line::ShapeYX0) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(Fv + Flm) * 
                  lsm[i].compute_dX0(t[ip], t0, LineShape::Variable::Y);
              } else if (deriv == Jacobian::Line::ShapeGX1) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(Fv - Flm) * 
                  lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::G);
              } else if (deriv == Jacobian::Line::ShapeYX1) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(Fv + Flm) * 
                  lsm[i].compute_dX1(t[ip], t0, LineShape::Variable::Y);
              } else if (deriv == Jacobian::Line::ShapeGX2) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(Fv - Flm) * 
                  lsm[i].compute_dX2(t[ip], t0, LineShape::Variable::G);
              } else if (deriv == Jacobian::Line::ShapeYX2) {
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(Fv - Flm) * 
                  lsm[i].compute_dX2(t[ip], t0, LineShape::Variable::Y);
              } else if (deriv == Jacobian::Line::Center) {
                const Complex dFv = Fv / f0[i] - dw * invGD + dw * z / f0[i];
                const Complex dFlm = Complex(0, 1) * dm;
                dxsec[iq](j, ip) += ST * pow2(f[j]) * std::real(
                  Complex(1 + G, Y) * dFv +
                  Complex(1 + G, -Y) * dFlm);
              } else if (deriv == Jacobian::Line::Strength) {
                dxsec[iq](j, ip) += theta_3 * p[ip] * std::exp(-a2[i] * theta_m1) * pow2(f[j]) * abs.real();
              } else if (deriv == Jacobian::Line::SpecialParameter1) {
                dxsec[iq](j, ip) += -theta_m1 * ST * pow2(f[j]) * abs.real();
              }
            }
          }
        }
      }
    }
  }
}
