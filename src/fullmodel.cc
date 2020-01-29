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
 * @file   fullmodel.cc
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */

#include "fullmodel.h"

void makarov2020_o2_lines(ArrayOfMatrix& xsec,
                          ArrayOfArrayOfMatrix& dxsec,
                          const Vector& f,
                          const Vector& p,
                          const Vector& t,
                          const makarov2020_o2_lines_control& ctrl)
{
  if (ctrl.pressure >= 0 or 
    ctrl.f0 >= 0 or
    ctrl.intens >= 0 or
    ctrl.a2 >= 0 or
    ctrl.gamma >= 0 or
    ctrl.y0 >= 0 or
    ctrl.y1 >= 0 or
    ctrl.g0 >= 0 or
    ctrl.g1 >= 0 or
    ctrl.dv0 >= 0 or
    ctrl.dv1 >= 0 or
    ctrl.x >= 0) {
    throw std::runtime_error ("Not implemented");
  }
  
  using Constant::pi;
  using Constant::sqrt_pi;
  using Constant::inv_sqrt_pi;
  using Constant::pow2;
  using Constant::pow3;
  using Constant::log10_euler;
  
  constexpr Index n = 38 + 6;
  
  // Central frequency [GHz]
  constexpr std::array<Numeric, n> f0 = {
    118.750334, 56.264774, 62.486253, 58.446588, 60.306056,
    59.590983, 59.164204, 60.434778, 58.323877, 61.150562,
    57.612486, 61.800158, 56.968211, 62.411220, 56.363399,
    62.997984, 55.783815, 63.568526, 55.221384, 64.127775,
    54.671180, 64.678910, 54.130025, 65.224078, 53.595775,
    65.764779, 53.066934, 66.302096, 52.542418, 66.836834,
    52.021429, 67.369601, 51.503360, 67.900868, 50.987745,
    68.431006, 50.474214, 68.960312, 368.498246, 424.763020,
    487.249273, 715.392902, 773.839490, 834.145546};
  
  // Intensity [kHz / kPa]
  constexpr std::array<Numeric, n> intens = {
    940.3, 543.4, 1503.0, 1442.1, 2103.4, 2090.7, 2379.9,
    2438.0, 2363.7, 2479.5, 2120.1, 2275.9, 1746.6, 1915.4,
    1331.8, 1490.2, 945.3, 1078.0, 627.1, 728.7, 389.7, 461.3,
    227.3, 274.0, 124.6, 153.0, 64.29, 80.40, 31.24, 39.80,
    14.32, 18.56, 6.193, 8.172, 2.529, 3.397, 0.975, 1.334,
    67.4, 637.7, 237.4, 98.1, 572.3, 183.1};
  
  // Temperature intensity modifier
  constexpr std::array<Numeric, n> a2 = {
    0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.386,
    0.621, 0.621, 0.910, 0.910, 1.255, 1.255, 1.654, 1.654,
    2.109, 2.108, 2.618, 2.617, 3.182, 3.181, 3.800, 3.800,
    4.474, 4.473, 5.201, 5.200, 5.983, 5.982, 6.819, 6.818,
    7.709, 7.708, 8.653, 8.652, 9.651, 9.650, 0.048, 0.044,
    0.049, 0.145, 0.141, 0.145};
  
  // Pressure broadening coefficient
  constexpr std::array<Numeric, n> gamma = {
    1.685, 1.703, 1.513, 1.495, 1.433, 1.408, 1.353, 1.353,
    1.303, 1.319, 1.262, 1.265, 1.238, 1.217, 1.207, 1.207,
    1.137, 1.137, 1.101, 1.101, 1.037, 1.038, 0.996, 0.996,
    0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811,
    0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 1.64, 1.64, 1.60,
    1.60, 1.62, 1.47};
  
  // First order line mixing first coefficient
  constexpr std::array<Numeric, n> y0 = {
    -0.041, 0.277, -0.372, 0.559, -0.573, 0.618, -0.366, 0.278,
    -0.089, -0.021, 0.060, -0.152, 0.216, -0.293, 0.373, -0.436,
    0.491, -0.542, 0.571, -0.613, 0.636, -0.670, 0.690, -0.718,
    0.740, -0.763, 0.788, -0.807, 0.834, -0.849, 0.876, -0.887,
    0.915, -0.922, 0.950, -0.955, 0.987, -0.988, 0, 0, 0, 0, 0, 0};
  
  // First order line mixing second coefficient
  constexpr std::array<Numeric, n> y1 = {
    0, 0.124, -0.002, 0.008, 0.045, -0.093, 0.264, -0.351, 0.359,
    -0.416, 0.326, -0.353, 0.484, -0.503, 0.579, -0.590, 0.616,
    -0.619, 0.611, -0.609, 0.574, -0.568, 0.574, -0.566, 0.60,
    -0.59, 0.63, -0.62, 0.64, -0.63, 0.65, -0.64, 0.65, -0.64,
    0.65, -0.64, 0.64, -0.62, 0, 0, 0, 0, 0, 0};
  
  // Second order line mixing strength adjustment first coefficient
  constexpr std::array<Numeric, n> g0 = {
    -0.000695, -0.090,  -0.103, -0.239, -0.172, -0.171, 0.028,
    0.150, 0.132, 0.170, 0.087, 0.069, 0.083, 0.067, 0.007,
    0.016, -0.021, -0.066, -0.095, -0.115, -0.118, -0.140,
    -0.173, -0.186, -0.217, -0.227, -0.234, -0.242, -0.266,
    -0.272, -0.301, -0.304, -0.334, -0.333, -0.361, -0.358,
    -0.348, -0.344, 0, 0, 0, 0, 0, 0};
  
  // Second order line mixing strength adjustment second coefficient
  constexpr std::array<Numeric, n> g1 = {
    0, -0.045, 0.007, 0.033, 0.081, 0.162, 0.179, 0.225, 0.054,
    0.003, 0.0004, -0.047, -0.034, -0.071, -0.180, -0.210, -0.285,
    -0.323, -0.363, -0.380, -0.378, -0.387, -0.392, -0.394, -0.424,
    -0.422, -0.465, -0.46, -0.51, -0.50, -0.55, -0.54, -0.58, -0.56,
    -0.62, -0.59, -0.68, -0.65, 0, 0, 0, 0, 0, 0};
  
  // Second order line mixing frequency adjustment first coefficient
  constexpr std::array<Numeric, n> dv0 = {
    -0.00028, 0.00597, -0.0195, 0.032, -0.0475, 0.0541, -0.0232,
    0.0154, 0.0007, -0.0084, -0.0025, -0.0014, -0.0004, -0.0020,
    0.005, -0.0066, 0.0072, -0.008, 0.0064, -0.0070, 0.0056, -0.0060,
    0.0047, -0.0049, 0.0040, -0.0041, 0.0036, -0.0037, 0.0033, -0.0034,
    0.0032, -0.0032, 0.0030, -0.0030, 0.0028, -0.0029, 0.0029, -0.0029,
    0, 0, 0, 0, 0, 0};
  
  // Second order line mixing frequency adjustment second coefficient
  constexpr std::array<Numeric, n> dv1 = {
    -0.00039, 0.009, -0.012, 0.016, -0.027, 0.029, 0.006, -0.015,
    0.010, -0.014, -0.013, 0.013, 0.004, -0.005, 0.010, -0.010, 0.010,
    -0.011, 0.008, -0.009, 0.003, -0.003, 0.0009, -0.0009, 0.0017,
    -0.0016, 0.0024, -0.0023, 0.0024, -0.0024, 0.0024, -0.0020, 0.0017,
    -0.0016, 0.0013, -0.0012, 0.0005, -0.0004, 0, 0, 0, 0, 0, 0};
  
  /*
   Quantum numbers kept so that Zeeman effect can be
   easily implemented in the future
   
  // N of upper level
  constexpr std::array<Index ,n> Np = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 1, 1, 1, 3, 3, 3};
  
  // N of lower level
  constexpr std::array<Index ,n> Npp = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 3, 3, 3, 5, 5, 5};
  
  // J of upper level
  constexpr std::array<Index ,n> Jp = {
    1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
    11, 11, 13, 13, 15, 15, 17, 17,
    19, 19,  21, 21, 23, 23, 25, 25,
    27, 27, 29, 29, 31, 31, 33, 33,
    35, 35, 37, 37, 1, 2, 2, 3, 4, 4};
  
  // J of upper level
  constexpr std::array<Index ,n> Jpp = {
    0, 2, 2, 4, 4, 6, 6, 8, 8, 10,
    10, 12, 12, 14, 14, 16, 16, 18,
    18, 20,  20, 22, 22, 24, 24, 26,
    26, 28, 28, 30, 30, 32, 32, 34,
    34, 36, 36, 38, 2, 2, 3, 4, 4, 5};
  */
  
  // Exponent for temperature
  constexpr Numeric x = 0.754;
  
  // Temperature dependency
  constexpr Numeric t0 = 300;
  
  // Conversion to per meter absorption from per kilometer in decibel
  constexpr Numeric conversion = 0.1820 * 1e-3 / (0.2085 * 10.0 * log10_euler);
  
  // Cancel when there are no lines before any calculations are performed
  if (ctrl.pos_in_xsec < 0) return;
  
  // Per pressure level
  for (Index ip=0; ip<p.nelem(); ip++) {
    // Model implementation
    const Numeric theta = t0 / t[ip];
    const Numeric theta_m1 = theta - 1;
    const Numeric theta_3 = pow3(theta);
    const Numeric theta_x = std::pow(theta, x);
    const Numeric theta_2x = pow2(theta_x);
    const Numeric GD_div_F0 = Linefunctions::DopplerConstant(t[ip], SpeciesTag("O2-66").SpeciesMass());
    
    for (Index i=0; i<n; i++) {
      const Numeric invGD = 1 / (GD_div_F0 * f0[i]*1e9);
      const Numeric fac = sqrt_pi * invGD;
      const Numeric G0 = 1e9 * (1e-5 * p[ip] * gamma[i]) * theta_x;
      const Numeric Y = 1e-5 * p[ip] * (y0[i] + y1[i] * theta_m1) * theta_x;
      const Numeric G = pow2(1e-5 * p[ip]) * (g0[i] + g1[i] * theta_m1) * theta_2x;
      const Numeric DV = 1e9 * pow2(1e-5 * p[ip]) * (dv0[i] + dv1[i] * theta_m1) * theta_2x;
      const Numeric ST = 1e-6 * theta_3 * (1e-3 * p[ip] * intens[i])/(1e9*f0[i]) * std::exp(-a2[i] * theta_m1);
      
      for (Index j=0; j<f.nelem(); j++) {
        const Complex z = Complex(f0[i]*1e9 + DV - f[j], G0) * invGD;
        const Complex Fv = fac * Faddeeva::w(z);
        const Complex Flm = 1 / Complex(G0, f[j] + f0[i]*1e9 + DV);
        
        xsec[ctrl.pos_in_xsec](j, ip) += std::real(conversion * ST * pow2(f[j]) *
        (
          /* around line center */
          Complex(1 + G, Y) * Fv +
          /* mirrored line far from line center */
          Complex(1 + G, -Y) * Flm
        ));
        
        if (ctrl.temperature >= 0) {
          const Numeric dinvGD = - invGD * Linefunctions::dDopplerConstant_dT(t[ip], GD_div_F0);
          const Numeric dG0 = -(x/t[ip]) * G0;
          const Numeric dY = -((y1[i]*t0 + x*(y0[i]*t - y1[i]*(t[ip]-t0)))/(t[ip] * (y0[i]*t[ip] - y1[i]*(t[ip]-t0)))) * Y;
          const Numeric dG = -((g1[i]*t0 + 2*x*(g0[i]*t - g1[i]*(t[ip]-t0)))/(t[ip] * (g0[i]*t[ip] - g1[i]*(t[ip]-t0)))) * G;
          const Numeric dDV = -((dv1[i]*t0 + 2*x*(dv0[i]*t - dv1[i]*(t[ip]-t0)))/(t[ip] * (dv0[i]*t[ip] - dv1[i]*(t[ip]-t0)))) * DV;
          const Numeric dST = (a2[i]*t0 - 3*t[ip]) / pow2(t[ip]) * ST;
          
          const Complex dFv = 2 * (Complex(0, fac * inv_sqrt_pi) - z * Fv) * (invGD * Complex(dDV, dG0) - dinvGD) + Fv * dinvGD;
          const Complex dFlm = - pi * pow2(Flm) * Complex(dG0, dDV);
          dxsec[ctrl.temperature][ctrl.pos_in_xsec](j, ip) += std::real(conversion * ST * pow2(f[j]) * 
          (
            /* around line center */
            Complex(1 + G, Y) * dFv + Complex(dG, dY) * Fv +
            /* mirrored line far from line center */
            Complex(1 + G, -Y) * dFlm + Complex(G, -dY) * Flm
          )) +
          std::real(conversion * dST * pow2(f[j]) *
          (
            /* around line center */
            Complex(1 + G, Y) * Fv +
            /* mirrored line far from line center */
            Complex(1 + G, -Y) * Flm
          ));
        }
        
        if (ctrl.frequency >= 0) {
          const Complex dFv = 2 * (Complex(0, fac * inv_sqrt_pi) - z * Fv) * invGD;
          const Complex dFlm = Complex(0, pi) * pow2(Flm);
          dxsec[ctrl.temperature][ctrl.pos_in_xsec](j, ip) += std::real(conversion * ST * pow2(f[j]) * 
          (
            /* around line center */
            Complex(1 + G, Y) * dFv +
            /* mirrored line far from line center */
            Complex(1 + G, -Y) * dFlm
          )) + 
          std::real(conversion * ST * 2 * f[j] *
          (
            /* around line center */
            Complex(1 + G, Y) * Fv +
            /* mirrored line far from line center */
            Complex(1 + G, -Y) * Flm
          ));
        }
        
      }
    }
  }
}
