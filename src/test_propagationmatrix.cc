/* Copyright (C) 2015 Richard Larsson <ric.larsson@gmail.com>
 * 
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation; either version 2, or (at your option) any
 *  later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 *  USA. */

/*!
 * \file   test_propagationmatrix.cc
 * \author <ric.larsson@gmail.com>
 * \date   2015-10-20
 * 
 * \brief  Test Propagation Matrix Internal Partial Derivatives and PropagationMatrix
 */

#include <random>
#include "absorption.h"
#include "arts.h"
#include "global_data.h"
#include "lineshapemodel.h"
#include "linefunctions.h"
#include "linescaling.h"
#include "transmissionmatrix.h"
#include "zeeman.h"
#include "zeemandata.h"
#include <Faddeeva/Faddeeva.hh>
#include "continua.h"

// void define_species_data();
// void define_species_map();

void test_matrix_buildup() {
  const Numeric k11 = 1;
  const Numeric k12 = -0.51;
  const Numeric k13 = -0.21;
  const Numeric k14 = 0.31;
  const Numeric k23 = -0.1;
  const Numeric k24 = -0.99;
  const Numeric k34 = 2;

  const Numeric r = 0.5;

  const Numeric a = -k11 * r;
  const Numeric b = -k12 * r;
  const Numeric c = -k13 * r;
  const Numeric d = -k14 * r;
  const Numeric u = -k23 * r;
  const Numeric v = -k24 * r;
  const Numeric w = -k34 * r;

  const Numeric b2 = b * b, c2 = c * c, d2 = d * d, u2 = u * u, v2 = v * v,
                w2 = w * w;

  const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;

  Numeric Const1;
  Const1 = b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2);
  Const1 += c2 * (c2 * 0.5 + d2 - u2 + v2 - w2);
  Const1 += d2 * (d2 * 0.5 + u2 - v2 - w2);
  Const1 += u2 * (u2 * 0.5 + v2 + w2);
  Const1 += v2 * (v2 * 0.5 + w2);
  Const1 *= 2;
  Const1 += 8 * (b * d * u * w - b * c * v * w - c * d * u * v);
  Const1 += w2 * w2;

  if (Const1 > 0.0)
    Const1 = sqrt(Const1);
  else
    Const1 = 0.0;

  const Complex sqrt_BpA = sqrt(Complex(Const2 + Const1, 0.0));
  const Complex sqrt_BmA = sqrt(Complex(Const2 - Const1, 0.0));
  const Numeric x = sqrt_BpA.real() * sqrt(0.5);
  const Numeric y = sqrt_BmA.imag() * sqrt(0.5);
  const Numeric x2 = x * x;
  const Numeric y2 = y * y;
  const Numeric cos_y = cos(y);
  const Numeric sin_y = sin(y);
  const Numeric cosh_x = cosh(x);
  const Numeric sinh_x = sinh(x);
  const Numeric x2y2 = x2 + y2;
  const Numeric inv_x2y2 = 1.0 / x2y2;

  std::cout << x << " " << y << " " << Const1 << " " << Const2 << "\n";

  Numeric C0, C1, C2, C3;
  Numeric inv_y = 0.0, inv_x = 0.0;  // Init'd to remove warnings

  // X and Y cannot both be zero
  if (x == 0.0) {
    inv_y = 1.0 / y;
    C0 = 1.0;
    C1 = 1.0;
    C2 = (1.0 - cos_y) * inv_x2y2;
    C3 = (1.0 - sin_y * inv_y) * inv_x2y2;
  } else if (y == 0.0) {
    inv_x = 1.0 / x;
    C0 = 1.0;
    C1 = 1.0;
    C2 = (cosh_x - 1.0) * inv_x2y2;
    C3 = (sinh_x * inv_x - 1.0) * inv_x2y2;
  } else {
    inv_x = 1.0 / x;
    inv_y = 1.0 / y;

    C0 = (cos_y * x2 + cosh_x * y2) * inv_x2y2;
    C1 = (sin_y * x2 * inv_y + sinh_x * y2 * inv_x) * inv_x2y2;
    C2 = (cosh_x - cos_y) * inv_x2y2;
    C3 = (sinh_x * inv_x - sin_y * inv_y) * inv_x2y2;
  }

  std::cout << C0 << " " << C1 << " " << C2 << " " << C3 << "\n";

  Matrix F(4, 4, 0), A(4, 4, 0);

  MatrixViewMap eigF = MapToEigen(F);
  Eigen::Matrix4d eigA;
  eigA << 0, b, c, d, b, 0, u, v, c, -u, 0, w, d, -v, -w, 0;

  eigF = C1 * eigA + C2 * eigA * eigA + C3 * eigA * eigA * eigA;
  eigF(0, 0) += C0;
  eigF(1, 1) += C0;
  eigF(2, 2) += C0;
  eigF(3, 3) += C0;
  eigF *= exp(a);

  std::cout << F << "\n";
}

void test_transmissionmatrix() {
  // Initializes as unity matrices
  TransmissionMatrix a(2, 4);
  std::cout << "Initialized TransmissionMatrix(2, 4):\n" << a << "\n";

  // Set a single input
  Eigen::Matrix4d A;
  A << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
  std::cout << "New Matrix:\n" << A << "\n\n";
  a.Mat4(0) = A;
  std::cout << "Updated TransmissionMatrix Position 1 wit New Matrix:\n"
            << a << "\n";

  // The stream can also set the values
  String S =
      "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 125 26 27 28 29 30 31 32";
  std::cout << "Stream:\n" << S << "\n\n";
  std::istringstream astream(S);
  astream >> a;
  std::cout << "Streamed into TransmissionMatrix:\n" << a << "\n";

  // Initialize empty
  RadiationVector b(3, 3);
  std::cout << "Initialized RadiationVector(3, 3)\n" << b << "\n";

  // Set is not defined but add is
  Eigen::Vector3d B;
  B << 1, 2, 3;
  std::cout << "New Vector:\n" << B << "\n\n";  // nb. not transposed
  b.Vec3(1).noalias() += B;
  std::cout << "Updated RadiationVector Position 1 with New Vector:\n"
            << b << "\n";

  // The stream can also set the values
  String T = "1 2 3 4 5 6 7 8 90";
  std::cout << "Stream:\n" << T << "\n\n";
  std::istringstream bstream(T);
  bstream >> b;
  std::cout << "Streamed into RadiationVector:\n" << b << "\n";
}

void test_r_deriv_propagationmatrix() {
  auto f = [](Numeric x) { return 0.1 * x; };
  auto df = [](Numeric x) { return 0.1 + 0 * x; };

  Index nstokes = 4;
  const Numeric x1 = 30;
  const Numeric x2 = -0.1;
  const Numeric r_normal = 1000;
  const Numeric r_extra1 = r_normal + f(x1);
  const Numeric r_extra2 = r_normal + f(x2);

  PropagationMatrix a(1, nstokes);
  a.Kjj() = 10 + Numeric(rand() % 1000) / 100;
  if (nstokes > 1) {
    a.K12() = 2 + Numeric(rand() % 1000) / 100;
    if (nstokes > 2) {
      a.K13() = 3 + Numeric(rand() % 1000) / 100;
      a.K23() = -1 + Numeric(rand() % 1000) / 100;
      if (nstokes > 3) {
        a.K14() = 5 + Numeric(rand() % 1000) / 100;
        a.K24() = -3 + Numeric(rand() % 1000) / 100;
        a.K34() = -2 + Numeric(rand() % 1000) / 100;
      }
    }
  }
  a.GetData() *= 1e-5;

  PropagationMatrix b(1, nstokes);
  b.Kjj() = 5 + Numeric(rand() % 1000) / 100;
  if (nstokes > 1) {
    b.K12() = -1 + Numeric(rand() % 1000) / 100;
    if (nstokes > 2) {
      b.K13() = -3 + Numeric(rand() % 1000) / 100;
      b.K23() = 4 + Numeric(rand() % 1000) / 100;
      if (nstokes > 3) {
        b.K14() = 2 + Numeric(rand() % 1000) / 100;
        b.K24() = -1 + Numeric(rand() % 1000) / 100;
        b.K34() = 3 + Numeric(rand() % 1000) / 100;
      }
    }
  }
  b.GetData() *= 5e-6;

  ArrayOfPropagationMatrix da(1, PropagationMatrix(1, nstokes));
  da[0].GetData() = 0;
  Tensor3 T_normal(1, nstokes, nstokes), T_extra(1, nstokes, nstokes);
  Tensor4 dT1(1, 1, nstokes, nstokes), dT2(1, 1, nstokes, nstokes);

  compute_transmission_matrix_and_derivative(
      T_normal, dT1, dT2, r_normal, a, b, da, da, df(x1), df(x2), 0);

  std::cout << "Transmission at r=" << r_normal << ":\n"
            << MapToEigen(T_normal(0, joker, joker)) << "\n"
            << "\n";
  std::cout << "First derivative:\n"
            << MapToEigen(dT1(0, 0, joker, joker)) << "\n"
            << "\n";
  std::cout << "Second derivative:\n"
            << MapToEigen(dT2(0, 0, joker, joker)) << "\n"
            << "\n";

  compute_transmission_matrix(T_extra, r_extra1, a, b);

  std::cout << "Transmission at perturbed r1=" << r_extra1 << ":\n"
            << MapToEigen(T_extra(0, joker, joker)) << "\n"
            << "\n";
  T_extra -= T_normal;
  T_extra /= x1;
  std::cout << "First derivative perturbed:\n"
            << MapToEigen(T_extra(0, joker, joker)) << "\n"
            << "\n";
  T_extra /= dT1(0, joker, joker, joker);
  std::cout << "First derivative perturbed relative:\n"
            << MapToEigen(T_extra(0, joker, joker)) << "\n"
            << "\n";

  compute_transmission_matrix(T_extra, r_extra2, a, b);

  std::cout << "Transmission at perturbed r2=" << r_extra2 << ":\n"
            << MapToEigen(T_extra(0, joker, joker)) << "\n"
            << "\n";
  T_extra -= T_normal;
  T_extra /= x2;
  std::cout << "Second derivative perturbed:\n"
            << MapToEigen(T_extra(0, joker, joker)) << "\n"
            << "\n";
  T_extra /= dT2(0, joker, joker, joker);
  std::cout << "Second derivative perturbed relative:\n"
            << MapToEigen(T_extra(0, joker, joker)) << "\n"
            << "\n";
}

void test_transmat_from_propmat() {
  const Numeric a = 2;
  const Numeric b = 3;
  const Numeric c = 4;
  const Numeric d = 1;
  const Numeric u = 5;
  const Numeric v = 1;
  const Numeric w = 5;

  PropagationMatrix test1(1, 1);
  PropagationMatrix test2(1, 2);
  PropagationMatrix test3(1, 3);
  PropagationMatrix test4(1, 4);
  test1.Kjj() = a;
  test2.Kjj() = a;
  test3.Kjj() = a;
  test4.Kjj() = a;
  test2.K12() = b;
  test3.K12() = b;
  test4.K12() = b;
  test3.K13() = c;
  test4.K13() = c;
  test3.K23() = u;
  test4.K23() = u;
  test4.K14() = d;
  test4.K24() = v;
  test4.K34() = w;

  std::cout << test1 << "\n\n"
            << test2 << "\n\n"
            << test3 << "\n\n"
            << test4 << "\n\n";

  TransmissionMatrix ans1(1, 1);
  TransmissionMatrix ans2(1, 2);
  TransmissionMatrix ans3(1, 3);
  TransmissionMatrix ans4(1, 4);

  ArrayOfTransmissionMatrix empty(0);

  stepwise_transmission(ans1,
                        empty,
                        empty,
                        test1,
                        test1,
                        ArrayOfPropagationMatrix(0),
                        ArrayOfPropagationMatrix(0),
                        1,
                        0,
                        0,
                        -1);

  stepwise_transmission(ans2,
                        empty,
                        empty,
                        test2,
                        test2,
                        ArrayOfPropagationMatrix(0),
                        ArrayOfPropagationMatrix(0),
                        1,
                        0,
                        0,
                        -1);

  stepwise_transmission(ans3,
                        empty,
                        empty,
                        test3,
                        test3,
                        ArrayOfPropagationMatrix(0),
                        ArrayOfPropagationMatrix(0),
                        1,
                        0,
                        0,
                        -1);

  stepwise_transmission(ans4,
                        empty,
                        empty,
                        test4,
                        test4,
                        ArrayOfPropagationMatrix(0),
                        ArrayOfPropagationMatrix(0),
                        1,
                        0,
                        0,
                        -1);

  std::cout << ans1 << "\n\n"
            << ans2 << "\n\n"
            << ans3 << "\n\n"
            << ans4 << "\n\n";
}

void test_transmat_to_cumulativetransmat() {
  int i = 1;
  ArrayOfPropagationMatrix propmats(5, PropagationMatrix(1, 4));
  for (auto& pm : propmats) {
    pm.K12() = i;
    i++;
    pm.K13() = i;
    i++;
    pm.K23() = i;
    i++;
    pm.K14() = i;
    i++;
    pm.K24() = i;
    i++;
    pm.K34() = i;
    i++;
    pm.Kjj() = 2 * i;
    i++;
    pm.GetData() *= 1e-2;
  }

  std::cout << "Propmats:\n" << propmats << "\n\n";

  ArrayOfTransmissionMatrix layers(5, TransmissionMatrix(1, 4));
  ArrayOfTransmissionMatrix empty(0);
  for (i = 0; i < 4; i++) {
    stepwise_transmission(layers[i + 1],
                          empty,
                          empty,
                          propmats[i],
                          propmats[i + 1],
                          ArrayOfPropagationMatrix(0),
                          ArrayOfPropagationMatrix(0),
                          1,
                          0,
                          0,
                          -1);
  }

  std::cout << "Layers:\n" << layers << "\n\n";

  ArrayOfTransmissionMatrix cumulative_forward =
      cumulative_transmission(layers, CumulativeTransmission::Forward);
  ArrayOfTransmissionMatrix cumulative_reflect =
      cumulative_transmission(layers, CumulativeTransmission::Reverse);

  std::cout << "Forward accumulation:\n" << cumulative_forward << "\n\n";
  std::cout << "Reflect accumulation:\n" << cumulative_reflect << "\n\n";
}

void test_sinc_likes_0limit() {
  Numeric start = 1.0;
  Numeric end = 1e-7;
  int n = 10000;

  Vector x(n), sx(n), shx(n), cx(n), chx(n);
  nlogspace(x, start, end, n);

  for (int i = 0; i < n; i++) {
    sx[i] = std::sin(x[i]);
    cx[i] = std::cos(x[i]);
    shx[i] = std::sinh(x[i]);
    chx[i] = std::cosh(x[i]);
  }

  std::cout
      << std::scientific << std::setprecision(15)
      << "x\tabs(sx/x-1)\tabs(shx/x-1)\tabs((chx-cx)/2x^2-1/2)\tabs((shx/x-sx/x)/2x^2-1/6)\n";

  for (int i = 0; i < n; i++)
    std::cout << x[i] << '\t' << std::abs(sx[i] / x[i] - 1.0) << '\t'
              << std::abs(shx[i] / x[i] - 1.0) << '\t'
              << std::abs((chx[i] - cx[i]) / (2 * x[i] * x[i]) - 0.5) << '\t'
              << std::abs((shx[i] / x[i] - sx[i] / x[i]) / (2 * x[i] * x[i]) -
                          1.0 / 6.0)
              << '\n';
}

void test_zeeman() {
  define_species_data();
  define_species_map();

  auto o266 = SpeciesTag("O2-66");
  auto o268 = SpeciesTag("O2-68");

  Numeric g;
  QuantumNumbers qn;
  qn.Set(QuantumNumberType::Hund, Index(Hund::CaseB));
  qn.Set(QuantumNumberType::Lambda, 0);
  qn.Set(QuantumNumberType::v1, 0);
  qn.Set(QuantumNumberType::S, 1);

  std::cout << "Table from Larsson, Lankhaar, Eriksson (2019)\n";
  for (Index i = 1; i < 51; i++) {
    qn.Set(QuantumNumberType::J, i);

    qn.Set(QuantumNumberType::N, qn[QuantumNumberType::J] - 1);
    std::cout << i << "_" << i - 1;

    g = Zeeman::GetAdvancedModel(
            QuantumIdentifier(o266.Species(), o266.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    g = Zeeman::GetAdvancedModel(
            QuantumIdentifier(o268.Species(), o268.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    g = Zeeman::GetSimpleModel(
            QuantumIdentifier(o266.Species(), o266.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    qn.Set(QuantumNumberType::N, qn[QuantumNumberType::J]);
    std::cout << '\t' << i << "_" << i;

    g = Zeeman::GetAdvancedModel(
            QuantumIdentifier(o266.Species(), o266.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    g = Zeeman::GetAdvancedModel(
            QuantumIdentifier(o268.Species(), o268.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    g = Zeeman::GetSimpleModel(
            QuantumIdentifier(o266.Species(), o266.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    qn.Set(QuantumNumberType::N, qn[QuantumNumberType::J] + 1);
    std::cout << '\t' << i << "_" << i + 1;

    g = Zeeman::GetAdvancedModel(
            QuantumIdentifier(o266.Species(), o266.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    g = Zeeman::GetAdvancedModel(
            QuantumIdentifier(o268.Species(), o268.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    g = Zeeman::GetSimpleModel(
            QuantumIdentifier(o266.Species(), o266.Isotopologue(), qn, qn))
            .gl();
    std::cout << '\t' << g;

    std::cout << '\n';
  }
}

constexpr bool test_quantum_numbers(const QuantumNumbers qns, const Index i)
{
  return (i > 0) ? (qns[i].isUndefined() ? test_quantum_numbers(qns, i-1) : false) : true;
}

void test_quantum()
{
  static_assert(test_quantum_numbers(QuantumNumbers(), Index(QuantumNumberType::FINAL_ENTRY) - 1),
                "Bad last entry in QuantumNumbers.  Did you recently expand the list?");
}


void test_mpm20()
{
//   using Constant::inv_sqrt_pi;
  using Constant::pow2;
  using Constant::pow3;
  using Constant::log10_euler;
  
  define_species_data();
  define_species_map();
  
  constexpr Index n = 38 + 6;
  
  // Central frequency [GHz]
  constexpr std::array<Numeric, n> f0 = {118.750334, 56.264774, 62.486253, 58.446588, 60.306056,
                                         59.590983, 59.164204, 60.434778, 58.323877, 61.150562,
                                         57.612486, 61.800158, 56.968211, 62.411220, 56.363399,
                                         62.997984, 55.783815, 63.568526, 55.221384, 64.127775,
                                         54.671180, 64.678910, 54.130025, 65.224078, 53.595775,
                                         65.764779, 53.066934, 66.302096, 52.542418, 66.836834,
                                         52.021429, 67.369601, 51.503360, 67.900868, 50.987745,
                                         68.431006, 50.474214, 68.960312, 368.498246, 424.763020,
                                         487.249273, 715.392902, 773.839490, 834.145546};
  
  // Intensity [kHz / kPa]
  constexpr std::array<Numeric, n> intens = {940.3, 543.4, 1503.0, 1442.1, 2103.4, 2090.7, 2379.9,
                                             2438.0, 2363.7, 2479.5, 2120.1, 2275.9, 1746.6, 1915.4,
                                             1331.8, 1490.2, 945.3, 1078.0, 627.1, 728.7, 389.7, 461.3,
                                             227.3, 274.0, 124.6, 153.0, 64.29, 80.40, 31.24, 39.80,
                                             14.32, 18.56, 6.193, 8.172, 2.529, 3.397, 0.975, 1.334,
                                             67.4, 637.7, 237.4, 98.1, 572.3, 183.1};
  
  // Temperature intensity modifier
  constexpr std::array<Numeric, n> a2 = {0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.386,
                                         0.621, 0.621, 0.910, 0.910, 1.255, 1.255, 1.654, 1.654,
                                         2.109, 2.108, 2.618, 2.617, 3.182, 3.181, 3.800, 3.800,
                                         4.474, 4.473, 5.201, 5.200, 5.983, 5.982, 6.819, 6.818,
                                         7.709, 7.708, 8.653, 8.652, 9.651, 9.650, 0.048, 0.044,
                                         0.049, 0.145, 0.141, 0.145};
  
  // Pressure broadening coefficient
  constexpr std::array<Numeric, n> gamma = {1.685, 1.703, 1.513, 1.495, 1.433, 1.408, 1.353, 1.353,
                                            1.303, 1.319, 1.262, 1.265, 1.238, 1.217, 1.207, 1.207,
                                            1.137, 1.137, 1.101, 1.101, 1.037, 1.038, 0.996, 0.996,
                                            0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811,
                                            0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 1.64, 1.64, 1.60,
                                            1.60, 1.62, 1.47};
  
  constexpr std::array<Numeric, n> y0 = {-0.041, 0.277, -0.372, 0.559, -0.573, 0.618, -0.366, 0.278,
                                         -0.089, -0.021, 0.060, -0.152, 0.216, -0.293, 0.373, -0.436,
                                         0.491, -0.542, 0.571, -0.613, 0.636, -0.670, 0.690, -0.718,
                                         0.740, -0.763, 0.788, -0.807, 0.834, -0.849, 0.876, -0.887,
                                         0.915, -0.922, 0.950, -0.955, 0.987, -0.988, 0, 0, 0, 0, 0, 0};
  
  constexpr std::array<Numeric, n> y1 = {0, 0.124, -0.002, 0.008, 0.045, -0.093, 0.264, -0.351, 0.359,
                                         -0.416, 0.326, -0.353, 0.484, -0.503, 0.579, -0.590, 0.616,
                                         -0.619, 0.611, -0.609, 0.574, -0.568, 0.574, -0.566, 0.60,
                                         -0.59, 0.63, -0.62, 0.64, -0.63, 0.65, -0.64, 0.65, -0.64,
                                         0.65, -0.64, 0.64, -0.62, 0, 0, 0, 0, 0, 0};
  
  constexpr std::array<Numeric, n> g0 = {-0.000695, -0.090,  -0.103, -0.239, -0.172, -0.171, 0.028,
                                         0.150, 0.132, 0.170, 0.087, 0.069, 0.083, 0.067, 0.007,
                                         0.016, -0.021, -0.066, -0.095, -0.115, -0.118, -0.140,
                                         -0.173, -0.186, -0.217, -0.227, -0.234, -0.242, -0.266,
                                         -0.272, -0.301, -0.304, -0.334, -0.333, -0.361, -0.358,
                                         -0.348, -0.344, 0, 0, 0, 0, 0, 0};
  
  constexpr std::array<Numeric, n> g1 = {0, -0.045, 0.007, 0.033, 0.081, 0.162, 0.179, 0.225, 0.054,
                                         0.003, 0.0004, -0.047, -0.034, -0.071, -0.180, -0.210, -0.285,
                                         -0.323, -0.363, -0.380, -0.378, -0.387, -0.392, -0.394, -0.424,
                                         -0.422, -0.465, -0.46, -0.51, -0.50, -0.55, -0.54, -0.58, -0.56,
                                         -0.62, -0.59, -0.68, -0.65, 0, 0, 0, 0, 0, 0};
  
  constexpr std::array<Numeric, n> dv0 = {-0.00028, 0.00597, -0.0195, 0.032, -0.0475, 0.0541, -0.0232,
                                          0.0154, 0.0007, -0.0084, -0.0025, -0.0014, -0.0004, -0.0020,
                                          0.005, -0.0066, 0.0072, -0.008, 0.0064, -0.0070, 0.0056, -0.0060,
                                          0.0047, -0.0049, 0.0040, -0.0041, 0.0036, -0.0037, 0.0033, -0.0034,
                                          0.0032, -0.0032, 0.0030, -0.0030, 0.0028, -0.0029, 0.0029, -0.0029,
                                          0, 0, 0, 0, 0, 0};
  
  constexpr std::array<Numeric, n> dv1 = {-0.00039, 0.009, -0.012, 0.016, -0.027, 0.029, 0.006, -0.015,
                                          0.010, -0.014, -0.013, 0.013, 0.004, -0.005, 0.010, -0.010, 0.010,
                                          -0.011, 0.008, -0.009, 0.003, -0.003, 0.0009, -0.0009, 0.0017,
                                          -0.0016, 0.0024, -0.0023, 0.0024, -0.0024, 0.0024, -0.0020, 0.0017,
                                          -0.0016, 0.0013, -0.0012, 0.0005, -0.0004, 0, 0, 0, 0, 0, 0};
  
  /*
  constexpr std::array<Index ,n> Np = {1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
                                       11, 11, 13, 13, 15, 15, 17, 17,
                                       19, 19,  21, 21, 23, 23, 25, 25,
                                       27, 27, 29, 29, 31, 31, 33, 33,
                                       35, 35, 37, 37, 1, 1, 1, 3, 3, 3};
  
  constexpr std::array<Index ,n> Npp = {1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
                                        11, 11, 13, 13, 15, 15, 17, 17,
                                        19, 19,  21, 21, 23, 23, 25, 25,
                                        27, 27, 29, 29, 31, 31, 33, 33,
                                        35, 35, 37, 37, 3, 3, 3, 5, 5, 5};
  
  constexpr std::array<Index ,n> Jp = {1, 1, 3, 3, 5, 5, 7, 7, 9, 9,
                                       11, 11, 13, 13, 15, 15, 17, 17,
                                       19, 19,  21, 21, 23, 23, 25, 25,
                                       27, 27, 29, 29, 31, 31, 33, 33,
                                       35, 35, 37, 37, 1, 2, 2, 3, 4, 4};
  
  constexpr std::array<Index ,n> Jpp = {0, 2, 2, 4, 4, 6, 6, 8, 8, 10,
                                        10, 12, 12, 14, 14, 16, 16, 18,
                                        18, 20,  20, 22, 22, 24, 24, 26,
                                        26, 28, 28, 30, 30, 32, 32, 34,
                                        34, 36, 36, 38, 2, 2, 3, 4, 4, 5};
  */
  
  // Model constants
  constexpr Numeric x = 0.754;
  constexpr Numeric t0 = 300;
  constexpr Numeric conversion = 0.1820 * 1.00000e-3 / (10.0 * log10_euler);
  
  // Test constants
  constexpr Index nf = 501;
  constexpr Numeric fstart = 25e9;
  constexpr Numeric fend = 165e9;
  constexpr Numeric t = 296;
  constexpr Numeric p = 1e5;
  Vector f(nf);
  nlinspace(f, fstart, fend, nf);
  ComplexVector xsec(nf);
  
  // Model implementation
  const Numeric theta = t0 / t;
  const Numeric theta_m1 = theta - 1;
  const Numeric theta_3 = pow3(theta);
  const Numeric theta_x = std::pow(theta, x);
  const Numeric theta_2x = pow2(theta_x);
//   const Numeric GD_div_F0 = Linefunctions::DopplerConstant(t, SpeciesTag("O2-66").SpeciesMass());
  
  xsec = 0;
  for (Index i=0; i<n; i++) {
//     const Numeric invGD = 1 / (GD_div_F0 * f0[i]*1e9);
//     const Numeric fac = inv_sqrt_pi * invGD;
    const Numeric G0 = 1e9 * 1e-5 * p * gamma[i] * theta_x;
    const Numeric Y = 1e-5 * p * (y0[i] + y1[i] * theta_m1) * theta_x;
    const Numeric G = pow2(1e-5 * p) * (g0[i] + g1[i] * theta_m1) * theta_2x;
    const Numeric DV = 1e9 * pow2(1e-5 * p) * (dv0[i] + dv1[i] * theta_m1) * theta_2x;
    const Numeric ST = 1e-6 * p * theta_3 * intens[i]/(1e9*f0[i]) * std::exp(-a2[i] * theta_m1);
    
    std::cout<<f0[i]<<" "<<ST<<"\n";
    
    for (Index j=0; j<f.nelem(); j++) {
      xsec[j] += ST * pow2(f[j]) * (
//         (G0 * (1 + G) + Y * (f[j] - f0[i]*1e9 - DV)) /
//         (pow2(f[j] - f0[i]*1e9 - DV) + pow2(G0)) +
        fac * Complex(1 + G, Y) *
      Faddeeva::w(Complex(f0[i]*1e9 + DV - f[j], G0) * invGD) +
      (G0 * (1 + G) - Y * (f[j] + f0[i]*1e9 + DV)) /
      (pow2(f[j] + f0[i]*1e9 + DV) + pow2(G0)));
    }
  }
  
  std::cout<<'\n';
  Matrix pxsec(nf, 1, 0);
//   MPM87O2AbsModel(pxsec, 0, 1, 1, 1, "user", f, Vector(1, p), Vector(1, t), Vector(1, 0), Vector(1, 1), Verbosity());
  PWR93O2AbsModel(pxsec, 0, 1, 1, 1, "user", "PWR98", f, Vector(1, p), Vector(1, t), Vector(1, 0), Vector(1, 1), Verbosity());
  
  xsec *= conversion/200;
  std::cout<< "xr = np.array([";
  for (auto xi: xsec)
    std::cout<<xi.real() << ',' << ' ';
  std::cout << ']' << ')' << '\n';
  
  std::cout<< "xr2 = np.array([";
  for (Index i=0; i<nf; i++)
    std::cout<<pxsec(i, 0) << ',' << ' ';
  std::cout << ']' << ')' << '\n';
}

int main() {
  /*test_speed_of_pressurebroadening();
    test_transmissionmatrix();
    test_r_deriv_propagationmatrix();
    test_transmat_from_propmat();
    test_transmat_to_cumulativetransmat();
    test_sinc_likes_0limit();*/
//   test_zeeman();
test_mpm20();
  return 0;
}
