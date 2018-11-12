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

#include "absorption.h"
#include "linescaling.h"
#include "linefunctions.h"
#include "global_data.h"
#include "arts.h"
#include "zeeman.h"
#include "linefunctiondata.h"
#include "transmissionmatrix.h"

// void define_species_data();
// void define_species_map();

void test_matrix_buildup()
{
  const Numeric k11 = 1;
  const Numeric k12 = -0.51;
  const Numeric k13 = -0.21;
  const Numeric k14 = 0.31;
  const Numeric k23 = -0.1;
  const Numeric k24 = -0.99;
  const Numeric k34 = 2;
  
  const Numeric r=0.5;
  
  const Numeric a = - k11*r;
  const Numeric b = - k12*r;
  const Numeric c = - k13*r;
  const Numeric d = - k14*r;
  const Numeric u = - k23*r;
  const Numeric v = - k24*r;
  const Numeric w = - k34*r;
  
  const Numeric b2 = b * b, c2 = c * c,
  d2 = d * d, u2 = u * u,
  v2 = v * v, w2 = w * w;
  
  const Numeric Const2 = b2 + c2 + d2 - u2 - v2 - w2;
  
  Numeric Const1;
  Const1  = b2 * (b2 * 0.5 + c2 + d2 - u2 - v2 + w2);
  Const1 += c2 * (c2 * 0.5 + d2 - u2 + v2 - w2);
  Const1 += d2 * (d2 * 0.5 + u2 - v2 - w2);
  Const1 += u2 * (u2 * 0.5 + v2 + w2);
  Const1 += v2 * (v2 * 0.5 + w2);
  Const1 *= 2;
  Const1 += 8 * (b * d * u * w - b * c * v * w - c * d * u * v);
  Const1 += w2 * w2;
  
  if(Const1 > 0.0)
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
  
  std::cout<<x<<" "<<y<<" "<<Const1<<" "<<Const2<<"\n";
  
  Numeric C0, C1, C2, C3;
  Numeric inv_y = 0.0, inv_x = 0.0;  // Init'd to remove warnings
  
  // X and Y cannot both be zero
  if(x == 0.0)
  {
    inv_y = 1.0 / y;
    C0 = 1.0;
    C1 = 1.0;
    C2 = (1.0 - cos_y) * inv_x2y2;
    C3 = (1.0 - sin_y*inv_y) * inv_x2y2;
  }
  else if(y == 0.0)
  {
    inv_x = 1.0 / x;
    C0 = 1.0;
    C1 = 1.0;
    C2 = (cosh_x - 1.0) * inv_x2y2;
    C3 = (sinh_x*inv_x - 1.0) * inv_x2y2;
  }
  else
  {
    inv_x = 1.0 / x;
    inv_y = 1.0 / y;
    
    C0 = (cos_y*x2 + cosh_x*y2) * inv_x2y2;
    C1 = (sin_y*x2*inv_y + sinh_x*y2*inv_x) * inv_x2y2;
    C2 = (cosh_x - cos_y) * inv_x2y2;
    C3 = (sinh_x*inv_x - sin_y*inv_y) * inv_x2y2;
  }
  
  std::cout<<C0<<" "<<C1<<" "<<C2<<" "<<C3<<"\n";
  
  Matrix F(4,4,0), A(4,4,0);
  
  MatrixViewMap eigF = MapToEigen(F);
  Eigen::Matrix4d eigA;
  eigA << 0,  b,  c, d, 
          b,  0,  u, v, 
          c ,-u,  0, w, 
          d, -v, -w, 0;
  
  
  eigF = C1 * eigA + C2 * eigA*eigA + C3 * eigA*eigA*eigA;
  eigF(0, 0) += C0; eigF(1, 1) += C0; eigF(2, 2) += C0; eigF(3, 3) += C0;
  eigF *= exp(a);
  
  std::cout<<F<<"\n";
}

void test_zeeman()
{
  const Vector rtp_mag = {20e-6, 20e-6, 0};
  
  for(Numeric za = 0; za <= 180; za += 15) {
//     std::cout << "arts_eta.append([";
    for(Numeric aa = 0; aa <= 360; aa += 15) {
      const Vector rtp_los = {za, aa};
      std::cout << "ZA.append(" << za << "),"
                << " AA.append(" << aa << ")," 
                << " THETA.append(" << RAD2DEG * zeeman_magnetic_theta(rtp_mag[0], rtp_mag[1], rtp_mag[2], DEG2RAD*za, DEG2RAD*aa) << "),"
                << " ETA.append(" << RAD2DEG * zeeman_magnetic_eta(rtp_mag[0], rtp_mag[1], rtp_mag[2], DEG2RAD*za, DEG2RAD*aa) << "),"
                << " deriv.append(" << RAD2DEG * zeeman_magnetic_deta_du(rtp_mag[0], rtp_mag[1], rtp_mag[2], DEG2RAD*za, DEG2RAD*aa) << ")\n";
    }
//     std::cout << "])\n";
  }
}




void test_linefunctionsdata()
{
  define_species_data();
  define_species_map();
  
  std::tuple<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric> X;
  String s = "VP LM1 4 SELF T1 16000 0.7 T1 100 1.3 T4 0.7e-4 0.5e-6 0.7 CO2 T1 16001 0.71 T1 101 1.31 T1 0.9e-4 0.7 H2O T1 16001.1 0.711 T1 101.1 1.311 T1 0.4e-4 0.7 AIR T1 18002 0.72 T1 102 1.32 T1 0.1e-4 0.7 THIS-IS-NOT-READ-BY-CIN";
  istringstream x1(s);
  LineFunctionData test;
  x1 >> test;
  
  std::cout<<s<<"\n";
  std::cout<<test<<"\n";
  
  String n;
  x1 >> n;
  std::cout << n << "\n";
  
  Numeric G0, D0, G2, D2, FVC, ETA, Y, G, DV;
  const ArrayOfArrayOfSpeciesTag aspt = {{SpeciesTag("CO2")}, {SpeciesTag("H2O")}, {SpeciesTag("H2O2")}};
  const Vector vmrs = {0.2, 0.3, 0.2};
  X = test.GetParams (296., 246., 2., 0.2, vmrs, aspt);
  G0=std::get<0>(X);
  D0=std::get<1>(X);
  G2=std::get<2>(X);
  D2=std::get<3>(X);
  FVC=std::get<4>(X);
  ETA=std::get<5>(X);
  Y=std::get<6>(X);
  G=std::get<7>(X);
  DV=std::get<8>(X);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  X = test.GetParams (296., 247., 2., 0.2, vmrs, aspt);
  G0=std::get<0>(X);
  D0=std::get<1>(X);
  G2=std::get<2>(X);
  D2=std::get<3>(X);
  FVC=std::get<4>(X);
  ETA=std::get<5>(X);
  Y=std::get<6>(X);
  G=std::get<7>(X);
  DV=std::get<8>(X);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  s = "VP # 1 AIR T1 16000 0.7 T1 100 1.3";
  istringstream x2(s);
  x2 >> test;
  std::cout << s << "\n";
  std::cout << test << "\n";
  X = test.GetParams (296., 246., 2., 0.2, vmrs, aspt);
  G0=std::get<0>(X);
  D0=std::get<1>(X);
  G2=std::get<2>(X);
  D2=std::get<3>(X);
  FVC=std::get<4>(X);
  ETA=std::get<5>(X);
  Y=std::get<6>(X);
  G=std::get<7>(X);
  DV=std::get<8>(X);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  X = test.GetParams (296., 247., 2., 0.2, vmrs, aspt);
  G0=std::get<0>(X);
  D0=std::get<1>(X);
  G2=std::get<2>(X);
  D2=std::get<3>(X);
  FVC=std::get<4>(X);
  ETA=std::get<5>(X);
  Y=std::get<6>(X);
  G=std::get<7>(X);
  DV=std::get<8>(X);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  s = "VP # 1 H2O2 T1 16000 0.7 T1 100 1.3";
  istringstream x3(s);
  x3 >> test;
  std::cout << s << "\n";
  std::cout << test << "\n";
  X = test.GetParams(296., 246., 2., 0.2, vmrs, aspt);
  G0=std::get<0>(X);
  D0=std::get<1>(X);
  G2=std::get<2>(X);
  D2=std::get<3>(X);
  FVC=std::get<4>(X);
  ETA=std::get<5>(X);
  Y=std::get<6>(X);
  G=std::get<7>(X);
  DV=std::get<8>(X);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  X = test.GetParams(296., 247., 2., 0.2, vmrs, aspt);
  G0=std::get<0>(X);
  D0=std::get<1>(X);
  G2=std::get<2>(X);
  D2=std::get<3>(X);
  FVC=std::get<4>(X);
  ETA=std::get<5>(X);
  Y=std::get<6>(X);
  G=std::get<7>(X);
  DV=std::get<8>(X);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
}


void test_speed_of_pressurebroadening()
{
  constexpr bool test_new=true;
  
  define_species_data();
  define_species_map();
  
  const Index N=5000000;
  const Numeric T0=100, T1=300, dT=(T1-T0)/N;
  const Numeric P = 133.33;
  const ArrayOfArrayOfSpeciesTag aspt = {{SpeciesTag("CO2")}, {SpeciesTag("H2O")}, {SpeciesTag("H2O2")}, {SpeciesTag("H2")}, {SpeciesTag("He")}};
  const Vector vmrs = {0.05, 0.3, 0.2, 0.1, 0.2};
  
  Numeric T = T0;
  Vector G0(N), D0(N), G2(N), D2(N), FVC(N), ETA(N), Y(N), G(N), DV(N);
  std::tuple<Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric, Numeric> X;
  
  // New line shape
  LineFunctionData lf;
  String new_lf = "VP LM1 2" 
  "SELF T1 16000 0.70 T5 102 0.72 T1 0.1e-4 0.7" 
  "AIR  T1 18002 0.72 T5 102 0.72 T1 0.1e-4 0.7";
  istringstream stream_new_lf(new_lf);
  stream_new_lf >> lf;
  
  // Old line shape
  PressureBroadeningData pb;
  pb.SetAirBroadeningFromCatalog(16000, 0.70, 18002, 0.72, 102, 0, 0, 0, 0, 0);
  LineMixingData lm;
  lm.StorageTag2SetType("L1");
  lm.SetDataFromVectorWithKnownType(Vector({296, 0.1e-4 ,0.7}));
  
  std::cout << lf << "\n";
  std::cout << LineFunctionData(pb, lm, "CO2", 296) << "\n";
  
  Index i=0;
  if(test_new) {
    for(i=0; i<N; i++) {
      X = lf.GetParams (296., T, P, vmrs[0], vmrs, aspt);
      G0[i]=std::get<0>(X);
      D0[i]=std::get<1>(X);
      G2[i]=std::get<2>(X);
      D2[i]=std::get<3>(X);
      FVC[i]=std::get<4>(X);
      ETA[i]=std::get<5>(X);
      Y[i]=std::get<6>(X);
      G[i]=std::get<7>(X);
      DV[i]=std::get<8>(X);
      T += dT;
    }
  }
  else {
    for(i=0; i<N; i++) {
      pb.GetPressureBroadeningParams(G0[i], G2[i], ETA[i], D0[i], D2[i], FVC[i],
                                     T, 296, P, P*vmrs[0],  0, 1, {}, vmrs);
      lm.GetLineMixingParams(Y[i], G[i], DV[i], T, P, 0.0, 1);
      T += dT;
    }
  }
  
  if(test_new)
    printf("V0 = np.array([\n");
  else
    printf("V1 = np.array([\n");
  
  for(i=N-2; i<N; i++)
    printf("[%f, %f, %f, %f, %f, %f, %f, %f, %f], \n", 
             G0[i], D0[i], G2[i], D2[i], FVC[i], ETA[i], Y[i], G[i], DV[i]);
  printf("]).T\n");
}


void test_transmissionmatrix()
{
  // Initializes as unity matrices
  TransmissionMatrix a(2, 4);
  std::cout << "Initialized TransmissionMatrix(2, 4):\n" << a << "\n";
  
  // Set a single input
  Eigen::Matrix4d A;
  A << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16;
  std::cout << "New Matrix:\n" << A << "\n\n";
  a.set4(A, 1);
  std::cout << "Updated TransmissionMatrix Position 1 wit New Matrix:\n" << a << "\n";
  
  // The stream can also set the values
  String S = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 125 26 27 28 29 30 31 32";
  std::cout << "Stream:\n" << S << "\n\n";
  std::istringstream astream(S);
  astream >> a;
  std::cout << "Streamed into TransmissionMatrix:\n" << a << "\n";
  
  // Initialize empty
  RadiationVector b(3, 3);
  std::cout << "Initialized RadiationVector(3, 3)\n" << b << "\n";
  
  // Set is not defined but add is
  Eigen::Vector3d B;
  B << 1,2,3;
  std::cout << "New Vector:\n" << B << "\n\n";  // nb. not transposed
  b.add3(B, 1);
  std::cout << "Updated RadiationVector Position 1 with New Vector:\n" << b << "\n";
  
  // The stream can also set the values
  String T = "1 2 3 4 5 6 7 8 90";
  std::cout << "Stream:\n" << T << "\n\n";
  std::istringstream bstream(T);
  bstream >> b;
  std::cout << "Streamed into RadiationVector:\n" << b << "\n";
}


int main()
{
    std::cout<<"Testing Propmat Partials\n";
//     test_new_lineshapes();
//     test_zeeman();
//     test_linefunctionsdata();
//     test_speed_of_pressurebroadening();
    test_transmissionmatrix();
    return 0;
}
