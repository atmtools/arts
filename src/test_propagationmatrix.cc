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


void test_new_lineshapes()
{
  define_species_data();
  define_species_map();
  
  const Index this_species=0;
  const Index water_species=-1;
  const ArrayOfIndex broadening_species={-1,-2,-1,-1,-1,-1};
  const Numeric T=250;
  const Numeric P=10;
  const Numeric H=50e-6;
  GriddedField1 gf1; gf1.data={4.016432e-01, 7.315888e-01, -3.313678e-05, 6.642877e-08};  // O2-66 partition function
  const ArrayOfGriddedField1 part_data(1, gf1);
  
  const ArrayOfNumeric empty_aon;
  const SpeciesTag st("O2-66");
  std::cout<< st<<"\n";
  LineRecord line(st.Species(),st.Isotopologue(),100e9,0.0,4.03935532732085e-19,296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
  line.SetVersionToLatest();
  
  std::cout<<line<<"\n\n";
  
  const Numeric QT=single_partition_function(T, SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF, part_data);
  const Numeric dQTdT=dsingle_partition_function_dT(QT, T, 0.1, SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF, part_data);
  const Numeric QT0=single_partition_function(line.Ti0(), SpeciesAuxData::AT_PARTITIONFUNCTION_COEFF, part_data);
  
  RetrievalQuantity rq;
  rq.PropType(JacPropMatType::Frequency);
  std::cout<<"PROPTYPE: " << (Index(rq.PropMatType()) == Index(JacPropMatType::Frequency)) << "\n";
  rq.Perturbation(0.1);
  const ArrayOfRetrievalQuantity jacs(1, rq);
  const ArrayOfIndex i_jacs = equivlent_propmattype_indexes(jacs);
  std::cout<<i_jacs<<"\n";
  
  const Linefunctions::SingleLevelLineData leveldata(line, Vector(1, 1.0), Vector(0), T, P, H,0.0,1.0, QT, dQTdT, QT0, broadening_species, this_species, water_species, 0, jacs, i_jacs);
  std::cout<<leveldata<<"\n";
  
  /*for(Numeric f=99e9; f<101e9; f+=100e4)*/ {
    Numeric f = 101e9;
    
    Complex F, N;
    ComplexVector dF(1,Complex(0,0)), dN(1);
    Linefunctions::set_lineshape_from_level_line_data(F, N, dF, dN, f, T, leveldata, line, jacs, i_jacs);
    std::cout<<"NEW: " << F << " " << dF << "\n";
    ComplexVector Fv(1), Nv(0);
    ComplexMatrix dFm(1, 1), dNm(0, 1);
    Vector fv(1, f);
    Range r(0, 1);
    Linefunctions::set_cross_section_for_single_line(Fv, dFm, Nv, dNm, r, jacs, i_jacs,
                                                     line, fv, Vector(1,1), Vector(0), P, T, 
                                                     Linefunctions::DopplerConstant(T, line.IsotopologueData().Mass()), P,
                                                     1.0, H, Linefunctions::dDopplerConstant_dT(T, line.IsotopologueData().Mass()),
                                                     0.0, QT, dQTdT, QT0, broadening_species, this_species, water_species, 0, Verbosity());
    std::cout<<F<<" "<<Fv[0]<<" "<<1.0-Fv[0]/F << "\n";
    std::cout<<dF[0]<<" "<<dFm(0,0)<<" "<<1.0-dFm(0,0)/dF[0] << "\n";
  }
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
  const ArrayOfSpeciesTag aspt = {SpeciesTag("CO2"), SpeciesTag("H2O"), SpeciesTag("O3")};
  const Vector vmrs = {0.2, 0.3, 0.2};
  test.GetParams (G0, D0, G2, D2, FVC, ETA, Y, G, DV,
                  296., 246., 2., 0.2, vmrs, aspt);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  s = "VP # 1 AIR T1 16000 0.7 T1 100 1.3";
  istringstream x2(s);
  x2 >> test;
  std::cout << s << "\n";
  std::cout << test << "\n";
  test.GetParams (G0, D0, G2, D2, FVC, ETA, Y, G, DV,
                  296., 246., 2., 0.2, vmrs, aspt);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
  
  s = "VP # 1 H2O2 T1 16000 0.7 T1 100 1.3";
  istringstream x3(s);
  x3 >> test;
  std::cout << s << "\n";
  std::cout << test << "\n";
  test.GetParams (G0, D0, G2, D2, FVC, ETA, Y, G, DV,
                  296., 246., 2., 0.2, vmrs, aspt);
  std::cout << G0 << " " << D0 << " " << G2 << " " << D2 << " " << FVC << " " << ETA << " " << Y << " " << G << " " << DV << "\n";
}

int main()
{
    std::cout<<"Testing Propmat Partials\n";
//     test_new_lineshapes();
//     test_zeeman();
    test_linefunctionsdata();
    return 0;
}
