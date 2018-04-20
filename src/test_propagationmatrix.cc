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
#include "jacobian.h"
#include "linescaling.h"
#include "lineshapes.cc"
#include "lin_alg.h"
#include "rte.h"
#include "propagationmatrix.h"
#include "linefunctions.h"
#include "wigner_functions.h"

void test_pressurebroadening()
{
    const ArrayOfNumeric empty_aon;
    
    const LineRecord line(species_index_from_species_name("O2"),0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    const Numeric vmr=0.2, dT = 1e-1;
    
    Numeric tmp1,tmp2,tmp3,tmp4, g_0, df_0, dg, dg2, de, ddf, ddf2, dfvc, g_1, df_1;
    ArrayOfIndex empty1;
    const Vector vmrs(1,vmr);
    Verbosity none;
    
    for(Numeric p =-2;p<6;p++)
    {
        const Numeric P = pow(10.,p);
        for(Numeric T = 150; T<350; T+=20)
        {
            line.PressureBroadening().GetPressureBroadeningParams(g_0,tmp1,tmp2,df_0,tmp3,tmp4,
                                                                  line.Ti0()/T, P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams(g_1,tmp1,tmp2,df_1,tmp3,tmp4,
                                                                  line.Ti0()/(T+dT), P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams_dT(dg, dg2, de, ddf, ddf2, dfvc, T, 
                                                                     line.Ti0(),P,
                                                                     vmr*P,0,-1,
                                                                     empty1,
                                                                     vmrs,none);
            //std::cout<<((g_1-g_0)/dT-dg)/dg<< "\n";
            /* This gives relative errors of (minus) 0.001 at 150 K and 0.0001 at 300 K.  Constant with pressure. */
        }
    }
}

Numeric test2(Numeric x)
{
    return x*x;
}
void test_partitionfunction()
{
    
    const ArrayOfNumeric empty_aon;
    const LineRecord line(species_index_from_species_name("O2"),0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    Numeric qt_cache, qt_cache_d, qref_cache, qref_cache_d, dqt;
    
    const Numeric dT =1e-1;

    // Numeric x = 7.315888e-01;
    Vector part_data;
    part_data={ 4.016432e-01,  7.315888e-01,  -3.313678e-05,  6.642877e-08 };//O2 66 partition function

    
    
    for(Numeric T = 150; T<350; T+=20)
    {
        CalculatePartitionFctFromCoeff(qref_cache,qt_cache,line.Ti0(),T,part_data);
        CalculatePartitionFctFromCoeff(qref_cache_d,qt_cache_d,line.Ti0(),T+dT,part_data);
        
        CalculatePartitionFctFromCoeff_dT(dqt, T, part_data);
        
        //std::cout<<(((qref_cache_d/qt_cache_d - qref_cache/qt_cache) /dT) - 
        //-qref_cache/qt_cache/qt_cache * dqt) /(-qref_cache/qt_cache/qt_cache * dqt) <<"\n";
        /* This gives relative errors of less than (minus) 0.001 */
    }
}

void test_K1_and_K2()
{
    // Physical constants
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const ArrayOfNumeric empty_aon;
    const LineRecord line(species_index_from_species_name("O2"),0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    const Numeric dT = 0.1;
    
    for(Numeric T = 150; T<350; T+=20)
    {
        // Following Futbolin's division into two parts for the Boltzmann ratio because
        // gamma is also used for the NLTE part later on
        const Numeric gamma = exp( - PLANCK_CONST * line.F() / ( BOLTZMAN_CONST * T ) );
        const Numeric gamma_ref = exp( - PLANCK_CONST * line.F() / ( BOLTZMAN_CONST * line.Ti0()) );
        // Stimulated emission
        const Numeric K2 = (1.-gamma)/(1.-gamma_ref);
        // Boltzmann level
        const Numeric K1 = exp( line.Elow() / BOLTZMAN_CONST * (T-line.Ti0())/(T*line.Ti0()) );
        
        // Following Futbolin's division into two parts for the Boltzmann ratio because
        // gamma is also used for the NLTE part later on
        const Numeric gamma_d = exp( - PLANCK_CONST * line.F() / ( BOLTZMAN_CONST * (T+dT) ) );
        // Stimulated emission
        const Numeric K2_d = (1.-gamma_d)/(1.-gamma_ref);
        // Boltzmann level
        const Numeric K1_d = exp( line.Elow() / BOLTZMAN_CONST * ((T+dT)-line.Ti0())/((T+dT)*line.Ti0()) );
        
        const Numeric dK1 =  line.Elow()/BOLTZMAN_CONST/T/T * K1;
        const Numeric dK2 = -line.F()*PLANCK_CONST/BOLTZMAN_CONST/T/T * (gamma/(1.0-gamma_ref));
        
        //std::cout<< (((K2_d-K2)/dT)- dK2)/dK2 <<"\n";
        /* This gives relative errors of just above (minus) 0.001, 
           and is fairly constant with Temperature */
        
        //std::cout<<((K1_d-K1)/dT-dK1)/dK1<<"\n";
        /* This gives relative errors of just above (minus) 0.0001, 
           and is fairly constant with Temperature */
        
        if(!(abs(((K1_d-K1)/dT-dK1)/dK1)<0.001))
            std::cout<<"K1 partial checks out bad.\n";
        if(!(abs((((K2_d-K2)/dT)- dK2)/dK2)<0.01))
            std::cout<<"K1 partial checks out bad.\n";
    }
}

void test_lineshape()
{
    // Physical constants
    extern const Numeric SPEED_OF_LIGHT;
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric AVOGADROS_NUMB;
    
    static const Numeric doppler_const = sqrt(2.0 * BOLTZMAN_CONST *
    AVOGADROS_NUMB) / SPEED_OF_LIGHT; 
    const ArrayOfNumeric empty_aon;
    
    const LineRecord line(species_index_from_species_name("O2"),0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    const Numeric vmr=0.2, dT = 1e-1;
    
    Numeric tmp1,tmp2,tmp3,tmp4, g, df, g_d, df_d, dg, dg2, de, ddf, ddf2, dfvc;
    ArrayOfIndex empty1;
    const Vector vmrs(1,vmr);
    Verbosity none;
    Vector empty_vector;
    
    Vector f_grid(5);
    f_grid[0]=line.F()-1.0e6;
    f_grid[1]=line.F()-50.0e3;
    f_grid[2]=line.F();
    f_grid[3]=line.F()+50.0e3;
    f_grid[4]=line.F()+1.0e6;
    
    Vector ls(5,0.0),ls_d(5,0.0),ls_phase(5,0.0),dls_dx(5,0.0),dls_dy(5,0.0),dx_dT(5,0.0);
    
    for(Numeric p =-2;p<6;p++)
    {
        const Numeric P = pow(10., p);
        //std::cout<<"\nPressure: "<<P<<" Pa. dT: "<<dT<<" K\n";
        for(Numeric T = 150; T<350; T+=20)
        {
            line.PressureBroadening().GetPressureBroadeningParams(g,tmp1,tmp2,df,tmp3,tmp4,
                                                                  line.Ti0()/T, P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams(g_d,tmp1,tmp2,df_d,tmp3,tmp4,
                                                                  line.Ti0()/(T+dT), P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams_dT(dg, dg2, de, ddf, ddf2, dfvc, T, 
                                                                     line.Ti0(),P,
                                                                     vmr*P,0,-1,
                                                                     empty1,
                                                                     vmrs,none);
            
            const Numeric sigma = line.F() * doppler_const 
            *sqrt( T / 31.989830);//mass is O2-66
            
            const Numeric dsigma_dT = line.F() * doppler_const 
            *sqrt( T / 31.989830)/T/2.0;//mass is O2-66
            
            const Numeric sigma_d = line.F() * doppler_const 
            *sqrt( (T+dT) / 31.989830);//mass is O2-66
            
            //std::cout<<((sigma_d-sigma)/dT - dsigma_dT)/dsigma_dT<<"\n";
            /* This gives relative errors of just above (minus) 0.0001, 
               for 150 K and much less at higher temperatures*/
            
            faddeeva_algorithm_916(    ls,
                                       ls_phase,
                                       dls_dx,
                                       empty_vector,
                                       dls_dy,
                                       empty_vector,
                                       empty_vector,
                                       line.F(),
                                       g,
                                       0.0,0.0,0.0,0.0,
                                       sigma,
                                       0.0,
                                       f_grid,
                                       false,
                                       true);
            
            faddeeva_algorithm_916(    ls_d,
                                       empty_vector,
                                       empty_vector,
                                       empty_vector,
                                       empty_vector,
                                       empty_vector,
                                       empty_vector,
                                       line.F(),
                                       g_d,
                                       0.0,0.0,0.0,0.0,
                                       sigma_d,
                                       0.0,
                                       f_grid,
                                       false,
                                       false);
            
            Numeric dy_dT,dnorm_dT;
            w_x_plus_iy_dT(dx_dT,
                           dy_dT,
                           dnorm_dT,
                           f_grid,     
                           line.F(),    // Note that this is NOT the line center, but the shifted center (normally, and in our case the shift is nil)
                           sigma, 
                           df,    
                           0.0,    
                           dsigma_dT, 
                           g,
                           dg);
            
//             std::cout<<"Temperature: "<<T<<" K\n";
//             for(Index iv=0;iv<5; iv++)
//             {
//                 std::cout<<((ls_d[iv]-ls[iv])/dT - dx_dT[iv]*dls_dx[iv] 
//                            - dy_dT*dls_dy[iv] - dnorm_dT*ls[iv])/
//                            (dx_dT[iv]*dls_dx[iv] + dy_dT*dls_dy[iv] + dnorm_dT*ls[iv])<<" ";
//             }
//            std::cout<<"\n";
            /* This gives relative errors of just above (minus) 4% at 0.01 Pa and 230 K
               near the Doppler half-width but much less at line center and far from 
               the line, where errors are generally less than 1%.  For pressures above 
               1 Pa, the errors are on magnitude of less than (minus) 0.001.*/
        }
    }
}


Complex test_faddeeva()
{
  Complex F;
  
  Vector frequencies; 
  linspace(frequencies, 0, 50, 1);
  
  Vector gamma;
  linspace(gamma, 0, 1000000, 1);
  
  // On tested architecture with normal mode
  // About 76 nanoseconds per call to w(z)
  // About 32 nanoseconds per call to lorentz(z)
  for(Numeric& g : gamma)
  {
    for(Numeric& f : frequencies)
    {
//       const Complex zv = Complex(f, g);
//       Faddeeva::w(zv); 
      
      const Complex zl = Complex(g, f);
      0.3183098861837907 / zl;
    }
  }
  
  return F;
}


void test_the_class()
{
  PropagationMatrix a(50, 4);
  Matrix tmp(4, 4, 0.);
  
  a.Kjj() += Vector(50, 1.2);
  a.K12() += Vector(50, 1.7);
  a.K13() += Vector(50, 1.1);
  a.K23() += Vector(50, 1.3);
  
  a.MatrixAtPosition(tmp, 1);
  std::cout << tmp << "\n";
  
  Tensor3 b(50, 4, 4);
  compute_transmission_matrix(b, 1.0, a, a);
  
  std::cout << b(0, joker, joker) << "\n";
}


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


// void test_new_lineshapes()
// {
//   const ArrayOfNumeric empty_aon;
//   
//   const LineRecord line(species_index_from_species_name("O2"),0,61150556350.7454,0.0,4.03935532732085e-19,
//                         296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
//                         0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
//   Vector part_data;
//   part_data = { 4.016432e-01,  7.315888e-01,  -3.313678e-05,  6.642877e-08 };//O2 66 partition function
//   
//   const QuantumIdentifier QI = line.QuantumIdentity();
//   
//   Vector f_grid(5);
//   f_grid[0]=line.F()-50.0e9;
//   f_grid[1]=line.F()-50.0e8;
//   f_grid[2]=line.F();
//   f_grid[3]=line.F()+50.0e8;
//   f_grid[4]=line.F()+500.0e9;
//   
//   ComplexVector F(5), F2(5);
//   ArrayOfComplexVector dF(1);
//   dF[0].resize(5);
//   
//   PropmatPartialsData ppd;
//   ppd.SetOnlyTemperatureTrue();
//   
//   const Numeric pressure = 1e5;
//   const Numeric vmr = 0.2;
//   Numeric G0, G2, eta, L0, L2, FVC, dG0, dG2, de, dL0, dL2, dFVC, T, gd_div_f0, dgd_div_f0_dT, qt, qt0, dqt, K1, K2, dK1, dK2;
//   const Numeric dT = 1/10.;  
//   
//   T = 300;
//   gd_div_f0 = Linefunctions::DopplerConstant(T, 32.0);
//   dgd_div_f0_dT = Linefunctions::dDopplerConstant_dT(T, 32.0);
//   line.PressureBroadening().GetPressureBroadeningParams(G0, G2, eta, L0, L2, FVC, line.Ti0()/T, pressure,
//                                                    vmr*pressure, -1, -1, ArrayOfIndex(), 
//                                                    Vector(), Verbosity());
//   line.PressureBroadening().GetPressureBroadeningParams_dT(dG0, dG2, de, dL0, dL2, dFVC, T, line.Ti0(), pressure,
//                                                       vmr*pressure, -1, -1, ArrayOfIndex(), 
//                                                       Vector(), Verbosity());
//   CalculatePartitionFctFromCoeff(qt0, qt, line.Ti0(), T, part_data);
//   CalculatePartitionFctFromCoeff_dT(dqt, T, part_data);
//   K1 = boltzman_ratio(T, line.Ti0(), line.Elow());
//   K2 = stimulated_relative_emission(stimulated_emission(T, line.F()), 
//                                     stimulated_emission(line.Ti0(), line.F()));
//   dK1 = dboltzman_ratio_dT(K1, T, line.Elow());
//   dK2 = dstimulated_relative_emission_dT(stimulated_emission(T, line.F()), 
//                                          stimulated_emission(line.Ti0(), line.F()),
//                                          line.F(), T);
//   
//   Linefunctions::set_faddeeva_algorithm916(F, dF, f_grid, 0, 0, line.F(), gd_div_f0, G0, L0, 0, ppd, QI, dgd_div_f0_dT, dG0, dL0, 0);
//   
//   Linefunctions::set_lorentz(F, dF,f_grid, 0., 0., line.F(), G0, L0, 0.0, ppd, QI, dG0, dL0, 0.0);
//   
//   std::cout<<dF[0].real()<<"\n";
//   std::cout<< K2 << " " << dK1 << " " << dK2<<"\n";
//   
//   T += dT;
//   gd_div_f0 = Linefunctions::DopplerConstant(T, 32.0);
//   line.PressureBroadening().GetPressureBroadeningParams(G0, G2, eta, L0, L2, FVC, line.Ti0()/T, pressure,
//                                                         vmr*pressure, -1, -1, ArrayOfIndex(), 
//                                                         Vector(), Verbosity());
//   
//   Linefunctions::set_faddeeva_algorithm916(F2, dF, f_grid, 0, 0, line.F(), gd_div_f0, G0, L0, 0);
//   
//   Linefunctions::set_lorentz(F2, dF,f_grid, 0., 0., line.F(), G0, L0, 0.0);
//   
//   for(Index i = 0; i < F.nelem(); i++)
//     std::cout<<1/dT*(F2[i].real()-F[i].real())<<" ";
//   std::cout<<"\n";
//   for(Index i = 0; i < F.nelem(); i++)
//     std::cout<<1/dT*(F2[i].real()-F[i].real())/dF[0][i].real()<<" ";
//   std::cout<<"\n";
//   
// }


void test_erfc()
{
  const Numeric xstart = -300;
  const Numeric xend   =  300;
  const Numeric dx     =  1;
  Numeric x            =  xstart;
  
  while(x < xend)
  {
    std::cout<<x<<" "<<Faddeeva::w(Complex(x, 0))<<"\n";
    x += dx;
  }
}


void test_funs_compression()
{
  Index n = 1 << 5;
  Vector v;
  nlinspace(v, 0, Numeric(n), n+1);
  Index i = 0;
  Range r(joker);
  
  Vector f1(n+1, 0), f2(n+1, 0);
  
  std::cout << std::endl << "Partial" << std::endl;
  while(1 << i <= n)
  {
    r = Linefunctions::binary_range(v, i, false);
    std::cout << v[r] << std::endl;
    i++;
    
    for(Index j = 0; j < r.get_extent(); j++)
      f1[r(j)] = 5 * v[r(j)];
  }
  
  i = 0;
  std::cout << std::endl << "Full" << std::endl;
  while(1 << i <= n)
  {
    r = Linefunctions::binary_range(v, i, true);
    std::cout << v[r] << std::endl;
    i++;
    
    for(Index j = 0; j < r.get_extent(); j++)
      f2[r(j)] = 5 * v[r(j)];
  }
  
  std::cout << std::endl << f1 << std::endl << f2 << std::endl;
}

/*
void test_transmission(int i)
{
  Vector v4(4), o4(4);
  v4 = 1;
  Eigen::Vector4d ev4, eo4;
  ev4 << 1,1,1,1;
  std::cout<<v4<<"\n"<<ev4<<"\n";
  Tensor3 t3(i, 4, 4);
  for(int j=0; j<i; j++) for(int m=0; m<4; m++) for(int n=0; n<4; n++) t3(j, m, n) = (1.0 + j*m + j*n + j + m + n)/i/i;
  
  TransmissionMatrix et3(t3);
  ArrayOfMatrix4f aom4f = et3.GetData4();
  Index count = 0;
  
  for(int j=0; j<i; j++){
    count++;
    if(i%2)
      mult(v4, t3(j, joker, joker), o4);
    else
      mult(o4, t3(j, joker, joker), v4);
  }
  
  for(int j=0; j<i; j++){
    count++;
    if(i%2)
      ev4.noalias() = aom4f[j] * eo4;
    else
      eo4.noalias() = aom4f[j] * ev4;
  }
  
  std::cout<<t3(i-1,joker,joker)<<"\n";
  std::cout<<aom4f[i-1]<<"\n";
  std::cout<<o4<<"\n"<<eo4<<"\n";
  std::cout<<count<<"\n";
}*/


int main()
{
    std::cout<<"Testing Propmat Partials\n";
//     test_pressurebroadening();
//     test_partitionfunction();
//     test_K1_and_K2();
//     test_lineshape();
//     test_the_class();
//    test_new_lineshapes();
//    test_erfc();
     test_funs_compression();
    return 0;
}
