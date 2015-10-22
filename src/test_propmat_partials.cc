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
 * \file   test_promat_partials.cc
 * \author <ric.larsson@gmail.com>
 * \date   2015-10-20
 * 
 * \brief  Test Propmat Internal Partial Derivatives
 */

#include "absorption.h"
#include "jacobian.h"
#include "linescaling.h"
#include "lineshapes.cc"
#include "lin_alg.h"
#include "rte.h"

void test_pressurebroadening()
{
    const ArrayOfNumeric empty_aon;
    
    const LineRecord line(6,0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    const Numeric vmr=0.2, dT = 1e-1;
    
    Numeric tmp1,tmp2,tmp3,tmp4, g_0, df_0, dg, ddf, g_1, df_1;
    ArrayOfIndex empty1;
    const Vector vmrs(1,vmr);
    Verbosity none;
    
    for(Numeric p =-2;p<6;p++)
    {
        const Numeric P = pow10(p);
        for(Numeric T = 150; T<350; T+=20)
        {
            line.PressureBroadening().GetPressureBroadeningParams(g_0,tmp1,tmp2,df_0,tmp3,tmp4,
                                                                  line.Ti0()/T, P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams(g_1,tmp1,tmp2,df_1,tmp3,tmp4,
                                                                  line.Ti0()/(T+dT), P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams_dT(dg,ddf, T, 
                                                                     line.Ti0(),P,
                                                                     vmr*P,0,-1,
                                                                     empty1,
                                                                     vmrs,none);
            //std::cout<<((g_1-g_0)/dT-dg)/dg<< "\n";
            /* This gives relative errors of (minus) 0.001 at 150 K and 0.0001 at 300 K.  Constant with pressure. */
        }
    }
}

void test_partitionfunction()
{
    
    const ArrayOfNumeric empty_aon;
    const LineRecord line(species_index_from_species_name("O2"),0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    Numeric qt_cache, qt_cache_d, qref_cache, qref_cache_d, dqt;
    
    const Numeric dT =1e-1;
    
    Vector part_data;
    part_data=MakeArray<Numeric>( 4.016432e-01,  7.315888e-01,  -3.313678e-05,  6.642877e-08 );//O2 66 partition function
    
    
    
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
    
    const LineRecord line(6,0,61150556350.7454,0.0,4.03935532732085e-19,
                          296.0,2.5505950629926e-21,13255.072408981,13047.9619025907,0.8,0.8,
                          0.0,empty_aon,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    
    const Numeric vmr=0.2, dT = 1e-1;
    
    Numeric tmp1,tmp2,tmp3,tmp4, g, df, g_d, df_d, dg, ddf;
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
        const Numeric P = pow10(p);
        //std::cout<<"\nPressure: "<<P<<" Pa. dT: "<<dT<<" K\n";
        for(Numeric T = 150; T<350; T+=20)
        {
            line.PressureBroadening().GetPressureBroadeningParams(g,tmp1,tmp2,df,tmp3,tmp4,
                                                                  line.Ti0()/T, P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams(g_d,tmp1,tmp2,df_d,tmp3,tmp4,
                                                                  line.Ti0()/(T+dT), P, vmr*P,0,-1,empty1,vmrs,none);
            
            line.PressureBroadening().GetPressureBroadeningParams_dT(dg,ddf, T, 
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

void test_matrixexp()
{
   // Startup
   Matrix tmp(4,4);
   Numeric length;
   
   // Below are examples from one run.
   
   // length for ip
   length=8407.36;
   //dtdx for ip
   tmp(0,0)=1.06286e-08; tmp(0,1)=3.81257e-10; tmp(0,2)=7.23907e-10; tmp(0,3)=-5.71387e-09; 
   tmp(1,0)=3.81268e-10; tmp(1,1)=1.06286e-08; tmp(1,2)=1.0299e-08; tmp(1,3)=4.25483e-09; 
   tmp(2,0)=7.23901e-10; tmp(2,1)=-1.0299e-08; tmp(2,2)=1.06286e-08; tmp(2,3)=-2.24089e-09; 
   tmp(3,0)=-5.71387e-09; tmp(3,1)=-4.25481e-09; tmp(3,2)=2.24093e-09; tmp(3,3)=1.06286e-08; 
   const Matrix dtdx=tmp;
   //ppath_ext for ip+1
   tmp(0,0)=6.57435e-10; tmp(0,1)=2.38216e-11; tmp(0,2)=4.52381e-11; tmp(0,3)=-3.49977e-10; 
   tmp(1,0)=2.38216e-11; tmp(1,1)=6.57435e-10; tmp(1,2)=6.81735e-10; tmp(1,3)=2.71398e-10; 
   tmp(2,0)=4.52381e-11; tmp(2,1)=-6.81735e-10; tmp(2,2)=6.57435e-10; tmp(2,3)=-1.42913e-10; 
   tmp(3,0)=-3.49977e-10; tmp(3,1)=-2.71398e-10; tmp(3,2)=1.42913e-10; tmp(3,3)=6.57435e-10; 
   const Matrix ppath_ext_1=tmp;
   //ppath_ext for ip
   tmp(0,0)=2.08609e-10; tmp(0,1)=7.33668e-12; tmp(0,2)=1.39302e-11; tmp(0,3)=-1.11706e-10; 
   tmp(1,0)=7.33668e-12; tmp(1,1)=2.08609e-10; tmp(1,2)=2.13876e-10; tmp(1,3)=8.65063e-11; 
   tmp(2,0)=1.39302e-11; tmp(2,1)=-2.13876e-10; tmp(2,2)=2.08609e-10; tmp(2,3)=-4.55607e-11; 
   tmp(3,0)=-1.11706e-10; tmp(3,1)=-8.65063e-11; tmp(3,2)=4.55607e-11; tmp(3,3)=2.08609e-10; 
   const Matrix ppath_ext_0=tmp
   ;//dppath_ext_dx for ip
   tmp(0,0)=2.08609e-10; tmp(0,1)=7.33668e-12; tmp(0,2)=1.39302e-11; tmp(0,3)=-1.11706e-10; 
   tmp(1,0)=7.33668e-12; tmp(1,1)=2.08609e-10; tmp(1,2)=0; tmp(1,3)=0; 
   tmp(2,0)=1.39302e-11; tmp(2,1)=0; tmp(2,2)=2.08609e-10; tmp(2,3)=0; 
   tmp(3,0)=-1.11706e-10; tmp(3,1)=0; tmp(3,2)=0; tmp(3,3)=2.08609e-10; 
   const Matrix ppath_ext_0_test=tmp
   ;//dppath_ext_dx for ip
   tmp(0,0)=-2.53045e-12; tmp(0,1)=-9.07729e-14; tmp(0,2)=-1.72351e-13; tmp(0,3)=1.36037e-12; 
   tmp(1,0)=-9.07729e-14; tmp(1,1)=-2.53045e-12; tmp(1,2)=-2.59077e-12; tmp(1,3)=-1.2184e-12; 
   tmp(2,0)=-1.72351e-13; tmp(2,1)=2.59077e-12; tmp(2,2)=-2.53045e-12; tmp(2,3)=6.41701e-13; 
   tmp(3,0)=1.36037e-12; tmp(3,1)=1.2184e-12; tmp(3,2)=-6.41701e-13; tmp(3,3)=-2.53045e-12; 
   const Matrix dppath_ext_dx=tmp;
   tmp(0,0)=-2.53045e-12; tmp(0,1)=-9.07729e-14; tmp(0,2)=-1.72351e-13; tmp(0,3)=1.36037e-12; 
   tmp(1,0)=-9.07729e-14; tmp(1,1)=-2.53045e-12; tmp(1,2)=0; tmp(1,3)=0; 
   tmp(2,0)=-1.72351e-13; tmp(2,1)=0; tmp(2,2)=-2.53045e-12; tmp(2,3)=0; 
   tmp(3,0)=1.36037e-12; tmp(3,1)=0; tmp(3,2)=0; tmp(3,3)=-2.53045e-12; 
   const Matrix dppath_ext_dx_test=tmp;
   //ppath_ext_dt for ip
   tmp(0,0)=2.08356e-10; tmp(0,1)=7.32761e-12; tmp(0,2)=1.3913e-11; tmp(0,3)=-1.1157e-10; 
   tmp(1,0)=7.32761e-12; tmp(1,1)=2.08356e-10; tmp(1,2)=2.13631e-10; tmp(1,3)=8.64051e-11; 
   tmp(2,0)=1.3913e-11; tmp(2,1)=-2.13631e-10; tmp(2,2)=2.08356e-10; tmp(2,3)=-4.55074e-11; 
   tmp(3,0)=-1.1157e-10; tmp(3,1)=-8.64051e-11; tmp(3,2)=4.55074e-11; tmp(3,3)=2.08356e-10; 
   const Matrix ppath_ext_dt=tmp;
   //trans_partial for ip
   tmp(0,0)=0.999996; tmp(0,1)=-1.30978e-07; tmp(0,2)=-2.48724e-07; tmp(0,3)=1.94076e-06; 
   tmp(1,0)=-1.3098e-07; tmp(1,1)=0.999996; tmp(1,2)=-3.76485e-06; tmp(1,3)=-1.50451e-06; 
   tmp(2,0)=-2.48723e-07; tmp(2,1)=3.76485e-06; tmp(2,2)=0.999996; tmp(2,3)=7.92277e-07; 
   tmp(3,0)=1.94076e-06; tmp(3,1)=1.50451e-06; tmp(3,2)=-7.92283e-07; tmp(3,3)=0.999996; 
   const Matrix trans_partial=tmp;
   
   
   Matrix T(4,4),dT(4,4), A(4,4),dA(4,4);
   
   A=ppath_ext_0;
   A*=-length;
   dA=dppath_ext_dx;
   dA*=-0.5*length;
   matrix_exp_dmatrix_exp(T,dT,A,dA,6);
   
   tmp=dtdx;
   tmp-=dT;
   tmp/=dT;
   
//    std::cout<<"test:\n";
//    std::cout<<tmp<<"\n\ndtdx:\n";
//    Really large relative errors of almost 20%.  Still believable results since relative change is less than 1e-4 for normal cases.
}

int main()
{
    std::cout<<"Testing Propmat Partials\n";
    test_pressurebroadening();
    test_partitionfunction();
    test_K1_and_K2();
    test_lineshape();
    test_matrixexp();
    return 0;
}