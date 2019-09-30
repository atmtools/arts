/************************************************************************
 * $Id: disotest.c 2966 2013-07-24 08:58:48Z svn-kylling $
 ************************************************************************/

/*
 *   Copyright (c) 2011 by Timothy E. Dowling
 *   
 *   This file is part of cdisort.
 *
 *   cdisort is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   cdisort is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with cdisort.  If not, see <http://www.gnu.org/licenses/>.
 */
 
/*
 * disotest()
 *
 * Runs test problems for disort() and checks answers.
 * These problems test almost all logical branches in disort().
 *
 * This file contains the following functions:
 *   main()
 *   disort_test01,disort_test02,...,disort_test13
 *   print_test()
 *
 * It is recommended that you use the code below as a template for creating your own
 * calls to disort(), rather than starting from scratch.  This will prevent mistakes
 * and ensure that every input argument gets a value.  Note in particular how getmom()
 * is sometimes called to fill an array section of PMOM(k,lc) (for one layer); several
 * people have done this incorrectly in attempting to write it ab initio. 
 *
 * Routines called :
 *   c_disort():   The discrete ordinates radiative transfer program
 *   c_getmom():   Sets phase function Legendre coefficients
 *   print_test(): Prints fluxes and intensities and their ratios to the correct values
 *   c_ratio():    Ratio (calculated/correct) with under/overflow protection, or (1.+calculated) if correct == 0.
 *
 *-----------------------REFERENCES (cited in code using the acronyms shown)------------------------------------
 *
 *    DGIS: Devaux C, Grandjean P, Ishiguro Y, Siewert CE, 1979, 
 *              On Multi-Region Problems in Radiative Transfer, Astrophys. Space Sci. 62, 225-233
 *      GS: Garcia RDM, Siewert CE, 1985, Benchmark Results in Radiative Transfer,
 *              Transport Theory and Statistical Physics 14, 437-483
 *      KS: Kylling A, Stamnes K, 1992, Efficient yet accurate solution of the linear transport
 *              equation in the presence of internal sources: The exponential-linear-in-depth
 *              approximation, J. Comp. Phys., 102, 265-276
 *       L: Lenoble J, ed, 1985:  Radiative Transfer in Absorbing
 *              and Scattering Atmospheres: Standard Computational Procedures, Deepak Publishing, Hampton, Virginia
 *      NT: Nakajima T, Tanaka M, 1988,  Algorithms for Radiative Intensity Calculations in 
 *              Moderately Thick Atmospheres Using a Truncation Approximation, J.Q.S.R.T. 40, 51-69
 *      OS: Ozisik M, Shouman S, 1980,  Source Function Expansion Method for Radiative Transfer in a Two-Layer
 *              Slab, J.Q.S.R.T. 24, 441-449
 *      SS: Stamnes K, Swanson R, 1981,  A New Look at the Discrete Ordinate Method for Radiative Transfer
 *              Calculations in Anisotropically Scattering Atmospheres, J. Atmos. Sci. 38, 387-399
 *      SD: Stamnes K, Dale H, 1981, A New Look at the Discrete Ordinate Method for Radiative Transfer
 *              Calculations in Anisotropically Scattering Atmospheres. II: Intensity Computations,
 *              J. Atmos. Sci. 38, 2696-2706
 *      S1: Stamnes K, 1982, On the Computation of Angular Distributions of Radiation in Planetary
 *              Atmospheres, J.Q.S.R.T. 28, 47-51
 *      S2: Stamnes K, 1982, Reflection and Transmission by a Vertically Inhomogeneous Planetary Atmosphere,
 *              Planet. Space Sci. 30, 727-732
 *      SC: Stamnes K, Conklin P, 1984, A New Multi-Layer Discrete Ordinate Approach to Radiative Transfer
 *              in Vertically Inhomogeneous Atmospheres, J.Q.S.R.T. 31, 273-282
 *      SW: Sweigart A, 1970, Radiative Transfer in Atmospheres Scattering According to the Rayleigh Phase Function
 *              with Absorption, The Astrophysical Journal Supplement Series 22, 1-80
 *    STWJ: Stamnes K, Tsay SC, Wiscombe W, Jayaweera K, 1988, A Numerically Stable Algorithm for
 *              Discrete-Ordinate-Method Radiative Transfer in Multiple Scattering and Emitting Layered Media,
 *              Appl. Opt. 27, 2502-2509
 *    STWL: Stamnes K, Tsay SC, Wiscombe W, Laszlo I: A General-Purpose Numerically Stable Computer
 *              Code for Discrete-Ordinate-Method Radiative Transfer in Scattering and Emitting Layered Media,
 *              DISORT Report v1.1 (2000)
 * VH1,VH2: Van de Hulst, H.C., 1980: Multiple Light Scattering, Tables, Formulas and Applications, Volumes 1 and 2,
 *              Academic Press, New York.
 *       W: Wiscombe, W., 1977:  The Delta-M Method: Rapid Yet Accurate Radiative Flux Calculations, J. Atmos. Sci.
 *              34, 1408-1422
 *-----------------------------------------------------------------------------------------------------------------
 *
 * Rewritten in C by T. Dowling, Summer 2010.
 */

#include "cdisort.h"

/* 
 * Disort-specific shift macros.
 * Using unit-offset shift macros to match Fortran version
 */
#undef  DTAUC
#define DTAUC(lc)  ds.dtauc[lc-1]
#undef  PHI
#define PHI(j)     ds.phi[j-1]
#undef  PMOM
#define PMOM(k,lc) ds.pmom[k+(lc-1)*(ds.nmom_nstr+1)]
#undef  SSALB
#define SSALB(lc)  ds.ssalb[lc-1]
#undef  TEMPER
#define TEMPER(lc) ds.temper[lc]
#undef  UMU
#define UMU(iu)    ds.umu[iu-1]
#undef  UTAU
#define UTAU(lu)   ds.utau[lu-1]

/*
 * Disotest-specific shift macros
 */
#undef  GOODUU
#define GOODUU(iu,lu,j) good.uu[iu-1+(lu-1+(j-1)*ds.ntau)*ds.numu]

/*========================== main() ======================================*/

int main(int argc, char* argv[])
{
  /*
   * Comment out any tests you do not wish to run.
   * See the individual functions for a description of each test.
   */

  disort_test01();
  disort_test02();
  disort_test03();
  disort_test04();
  disort_test05();
  disort_test06(); 
  disort_test07();
  disort_test08();
  disort_test09(); 
  disort_test10();
  disort_test11();
  disort_test12();
  disort_test13(); 
  /*
   * Test 14 compares a 4-stream disort() to twostr(), and is not part of
   * the original Fortran distribution.
   */ 
  disort_test14();

  return 0;
}

/*========================== end of main() ===============================*/

/*========================== disort_test01() =============================*/

/**********************************************************************
 ****  Test Problem 1:  Isotropic Scattering                       ****
 ****  (Compare to Ref. VH1, Table 12)                             ****
 **********************************************************************/

void disort_test01(void)
{

  const int
    ncase = 6;
  char
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas;
  extern void
    c_disort(),c_errmsg(),c_getmom();

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.planck = FALSE;
  ds.flag.onlyfl = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;

  ds.nlyr   = 1;
  ds.nstr   = 16;
  ds.nphase = ds.nstr;
  ds.nmom   = ds.nstr;
  ds.nphi   = 1;
  ds.numu   = 6;
  ds.ntau   = 2;

  ds.flag.brdf_type = BRDF_NONE;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&good);
  
  c_getmom(ISOTROPIC,0.,ds.nmom,ds.pmom);

  UTAU(1)   =  0.;
  UMU(1)    = -1.;
  UMU(2)    = -0.5;
  UMU(3)    = -0.1;
  UMU(4)    =   .1;
  UMU(5)    =   .5;
  UMU(6)    =  1.;
  PHI(1)    =  0.;

  ds.bc.umu0   =   .1;
  ds.bc.phi0   =  0.;
  ds.bc.albedo =  0.;
  ds.bc.fluor  =  0.;

  for (icas = 1; icas <= ncase; icas++) {
    switch(icas) {
      case 1:
        UTAU(2)  = .03125;
        SSALB(1) = .2;

        ds.bc.fbeam = M_PI/ds.bc.umu0;
        ds.bc.fisot = 0.0;

        /* Correct answers */ 
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=2.29844;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =7.94108E-02;
        good.rad[0].flup  =7.99451E-02,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =2.54067E+01,good.rad[1].dfdt  =1.86531E+01;
        GOODUU(1,1,1) =0.,        GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=1.17771E-01,GOODUU(5,1,1)=2.64170E-02,GOODUU(6,1,1)=1.34041E-02,
        GOODUU(1,2,1)=1.33826E-02,GOODUU(2,2,1)=2.63324E-02,GOODUU(3,2,1)=1.15898E-01,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
      case 2:
        UTAU(2)  =  .03125;
        SSALB(1) = 1.;

        ds.bc.fbeam = M_PI/ds.bc.umu0;
        ds.bc.fisot = 0.0;

        /* Correct answers */
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=2.29844;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =4.20233E-01;
        good.rad[0].flup  =4.22922E-01,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =0.,         good.rad[1].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=6.22884E-01,GOODUU(5,1,1)=1.39763E-01,GOODUU(6,1,1)=7.09192E-02,
        GOODUU(1,2,1)=7.08109E-02,GOODUU(2,2,1)=1.39337E-01,GOODUU(3,2,1)=6.13458E-01,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
      case 3:
        UTAU(2)  = .03125;
        SSALB(1) = .99;

        ds.bc.fbeam = 0.;
        ds.bc.fisot = 1.;

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.;
        good.rad[0].rfldn =3.14159,    good.rad[1].rfldn =3.04897;
        good.rad[0].flup  =9.06556E-02,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =6.66870E-02,good.rad[1].dfdt  =5.88936E-02;
        GOODUU(1,1,1)=1.,         GOODUU(2,1,1)=1.,         GOODUU(3,1,1)=1.,         GOODUU(4,1,1)=1.33177E-01,GOODUU(5,1,1)=2.99879E-02,GOODUU(6,1,1)=1.52233E-02,
        GOODUU(1,2,1)=9.84447E-01,GOODUU(2,2,1)=9.69363E-01,GOODUU(3,2,1)=8.63946E-01,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
      case 4:
        UTAU(2)  = 32.;
        SSALB(1) =   .2;

        ds.bc.fbeam = M_PI/ds.bc.umu0;
        ds.bc.fisot = 0.0;

        /* Correct answers */
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=0.;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =0.;
        good.rad[0].flup  =2.59686E-01,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =2.57766E+01,good.rad[1].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=2.62972E-01,GOODUU(5,1,1)=9.06967E-02,GOODUU(6,1,1)=5.02853E-02,
        GOODUU(1,2,1)=1.22980E-15,GOODUU(2,2,1)=1.30698E-17,GOODUU(3,2,1)=6.88840E-18,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
      case 5:
        UTAU(2)  = 32.;
        SSALB(1) =  1.;

        ds.bc.fbeam = M_PI/ds.bc.umu0;
        ds.bc.fisot = 0.0;

        /* Correct answers */
        good.rad[0].rfldir=3.14159,good.rad[1].rfldir=0.;
        good.rad[0].rfldn =0.,     good.rad[1].rfldn =6.76954E-02;
        good.rad[0].flup  =3.07390,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =0.,     good.rad[1].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=1.93321E+00,GOODUU(5,1,1)=1.02732E+00,GOODUU(6,1,1)=7.97199E-01,
        GOODUU(1,2,1)=2.71316E-02,GOODUU(2,2,1)=1.87805E-02,GOODUU(3,2,1)=1.16385E-02,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
      case 6:
        UTAU(2)  = 32.;
        SSALB(1) =   .99;

        ds.bc.fbeam = 0.;
        ds.bc.fisot = 1.;

        /* Correct answers */ 
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.;
        good.rad[0].rfldn =3.14159,    good.rad[1].rfldn =4.60048E-03;
        good.rad[0].flup  =2.49618,    good.rad[1].flup  =0.;
        good.rad[0].dfdt  =1.14239E-01,good.rad[1].dfdt  =7.93633E-05;
        GOODUU(1,1,1)=1.,         GOODUU(2,1,1)=1.,         GOODUU(3,1,1)=1.,         GOODUU(4,1,1)=8.77510E-01,GOODUU(5,1,1)=8.15136E-01,GOODUU(6,1,1)=7.52715E-01,
        GOODUU(1,2,1)=1.86840E-03,GOODUU(2,2,1)=1.26492E-03,GOODUU(3,2,1)=7.79280E-04,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
    }

    DTAUC(1) = UTAU(2);

    sprintf(ds.header,"Test Case No. 1%.1s:  Isotropic Scattering, Ref. VH1, Table 12:  b =%9.5f, a =%5.2f",abc+icas-1,UTAU(2),SSALB(1));

    c_disort_out_alloc(&ds,&out);
    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);
    c_disort_out_free(&ds,&out);
  }

  /* Free allocated memory. */
  c_disort_out_free(&ds,&good);
  c_disort_state_free(&ds);

  return;
}

/*========================== end of disort_test01() ======================*/

/*========================== disort_test02() =============================*/

/**********************************************************************
 ****  Test Problem 2:  Rayleigh Scattering, Beam Source           ****
 ****  (Compare To Ref. SW, Table 1)                               ****
 **********************************************************************/

void disort_test02(void)
{
  char
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas,iod,iss;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.planck = FALSE;
  ds.flag.onlyfl = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;


  ds.nstr   = 16;
  ds.nlyr   = 1;
  ds.nphase = ds.nstr;
  ds.nmom   = ds.nstr;
  ds.ntau   = 2;
  ds.numu   = 6;
  ds.nphi   = 1;

  ds.flag.brdf_type = BRDF_NONE;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&good);

  c_getmom(RAYLEIGH,0.,ds.nmom,ds.pmom);

  UTAU(1)   = 0.;

  UMU(1)    = -0.981986;
  UMU(2)    = -0.538263;
  UMU(3)    = -0.018014;
  UMU(4)    = 0.018014;
  UMU(5)    = 0.538263;
  UMU(6)    = 0.981986;

  PHI(1)    = 0.0;

  ds.bc.fbeam  = M_PI;
  ds.bc.umu0   = 0.080442;
  ds.bc.phi0   = 0.;
  ds.bc.fisot  = 0.;
  ds.bc.albedo = 0.;
  ds.bc.fluor  = 0.;

  icas = 0;
  for (iod = 1; iod <= 2; iod++) {
    if (iod == 1) {
      UTAU(2) = 0.2;
    }
    else {
      UTAU(2) = 5.0;
    }
    DTAUC(1) = UTAU(2);
    for (iss = 1; iss <= 2; iss++) {
      if(iss == 1) {
        SSALB(1) = 0.5;
      }
      else {
        SSALB(1) = 1.0;
      }
                
      icas++;

      switch(icas) {
        case 1:
          /* Correct answers */
          good.rad[0].rfldir=2.52716E-01,good.rad[1].rfldir=2.10311E-02;
          good.rad[0].rfldn =0.,         good.rad[1].rfldn =4.41791E-02;
          good.rad[0].flup  =5.35063E-02,good.rad[1].flup  =0.;
          good.rad[0].dfdt  =1.66570E+00,good.rad[1].dfdt  =1.89848E-01;
          GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=1.61796E-01,GOODUU(5,1,1)=2.11501E-02,GOODUU(6,1,1)=7.86713E-03,
          GOODUU(1,2,1)=7.71897E-03,GOODUU(2,2,1)=2.00778E-02,GOODUU(3,2,1)=2.57685E-02,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
        break;
        case 2:
          /* Correct answers */
          good.rad[0].rfldir=2.52716E-01,good.rad[1].rfldir=2.10311E-02;
          good.rad[0].rfldn =0.,         good.rad[1].rfldn =1.06123E-01;
          good.rad[0].flup  =1.25561E-01,good.rad[1].flup  =0.;
          good.rad[0].dfdt  =0.,         good.rad[1].dfdt  =0.;
          GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=3.47678E-01,GOODUU(5,1,1)=4.87120E-02,GOODUU(6,1,1)=1.89387E-02,
          GOODUU(1,2,1)=1.86027E-02,GOODUU(2,2,1)=4.64061E-02,GOODUU(3,2,1)=6.77603E-02,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
        break;
        case 3:
          /* Correct answers */
          good.rad[0].rfldir=2.52716E-01,good.rad[1].rfldir=2.56077E-28;
          good.rad[0].rfldn =0.,         good.rad[1].rfldn =2.51683E-04;
          good.rad[0].flup  =6.24730E-02,good.rad[1].flup  =0.;
          good.rad[0].dfdt  =1.67462E+00,good.rad[1].dfdt  =1.75464E-04;
          GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=1.62566E-01,GOODUU(5,1,1)=2.45786E-02,GOODUU(6,1,1)=1.01498E-02,
          GOODUU(1,2,1)=1.70004E-04,GOODUU(2,2,1)=3.97168E-05,GOODUU(3,2,1)=1.32472E-05,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
        break;
        case 4:
          /* Correct answers */
          good.rad[0].rfldir=2.52716E-01,good.rad[1].rfldir=0.;
          good.rad[0].rfldn =0.,         good.rad[1].rfldn =2.68008E-02;
          good.rad[0].flup  =2.25915E-01,good.rad[1].flup  =0.;
          good.rad[0].dfdt  =0.,         good.rad[1].dfdt  =0.;
          GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=3.64010E-01,GOODUU(5,1,1)=8.26993E-02,GOODUU(6,1,1)=4.92370E-02,
          GOODUU(1,2,1)=1.05950E-02,GOODUU(2,2,1)=7.69337E-03,GOODUU(3,2,1)=3.79276E-03,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
        break;
      }

      sprintf(ds.header,"Test Case No. 2%.1s, Rayleigh Scattering, Ref. SW, Table 1:  tau =%5.2f, mu0 =%9.6f, ss-albedo =%4.2f",
                        abc+icas-1,UTAU(2),ds.bc.umu0,SSALB(1));

      c_disort_out_alloc(&ds,&out);

      c_disort(&ds,&out);

      print_test(&ds,&out,&ds,&good);

      c_disort_out_free(&ds,&out);
    }
  }

  /* Free allocated memory */
  c_disort_out_free(&ds,&good);
  c_disort_state_free(&ds);

  return;
}

/*========================== end of disort_test02() ======================*/

/*========================== disort_test03() =============================*/

/**********************************************************************
 ****  Test Problem 3:  Henyey-Greenstein Scattering               ****
 ****  (Compare To Ref. VH2, Table 37)                             ****
 **********************************************************************/

void disort_test03(void)
{
  const int
    ncase = 2;
  char
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  double
    gg;
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.onlyfl = FALSE;
  ds.flag.planck = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;


  ds.nstr   = 16;
  ds.nlyr   = 1;
  ds.nphase = ds.nstr;
  ds.nmom   = 32;
  ds.ntau   = 2;
  ds.numu   = 6;
  ds.nphi   = 1;

  ds.flag.brdf_type = BRDF_NONE;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&good);

  gg = .75;
  c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

  UTAU(1)   = 0.;

  UMU(1)    = -1.0;
  UMU(2)    = -0.5;
  UMU(3)    = -0.1;
  UMU(4)    = 0.1;
  UMU(5)    = 0.5;
  UMU(6)    = 1.;

  SSALB(1)  = 1.;

  PHI(1)    = 0.;

  ds.bc.umu0   = 1.;
  ds.bc.phi0   = 0.;
  ds.bc.fbeam  = M_PI/ds.bc.umu0;
  ds.bc.fisot  = 0.;
  ds.bc.albedo = 0.;
  ds.bc.fluor  = 0.;

  for (icas = 1; icas <= ncase; icas++) {
    switch(icas) {
      case 1:
        UTAU(2) = 1.0;

        /* Correct answers */
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=1.15573;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =1.73849;
        good.rad[0].flup  =2.47374E-01,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =0.,         good.rad[1].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=1.51159E-01,GOODUU(5,1,1)=1.01103E-01,GOODUU(6,1,1)=3.95460E-02,
        GOODUU(1,2,1)=3.05855E+00,GOODUU(2,2,1)=2.66648E-01,GOODUU(3,2,1)=2.13750E-01,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
      case 2:
        UTAU(2) = 8.0;

        /* Correct answers */
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=1.05389E-03;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =1.54958;
        good.rad[0].flup  =1.59096E+00,good.rad[1].flup  =0.;
        good.rad[0].dfdt  =0.,         good.rad[1].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=3.79740E-01,GOODUU(5,1,1)=5.19598E-01,GOODUU(6,1,1)=4.93302E-01,
        GOODUU(1,2,1)=6.69581E-01,GOODUU(2,2,1)=4.22350E-01,GOODUU(3,2,1)=2.36362E-01,GOODUU(4,2,1)=0.,         GOODUU(5,2,1)=0.,         GOODUU(6,2,1)=0.;
      break;
    }

    DTAUC(1) = UTAU(2);

    sprintf(ds.header,"Test Case No. 3%.1s, Henyey-Greenstein Scattering, Ref. VH2, Table 37, g = %3.2f, b =%9.5f, a =%5.2f",
                      abc+icas-1,gg,UTAU(2),SSALB(1));

    c_disort_out_alloc(&ds,&out);

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    c_disort_out_free(&ds,&out);
  }

  /* Free allocated memory */
  c_disort_out_free(&ds,&good);
  c_disort_state_free(&ds);

  return;
}

/*========================== end of disort_test03() ======================*/

/*========================== disort_test04() =============================*/

/**********************************************************************
 ****  Test Problem 4:  Haze-L Scattering, Beam Source             ****
 ****  (Compare to Ref. GS, Tables 12-16)                          ****
 **********************************************************************/

void disort_test04(void)
{
  const int
    ncase = 3;
  char
    title[128],
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.planck = FALSE;
  ds.flag.onlyfl = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;

  ds.nstr   = 32;
  ds.nphase = ds.nstr;
  ds.nlyr   = 1;
  ds.nmom   = ds.nstr;
  ds.ntau   = 3;
  ds.numu   = 6;

  ds.bc.fbeam  = M_PI;
  ds.bc.phi0   = 0.;
  ds.bc.fisot  = 0.;
  ds.bc.albedo = 0.;
  ds.bc.fluor  = 0.;

  ds.flag.brdf_type = BRDF_NONE;

  for (icas = 1; icas <= ncase; icas++) {
    sprintf(title,"Test Case No. 4%.1s, Haze-L Scattering, Ref. GS, Table ",abc+icas-1);
    switch(icas) {
      case 1:
        ds.nphi = 1;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        c_getmom(HAZE_GARCIA_SIEWERT,0.,ds.nmom,ds.pmom);

        DTAUC(1)  = 1.;

        UTAU(1)   = 0.;
        UTAU(2)   = 0.5;
        UTAU(3)   = 1.;

        UMU(1)    = -1.;
        UMU(2)    = -0.5;
        UMU(3)    = -0.1;
        UMU(4)    = 0.1;
        UMU(5)    = 0.5;
        UMU(6)    = 1.;

        SSALB(1) = 1.0;
        PHI(1)   = 0.0;

        ds.bc.umu0  = 1.0;
        sprintf(ds.header,"%s 12",title);

        /* Correct answers */
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=1.90547,    good.rad[2].rfldir=1.15573;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =1.17401,    good.rad[2].rfldn =1.81264;
        good.rad[0].flup  =1.73223E-01,good.rad[1].flup  =1.11113E-01,good.rad[2].flup  =0.;
        good.rad[0].dfdt  =0.,         good.rad[1].dfdt  =0.,         good.rad[2].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=9.26837E-02,GOODUU(5,1,1)=6.59569E-02,GOODUU(6,1,1)=3.64755E-02,
        GOODUU(1,2,1)=2.51608E+00,GOODUU(2,2,1)=1.19287E-01,GOODUU(3,2,1)=1.34962E-01,GOODUU(4,2,1)=1.23887E-01,GOODUU(5,2,1)=4.02058E-02,GOODUU(6,2,1)=1.77746E-02,
        GOODUU(1,3,1)=3.37302E+00,GOODUU(2,3,1)=2.19835E-01,GOODUU(3,3,1)=1.56893E-01,GOODUU(4,3,1)=0.,         GOODUU(5,3,1)=0.,         GOODUU(6,3,1)=0.;
      break;
      case 2:
        ds.nphi = 1;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        c_getmom(HAZE_GARCIA_SIEWERT,0.,ds.nmom,ds.pmom);

        DTAUC(1) = 1.;

        UTAU(1)  = 0.;
        UTAU(2)  = 0.5;
        UTAU(3)  = 1.;

        UMU(1)   = -1.;
        UMU(2)   = -0.5;
        UMU(3)   = -0.1;
        UMU(4)   = 0.1;
        UMU(5)   = 0.5;
        UMU(6)   = 1.;

        SSALB(1) = 0.9;
        PHI(1)   = 0.0;

        ds.bc.umu0  = 1.0;
        sprintf(ds.header,"%s 13",title);

        /* Correct answers */
        good.rad[0].rfldir=3.14159,    good.rad[1].rfldir=1.90547,    good.rad[2].rfldir=1.15573;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =1.01517,    good.rad[2].rfldn =1.51554;
        good.rad[0].flup  =1.23665E-01,good.rad[1].flup  =7.88690E-02,good.rad[2].flup  =0.;
        good.rad[0].dfdt  =3.43724E-01,good.rad[1].dfdt  =3.52390E-01,good.rad[2].dfdt  =3.19450E-01;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=6.53056E-02,GOODUU(5,1,1)=4.55144E-02,GOODUU(6,1,1)=2.82693E-02,
        GOODUU(1,2,1)=2.24258E+00,GOODUU(2,2,1)=9.66049E-02,GOODUU(3,2,1)=9.61335E-02,GOODUU(4,2,1)=8.43278E-02,GOODUU(5,2,1)=2.79473E-02,GOODUU(6,2,1)=1.38835E-02,
        GOODUU(1,3,1)=2.97057E+00,GOODUU(2,3,1)=1.67698E-01,GOODUU(3,3,1)=1.08115E-01,GOODUU(4,3,1)=0.,         GOODUU(5,3,1)=0.,         GOODUU(6,3,1)=0.;
      break;
      case 3:
        ds.nphi = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        c_getmom(HAZE_GARCIA_SIEWERT,0.,ds.nmom,ds.pmom);

        DTAUC(1) = 1.;

        UTAU(1)  = 0.;
        UTAU(2)  = 0.5;
        UTAU(3)  = 1.;

        UMU(1)   = -1.;
        UMU(2)   = -0.5;
        UMU(3)   = -0.1;
        UMU(4)   = 0.1;
        UMU(5)   = 0.5;
        UMU(6)   = 1.;

        SSALB(1) =   0.9;

        PHI(1)   =   0.0;
        PHI(2)   =  90.0;
        PHI(3)   = 180.0;

        ds.bc.umu0  =   0.5;
        sprintf(ds.header,"%s 14-16",title);

        /* Correct answers */
        good.rad[0].rfldir=1.57080,    good.rad[1].rfldir=5.77864E-01,good.rad[2].rfldir=2.12584E-01;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =7.02764E-01,good.rad[2].rfldn =8.03294E-01;
        good.rad[0].flup  =2.25487E-01,good.rad[1].flup  =1.23848E-01,good.rad[2].flup  =0.;
        good.rad[0].dfdt  =3.85003E-01,good.rad[1].dfdt  =3.37317E-01,good.rad[2].dfdt  =2.16403E-01;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=8.70812E-01,GOODUU(5,1,1)=2.24960E-01,GOODUU(6,1,1)=2.27572E-02,
        GOODUU(1,2,1)=4.77016E-02,GOODUU(2,2,1)=3.02631E+00,GOODUU(3,2,1)=1.41195E+00,GOODUU(4,2,1)=6.97692E-01,GOODUU(5,2,1)=1.09130E-01,GOODUU(6,2,1)=9.32861E-03,
        GOODUU(1,3,1)=8.38488E-02,GOODUU(2,3,1)=2.70538E+00,GOODUU(3,3,1)=8.76523E-01,GOODUU(4,3,1)=0.,         GOODUU(5,3,1)=0.,         GOODUU(6,3,1)=0.,
        GOODUU(1,1,2)=0.,         GOODUU(2,1,2)=0.,         GOODUU(3,1,2)=0.,         GOODUU(4,1,2)=8.88117E-02,GOODUU(5,1,2)=5.77411E-02,GOODUU(6,1,2)=2.27572E-02,
        GOODUU(1,2,2)=4.77016E-02,GOODUU(2,2,2)=5.80971E-02,GOODUU(3,2,2)=1.04502E-01,GOODUU(4,2,2)=9.16071E-02,GOODUU(5,2,2)=2.95842E-02,GOODUU(6,2,2)=9.32861E-03,
        GOODUU(1,3,2)=8.38488E-02,GOODUU(2,3,2)=9.42187E-02,GOODUU(3,3,2)=8.95457E-02,GOODUU(4,3,2)=0.,         GOODUU(5,3,2)=0.,         GOODUU(6,3,2)=0.,
        GOODUU(1,1,3)=0.,         GOODUU(2,1,3)=0.,         GOODUU(3,1,3)=0.,         GOODUU(4,1,3)=6.98247E-02,GOODUU(5,1,3)=5.02877E-02,GOODUU(6,1,3)=2.27572E-02,
        GOODUU(1,2,3)=4.77016E-02,GOODUU(2,2,3)=2.58544E-02,GOODUU(3,2,3)=6.25954E-02,GOODUU(4,2,3)=5.91273E-02,GOODUU(5,2,3)=2.47702E-02,GOODUU(6,2,3)=9.32861E-03,
        GOODUU(1,3,3)=8.38488E-02,GOODUU(2,3,3)=3.99383E-02,GOODUU(3,3,3)=4.67155E-02,GOODUU(4,3,3)=0.,         GOODUU(5,3,3)=0.,         GOODUU(6,3,3)=0.;
      break;
    }

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    /* Free allocated memory */
    c_disort_out_free(&ds,&good);
    c_disort_out_free(&ds,&out);
    c_disort_state_free(&ds);
  }
}

/*========================== end of disort_test04() ======================*/

/*========================== disort_test05() =============================*/

/**********************************************************************
 ****  Test Problem 5:  Cloud C.1 Scattering, Beam Source          ****
 ****  (Compare to Ref. GS, Tables 19-20)                          ****
 **********************************************************************/

void disort_test05(void)
{
  const int
    ncase = 2;
  char
    title[128],
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd   = GENERAL_BC;
  ds.flag.usrtau  = TRUE;
  ds.flag.usrang  = TRUE;
  ds.flag.lamber  = TRUE;
  ds.flag.planck  = FALSE;
  ds.flag.onlyfl  = FALSE;
  ds.flag.quiet   = TRUE;
  ds.flag.spher   = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;


  ds.nstr   =  48;
  ds.nphase = ds.nstr;
  ds.nlyr   =   1;
  ds.nmom   = 299;
  ds.ntau   =   3;
  ds.numu   =   6;
  ds.nphi   =   1;

  ds.flag.brdf_type = BRDF_NONE;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&good);

  c_getmom(CLOUD_GARCIA_SIEWERT,0.,ds.nmom,ds.pmom);

  DTAUC(1)  = 64.;

  UMU(1)    = -1.;
  UMU(2)    = -0.5;
  UMU(3)    = -0.1;
  UMU(4)    =  0.1;
  UMU(5)    =  0.5;
  UMU(6)    =  1.;

  PHI(1)    =  0.;

  ds.bc.fbeam  =  M_PI;
  ds.bc.umu0   =  1.;
  ds.bc.phi0   =  0.;
  ds.bc.fisot  =  0.;
  ds.bc.albedo =  0.;
  ds.bc.fluor  =  0.;

  for (icas = 1; icas <= ncase; icas++) {
    sprintf(title,"Test Case No. 5%.1s, Cloud C.1 Scattering, Ref. GS, Table ",abc+icas-1);
    switch(icas) {
      case 1:
        UTAU(1)  =  0.0;
        UTAU(2)  = 32.0;
        UTAU(3)  = 64.0;
        SSALB(1) =  1.0;
        sprintf(ds.header,"%s 19",title);

        /* Correct answers */
        good.rad[0].rfldir=3.14159,good.rad[1].rfldir=3.97856E-14,good.rad[2].rfldir=5.03852E-28;
        good.rad[0].rfldn =0.,     good.rad[1].rfldn =2.24768,    good.rad[2].rfldn =4.79851E-01;
        good.rad[0].flup  =2.66174,good.rad[1].flup  =1.76783,    good.rad[2].flup  =0.;
        good.rad[0].dfdt  =0.,     good.rad[1].dfdt  =0.,         good.rad[2].dfdt  =0.;
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=0.,         GOODUU(3,1,1)=0.,         GOODUU(4,1,1)=4.58927E-01,GOODUU(5,1,1)=7.72983E-01,GOODUU(6,1,1)=1.07196E+00,
        GOODUU(1,2,1)=7.53662E-01,GOODUU(2,2,1)=6.96362E-01,GOODUU(3,2,1)=6.50541E-01,GOODUU(4,2,1)=6.27631E-01,GOODUU(5,2,1)=5.81809E-01,GOODUU(6,2,1)=5.24532E-01,
        GOODUU(1,3,1)=1.95230E-01,GOODUU(2,3,1)=1.31990E-01,GOODUU(3,3,1)=7.20655E-02,GOODUU(4,3,1)=0.,         GOODUU(5,3,1)=0.,         GOODUU(6,3,1)=0.;
      break;
      case 2:
        UTAU(1)  =  3.2;
        UTAU(2)  = 12.8;
        UTAU(3)  = 48.0;
        SSALB(1) =  0.9;
        sprintf(ds.header,"%s 20",title);

        /* Correct answers */
        good.rad[0].rfldir=1.28058E-01,good.rad[1].rfldir=8.67322E-06,good.rad[2].rfldir=4.47729E-21;
        good.rad[0].rfldn =1.74767,    good.rad[1].rfldn =2.33975E-01,good.rad[2].rfldn =6.38345E-05;
        good.rad[0].flup  =2.70485E-01,good.rad[1].flup  =3.74252E-02,good.rad[2].flup  =1.02904E-05;
        good.rad[0].dfdt  =3.10129E-01,good.rad[1].dfdt  =4.52671E-02,good.rad[2].dfdt  =1.25021E-05;
        GOODUU(1,1,1)=6.79623E+01,GOODUU(2,1,1)=2.21027E-01,GOODUU(3,1,1)=1.36619E-01,GOODUU(4,1,1)=1.14084E-01,GOODUU(5,1,1)=8.73870E-02,GOODUU(6,1,1)=8.81626E-02,
        GOODUU(1,2,1)=2.05706E-01,GOODUU(2,2,1)=4.92736E-02,GOODUU(3,2,1)=2.65449E-02,GOODUU(4,2,1)=2.02154E-02,GOODUU(5,2,1)=1.29661E-02,GOODUU(6,2,1)=9.51334E-03,
        GOODUU(1,3,1)=3.41286E-05,GOODUU(2,3,1)=1.39916E-05,GOODUU(3,3,1)=7.47039E-06,GOODUU(4,3,1)=5.65602E-06,GOODUU(5,3,1)=3.58245E-06,GOODUU(6,3,1)=2.57858E-06;
        break;
      break;
    }

    c_disort_out_alloc(&ds,&out);

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    c_disort_out_free(&ds,&out);
  }

  /* Free allocated memory */
  c_disort_out_free(&ds,&good);
  c_disort_state_free(&ds);

  return;
}

/*========================== end of disort_test05() ======================*/

/*========================== disort_test06() =============================*/

/***********************************************************************
 ****  Test Problem 6:  No Scattering, Increasingly Complex Sources ****
 ***********************************************************************/

void disort_test06(void)
{
  const int
    ncase = 8;
  char
    title[128],
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.onlyfl = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;

  ds.nstr   = 16;
  ds.nphase = ds.nstr;
  ds.nlyr   =  1;
  ds.numu   =  4;
  ds.nphi   =  1;
  ds.nmom   =  0;

  ds.wvnmlo     =     0.;
  ds.wvnmhi     = 50000.;

  ds.bc.fbeam   = 200.;
  ds.bc.umu0    =   0.5;
  ds.bc.phi0    =   0.;
  ds.bc.fisot   =   0.;
  ds.bc.temis   =   1.;
  ds.bc.fluor   =   0.;

  ds.flag.brdf_type = BRDF_NONE;

  for (icas = 1; icas <= ncase; icas++) {
    sprintf(title,"Test Case No. 6%.1s: No Scattering; Source = Beam",abc+icas-1);

    switch(icas) {
      case 1:
        /*
         * Transparent medium, beam source
         */
        ds.flag.lamber = TRUE;
        ds.flag.planck = FALSE;
	
        ds.ntau = 2;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1)  = 0.;

        UMU(1)    = -1.;
        UMU(2)    = -0.1;
        UMU(3)    =  0.1;
        UMU(4)    =  1.;

        PHI(1)    = 90.;

        UTAU(1)   = 0.0;
        UTAU(2)   = 0.0;

        DTAUC(1)  = 0.;
        ds.bc.albedo = 0.;

        sprintf(ds.header,"%s; Bottom Albedo = %.0f",title,ds.bc.albedo);

        /* Correct answers */
        good.rad[0].rfldir=100.,good.rad[1].rfldir=100.,
        good.rad[0].rfldn =0.,  good.rad[1].rfldn =0.,
        good.rad[0].flup  =0.,  good.rad[1].flup  =0.,
        good.rad[0].dfdt  =200.,good.rad[1].dfdt  =200.,
        GOODUU(1,1,1)=0.,GOODUU(2,1,1)=0.,GOODUU(3,1,1)=0.,GOODUU(4,1,1)=0.,
        GOODUU(1,2,1)=0.,GOODUU(2,2,1)=0.,GOODUU(3,2,1)=0.,GOODUU(4,2,1)=0.;
      break;
      case 2:
        /*
         * Add some optical depth
         */
        ds.flag.lamber = TRUE;
        ds.flag.planck = FALSE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1) = 0.;

        UMU(1)   = -1.;
        UMU(2)   = -0.1;
        UMU(3)   =  0.1;
        UMU(4)   =  1.;

        PHI(1)   = 90.;

        UTAU(1)  = 0.;
        UTAU(2)  = 0.5;
        UTAU(3)  = 1.;

        DTAUC(1) = 1.0;

        sprintf(ds.header,"%s; Bottom Albedo = %.0f",title,ds.bc.albedo);

        /* Correct answers */
        good.rad[0].rfldir=100.,good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn =0.,  good.rad[1].rfldn =0.,         good.rad[2].rfldn =0.;
        good.rad[0].flup  =0.,  good.rad[1].flup  =0.,         good.rad[2].flup  =0.;
        good.rad[0].dfdt  =200.,good.rad[1].dfdt  =7.35759E+01,good.rad[2].dfdt  =2.70671E+01;
        GOODUU(1,1,1)=0.,GOODUU(2,1,1)=0.,GOODUU(3,1,1)=0.,GOODUU(4,1,1)=0.,
        GOODUU(1,2,1)=0.,GOODUU(2,2,1)=0.,GOODUU(3,2,1)=0.,GOODUU(4,2,1)=0.,
        GOODUU(1,3,1)=0.,GOODUU(2,3,1)=0.,GOODUU(3,3,1)=0.,GOODUU(4,3,1)=0.;
      break;
      case 3:
        /*
         * Add some isotropic reflection
         */
        ds.flag.lamber = TRUE;
        ds.flag.planck = FALSE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1)  = 0.;

        UMU(1)    = -1.;
        UMU(2)    = -0.1;
        UMU(3)    =  0.1;
        UMU(4)    =  1.;

        PHI(1)    = 90.;

        UTAU(1)   = 0.;
        UTAU(2)   = 0.5;
        UTAU(3)   = 1.;

        DTAUC(1) = 1.0;

	ds.flag.brdf_type = BRDF_HAPKE;
        ds.bc.albedo = .5;

        sprintf(ds.header,"%s; Bottom Albedo=%3.1f Lambert",title,ds.bc.albedo);

        /* Correct answers */
        good.rad[0].rfldir=100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =0.,         good.rad[2].rfldn =0.;
        good.rad[0].flup  =1.48450E+00,good.rad[1].flup  =2.99914E+00,good.rad[2].flup  =6.76676E+00;
        good.rad[0].dfdt  =2.02010E+02,good.rad[1].dfdt  =7.79962E+01,good.rad[2].dfdt  = 4.06006E+01;
        GOODUU(1,1,1)=0.,GOODUU(2,1,1)=0.,GOODUU(3,1,1)=9.77882E-05,GOODUU(4,1,1)=7.92386E-01,
        GOODUU(1,2,1)=0.,GOODUU(2,2,1)=0.,GOODUU(3,2,1)=1.45131E-02,GOODUU(4,2,1)=1.30642E+00,
        GOODUU(1,3,1)=0.,GOODUU(2,3,1)=0.,GOODUU(3,3,1)=2.15393E+00,GOODUU(4,3,1)=2.15393E+00;
      break;
      case 4:
        /*
         * Use non-isotropic reflection
         */
        ds.flag.lamber = FALSE;
        ds.flag.planck = FALSE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1) = 0.;

        UMU(1)   = -1.;
        UMU(2)   = -0.1;
        UMU(3)   =  0.1;
        UMU(4)   =  1.;

        PHI(1)   = 90.;

        UTAU(1)  = 0.;
        UTAU(2)  = 0.5;
        UTAU(3)  = 1.;

        DTAUC(1) = 1.;

        sprintf(ds.header,"%s; Bottom Albedo = Non-Lambert",title);

        /* Correct answers */
        good.rad[0].rfldir=100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =0.,         good.rad[2].rfldn =0.;
        good.rad[0].flup  =6.70783E-01,good.rad[1].flup  =1.39084E+00,good.rad[2].flup  =3.31655E+00;
        good.rad[0].dfdt  =2.00936E+02,good.rad[1].dfdt  =7.57187E+01,good.rad[2].dfdt  =3.45317E+01;
        GOODUU(1,1,1)=0.,GOODUU(2,1,1)=0.,GOODUU(3,1,1)=6.80068E-05,GOODUU(4,1,1)=3.15441E-01,
        GOODUU(1,2,1)=0.,GOODUU(2,2,1)=0.,GOODUU(3,2,1)=1.00931E-02,GOODUU(4,2,1)=5.20074E-01,
        GOODUU(1,3,1)=0.,GOODUU(2,3,1)=0.,GOODUU(3,3,1)=1.49795E+00,GOODUU(4,3,1)=8.57458E-01;
      break;
      case 5:
        /*
         * Add some bottom-boundary emission
         */
        ds.flag.lamber = FALSE;
        ds.flag.planck = TRUE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1)  = 0.;

        UMU(1)    = -1.;
        UMU(2)    = -0.1;
        UMU(3)    =  0.1;
        UMU(4)    =  1.;

        PHI(1)    = 90.;

        UTAU(1)   = 0.;
        UTAU(2)   = 0.5;
        UTAU(3)   = 1.;

        DTAUC(1)  =   1.;

        TEMPER(0) =   0.;
        TEMPER(1) =   0.;

        ds.bc.ttemp  =   0.;
        ds.bc.btemp  = 300.;
        sprintf(ds.header,"%s, Bottom Emission; Bott Alb = Non-Lambert",title);

        /* Correct answers */
        good.rad[0].rfldir=100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn =0.,         good.rad[1].rfldn =0.,         good.rad[2].rfldn =0.;
        good.rad[0].flup  =7.95458E+01,good.rad[1].flup  =1.59902E+02,good.rad[2].flup  =3.56410E+02;
        good.rad[0].dfdt  =3.07079E+02,good.rad[1].dfdt  =3.07108E+02,good.rad[2].dfdt  =7.17467E+02;
        GOODUU(1,1,1)=0.,GOODUU(2,1,1)=0.,GOODUU(3,1,1)=4.53789E-03,GOODUU(4,1,1)=4.33773E+01,
        GOODUU(1,2,1)=0.,GOODUU(2,2,1)=0.,GOODUU(3,2,1)=6.73483E-01,GOODUU(4,2,1)=7.15170E+01,
        GOODUU(1,3,1)=0.,GOODUU(2,3,1)=0.,GOODUU(3,3,1)=9.99537E+01,GOODUU(4,3,1)=1.17912E+02;
      break;
      case 6:
        /*
         * Add some top-boundary diffuse incidence (prescribed + emitted)
         */
        ds.flag.lamber = FALSE;
        ds.flag.planck = TRUE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1)  = 0.;

        UMU(1)    = -1.;
        UMU(2)    = -0.1;
        UMU(3)    =  0.1;
        UMU(4)    =  1.;

        PHI(1)    = 90.;

        UTAU(1)   = 0.;
        UTAU(2)   = 0.5;
        UTAU(3)   = 1.;

        DTAUC(1)  =  1.;

        TEMPER(0) =  0.;
        TEMPER(1) =  0.;

        ds.bc.fisot  = 100./M_PI;
        ds.bc.ttemp  = 250.;
        ds.bc.btemp  = 300.;
        sprintf(ds.header,"%s, Bottom+Top Emission; Bott Alb = Non-Lambert",title);

        /* Correct answers */
        good.rad[0].rfldir=100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn =3.21497E+02,good.rad[1].rfldn =1.42493E+02,good.rad[2].rfldn =7.05305E+01;
        good.rad[0].flup  =8.27917E+01,good.rad[1].flup  =1.66532E+02,good.rad[2].flup  =3.71743E+02;
        good.rad[0].dfdt  =9.54523E+02,good.rad[1].dfdt  =5.27085E+02,good.rad[2].dfdt  =8.45341E+02;
        GOODUU(1,1,1)=1.02336E+02,GOODUU(2,1,1)=1.02336E+02,GOODUU(3,1,1)=4.80531E-03,GOODUU(4,1,1)=4.50168E+01,
        GOODUU(1,2,1)=6.20697E+01,GOODUU(2,2,1)=6.89532E-01,GOODUU(3,2,1)=7.13172E-01,GOODUU(4,2,1)=7.42191E+01,
        GOODUU(1,3,1)=3.76472E+01,GOODUU(2,3,1)=4.64603E-03,GOODUU(3,3,1)=1.05844E+02,GOODUU(4,3,1)=1.22368E+02;
      break;
      case 7:
        /*
         * Add some internal emission
         */
        ds.flag.lamber = FALSE;
        ds.flag.planck = TRUE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1)  = 0.;

        UMU(1)    =  -1.;
        UMU(2)    =  -0.1;
        UMU(3)    =   0.1;
        UMU(4)    =   1.;

        PHI(1)    =  90.;

        UTAU(1)   = 0.;
        UTAU(2)   = 0.5;
        UTAU(3)   = 1.;

        DTAUC(1)  =   1.;

        TEMPER(0) = 250.;
        TEMPER(1) = 300.;

        ds.bc.ttemp  = 250.;
        ds.bc.btemp  = 300.;
        sprintf(ds.header,"%s, Bottom+Top+Internal Emission; Bott Alb = Non-Lambert",title);

        /* Correct answers */
        good.rad[0].rfldir=100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir= 1.35335E+01;
        good.rad[0].rfldn =3.21497E+02,good.rad[1].rfldn =3.04775E+02,good.rad[2].rfldn = 3.63632E+02;
        good.rad[0].flup  =3.35292E+02,good.rad[1].flup  =4.12540E+02,good.rad[2].flup  = 4.41125E+02;
        good.rad[0].dfdt  =5.80394E+02,good.rad[1].dfdt  =1.27117E+02,good.rad[2].dfdt  =-1.68003E+02;
        GOODUU(1,1,1)=1.02336E+02,GOODUU(2,1,1)=1.02336E+02,GOODUU(3,1,1)=7.80733E+01,GOODUU(4,1,1)=1.16430E+02,
        GOODUU(1,2,1)=9.78748E+01,GOODUU(2,2,1)=1.01048E+02,GOODUU(3,2,1)=1.15819E+02,GOODUU(4,2,1)=1.34966E+02,
        GOODUU(1,3,1)=1.10061E+02,GOODUU(2,3,1)=1.38631E+02,GOODUU(3,3,1)=1.38695E+02,GOODUU(4,3,1)=1.40974E+02;
      break;
      case 8:
        /*
         * Increase the optical depth
         */
        ds.flag.lamber = FALSE;
        ds.flag.planck = TRUE;

        ds.ntau = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        SSALB(1)  = 0.;

        UMU(1)    = -1.;
        UMU(2)    = -0.1;
        UMU(3)    =  0.1;
        UMU(4)    =  1.;

        PHI(1)    = 90.;

        DTAUC(1)  =  10.0;

        TEMPER(0) = 250.0;
        TEMPER(1) = 300.0;

        UTAU(1)   =   0.0;
        UTAU(2)   =   1.0;
        UTAU(3)   =  10.0;

        ds.bc.ttemp  = 250.0;
        ds.bc.btemp  = 300.0;
        sprintf(ds.header,"%s, Bottom+Top+Internal Emission; Bott Alb = Non-Lambert",title);

        /* Correct answers */
        good.rad[0].rfldir=100.,       good.rad[1].rfldir=1.35335E+01,good.rad[2].rfldir= 2.06115E-07;
        good.rad[0].rfldn =3.21497E+02,good.rad[1].rfldn =2.55455E+02,good.rad[2].rfldn = 4.43444E+02;
        good.rad[0].flup  =2.37350E+02,good.rad[1].flup  =2.61130E+02,good.rad[2].flup  = 4.55861E+02;
        good.rad[0].dfdt  =4.23780E+02,good.rad[1].dfdt  =6.19828E+01,good.rad[2].dfdt  =-3.11658E+01;
        GOODUU(1,1,1)=1.02336E+02,GOODUU(2,1,1)=1.02336E+02,GOODUU(3,1,1)=7.12616E+01,GOODUU(4,1,1)=7.80736E+01,
        GOODUU(1,2,1)=8.49992E+01,GOODUU(2,2,1)=7.73186E+01,GOODUU(3,2,1)=7.88310E+01,GOODUU(4,2,1)=8.56423E+01,
        GOODUU(1,3,1)=1.38631E+02,GOODUU(2,3,1)=1.45441E+02,GOODUU(3,3,1)=1.44792E+02,GOODUU(4,3,1)=1.45163E+02;
      break;
    }

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    /* Free allocated memory */
    c_disort_out_free(&ds,&good);
    c_disort_out_free(&ds,&out);
    c_disort_state_free(&ds);
  }
  
  return;
}

/*========================== end of disort_test06() ======================*/

/*========================== disort_test07() =============================*/

/**********************************************************************
 ****  Test Problem 7:  Absorption + Scattering + All Possible     ****
 ****  Sources, Various Surface Reflectivities ( One Layer )       ****
 **** (Compare 7a,f Fluxes and Intensities to Ref. KS, Tables I-II ****
 **********************************************************************/

void disort_test07(void)
{
  const int
    ncase = 5;
  char
    title[128],
    abc[] = "abcdefghijklmnopqrstuvwxyz";
  double
    gg;
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;

  ds.nlyr        = 1;

  ds.flag.brdf_type = BRDF_NONE;

  for (icas = 1; icas <= ncase; icas++) {
    sprintf(title,"Test Case No. 7%.1s",abc+icas-1);
    switch(icas) {
      case 1:
        ds.flag.ibcnd  = GENERAL_BC;
        ds.flag.lamber = TRUE;
        ds.flag.onlyfl = TRUE;
        ds.flag.planck = TRUE;
	ds.flag.usrang = FALSE;
        ds.flag.quiet  = TRUE;
	ds.flag.intensity_correction = TRUE;
        ds.flag.old_intensity_correction = TRUE;

        ds.nstr   = 16;
	ds.nphase = ds.nstr;
        ds.nmom   = 16;
        ds.nphi   =  1;
        ds.ntau   =  2;
        ds.numu   =  2;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        DTAUC(1) = 1.;
        SSALB(1) = 0.1;

        gg = .05;
        c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

        TEMPER(0)   = 200.;
        TEMPER(1)   = 300.;

        ds.wvnmlo      = 300.;
        ds.wvnmhi      = 800.;

        UTAU(1)     =   0.;
        UTAU(2)     =   1.;

        UMU(1)      =  -1.;
        UMU(2)      =   1.;

        PHI(1)      =   0.;

        ds.bc.umu0     =   0.5;
        ds.bc.phi0     =   0.;
        ds.bc.fbeam    =   0.;
        ds.bc.fisot    =   0.;
        ds.bc.albedo   =   0.;
        ds.bc.ttemp    =   0.;
        ds.bc.btemp    =   0.;
        ds.bc.temis    =   1.;
	ds.bc.fluor    =   0.;
        sprintf(ds.header,"%s: Absorption + Scattering, Internal Thermal Sources; Ref. KS, Table I, tau = %3.1f, a = %3.1f, g = %4.2f",
               title,DTAUC(1),SSALB(1),gg);

        /* Correct answers */
        good.rad[0].rfldir= 0.,         good.rad[1].rfldir= 0.,
        good.rad[0].rfldn = 0.,         good.rad[1].rfldn = 1.21204E+02,
        good.rad[0].flup  = 8.62936E+01,good.rad[1].flup  = 0.,
        good.rad[0].dfdt  =-5.13731E+01,good.rad[1].dfdt  =-5.41036E+02;
      break;
      case 2:
        ds.flag.ibcnd  = GENERAL_BC;
        ds.flag.usrang = TRUE;
        ds.flag.lamber = TRUE;
        ds.flag.planck = TRUE;
        ds.flag.onlyfl = FALSE;

        ds.ntau = 2;
        ds.numu = 2;
        ds.nphi = 1;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        DTAUC(1) = 100.0;
        SSALB(1) =    .95;

        gg = .75;
        c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

        TEMPER(0)   =  200.;
        TEMPER(1)   =  300.;

        ds.wvnmlo      = 2702.99;
        ds.wvnmhi      = 2703.01;

        UTAU(1)     = 0.;
        UTAU(2)     = 100.;

        UMU(1)      = -1.;
        UMU(2)      =  1.;

        PHI(1)      = 0.;

        ds.bc.umu0     = 0.5;
        ds.bc.phi0     = 0.;
        ds.bc.fbeam    = 0.;
        ds.bc.fisot    = 0.;
        ds.bc.albedo   = 0.;
        ds.bc.ttemp    = 0.;
        ds.bc.btemp    = 0.;
        ds.bc.temis    = 1.;
	ds.bc.fluor    =   0.;
        sprintf(ds.header,"%s: Absorption + Scattering, Internal Thermal Sources; Ref. KS, Table II, tau = %5.1f, a = %4.2f, g = %4.2f",
                title,DTAUC(1),SSALB(1),gg);

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir= 0.,
        good.rad[0].rfldn =0.,         good.rad[1].rfldn = 2.07786E-05,
        good.rad[0].flup  =1.10949E-06,good.rad[1].flup  = 0.,
        good.rad[0].dfdt  =8.23219E-08,good.rad[1].dfdt  =-5.06461E-06,
        GOODUU(1,1,1)=0.,         GOODUU(2,1,1)=4.65744E-07,
        GOODUU(1,2,1)=7.52311E-06,GOODUU(2,2,1)=0.;
      break;
      case 3:
        ds.flag.ibcnd  = GENERAL_BC;
        ds.flag.usrang = TRUE;
        ds.flag.planck = TRUE;
        ds.flag.onlyfl = FALSE;
        ds.flag.lamber = TRUE;
        ds.flag.quiet  = TRUE;
	ds.flag.intensity_correction = TRUE;
        ds.flag.old_intensity_correction = TRUE;

        ds.nstr   =  12;
	ds.nphase = ds.nstr;
        ds.nmom   =  12;
        ds.ntau   =   3;
        ds.numu   =   4;
        ds.nphi   =   2;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        TEMPER(0) = 300.0;
        TEMPER(1) = 200.0;

        gg = .8;
        c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

        DTAUC(1)    =     1.;
        SSALB(1)    =     0.5;

        ds.wvnmlo   =     0.;
        ds.wvnmhi   = 50000.;

        UTAU(1)     =     0.;
        UTAU(2)     =     0.5;
        UTAU(3)     =     1.;

        UMU(1)      =    -1.;
        UMU(2)      =    -0.1;
        UMU(3)      =     0.1;
        UMU(4)      =     1.;

        PHI(1)      =     0.;
        PHI(2)      =    90.;

        ds.bc.fbeam    =   200.;
        ds.bc.umu0     =     0.5;
        ds.bc.phi0     =     0.;
        ds.bc.fisot    =   100.;
        ds.bc.ttemp    =   100.;
        ds.bc.btemp    =   320.;
        ds.bc.temis    =     1.;
        ds.bc.albedo   =     0.;
	ds.bc.fluor    =   0.;
        sprintf(ds.header,"%s: Absorption + Henyey-Greenstein Scattering, All Sources, Bottom Albedo = %.0f",title,ds.bc.albedo);

        /* Correct answers */
        good.rad[0].rfldir= 100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn = 3.19830E+02,good.rad[1].rfldn =3.54099E+02,good.rad[2].rfldn =3.01334E+02;
        good.rad[0].flup  = 4.29572E+02,good.rad[1].flup  =4.47018E+02,good.rad[2].flup  =5.94576E+02;
        good.rad[0].dfdt  =-8.04270E+01,good.rad[1].dfdt  =2.51589E+02,good.rad[2].dfdt  =7.15964E+02;
        GOODUU(1,1,1)=1.01805E+02,GOODUU(2,1,1)=1.01805E+02,GOODUU(3,1,1)=1.46775E+02,GOODUU(4,1,1)=1.49033E+02,
        GOODUU(1,2,1)=1.06583E+02,GOODUU(2,2,1)=1.28565E+02,GOODUU(3,2,1)=1.04464E+02,GOODUU(4,2,1)=1.59054E+02,
        GOODUU(1,3,1)=9.66519E+01,GOODUU(2,3,1)=8.65854E+01,GOODUU(3,3,1)=1.89259E+02,GOODUU(4,3,1)=1.89259E+02,
        GOODUU(1,1,2)=1.01805E+02,GOODUU(2,1,2)=1.01805E+02,GOODUU(3,1,2)=1.29641E+02,GOODUU(4,1,2)=1.49033E+02,
        GOODUU(1,2,2)=1.06583E+02,GOODUU(2,2,2)=1.06408E+02,GOODUU(3,2,2)=9.48418E+01,GOODUU(4,2,2)=1.59054E+02,
        GOODUU(1,3,2)=9.66519E+01,GOODUU(2,3,2)=7.49310E+01,GOODUU(3,3,2)=1.89259E+02,GOODUU(4,3,2)=1.89259E+02;
      break;
      case 4:
        ds.flag.lamber = TRUE;
        ds.flag.ibcnd  = GENERAL_BC;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        TEMPER(0) = 300.0;
        TEMPER(1) = 200.0;

        gg = .8;
        c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

        DTAUC(1)    =     1.;
        SSALB(1)    =     0.5;

        ds.wvnmlo   =     0.;
        ds.wvnmhi   = 50000.;

        UTAU(1)     =     0.;
        UTAU(2)     =     0.5;
        UTAU(3)     =     1.;

        UMU(1)      =    -1.;
        UMU(2)      =    -0.1;
        UMU(3)      =     0.1;
        UMU(4)      =     1.;

        PHI(1)      =     0.;
        PHI(2)      =    90.;

        ds.bc.albedo = 1.0;
	ds.bc.fluor  =   0.;
	ds.flag.brdf_type = BRDF_HAPKE;

        sprintf(ds.header,"%s: Absorption + Henyey-Greenstein Scattering, All Sources, Bottom Albedo = %.0f",title,ds.bc.albedo);

        /* Correct answers */
        good.rad[0].rfldir= 100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn = 3.19830E+02,good.rad[1].rfldn =3.50555E+02,good.rad[2].rfldn =2.92063E+02;
        good.rad[0].flup  = 3.12563E+02,good.rad[1].flup  =2.68126E+02,good.rad[2].flup  =3.05596E+02;
        good.rad[0].dfdt  =-1.68356E+02,good.rad[1].dfdt  =1.01251E+02,good.rad[2].dfdt  =4.09326E+02;
        GOODUU(1,1,1)=1.01805E+02,GOODUU(2,1,1)=1.01805E+02,GOODUU(3,1,1)=1.40977E+02,GOODUU(4,1,1)=9.62764E+01,
        GOODUU(1,2,1)=1.06203E+02,GOODUU(2,2,1)=1.23126E+02,GOODUU(3,2,1)=9.19545E+01,GOODUU(4,2,1)=8.89528E+01,
        GOODUU(1,3,1)=9.56010E+01,GOODUU(2,3,1)=7.25576E+01,GOODUU(3,3,1)=9.72743E+01,GOODUU(4,3,1)=9.72743E+01,
        GOODUU(1,1,2)=1.01805E+02,GOODUU(2,1,2)=1.01805E+02,GOODUU(3,1,2)=1.23843E+02,GOODUU(4,1,2)=9.62764E+01,
        GOODUU(1,2,2)=1.06203E+02,GOODUU(2,2,2)=1.00969E+02,GOODUU(3,2,2)=8.23318E+01,GOODUU(4,2,2)=8.89528E+01,
        GOODUU(1,3,2)=9.56010E+01,GOODUU(2,3,2)=6.09031E+01,GOODUU(3,3,2)=9.72743E+01,GOODUU(4,3,2)=9.72743E+01;
      break;
      case 5:
        ds.flag.ibcnd  = GENERAL_BC;
        ds.flag.lamber = FALSE;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        TEMPER(0) = 300.0;
        TEMPER(1) = 200.0;

        gg = .8;
        c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

        DTAUC(1)    =     1.;
        SSALB(1)    =     0.5;

        ds.wvnmlo   =     0.;
        ds.wvnmhi   = 50000.;

        UTAU(1)     =     0.;
        UTAU(2)     =     0.5;
        UTAU(3)     =     1.;

        UMU(1)      =    -1.;
        UMU(2)      =    -0.1;
        UMU(3)      =     0.1;
        UMU(4)      =     1.;

        PHI(1)      =     0.;
        PHI(2)      =    90.;

        sprintf(ds.header,"%s: Absorption + Henyey-Greenstein Scattering, All Sources, Bottom Albedo = BDR Function",title);

        /* Correct answers */
        good.rad[0].rfldir= 100.,       good.rad[1].rfldir=3.67879E+01,good.rad[2].rfldir=1.35335E+01;
        good.rad[0].rfldn = 3.19830E+02,good.rad[1].rfldn =3.53275E+02,good.rad[2].rfldn =2.99002E+02;
        good.rad[0].flup  = 4.04300E+02,good.rad[1].flup  =4.07843E+02,good.rad[2].flup  =5.29248E+02;
        good.rad[0].dfdt  =-9.98568E+01,good.rad[1].dfdt  =2.17387E+02,good.rad[2].dfdt  =6.38461E+02;
        GOODUU(1,1,1)=1.01805E+02,GOODUU(2,1,1)=1.01805E+02,GOODUU(3,1,1)=1.45448E+02,GOODUU(4,1,1)=1.38554E+02,
        GOODUU(1,2,1)=1.06496E+02,GOODUU(2,2,1)=1.27296E+02,GOODUU(3,2,1)=1.01395E+02,GOODUU(4,2,1)=1.45229E+02,
        GOODUU(1,3,1)=9.63993E+01,GOODUU(2,3,1)=8.29009E+01,GOODUU(3,3,1)=1.60734E+02,GOODUU(4,3,1)=1.71307E+02,
        GOODUU(1,1,2)=1.01805E+02,GOODUU(2,1,2)=1.01805E+02,GOODUU(3,1,2)=1.28281E+02,GOODUU(4,1,2)=1.38554E+02,
        GOODUU(1,2,2)=1.06496E+02,GOODUU(2,2,2)=1.05111E+02,GOODUU(3,2,2)=9.16726E+01,GOODUU(4,2,2)=1.45229E+02,
        GOODUU(1,3,2)=9.63993E+01,GOODUU(2,3,2)=7.11248E+01,GOODUU(3,3,2)=1.59286E+02,GOODUU(4,3,2)=1.71307E+02;
      break;
    }

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    /* Free allocated memory */
    c_disort_out_free(&ds,&good);
    c_disort_out_free(&ds,&out);
    c_disort_state_free(&ds);
  }

  return;
}

/*========================== end of disort_test07() ======================*/

/*========================== disort_test08() =============================*/

/**********************************************************************
 ****  Test Problem 8:  Absorbing/Isotropic-Scattering Medium      ****
 ****  With Two Computational Layers                               ****
 **** (Compare Fluxes To Ref. OS, Table 1)                         ****
 **********************************************************************/

void disort_test08(void)
{
  const int
    ncase = 3;
  disort_state
    ds;
  disort_output
    out,good;
  register int
    icas,lc;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.planck = FALSE;
  ds.flag.onlyfl = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;


  ds.nstr   = 8;
  ds.nphase = ds.nstr;
  ds.nlyr   = 2;
  ds.nmom   = ds.nstr;
  ds.ntau   = 3;
  ds.numu   = 4;
  ds.nphi   = 1;

  ds.flag.brdf_type = BRDF_NONE;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&good);

  for (lc = 1; lc <= ds.nlyr; lc++) {
    c_getmom(ISOTROPIC,0.,ds.nmom,&PMOM(0,lc));
  }

  UMU(1)    = -1.;
  UMU(2)    = -0.2;
  UMU(3)    =  0.2;
  UMU(4)    =  1.;

  PHI(1)    = 60.;

  ds.bc.fbeam  =  0.;
  ds.bc.fisot  =  1./M_PI;
  ds.bc.albedo =  0.;
  ds.bc.phi0   = 0.0;
  ds.bc.umu0   = 0.5;
  ds.bc.fluor  =   0.;

  for (icas = 1; icas <= ncase; icas++) {
    switch(icas) {
      case 1:
        DTAUC(1) = .25;
        DTAUC(2) = .25;
        SSALB(1) = .5;
        SSALB(2) = .3;
        UTAU(1)  = .0;
        UTAU(2)  = .25;
        UTAU(3)  = .5;
        sprintf(ds.header,"Test Case No. 8a:  Ref. OS, Table 1, Line 4 (Two Inhomogeneous Layers)");

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.,         good.rad[2].rfldir=0.;
        good.rad[0].rfldn =1.,         good.rad[1].rfldn =7.22235E-01,good.rad[2].rfldn =5.13132E-01;
        good.rad[0].flup  =9.29633E-02,good.rad[1].flup  =2.78952E-02,good.rad[2].flup  =0.;
        good.rad[0].dfdt  =1.12474E+00,good.rad[1].dfdt  =6.51821E-01,good.rad[2].dfdt  =5.63361E-01;
        GOODUU(1,1,1)=3.18310E-01,GOODUU(2,1,1)=3.18310E-01,GOODUU(3,1,1)=5.62566E-02,GOODUU(4,1,1)=1.94423E-02,
        GOODUU(1,2,1)=2.62711E-01,GOODUU(2,2,1)=1.36952E-01,GOODUU(3,2,1)=1.84909E-02,GOODUU(4,2,1)=5.52188E-03,
        GOODUU(1,3,1)=2.10014E-01,GOODUU(2,3,1)=5.60376E-02,GOODUU(3,3,1)=0.,         GOODUU(4,3,1)=0.;
      break;
      case 2:
        DTAUC(1) = .25;
        DTAUC(2) = .25;
        SSALB(1) = .8;
        SSALB(2) = .95;
        UTAU(1)  = .0;
        UTAU(2)  = .25;
        UTAU(3)  = .5;
        sprintf(ds.header,"Test Case No. 8b:  Ref. OS, Table 1, Line 1 (Two Inhomogeneous Layers)");

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.,         good.rad[2].rfldir=0.;
        good.rad[0].rfldn =1.,         good.rad[1].rfldn =7.95332E-01,good.rad[2].rfldn =6.50417E-01;
        good.rad[0].flup  =2.25136E-01,good.rad[1].flup  =1.26349E-01,good.rad[2].flup  =0.;
        good.rad[0].dfdt  =5.12692E-01,good.rad[1].dfdt  =3.56655E-01,good.rad[2].dfdt  =5.68095E-02;
        GOODUU(1,1,1)=3.18310E-01,GOODUU(2,1,1)=3.18310E-01,GOODUU(3,1,1)=1.23687E-01,GOODUU(4,1,1)=4.95581E-02,
        GOODUU(1,2,1)=2.77499E-01,GOODUU(2,2,1)=1.83950E-01,GOODUU(3,2,1)=8.35695E-02,GOODUU(4,2,1)=2.50575E-02,
        GOODUU(1,3,1)=2.40731E-01,GOODUU(2,3,1)=1.29291E-01,GOODUU(3,3,1)=0.,         GOODUU(4,3,1)=0.;
      break;
      case 3:
        DTAUC(1) = 1.;
        DTAUC(2) = 2.;
        SSALB(1) =  .8;
        SSALB(2) =  .95;
        UTAU(1)  = 0.;
        UTAU(2)  = 1.;
        UTAU(3)  = 3.;
        sprintf(ds.header,"Test Case No. 8c:  Ref. OS, Table 1, Line 13 (Two Inhomogeneous Layers)");

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.,         good.rad[2].rfldir=0.;
        good.rad[0].rfldn =1.,         good.rad[1].rfldn =4.86157E-01,good.rad[2].rfldn =1.59984E-01;
        good.rad[0].flup  =3.78578E-01,good.rad[1].flup  =2.43397E-01,good.rad[2].flup  =0.;
        good.rad[0].dfdt  =5.65095E-01,good.rad[1].dfdt  =2.76697E-01,good.rad[2].dfdt  =1.35679E-02;
        GOODUU(1,1,1)=3.18310E-01,GOODUU(2,1,1)=3.18310E-01,GOODUU(3,1,1)=1.49335E-01,GOODUU(4,1,1)=1.04766E-01,
        GOODUU(1,2,1)=1.89020E-01,GOODUU(2,2,1)=9.88158E-02,GOODUU(3,2,1)=9.65192E-02,GOODUU(4,2,1)=6.54445E-02,
        GOODUU(1,3,1)=6.84762E-02,GOODUU(2,3,1)=2.96698E-02,GOODUU(3,3,1)=0.,         GOODUU(4,3,1)=0.;
      break;
    }

    c_disort_out_alloc(&ds,&out);

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    c_disort_out_free(&ds,&out);
  }

  /* Free allocated memory */
  c_disort_out_free(&ds,&good);
  c_disort_state_free(&ds);

  return;
}

/*========================== end of disort_test08() ======================*/

/*========================== disort_test09() =============================*/

/**********************************************************************
 ****  Test Problem 9:  General Emitting/Absorbing/Scattering      ****
 ****  Medium with Every Computational Layer Different.            ****
 **** (Compare 9a,b Fluxes to Ref. DGIS, Tables VI-VII, beta = 0)  ****
 **********************************************************************/

void disort_test09(void)
{
  register int
    icas,lc,k;
  const int
    ncase = 1;
  double
    gg;
  disort_state
    ds;
  disort_output
    out,good;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.usrang = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.onlyfl = FALSE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;

  ds.nstr   = 8;
  ds.nphase = ds.nstr;
  ds.nlyr   = 6;
  ds.nmom   = 8;
  ds.ntau   = 5;
  ds.numu   = 4;
  ds.nphi   = 1;

  ds.bc.fbeam = 0.;
  ds.bc.fisot = 1./M_PI;
  ds.bc.phi0  = 0.0;
  ds.bc.umu0  = 0.5;
  ds.bc.fluor =   0.;

  ds.flag.brdf_type = BRDF_NONE;

  for (icas = 1; icas <= ncase; icas++) {
    switch(icas) {
      case 1:
        ds.flag.planck = FALSE;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        for (lc = 1; lc <= ds.nlyr; lc++) {
          DTAUC(lc) = (double)lc;
          SSALB(lc) = 0.6+(double)lc*0.05;
        }

        UTAU(1)  = 0.;
        UTAU(2)  = 1.05;
        UTAU(3)  = 2.1;
        UTAU(4)  = 6.;
        UTAU(5)  = 21.;

        UMU(1)   = -1.;
        UMU(2)   = -0.2;
        UMU(3)   =  0.2;
        UMU(4)   =  1.;

        PHI(1)   = 60.;

        for (lc = 1; lc <= ds.nlyr; lc++) {
          c_getmom(ISOTROPIC,0.,ds.nmom,&PMOM(0,lc));
        }

        ds.bc.albedo  = 0.;
        sprintf(ds.header,"Test Case No. 9a:  Ref. DGIS, Tables VI-VII, beta=l=0 (multiple inhomogeneous layers)");

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.,         good.rad[2].rfldir=0.,         good.rad[3].rfldir=0.,         good.rad[4].rfldir=0.;
        good.rad[0].rfldn =1.,         good.rad[1].rfldn =3.55151E-01,good.rad[2].rfldn =1.44265E-01,good.rad[3].rfldn =6.71445E-03,good.rad[4].rfldn =6.16968E-07;
        good.rad[0].flup  =2.27973E-01,good.rad[1].flup  =8.75098E-02,good.rad[2].flup  =3.61819E-02,good.rad[3].flup  =2.19291E-03,good.rad[4].flup  =0.;
        good.rad[0].dfdt  =8.82116E-01,good.rad[1].dfdt  =2.32366E-01,good.rad[2].dfdt  =9.33443E-02,good.rad[3].dfdt  =3.92782E-03,good.rad[4].dfdt  =1.02500E-07;
        GOODUU(1,1,1)=3.18310E-01,GOODUU(2,1,1)=3.18310E-01,GOODUU(3,1,1)=9.98915E-02,GOODUU(4,1,1)=5.91345E-02,
        GOODUU(1,2,1)=1.53507E-01,GOODUU(2,2,1)=5.09531E-02,GOODUU(3,2,1)=3.67006E-02,GOODUU(4,2,1)=2.31903E-02,
        GOODUU(1,3,1)=7.06614E-02,GOODUU(2,3,1)=2.09119E-02,GOODUU(3,3,1)=1.48545E-02,GOODUU(4,3,1)=9.72307E-03,
        GOODUU(1,4,1)=3.72784E-03,GOODUU(2,4,1)=1.08815E-03,GOODUU(3,4,1)=8.83316E-04,GOODUU(4,4,1)=5.94743E-04,
        GOODUU(1,5,1)=2.87656E-07,GOODUU(2,5,1)=1.05921E-07,GOODUU(3,5,1)=0.,         GOODUU(4,5,1)=0.;
      break;
      case 2:
        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        for (lc = 1; lc <= ds.nlyr; lc++) {
          DTAUC(lc) = (double)lc;
          SSALB(lc) = 0.6+(double)lc*0.05;
        }

        UTAU(1)  = 0.;
        UTAU(2)  = 1.05;
        UTAU(3)  = 2.1;
        UTAU(4)  = 6.;
        UTAU(5)  = 21.;

        UMU(1)   = -1.;
        UMU(2)   = -0.2;
        UMU(3)   =  0.2;
        UMU(4)   =  1.;

        PHI(1)   = 60.;

        PMOM(0,1) = 1.0;
        PMOM(1,1) = 2.00916/3.;
        PMOM(2,1) = 1.56339/5.;
        PMOM(3,1) = 0.67407/7.;
        PMOM(4,1) = 0.22215/9.;
        PMOM(5,1) = 0.04725/11.;
        PMOM(6,1) = 0.00671/13.;
        PMOM(7,1) = 0.00068/15.;
        PMOM(8,1) = 0.00005/17.;
        for (lc = 2; lc <= ds.nlyr; lc++) {
          for (k = 0; k <= ds.nmom; k++) {
            PMOM(k,lc) = PMOM(k,1);
          }
        }

        sprintf(ds.header,"Test Case No. 9b:  Ref. DGIS, Tables VI-VII, beta=0,l=%1d (multiple inhomogeneous layers)",ds.nmom);

        /* Correct answers */
        good.rad[0].rfldir=0.,         good.rad[1].rfldir=0.,         good.rad[2].rfldir=0.,         good.rad[3].rfldir=0.,         good.rad[4].rfldir=0.;
        good.rad[0].rfldn =1.,         good.rad[1].rfldn =4.52357E-01,good.rad[2].rfldn =2.36473E-01,good.rad[3].rfldn =2.76475E-02,good.rad[4].rfldn =7.41853E-05;
        good.rad[0].flup  =1.00079E-01,good.rad[1].flup  =4.52014E-02,good.rad[2].flup  =2.41941E-02,good.rad[3].flup  =4.16016E-03,good.rad[4].flup  =0.;
        good.rad[0].dfdt  =8.04577E-01,good.rad[1].dfdt  =2.55330E-01,good.rad[2].dfdt  =1.30976E-01,good.rad[3].dfdt  =1.36227E-02,good.rad[4].dfdt  =1.22022E-05;
        GOODUU(1,1,1)=3.18310E-01,GOODUU(2,1,1)=3.18310E-01,GOODUU(3,1,1)=7.39198E-02,GOODUU(4,1,1)=1.32768E-02,
        GOODUU(1,2,1)=1.96609E-01,GOODUU(2,2,1)=5.92369E-02,GOODUU(3,2,1)=3.00230E-02,GOODUU(4,2,1)=7.05566E-03,
        GOODUU(1,3,1)=1.15478E-01,GOODUU(2,3,1)=3.01809E-02,GOODUU(3,3,1)=1.52672E-02,GOODUU(4,3,1)=4.06932E-03,
        GOODUU(1,4,1)=1.46177E-02,GOODUU(2,4,1)=3.85590E-03,GOODUU(3,4,1)=2.38301E-03,GOODUU(4,4,1)=7.77890E-04,
        GOODUU(1,5,1)=3.37742E-05,GOODUU(2,5,1)=1.20858E-05,GOODUU(3,5,1)=0.,         GOODUU(4,5,1)=0.;
      break;
      case 3:
        ds.flag.planck = TRUE;

        ds.nphi = 3;

        /* Allocate memory */
        c_disort_state_alloc(&ds);
        c_disort_out_alloc(&ds,&out);
        c_disort_out_alloc(&ds,&good);

        for (lc = 1; lc <= ds.nlyr; lc++) {
          DTAUC(lc) = (double)lc;
          SSALB(lc) = 0.6+(double)lc*0.05;
        }

        UTAU(1)  = 0.;
        UTAU(2)  = 1.05;
        UTAU(3)  = 2.1;
        UTAU(4)  = 6.;
        UTAU(5)  = 21.;

        UMU(1)   = -1.;
        UMU(2)   = -0.2;
        UMU(3)   =  0.2;
        UMU(4)   =  1.;

        TEMPER(0) = 600.0;
        for (lc = 1; lc <= ds.nlyr; lc++) {
          gg = (double)lc/7.;
          c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,&PMOM(0,lc));
          ds.temper[lc] = 600.+(double)lc*10.;
        }

        PHI(1)     =  60.0;
        PHI(2)     = 120.0;
        PHI(3)     = 180.0;

        ds.wvnmlo     = 999.0;
        ds.wvnmhi     = 1000.0;

        ds.bc.fbeam   = M_PI;
        ds.bc.fisot   = 1.0;
        ds.bc.albedo  = 0.5;
        ds.bc.btemp   = 700.0;
        ds.bc.ttemp   = 550.0;
        ds.bc.temis   = 1.0;
        sprintf(ds.header,"Test Case No. 9c:  Generalization of 9A to include all possible complexity");

        /* Correct answers */
        good.rad[0].rfldir=1.57080E+00,good.rad[1].rfldir=1.92354E-01,good.rad[2].rfldir=2.35550E-02;good.rad[3].rfldir=9.65131E-06,good.rad[4].rfldir=9.03133E-19;
        good.rad[0].rfldn =6.09217E+00,good.rad[1].rfldn =4.97279E+00,good.rad[2].rfldn =4.46616E+00;good.rad[3].rfldn =4.22731E+00,good.rad[4].rfldn =4.73767E+00;
        good.rad[0].flup  =4.68414E+00,good.rad[1].flup  =4.24381E+00,good.rad[2].flup  =4.16941E+00;good.rad[3].flup  =4.30667E+00,good.rad[4].flup  =5.11524E+00;
        good.rad[0].dfdt  =3.49563E+00,good.rad[1].dfdt  =8.81206E-01,good.rad[2].dfdt  =3.50053E-01;good.rad[3].dfdt  =1.93471E-02,good.rad[4].dfdt  =7.15349E-02;
        GOODUU(1,1,1)=1.93920,GOODUU(2,1,1)=1.93920,GOODUU(3,1,1)=1.61855,GOODUU(4,1,1)=1.43872,
        GOODUU(1,2,1)=1.66764,GOODUU(2,2,1)=1.44453,GOODUU(3,2,1)=1.38339,GOODUU(4,2,1)=1.33890,
        GOODUU(1,3,1)=1.48511,GOODUU(2,3,1)=1.35009,GOODUU(3,3,1)=1.33079,GOODUU(4,3,1)=1.32794,
        GOODUU(1,4,1)=1.34514,GOODUU(2,4,1)=1.35131,GOODUU(3,4,1)=1.35980,GOODUU(4,4,1)=1.37918,
        GOODUU(1,5,1)=1.48927,GOODUU(2,5,1)=1.54270,GOODUU(3,5,1)=1.62823,GOODUU(4,5,1)=1.62823,
        GOODUU(1,1,2)=1.93920,GOODUU(2,1,2)=1.93920,GOODUU(3,1,2)=1.57895,GOODUU(4,1,2)=1.43872,
        GOODUU(1,2,2)=1.66764,GOODUU(2,2,2)=1.42925,GOODUU(3,2,2)=1.37317,GOODUU(4,2,2)=1.33890,
        GOODUU(1,3,2)=1.48511,GOODUU(2,3,2)=1.34587,GOODUU(3,3,2)=1.32921,GOODUU(4,3,2)=1.32794,
        GOODUU(1,4,2)=1.34514,GOODUU(2,4,2)=1.35129,GOODUU(3,4,2)=1.35979,GOODUU(4,4,2)=1.37918,
        GOODUU(1,5,2)=1.48927,GOODUU(2,5,2)=1.54270,GOODUU(3,5,2)=1.62823,GOODUU(4,5,2)=1.62823,
        GOODUU(1,1,3)=1.93920,GOODUU(2,1,3)=1.93920,GOODUU(3,1,3)=1.56559,GOODUU(4,1,3)=1.43872,
        GOODUU(1,2,3)=1.66764,GOODUU(2,2,3)=1.42444,GOODUU(3,2,3)=1.37034,GOODUU(4,2,3)=1.33890,
        GOODUU(1,3,3)=1.48511,GOODUU(2,3,3)=1.34469,GOODUU(3,3,3)=1.32873,GOODUU(4,3,3)=1.32794,
        GOODUU(1,4,3)=1.34514,GOODUU(2,4,3)=1.35128,GOODUU(3,4,3)=1.35979,GOODUU(4,4,3)=1.37918,
        GOODUU(1,5,3)=1.48927,GOODUU(2,5,3)=1.54270,GOODUU(3,5,3)=1.62823,GOODUU(4,5,3)=1.62823;
      break;
    }

    c_disort(&ds,&out);

    print_test(&ds,&out,&ds,&good);

    /* Free allocated memory */
    c_disort_out_free(&ds,&good);
    c_disort_out_free(&ds,&out);
    c_disort_state_free(&ds);
  }

  return;
}

/*========================== end of disort_test09() ======================*/

/*========================== disort_test10() =============================*/

/**********************************************************************
 ****  Test Problem 10: Compare ds.flag.usrang = TRUE vs. FALSE    ****
 ****  Take Problem 9c (our most general case) but only 4 Streams  ****
 **********************************************************************/

void disort_test10(void)
{
  register int
    lc;
  double
    gg;
  disort_state
    ds_good,ds_out;
  disort_output
    good,out;

  ds_good.accur = ds_out.accur = 0.;
  ds_good.flag.prnt[0]=ds_out.flag.prnt[0]=TRUE;
  ds_good.flag.prnt[3]=ds_out.flag.prnt[3]=FALSE;
  ds_good.flag.prnt[4]=ds_out.flag.prnt[4]=TRUE;

  ds_good.flag.ibcnd =ds_out.flag.ibcnd  = GENERAL_BC;
  ds_good.flag.usrtau=ds_out.flag.usrtau = TRUE;
  ds_good.flag.lamber=ds_out.flag.lamber = TRUE;
  ds_good.flag.planck=ds_out.flag.planck = TRUE;
  ds_good.flag.onlyfl=ds_out.flag.onlyfl = FALSE;
  ds_good.flag.quiet =ds_out.flag.quiet  = TRUE;
  ds_good.flag.spher =ds_out.flag.spher  = FALSE;
  ds_good.flag.general_source          =ds_out.flag.general_source           = FALSE;
  ds_good.flag.output_uum              =ds_out.flag.output_uum           = FALSE;
  ds_good.flag.intensity_correction=ds_out.flag.intensity_correction = TRUE;
  ds_good.flag.old_intensity_correction=ds_out.flag.old_intensity_correction = TRUE;


  ds_good.nstr=ds_out.nstr = 4;
  ds_good.nphase=ds_out.nphase = ds_out.nstr;
  ds_good.nlyr=ds_out.nlyr = 6;
  ds_good.nmom=ds_out.nmom = ds_good.nstr;
  ds_good.ntau=ds_out.ntau = 3;
  ds_good.nphi=ds_out.nphi = 2;

  ds_good.bc.fbeam =ds_out.bc.fbeam    = M_PI;
  ds_good.bc.umu0  =ds_out.bc.umu0     = 0.5;
  ds_good.bc.phi0  =ds_out.bc.phi0     = 0.;
  ds_good.bc.fisot =ds_out.bc.fisot    = 1.;
  ds_good.bc.albedo=ds_out.bc.albedo   = 0.5;
  ds_good.bc.btemp =ds_out.bc.btemp    =  700.;
  ds_good.bc.ttemp =ds_out.bc.ttemp    =  550.;
  ds_good.bc.temis =ds_out.bc.temis    =    1.;
  ds_good.bc.fluor =ds_out.bc.fluor    = 0.;

  ds_good.wvnmlo=ds_out.wvnmlo =  999.;
  ds_good.wvnmhi=ds_out.wvnmhi = 1000.;

  ds_good.flag.brdf_type = BRDF_NONE;

  /* 
   * Case 1
   */
  ds_good.flag.prnt[1]= TRUE;
  ds_good.flag.prnt[2]= TRUE;

  ds_good.flag.usrang = TRUE;
  ds_good.numu        =  4;

  /* Allocate memory */
  c_disort_state_alloc(&ds_good);
  c_disort_out_alloc(&ds_good,&good);

  ds_good.temper[0] = 600.;
  for (lc = 1; lc <= ds_good.nlyr; lc++) {
    ds_good.dtauc[lc-1] = (double)lc;
    ds_good.ssalb[lc-1] = .6+(double)lc*.05;

    gg = (double)lc/(ds_good.nlyr+1);
    c_getmom(HENYEY_GREENSTEIN,gg,ds_good.nmom,&ds_good.pmom[0+(lc-1)*(ds_good.nmom_nstr+1)]);
    ds_good.temper[lc] = 600.+(double)lc*10.;
  }

  ds_good.utau[1-1] =  0.;
  ds_good.utau[2-1] =  2.1;
  ds_good.utau[3-1] = 21.;

  ds_good.phi[1-1] =  60.;
  ds_good.phi[2-1] = 120.;

  ds_good.umu[1-1] = -0.788675129;
  ds_good.umu[2-1] = -0.211324871;
  ds_good.umu[3-1] =  0.211324871;
  ds_good.umu[4-1] =  0.788675129;
  sprintf(ds_good.header,"Test Case No. 10a:  like 9c, ds.flag.usrang = TRUE");

  c_disort(&ds_good,&good);

  /* 
   * Case 2
   */
  ds_out.flag.prnt[1] = FALSE;
  ds_out.flag.prnt[2] = FALSE;

  ds_out.flag.brdf_type = BRDF_NONE;

  ds_out.flag.usrang  = FALSE;
  ds_out.numu         = 0;

  /* Allocate memory */
  c_disort_state_alloc(&ds_out);
  c_disort_out_alloc(&ds_out,&out);

  ds_out.temper[0] = 600.;
  for (lc = 1; lc <= ds_out.nlyr; lc++) {
    ds_out.dtauc[lc-1] = (double)lc;
    ds_out.ssalb[lc-1] = .6+(double)lc*.05;

    gg = (double)lc/(ds_out.nlyr+1);
    c_getmom(HENYEY_GREENSTEIN,gg,ds_out.nmom,&ds_out.pmom[0+(lc-1)*(ds_out.nmom_nstr+1)]);
    ds_out.temper[lc] = 600.+(double)lc*10.;
  }

  ds_out.utau[1-1] =  0.;
  ds_out.utau[2-1] =  2.1;
  ds_out.utau[3-1] = 21.;

  ds_out.phi[1-1] =  60.;
  ds_out.phi[2-1] = 120.;

  sprintf(ds_out.header,"Test Case No. 10b:  like 9C, ds.flag.usrang = FALSE");

  c_disort(&ds_out,&out);

  print_test(&ds_out,&out,&ds_good,&good);

  /* Free allocated memory */
  c_disort_out_free(&ds_good,&good);
  c_disort_out_free(&ds_out,&out);
  c_disort_state_free(&ds_good);
  c_disort_state_free(&ds_out);
  
  return;
}

/*========================== end of disort_test10() ======================*/

/*========================== disort_test11() =============================*/

/**********************************************************************
 ****  Test Problem 11: Single-Layer vs. Multiple Layers           ****
 ****  11a: Results at user levels for one computational layer     ****
 ****  11b: Single layer of 11a subdivided into multiple           ****
 ****       computational layers at the 11a user levels            ****
 **********************************************************************/

void disort_test11(void)
{
  register int
    lc;
  disort_state
    ds_good,ds_out;
  disort_output
    good,out;

  ds_good.accur=ds_out.accur = 0.;
  ds_good.flag.prnt[0]=ds_out.flag.prnt[0]=TRUE;
  ds_good.flag.prnt[1]=ds_out.flag.prnt[1]=FALSE;
  ds_good.flag.prnt[2]=ds_out.flag.prnt[2]=FALSE;
  ds_good.flag.prnt[3]=ds_out.flag.prnt[3]=FALSE;
  ds_good.flag.prnt[4]=ds_out.flag.prnt[4]=TRUE;

  ds_good.flag.ibcnd =ds_out.flag.ibcnd  = GENERAL_BC;
  ds_good.flag.usrang=ds_out.flag.usrang = TRUE;
  ds_good.flag.lamber=ds_out.flag.lamber = TRUE;
  ds_good.flag.planck=ds_out.flag.planck = FALSE;
  ds_good.flag.onlyfl=ds_out.flag.onlyfl = FALSE;
  ds_good.flag.quiet=ds_out.flag.quiet   = TRUE;
  ds_good.flag.spher =ds_out.flag.spher  = FALSE;
  ds_good.flag.general_source          =ds_out.flag.general_source           = FALSE;
  ds_good.flag.output_uum              =ds_out.flag.output_uum           = FALSE;
  ds_good.flag.intensity_correction=ds_out.flag.intensity_correction = TRUE;
  ds_good.flag.old_intensity_correction=ds_out.flag.old_intensity_correction = TRUE;


  ds_good.nstr=ds_out.nstr = 16;
  ds_good.nphase = ds_out.nstr;
  ds_good.nmom=ds_out.nmom = ds_good.nstr;
  ds_good.numu=ds_out.numu = 4;
  ds_good.nphi=ds_out.nphi = 2;

  ds_good.bc.fbeam =ds_out.bc.fbeam  = 1.;
  ds_good.bc.umu0  =ds_out.bc.umu0   = 0.5;
  ds_good.bc.phi0  =ds_out.bc.phi0   = 0.;
  ds_good.bc.fisot =ds_out.bc.fisot  = 0.5/M_PI;
  ds_good.bc.albedo=ds_out.bc.albedo = 0.5;
  ds_good.bc.fluor =ds_out.bc.fluor    = 0.;

  ds_good.flag.brdf_type = BRDF_NONE;

  /* 
   * Case 1
   */
  ds_good.flag.prnt[1] = TRUE, ds_good.flag.prnt[2] = TRUE;
  ds_good.flag.usrtau  = TRUE;

  ds_good.nlyr = 1;
  ds_good.ntau = 4;

  /* Allocate memory */
  c_disort_state_alloc(&ds_good);
  c_disort_out_alloc(&ds_good,&good);

  ds_good.umu[1-1] = -1.0;
  ds_good.umu[2-1] = -0.1;
  ds_good.umu[3-1] =  0.1;
  ds_good.umu[4-1] =  1.;

  ds_good.phi[1-1] =  0.;
  ds_good.phi[2-1] = 90.;

  ds_good.dtauc[1-1] =  1.;
  ds_good.ssalb[1-1] =   .9;

  c_getmom(ISOTROPIC,0.,ds_good.nmom,ds_good.pmom);

  ds_good.utau[1-1] = 0.;
  ds_good.utau[2-1] =  .05;
  ds_good.utau[3-1] =  .5;
  ds_good.utau[4-1] = 1.;

  sprintf(ds_good.header,"Test Case No. 11a: One Isotropic-Scattering Layer");

  c_disort(&ds_good,&good);

  /* 
   * Case 2
   */
  ds_out.flag.prnt[1] = FALSE, ds_out.flag.prnt[2] = FALSE;
  ds_out.flag.usrtau  = FALSE;

  ds_out.flag.brdf_type = BRDF_NONE;

  ds_out.nlyr = 3;

  /* Allocate memory */
  c_disort_state_alloc(&ds_out);
  c_disort_out_alloc(&ds_out,&out);

  ds_out.umu[1-1] = -1.0;
  ds_out.umu[2-1] = -0.1;
  ds_out.umu[3-1] =  0.1;
  ds_out.umu[4-1] =  1.;

  ds_out.phi[1-1] =  0.;
  ds_out.phi[2-1] = 90.;

  for (lc = 1; lc <= ds_out.nlyr; lc++) {
    ds_out.dtauc[lc-1] = ds_good.utau[lc+1-1]-ds_good.utau[lc-1];
    ds_out.ssalb[lc-1] = 0.9;
    c_getmom(ISOTROPIC,0.,ds_out.nmom,&ds_out.pmom[0+(lc-1)*(ds_out.nmom_nstr+1)]);
  }
  
  sprintf(ds_out.header,"Test Case No. 11b: Same as 11a but treated as multiple layers");

  c_disort(&ds_out,&out);

  print_test(&ds_out,&out,&ds_good,&good);

  /* Free allocated memory */
  c_disort_out_free(&ds_good,&good);
  c_disort_out_free(&ds_out,&out);
  c_disort_state_free(&ds_good);
  c_disort_state_free(&ds_out);
  
  return;
}

/*========================== end of disort_test11() ======================*/

/*========================== disort_test12() =============================*/

/***********************************************************************
 ****  Test Problem 12: Test Absorption-Optical-Depth Shortcut      ****
 ****  compares cases where the DISORT shortcut for absorption      ****
 ****  optical depth > 10 is not used (12a), then is used (12b)     ****
 ****  (this shortcut is only employed when ds.flag.planck = FALSE) ****
 ***********************************************************************/

void disort_test12(void)
{
  register int
    lc;
  double
    gg;
  disort_state
    ds_good,ds_out;
  disort_output
    good,out;

  ds_good.accur=ds_out.accur = 0.;
  ds_good.flag.prnt[0]=ds_out.flag.prnt[0] = TRUE;
  ds_good.flag.prnt[1]=ds_out.flag.prnt[1] = FALSE;
  ds_good.flag.prnt[2]=ds_out.flag.prnt[2] = FALSE;
  ds_good.flag.prnt[3]=ds_out.flag.prnt[3] = FALSE;
  ds_good.flag.prnt[4]=ds_out.flag.prnt[4] = TRUE;

  ds_good.flag.ibcnd     = ds_out.flag.ibcnd     = GENERAL_BC;
  ds_good.flag.usrang    = ds_out.flag.usrang    = TRUE;
  ds_good.flag.lamber    = ds_out.flag.lamber    = TRUE;
  ds_good.flag.planck    = ds_out.flag.planck    = FALSE;
  ds_good.flag.onlyfl    = ds_out.flag.onlyfl    = FALSE;
  ds_good.flag.quiet     = ds_out.flag.quiet     = TRUE;
  ds_good.flag.spher     = ds_out.flag.spher     = FALSE;
  ds_good.flag.brdf_type = ds_out.flag.brdf_type = BRDF_NONE;
  ds_good.flag.general_source          =ds_out.flag.general_source           = FALSE;
  ds_good.flag.output_uum              =ds_out.flag.output_uum           = FALSE;
  ds_good.flag.intensity_correction=ds_out.flag.intensity_correction = TRUE;
  ds_good.flag.old_intensity_correction=ds_out.flag.old_intensity_correction = TRUE;

  ds_good.nstr=ds_out.nstr = 20;
  ds_good.nphase = ds_out.nstr;
  ds_good.nmom   = ds_out.nmom = ds_good.nstr;
  ds_good.numu   = ds_out.numu = 4;
  ds_good.nphi   = ds_out.nphi = 1;

  /*
   * Case 1
   */
  ds_good.flag.prnt[1] = TRUE, ds_good.flag.prnt[2] = TRUE;
  ds_good.flag.usrtau  = TRUE;

  ds_good.nlyr = 1;
  ds_good.ntau = 4;

  /* Allocate memory */
  c_disort_state_alloc(&ds_good);
  c_disort_out_alloc(&ds_good,&good);

  ds_good.utau[1-1] =  0.;
  ds_good.utau[2-1] = 10.;
  ds_good.utau[3-1] = 19.9;
  ds_good.utau[4-1] = 20.1;

  ds_good.umu[1-1] = -1.;
  ds_good.umu[2-1] = -0.1;
  ds_good.umu[3-1] =  0.1;
  ds_good.umu[4-1] = 1.;

  ds_good.phi[1-1] = 0.;

  ds_good.dtauc[1-1] = 20.1;
  ds_good.ssalb[1-1] =   .5;

  gg = .9;
  c_getmom(HENYEY_GREENSTEIN,gg,ds_good.nmom,ds_good.pmom);

  ds_good.bc.fbeam    = 1.;
  ds_good.bc.umu0     = 1.;
  ds_good.bc.phi0     = 0.;
  ds_good.bc.fisot    = 0.;
  ds_good.bc.albedo   = 1.;
  ds_good.bc.fluor    = 0.;


  sprintf(ds_good.header,"Test Case No. 12a:  Overhead Beam Striking Absorbing/Scattering Medium");

  c_disort(&ds_good,&good);

  /*
   * Case 2
   */
  ds_out.flag.prnt[1] = FALSE, ds_out.flag.prnt[2] = FALSE;
  ds_out.flag.usrtau  = FALSE;

  ds_out.nlyr = ds_good.ntau-1;

  /* Allocate memory */
  c_disort_state_alloc(&ds_out);
  c_disort_out_alloc(&ds_out,&out);

  ds_out.umu[1-1] = -1.;
  ds_out.umu[2-1] = -0.1;
  ds_out.umu[3-1] =  0.1;
  ds_out.umu[4-1] =  1.;

  ds_out.phi[1-1] = 0.;

  gg = .9;
  for (lc = 1; lc <= ds_out.nlyr; lc++) {
    ds_out.dtauc[lc-1] = ds_good.utau[lc+1-1]-ds_good.utau[lc-1];
    ds_out.ssalb[lc-1] = .5;
    c_getmom(HENYEY_GREENSTEIN,gg,ds_out.nmom,&ds_out.pmom[0+(lc-1)*(ds_out.nmom_nstr+1)]);
  }
  ds_out.bc.fbeam    = 1.;
  ds_out.bc.umu0     = 1.;
  ds_out.bc.phi0     = 0.;
  ds_out.bc.fisot    = 0.;
  ds_out.bc.albedo   = 1.;
  ds_out.bc.fluor    = 0.;

  sprintf(ds_out.header,"Test Case No. 12b: Same as 12a but uses shortcut for absorption optical depth > 10");

  c_disort(&ds_out,&out);

  print_test(&ds_out,&out,&ds_good,&good);

  /* Free allocated memory */
  c_disort_out_free(&ds_good,&good);
  c_disort_out_free(&ds_out,&out);
  c_disort_state_free(&ds_good);
  c_disort_state_free(&ds_out);
  
  return;
}

/*========================== end of disort_test12() ======================*/

/*========================== disort_test13() =============================*/

/**********************************************************************
 ****  Test Problem 13: Test shortcut for flux albedo, transmission ***
 **** ( shortcut gives flux albedo, transmission of entire medium  ****
 ****   as a function of sun angle )                               ****
 ****  13a,c = Shortcut;  13b,d = Brute Force Method               ****
 **********************************************************************/

void disort_test13(void)
{
  register int
    lc;
  double
    gg;
  disort_state
    ds;
  disort_output
    out;

  ds.accur = 0.;
  ds.flag.prnt[0]=TRUE, ds.flag.prnt[1]=FALSE, ds.flag.prnt[2]=FALSE, ds.flag.prnt[3]=FALSE, ds.flag.prnt[4]=TRUE;
  ds.flag.quiet  = TRUE;
  ds.flag.spher  = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.intensity_correction = TRUE;
  ds.flag.old_intensity_correction = TRUE;

  ds.nstr      = 16;
  ds.nphase    = ds.nstr;
  ds.nmom      = ds.nstr;
  ds.nphi      = 0;
  ds.bc.phi0   = 0.0;
  ds.bc.albedo = 0.5;
  ds.bc.fluor  = 0.;

  ds.flag.brdf_type = BRDF_NONE;

  /*
   * Case 1
   */
  ds.flag.prnt[1] = FALSE, ds.flag.prnt[3] = TRUE;

  ds.flag.ibcnd  = SPECIAL_BC;
  ds.flag.usrang = TRUE;
  ds.flag.onlyfl = FALSE;

  ds.nlyr = 1;
  ds.numu = 1;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&out);

  DTAUC(1) = 1.;
  SSALB(1) =  .99;

  UMU(1) = 0.5;

  gg = .8;
  c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

  ds.bc.fbeam  =  0.;

  sprintf(ds.header,"Test Case No. 13a:  Albedo and Transmissivity from Shortcut, Single Layer");

  c_disort(&ds,&out);

  /* Free allocated memory */
  c_disort_out_free(&ds,&out);
  c_disort_state_free(&ds);

  /*
   * Case 2
   */
  ds.flag.prnt[1] = TRUE, ds.flag.prnt[3] = FALSE;

  ds.flag.ibcnd  = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.lamber = TRUE;
  ds.flag.planck = FALSE;
  ds.flag.onlyfl = TRUE;
  ds.flag.usrang = FALSE;

  ds.numu = 1;
  ds.ntau = 2;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&out);

  DTAUC(1) = 1.;
  SSALB(1) =  .99;

  UMU(1) = 0.5;

  gg = .8;
  c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,ds.pmom);

  UTAU(1) = 0.;
  UTAU(2) = 1.;

  ds.bc.umu0  = 0.5;
  ds.bc.fbeam = 1./ds.bc.umu0;
  ds.bc.fisot = 0.;

  sprintf(ds.header,"Test Case No. 13b:  Albedo and Transmissivity by Regular Method, Single Layer");

  c_disort(&ds,&out);

  /* Free allocated memory */
  c_disort_out_free(&ds,&out);
  c_disort_state_free(&ds);

  /*
   * Case 3
   */
  ds.flag.prnt[1] = FALSE, ds.flag.prnt[3] = TRUE;

  ds.flag.ibcnd = SPECIAL_BC;
  ds.flag.onlyfl = FALSE;
  ds.flag.usrang = TRUE;

  ds.numu = 1;
  ds.nlyr = 2;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&out);

  for (lc = 1; lc <= ds.nlyr; lc++) {
    DTAUC(lc) = 1./ds.nlyr;
    gg = .8;
    c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,&PMOM(0,lc));
  }
  SSALB(1) = .99;
  SSALB(2) = .50;

  UMU(1)   = .5;

  sprintf(ds.header,"Test Case No. 13c:  Albedo and Transmissivity from Shortcut, Multiple Layer");

  c_disort(&ds,&out);

  /* Free allocated memory */
  c_disort_out_free(&ds,&out);
  c_disort_state_free(&ds);

  /*
   * Case 4
   */
  ds.flag.prnt[1] = TRUE, ds.flag.prnt[3] = FALSE;

  ds.flag.ibcnd = GENERAL_BC;
  ds.flag.usrtau = TRUE;
  ds.flag.onlyfl = TRUE;
  ds.flag.usrang = FALSE;

  ds.numu = 1;
  ds.ntau = 2;

  /* Allocate memory */
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds,&out);

  for (lc = 1; lc <= ds.nlyr; lc++) {
    DTAUC(lc) = 1./ds.nlyr;
    gg = .8;
    c_getmom(HENYEY_GREENSTEIN,gg,ds.nmom,&PMOM(0,lc));
  }
  SSALB(1) = .99;
  SSALB(2) = .50;

  UMU(1)   = .5;

  UTAU(1) = 0.;
  UTAU(2) = 1.;

  ds.bc.umu0  = 0.5;
  ds.bc.fbeam = 1./ds.bc.umu0;
  ds.bc.fisot = 0.;

  sprintf(ds.header,"Test Case No. 13d:  Albedo and Transmissivity by Regular Method, Multiple Layer");

  c_disort(&ds,&out);

  /* Free allocated memory */
  c_disort_out_free(&ds,&out);
  c_disort_state_free(&ds);

  return;
}

/*========================== end of disort_test13() ======================*/

/*========================== disort_test14() =============================*/

/**********************************************************************
 ****  Test Problem 14: Like Test 10b, but compare twostr() to     ****
 ****  disort() with 4 streams.                                    ****
 **********************************************************************/

void disort_test14(void)
{
  register int
    lc;
  double
     radius,
    *gg;
  int
    deltam,
    ierror[TWOSTR_NERR];
  disort_state
    ds_ds,twostr_ds;
  disort_output
    ds_out,twostr_out;

  ds_ds.accur = 0.;
  ds_ds.flag.prnt[0]=TRUE;
  ds_ds.flag.prnt[1]=TRUE;
  ds_ds.flag.prnt[2]=FALSE;
  ds_ds.flag.prnt[3]=FALSE;
  ds_ds.flag.prnt[4]=TRUE;
  ds_ds.flag.quiet  = TRUE;
  ds_ds.flag.spher  = FALSE;
  ds_ds.flag.general_source           = FALSE;
  ds_ds.flag.output_uum = FALSE;
  ds_ds.flag.intensity_correction = TRUE;
  ds_ds.flag.old_intensity_correction = TRUE;


  ds_ds.flag.ibcnd  = GENERAL_BC;
  ds_ds.flag.usrtau = TRUE;
  ds_ds.flag.lamber = TRUE;
  ds_ds.flag.planck = TRUE;
  ds_ds.flag.onlyfl = TRUE;

  ds_ds.nstr = 4;
  ds_ds.nlyr = 6;
  ds_ds.nphase = ds_ds.nstr;
  ds_ds.nmom = ds_ds.nstr;
  ds_ds.ntau = 3;
  ds_ds.nphi = 2;

  ds_ds.bc.fbeam  = twostr_ds.bc.fbeam    = M_PI;
  ds_ds.bc.umu0   = twostr_ds.bc.umu0     =    0.5;
  ds_ds.bc.phi0   = twostr_ds.bc.phi0     =    0.;
  ds_ds.bc.fisot  = twostr_ds.bc.fisot    =    1.;
  ds_ds.bc.albedo = twostr_ds.bc.albedo   =    0.5;
  ds_ds.bc.fluor  = twostr_ds.bc.fluor    =    0.;
  ds_ds.bc.btemp  = twostr_ds.bc.btemp    =  700.;
  ds_ds.bc.ttemp  = twostr_ds.bc.ttemp    =  550.;
  ds_ds.bc.temis  = twostr_ds.bc.temis    =    1.;
  ds_ds.wvnmlo    = twostr_ds.wvnmlo      =  999.;
  ds_ds.wvnmhi    = twostr_ds.wvnmhi      = 1000.;

  ds_ds.flag.brdf_type = BRDF_NONE;

  /* 
   * Case 1: disort()
   */
  ds_ds.flag.usrang = FALSE;
  ds_ds.numu        =  0;

  /* Allocate memory */
  c_disort_state_alloc(&ds_ds);
  c_disort_out_alloc(&ds_ds,&ds_out);

  gg = c_dbl_vector(0,ds_ds.nlyr-1,"gg");

  ds_ds.temper[0] = 600.;
  for (lc = 1; lc <= ds_ds.nlyr; lc++) {
    ds_ds.dtauc[lc-1] = (double)lc;
    ds_ds.ssalb[lc-1] = .6+(double)lc*.05;

    gg[lc-1] = (double)lc/(ds_ds.nlyr+1);
    c_getmom(HENYEY_GREENSTEIN,gg[lc-1],ds_ds.nmom,&ds_ds.pmom[0+(lc-1)*(ds_ds.nmom_nstr+1)]);
    ds_ds.temper[lc] = 600.+(double)lc*10.;
  }

  free(gg);

  ds_ds.utau[1-1] =  0.;
  ds_ds.utau[2-1] =  2.1;
  ds_ds.utau[3-1] = 21.;

  ds_ds.phi[1-1] =  60.;
  ds_ds.phi[2-1] = 120.;

  sprintf(ds_ds.header,"Test Case No. 14a: disort() as in 10b");

  c_disort(&ds_ds,&ds_out);

  /* 
   * Case 2: twostr()
   */
  twostr_ds.flag.prnt[0] = TRUE;
  twostr_ds.flag.prnt[1] = FALSE;

  twostr_ds.flag.usrtau = TRUE;
  twostr_ds.flag.planck = TRUE;
  twostr_ds.flag.quiet  = FALSE;

  twostr_ds.ntau = 3;
  twostr_ds.nlyr = 6;

  deltam = TRUE;
  /*
   * If flag.spher = TRUE, also need to set radius and zd[].
   */
  twostr_ds.flag.spher = TRUE;
  radius = 6378.;

  /* Allocate memory */
  c_twostr_state_alloc(&twostr_ds);
  c_twostr_out_alloc(&twostr_ds,&twostr_out);
  /*
   * gg is needed as an input vector for twostr()
   */
  gg = c_dbl_vector(0,twostr_ds.nlyr-1,"gg");

  twostr_ds.temper[0] = 600.;
  if ( twostr_ds.flag.spher )   twostr_ds.zd[0]     = 10.*(double)(twostr_ds.nlyr+1);
  for (lc = 1; lc <= twostr_ds.nlyr; lc++) {
    twostr_ds.dtauc[lc-1] = (double)lc;
    twostr_ds.ssalb[lc-1] = .6+(double)lc*.05;

    gg[lc-1]             = (double)lc/(twostr_ds.nlyr+1);
    twostr_ds.temper[lc] = 600.+(double)lc*10.;
    if ( twostr_ds.flag.spher ) twostr_ds.zd[lc]     = 10.*(double)(twostr_ds.nlyr-lc+1);
  }

  twostr_ds.utau[1-1] =  0.;
  twostr_ds.utau[2-1] =  2.1;
  twostr_ds.utau[3-1] = 21.;

  sprintf(twostr_ds.header,"Test Case No. 14b: twostr() instead of disort()");

  c_twostr(&twostr_ds,&twostr_out,deltam,gg,ierror,radius);

  print_test(&twostr_ds,&twostr_out,&ds_ds,&ds_out);

  /* Free allocated memory */
  c_disort_out_free(&ds_ds,&ds_out);
  c_twostr_out_free(&twostr_ds,&twostr_out);
  c_disort_state_free(&ds_ds);
  c_twostr_state_free(&twostr_ds);
  free(gg);

  return;
}

/*========================== end of disort_test14() ======================*/

/*========================== print_test() ================================*/

/*
 * Print DISORT results and, directly beneath them, their ratios to the 
 * correct answers, calc/good, or 1.+calc in the case good = 0.;
 * print number of non-unit ratios that occur but try
 * to count just the cases where there is a real disagreement and not
 * those where flux or intensity are down at their noise level (defined as
 * 10^(-6) times their maximum value).  d(flux)/d(tau) is treated the
 * same as fluxes in this noise estimation even though it is a different
 * type of quantity (although with flux units).
 *
 * Correct input values are "good"; calculated input values are "calc".
 *
 * Fortran name: prtfin().
 */

#undef  BAD_RATIO
#define BAD_RATIO(r) ((r) < 0.99 || (r) > 1.01)

/* Using unit-offset shift macros to match Fortran version */

/* Disort-specific shift macros */
#undef  UMU
#define UMU(iu)    ds_good->umu[iu-1]
#undef  UTAU
#define UTAU(lu)   ds_good->utau[lu-1]

/* Disotest-specific shift macros */
#undef  GOODUU
#define GOODUU(iu,lu,j) good->uu[iu-1+(lu-1+(j-1)*ds_good->ntau)*ds_good->numu]
#undef  CALCUU
#define CALCUU(iu,lu,j) calc->uu[iu-1+(lu-1+(j-1)*ds_calc->ntau)*ds_calc->numu]

void print_test(disort_state  *ds_calc,
                disort_output *calc,
                disort_state  *ds_good,
                disort_output *good)
{
  register int
    iu,j,lu,numbad;
  extern void
    c_errmsg();
  extern double
    c_ratio();
  double
    flxmax,umax,fnoise,unoise,
    rat1,rat2,rat3,rat4,
    ratv[ds_good->nphi];

  flxmax = 0.0;
  for (lu = 0; lu < ds_good->ntau; lu++) {
    flxmax = MAX(MAX(MAX(flxmax,good->rad[lu].rfldir),good->rad[lu].rfldn),good->rad[lu].flup);
  }

  fnoise = 1.e-6*flxmax;
  if (flxmax <= 0.) {
    c_errmsg("print_test()--all fluxes zero or negative",DS_WARNING);
  }
  if (fnoise <= 0.) {
    c_errmsg("print_test()--all fluxes near underflowing",DS_WARNING);
  }

  numbad = 0;

  fprintf(stdout,"\n\n                  <-------------- FLUXES -------------->\n"
                 "    Optical       Downward       Downward         Upward    d(Net Flux)\n"
                 "      Depth         Direct        Diffuse        Diffuse    / d(Op Dep)\n");

  for (lu = 1; lu <= ds_good->ntau; lu++) {
    fprintf(stdout,"%11.4f%15.4e%15.4e%15.4e%15.4e\n",
                   UTAU(lu),calc->rad[lu-1].rfldir,calc->rad[lu-1].rfldn,calc->rad[lu-1].flup,calc->rad[lu-1].dfdt);

    fprintf(stdout,"%11.4f%15.4e%15.4e%15.4e%15.4e\n",
                   UTAU(lu),good->rad[lu-1].rfldir,good->rad[lu-1].rfldn,good->rad[lu-1].flup,good->rad[lu-1].dfdt);

    rat1 = c_ratio(calc->rad[lu-1].rfldir,good->rad[lu-1].rfldir);
    rat2 = c_ratio(calc->rad[lu-1].rfldn, good->rad[lu-1].rfldn);
    rat3 = c_ratio(calc->rad[lu-1].flup,  good->rad[lu-1].flup);
    rat4 = c_ratio(calc->rad[lu-1].dfdt,  good->rad[lu-1].dfdt);

    fprintf(stdout,"               (%9.4f)    (%9.4f)    (%9.4f)    (%9.4f)\n",rat1,rat2,rat3,rat4);

    /*
     * NOTE: In the original Fortran, for a/b, ratio() returns a huge number if b == 0., hence
     *       there is an extra conditional of the form fabs(output) > fnoise so that nearly-zero output
     *       will not be counted as bad. In contrast, this C version has ratio() returning 1.+a when b == 0.,
     *       hence the conditional involving fnoise is removed (the same applies to unoise below).
     */
    if(BAD_RATIO(rat1)) numbad++;
    if(BAD_RATIO(rat2)) numbad++;
    if(BAD_RATIO(rat3)) numbad++;
    if(BAD_RATIO(rat4)) numbad++;
  }

  if (!ds_good->flag.onlyfl) {
    /*
     * Print intensities
     */
    umax = 0.;
    for (j = 1; j <= ds_good->nphi; j++) {
      for (lu = 1; lu <= ds_good->ntau; lu++) {
        for (iu = 1; iu <= ds_good->numu; iu++) {
          umax = MAX(umax,GOODUU(iu,lu,j));
        }
      }
    }

    unoise = 1.e-6*umax;

    if (umax <= 0.) {
      c_errmsg("print_test()--all intensities zero or negative",DS_WARNING);
    }

    if (unoise <= 0.) {
      c_errmsg("print_test()--all intensities near underflowing",DS_WARNING);
    }

    fprintf(stdout,"\n\n ********  I N T E N S I T I E S  *********"
                   "\n\n             Polar   Azimuthal Angles (Degrees)"
                   "\n   Optical   Angle"
                   "\n     Depth  Cosine");
    for (j = 0; j < ds_good->nphi; j++) {
      fprintf(stdout,"%10.1f    ",ds_good->phi[j]);
    }
    fprintf(stdout,"\n");

    for (lu = 1; lu <= ds_good->ntau; lu++) {
      for (iu = 1; iu <= ds_good->numu; iu++) {
        if (iu == 1) {
          fprintf(stdout,"\n%10.3f%8.3f",UTAU(lu),UMU(iu));
          for (j = 1; j <= ds_good->nphi; j++) {
            fprintf(stdout,"%14.4e",CALCUU(iu,lu,j));
          }
          fprintf(stdout,"\n");
        }
        if (iu > 1) {
          fprintf(stdout,"          %8.3f",UMU(iu));
          for(j = 1; j <= ds_good->nphi; j++) {
            fprintf(stdout,"%14.4e",CALCUU(iu,lu,j));
          }
          fprintf(stdout,"\n");
        }
        for (j = 1; j <= ds_good->nphi; j++) {
          ratv[j-1] = c_ratio(CALCUU(iu,lu,j),GOODUU(iu,lu,j));
          /*
           * NOTE: This C version has the conditional fabs(output) > unoise removed;
           *       see note above regarding fnoise.
           */
          if(BAD_RATIO(ratv[j-1])) {
            numbad++;
          }
        }
        fprintf(stdout,"                  ");
        for (j = 1; j <= ds_good->nphi; j++) {
          fprintf(stdout,"   (%9.4f)",ratv[j-1]);
        }
        fprintf(stdout,"\n");
      }
    }
  }

  if (numbad > 0) {
    if (numbad == 1) {
      fprintf(stdout,"\n\n =============================================\n"
                     " ====  %4d  SERIOUSLY NON-UNIT RATIO     ====\n"
                     " =============================================\n",numbad);
    }
    else {
      fprintf(stdout,"\n\n =============================================\n"
                     " ====  %4d  SERIOUSLY NON-UNIT RATIOS    ====\n"
                     " =============================================\n",numbad);
    }
  }

  return;
}

/*========================== end of print_test() =========================*/

#undef GOODUU
#undef PMOM

/* * * * * * * * * * * * * * end of disotest.c  * * * * * * * * * * * * * */

