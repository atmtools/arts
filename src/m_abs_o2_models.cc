/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>
                                                                                
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.
                                                                                
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
                                                                                
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

/*! 
  \file   m_abs_o2_models.cc
  \author Oliver Lemke, Thomas Kuhn, Nikolay Koulev, Axel von Engeln
  \date   2004-01-14
  
  \brief  This file has functions for calculating the propagation
          tensor in the case of Zeeman effect for oxygen.
 
*/

// =================================================================================

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "abs_o2_models.h"
#include "arts.h"
#include "auto_md.h"
#include "complex.h"
#include "exceptions.h"
#include "math_funcs.h"
#include "matpackIII.h"
#include "messages.h"
#include "mystring.h"
#include "xml_io.h"

// =================================================================================

//-----------------------------------------------------------------------
//! HUMLIK_Faddeeva_Wells calculates the complex Faddeeva function (internal function)
/*! 
  This function calculates the complex Faddeeva function according to
  R. J. Wells, Rapid approximation to the Voigt/Faddeeva function and its 
  derivatives, JQSRT, vol.62, pp.29-48, 1999.
  This function must be called for each spectral line center frequency once and 
  calculates the complex Faddeeva function for the entire input frequency grid.

  -------------------------------------------------------------------
  DESCRIPTION OF THE INPUT VECTOR X:
  
                         F_i - F_o
  X_i =  SQRT(ln(2)) * --------------
                          GAMMA_D
       
  DESCRIPTION OF THE INPUT PARAMETER Y:
  
                          GAMMA_L 
  Y   =  SQRT(ln(2)) * --------------
                          GAMMA_D
  
  DESCRIPTION:
  F_i     :  ELEMENT OF THE INPUT FREQUENCY GRID OF CALCULATION  [Hz]
  F_o     :  LINE CENTER FREQUENCY OF SPECTRAL LINE              [Hz]
  GAMMA_D :  DOPPLER HALF WIDTH OF SPECTRAL LINE                 [Hz]
  GAMMA_L :  PRESSURE BROADENING HALF WIDTH OF SPECTRAL LINE     [Hz]
  SQRT(ln(2)) = 0.83255E0 
  -------------------------------------------------------------------
  
  \retval K      real part of the complex Faddeeva function (equal to the Voigt function)
  \retval L      imaginary part of the complex Faddeeva function (used for line mixing)
 
  \param X       X[i] =  SQRT(ln(2)) * ( (f_grid[i] - F_o)) / gamma_Doppler )
                 where F_o is the line center frequency of the specific line in question
  \param Y       Y    =  SQRT(ln(2)) * (gamma_Lorentz / gamma_Doppler)

  \author Thomas Kuhn
  \date   2003-11-30
*/

void HUMLIK_Faddeeva_Wells(  // Input
    ConstVectorView X,
    Numeric Y,
    // Output
    Vector& K,
    Vector& L) {
  /*
    To calculate the Faddeeva function with relative error less than 10^(-R).
    R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
    subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
  */

  const Index N = X.nelem(); /* Number of frequencies in the frequency grid */

  K.resize(
      N); /* resize the output vectors according to the size input vector X */
  L.resize(
      N); /* resize the output vectors according to the size input vector X */

  const Numeric R0 = 146.7; /*  Region boundaries  */
  const Numeric R1 = 14.67; /*  for R=4            */

  /* Constants */
  const Numeric RRTPI = 0.56418958; /* 1/SQRT(pi)          */
  const Numeric Y0 = 1.5;           /* for CPF12 algorithm */
  const Numeric Y0PY0 = Y0 + Y0;    /* for CPF12 algorithm */
  const Numeric Y0Q = Y0 * Y0;      /* for CPF12 algorithm */
  const Numeric C[] = {1.0117281,
                       -0.75197147,
                       0.012557727,
                       0.010022008,
                       -0.00024206814,
                       0.00000050084806};
  const Numeric S[] = {1.393237,
                       0.23115241,
                       -0.15535147,
                       0.0062183662,
                       0.000091908299,
                       -0.00000062752596};
  const Numeric T[] = {
      0.31424038, 0.94778839, 1.5976826, 2.2795071, 3.0206370, 3.8897249};

  /* Local variables */
  int RG1, RG2, RG3;                         /* y polynomial flags        */
  Numeric ABX, XQ, YQ, YRRTPI;               /* |x|, x^2, y^2, y/SQRT(pi) */
  Numeric XLIM0, XLIM1, XLIM2, XLIM3, XLIM4; /* |x| on region boundaries  */
  /* W4 temporary variables */
  Numeric A0 = 0.0, D0 = 0.0, D2 = 0.0, E0 = 0.0, E2 = 0.0, E4 = 0.0, H0 = 0.0,
          H2 = 0.0, H4 = 0.0, H6 = 0.0;
  Numeric P0 = 0.0, P2 = 0.0, P4 = 0.0, P6 = 0.0, P8 = 0.0, Z0 = 0.0, Z2 = 0.0,
          Z4 = 0.0, Z6 = 0.0, Z8 = 0.0;
  Numeric B1 = 0.0, F1 = 0.0, F3 = 0.0, F5 = 0.0, Q1 = 0.0, Q3 = 0.0, Q5 = 0.0,
          Q7 = 0.0;
  Numeric XP[5], XM[5], YP[5], YM[5]; /* CPF12 temporary values    */
  Numeric MQ[5], PQ[5], MF[5], PF[5];
  Numeric D, YF, YPY0, YPY0Q;

  /***** Start of executable code *****************************************/

  RG1 = 1; /* Set flags  */
  RG2 = 1;
  RG3 = 1;
  YQ = Y * Y;         /* y^2        */
  YRRTPI = Y * RRTPI; /* y/SQRT(pi) */

  /* Region boundaries when both K and L are required or when R<>4 */
  XLIM0 = R0 - Y;
  XLIM1 = R1 - Y;
  XLIM3 = 3.097 * Y - 0.45;

  /* For speed the following 3 lines should replace the 3 above if R=4 and L is not required */
  /* XLIM0 = 15100.0 + Y*(40.0 + Y*3.6)                                                      */
  /* XLIM1 = 164.0 - Y*(4.3 + Y*1.8)                                                         */
  /* XLIM3 = 5.76*YQ                                                                         */

  XLIM2 = 6.8 - Y;
  XLIM4 = 18.1 * Y + 1.65;
  if (Y <= 0.000001) /* When y<10^-6       */
  {
    XLIM1 = XLIM0; /* avoid W4 algorithm */
    XLIM2 = XLIM0;
  };

  /*..... */

  for (int I = 0; I < N; ++I) /* Loop over all points */
  {
    ABX = fabs(X[I]); /*! |x|  */
    XQ = ABX * ABX;   /* ! x^2 */
    if (ABX > XLIM0)  /* ! Region 0 algorithm */
    {
      K[I] = YRRTPI / (XQ + YQ);
      L[I] = K[I] * X[I] / Y;
    } else if (ABX > XLIM1) /* ! Humlicek W4 Region 1 */
    {
      if (RG1 != 0) /* ! First point in Region 1 */
      {
        RG1 = 0;
        A0 = YQ + 0.5; /* ! Region 1 y-dependents */
        D0 = A0 * A0;
        D2 = YQ + YQ - 1.0;
        B1 = YQ - 0.5;
      };
      D = RRTPI / (D0 + XQ * (D2 + XQ));
      K[I] = D * Y * (A0 + XQ);
      L[I] = D * X[I] * (B1 + XQ);
    } else if (ABX > XLIM2) /* ! Humlicek W4 Region 2 */
    {
      if (RG2 != 0) /* ! First point in Region 2 */
      {
        RG2 = 0;
        H0 = 0.5625 +
             YQ * (4.5 +
                   YQ * (10.5 + YQ * (6.0 + YQ))); /* ! Region 2 y-dependents */
        H2 = -4.5 + YQ * (9.0 + YQ * (6.0 + YQ * 4.0));
        H4 = 10.5 - YQ * (6.0 - YQ * 6.0);
        H6 = -6.0 + YQ * 4.0;
        E0 = 1.875 + YQ * (8.25 + YQ * (5.5 + YQ));
        E2 = 5.25 + YQ * (1.0 + YQ * 3.0);
        E4 = 0.75 * H6;
        F1 = -1.875 + YQ * (5.25 + YQ * (4.5 + YQ));
        F3 = 8.25 - YQ * (1.0 - YQ * 3.0);
        F5 = -5.5 + YQ * 3.0;
      };
      D = RRTPI / (H0 + XQ * (H2 + XQ * (H4 + XQ * (H6 + XQ))));
      K[I] = D * Y * (E0 + XQ * (E2 + XQ * (E4 + XQ)));
      L[I] = D * X[I] * (F1 + XQ * (F3 + XQ * (F5 + XQ)));
    } else if (ABX < XLIM3) /* ! Humlicek W4 Region 3 */
    {
      if (RG3 != 0) /* ! First point in Region 3 */
      {
        RG3 = 0;
        Z0 =
            272.1014 +
            Y * (1280.829 +
                 Y * (2802.870 +
                      Y * (3764.966 +
                           Y * (3447.629 +
                                Y * (2256.981 +
                                     Y * (1074.409 +
                                          Y * (369.1989 +
                                               Y * (88.26741 +
                                                    Y * (13.39880 +
                                                         Y))))))))); /*! Region 3 y-dependents*/
        Z2 = 211.678 +
             Y * (902.3066 +
                  Y * (1758.336 +
                       Y * (2037.310 +
                            Y * (1549.675 +
                                 Y * (793.4273 +
                                      Y * (266.2987 +
                                           Y * (53.59518 + Y * 5.0)))))));
        Z4 = 78.86585 +
             Y * (308.1852 +
                  Y * (497.3014 +
                       Y * (479.2576 +
                            Y * (269.2916 + Y * (80.39278 + Y * 10.0)))));
        Z6 = 22.03523 +
             Y * (55.02933 + Y * (92.75679 + Y * (53.59518 + Y * 10.0)));
        Z8 = 1.496460 + Y * (13.39880 + Y * 5.0);
        P0 = 153.5168 +
             Y * (549.3954 +
                  Y * (919.4955 +
                       Y * (946.8970 +
                            Y * (662.8097 +
                                 Y * (328.2151 +
                                      Y * (115.3772 +
                                           Y * (27.93941 +
                                                Y * (4.264678 +
                                                     Y * 0.3183291))))))));
        P2 = -34.16955 +
             Y * (-1.322256 +
                  Y * (124.5975 +
                       Y * (189.7730 +
                            Y * (139.4665 +
                                 Y * (56.81652 +
                                      Y * (12.79458 + Y * 1.2733163))))));
        P4 = 2.584042 +
             Y * (10.46332 +
                  Y * (24.01655 +
                       Y * (29.81482 + Y * (12.79568 + Y * 1.9099744))));
        P6 = -0.07272979 + Y * (0.9377051 + Y * (4.266322 + Y * 1.273316));
        P8 = 0.0005480304 + Y * 0.3183291;
        Q1 = 173.2355 +
             Y * (508.2585 +
                  Y * (685.8378 +
                       Y * (557.5178 +
                            Y * (301.3208 +
                                 Y * (111.0528 +
                                      Y * (27.62940 +
                                           Y * (4.264130 + Y * 0.3183291)))))));
        Q3 = 18.97431 +
             Y * (100.7375 +
                  Y * (160.4013 +
                       Y * (130.8905 +
                            Y * (55.88650 + Y * (12.79239 + Y * 1.273316)))));
        Q5 = 7.985877 +
             Y * (19.83766 + Y * (28.88480 + Y * (12.79239 + Y * 1.909974)));
        Q7 = 0.6276985 + Y * (4.264130 + Y * 1.273316);
      };
      D = 1.7724538 / (Z0 + XQ * (Z2 + XQ * (Z4 + XQ * (Z6 + XQ * (Z8 + XQ)))));
      K[I] = D * (P0 + XQ * (P2 + XQ * (P4 + XQ * (P6 + XQ * P8))));
      L[I] =
          D * X[I] * (Q1 + XQ * (Q3 + XQ * (Q5 + XQ * (Q7 + XQ * 0.3183291))));
    } else /* ! Humlicek CPF12 algorithm */
    {
      YPY0 = Y + Y0;
      YPY0Q = YPY0 * YPY0;
      K[I] = 0.000;
      L[I] = 0.000;
      for (int J = 0; J <= 5; ++J) {
        D = X[I] - T[J];
        MQ[J] = D * D;
        MF[J] = 1.0 / (MQ[J] + YPY0Q);
        XM[J] = MF[J] * D;
        YM[J] = MF[J] * YPY0;
        D = X[I] + T[J];
        PQ[J] = D * D;
        PF[J] = 1.0 / (PQ[J] + YPY0Q);
        XP[J] = PF[J] * D;
        YP[J] = PF[J] * YPY0;
        L[I] = L[I] + C[J] * (XM[J] + XP[J]) + S[J] * (YM[J] - YP[J]);
      };
      if (ABX <= XLIM4) /* ! Humlicek CPF12 Region I */
      {
        for (int J = 0; J <= 5; ++J) {
          K[I] = K[I] + C[J] * (YM[J] + YP[J]) - S[J] * (XM[J] - XP[J]);
        }
      } else /* ! Humlicek CPF12 Region II */
      {
        YF = Y + Y0PY0;
        for (int J = 0; J <= 5; ++J) {
          K[I] = K[I] +
                 (C[J] * (MQ[J] * MF[J] - Y0 * YM[J]) + S[J] * YF * XM[J]) /
                     (MQ[J] + Y0Q) +
                 (C[J] * (PQ[J] * PF[J] - Y0 * YP[J]) - S[J] * YF * XP[J]) /
                     (PQ[J] + Y0Q);
        };
        K[I] = Y * K[I] + exp(-XQ);
      };
    };
  };

} /* end of HUMLIK_Faddeeva_Wells */

// =================================================================================

void CEF(  // Input
    const Complex zeta,
    const Numeric Gamma_D,
    //oUTPUT
    Numeric& Re_CEF,
    Numeric& Im_CEF) {
  // The special line shape, Complex Error Function.
  // see P. W. Rosenkranz Chapter 2 of the Janssen book

  extern const Numeric PI;

  Complex retval;

  retval = (1.000e0 / sqrt(PI) / Gamma_D) *
           (122.60793178 * pow(zeta, Numeric(0.)) +
            214.38238869 * pow(zeta, Numeric(1.)) +
            181.92853309 * pow(zeta, Numeric(2.)) +
            93.15558046 * pow(zeta, Numeric(3.)) +
            30.18014220 * pow(zeta, Numeric(4.)) +
            5.91262621 * pow(zeta, Numeric(5.)) +
            0.56418958 * pow(zeta, Numeric(6.))) /
           (122.60793178 * pow(zeta, Numeric(0.)) +
            352.73062511 * pow(zeta, Numeric(1.)) +
            457.33447878 * pow(zeta, Numeric(2.)) +
            348.70391772 * pow(zeta, Numeric(3.)) +
            170.35400182 * pow(zeta, Numeric(4.)) +
            53.99290691 * pow(zeta, Numeric(5.)) +
            10.47985711 * pow(zeta, Numeric(6.)) +
            1.00000000 * pow(zeta, Numeric(7.)));

  Re_CEF = real(retval);
  Im_CEF = imag(retval);

  return;
}

// ########################################################################################

/*! 
  calculate Zeeman intensity and frequency shift for
  all splitted lines of one transition

   \param xi   Output: intensity factors
   \param eta  Output: frequency shifts
   \param N_r  Input: rotational quantum number, 
                      N+ transitions are positive, N- negative
   \param AN_r Input: total number of splitted lines

  \author Axel von Engeln
  \date   2004-01-07
*/
void Zeeman_o2_splitting_factors(  //Output
    Matrix& xi,
    Matrix& eta,
    //Input
    const Index N_r,
    const Index AN_r) {
  // rotational angular momentum quantum number
  Numeric N = abs(N_r);

  // check if N+ or N- transition
  if (N_r > 0) {
    // N+ transition

    // cycle through all magnetic quantum number (eq. split lines)
    for (Index M = -AN_r; M <= AN_r; M++) {
      // do the 3 polarizations

      // common divisor of sigma (sdiv) and pi (pdiv) lines
      Numeric sdiv = 4 * (N + 1) * (2 * N + 1) * (2 * N + 3);
      Numeric pdiv = N * (N + 1);

      // Sigma- intensity
      xi(0, M + AN_r) = 3 * (N - M + 1) * (N - M + 2) / sdiv;

      // Sigma- frequency shift
      eta(0, M + AN_r) = (M * (N - 1) - N) / pdiv;

      // PI intensity
      xi(1, M + AN_r) = 3 * ((N + 1) * (N + 1) - M * M) /
                        ((N + 1) * (2 * N + 1) * (2 * N + 3));

      // PI frequency shift
      eta(1, M + AN_r) = M * (N - 1) / pdiv;

      // Sigma+ intensity
      xi(2, M + AN_r) = 3 * (N + M + 1) * (N + M + 2) / sdiv;

      // Sigma+ frequency shift
      eta(2, M + AN_r) = (M * (N - 1) + N) / pdiv;
    }
  } else {
    // N- transition

    // cycle through all magnetic quantum number (eq. split lines)
    for (Index M = -AN_r; M <= AN_r; M++) {
      // do the 3 polarizations

      // common divisor of sigma (sdiv) and pi (pdiv) lines
      Numeric sdiv = 4 * N * (2 * N - 1) * (2 * N + 1);
      Numeric pdiv = N * (N + 1);

      // Sigma- intensity
      xi(0, M + AN_r) = 3 * (N + M) * (N + M - 1) / sdiv;

      // Sigma- frequency shift
      eta(0, M + AN_r) = (M * (N + 2) - (N + 1)) / pdiv;

      // PI intensity
      xi(1, M + AN_r) = 3 * (N * N - M * M) / (N * (2 * N - 1) * (2 * N + 1));

      // PI frequency shift
      eta(1, M + AN_r) = M * (N + 2) / pdiv;

      // Sigma+ intensity
      xi(2, M + AN_r) = 3 * (N - M) * (N - M - 1) / sdiv;

      // Sigma+ frequency shift
      eta(2, M + AN_r) = (M * (N + 2) + (N + 1)) / pdiv;
    }
  }

  return;
}

// ########################################################################################

//! Zeeman_o2_line_splitting calculates the extinction matrix and absortion vector

/*! 
The function calculates the extinction matrix and the absorption vector  
at a given frequency for a specified Zeeman transtition of O2. Both quantities 
are needed for the solution of the vector RTE if the Zeeman effect is significant.

The general structure of the output quantities is as follows:

              |A  B  0  C|
ext_mat_tmp = |B  A -D  0|
              |0  D  A  E|
              |C  0 -E  A|

and

                  |A |
abs_vec_tmp = 2 x |B |
                  |0 |
                  |C |

where the components A, B, C, D, and E are individually calculated in
Zeeman_o2_line_splitting.

!!! WARNING: This function is not threadsafe! It uses non-const static vars (OLE) !!!

   \retval ext_mat_tmp     Temporary 4x4 extinction matrix, [1/m]
   \retval abs_vec_tmp     Temporary 1x4 absorption vector, [1/m]

   \param B_field         Magnitude of Earth's magnetic field along LOS, [T]
   \param phi             Angle between the mag. field and LOS, [rad]  
   \param f_c             Spectral line (unsplit) center frequencies, [GHz] 
   \param S               Spectral line (unsplit) intensities, [cm� Hz]
   \param gamma_L         Spectral line (unsplit) pressure broadening, [GHz]  
   \param gamma_D         Spectral line (unsplit) Doppler broadening, [GHz]
   \param N_r             Rotational quantum number notation of the Zeeman transition
   \param f_grid_point    Frequency at which the calculation is done, [GHz]

   \author Nikolay Koulev, Thomas Kuhn
   \date 2003-11-30

*/
void Zeeman_o2_line_splitting(  // Output:
    MatrixView ext_mat_tmp,
    VectorView abs_vec_tmp,
    // Input:
    const Numeric B_field,
    const Numeric phi,
    const Numeric f_c,
    const Numeric S,
    const Numeric gamma_L,
    const Numeric gamma_D,
    const Index N_r,
    const Numeric f_grid_point)

{
  // Check of the rotational quantum number notation of a given Zeeman transition.
  // This parameter has integer values in the range [-39,-1]&[1,39]  and
  // [-33,-1]&[1,33] in the models of Liebe and Rosenkranz, respectively.
  assert((abs(N_r)) >= 1);
  assert((abs(N_r)) <= 39);

  // Propagation tensor G at a given frequency.
  // The 2x2 complex tensor is stored in ARTS in the form
  // of 4x2 real(!) matrix. The first index accounts for the
  // elements of the tensor taken rowwise, the second - either
  // for the real or the imaginary part, respecttively.
  //
  // |G_s_11 G_s_12|      |Re(G_s_11)  Im(G_s_11)|
  // |G_s_21 G_s_22|  ->  |Re(G_s_12)  Im(G_s_12)|
  //                      |Re(G_s_21)  Im(G_s_21)|
  //                      |Re(G_s_22)  Im(G_s_22)|
  Matrix G_s(4, 2, 0.);

  // Magnetic susceptibility tensor at a given frequency.
  // The 2x2 complex tensor is stored in ARTS in the form
  // of 4x2 real(!) matrix. The first index accounts for the
  // elements of the tensor taken rowwise, the second - either
  // for the real or the imaginary part, respecttively.
  //
  // |Chi_11 Chi_12|      |Re(Chi_11)  Im(Chi_11)|
  // |Chi_21 Chi_22|  ->  |Re(Chi_12)  Im(Chi_12)|
  //                      |Re(Chi_21)  Im(Chi_21)|
  //                      |Re(Chi_22)  Im(Chi_22)|
  Matrix Chi(4, 2, 0.);

  // The polarization  matrix accounts for the contributions of
  // the 3 different polarizations of the components of the Zeeman
  // split due to 3 different values of DeltaM.
  //
  // The polarization matrix is 2x2 and is differently defined for the 3
  // different values of DeltaM. For the sake of programming convenience a
  // Tensor3 P(i,j,k) is constructed through which elements
  // the polarization matrix is expressed and dealt with.
  // The first index points to the value of DeltaM, the other
  // 2 indexes are simply those of the matrix elements of the
  // polarisation matrix for this value of DeltaM.
  //
  //                                       |P(0,0,0)  P(0,0,1)|
  // DeltaM = -1 --> Polarization matrix = |                  |
  //                                       |P(0,1,0)  P(0,1,1)|
  //
  //                                       |P(1,0,0)  P(1,0,1)|
  // DeltaM =  0 --> Polarization matrix = |                  |
  //                                       |P(1,1,0)  P(1,1,1)|
  //
  //                                       |P(2,0,0)  P(2,0,1)|
  // DeltaM =  1 --> Polarization matrix = |                  |
  //                                       |P(2,1,0)  P(2,1,1)|
  //
  Tensor3 P(3, 2, 2);

  // Definitions of the elements of the Tensor3 P.
  // DeltaM==-1
  P(0, 0, 0) = 0.5 * (1.0 - cos(phi)) * (1.0 - cos(phi));
  P(0, 0, 1) = 0.5 * sin(phi) * sin(phi);
  P(0, 1, 0) = P(0, 0, 1);
  P(0, 1, 1) = 0.5 * (1.0 + cos(phi)) * (1.0 + cos(phi));
  // DeltaM==0
  P(1, 0, 0) = P(0, 0, 1);
  P(1, 0, 1) = -P(0, 0, 1);
  P(1, 1, 0) = -P(0, 0, 1);
  P(1, 1, 1) = P(0, 0, 1);
  // DeltaM==1
  P(2, 0, 0) = 0.5 * (1.0 + cos(phi)) * (1.0 + cos(phi));
  P(2, 0, 1) = P(0, 0, 1);
  P(2, 1, 0) = P(0, 0, 1);
  P(2, 1, 1) = 0.5 * (1.0 - cos(phi)) * (1.0 - cos(phi));

  // variables for complex error function / Faddeeva function
  Vector FX(1, 0.0);
  Vector FK(1, 0.0);
  Vector FL(1, 0.0);
  Numeric FY = 0.0;

  // total angular momentum N of the lower level
  Numeric AN_r = abs(N_r);

  // number of the split components for each of the 3 polarizations
  Numeric BN_r = 2 * AN_r;

  if (N_r > 0) {
    // N+ Zeeman transition (N -> N+1)
    BN_r++;
  } else {
    // N- Zeeman transition (N -> N-1)
    AN_r--;
    BN_r--;
  };

  // calculate intensity factors and frequency shifts for all split
  // lines only once
  static Index calzeefacs = 1;
  static Matrix eta(3, (Index)BN_r), xi(3, (Index)BN_r);
  if (calzeefacs == 1) {
    Zeeman_o2_splitting_factors(xi, eta, N_r, (Index)AN_r);
    calzeefacs = 0;
  }

  // Starting the loop over the 3 possible values of the change of the
  // magnetic quantum number M upon a Zeeman transition.
  for (Index k = 0; k < 3; k++) {
    // Complex refractive index.
    Numeric Re_N_s = 0.0;
    Numeric Im_N_s = 0.0;

    // Starting the loop over 2N+1 or 2N-1 values of
    // the magnetic quantum number M for a given transition
    for (Index j = 0; j < BN_r; j++) {
      // Center frequency of the individual Zeeman components
      Numeric f_z;
      f_z = f_c + 28.026 * eta(k, j) * B_field;  // [GHz]

      // Doppler broadening for Zeeman splitting.
      Numeric ALPHA = gamma_D * f_z;  // [GHz]

      // complex error function / Faddeeva function
      // (SQRT(ln(2)) = 0.8325546e0)
      FX[0] = 0.8325546 * abs(f_grid_point - f_z) / ALPHA;
      FY = 0.8325546 * gamma_L / ALPHA;
      HUMLIK_Faddeeva_Wells(FX, FY, FK, FL);
      // real part of complex refractive index
      Re_N_s += 0.5 * S * xi(k, j) * FL[0];
      // imaginary part of complex refractive index
      Im_N_s += -0.5 * S * xi(k, j) * FK[0];
    };

    // Calculating the halved value of magnetic susceptibility tensor
    // at given frequency as a product of the polarisation matrix and
    // 	the complex refractive index.
    Chi(0, 0) += P(k, 0, 0) * Re_N_s;
    Chi(0, 1) += P(k, 0, 0) * Im_N_s;

    Chi(1, 0) += P(k, 0, 1) * Re_N_s;
    Chi(1, 1) += P(k, 0, 1) * Im_N_s;

    Chi(2, 0) += P(k, 1, 0) * Re_N_s;
    Chi(2, 1) += P(k, 1, 0) * Im_N_s;

    Chi(3, 0) += P(k, 1, 1) * Re_N_s;
    Chi(3, 1) += P(k, 1, 1) * Im_N_s;
  };

  // Calculating the propagation tensor at a given frequency
  // through the complex magnetic susceptibility tensor.
  G_s(0, 0) = -Chi(0, 1);
  G_s(0, 1) = (1.0 + Chi(0, 0));

  G_s(1, 0) = -Chi(1, 1);
  G_s(1, 1) = Chi(1, 0);

  G_s(2, 0) = -Chi(2, 1);
  G_s(2, 1) = Chi(2, 0);

  G_s(3, 0) = -Chi(3, 1);
  G_s(3, 1) = (1.0 + Chi(3, 0));

  // calculate elements of extinction/absorption
  Numeric AA = G_s(0, 0) + G_s(3, 0);
  Numeric BB = G_s(1, 0) + G_s(2, 0);
  Numeric CC = G_s(0, 0) - G_s(3, 0);
  Numeric DD = G_s(0, 1) - G_s(3, 1);
  Numeric EE = G_s(1, 1) + G_s(1, 1);

  // Non-zero absorption vector components at a given frequency.
  // The zero components are already as 0 initialized.
  abs_vec_tmp[0] += AA;
  abs_vec_tmp[1] += BB;
  abs_vec_tmp[3] += CC;

  // First row of the non-zero extinction matrix components at a given frequency.
  // The zero components are already as 0 initialized.
  ext_mat_tmp(0, 0) += AA;
  ext_mat_tmp(0, 1) += BB;
  ext_mat_tmp(0, 3) += CC;

  // Second row
  ext_mat_tmp(1, 0) += BB;
  ext_mat_tmp(1, 1) += AA;
  ext_mat_tmp(1, 2) += -DD;

  // Third row
  ext_mat_tmp(2, 1) += DD;
  ext_mat_tmp(2, 2) += AA;
  ext_mat_tmp(2, 3) += EE;

  // Fourth row
  ext_mat_tmp(3, 0) += CC;
  ext_mat_tmp(3, 2) += -EE;
  ext_mat_tmp(3, 3) += AA;

  return;
}

// #################################################################################

// This is an internal function which calculates the O2 absorption in the
// pressure broadening (Van Vleck Weisskopf line shape). For the 60 GHz line complex
// the line mixing is included.
void PWRO2VVWMixing(      //output
    Matrix& ext_mat_tmp,  // extincton matrix [1/m]
    Vector& abs_vec_tmp,  // absorption vector [1/m]
    // Input
    Numeric STR,    // line inensity of the unsplit line [cm� Hz]
    Numeric Y,      // Line mixing parameter [dimensionless]
    Numeric DF,     // Pressure broadening [GHz]
    Numeric FO,     // line center frequency [GHz]
    Numeric ff,     // frequency on the frequency grid [GHz]
    Numeric vmro2,  // volume mixing ratio of O2 [1]
    Numeric p,      // atmospheric pressure [hPa]
    Numeric TH)     // dimensionless temperature 300K/T  [1]
{
  extern const Numeric PI;

  // call VVW function for line mixing but without Zeeman splitting.
  // This is sufficient for the troposhere where ther individual lines
  // are too wide that line mixing and pressure broadening dominate over
  // Zeeman splitting.
  Numeric SF1 = (DF + (ff - FO) * Y) / ((ff - FO) * (ff - FO) + DF * DF);
  Numeric SF2 = (DF - (ff + FO) * Y) / ((ff + FO) * (ff + FO) + DF * DF);

  // sum the line absorption part for a specific spectral line on the frequency grid
  Numeric SUM = STR * (SF1 + SF2) * (ff / FO) * (ff / FO);

  // O2 absorption [1/m] -----------------------------------------------------
  // Rosenkranz uses the factor 0.5034e12 in the calculation of the abs coeff.
  // This factor is the product of several terms:
  // 0.5034e12 = ISORATIO *   VMR   * (Hz/GHz) * (k_B*300K)^-1
  //           = 0.995262 * 0.20946 *   10^-9  * 2.414322e21(hPa*cm^2*km)^-1
  //             |---- 0.2085 ----|   |---- 2.414322e12(hPa*cm^2*km)^-1 ---|
  //             |---- 0.2085 ----|   |---- 2.414322e10( Pa*cm^2*km)^-1 ---|
  // O2ABS = 2.4143e12 * SUM * PRESDA * pow(TH, 3.0) / PI;
  // O2ABS = CONT + (2.414322e10 * SUM * p[i] * pow(TH, 3.0) / PI);
  // unit conversion x Nepers/km = y 1/m  --->  y = x * 1.000e-3
  // therefore 2.414322e10 --> 2.414322e7
  // absorption coefficient [1/m]
  Numeric abs_line = 2.414322e7 / PI * vmro2 * p * pow(TH, 3.0) * SUM;

  // fill I/O extinction matrix and absorption vector properly
  abs_vec_tmp[0] += abs_line;
  for (Index k = 0; k < 4; ++k) ext_mat_tmp(k, k) += abs_vec_tmp[0];

  return;
}

// #################################################################################

// This is an internal function which calculates the O2 absorption in the
// pressure broadening (Van Vleck Weisskopf line shape) as well as in the Doppler
// broadening regime. The line shape function is a Voigt function.
// For the 60 GHz line complex the line mixing is included.
void PWRO2VoigtMixing(    //output
    Matrix& ext_mat_tmp,  // extincton matrix [1/m]
    Vector& abs_vec_tmp,  // absorption vector [1/m]
    // Input
    const Numeric STR,      // line inensity of the unsplit line [cm� Hz]
    const Numeric Y,        // Line mixing parameter [dimensionless]
    const Numeric DF,       // Pressure broadening [GHz]
    const Numeric gamma_D,  // Doppler broadening of spectral line [GHz]
    const Numeric FO,       // line center frequency [GHz]
    const Numeric ff,       // frequency on the frequency grid [GHz]
    const Numeric vmro2,    // volume mixing ratio of O2 [1]
    const Numeric p,        // atmospheric pressure [hPa]
    const Numeric TH)       // dimensionless temperature 300K/T  [1]
{
  extern const Numeric PI;
  const Numeric lnpi = 0.4697186e0;

  // variables for complex error function / Faddeeva function
  Vector FX(1, 0.0);
  Numeric FY = 0.0;
  Vector FK(1, 0.0);
  Vector FL(1, 0.0);

  // complex error function / Faddeeva function
  // (SQRT(ln(2)) = 0.8325546e0)
  FX[0] = 0.8325546 * (ff - FO) / gamma_D;
  FY = 0.8325546 * DF / gamma_D;
  HUMLIK_Faddeeva_Wells(FX, FY, FK, FL);

  Numeric SF1 = (lnpi / gamma_D * FK[0]) + FL[0] * Y;

  // complex error function / Faddeeva function
  // (SQRT(ln(2)) = 0.8325546e0)
  FX[0] = 0.8325546 * (ff + FO) / gamma_D;
  FY = 0.8325546 * DF / gamma_D;
  HUMLIK_Faddeeva_Wells(FX, FY, FK, FL);

  Numeric SF2 = (lnpi / gamma_D * FK[0]) + FL[0] * Y;

  // sum the line absorption part for a specific spectral line on the frequency grid
  Numeric SUM = STR * (SF1 + SF2) * (ff / FO);

  // O2 absorption [1/m] -----------------------------------------------------
  // Rosenkranz uses the factor 0.5034e12 in the calculation of the abs coeff.
  // This factor is the product of several terms:
  // 0.5034e12 = ISORATIO *   VMR   * (Hz/GHz) * (k_B*300K)^-1
  //           = 0.995262 * 0.20946 *   10^-9  * 2.414322e21(hPa*cm^2*km)^-1
  //             |---- 0.2085 ----|   |---- 2.414322e12(hPa*cm^2*km)^-1 ---|
  //             |---- 0.2085 ----|   |---- 2.414322e10( Pa*cm^2*km)^-1 ---|
  // O2ABS = 2.4143e12 * SUM * PRESDA * pow(TH, 3.0) / PI;
  // O2ABS = CONT + (2.414322e10 * SUM * p[i] * pow(TH, 3.0) / PI);
  // unit conversion x Nepers/km = y 1/m  --->  y = x * 1.000e-3
  // therefore 2.414322e10 --> 2.414322e7
  // absorption coefficient [1/m]
  Numeric abs_line = 2.414322e7 / PI * vmro2 * p * pow(TH, 3.0) * SUM;
  //  cout << "%%%  abs_line=" << abs_line << "    " << SUM << "\n";
  // fill I/O extinction matrix and absorption vector properly
  abs_vec_tmp[0] += abs_line;
  for (Index k = 0; k < 4; ++k) ext_mat_tmp(k, k) += abs_vec_tmp[0];

  return;
}

// #################################################################################

// number of spectral lines in the different model versions
static const Index n_lines_PWR88 = 40;  // all O2 lines (range: 50-850 GHz)
static const Index n_lines_PWR93 = 40;  // all O2 lines (range: 50-850 GHz)
static const Index n_lines_PWR98 = 40;  // all O2 lines (range: 50-850 GHz)

// rotational quantum number of the spectral lines
static const Index QM93[] = {
    0,   -1,  +1,  -3,  +3,  -5,  +5,  -7,  +7,  -9,  +9,  -11, +11, -13,
    +13, -15, +15, -17, +17, -19, +19, -21, +21, -23, +23, -25, +25, -27,
    +27, -29, +29, -31, +31, -33, +33, 0,   0,   -0,  0,   0,   0};

// rotational quantum number of the spectral lines
static const Index QM98[] = {
    0,   -1,  +1,  -3,  +3,  -5,  +5,  -7,  +7,  -9,  +9,  -11, +11, -13,
    +13, -15, +15, -17, +17, -19, +19, -21, +21, -23, +23, -25, +25, -27,
    +27, -29, +29, -31, +31, -33, +33, 0,   0,   -0,  0,   0,   0};

// line center frequencies for the model version 1993
static const Numeric F93[] = {
    0,        118.7503, 56.2648,  62.4863,  58.4466,  60.3061, 59.5910,
    59.1642,  60.4348,  58.3239,  61.1506,  57.6125,  61.8002, 56.9682,
    62.4112,  56.3634,  62.9980,  55.7838,  63.5685,  55.2214, 64.1278,
    54.6712,  64.6789,  54.1300,  65.2241,  53.5957,  65.7648, 53.0669,
    66.3021,  52.5424,  66.8368,  52.0214,  67.3696,  51.5034, 67.9009,
    368.4984, 424.7631, 487.2494, 715.3932, 773.8397, 834.1453};

// line center frequencies for the model version 1998
static const Numeric F98[] = {
    0,        118.7503, 56.2648,  62.4863,  58.4466,  60.3061, 59.5910,
    59.1642,  60.4348,  58.3239,  61.1506,  57.6125,  61.8002, 56.9682,
    62.4112,  56.3634,  62.9980,  55.7838,  63.5685,  55.2214, 64.1278,
    54.6712,  64.6789,  54.1300,  65.2241,  53.5957,  65.7648, 53.0669,
    66.3021,  52.5424,  66.8368,  52.0214,  67.3696,  51.5034, 67.9009,
    368.4984, 424.7632, 487.2494, 715.3931, 773.8397, 834.1458};

// line strength (taken from HITRAN92) for the model version 1993 at T=300K [cm� * Hz]
static const Numeric S93[] = {
    0,          0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14, 0.3351E-14,
    0.3292E-14, 0.3721E-14, 0.3891E-14, 0.3640E-14, 0.4005E-14, 0.3227E-14,
    0.3715E-14, 0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14, 0.1391E-14,
    0.1808E-14, 0.9124E-15, 0.1230E-14, 0.5603E-15, 0.7842E-15, 0.3228E-15,
    0.4689E-15, 0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15, 0.4264E-16,
    0.6899E-16, 0.1924E-16, 0.3229E-16, 0.8191E-17, 0.1423E-16, 0.6460E-15,
    0.7047E-14, 0.3011E-14, 0.1826E-14, 0.1152E-13, 0.3971E-14};

// line strength (intensities in the submm range are updated according to HITRAN96)
// for the model version 1998 at T=300K [cm� * Hz]
static const Numeric S98[] = {
    0,          0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14, 0.3351E-14,
    0.3292E-14, 0.3721E-14, 0.3891E-14, 0.3640E-14, 0.4005E-14, 0.3227E-14,
    0.3715E-14, 0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14, 0.1391E-14,
    0.1808E-14, 0.9124E-15, 0.1230E-14, 0.5603E-15, 0.7842E-15, 0.3228E-15,
    0.4689E-15, 0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15, 0.4264E-16,
    0.6899E-16, 0.1924E-16, 0.3229E-16, 0.8191E-17, 0.1423E-16, 0.6494E-15,
    0.7083E-14, 0.3025E-14, 0.1835E-14, 0.1158E-13, 0.3993E-14};

// line mixing y parameter for the calculation of Y [1/bar]
static const Numeric Y93[] = {
    0,       -0.0233, 0.2408,  -0.3486, 0.5227,  -0.5430, 0.5877,
    -0.3970, 0.3237,  -0.1348, 0.0311,  0.0725,  -0.1663, 0.2832,
    -0.3629, 0.3970,  -0.4599, 0.4695,  -0.5199, 0.5187,  -0.5597,
    0.5903,  -0.6246, 0.6656,  -0.6942, 0.7086,  -0.7325, 0.7348,
    -0.7546, 0.7702,  -0.7864, 0.8083,  -0.8210, 0.8439,  -0.8529,
    0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000};

// y parameter for the calculation of Y [1/bar].
// These values are from P. W. Rosenkranz, Interference coefficients for the
// overlapping oxygen lines in air, JQSRT, 1988, Volume 39, 287-297.
static const Numeric Y88[] = {
    0,       -0.0244, 0.2772,  -0.4068, 0.6270,  -0.6183, 0.6766,
    -0.4119, 0.3290,  0.0317,  -0.1591, 0.1145,  -0.2068, 0.3398,
    -0.4158, 0.3922,  -0.4482, 0.4011,  -0.4442, 0.4339,  -0.4687,
    0.4783,  -0.5074, 0.5157,  -0.5403, 0.5400,  -0.5610, 0.5719,
    -0.5896, 0.6046,  -0.6194, 0.6347,  -0.6468, 0.6627,  -0.6718,
    0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000};

// temperature exponent of the line strength in [1]
/*static const Numeric BE[] = {	0,   0.009,   0.015,   0.083,   0.084, 
                                  0.212,	 0.212,	  0.391,   0.391, 
                                  0.626,	 0.626,	  0.915,   0.915, 
                                  1.260,	 1.260,	  1.660,   1.665,   
                                  2.119,	 2.115,	  2.624,   2.625, 
                                  3.194,	 3.194,	  3.814,   3.814, 
                                  4.484,	 4.484,	  5.224,   5.224, 
                                  6.004,	 6.004,	  6.844,   6.844, 
                                  7.744,	 7.744,	  0.048,   0.044, 
                                  0.049,	 0.145,	  0.141,   0.145};*/

// line width parameter [GHz/bar]
static const Numeric W300[] = {
    0,     1.630, 1.646, 1.468, 1.449, 1.382, 1.360, 1.319, 1.297, 1.266, 1.248,
    1.221, 1.207, 1.181, 1.171, 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 1.050,
    1.050, 1.020, 1.020, 1.000, 1.000, 0.970, 0.970, 0.940, 0.940, 0.920, 0.920,
    0.890, 0.890, 1.920, 1.920, 1.920, 1.810, 1.810, 1.810};

// v parameter for the calculation of Y [1/bar]
static const Numeric V[] = {0,       0.0079, -0.0978, 0.0844, -0.1273, 0.0699,
                            -0.0776, 0.2309, -0.2825, 0.0436, -0.0584, 0.6056,
                            -0.6619, 0.6451, -0.6759, 0.6547, -0.6675, 0.6135,
                            -0.6139, 0.2952, -0.2895, 0.2654, -0.2590, 0.3750,
                            -0.3680, 0.5085, -0.5002, 0.6206, -0.6091, 0.6526,
                            -0.6393, 0.6640, -0.6475, 0.6729, -0.6545, 0.0000,
                            0.0000,  0.0000, 0.0000,  0.0000, 0.0000};

//!   P. W. Rosenkranz oxygen absorption model
/*!
  
  - REFERENCES FOR EQUATIONS AND COEFFICIENTS:
    P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
    BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
    H.J. Liebe et al, JQSRT vol.48, pp.629-643, 1992.
    M.J. Schwartz, Ph.D. thesis, M.I.T., 1997.
  - SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
  - This version differs from Liebe's MPM92 in two significant respects:
    1. It uses the modification of the 1- line width temperature dependence
       recommended by Schwartz: (1/T).
    2. It uses the same temperature dependence (X) for submillimeter 
       line widths as in the 60 GHz band: (1/T)**0.8 
  
  history of the Rosenkranz absorption model:
  05-01-95  P. Rosenkranz: first version
  11-05-97  P. Rosenkranz: 1- line modification.
  12-16-98  P. Rosenkranz: updated submm freq's and intensities from HITRAN96
  
  \retval   ext_mat_tmp               Tensor3 of the Extinction Matrix [1/m] 
  \retval   abs_vec_tmp               Matrix of the Absorption Vector  [1/m]
  \param    geomag_strength           mag. field strength             [Tesla]
  \param    geomag_angle              mag. field orientation angle    [radians]
  \param    zeeman_o2_onoff           Zeeman splitting on or off
  \param    zeeman_o2_pressure_limit  Zeeman pressure limit           [Pa]
  \param    f_grid_point              frequency vector                [Hz]
  \param    p                         predefined pressure             [Pa]
  \param    t                         predefined temperature grid     [K] 
  \param    vmro2                     O2 vmr                          [1]
  \param    vmrh2o                    H2O vmr                         [1]
  \param    abs_model                 model version of the P. W. Rosenkranz oxygen absorption model
  \param    abs_user_parameters       scaling factor(s)

  \note     Except for  model 'user' the user defined input parameters CCin, CLin, CWin, and COin 
           are neglected (model dominates over parameters).<br>
           Allowed models:<br> 
           'Rosenkranz', 'RosenkranzLines', 'RosenkranzContinuum',
           'RosenkranzNoCoupling', and 'user'. <br> 
           For the parameter  version the following three string values are allowed:
           'PWR88', 'PWR93', 'PWR98'.<br> 
           See the user guide for detailed explanations.

           \remark   References:<br>  
           <ol>
           <li>
           P. W. Rosenkranz, Chapter 2: <i>Absorption of Microwaves by 
           Atmospheric Gases</i>, in M. A. Janssen (editor), 
           <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
           John Wiley \& Sons, Inc., 1993.
           </li>
           <li>
           P. W. Rosenkranz, <i>Interference coefficients for the 
           overlapping oxygen lines in air</i>, JQSRT, vol.39, pp.287-297, 1988.
           </li>
              <li>
               H.J. Liebe, P. W. Rosenkranz, G. A. Hufford, 
               <i>Atmospheric 60GHz Oxygen Spectrum: New Laboratory Measurements 
               and Line Parameters</i>, JQSRT, vol.48, pp.629-643, 1992.
              </li>
              <li>
               M.J. Schwartz, <i>Observation and Modeling of Atmospheric Oxygen 
               Millimeter-wave Transmittance</i>, Ph.D. Thesis, M.I.T. (1997).
              </li>
             </ol>


  \author Thomas Kuhn
  \date 2003-11-28
 */

Index absPWRO2Model(         // WS Output:
    Tensor3& ext_ten_zee,    // Tensor3 of the Extinction Matrix [1/m].
    MatrixView abs_mat_zee,  // Matrix of the Absorption Vector  [1/m].
    // WS Input:
    const Numeric geomag_strength,  // mag. field strength          [Tesla]
    const Numeric geomag_angle,     // mag. field orientation angle [radians]
    const Index zeeman_o2_onoff,    // Zeeman splitting on or off
    const Numeric zeeman_o2_pressure_limit,  // Zeeman pressure limit   [Pa]
    const Index zeeman_o2_line,         // Zeeman splitting line to calcualte
    ConstVectorView f_grid,             // frequency vector             [Hz]
    const Numeric p,                    // pressure                     [Pa]
    const Numeric t,                    // temperature                  [K]
    const Numeric vmro2,                // O2 vmr                       [1]
    const Numeric vmrh2o,               // H2O vmr                      [1]
    const String& abs_model,            // model selection string
    const Vector& abs_user_parameters)  // scaling factor(s)
{
  extern const Numeric PI;
  const Numeric Hz_to_GHz = 1.000000e-9;  // [GHz/Hz]
  const Numeric Pa_to_hPa = 1.000000e-2;  // [Pa/hPa]

  // convert p to hPa and use phPa within this routine
  Numeric phPa = p * Pa_to_hPa;

  // return value:  1=true, 0=false
  Index retval = 1;

  // isotopic ratio for O^16O^16
  //  const Numeric ISORATIO = 0.995262;

  // default scaling factors (=unity) for the original Rosenkranz models
  Numeric CC = 1.0;  // scaling factor of the O2 continuum
  Numeric CL = 1.0;  // scaling factor of the O2 line intensity
  Numeric CW = 1.0;  // scaling factor of the O2 line width parameter
  Numeric CO = 1.0;  // scaling factor of the O2 line mixing parameter

  // widths in MHz/mbar for the O2 continuum
  const Numeric WB300 = 0.56;  // [MHz/mbar]=[MHz/hPa]
  const Numeric X = 0.80;      // [1]
  //  const Numeric CONTFAC = 1.2300e-10; //

  // select the parameter set (!!model dominates values!!):
  Index n_lines = 0;
  int flag98 = 1;
  const Numeric* F;
  const Numeric* S300;
  const Numeric* Y300;
  const Index* QM;

  if (abs_model == "PWR88")  // ---------- model version 1988
  {
    // number of spectral lines in the model
    n_lines = n_lines_PWR88;
    // fill in the appropriate version parameters
    F = F93;     // line center frequencies from 1993 version
    S300 = S93;  // line intensity from 1993 version
    Y300 = Y88;  // line mixing from 1988 version
    QM = QM93;   // line quantum number
  } else if (abs_model == "PWR88scaling") {
    // check if the input parameter vector has appropriate size
    if (abs_user_parameters.nelem() != 4)
      throw runtime_error(
          "Method absPWRO2Model: Vector abs_user_parameters must have 4 element.");
    // reset the scaling factors
    CC = abs_user_parameters[0];
    CL = abs_user_parameters[1];
    CW = abs_user_parameters[2];
    CO = abs_user_parameters[3];
    // number of spectral lines in the model
    n_lines = n_lines_PWR88;
    // fill in the appropriate version parameters
    F = F93;                        // line center frequencies from 1993 version
    S300 = S93;                     // line intensity from 1993 version
    Y300 = Y88;                     // line mixing from 1988 version
    QM = QM93;                      // line quantum number
  } else if (abs_model == "PWR93")  // ---------- model version 1993
  {
    // number of spectral lines in the model
    n_lines = n_lines_PWR93;
    // fill in the appropriate version parameters
    F = F93;
    S300 = S93;
    Y300 = Y93;
    QM = QM93;  // line quantum number
  } else if (abs_model == "PWR93Lines") {
    // reset the scaling factors
    CC = 0.0;
    // number of spectral lines in the model
    n_lines = n_lines_PWR93;
    // fill in the appropriate version parameters
    F = F93;
    S300 = S93;
    Y300 = Y93;
    QM = QM93;  // line quantum number
  } else if (abs_model == "PWR93Continuum") {
    // reset the scaling factors
    CL = 0.0;
    CW = 0.0;
    CO = 0.0;
    // number of spectral lines in the model
    n_lines = n_lines_PWR93;
    // fill in the appropriate version parameters
    F = F93;
    S300 = S93;
    Y300 = Y93;
    QM = QM93;  // line quantum number
  } else if (abs_model == "PWR93NoCoupling") {
    // reset the scaling factors
    CO = 0.0;
    // number of spectral lines in the model
    n_lines = n_lines_PWR93;
    // fill in the appropriate version parameters
    F = F93;
    S300 = S93;
    Y300 = Y93;
    QM = QM93;  // line quantum number
  } else if (abs_model == "PWR93scaling") {
    // check if the input parameter vector has appropriate size
    if (abs_user_parameters.nelem() != 4)
      throw runtime_error(
          "Method absPWRO2Model: Vector abs_user_parameters must have 4 element.");
    // reset the scaling factors
    CC = abs_user_parameters[0];
    CL = abs_user_parameters[1];
    CW = abs_user_parameters[2];
    CO = abs_user_parameters[3];
    // number of spectral lines in the model
    n_lines = n_lines_PWR93;
    // fill in the appropriate version parameters
    F = F93;
    S300 = S93;
    Y300 = Y93;
    QM = QM93;                      // line quantum number
  } else if (abs_model == "PWR98")  // ---------- model version 1998
  {
    // set flag for 1998 version identification
    flag98 = 0;
    // number of spectral lines in the model
    n_lines = n_lines_PWR98;
    // fill in the appropriate version parameters
    F = F98;
    S300 = S98;
    Y300 = Y93;
    QM = QM98;  // line quantum number
  } else if (abs_model == "PWR98Lines") {
    // set flag for 1998 version identification
    flag98 = 0;
    // reset the scaling factors
    CC = 0.000;
    // number of spectral lines in the model
    n_lines = n_lines_PWR98;
    // fill in the appropriate version parameters
    F = F98;
    S300 = S98;
    Y300 = Y93;
    QM = QM98;  // line quantum number
  } else if (abs_model == "PWR98Continuum") {
    // set flag for 1998 version identification
    flag98 = 0;
    // reset the scaling factors
    CL = 0.0;
    CW = 0.0;
    CO = 0.0;
    // number of spectral lines in the model
    n_lines = n_lines_PWR98;
    // fill in the appropriate version parameters
    F = F98;
    S300 = S98;
    Y300 = Y93;
    QM = QM98;  // line quantum number
  } else if (abs_model == "PWR98NoCoupling") {
    // set flag for 1998 version identification
    flag98 = 0;
    // reset the scaling factors
    CO = 0.0;
    // number of spectral lines in the model
    n_lines = n_lines_PWR98;
    // fill in the appropriate version parameters
    F = F98;
    S300 = S98;
    Y300 = Y93;
    QM = QM98;  // line quantum number
  } else if (abs_model == "PWR98scaling") {
    // set flag for 1998 version identification
    flag98 = 0;
    // check if the input parameter vector has appropriate size
    if (abs_user_parameters.nelem() != 4)
      throw runtime_error(
          "Method absPWRO2Model: Vector abs_user_parameters must have 4 element.");
    // set the scaling factors
    CC = abs_user_parameters[0];
    CL = abs_user_parameters[1];
    CW = abs_user_parameters[2];
    CO = abs_user_parameters[3];
    // number of spectral lines in the model
    n_lines = n_lines_PWR98;
    // fill in the appropriate version parameters
    F = F98;
    S300 = S98;
    Y300 = Y93;
    QM = QM98;  // line quantum number
  } else {
    retval = 0;
    ostringstream os;
    os << "Method absPWRO2Model: wrong model values given.\n"
       << "Valid models are: 'Rosenkranz', 'RosenkranzLines', RosenkranzContinuum, "
       << "'RosenkranzNoCoupling', and 'user'" << '\n';
    throw runtime_error(os.str());
  }

  // range of lines to take into account for the line absorption part
  //  const Index first_line = 0;         // first line for calculation
  //  const Index last_line  = n_lines-1; // last line for calculation

  // const Index first_line = 3;         // first line for calculation
  // const Index last_line  = 3;        // last line for calculation

  // relative inverse temperature [1]
  Numeric TH = 300.0 / t;
  Numeric TH1 = (TH - 1.0);
  Numeric B = pow(TH, X);
  // partial pressure of H2O and dry air [hPa]
  Numeric PRESWV = phPa * vmrh2o;
  Numeric PRESDA = phPa - PRESWV;
  Numeric ACONST = .5034E12 * PRESDA * TH * TH * TH * vmro2 / 0.2085;
  Numeric DEN = 0.001 * (PRESDA * B + 1.1 * PRESWV * TH);  // [hPa]
  Numeric DENS = 0.001 * (PRESDA + 1.1 * PRESWV) * TH;     // [hPa]
  Numeric DFNR = WB300 * DEN;                              // [GHz]

  // --- Doppler broadening constant ----------------------------
  Numeric gamma_D = .760E-7 * sqrt(t);

  // continuum absorption [1/m/GHz]
  //  Numeric CCONT  = CC * CONTFAC * TH * TH * phPa;

  // add absorption into temporary array, necessary since some
  // elements might be set to zero later on if found negative
  Matrix abs_mat_tmp(f_grid.nelem(), 4, 0.0);
  Tensor3 ext_ten_tmp(f_grid.nelem(), 4, 4, 0.0);

  // loop over all spectral line center frequencies
  Numeric BFAC = 0.0;

  for (Index K = 1; K < n_lines + 1; K++) {
    if (K > 34) {
      Numeric KN = K + K - 35;
      BFAC = exp(-6.89526E-3 * KN * (KN + 1) * TH1);
    } else {
      Numeric K2 = K / 2;
      if (K2 * 2 != K) {
        BFAC = exp(-6.89526E-3 * K * (K + 1) * TH1);
      }
    }

    Numeric DF = CW * W300[K];
    // change in 118GHz line in 1998 version compared to older versions
    if ((flag98 == 0) && (abs((F[K] - 118.75)) < 0.10)) {
      // 118 line update in '98 version according to
      // Ph.D. Thesis of M. J. Schwartz, MIT, 1997
      DF *= DENS;  // [GHz]
    } else {
      // versions 1993 and 1988 have a simpler temperature dependence
      DF *= DEN;  // [GHz]
    };

    Numeric Y = .001 * phPa * B * (Y300[K] + V[K] * TH1);

    Numeric STR = S300[K] * BFAC * ACONST;

    // Zeeman effect on/off switch, only 1 line is calculated
    // with Zeeman! Shitty programming.....
    if ((zeeman_o2_onoff == 1) && (p < zeeman_o2_pressure_limit) &&
        (zeeman_o2_line == QM[K])) {
      for (Index J = 0; J < f_grid.nelem();
           J++)  // move into routine one day...
      {
        // cout << "***absPWRO2Model***  Zeeman,  QM=" << QM[l] << "\n";
        // call Zeeman splitting function for a single spectral line
        Zeeman_o2_line_splitting(ext_ten_tmp(J, joker, joker),
                                 abs_mat_tmp(J, joker),
                                 geomag_strength,
                                 geomag_angle,
                                 F[K],
                                 STR,
                                 DF,
                                 gamma_D,
                                 QM[K],
                                 f_grid[J] * Hz_to_GHz);
      }

    } else {
      // Doppler broadening
      Numeric ALPHA = gamma_D * F[K];

      /* RCD = ratio of collision to doppler broadening */
      Numeric RCD = DF / ALPHA;

      for (Index J = 0; J < f_grid.nelem(); J++) {
        Numeric ff = f_grid[J] * Hz_to_GHz;

        Numeric FD = ff - F[K];

        if ((RCD - 40) <= 0.0) {
          /* FD = F[K] - freq[J]; */
          Vector XD(1);
          Vector FK(1, 0.0), FL(1, 0.0);
          XD[0] = 0.8325546 * abs(-FD / ALPHA);
          /* call voigt function */
          HUMLIK_Faddeeva_Wells(XD, 0.8325546 * RCD, FK, FL);

          Numeric add_abs = STR * FK[0] * (ff / F[K]) * (ff / F[K]);

          abs_mat_tmp(J, 0) += add_abs;
          for (Index N = 0; N < 4; N++) ext_ten_tmp(J, N, N) += add_abs;
        } else {
          /* FD = freq[J] - F[K]; */
          Numeric FS = ff + F[K];
          Numeric SF1 = (DF + FD * Y) / (FD * FD + DF * DF);
          Numeric SF2 = (DF - FS * Y) / (FS * FS + DF * DF);

          Numeric add_abs = STR * (SF1 + SF2) * (ff / F[K]) * (ff / F[K]) / PI;

          abs_mat_tmp(J, 0) += add_abs;
          for (Index N = 0; N < 4; N++) ext_ten_tmp(J, N, N) += add_abs;
          /*	      printf("%e\n",ab_tmp[J]);*/
        }
      }
    }
  }

  // Check for negative line spectrum, add continuum, add to total:
  for (Index J = 0; J < f_grid.nelem(); J++) {
    if (abs_mat_tmp(J, 0) < 0.0) {
      abs_mat_tmp(J, joker) = 0;
      ext_ten_tmp(J, joker, joker) = 0.0;
    }
    Numeric ff = f_grid[J] * Hz_to_GHz;
    abs_mat_tmp(J, 0) +=
        1.6E-17 * ff * ff * DFNR / (TH * (ff * ff + DFNR * DFNR)) * ACONST / PI;
    abs_mat_zee(J, joker) += abs_mat_tmp(J, joker);
    ext_ten_zee(J, joker, joker) += ext_ten_tmp(J, joker, joker);
  }

  return retval;  // 1=true, 0=false
}

// #################################################################################

//!   oxygen absorption models
/*
  See arts -d absO2Model for detailed documentation.
  
  \retval   abs_coef            absorption/extinction of oxygen  [1/m]

  \param    f_grid         predefined frequency grid        [Hz]
  \param    abs_p          predefined pressure              [Pa]
  \param    abs_t          predefined temperature grid      [K] 
  \param    abs_vmr        H2O volume mixing ratio profile  [1]
  \param    model          model version of the P. W. Rosenkranz oxygen absorption model
  \param    abs_user_parameters  allows user defined input parameter set 
                                 (CCin, CLin, CWin, and COin)<br> or choice of 
                                 pre-defined parameters of specific models (see note below).

  \note     Except for  model 'user' the input parameter vector 'abs_user_parameters ' 
            must have exactly one or four emements which contains the overall scaling factor or the 
            scaling factors for the O2 continuum, O2 line intensities, O2 line broadening, 
            and O2 line mixing.<br>
            Allowed models are<br> 
            'PWR88'           : Rosenkranz 1993 model with line mixing parameters from 
            the JQSRT 1988 paper<br> 
            'PWR88scaling'    : possibility to scale the whole 'PWR88' model by a constant factor<br> 
            'PWR93'           : full Rosenkranz 1993 version as described in the Janssen book<br> 
            'PWR93Lines'      : only the line absorption of 'PWR93'<br> 
            'PWR93Continuum'  : only the continuum absorption of 'PWR93'<br> 
            'PWR93NoCoupling' : full 'PWR93' but without line mixing<br> 
            'PWR93scaling'    : possibility to scale the whole 'PWR93' model by a constant factor<br> 
            'PWR98'           : full modified Rosenkranz 1993 version 1998 as described in the F77 code<br> 
            'PWR98Lines'      : only the line absorption of 'PWR98'<br> 
            'PWR98Continuum'  : only the continuum absorption of 'PWR98'<br>
            'PWR98NoCoupling' : full 'PWR98' but without line mixing<br>
            'PWR98scaling'    : possibility to scale the whole 'PWR98' model by a constant factor<br>
            <br> 
            See the user guide for detailed explanations.

   \remark   dieeferent pre-defined absorption models can be selected by using
             the input variable 'abs_model'.<br>
             <br>  
             References:<br>  
             <ol>
              <li>
              P. W. Rosenkranz, Chapter 2: <i>Absorption of Microwaves by 
               Atmospheric Gases</i>, in M. A. Janssen (editor), 
               <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
               John Wiley \& Sons, Inc., 1993.
               </li>
               <li>
               P. W. Rosenkranz, <i>Interference coefficients for the 
               overlapping oxygen lines in air</i>, JQSRT, vol.39, pp.287-297, 1988.
              </li>
              <li>
               H.J. Liebe, P. W. Rosenkranz, G. A. Hufford, 
               <i>Atmospheric 60GHz Oxygen Spectrum: New Laboratory Measurements 
               and Line Parameters</i>, JQSRT, vol.48, pp.629-643, 1992.
              </li>
              <li>
               M.J. Schwartz, <i>Observation and Modeling of Atmospheric Oxygen 
               Millimeter-wave Transmittance</i>, Ph.D. Thesis, M.I.T. (1997).
              </li>
             </ol>


  \author Thomas Kuhn
  \date 2003-12-01
 */

// #################################################################################

//!   oxygen absorption models with possible Zeeman splitting of each spectral line
/*!
  See arts -d absO2ZeemanModel for detailed documentation.
  
  \retval   ext_mat_zee               tensor3 of the Extinction Matrix [1/m]
  \retval   abs_vec_zee               matrix of the Absorption Vector  [1/m]
  
  \param    geomag_los               geomagnetic field B, strength and angle
  \param    f_grid                   predefined frequency grid        [Hz]
  \param    zeeman_o2_onoff          Zeeman effect flag on or off          [1]
  \param    zeeman_o2_pressure_limit pressure limit from which on Zeeman effect 
                                    should be considered if zeeman_o2_onoff is true
  \param    zeeman_o2_line           line to calculate with Zeeman 
  \param    ppath_index              index of propagation path index
  \param    rtp_pressure             total atmospheric pressure       [Pa]
  \param    rtp_temperature          temperature                      [K]
  \param    rtp_vmr             list of species vmr              [1]
  \param    species_index            index of key species   
  \param    abs_model                model version of the P. W. Rosenkranz oxygen 
                                    absorption model
  \param    abs_user_parameters      allows user defined input parameter set 
                                    (CCin, CLin, CWin, and COin)<br> or choice of 
                                    pre-defined parameters of specific models (see note below).
  \param    stokes_dim               dimension of Stokes vector

  \note     Except for  model 'user' the input parameter vector 'abs_user_parameters ' 
            must have exactly one or four emements which contains the overall scaling factor or the 
            scaling factors for the O2 continuum, O2 line intensities, O2 line broadening, 
            and O2 line mixing.<br>
            Allowed models are<br> 
            'PWR88'           : Rosenkranz 1993 model with line mixing parameters from 
            the JQSRT 1988 paper<br> 
            'PWR88scaling'    : possibility to scale the whole 'PWR88' model by a constant factor<br> 
            'PWR93'           : full Rosenkranz 1993 version as described in the Janssen book<br> 
            'PWR93Lines'      : only the line absorption of 'PWR93'<br> 
            'PWR93Continuum'  : only the continuum absorption of 'PWR93'<br> 
            'PWR93NoCoupling' : full 'PWR93' but without line mixing<br> 
            'PWR93scaling'    : possibility to scale the whole 'PWR93' model by a constant factor<br> 
            'PWR98'           : full modified Rosenkranz 1993 version 1998 as described in the F77 code<br> 
            'PWR98Lines'      : only the line absorption of 'PWR98'<br> 
            'PWR98Continuum'  : only the continuum absorption of 'PWR98'<br>
            'PWR98NoCoupling' : full 'PWR98' but without line mixing<br>
            'PWR98scaling'    : possibility to scale the whole 'PWR98' model by a constant factor<br>
            <br> 
            See the user guide for detailed explanations.

   \remark   dieeferent pre-defined absorption models can be selected by using
             the input variable 'abs_model'.<br>
             <br>  
             References:<br>  
             <ol>
              <li>
              P. W. Rosenkranz, Chapter 2: <i>Absorption of Microwaves by 
               Atmospheric Gases</i>, in M. A. Janssen (editor), 
               <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
               John Wiley \& Sons, Inc., 1993.
              </li>
              <li>
              P. W. Rosenkranz, <i>Interference coefficients for the 
               overlapping oxygen lines in air</i>, JQSRT, vol.39, pp.287-297, 1988.
              </li>
              <li>
               H.J. Liebe, P. W. Rosenkranz, G. A. Hufford, 
               <i>Atmospheric 60GHz Oxygen Spectrum: New Laboratory Measurements 
               and Line Parameters</i>, JQSRT, vol.48, pp.629-643, 1992.
              </li>
              <li>
               M.J. Schwartz, <i>Observation and Modeling of Atmospheric Oxygen 
               Millimeter-wave Transmittance</i>, Ph.D. Thesis, M.I.T. (1997).
              </li>
             </ol>


  \author Thomas Kuhn
  \date 2003-12-01
 */
void absO2ZeemanModel(  // WS Output:
    // define internal absorption arrays
    Tensor3& ext_mat_zee,  // Tensor3 of the Extinction Matrix [1/m].
    Matrix& abs_vec_zee,   // Matrix of the Absorption Vector  [1/m].
    // WS Input:
    const Matrix& geomag_los,  // [Magnetic field B, angle between B and los]
    const Vector& f_grid,      // frequency vector                 [Hz]
    const Index& zeeman_o2_onoff,
    const Numeric& zeeman_o2_pressure_limit,
    const Index& zeeman_o2_line,
    const Index& ppath_index _U_,       // index of propagation path index
    const Numeric& rtp_pressure,        // total atmospheric pressure       [Pa]
    const Numeric& rtp_temperature,     // temperature                      [K]
    const Vector& rtp_vmr,              // list of species vmr              [1]
    const ArrayOfIndex& species_index,  // index of key species
    const String& abs_model,            // model selection string
    const Vector& abs_user_parameters,  // scaling factor(s)
    const Index& stokes_dim             // dimension of Stokes vector
) {
  // check dimension of vector RT calculation:
  /*
  if (stokes_dim != 4)
    {
      ostringstream os;
      os << "absO2ZeemanModel: error detected. Zeeman effect can only be "
         << "calculated for stokes_dim == 4.\n" 
         << "user defined stokes_dim is " << stokes_dim << "\n";
      throw runtime_error(os.str());
    };
  */

  // set for the moment the O2 VMR to a constant value of 0.21
  //const Numeric vmro2 = 0.20946e0;
  Numeric vmro2;
  vmro2 = rtp_vmr[species_index[SPECIES_INDEX_O2]];

  Numeric vmrh2o;
  vmrh2o = rtp_vmr[species_index[SPECIES_INDEX_H2O]];

  // resize and initialization the output variables
  ext_mat_zee.resize(f_grid.nelem(), stokes_dim, stokes_dim);
  abs_vec_zee.resize(f_grid.nelem(), stokes_dim);
  ext_mat_zee = 0.0;
  abs_vec_zee = 0.0;

  // Loop all frequencies (now done inside absPWRO2model
  //  for ( Index i=0; i<f_grid.nelem(); ++i )
  {
    // return index
    Index ret;

    // call the absorption function
    ret = absPWRO2Model(ext_mat_zee,
                        abs_vec_zee,
                        // FIXME: should be zero or direct numbers here...
                        geomag_los(0, 0),  // magnetic field strength
                        geomag_los(0, 1),  // angle propagation, magnetic field
                        zeeman_o2_onoff,
                        zeeman_o2_pressure_limit,
                        zeeman_o2_line,
                        f_grid,
                        rtp_pressure,
                        rtp_temperature,
                        vmro2,
                        vmrh2o,
                        abs_model,
                        abs_user_parameters);

    // output variable abs/ext (vector/matrix) should be passed
    // to the agenda for gas absorption
    if (ret != 1) {
      ostringstream os;
      os << "absO2ZeemanModel:\n"
         << "error in O2 absorption function *absPWRO2Model* detected\n"
         << "occured error code:" << ret << "\n"
         << "error codes:  1 = no error\n"
         << "              0 = wrong model name\n";
      throw runtime_error(os.str());
    };
  };

  return;
}

// #################################################################################

void ZeemanO2Settings(  // Output
    Index& zeeman_o2_onoff,
    Numeric& zeeman_o2_pressure_limit,
    Index& zeeman_o2_line,
    // Keywords
    const Index& ZeemanO2OnOff,
    const Numeric& ZeemanO2PressureLimit,
    const Index& ZeemanO2Line) {
  // set Zeeman O2 on / off
  zeeman_o2_onoff = ZeemanO2OnOff;

  // set upper pressure level for Zeeman
  zeeman_o2_pressure_limit = ZeemanO2PressureLimit;

  // line to calculate with Zeeman splitting
  zeeman_o2_line = ZeemanO2Line;

  return;
}

// #################################################################################

//! test_zeeman
/*!
   See the the online help (arts -d FUNCTION_NAME)


   \author Thomas Kuhn, Claudia Emde, Oliver Lemke, Nikolay Koulev
   \date   2003-12-10

*/
/*
void test_zeeman(
                 const Agenda& opt_prop_gas_agenda
                )
{
  opt_prop_gas_agenda.execute();
}
*/
