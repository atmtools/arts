/* Copyright (C) 2003 Oliver Lemke <olemke@uni-bremen.de>
                                                                                
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
  \author Oliver Lemke 
  \date   2003-11-28
  
  \brief  This file has functions for calculating the propagation
          tensor in the case of Zeeman effect for oxygen.
 
*/


// =================================================================================


#include <iostream>
#include <cmath>
#include <cstdlib>
#include "math_funcs.h"
#include "arts.h"
#include "matpackIII.h"
#include "mystring.h"
#include "xml_io.h"
#include "complex.h"
#include "exceptions.h"
#include "messages.h"
#include "abs_o2_models.h"
#include "auto_md.h"



// =================================================================================


// arts defined constants
extern const Numeric EULER_NUMBER;
extern const Numeric LOG10_EULER_NUMBER;
extern const Numeric NAT_LOG_TEN;
extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT; 
extern const Numeric BOLTZMAN_CONST;
const Numeric twoPI_c=2*PI/SPEED_OF_LIGHT;
const Complex complex_i = pow(Complex (-1., 0), 0.5);
 
// conversion from neper to decibel:
const Numeric Np_to_dB  = (10.000000 * LOG10_EULER_NUMBER); // [dB/Np]
// conversion from decibel to neper:
const Numeric dB_to_Np  = (1.000000 / Np_to_dB);            // [Np/dB]
// conversion from GHz to Hz:
const Numeric GHz_to_Hz = 1.000000e9;                       // [Hz/GHz]
// conversion from Hz to GHz:
const Numeric Hz_to_GHz = 1.000000e-9;                      // [GHz/Hz]
// conversion from kPa to Pa:
const Numeric kPa_to_Pa = 1.000000e3;                       // [kPa/Pa]
// conversion from Pa to kPa:
const Numeric Pa_to_kPa = 1.000000e-3;                      // [Pa/kPa]
// conversion from hPa to Pa (hPa = mbar):
const Numeric hPa_to_Pa = 1.000000e2;                       // [hPa/Pa]
// conversion from Pa to hPa (hPa = mbar):
const Numeric Pa_to_hPa = 1.000000e-2;                      // [Pa/hPa]



// =================================================================================



//-------------------------------------------------------------------------
//! Lagrange Interpolation (internal function).
/*! 
  This function calculates the Lagrange interpolation of four interpolation 
  points as described in 
  <a href="http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html">
  Lagrange Interpolating Polynomial</a>.<br>
  The input are the four x-axis values [x0,x1,x2,x3] and their associated 
  y-axis values [y0,y1,y2,y3]. Furthermore the x-axis point "a" at which the 
  interpolation should be calculated must be given as input. NOTE that 
  the relation x2 =< x < x3 MUST hold!

  \param x     x-vector with four elements [x0,x1,x2,x3]
  \param y     y-vector with four elements: yj = y(xj), j=0,1,2,3
  \param a     interpolation point on the x-axis with x1 =< a < x2 

  \return FIXME

  \author Thomas Kuhn
  \date   2003-11-25
*/

Numeric LagrangeInterpol4( ConstVectorView x,
			   ConstVectorView y,
			   const Numeric a)
{
  // lowermost grid spacing on x-axis
  const Numeric Dlimit = 1.00000e-15;

  // Check that dimensions of x and y vector agree
  const Index n_x = x.nelem();
  const Index n_y = y.nelem();
  if ( (n_x != 4) || (n_y != 4) )
    {
      ostringstream os;
      os << "The vectors x and y must all have the same length of 4 elements!\n"
	 << "Actual lengths:\n"
         << "x:" << n_x << ", " << "y:" << n_y << ".";
      throw runtime_error(os.str());
    }

  // assure that x1 =< a < x2 holds
  if ( (a < x[1]) || (a > x[2]) )
    {
      ostringstream os;
      os << "LagrangeInterpol4: the relation x[1] =< a < x[2] is not satisfied. " 
         << "No interpolation can be calculated.\n";
      throw runtime_error(os.str());
    };

  // calculate the Lagrange polynomial coefficients for a polynomial of the order of 3
  Numeric b[4];
  for (Index i=0 ; i < 4 ; ++i)
    {
      b[i] = 1.000e0;
      for (Index k=0 ; k < 4 ; ++k)
	{
	  if ( (k != i) && (fabs(x[i]-x[k]) > Dlimit) )  
	    b[i] = b[i] * ( (a-x[k]) / (x[i]-x[k]) );
	};
    };

  Numeric ya = 0.000e0;
  for (Index i=0 ; i < n_x ; ++i) ya = ya + b[i]*y[i];

  return ya;
}


// =================================================================================


//-----------------------------------------------------------------------
//! HUMLIK_Faddeeva_Wells calculates the complex Faddeeva function (internal function)
/*! 
  This function calculates the complex Faddeeva function according to
  R. J. Wells, Rapid approximation to the Voigt/Faddeeva function and its 
  derivatives, JQSRT, vol.62, pp.29-48, 1999.
  This function must be called for each spectral line center frequency once and 
  calculates the complex Faddeeva function for the entire input frequency grid.

  CTKS  -------------------------------------------------------------------
  CTKS  
  CTKS  DESCRIPTION OF THE INPUT VECTOR X:
  CTKS  
  CTKS                         F_i - F_o
  CTKS  X_i =  SQRT(ln(2)) * --------------
  CTKS                          GAMMA_D
  CTKS       
  CTKS  DESCRIPTION OF THE INPUT PARAMETER Y:
  CTKS  
  CTKS                          GAMMA_L 
  CTKS  Y   =  SQRT(ln(2)) * --------------
  CTKS                          GAMMA_D
  CTKS  
  CTKS  DESCRIPTION:
  CTKS  F_i     :  ELEMENT OF THE INPUT FREQUENCY GRID OF CALCULATION  [Hz]
  CTKS  F_o     :  LINE CENTER FREQUENCY OF SPECTRAL LINE              [Hz]
  CTKS  GAMMA_D :  DOPPLER HALF WIDTH OF SPECTRAL LINE                 [Hz]
  CTKS  GAMMA_L :  PRESSURE BROADENING HALF WIDTH OF SPECTRAL LINE     [Hz]
  CTKS  SQRT(ln(2)) = 0.83255E0 
  CTKS  -------------------------------------------------------------------
  
  \retval K      real part of the complex Faddeeva function (equal to the Voigt function)
  \retval L      imaginary part of the complex Faddeeva function (used for line mixing)

  \param N       number of elements of vectors X, K, and L
  \param X       X[i] =  SQRT(ln(2)) * ( (f_grid[i] - F_o)) / gamma_Doppler )
                 where F_o is the line center frequency of the specific line in question
  \param Y       Y    =  SQRT(ln(2)) * (gamma_Lorentz / gamma_Doppler)

  \author Thomas Kuhn
  \date   2003-11-30
*/

void HUMLIK_Faddeeva_Wells ( // Input
		       ConstVectorView X, 
		       Numeric Y, 
                       // Output
		       Vector& K, 
		       Vector& L )
{

  /*
    To calculate the Faddeeva function with relative error less than 10^(-R).
    R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
    subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
  */

  const Index N = X.nelem();	  /* Number of frequencies in the frequency grid */

  K.resize(N); /* resize the output vectors according to the size input vector X */
  L.resize(N); /* resize the output vectors according to the size input vector X */


  const Numeric R0 = 146.7;                              /*  Region boundaries  */
  const Numeric R1 = 14.67;                              /*  for R=4            */


  /* Constants */
  const Numeric RRTPI = 0.56418958;                      /* 1/SQRT(pi)          */
  const Numeric Y0    = 1.5;                             /* for CPF12 algorithm */
  const Numeric Y0PY0 = Y0+Y0;                           /* for CPF12 algorithm */
  const Numeric Y0Q   = Y0*Y0;                           /* for CPF12 algorithm */
  const Numeric C[]   = { 1.0117281,     -0.75197147,        0.012557727,
	                  0.010022008,   -0.00024206814,     0.00000050084806 };
  const Numeric S[]   = { 1.393237,       0.23115241,       -0.15535147,
		          0.0062183662,   0.000091908299,   -0.00000062752596 };
  const Numeric T[]   = { 0.31424038,     0.94778839,        1.5976826,
		          2.2795071,      3.0206370,         3.8897249 };


  /* Local variables */
  int RG1, RG2, RG3;                                     /* y polynomial flags        */
  Numeric ABX, XQ, YQ, YRRTPI;                           /* |x|, x^2, y^2, y/SQRT(pi) */
  Numeric XLIM0, XLIM1, XLIM2, XLIM3, XLIM4;             /* |x| on region boundaries  */
  /* W4 temporary variables */
  Numeric A0=0.0, D0=0.0, D2=0.0, E0=0.0, E2=0.0, E4=0.0, H0=0.0, H2=0.0, H4=0.0, H6=0.0;
  Numeric P0=0.0, P2=0.0, P4=0.0, P6=0.0, P8=0.0, Z0=0.0, Z2=0.0, Z4=0.0, Z6=0.0, Z8=0.0;
  Numeric B1=0.0, F1=0.0, F3=0.0, F5=0.0, Q1=0.0, Q3=0.0, Q5=0.0, Q7=0.0;
  Numeric XP[5], XM[5], YP[5], YM[5];                    /* CPF12 temporary values    */
  Numeric MQ[5], PQ[5], MF[5], PF[5];
  Numeric D, YF, YPY0, YPY0Q;  


  /***** Start of executable code *****************************************/


  RG1 = 1;               /* Set flags  */
  RG2 = 1;
  RG3 = 1;
  YQ  = Y*Y;             /* y^2        */
  YRRTPI = Y*RRTPI;      /* y/SQRT(pi) */


  /* Region boundaries when both K and L are required or when R<>4 */
  XLIM0 = R0 - Y ;
  XLIM1 = R1 - Y;
  XLIM3 = 3.097*Y - 0.45;

  /* For speed the following 3 lines should replace the 3 above if R=4 and L is not required */
  /* XLIM0 = 15100.0 + Y*(40.0 + Y*3.6)                                                      */
  /* XLIM1 = 164.0 - Y*(4.3 + Y*1.8)                                                         */
  /* XLIM3 = 5.76*YQ                                                                         */

  XLIM2 = 6.8 - Y;
  XLIM4 = 18.1*Y + 1.65;
  if ( Y <= 0.000001 )                                  /* When y<10^-6       */
    {
      XLIM1 = XLIM0;                                    /* avoid W4 algorithm */
      XLIM2 = XLIM0;
    };

/*..... */

  for( int I=0; I<N; ++I)                               /* Loop over all points */
    {                                                     
      ABX = fabs( X[I] );                                /*! |x|  */
      XQ  = ABX*ABX;                                    /* ! x^2 */
      if ( ABX > XLIM0 )                                /* ! Region 0 algorithm */
	{
	  K[I] = YRRTPI / (XQ + YQ);
	  L[I] = K[I]*X[I] / Y;
	}
      else if ( ABX > XLIM1 )                           /* ! Humlicek W4 Region 1 */
	{
	  if ( RG1 != 0 )                               /* ! First point in Region 1 */
	    {
	      RG1 = 0;
	      A0 = YQ + 0.5;                            /* ! Region 1 y-dependents */
	      D0 = A0*A0;
	      D2 = YQ + YQ - 1.0;
	      B1 = YQ - 0.5;
	    };
	  D = RRTPI / (D0 + XQ*(D2 + XQ));
	  K[I] = D*Y   *(A0 + XQ);
	  L[I] = D*X[I]*(B1 + XQ);
	}
      else if ( ABX > XLIM2 )                           /* ! Humlicek W4 Region 2 */ 
	{
	  if ( RG2 != 0 )                               /* ! First point in Region 2 */
	    {
	      RG2 = 0;
	      H0 =  0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ))); /* ! Region 2 y-dependents */
	      H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0));
	      H4 = 10.5    - YQ*(6.0 - YQ*  6.0);
	      H6 = -6.0    + YQ* 4.0;
	      E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ));
	      E2 =  5.25   + YQ*(1.0  + YQ* 3.0);
	      E4 =  0.75*H6;
	      F1 = -1.875  + YQ*(5.25 + YQ*(4.5 + YQ));
	      F3 =  8.25   - YQ*(1.0  - YQ* 3.0);
	      F5 = -5.5    + YQ* 3.0;
	    };
	  D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))));
	  K[I] = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)));
	  L[I] = D*X[I]*(F1 + XQ*(F3 + XQ*(F5 + XQ)));
	}
      else if ( ABX < XLIM3 )                           /* ! Humlicek W4 Region 3 */
	{
	  if ( RG3 != 0 )                               /* ! First point in Region 3 */
	    {
	      RG3 = 0;
	      Z0 = 272.1014 + Y*(1280.829 + Y*(2802.870 + Y*(3764.966 
			    + Y*(3447.629 + Y*(2256.981 + Y*(1074.409 
			    + Y*(369.1989 + Y*(88.26741 + Y*(13.39880 
                            + Y)))))))));               /*! Region 3 y-dependents*/
	      Z2 = 211.678  + Y*(902.3066 + Y*(1758.336 + Y*(2037.310
                            + Y*(1549.675 + Y*(793.4273 + Y*(266.2987
                            + Y*(53.59518 + Y*5.0)))))));
	      Z4 = 78.86585 + Y*(308.1852 + Y*(497.3014 + Y*(479.2576
			    + Y*(269.2916 + Y*(80.39278 + Y*10.0)))));
	      Z6 = 22.03523 + Y*(55.02933 + Y*(92.75679 + Y*(53.59518
                            + Y*10.0)));
	      Z8 = 1.496460 + Y*(13.39880 + Y*5.0);
	      P0 = 153.5168 + Y*(549.3954 + Y*(919.4955 + Y*(946.8970
                            + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941
                            + Y*(4.264678 + Y*0.3183291))))))));
	      P2 = -34.16955 + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730
                             + Y*(139.4665 + Y*(56.81652 + Y*(12.79458
                             + Y*1.2733163))))));
	      P4 = 2.584042 + Y*(10.46332 + Y*(24.01655 + Y*(29.81482
                            + Y*(12.79568 + Y*1.9099744))));
	      P6 = -0.07272979 + Y*(0.9377051+ Y*(4.266322 + Y*1.273316));
	      P8 = 0.0005480304 + Y*0.3183291;
	      Q1 = 173.2355 + Y*(508.2585 + Y*(685.8378 + Y*(557.5178
                            + Y*(301.3208 + Y*(111.0528 + Y*(27.62940
                            + Y*(4.264130 + Y*0.3183291)))))));
	      Q3 = 18.97431 + Y*(100.7375 + Y*(160.4013 + Y*(130.8905
                            + Y*(55.88650 + Y*(12.79239+Y*1.273316)))));
	      Q5 = 7.985877 + Y*(19.83766 + Y*(28.88480 + Y*(12.79239
			    + Y*1.909974)));
	      Q7 = 0.6276985 + Y*(4.264130 + Y*1.273316);
	    };
	  D = 1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))));
	  K[I] = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))));
	  L[I] = D*X[I]*(Q1 + XQ*(Q3 + XQ*(Q5 + XQ*(Q7 + XQ*0.3183291))));
	}
      else                                              /* ! Humlicek CPF12 algorithm */
	{
	  YPY0 = Y + Y0;
	  YPY0Q = YPY0*YPY0;
	  K[I] = 0.000;
	  L[I] = 0.000;
	  for(int J=0; J <= 5; ++J)
	    {
	      D = X[I] - T[J];
	      MQ[J] = D*D;
	      MF[J] = 1.0 / (MQ[J] + YPY0Q);
	      XM[J] = MF[J]*D;
	      YM[J] = MF[J]*YPY0;
	      D = X[I] + T[J];
	      PQ[J] = D*D;
	      PF[J] = 1.0 / (PQ[J] + YPY0Q);
	      XP[J] = PF[J]*D;
	      YP[J] = PF[J]*YPY0;
	      L[I]  = L[I] + C[J]*(XM[J]+XP[J]) + S[J]*(YM[J]-YP[J]);
	    };
	  if ( ABX <= XLIM4 )                           /* ! Humlicek CPF12 Region I */
	    {
	      for(int J=0; J <= 5; ++J)
		{
		  K[I] = K[I] + C[J]*(YM[J]+YP[J]) - S[J]*(XM[J]-XP[J]);
		}
	    }
	  else                                          /* ! Humlicek CPF12 Region II */
	    {
	      YF = Y + Y0PY0;
	      for(int J=0; J <= 5; ++J)
		{
		  K[I] = K[I]
		    + (C[J]*(MQ[J]*MF[J]-Y0*YM[J]) + S[J]*YF*XM[J]) / (MQ[J]+Y0Q)
		    + (C[J]*(PQ[J]*PF[J]-Y0*YP[J]) - S[J]*YF*XP[J]) / (PQ[J]+Y0Q);
		};
	      K[I] = Y*K[I] + exp( -XQ );
	    };
	};
    };

};  /* end of HUMLIK_Faddeeva_Wells */



// =================================================================================



//! oxygen absorption model including Zeeman splitting (workspace method)
/*!

   See arts -d ZeemanO2 for detailed documentation.

   \retval   abs                  absorption coefficient of 
                                  H2O self continuum according to CKD_MT 1.00   [1/m]
   \param    f_grid               predefined frequency grid                     [Hz]
   \param    abs_p                predefined pressure grid                      [Pa]
   \param    abs_t                predefined temperature grid                   [K] 
   \param    abs_vmr              H2O volume mixing ratio profile               [1]
   \param    model                allows modifications of the original 
                                  model by selecting model="user". To 
				  select the pre-defined original model use
				  model="CKDMT100"
   \param    abs_user_parameters  absorption scaling factor                     [1]

   \note     This oxygen absorption model includes Zeeman splitting of the 60 GHz band.


   \remark   References:<br>  
             <ol>
              <li>
	       P. W. Rosenkranz, Chapter 2: <i>Absorption of Microwaves by 
               Atmospheric Gases</i>, in M. A. Janssen (editor), 
               <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
               John Wiley & Sons, Inc., 1993.
              </li>
              <li>
               H.J. Liebe, P. W. Rosenkranz, G. A. Hufford, 
               <i>Atmospheric 60GHz Oxygen Spectrum: New Laboratory Measurements 
               and Line Parameters</i>, JQSRT vol.48, pp.629-643 (1992).
              </li>
              <li>
               M.J. Schwartz, <i>Observation and Modeling of Atmospheric Oxygen 
               Millimeter-wave Transmittance</i>, Ph.D. Thesis, M.I.T. (1997).
              </li>
             </ol>

   \author   Nikolay Koulev, Oliver Lemke, Thomas Kuhn
   \date     2003-11-28
*/ 
/* ********************************************************************************************
void ZeemanO2 (// WS Output:
	       Tensor3& ext_mat_zee, // Tensor3 of the Extinction Matrix for the specified frequency grid in [1/m]. 
	       Matrix& abs_vec_zee,  // Matrix of the Absorption Vector for the specified frequency grid in [1/m].
	       // WS Input:
	       const Vector& LineCenterFreq   // model line center frequencies
	       const Vector& LineIntensities  // model line intensities 
               const Vector& LineBroadening   // model line pressure broadening
	       const Numeric T                // atmospheric temperature
	       const Vector& f_grid,          // frequency grid.
	       const Numeric& N_r)            // Rotational quantum number indentification of the Zeeman line.

{
  
  // molecular weight of oxygen (main isotope) [g/mol]
  const Numeric MWO2 = 31.989830e0; 

  
  // Magnitude of Earth's magnetic field along LOS in [T].
  Numeric B_field;
  
  // Dummy value for the magnetic field in [T].
  B_field = 3.5e-5;

  // Angle between the directions of Earth's magnetic field and LOS.    
  Numeric phi;

  // Dummy value for phi.
  phi = 90./57.2957914331333;
  

  // Pressure in [Pa]- just provisionary now. 
  //Later will bw replaced by p_grid.
  Numeric p;

  // Test value for the pressure in [Pa] at 95 km 'midlatitude-summer'.
  p = 6.200000e-02;

  // Temperature in [K].
  Numeric T;

  // Test value for the temperature in [K] at 95 km for 'midlatitude-summer'.
  T = 1.783000e+02; 
  //-------------------------------------------------

    


    //------------------------------------------------------------
    // Definitions of the line intensity and broadening parameters,
    // wavenumber and frequency shift.
    //--------------------------------

    // Doppler broadening line width in [Hz].
    Numeric gamma;
    
    
    // Pressure broadening line width [Hz].
    Numeric agam;
    
    
    // Wavenumber.
    Numeric k_0;

    // Intensity factor, or relative intensity, of the 
    // individual components of the Zeeman split relative 
    // to the the intensity of the unsplit line.
    Numeric xi;
    
    
    // This the Zeeman frequency shift of the individual 
    // components of the Zeeman split.
    Numeric eta;
    //--------------------------------------------------





    //---------------------------------------------------------------
    //  Definitions & resizements of the tensor and array quantities.
    //---------------------------------------------------------------
    
    // Magnetic susceptibility tensor at a given frequency. 
    // It's in the form of an array because complex matrices
    // or tensors cannot be handled by ARTS.
    Array<Complex> Chi(f_grid.nelem()*2*2,0);
    
    // Complex refractive index.
    Array<Complex> N_s(f_grid.nelem(),0);
    
    // Propagation tensor at a given frequency. It's in the form of an array 
    // because complex matrices or tensors cannot be handled by ARTS. 
    Array<Complex> G_s(f_grid.nelem()*2*2);
      
    // Resize of the tensor consisting of the Zeeman 
    // Extinction matrix and the frequency grid.
    ext_mat_zee.resize(f_grid.nelem(),4,4);
    
    // Resize of the matrix consisting of the Z
    //Z eeman absorption vector and the frequency grid.
    abs_vec_zee.resize(f_grid.nelem(),4);      

    // Polarization  matrix accounting for the contributions of 
    // the 3 different polarizations of the components of the Zeeman 
    // split due to 3 different values of DeltaM.
    Matrix P;
    P.resize(2,2);
    //---------------------------------------------------------------
    
    
    
    
    
    // Matrix variable for the content of the zeeman_intensity_coeff.xml.
    Matrix N;
    
    // Reading in content from zeeman_intensity_coeff.xml.
    xml_read_from_file ("zeeman_intensity_coeff.xml", N);
    
    
    
     if (N_r>0)
      { 
	
	// Assigning a value of the total angular momentum N of the lower  
	// level for a N+ Zeeman transition (N -> N+1).
	AN_r= abs(N_r);
	
	// Assigning a value of the number of the split components for each 
	// of the three possible polarizations of a N+ Zeeman transition.
	BN_r= 2*abs(N_r)+1;
	
	//! Resize of the frequency shift matrix 
	// for a N+ Zeeman transition.
	f_z_mat.resize(2*Index(abs(N_r))+1,3);
	
	//! Resize of the relative intensity matrix 
	// for a N+ Zeeman transition.
	xi_mat.resize(2*Index(abs(N_r))+1,3);
	
      }
    
    else if (N_r<0)
      
      {
	// Assigning a value of the total angular momentum N of the lower 
	// level for a N- Zeeman transition (N -> N-1).
	AN_r= abs(N_r)-1;
	
	// Assigning a value of the number of the split components for each
	// of the three possible polarizations of a N- Zeeman transition.
	BN_r= 2*abs(N_r)-1;
	
	
	// Resize of the frequency shift matrix 
	// for a N- Zeeman transition.
	f_z_mat.resize(2*Index(abs(N_r))-1,3);
	
	// Resize of the relative intensity matrix 
	// for a N- Zeeman transition.
	xi_mat.resize(2*Index(abs(N_r))-1,3);
      }
    
    
    
    
    
    
    
    // Starting the loop over the frequency grid.
    for (Index f_grid_index=0; f_grid_index < f_grid.nelem(); f_grid_index++)
      {
	
	//  Wavenumber claculation.
	k_0=twoPI_c*f_grid[f_grid_index];
	
	
	// Starting the loop over the 3 possible values of the change of the 
	// magnetic quantum number M upon a Zeeman transition.	  
	for (Index k=0;k<3;k++)
	  { 
	    Index DeltaM;
	    DeltaM = k-1;
	    
	    
	    
	    // Starting the loop over 2N+1 or 2N-1 values of
	    // the magnetic quantum number M for a given transition 
	    for (int j=0 ; j<BN_r;j++)
	      { 
		Numeric M;
		M = j-AN_r;
		
		
		
		// Pi polarization for a N+ Zeeman transition.
		if (N_r>0 && DeltaM==0) //
		  {
		    // Matrix representation of the relative intensity value 
		    xi_mat(j,k) = 3*((abs(N_r)+1)*(abs(N_r)+1)-M*M)
		      /((abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		    
		   		   		    
		    // Zeeman frequency shift value
		    eta = M*(abs(N_r)-1)/(abs(N_r)*(abs(N_r)+1));
		   
		
		  }
		
		
		// Pi polarization for a N- Zeeman transition.
		else if (N_r<0 && DeltaM==0)
		  
		  {
		    // Matrix representation of the relative intensity values 
		    xi_mat(j,k)= 3*(abs(N_r)*abs(N_r)-M*M)
		      /(abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		    
		    
		    // Zeeman frequency shift values
		    eta = M*(abs(N_r)+2)/(abs(N_r)*(abs(N_r)+1));
		  		    
		  }
		
		// Sigma+ polarization for a N+ Zeeman transition.
		else if (N_r>0 && DeltaM==1)
		  
		  {
		    // Matrix representation of the relative intensity
		    xi_mat(j,k) = 3*(abs(N_r)+M+1)*(abs(N_r)+M+2)/
		      (4*(abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		    
		    // Zeeman frequency shift value 
		    eta = (M*(abs(N_r)-1)+abs(N_r))/(abs(N_r)*(abs(N_r)+1));
						
		  }
		
		// Sigma- polarization for a N+ Zeeman transition.
		else if (N_r>0 && DeltaM==-1)
		  
		  {	
		    // Matrix representation of the relative intensity values
		    xi_mat(j,k) = 3*(abs(N_r)-M+1)*(abs(N_r)-M+2)/
		      (4*(abs(N_r)+1)*(2*abs(N_r)+1)*(2*abs(N_r)+3));
		   		    
		    // Zeeman frequency shift value
		    eta = (M*(abs(N_r)-1)-abs(N_r))/(abs(N_r)*(abs(N_r)+1));
						
		  }
		
		//! Sigma+ polarization for a N- Zeeman transition.
		else if (N_r<0 && DeltaM==1)
		  
		  {
		    // Matrix representation of the relative intensity values
		    xi_mat(j,k) = 3*(abs(N_r)-M)*(abs(N_r)-M-1)/
		      (4*abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		    
		    // Zeeman frequency shift value 
		    eta = (M*(abs(N_r)+2)+(abs(N_r)+1))/(abs(N_r)*(abs(N_r)+1));
		    
		  }
		
		// Sigma- polarization for a N- Zeeman transition.
		else if (N_r<0 && DeltaM==-1)
		  
		  {
		    // Matrix representation of the relative intensity values
		    xi_mat(j,k)= 3*(abs(N_r)+M)*(abs(N_r)+M-1)/
		      (4*abs(N_r)*(2*abs(N_r)-1)*(2*abs(N_r)+1));
		    
		   
		    
		    // Zeeman frequency shift value 
		    eta = (M*(abs(N_r)+2)-(abs(N_r)+1))/(abs(N_r)*(abs(N_r)+1));
		    
		  }
		
		else 
		  {
		    cout << "Error N_r or DeltaM" << endl;
		  }
		
		
		
		// Center frequency of the of the individual 
		// components of the Zeeman split.
		Numeric f_z;
		f_z = f_c + 28.03*1e+9*eta*B_field;
		f_z_mat(j,k)=f_z;
		
	
		
		// Line strength of the unsplit line in [Hz].
		Numeric S;
		Numeric theta = 300./T; 
		S = a1*1e-3*p*1e3*pow(theta,3)*exp(a2*(1-theta));
		
		// Numeric expressions for both the pressure (agam)
		// and Doppler (gamma) broadening line width 
		agam = a3*1e3* p * pow(theta,0.85);
		gamma = (1.096/sqrt(theta))*f_z*1e3;
		
		
		// Complex argument of the special line shape (CEF), 
		// used in the case of Zeeman splitting.
                Complex zeta;
		zeta = agam/gamma + complex_i*(f_grid[f_grid_index] - f_z)/gamma;
		
		
		
		// The special line shape, Complex Error Function.
		Complex CEF;
		CEF = (1/(sqrt(PI)*gamma))*(122.60793178*pow(zeta,0) 

		     + 214.38238869*pow(zeta,1) + 181.92853309*pow(zeta,2) 

		     + 93.15558046*pow(zeta,3) + 30.18014220*pow(zeta,4)

		     + 5.91262621*pow(zeta,5) + 0.56418958*pow(zeta,6))

		  /(122.60793178*pow(zeta,0) + 352.73062511*pow(zeta,1)
 
                     + 457.33447878*pow(zeta,2) + 348.70391772*pow(zeta,3) 

                     + 170.35400182*pow(zeta,4) + 53.99290691*pow(zeta,5)
 
		     + 10.47985711*pow(zeta,6) + 1.00000000*pow(zeta,7));
	      
		  
	      
	      
		  // The complex refractive index. 
		  N_s[f_grid_index] += S*xi*(-1)*complex_i*CEF;
	      
		 
	      
		  //!With Lorenzian.
		  //N_s[f_grid_index] = S*xi/(f_z - f_grid[f_grid_index] + 
		  //complex_i*agam);
	      
		  //!With Lorenzian - Whiting's approximation.
		  //N_s[f_grid_index] = S*xi/(f_z - f_grid[f_grid_index] + 
		  //complex_i*(agam/2 + pow(agam/2*agam/2 + gamma*gamma,0.5)));
	      
	      
	      
		}
	  
	  
	  
	  
	  
	  
	  
	  
	      // Polarization  matrix accounting for the contributions of 
	      // the 3 different polarizations of the components of the Zeeman 
	      // split due to 3 different values of DeltaM.
	      Matrix P;
	      P.resize(2,2);
	  
	      if (DeltaM==0)
		{
		  P(0,0) =  0.5*sin(phi)*sin(phi);
		  P(0,1) = -0.5*sin(phi)*sin(phi);
		  P(1,0) =  0.5*sin(phi)*sin(phi);
		  P(1,1) = -0.5*sin(phi)*sin(phi);
		}
	      else if (DeltaM==1)
		{
		  P(0,0) = 0.5*(1+cos(phi))*(1+cos(phi));
		  P(0,1) = 0.5*sin(phi)*sin(phi);
		  P(1,0) = 0.5*sin(phi)*sin(phi);
		  P(1,1) = 0.5*(1-cos(phi))*(1-cos(phi));;
		}
	      else if (DeltaM==-1)
		{
		  P(0,0) = 0.5*(1-cos(phi))*(1-cos(phi));
		  P(0,1) = 0.5*sin(phi)*sin(phi);
		  P(1,0) = 0.5*sin(phi)*sin(phi);
		  P(1,1) = 0.5*(1+cos(phi))*(1+cos(phi));
		}
	  
	  
	     
	  

	      // Calculating the halved value of magnetic susceptibility tensor 
	      // at given frequency. It's in the form of an array because 
	      // complex matrices or tensors cannot be handled by ARTS.		
	      Chi[f_grid_index*4+0] += P(0,0)*N_s[f_grid_index];
	      Chi[f_grid_index*4+1] += P(0,1)*N_s[f_grid_index];
	      Chi[f_grid_index*4+2] += P(1,0)*N_s[f_grid_index];
	      Chi[f_grid_index*4+3] += P(1,1)*N_s[f_grid_index];
	  
	  	  
	    }
      
      
      
      
	  
      
      
	  // Calculating the propagation tensor at a given frequency. 
	  // It's in the form of an array because complex matrices or 
	  // tensors cannot be handled by ARTS.
	  G_s[f_grid_index*4+0] = complex_i*k_0*(1. + Chi[f_grid_index*4+0]);
	  G_s[f_grid_index*4+1] = complex_i*k_0*Chi[f_grid_index*4+1];
	  G_s[f_grid_index*4+2] = complex_i*k_0*Chi[f_grid_index*4+2];
	  G_s[f_grid_index*4+3] = complex_i*k_0*(1. + Chi[f_grid_index*4+3]);
      
      
      
	  // First row of the extinction matrix at a given frequency.
	  ext_mat_zee(f_grid_index,0,0) = 
	    0.5*(2.*real(G_s[f_grid_index*4+0]) + 2.*real(G_s[f_grid_index*4+3]));

	  ext_mat_zee(f_grid_index,0,1) = 
	  - 0.5*(2.*real(G_s[f_grid_index*4+3]) - 2.*real(G_s[f_grid_index*4+0]));

	  ext_mat_zee(f_grid_index,0,2) = 
	    0.5*(2.*real(G_s[f_grid_index*4+1]) + 2.*real(G_s[f_grid_index*4+2]));

	  ext_mat_zee(f_grid_index,0,3) = 0.0;
      
	  // Second row of the extinction matrix at a given frequency.
	  ext_mat_zee(f_grid_index,1,0) = 
	  - 0.5*(2.*real(G_s[f_grid_index*4+3]) - 2.*real(G_s[f_grid_index*4+0]));

	  ext_mat_zee(f_grid_index,1,1) = 0.5*(2.*real(G_s[f_grid_index*4+0]) 
					     + 2.*real(G_s[f_grid_index*4+3]));
	  ext_mat_zee(f_grid_index,1,2) = 0.0;
	  ext_mat_zee(f_grid_index,1,3) = 0.5*(2.*imag(G_s[f_grid_index*4+1]) + 
					       2.*imag(G_s[f_grid_index*4+2]));
      
	  // Third row of the extinction matrix.
	  ext_mat_zee(f_grid_index,2,0) = 0.5*(2.*real(G_s[f_grid_index*4+1]) + 
					       2.*real(G_s[f_grid_index*4+2]));
	  ext_mat_zee(f_grid_index,2,1) = 0.0;
	  ext_mat_zee(f_grid_index,2,2) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 
					       2.*real(G_s[f_grid_index*4+3]));
	  ext_mat_zee(f_grid_index,2,3) = -0.5*(2.*imag(G_s[f_grid_index*4+3]) - 
						2.*imag(G_s[f_grid_index*4+0]));
      
	  // Fourth row of the extinction matrix at a given frequency.
	  ext_mat_zee(f_grid_index,3,0) = 0.0; 
	  ext_mat_zee(f_grid_index,3,1) = -0.5*(2.*imag(G_s[f_grid_index*4+1]) + 
						2.*imag(G_s[f_grid_index*4+2]));
	  ext_mat_zee(f_grid_index,3,2) = 0.5*(2.*imag(G_s[f_grid_index*4+3]) - 
					       2.*imag(G_s[f_grid_index*4+0]));
	  ext_mat_zee(f_grid_index,3,3) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 
					       2.*real(G_s[f_grid_index*4+3]));
      
	  // Absorption vector components at a given frequency.
	  abs_vec_zee(f_grid_index,0) = 0.5*(2.*real(G_s[f_grid_index*4+0]) + 
					     2.*real(G_s[f_grid_index*4+3]));
	  abs_vec_zee(f_grid_index,1) = -0.5*(2.*real(G_s[f_grid_index*4+3]) - 
					      2.*real(G_s[f_grid_index*4+0]));
	  abs_vec_zee(f_grid_index,2) = 0.5*(2.*real(G_s[f_grid_index*4+1]) + 
					     2.*real(G_s[f_grid_index*4+2]));
	  abs_vec_zee(f_grid_index,3) = 0.0;
      
	}


  
}

 ******************************************************************************************** */


void CEF( // Input
	 const Complex zeta, 
	 const Numeric Gamma_D,
	 //oUTPUT
	 Numeric& Re_CEF,
	 Numeric& Im_CEF )
{

  // The special line shape, Complex Error Function.
  // see P. W. Rosenkranz Chapter 2 of the Janssen book

  Complex retval;

  retval = (1.000e0 / sqrt(PI) /  Gamma_D) *
             (  122.60793178 * pow(zeta,0.) 
	      + 214.38238869 * pow(zeta,1.) + 181.92853309 * pow(zeta,2.) 
	      +  93.15558046 * pow(zeta,3.) +  30.18014220 * pow(zeta,4.)
	      +   5.91262621 * pow(zeta,5.) +   0.56418958 * pow(zeta,6.)) /
             (  122.60793178 * pow(zeta,0.) + 352.73062511 * pow(zeta,1.)
              + 457.33447878 * pow(zeta,2.) + 348.70391772 * pow(zeta,3.) 
              + 170.35400182 * pow(zeta,4.) +  53.99290691 * pow(zeta,5.)
              +  10.47985711 * pow(zeta,6.) +   1.00000000 * pow(zeta,7.));
    
  Re_CEF = real( retval );
  Im_CEF = imag( retval );

  return;
};



// ########################################################################################


void Zeeman_o2_splitting_factors( Numeric& xi,
				  Numeric& eta,
				  const Numeric N_r,
				  const Numeric DeltaM,
				  const Numeric M)
{


  // Pi polarization for a N+ Zeeman transition.
  if (N_r>0 && DeltaM==0) //
    {
      // Matrix representation of the relative intensity value 
      xi  = 3*((abs(Numeric(N_r))+1)*(abs(Numeric(N_r))+1)-M*M)
	/((abs(Numeric(N_r))+1)*(2*abs(Numeric(N_r))+1)*(2*abs(Numeric(N_r))+3));
      
      // Zeeman frequency shift value
      eta = M*(abs(Numeric(N_r))-1)/(abs(Numeric(N_r))*(abs(Numeric(N_r))+1));
    }
  // Pi polarization for a N- Zeeman transition.
  else if (N_r<0 && DeltaM==0)
    {
      // Matrix representation of the relative intensity values 
      xi  = 3*(abs(Numeric(N_r))*abs(Numeric(N_r))-M*M)
	/(abs(Numeric(N_r))*(2*abs(Numeric(N_r))-1)*(2*abs(Numeric(N_r))+1));
      
      // Zeeman frequency shift values
      eta = M*(abs(Numeric(N_r))+2)/(abs(Numeric(N_r))*(abs(Numeric(N_r))+1));
    }
  // Sigma+ polarization for a N+ Zeeman transition.
  else if (N_r>0 && DeltaM==1)
    {
      // Matrix representation of the relative intensity
      xi  = 3*(abs(Numeric(N_r))+M+1)*(abs(Numeric(N_r))+M+2)/
	(4*(abs(Numeric(N_r))+1)*(2*abs(Numeric(N_r))+1)*(2*abs(Numeric(N_r))+3));
      
      // Zeeman frequency shift value 
      eta = (M*(abs(Numeric(N_r))-1)+abs(Numeric(N_r)))/(abs(Numeric(N_r))*(abs(Numeric(N_r))+1));
    }
  // Sigma- polarization for a N+ Zeeman transition.
  else if (N_r>0 && DeltaM==-1)
    {	
      // Matrix representation of the relative intensity values
      xi  = 3*(abs(Numeric(N_r))-M+1)*(abs(Numeric(N_r))-M+2)/
	(4*(abs(Numeric(N_r))+1)*(2*abs(Numeric(N_r))+1)*(2*abs(Numeric(N_r))+3));
      
      // Zeeman frequency shift value
      eta = (M*(abs(Numeric(N_r))-1)-abs(Numeric(N_r)))/(abs(Numeric(N_r))*(abs(Numeric(N_r))+1));
    }
  // Sigma+ polarization for a N- Zeeman transition.
  else if (N_r<0 && DeltaM==1)
    {
      // Matrix representation of the relative intensity values
      xi  = 3*(abs(Numeric(N_r))-M)*(abs(Numeric(N_r))-M-1)/
	(4*abs(Numeric(N_r))*(2*abs(Numeric(N_r))-1)*(2*abs(Numeric(N_r))+1));
      
      // Zeeman frequency shift value 
      eta = (M*(abs(Numeric(N_r))+2)+(abs(Numeric(N_r))+1))/(abs(Numeric(N_r))*(abs(Numeric(N_r))+1));
    }
  // Sigma- polarization for a N- Zeeman transition.
  else if (N_r<0 && DeltaM==-1)
    {
      // Matrix representation of the relative intensity values
      xi  = 3*(abs(Numeric(N_r))+M)*(abs(Numeric(N_r))+M-1)/
	(4*abs(Numeric(N_r))*(2*abs(Numeric(N_r))-1)*(2*abs(Numeric(N_r))+1));
      
      // Zeeman frequency shift value 
      eta = (M*(abs(Numeric(N_r))+2)-(abs(Numeric(N_r))+1))/(abs(Numeric(N_r))*(abs(Numeric(N_r))+1));
    }
  else 
    {
      throw runtime_error("Zeeman_o2_splitting_factors: Error in N_r or DeltaM");
    };
  
  
  return;
};



// ########################################################################################



// Propagation tensor G_s at a given frequency is in the form of an array 
// because complex matrices or tensors cannot be handled by ARTS. 
void Zeeman_o2_line_splitting(// Output:
			      Matrix& ext_mat_tmp,           // temporary 4x4 extinction matrix
			      Vector& abs_vec_tmp,           // temporary 1x4 absorption vector
			      // Input:
			      const Numeric B_field,         // magnitude of Earth's magnetic field along LOS in [T].
			      const Numeric phi,             // angle between the mag. field and LOS [radians]  
			      const Numeric f_c,             // spectral line center frequencies
			      const Numeric S,               // spectral line intensities 
			      const Numeric gamma_L,         // spectral line pressure broadening
			      const Numeric gamma_D,         // spectral line  Doppler broadening 
			                                     // with missing central frequency
			      const Index   N_r,             // Rotational quantum number of O2 line.
			      const Numeric f_grid_point)    // frequency grid point.
  
{
  
  
  // check if qunatum number is correct
  // the models of Liebe and Rosenkranz have parameters up to +/-39 and +/-33, respectively.
  if ( (abs(Numeric(N_r)) < 1) || (abs(Numeric(N_r)) > 39) )
    {
      ostringstream os;
      os << "Zeeman_o2_line_splitting: wrong quantum number N_r given in input. \n"
	 << "valid range is +/- 1 <= " << N_r << " <= 39 \n";
      throw runtime_error(os.str());
    };



  // wavenumber claculation.
  Numeric k_0 = twoPI_c * f_grid_point;

  // Complex refractive index.
  Numeric Re_N_s=0.000e0;
  Numeric Im_N_s=0.000e0;
  
  // Propagation tensor G at a given frequency is in the form of an array 
  // because complex matrices or tensors cannot be handled by ARTS. 
  // Array<Complex> G_s(2*2);
  Matrix G_s;
  G_s.resize(4,2); // (4 Stokes components) x (real+imaginary parts)
  G_s = 0.000e0;

  // Magnetic susceptibility tensor at a given frequency. 
  // It's in the form of an array because complex matrices
  // or tensors cannot be handled by ARTS.
  Matrix Chi;
  Chi.resize(4,2); // (4 Stokes components) x (real+imaginary parts)
  Chi = 0.0000;
  
  
  // Polarization  matrix accounting for the contributions of 
  // the 3 different polarizations of the components of the Zeeman 
  // split due to 3 different values of DeltaM.
  Tensor3 P;
  P.resize(3,2,2);
  // Polarization  matrix accounting for the contributions of 
  // the 3 different polarizations of the components of the Zeeman 
  // split due to 3 different values of DeltaM.      
  // DeltaM==0
  P(0,0,0) =  0.5 *    sin(phi)  *    sin(phi);
  P(0,0,1) = -0.5 *    sin(phi)  *    sin(phi);
  P(0,1,0) =  0.5 *    sin(phi)  *    sin(phi);
  P(0,1,1) = -0.5 *    sin(phi)  *    sin(phi);
  // DeltaM==1
  P(1,0,0) =  0.5 * (1+cos(phi)) * (1+cos(phi));
  P(1,0,1) =  0.5 *    sin(phi)  *    sin(phi);
  P(1,1,0) =  0.5 *    sin(phi)  *    sin(phi);
  P(1,1,1) =  0.5 * (1-cos(phi)) * (1-cos(phi));
  // DeltaM==-1
  P(2,0,0) =  0.5 * (1-cos(phi)) * (1-cos(phi));
  P(2,0,1) =  0.5 *    sin(phi)  *    sin(phi);
  P(2,1,0) =  0.5 *    sin(phi)  *    sin(phi);
  P(2,1,1) =  0.5 * (1+cos(phi)) * (1+cos(phi));


  // variables for complex error function / Faddeeva function
  Vector  FX;
  FX.resize(1);
  FX=0.00e0;
  Numeric FY=0.00e0;
  Vector  FK;
  FK.resize(1);
  FK=0.00e0;
  Vector  FL; 
  FL.resize(1);
  FL=0.00e0;
  

  // related quantities to quantum number N_r
  Numeric AN_r = 0.00e0;
  Numeric BN_r = 0.00e0;
  if (N_r>0)
    { 
      // Assigning a value of the total angular momentum N of the lower  
      // level for a N+ Zeeman transition (N -> N+1).
      AN_r = abs(Numeric(N_r));
      
      // Assigning a value of the number of the split components for each 
      // of the three possible polarizations of a N+ Zeeman transition.
      BN_r = 2*abs(Numeric(N_r))+1;
    }
  else if (N_r<0)
    {
      // Assigning a value of the total angular momentum N of the lower 
      // level for a N- Zeeman transition (N -> N-1).
      AN_r = abs(Numeric(N_r))-1;
      
      // Assigning a value of the number of the split components for each
      // of the three possible polarizations of a N- Zeeman transition.
      BN_r = 2*abs(Numeric(N_r))-1;
    };
  
  
  // Starting the loop over the 3 possible values of the change of the 
  // magnetic quantum number M upon a Zeeman transition.	  
  for (Index k=0; k<3; k++)
    { 
      Index DeltaM = k-1;
      
      // Starting the loop over 2N+1 or 2N-1 values of
      // the magnetic quantum number M for a given transition 
      for (Index j=0; j<BN_r; j++)
	{ 
	  Numeric M = (Numeric)j - AN_r;
	  
	  // calculate the splitting factors for the splitted line intensity 
	  // and frequency shift.
	  // (1) intensity factor, or relative intensity, of the 
	  //     individual components of the Zeeman split relative 
	  //     to the the intensity of the unsplit line.
	  Numeric xi = 0.00e0;    
	  // (2) this the Zeeman frequency shift of the individual 
	  //     components of the Zeeman split.
	  Numeric eta = 0.00e0;
	  Zeeman_o2_splitting_factors( xi, eta, N_r, DeltaM, M);

	  // Center frequency of the individual Zeeman components
	  Numeric f_z;
	  f_z = f_c + 28.03e0 * eta * B_field;  // [GHz]
	  
	  // Complex argument of the special line shape (CEF), 
	  // used in the case of Zeeman splitting.
	  Numeric Gamma_D = gamma_D * f_z;  // [GHz]
	  
	  // complex error function / Faddeeva function
	  // (SQRT(ln(2)) = 0.8325546e0)
	  FX[0] = 0.8325546e0 * (f_grid_point - f_z) / Gamma_D;
	  FY    = 0.8325546e0 * gamma_L / Gamma_D;
	  HUMLIK_Faddeeva_Wells( FX, FY, FK, FL );
	  Re_N_s +=  S * xi * FL[0]; // real      part of complex refractive index
	  Im_N_s += -S * xi * FK[0]; // imaginary part of complex refractive index
	  cout << "1******************************************************1 \n";
	  cout << "f_z=" << f_z << ", f_grid_point="      << f_grid_point << "\n"  
               << "Gamma_D=" << Gamma_D << ", gamma_L="   << gamma_L      << "\n";
	  cout << "WEL:   X=" << FX    << ", Y=" << FY    << "\n";
	  cout << "WEL:   K=" << FK[0] << ", L=" << FL[0] << "\n";
	  cout << "1******************************************************1 \n";

	  // The complex refractive index.
	  // Complex zeta = Complex(gamma_L/Gamma_D, ((f_grid_point - f_z)/Gamma_D));
	  Complex zeta = Complex((gamma_L/Gamma_D),((f_grid_point - f_z)/Gamma_D) );
	  /*
	  cout << "f_z=" << f_z << ", f_grid_point=" << f_grid_point 
               << ", Re(zeta)=" << real(zeta) << ", Im(zeta)=" << imag(zeta) << "\n" 
               << "gamma_D=" << gamma_D << ", Gamma_D=" << Gamma_D 
               << ", gamma_L=" << gamma_L << "\n";
	  */

	  Numeric KK;
	  Numeric LL;
	  CEF( zeta, Gamma_D, KK, LL );
	  // cout << "CEF:   KK=" << KK << ", LL=" << LL << "\n";
	  // Dummy = S * xi * (-1) * complex_i*
	  // Re_N_s +=  S * xi * LL; // real      part of complex refractive index
	  //Im_N_s += -S * xi * KK; // imaginary part of complex refractive index

	  // Calculating the halved value of magnetic susceptibility tensor 
	  // at given frequency. It's in the form of an array because 
	  // complex matrices or tensors cannot be handled by ARTS.		
	  Chi(0,0) += P(k,0,0) * Re_N_s;
	  Chi(0,1) += P(k,0,0) * Im_N_s;

	  Chi(1,0) += P(k,0,1) * Re_N_s;
	  Chi(1,1) += P(k,0,1) * Im_N_s;
	  
	  Chi(2,0) += P(k,1,0) * Re_N_s;
	  Chi(2,1) += P(k,1,0) * Im_N_s;

	  Chi(3,0) += P(k,1,1) * Re_N_s;
	  Chi(3,1) += P(k,1,1) * Im_N_s;
	};
    };
  
  
  // Calculating the propagation tensor at a given frequency. 
  // It's in the form of an array because complex matrices or 
  // tensors cannot be handled by ARTS.
  G_s(0,0) = -k_0 * Chi(0,1);
  G_s(0,1) =  k_0 * (1.000e0 + Chi(0,0));

  G_s(1,0) = -k_0 * Chi(1,1);
  G_s(1,1) =  k_0 * Chi(1,0);

  G_s(2,0) = -k_0 * Chi(2,1);
  G_s(2,1) =  k_0 * Chi(2,0);

  G_s(3,0) = -k_0 * Chi(3,1);
  G_s(3,1) =  k_0 * (1.000e0 + Chi(3,0));
  

  // calculate elements of extinction/absorption
  Numeric AA = 2.000e0  * G_s(0,0);
  Numeric BB = G_s(1,0) + G_s(2,0);
  Numeric CC = G_s(0,0) - G_s(3,0);
  Numeric DD = G_s(0,1) - G_s(3,1);
  Numeric EE = G_s(1,1) + G_s(1,1);

  // Absorption vector components at a given frequency.
  abs_vec_tmp[0]   +=  2.000e0 * AA;
  abs_vec_tmp[1]   +=  2.000e0 * BB;
  abs_vec_tmp[3]   +=  2.000e0 * CC;

  // First row of the extinction matrix at a given frequency.
  ext_mat_tmp(0,0) += -AA;
  ext_mat_tmp(0,1) += -BB;
  ext_mat_tmp(0,3) += -CC;
  
  // Second row of the extinction matrix at a given frequency.
  ext_mat_tmp(1,0) += -BB;
  ext_mat_tmp(1,1) += -AA;
  ext_mat_tmp(1,2) +=  DD;
  
  // Third row of the extinction matrix.
  ext_mat_tmp(2,1) += -DD;
  ext_mat_tmp(2,2) += -AA;
  ext_mat_tmp(2,3) += -EE;
  
  // Fourth row of the extinction matrix at a given frequency.
  ext_mat_tmp(3,0) += -CC;
  ext_mat_tmp(3,2) +=  EE;
  ext_mat_tmp(3,3) += -AA;
  
  return;

};



// #################################################################################


void PWRO2Mixing(//output
		 Matrix& ext_mat_tmp, 
		 Vector& abs_vec_tmp,
		 // Input
		 Numeric STR,
		 Numeric Y,
		 Numeric DF,
		 Numeric FO,
		 Numeric ff,
		 Numeric vmro2,
		 Numeric p,
                 Numeric TH)
{

  // call VVW function for line mixing but without Zeeman splitting.
  // This is sufficient for the troposhere where ther individual lines 
  // are too wide that line mixing and pressure broadening dominate over 
  // Zeeman splitting.
  Numeric SF1      = ( DF + (ff-FO)*Y ) / ( (ff-FO)*(ff-FO) + DF*DF );
  Numeric SF2      = ( DF - (ff+FO)*Y ) / ( (ff+FO)*(ff+FO) + DF*DF );

  // sum the line absorption part for a specific spectral line on the frequency grid
  Numeric SUM      = STR * (SF1+SF2) * (ff/FO) * (ff/FO);

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
  Numeric abs_line = 2.414322e7 / PI * vmro2   * p * pow(TH, (Numeric)3.0) * SUM;

  // fill I/O extinction matrix and absorption vector properly
  abs_vec_tmp[0]   +=  abs_line;
  for (Index k=0; k<4; ++k) ext_mat_tmp(k,k) +=  abs_vec_tmp[0];
    
  return;
}



// #################################################################################



void PWRO2VoigtMixing(//output
		      Matrix& ext_mat_tmp, 
		      Vector& abs_vec_tmp,
		      // Input
		      const Numeric STR,
		      const Numeric Y,
		      const Numeric DF,
		      const Numeric gamma_D,
		      const Numeric FO,
		      const Numeric ff,
		      const Numeric vmro2,
		      const Numeric p,
		      const Numeric TH)
{


  // coonstant = sqrt(ln(2) / pi)
  const Numeric lnpi = 0.4697186e0;

  // variables for complex error function / Faddeeva function
  Vector  FX;
  FX.resize(1);
  FX=0.00e0;
  Numeric FY=0.00e0;
  Vector  FK;
  FK.resize(1);
  FK=0.00e0;
  Vector  FL; 
  FL.resize(1);
  FL=0.00e0;


  // complex error function / Faddeeva function
  // (SQRT(ln(2)) = 0.8325546e0)
  FX[0] = 0.8325546e0 * (ff - FO) / gamma_D;
  FY    = 0.8325546e0 * DF / gamma_D;
  HUMLIK_Faddeeva_Wells( FX, FY, FK, FL );

  Numeric SF1 = (lnpi / gamma_D * FK[0]) + FL[0]*Y;


  // complex error function / Faddeeva function
  // (SQRT(ln(2)) = 0.8325546e0)
  FX[0] = 0.8325546e0 * (ff + FO) / gamma_D;
  FY    = 0.8325546e0 * DF / gamma_D;
  HUMLIK_Faddeeva_Wells( FX, FY, FK, FL );

  Numeric SF2 = (lnpi / gamma_D * FK[0]) + FL[0]*Y;

  /*
  cout << "2******************************************************2 \n";
  cout << "FO="       << FO      << ", ff ="  << ff    << "\n"  
       << "gamma_D="  << gamma_D << ", DF ="  << DF    << "\n";
  cout << "WEL:   X=" << FX[0]   << ", Y  ="  << FY    << "\n";
  cout << "WEL:   K=" << FK[0]   << ", L  ="  << FL[0] << "\n";
  cout << "WEL: SF1=" << SF1     << ", SF2= " << SF2   << "\n";
  cout << "2******************************************************2 \n";
  */

  // sum the line absorption part for a specific spectral line on the frequency grid
  Numeric SUM      = STR * (SF1+SF2) * (ff/FO);

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
  Numeric abs_line = 2.414322e7 / PI * vmro2   * p * pow(TH, (Numeric)3.0) * SUM;

  // fill I/O extinction matrix and absorption vector properly
  abs_vec_tmp[0]   +=  abs_line;
  for (Index k=0; k<4; ++k) ext_mat_tmp(k,k) +=  abs_vec_tmp[0];
    
  return;
}



// #################################################################################



//!   P. W. Rosenkranz oxygen absorption model
/*!
  See arts -d absPWRO2Model for detailed documentation.
  
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
  
  \retval   abs                  absorption coefficient of oxygen [1/m]
  \param    f_grid               predefined frequency grid        [Hz]
  \param    abs_p                predefined pressure              [Pa]
  \param    abs_t                predefined temperature grid      [K] 
  \param    abs_vmr              H2O volume mixing ratio profile  [1]
  \param    model                model version of the P. W. Rosenkranz oxygen absorption model
  \param    abs_user_parameters  scaling factor(s)

  \note     Except for  model 'user' the input parameters CCin, CLin, CWin, and COin 
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
               John Wiley & Sons, Inc., 1993.
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

Index absPWRO2Model(// WS Output:
		    Matrix&         ext_mat_tmp,          // Tensor3 of the Extinction Matrix [1/m]. 
		    Vector&         abs_vec_tmp,          // Matrix of the Absorption Vector  [1/m].
		    // WS Input:
		    const Numeric   geomag_strength,      // mag. field strength          [Gauss]
		    const Numeric   geomag_angle,         // mag. field orientation angle [radians]
		    const Index     zeeman_o2_onoff,           // Zeeman splitting on or off
                    const Numeric   zeeman_o2_pressure_limit,  // Zeeman pressure limit   [Pa]
		    const Numeric   f_grid_point,         // frequency vector             [Hz]
		    const Numeric   p,                    // pressure                     [Pa]
		    const Numeric   t,                    // temperature                  [K]
		    const Numeric   vmro2,                // O2 vmr                       [1]
		    const Numeric   vmrh2o,               // H2O vmr                      [1]
		    const String&   abs_model,            // model selection string
		    const Vector&   abs_user_parameters ) // scaling factor(s)
{


  // return value:  1=true, 0=false
  Index retval = 1;


  // molecular weight of oxygen (main isotope) [g/mol]
  // const Numeric MWO2 = 31.989830e0; 


  // isotopic ratio for O^18O^18
  const Numeric ISORATIO = 0.995262;


  // number of spectral lines in the different model versions
  const Index n_lines_PWR88 = 40; // all O2 lines (range: 50-850 GHz)
  const Index n_lines_PWR93 = 40; // all O2 lines (range: 50-850 GHz)
  const Index n_lines_PWR98 = 40; // all O2 lines (range: 50-850 GHz)


  // rotational quantum number of the spectral lines
  const Index QM93[n_lines_PWR93]  = {  -1,       +1,       -3,       +3,
                                        -5,       +5,       -7,       +7,
                                        -9,       +9,       -11,      +11,
                                        -13,      +13,      -15,      +15,
                                        -17,      +17,      -19,      +19,
                                        -21,      +21,      -23,      +23,
                                        -25,      +25,      -27,      +27,
                                        -29,      +29,      -31,      +31,
                                        -33,      +33,       0,        0,
                                        -0,        0,        0,        0};

  // rotational quantum number of the spectral lines
  const Index QM98[n_lines_PWR98]  = {  -1,       +1,       -3,       +3,
                                        -5,       +5,       -7,       +7,
                                        -9,       +9,       -11,      +11,
                                        -13,      +13,      -15,      +15,
                                        -17,      +17,      -19,      +19,
                                        -21,      +21,      -23,      +23,
                                        -25,      +25,      -27,      +27,
                                        -29,      +29,      -31,      +31,
                                        -33,      +33,       0,        0,
                                        -0,        0,        0,        0};

  // line center frequencies for the model version 1993
  const Numeric F93[n_lines_PWR93] = { 118.7503,  56.2648,  62.4863,  58.4466, 
				        60.3061,  59.5910,  59.1642,  60.4348, 
				        58.3239,  61.1506,  57.6125,  61.8002, 
				        56.9682,  62.4112,  56.3634,  62.9980, 
				        55.7838,  63.5685,  55.2214,  64.1278, 
				        54.6712,  64.6789,  54.1300,  65.2241, 
				        53.5957,  65.7648,  53.0669,  66.3021, 
				        52.5424,  66.8368,  52.0214,  67.3696, 
				        51.5034,  67.9009, 368.4984, 424.7631, 
				       487.2494, 715.3932, 773.8397, 834.1453};

  // line center frequencies for the model version 1998
  const Numeric F98[n_lines_PWR98] = { 118.7503,  56.2648,  62.4863,  58.4466,  
                                        60.3061,  59.5910,  59.1642,  60.4348,  
				        58.3239,  61.1506,  57.6125,  61.8002,
				        56.9682,  62.4112,  56.3634,  62.9980,  
				        55.7838,  63.5685,  55.2214,  64.1278,
                                        54.6712,  64.6789,  54.1300,  65.2241,
				        53.5957,  65.7648,  53.0669,  66.3021,
                                        52.5424,  66.8368,  52.0214,  67.3696,  
				        51.5034,  67.9009, 368.4984, 424.7632,
				       487.2494, 715.3931, 773.8397, 834.1458};


  // line strength (taken from HITRAN92) for the model version 1993 at T=300K [cm² * Hz]
  const Numeric S93[n_lines_PWR93] = {  0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14,
				        0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14,
				        0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,
				        0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14,
				        0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14,
				        0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,
				        0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15,
				        0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16,
				        0.8191E-17, 0.1423E-16, 0.6460E-15, 0.7047E-14, 
				        0.3011E-14, 0.1826E-14, 0.1152E-13, 0.3971E-14};

  // line strength (intensities in the submm range are updated according to HITRAN96)
  // for the model version 1998 at T=300K [cm² * Hz]
  const Numeric S98[n_lines_PWR98] = {  0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14,
				        0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14,
				        0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,
				        0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14,
				        0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14,
				        0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,
				        0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15,
				        0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16,
				        0.8191E-17, 0.1423E-16, 0.6494E-15, 0.7083E-14, 
				        0.3025E-14, 0.1835E-14, 0.1158E-13, 0.3993E-14};

  // line mixing y parameter for the calculation of Y [1/bar]
  const Numeric Y93[n_lines_PWR98] = { -0.0233,  0.2408, -0.3486,  0.5227,
				       -0.5430,  0.5877, -0.3970,  0.3237, 
				       -0.1348,  0.0311,  0.0725, -0.1663, 
					0.2832, -0.3629,  0.3970, -0.4599,
					0.4695, -0.5199,  0.5187, -0.5597,
					0.5903, -0.6246,  0.6656, -0.6942,
					0.7086, -0.7325,  0.7348, -0.7546,
					0.7702, -0.7864,  0.8083, -0.8210,
					0.8439, -0.8529,  0.0000,  0.0000,
					0.0000,  0.0000,  0.0000,  0.0000};
  
  // y parameter for the calculation of Y [1/bar].
  // These values are from P. W. Rosenkranz, Interference coefficients for the 
  // overlapping oxygen lines in air, JQSRT, 1988, Volume 39, 287-297.
  const Numeric Y88[n_lines_PWR88]  = {-0.0244,  0.2772, -0.4068,  0.6270,
				       -0.6183,  0.6766, -0.4119,  0.3290, 
				        0.0317, -0.1591,  0.1145, -0.2068, 
				        0.3398, -0.4158,  0.3922, -0.4482,
				        0.4011, -0.4442,  0.4339, -0.4687,
				        0.4783, -0.5074,  0.5157, -0.5403,
				        0.5400, -0.5610,  0.5719, -0.5896,
				        0.6046, -0.6194,  0.6347, -0.6468, 
				        0.6627, -0.6718,  0.0000,  0.0000,
				        0.0000,  0.0000,  0.0000,  0.0000};

  // temperature exponent of the line strength in [1]
  const Numeric BE[n_lines_PWR93] = {   0.009,   0.015,   0.083,   0.084, 
				        0.212,   0.212,   0.391,   0.391, 
			                0.626,   0.626,   0.915,   0.915, 
			        	1.260,   1.260,   1.660,   1.665,   
                                        2.119,   2.115,   2.624,   2.625, 
                                        3.194,   3.194,   3.814,   3.814, 
                                        4.484,   4.484,   5.224,   5.224, 
                                        6.004,   6.004,   6.844,   6.844, 
                                        7.744,   7.744,   0.048,   0.044, 
	      		                0.049,   0.145,   0.141,   0.145};

  // line width parameter [GHz/bar]
  const Numeric W300[n_lines_PWR93] = { 1.630,   1.646,   1.468,   1.449, 
				        1.382,   1.360,   1.319,   1.297, 
			                1.266,   1.248,   1.221,   1.207, 
			                1.181,   1.171,   1.144,   1.139, 
			                1.110,   1.108,   1.079,   1.078, 
			                1.050,   1.050,   1.020,   1.020, 
			                1.000,   1.000,   0.970,   0.970,
			                0.940,   0.940,   0.920,   0.920,
			                0.890,   0.890,   1.920,   1.920, 
			                1.920,   1.810,   1.810,   1.810};


  // v parameter for the calculation of Y [1/bar]
  const Numeric V[n_lines_PWR93] ={     0.0079, -0.0978,  0.0844, -0.1273,
				        0.0699, -0.0776,  0.2309, -0.2825, 
			      	        0.0436, -0.0584,  0.6056, -0.6619, 
				        0.6451, -0.6759,  0.6547, -0.6675,
				        0.6135, -0.6139,  0.2952, -0.2895, 
				        0.2654, -0.2590,  0.3750, -0.3680, 
				        0.5085, -0.5002,  0.6206, -0.6091,
				        0.6526, -0.6393,  0.6640, -0.6475,
				       0.6729, -0.6545,  0.0000,  0.0000,
				       0.0000,  0.0000,  0.0000,  0.0000};

  // standard scaling factors (=unity) for the original Rosenkranz models 
  const Numeric CC_PWR = 1.00000e0; // scaling factor of the O2 continuum
  const Numeric CL_PWR = 1.00000e0; // scaling factor of the O2 line intensity
  const Numeric CW_PWR = 1.00000e0; // scaling factor of the O2 line width parameter
  const Numeric CO_PWR = 1.00000e0; // scaling factor of the O2 line mixing parameter


  // widths in MHz/mbar for the O2 continuum
  const Numeric WB300 = 0.56; // [MHz/mbar]=[MHz/hPa]
  const Numeric X     = 0.80; // [1]
  const Numeric CONTFAC = 1.2300e-10;


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  Index n_lines = 0;
  int   flag98  = 1;
  const Numeric *F;
  const Numeric *S300;
  const Numeric *Y300;
  const Index *QM;


  if ( abs_model == "PWR88" )       // ---------- model version 1988
    {
      // set the scaling factors
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = CO_PWR; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR88;
      // fill in the appropriate version parameters
      F    = F93;  // line center frequencies from 1993 version
      S300 = S93;  // line intensity from 1993 version
      Y300 = Y88;  // line mixing from 1988 version
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR88scaling" )
    {
      // check if the input parameter vector has appropriate size
      if ( abs_user_parameters.nelem() != 1 )
	throw runtime_error("Method absPWRO2Model: Vector abs_user_parameters must have 1 element.");
      // set the scaling factors
      CC = abs_user_parameters[0]; // scaling factor of the O2 continuum
      CL = abs_user_parameters[1]; // scaling factor of the O2 line intensity
      CW = abs_user_parameters[2]; // scaling factor of the O2 line width parameter
      CO = abs_user_parameters[3]; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR88;
      // fill in the appropriate version parameters
      F    = F93;  // line center frequencies from 1993 version
      S300 = S93;  // line intensity from 1993 version
      Y300 = Y88;  // line mixing from 1988 version
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR93" )  // ---------- model version 1993
    {
      // set the scaling factors
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = CO_PWR; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR93;
      // fill in the appropriate version parameters
      F    = F93;
      S300 = S93;
      Y300 = Y93;
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR93Lines" )
    {
      // set the scaling factors
      CC = 0.000e0; // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = CO_PWR; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR93;
      // fill in the appropriate version parameters
      F    = F93;
      S300 = S93;
      Y300 = Y93;
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR93Continuum" )
    {
      // set the scaling factors
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = 0.000e0; // scaling factor of the O2 line intensity
      CW = 0.000e0; // scaling factor of the O2 line width parameter
      CO = 0.000e0; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR93;
      // fill in the appropriate version parameters
      F    = F93;
      S300 = S93;
      Y300 = Y93;
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR93NoCoupling" )
    {
      // set the scaling factors
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = 0.000e0; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR93;
      // fill in the appropriate version parameters
      F    = F93;
      S300 = S93;
      Y300 = Y93;
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR93scaling" )
    {
      // check if the input parameter vector has appropriate size
      if ( abs_user_parameters.nelem() != 1 )
	throw runtime_error("Method absPWRO2Model: Vector abs_user_parameters must have 1 element.");
      // set the scaling factors
      CC = abs_user_parameters[0]; // scaling factor of the O2 continuum
      CL = abs_user_parameters[1]; // scaling factor of the O2 line intensity
      CW = abs_user_parameters[2]; // scaling factor of the O2 line width parameter
      CO = abs_user_parameters[3]; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR93;
      // fill in the appropriate version parameters
      F    = F93;
      S300 = S93;
      Y300 = Y93;
      QM   = QM93; // line quantum number
    }
  else if ( abs_model == "PWR98" )  // ---------- model version 1998
    {
      // set flag for 1998 version identification
      flag98 = 0;
      // set the scaling factors
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = CO_PWR; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR98;
      // fill in the appropriate version parameters
      F    = F98;
      S300 = S98;
      Y300 = Y93;
      QM   = QM98; // line quantum number
    }
  else if ( abs_model == "PWR98Lines" )
    {
      // set flag for 1998 version identification
      flag98 = 0;
      // set the scaling factors
      CC = 0.000;  // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = CO_PWR; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR98;
      // fill in the appropriate version parameters
      F    = F98;
      S300 = S98;
      Y300 = Y93;
      QM   = QM98; // line quantum number
    }
  else if ( abs_model == "PWR98Continuum" )
    {
      // set flag for 1998 version identification
      flag98 = 0;
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = 0.000;  // scaling factor of the O2 line intensity
      CW = 0.000;  // scaling factor of the O2 line width parameter
      CO = 0.000;  // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR98;
      // fill in the appropriate version parameters
      F    = F98;
      S300 = S98;
      Y300 = Y93;
      QM   = QM98; // line quantum number
    }
  else if ( abs_model == "PWR98NoCoupling" )
    {
      // set flag for 1998 version identification
      flag98 = 0;
      // set the scaling factors
      CC = CC_PWR; // scaling factor of the O2 continuum
      CL = CL_PWR; // scaling factor of the O2 line intensity
      CW = CW_PWR; // scaling factor of the O2 line width parameter
      CO = 0.000;  // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR98;
      // fill in the appropriate version parameters
      F    = F98;
      S300 = S98;
      Y300 = Y93;
      QM   = QM98; // line quantum number
    }
  else if ( abs_model == "PWR98scaling" )
    {
      // set flag for 1998 version identification
      flag98 = 0;
      // check if the input parameter vector has appropriate size
      if ( abs_user_parameters.nelem() != 1 )
	throw runtime_error("Method absPWRO2Model: Vector abs_user_parameters must have 1 element.");
      // set the scaling factors 
      CC = abs_user_parameters[0]; // scaling factor of the O2 continuum
      CL = abs_user_parameters[2]; // scaling factor of the O2 line intensity
      CW = abs_user_parameters[3]; // scaling factor of the O2 line width parameter
      CO = abs_user_parameters[4]; // scaling factor of the O2 line mixing parameter
      // number of spectral lines in the model
      n_lines = n_lines_PWR98;
      // fill in the appropriate version parameters
      F    = F98;
      S300 = S98;
      Y300 = Y93;
      QM   = QM98; // line quantum number
    }
  else
    {
      ostringstream os;
      os << "Method absPWRO2Model: wrong model values given.\n"
	 << "Valid models are: 'Rosenkranz', 'RosenkranzLines', RosenkranzContinuum, " 
         << "'RosenkranzNoCoupling', and 'user'" << '\n';
      throw runtime_error(os.str());
    }


  // initialization of abs/ext arrays with zero
  abs_vec_tmp = 0.000e0;
  ext_mat_tmp = 0.000e0;


  // range of lines to take into account for the line absorption part
  //const Index first_line = 0;         // first line for calculation
  //const Index last_line  = n_lines-1; // last line for calculation
  
  const Index first_line = 9;         // first line for calculation
  const Index last_line  = 9;        // last line for calculation
  

  // relative inverse temperature [1]
  Numeric TH     = 3.0000e2 / t;
  Numeric TH1    = (TH-1.000e0);
  Numeric B      = pow(TH, X);
  // partial pressure of H2O and dry air [hPa]
  Numeric PRESWV = Pa_to_hPa * (p * vmrh2o);
  Numeric PRESDA = Pa_to_hPa * (p * (1.000e0 - vmrh2o));
  Numeric DEN    = 0.001*(PRESDA*B + 1.1*PRESWV*TH); // [hPa]
  Numeric DENS   = 0.001*(PRESDA   + 1.1*PRESWV)*TH; // [hPa]
  Numeric DFNR   = WB300*DEN;                        // [GHz]
  

  // continuum absorption [1/m/GHz]
  Numeric CCONT  = CC * CONTFAC * pow( TH, (Numeric)2.00e0 ) * p;
  

  // input frequency in [GHz]
  Numeric ff   = Hz_to_GHz * f_grid_point; 
      

  // O2 continuum absorption [1/m]
  Numeric abscont = CCONT * (ff * ff * DFNR / (ff*ff + DFNR*DFNR));


  // initialize the absorption vector with the O2 continuum
  abs_vec_tmp[0] =  abscont;


  // loop over the spectral line center frequencies
  for ( Index l=first_line; l<=last_line; ++l )
    {
      
      // --- line pressure broadening [GHz] ------------------------------
      Numeric DF = 0.000e0;
      // change in 118GHz line in 1998 version compared to older versions
      if ( (flag98 == 0) && (fabs((F[l]-118.75e0)) < 0.10e0) )
	{
	  // 118 line update in '98 version according to 
	  // Ph.D. Thesis of M. J. Schwartz, MIT, 1997
	  DF = CW * W300[l] * DENS; // [GHz]
	}
      else
	{
	  // versions 1993 and 1988 have a simpler temperature dependence
	  DF = CW * W300[l] * DEN; // [GHz]
	};
      
      // --- Doppler line broadening part [1] ----------------------------
      Numeric gamma_D = 1.09600e-6 / sqrt(TH);
      
      // --- line mixing parameter [1] -----------------------------------
      Numeric Y       = CO * 0.001e0 * 0.010e0 * p * B * ( Y300[l] + V[l]*TH1 );
      
      // --- line strength parameter [1] ---------------------------------
      Numeric STR     = CL * S300[l] * exp(-BE[l] * TH1);
      // line intensity parameter in units of the MPM model
      // For conversion factors see Peter Rayer, "NWP SAF Microwave 
      // transmittance models for RTTOV", version 1, 12 Feb 2001.
      Numeric A1      = STR / 0.041907e0 / 3.141592e0 / F[l] * 0.50326e12; 
      
      // Zeeman effect on/off switch
      if ( (zeeman_o2_onoff == 1)           && 
           (p < zeeman_o2_pressure_limit)   &&
           ( abs(QM[l]) > 0) )
	{
	  cout << "***absPWRO2Model***  Zeeman,  QM=" << QM[l] << "\n";
	  // call Zeeman splitting function for a single spectral line
	  Zeeman_o2_line_splitting(ext_mat_tmp,     abs_vec_tmp,
				   geomag_strength, geomag_angle, 
				   F[l],            A1, 
				   DF,              gamma_D,
				   QM[l],
				   ff);
	}
      else
	{
	  // call the original P. W. Rosenkranz model part with VVW line shape 
	  // and line mixing. This part is for the troposphere (high pressure) 
	  cout << "----------------------------------------------------------------------------\n";
	  cout << "***absPWRO2Model***  VVW-mixing \n";
	  PWRO2Mixing(ext_mat_tmp, abs_vec_tmp,	
                      STR, Y, DF, F[l], ff,
		      (ISORATIO * vmro2), p,  TH);

	  cout << "%%%  PWRO2Mixing      abs=" << abs_vec_tmp[0] << "\n";

	  PWRO2VoigtMixing(ext_mat_tmp, abs_vec_tmp,	
			   STR, Y, DF, (gamma_D*F[l]), 
                           F[l], f_grid_point,
			   (ISORATIO * vmro2), p,  TH);

	  cout << "%%%  PWRO2VoigtMixing abs=" << abs_vec_tmp[0] << "\n";
	};
    };
  
  /*
  cout << "ext_mat_tmp(0,0)=" << ext_mat_tmp(0,0) << "\n"
       << "ext_mat_tmp(0,1)=" << ext_mat_tmp(0,1) << "\n"
       << "ext_mat_tmp(0,2)=" << ext_mat_tmp(0,2) << "\n"
       << "ext_mat_tmp(0,3)=" << ext_mat_tmp(0,3) << "\n";
  */

  return retval;  // 1=true, 0=false
};



// #################################################################################



//!   oxygen absorption models
/*!
  See arts -d absO2Model for detailed documentation.
  
  \retval   abs            absorption/extinction of oxygen  [1/m]

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
               John Wiley & Sons, Inc., 1993.
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
  
  \retval   abs            absorption coefficient of oxygen [1/m]
  
  \param    geomag_los     geomagnetic field B, strength and angle
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
               John Wiley & Sons, Inc., 1993.
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
void absO2ZeemanModel(// WS Output:
		      // define internal absorption arrays
		      Tensor3&         ext_mat_zee,          // Tensor3 of the Extinction Matrix [1/m]. 
		      Matrix&          abs_vec_zee,          // Matrix of the Absorption Vector  [1/m].
		      // WS Input:
		      const Matrix&    geomag_los,           // [Magnetic field B, angle between B and los]
		      const Vector&    f_grid,               // frequency vector                 [Hz]
		      const Index&     f_index,              // frequency grid point index
		      const Index&     zeeman_o2_onoff, 
		      const Numeric&   zeeman_o2_pressure_limit,
		      const Index&     ppath_index,          // index of propagation path index
		      const Numeric&   rte_pressure,         // total atmospheric pressure       [Pa]
		      const Numeric&   rte_temperature,      // temperature                      [K]
		      const Vector&    rte_vmr_list,         // list of species vmr              [1]
		      const ArrayOfIndex& species_index,     // index of key species   
		      const String&    abs_model,            // model selection string
		      const Vector&    abs_user_parameters,  // scaling factor(s)
		      const Index&     stokes_dim            // dimension of Stokes vector
		     )
{



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
  vmro2 = rte_vmr_list[ species_index[SPECIES_INDEX_O2] ];

  Numeric vmrh2o;
  vmrh2o =  rte_vmr_list[ species_index[SPECIES_INDEX_H2O] ];

  // set the lower pressure level where Zeeman splitting is effective.
  // An average altitude of 30km (see Rosenkranz F77 code description) is assumed.
  // This yields a pressure level of approximately 13 hPa (MLS atmosphere).
  bool zeeman = false;  // 1=true, 0=false
  if (rte_pressure <= 1.000e3)
    zeeman = true;

  // resize the output variables
  ext_mat_zee.resize(f_grid.nelem(), stokes_dim, stokes_dim);
  abs_vec_zee.resize(f_grid.nelem(), stokes_dim);
  

  // f_index is -1 if all frequencies should be calculated
  // the other case f_index != -1 has to be implemented later
  // number of elements in frequency grid
  Index n_f=0;
  if(f_index == -1)
     n_f = f_grid.nelem();

   // Loop pressure/temperature grid points:
  for ( Index i=0; i<n_f; ++i )
    {

      // Numeric geomag_strength = 3.50000e-5;
      // Numeric geomag_angle    = (90.0000000e0/57.29579e0);
      Numeric geomag_strength = geomag_los(ppath_index,0);
      Numeric geomag_angle    = geomag_los(ppath_index,1);


      // set the frequency where the calculation is performed
      Numeric f_grid_point = f_grid[i]; // [Hz]


      // define output variables extinction matric and absorption vector
      Matrix ext_mat_tmp;
      Vector abs_vec_tmp;
      ext_mat_tmp.resize(4,4);
      abs_vec_tmp.resize(4);


      // define flag ok:  1=true, 0=false
      Index ok=1;


      // call the absorption function
      ok = absPWRO2Model( ext_mat_tmp,
			  abs_vec_tmp, 
			  geomag_strength,
			  geomag_angle,
			  zeeman_o2_onoff,
                          zeeman_o2_pressure_limit,
			  f_grid_point,
			  rte_pressure,
			  rte_temperature, 
			  vmro2,
			  vmrh2o,
			  abs_model, 
			  abs_user_parameters);
      
      // output variable abs/ext (vector/matrix) should be passed 
      // to the agenda for gas absorption
      if (ok == 1) 
 	{
	  for ( Index k=0; k<stokes_dim; ++k ) 
 	    {
	      abs_vec_zee(i,k) += abs_vec_tmp[k];
	      for ( Index l=0; l<stokes_dim; ++l ) 
		ext_mat_zee(i,k,l) += ext_mat_tmp(k,l);
	    }
	}
      else
	{
	  cout << "absO2ZeemanModel: ERROR in O2 calculation!!\n";
	};
    };


  return;
}



// #################################################################################


void ZeemanO2Settings( // Output
		       Index&         zeeman_o2_onoff, 
		       Numeric&       zeeman_o2_pressure_limit,
		       // Keywords
		       const Index&   ZeemanO2OnOff,
                       const Numeric& ZeemanO2PressureLimit )
{

  // set Zeeman O2 on / off
  zeeman_o2_onoff          = ZeemanO2OnOff;

  // set upper pressure level for Zeeman
  zeeman_o2_pressure_limit = ZeemanO2PressureLimit;

  return;
};




// #################################################################################



//! test_zeeman
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Thomas Kuhn
   \date   2003-12-08
*/
void test_zeeman(
		 const Agenda& opt_prop_gas_agenda
		 )
{
  opt_prop_gas_agenda.execute();
}
