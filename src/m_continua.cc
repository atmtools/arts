// #################################################################################
/* Copyright (C) 2001,2002,2003 
   Thomas Kuhn    <tkuhn@uni-bremen.de>
   Stefan Buehler <sbuehler@uni-bremen.de>

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

// FIXME: Port doxygen file header and documentation from arts-1-0 version to here. 

#include "matpackI.h"
#include "mystring.h"
#include "exceptions.h"
#include "messages.h"
#include "continua.h"

// #################################################################################

// global constants as defined in constants.cc

extern const Numeric EULER_NUMBER;
extern const Numeric LOG10_EULER_NUMBER;
extern const Numeric NAT_LOG_TEN;
extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric BOLTZMAN_CONST;


// numerical constants specific defined for the file continua.cc

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

// MPM pre-factor for unit setting:
const Numeric dB_m_Hz   = 0.1820427855916028e-06; // [dB/m/Hz]   (4 * pi / c) * 10 * log(e)
const Numeric dB_km_GHz = 0.1820427855916028e+06; // [dB/km/GHz] (4 * pi / c) * 10 * log(e)


// absorption unit conversions

// conversion from dB/km to Np/km for absorption units:
const Numeric dB_km_to_Np_km = dB_to_Np;
// conversion from dB/km to Np/m for absorption units:
const Numeric dB_km_to_Np_m  = (1.00000e-3 / (10.0 * LOG10_EULER_NUMBER));
// conversion from dB/km to 1/m for absorption units:
const Numeric dB_km_to_1_m   = (1.00000e-3 / (10.0 * LOG10_EULER_NUMBER));


// lower limit for absorption calculation due to underflow error:

const Numeric VMRCalcLimit = 1.000e-25;

// #################################################################################


//---------------------------------------------------
//! Calculate water vapor absorption.
/*! 
  See arts -d absMPM02_H2O for detailed documentation.

  Allowed options for abs_model:<br>
  "MPM02"          - Calculate lines and continuum.<br>
  "MPM02Lines"     - Calculate only lines.<br>
  "MPM02Continuum" - Calculate only continuum.<br>
  "user"           - Use parameters given by abs_user_parameters,
		     instead of the predefined settings.

  Meaning of abs_user_parameters:<br>
  Only used if abs_model=="user". In that case, abs_user_parameters must have 3 elements:<br>
  1. Continuum scaling factor<br>
  2. Line strength scaling factor<br>
  3. Line broadening scaling factor<br>
  Setting all scaling factors to 1 gives the same behavior as abs_model=="MPM02".

  \retval abs Absorption coefficients [1/m], dimension: [ f_grid, abs_p (=abs_t) ].

  \param f_grid Frequency grid [Hz].
  \param abs_p  List of pressures [Pa].
  \param abs_t  List of temperatures [K]. Must have same length as abs_p!
  \param abs_vmr List of volume mixing ratios [absolute number]. Must have same length as abs_p!
  \param abs_model String specifying the model to use. 
  \param abs_user_parameters Parameters for abs_model="user".

  \author Thomas Kuhn, Stefan Buehler (new interface)
  \date   2002-05-06,  2003-11-16
*/
void absMPM02_H2O(// WS Output:
                  Matrix& abs,
                  // WS Input:
                  const Vector& f_grid,
                  const Vector& abs_p,
                  const Vector& abs_t,
                  const Vector& abs_vmr,
                  const String& abs_model,
                  const Vector& abs_user_parameters)
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21)
  const Numeric CC_MPM02 = 1.00000;
  const Numeric CL_MPM02 = 1.00000;
  const Numeric CW_MPM02 = 1.00000;
  // ---------------------------------------------------------------------------------------

  // Vector abs_user_parameters must generally be empty. But if abs_model is
  // "user", then it must have exactly 3 elements.
  if ( abs_model == "user" )
    {
      if ( abs_user_parameters.nelem() != 3 )
        throw runtime_error("Vector abs_user_parameters must have 3 elements.");
    }
  else
    {
      if ( abs_user_parameters.nelem() != 0 )
        throw runtime_error("Vector abs_user_parameters must be empty.");
    }   

  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  // number of lines of Liebe line catalog (0-33 lines, 34 cont. pseudo line)
  Index i_first = 0;
  Index i_last  = 34;
  if ( abs_model == "MPM02" )
    {
      CC      = CC_MPM02;
      CL      = CL_MPM02;
      CW      = CW_MPM02;
      i_first = 0;
      i_last  = 34;
    }
  else if ( abs_model == "MPM02Lines" )
    {
      CC      = 0.000;
      CL      = CL_MPM02;
      CW      = CW_MPM02;
      i_first = 0;
      i_last  = 33;
    }
  else if ( abs_model == "MPM02Continuum" )
    {
      CC      = CC_MPM02;
      CL      = 0.000;
      CW      = 0.000;
      i_first = 34;
      i_last  = 34;
    }
  else if ( abs_model == "user" )
    {
      CC      = abs_user_parameters[1];
      CL      = abs_user_parameters[2];
      CW      = abs_user_parameters[3];
      i_first = 0;
      i_last  = 34;
    }
  else
    {
      ostringstream os;
      os << "Wrong abs_model string given.\n"
	 << "Valid abs_models are: 'MPM02', 'MPM02Lines', 'MPM02Continuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "  H2O-MPM02: (abs_model=" << abs_model << ") parameter values in use:\n" 
	<< "  CC = " << CC << "\n"
	<< "  CL = " << CL << "\n"
	<< "  CW = " << CW << "\n";
  
  
  const Index n_p = abs_p.nelem();	// Number of pressure levels
  const Index n_f = f_grid.nelem();	// Number of frequencies

  // Check that dimensions of abs_p, abs_t, and abs_vmr agree:
  if ( n_p != abs_t.nelem() ||
       n_p != abs_vmr.nelem() )
    {
      ostringstream os;
      os << "The vectors abs_p, abs_t, and abs_vmr must all have the same length!\n"
	 << "Actual lengths: "
         << n_p << ", "
         << abs_t.nelem() << ", "
         << abs_vmr.nelem() << ".";
      throw runtime_error(os.str());
    }

  // Make abs the right dimension and inititalize to zero:
  abs.resize(n_f,n_p);

  // FIXME: Now should come the real work...
  
  // Assign a dummy value to abs for testing:
  abs = 99;

  return;
}


// =================================================================================




/*! 
  See arts -d LagrangeInterpol4 for detailed documentation.

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

  \author Thomas Kuhn, Stefan Buehler (new interface)
  \date   2003-11-25
*/

Numeric LagrangeInterpol4( const Vector& x,
			   const Vector& y,
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
	  if ( (k != i) && (abs(x[i]-x[k]) > Dlimit) )  
	    b[i] = b[i] * ( (a-x[k]) / (x[i]-x[k]) );
	};
    };

  Numeric ya = 0.000e0;
  for (Index i=0 ; i < n_x ; ++i) ya = ya + b[i]*y[i];

  return ya;
}


// =================================================================================


/*! 
  This function calculates the CKD 4-point interpolation from the CKD 
  internal frequency grid to a point on the user defined frequency grid.

  In all the CKD routines in which the interpolation is important 
  a four point interpolation scheme with continuous derivatives at the 
  interior two points. This means that for interpolations between X_2 
  and X_3 with the value set Y_1, Y_2, Y_3, and Y_4, the four term 
  polynomial goes through values Y_2 and Y_3 and satisfies the derivative 
  constraint {(Y_3-Y_1) / (X_3-X_1)} at X_2 and {(Y_4-Y_2) / (X_4-X_2)} 
  at X_3. This should be marginally better than the standard for-point 
  Lagrange interpolation scheme.

  \retval xint   interpolated absorption at user defined frequency grid point [1/cm]

  \param V1A     first frequency on CKD internal frequency grid        [cm-1]
  \param DVA     frequency spacing on CKD internal frequency grid      [cm-1]
  \param A       absorption values at the CKD internal frequency grid  [1/cm]
  \param VI      point on the user defined frequency grid              [cm-1]

  \author Thomas Kuhn
  \date   2003-11-27
*/

Numeric XINT_FUN( const Numeric V1A,
                  const Numeric DVA,
		  const Vector& A,
                  const Numeric VI)
{

  // ----------------------------------------------------------------------
  //     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED               
  //     FROM V1A TO V2A IN INCREMENTS OF DVA INTO XINT
  // ----------------------------------------------------------------------

  /*
  */

  const Numeric ONEPL  = 1.001;

  Numeric RECDVA = 1.00e0/DVA;
    
  int J      = (int) ((VI-V1A)*RECDVA + ONEPL) ; 
  Numeric VJ = V1A + DVA * (Numeric)(J-1);    
  Numeric P  = RECDVA * (VI-VJ);        
  Numeric C  = (3.000e0-2.000e0*P) * P * P;         
  Numeric B  = 0.500e0 * P * (1.000e0-P);          
  Numeric B1 = B * (1.000e0-P);              
  Numeric B2 = B * P;                   
  
  Numeric xint = 0.0000e0;
  xint = -A[J-1] * B1             +
          A[J]   * (1.00e0-C+B2)  + 
          A[J+1] * (C+B1)         - 
          A[J+2] * B2;
  
  return xint;
}


// =================================================================================


/*! 
  This function calculates the CKD radiation field term according to
  RADFN = VI * tanh(XKT*VI/2)
  for special ranges special approximations are made. It is assured that 
  RADFN = 0.00 for VI=0.00. 

  \retval RADFN  radiation field term at input frequency VI       [cm-1]

  \param VI      point on the user defined frequency grid         [cm-1]
  \param XKT     XKT = (T*k_B) / (h*c)                            [cm-1]

  \author Thomas Kuhn
  \date   2003-11-27
*/

Numeric RADFN_FUN (const Numeric VI, 
		   const Numeric XKT)
{
// ---------------------------------------------------------------------- B18060
//              LAST MODIFICATION:    12 AUGUST 1991                      B17940
//                                                                        B17950
//                 IMPLEMENTATION:    R.D. WORSHAM                        B17960
//                                                                        B17970
//            ALGORITHM REVISIONS:    S.A. CLOUGH                         B17980
//                                    R.D. WORSHAM                        B17990
//                                    J.L. MONCET                         B18000
//                                                                        B18010
//                                                                        B18020
//                    ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC.         B18030
//                    840 MEMORIAL DRIVE,  CAMBRIDGE, MA   02139          B18040
//                                                                        B18050
//                                                                        B18070
//              WORK SUPPORTED BY:    THE ARM PROGRAM                     B18080
//                                    OFFICE OF ENERGY RESEARCH           B18090
//                                    DEPARTMENT OF ENERGY                B18100
//                                                                        B18110
//                                                                        B18120
//     SOURCE OF ORIGINAL ROUTINE:    AFGL LINE-BY-LINE MODEL             B18130
//                                                                        B18140
//                                            FASCOD3                     B18150
//                                                                        B18160
// ---------------------------------------------------------------------- B18060
//                                                                        B18170
//  IN THE SMALL XVIOKT REGION 0.5 IS REQUIRED

  Numeric XVI   = VI;                            
  Numeric RADFN = 0.00e0;                                
									   
  if (XKT > 0.000e0)
    {                                    
      Numeric XVIOKT = XVI/XKT;                 
      
      if (XVIOKT <= 0.01e0)
	{
	  RADFN = 0.500e0 * XVIOKT * XVI;
	}
      else if (XVIOKT <= 10.0e0)
	{
	  Numeric EXPVKT = exp(-XVIOKT);
	  RADFN = XVI * (1.00e0-EXPVKT) / (1.00e0+EXPVKT);
	}
      else
	{
	  RADFN = XVI;                         
	}                                  
    }
  else
    {
      RADFN = XVI;
    }
  
  return RADFN;
}


// =================================================================================


// CKD version MT 1.00 H2O self continuum absorption model
/**

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

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue, 
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author   Thomas Kuhn
   \date     2003-11-21
*/ 
void abs_CKDMT_H2O_H2O(// WS Output:
		       Matrix&            abs,
		       // WS Input:
		       const Vector&      f_grid,
		       const Vector&      abs_p,
		       const Vector&      abs_t,
		       const Vector&      abs_vmr,
		       const String&      abs_model,
		       const Vector&      abs_user_parameters)
{

  // check the model name about consistency
  if ((abs_model != "user") &&  (abs_model != "CKDMT100"))
    {
      ostringstream os;
      os << "ERROR in method abs_CKDMT_H2O_H2O:\n"  
	 << "CKD_MT1.00 H2O self continuum detected wrong model name >>" << abs_model << "<<.\n"
	 << "VALID model names are >>user<< and >>CKDMT100<<\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the self H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( abs_model == "user" )
    {
      if ( abs_user_parameters.nelem() != 1 )
        throw runtime_error("Vector abs_user_parameters must have 1 element.");
      // input scaling factor of calculated absorption
      ScalingFac = abs_user_parameters[0];
    }
  if ( abs_model == "CKDMT100" )
    {
      if ( abs_user_parameters.nelem() != 0 )
	throw runtime_error("Vector abs_user_parameters must be empty.");
      // input scaling factor of calculated absorption
      ScalingFac = 1.0000e0;
    }

  const Index n_p = abs_p.nelem();	// Number of pressure levels
  const Index n_f = f_grid.nelem();	// Number of frequencies


  // Check that dimensions of abs_p, abs_t, and abs_vmr agree:
  if ( n_p != abs_t.nelem() ||
       n_p != abs_vmr.nelem() )
    {
      ostringstream os;
      os << "The vectors abs_p, abs_t, and abs_vmr must all have the same length!\n"
	 << "Actual lengths: "
         << n_p << ", "
         << abs_t.nelem() << ", "
         << abs_vmr.nelem() << ".";
      throw runtime_error(os.str());
    }

  // Make abs the right dimension and inititalize to zero:
  abs.resize(n_f,n_p);


  // empirical correction factors
  const Numeric XFACREV[15] = 
    {1.003, 1.009, 1.015, 1.023, 1.029,1.033,
     1.037, 1.039, 1.040, 1.046, 1.036,1.027,
     1.01,  1.002, 1.00};

  // wavenumber range where CKD H2O self continuum is valid
  const Numeric VABS_min = -2.000e1; // [cm^-1]
  const Numeric VABS_max =  2.000e4; // [cm^-1]

  // It is assumed here that f_grid is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_grid. n_f_new < n_f
  Numeric V1ABS = f_grid[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_grid[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD_MT 1.00 H2O self continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << SL296_ckd_mt_100_v1 << "<->" << SL296_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ------------------- subroutine SL296/SL260 ----------------------------

  if (SL296_ckd_mt_100_v1 != SL260_ckd_mt_100_v1)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 H2O self continuum:\n"
	 << "parameter V1 not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (SL296_ckd_mt_100_v2 != SL260_ckd_mt_100_v2)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 H2O self continuum:\n"
	 << "parameter V2 not the same for different ref. temperatures.\n";
	throw runtime_error(os.str());
    }
  if (SL296_ckd_mt_100_dv != SL260_ckd_mt_100_dv)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 H2O self continuum:\n"
	 << "parameter DV not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (SL296_ckd_mt_100_npt != SL260_ckd_mt_100_npt)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 H2O self continuum:\n"
	 << "parameter NPT not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  
  // retrieve the appropriate array sequence of the self continuum
  // arrays of the CKD model.
  Numeric DVC = SL296_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-SL296_ckd_mt_100_v1) / SL296_ckd_mt_100_dv);
  if (V1C < SL296_ckd_mt_100_v1) I1 = -1;
  V1C = SL296_ckd_mt_100_v1 + (SL296_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-SL296_ckd_mt_100_v1) / SL296_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > SL296_ckd_mt_100_npt) NPTC = SL296_ckd_mt_100_npt+1;

  V2C = V1C + SL296_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      ostringstream os;
      out3 << "WARNING:\n"  
	   << "  CKD_MT 1.00 H2O self continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  // continuum coefficient vectors
  Vector SH2OT0; // [cm^3/molecules]
  Vector SH2OT1; // [cm^3/molecules]
  SH2OT0.resize(NPTC+addF77fields);
  SH2OT1.resize(NPTC+addF77fields);

  // loop over all CKD internal frequency grid points selected for 
  // the user defined frequency range
  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > SL296_ckd_mt_100_npt) ) 
	{
	  SH2OT0[J] = 0.0e0;   // at T=296 K
	  SH2OT1[J] = 0.0e0;   // at T=260 K
	}
      else
	{
	  SH2OT0[J] = SL296_ckd_mt_100[I];    // at T=296 K
	  SH2OT1[J] = SL260_ckd_mt_100[I];    // at T=260 K 
	}
    }

  // ------------------- subroutine SL296/SL260 ----------------------------
  
  // Loop over user defined pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {

      // XKT = (T*k_B) / (h*c) [1/cm]
      Numeric XKT  = abs_t[i] / 1.4387752e0;
      // atmospheric temperature factor for self cont. inter/extrapolation
      Numeric Tfac = (abs_t[i]-296.000e0) / (260.000e0-296.000e0); // [1]
      // ((CKD molecular cross section) * (number density) * (radiation field)) [1/cm]
      Vector k;
      k.resize(NPTC+addF77fields);
      k[0] = 0.00e0; // not used array field

      // loop over the CKD internal frequency grid
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  k[J]         = 0.00e0;
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SH2O = 0.0e0;
	  Numeric SFAC = 0.00e0;  
	  if (SH2OT0[J] > 0.0e0)
	    {
	      SH2O = SH2OT0[J] * pow( (SH2OT1[J]/SH2OT0[J]), Tfac ); 
	      SFAC = 1.00e0;

	      if ( (VJ >= 820.0e0) && (VJ <= 960.0e0) )
		{
		  int JFAC = (int)((VJ - 820.0e0)/10.0e0 + 0.00001e0);
		  if ( (JFAC >= 0) && (JFAC <=14) )
		    SFAC = XFACREV[JFAC];
		}

	      SH2O = SFAC * SH2O;
	    }

	  k[J] =  RADFN_FUN(VJ,XKT)                         *
                  (abs_vmr[i] * abs_p[i] / abs_t[i])        * 
                  (abs_vmr[i] * abs_p[i] / abs_t[i])        * 
	          (2.9600e2   / 1.0130e5 / BOLTZMAN_CONST)  * 
	          1.000e-6                                  * 
                  (SH2O*1.000e-20);
	};

      // Loop over user defined frequency grid.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_grid[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V >= 0.000e0) && (V < SL296_ckd_mt_100_v2) )
	    {
	      // arts absorption coefficient [1/m]
	      // improved Lagrange polynomial interpolation scheme used to
	      // interpolate k from the DKC internal frequency grid on the 
	      // user defined frequency grid.
	      // The factor 100 comes from the conversion from 1/cm to 1/m for
	      // the absorption coefficient.
	      abs(s,i) = ScalingFac * 1.000e2 * XINT_FUN(V1C,DVC,k,V);
	    };
	};
    };
  
}

// =================================================================================

// CKD version MT 1.00 H2O foreign continuum absorption model
/**

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

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue, 
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author   Thomas Kuhn
   \date     2003-11-21
*/ 
void abs_CKDMT_H2O_AIR(// WS Output:
		       Matrix&           abs,
		       // WS Input:
		       const Vector&     f_grid,
		       const Vector&     abs_p,
		       const Vector&     abs_t,
		       const Vector&     abs_vmr,
		       const String&     abs_model,
		       const Vector&     abs_user_parameters)
{


  // check the model name about consistency
  if ((abs_model != "user") &&  (abs_model != "CKDMT100"))
    {
      ostringstream os;
      os << "ERROR in method abs_CKDMT_H2O_H2O:\n"  
	 << "CKD_MT1.00 H2O self continuum detected wrong model name >>" << abs_model << "<<.\n"
	 << "VALID model names are >>user<< and >>CKDMT100<<\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the foreign H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( abs_model == "user" )
    {
      if ( abs_user_parameters.nelem() != 1 )
	throw runtime_error("Vector abs_user_parameters must have 1 element.");
      // input scaling factor of calculated absorption
      ScalingFac = abs_user_parameters[0];
    }
  if ( abs_model == "CKDMT100" )
    {
      if ( abs_user_parameters.nelem() != 0 )
	throw runtime_error("Vector abs_user_parameters must be empty.");
      // input scaling factor of calculated absorption
      ScalingFac = 1.0000e0;
    }

  const Index n_p = abs_p.nelem();	// Number of pressure levels
  const Index n_f = f_grid.nelem();	// Number of frequencies

  // Check that dimensions of abs_p, abs_t, and abs_vmr agree:
  if ( n_p != abs_t.nelem() ||
       n_p != abs_vmr.nelem() )
    {
      ostringstream os;
      os << "The vectors abs_p, abs_t, and abs_vmr must all have the same length!\n"
	 << "Actual lengths: "
         << n_p << ", "
         << abs_t.nelem() << ", "
         << abs_vmr.nelem() << ".";
      throw runtime_error(os.str());
    }


  // Make abs the right dimension and inititalize to zero:
  abs.resize(n_f,n_p);


  // wavenumber range where CKD H2O foreign continuum is valid
  const Numeric VABS_min = -2.000e1; // [cm^-1]
  const Numeric VABS_max =  2.000e4; // [cm^-1]


  // It is assumed here that f_grid is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_grid. n_f_new < n_f
  Numeric V1ABS = f_grid[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_grid[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD_MT 1.00 H2O foreign continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << FH2O_ckd_mt_100_v1 << "<->" << FH2O_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ---------------------- subroutine FRN296 ------------------------------

  // retrieve the appropriate array sequence of the foreign continuum
  // arrays of the CKD model.
  Numeric DVC = FH2O_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-FH2O_ckd_mt_100_v1) / FH2O_ckd_mt_100_dv);
  if (V1C < FH2O_ckd_mt_100_v1) I1 = -1;
  V1C = FH2O_ckd_mt_100_v1 + (FH2O_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-FH2O_ckd_mt_100_v1) / FH2O_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > FH2O_ckd_mt_100_npt) NPTC = FH2O_ckd_mt_100_npt+1;

  V2C = V1C + FH2O_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n" 
	   << "  CKD_MT 1.00 H2O foreign continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  // continuum coefficient vector
  Vector FH2OT0; // [cm^3/molecules]
  FH2OT0.resize(NPTC+addF77fields);

  // loop over all CKD internal frequency grid points selected for 
  // the user defined frequency range
  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > FH2O_ckd_mt_100_npt) ) 
	{
	  FH2OT0[J] = 0.0e0;
	}
      else
	{
	  FH2OT0[J] = FH2O_ckd_mt_100[I];
	};
    };

  // ---------------------- subroutine FRN296 ------------------------------
  
  


  // Loop over user defined pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      // XKT = (T*k_B) / (h*c) [1/cm]
      Numeric XKT = abs_t[i] / 1.4387752e0;
      
      // ((CKD molecular cross section) * (number density) * (radiation field)) [1/cm]
      Vector k;
      k.resize(NPTC+addF77fields); // [1/cm]
      k[0] = 0.00e0; // not used array field

      // loop over the CKD internal frequency grid
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric FH2O = FH2OT0[J];
	  k[J]         =  RADFN_FUN(VJ,XKT)                          *
                       (abs_vmr[i] * abs_p[i] / abs_t[i])            * 
	               ((1.0000e0-abs_vmr[i]) * abs_p[i] / abs_t[i]) * 
	               (2.9600e2   / 1.0130e5 / BOLTZMAN_CONST)      * 
	               1.000e-6                                      * 
                       (FH2O*1.000e-20);
	};
      
      // Loop over user defined frequency grid.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_grid[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V >= 0.000e0) && (V < VABS_max) )
	    {
	      // arts absorption coefficient [1/m]
	      // improved Lagrange polynomial interpolation scheme used to
	      // interpolate k from the DKC internal frequency grid on the 
	      // user defined frequency grid.
	      // The factor 100 comes from the conversion from 1/cm to 1/m for
	      // the absorption coefficient.
	      abs(s,i) = ScalingFac * 1.000e2 * XINT_FUN(V1C,DVC,k,V);
	    };
	};
    };
}
