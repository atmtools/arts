/* Copyright (C) 2001 Thomas Kuhn    <tkuhn@uni-bremen.de>
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

/**
   \file   continua.cc


   \param xsec Output:  Absorption cross section, defined such that the
                 absorption coefficient alpha is:<br>
                 alpha [1/m] = xsec * VMR.<br>
                 The functions adds to xsec, rather than replacing the
                 previous content. <br>
   <br>
   <br>
   <H3>The following full water vapor models are implemented:</H3><br>
   <ol>
   <li><b>H2O-MPM87</b> absorption model (line and continuum) according to <br>
          H. J. Liebe,<br> 
          <i>A contribution to modeling atmospheric millimeter-wave properties</i>,<br>
          Frequenz,  41, 1987, 31-36<br>
          and<br>
          H. J. Liebe and D. H. Layton,<br>
          <i>Millimeter-wave properties of the atmosphere:
          Laboratory studies and propagation modeling</i>,<br>
          U.S. Dept. of Commerce, National Telecommunications and Information
          Administration, Institute for Communication Sciences,<br>
          325 Broadway, Boulder, CO 80303-3328, report 87224.
   </li>
   <li><b>H2O-MPM89</b> absorption model (line and continuum) according to <br>
          H. J. Liebe,<br> Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631.
   </li>
   <li><b>H2O-MPM93</b> absorption model (line and continuum) according to <br>
          H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
          <i>Propagation modeling of moist air and suspended water/ice
          particles at frequencies below 1000 GHz</i>,<br>
          AGARD 52nd Specialists Meeting of the Electromagnetic Wave
          Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21 
          <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   <li><b>H2O-CP98</b> absorption model (line and continuum) according to <br>
          S. L. Cruz-Pol et al.,<br> Radio Science, 33(5), 1319, 1998
          <a href="http://ece.uprm.edu/~pol/Atmosphere.html">(WWW access)</a>.
   </li>
   <li><b>H2O-PWR98</b> absorption model (line and continuum) according to <br>
          P. W. Rosenkranz,<br> 
          Radio Science, 33(4),  919, 1998 and<br>
          Radio Science, 34(4), 1025, 1999 
          <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   </ol>
   <br>
   <br>
   <br>
   <br>
   <H3>The following full oxygen models are implemented:</H3><br>
   <ol>
   <li><b>O2-MPM93</b> absorption model (line and continuum) according to <br>
          H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
          <i>Propagation modeling of moist air and suspended water/ice
          particles at frequencies below 1000 GHz</i>,<br>
          AGARD 52nd Specialists Meeting of the Electromagnetic Wave
          Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21 
          <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   <li><b>O2-PWR93</b> absorption model (line and continuum) according to <br>
          P. W. Rosenkranz,<br> Chapter 2, in M. A. Janssen, <br>
          <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
          John Wiley & Sons, Inc., 1993
          <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   </ol>
   <br>
   <br>
   <br>
   <br>
   <H3>The following continuum parameterizations are implemented:</H3><br>
   <ol>
   <li><b>H2O-H2O (H2O-SelfContStandardType)</b>:<br> 
         P. W. Rosenkranz,<br> 
         Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and <br>
         Radio Science, Vol. 34, No 4, Page 1025, 1999 
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   <li><b>H2O-air (H2O-ForeignContStandardType)</b>: <br>
         P. W. Rosenkranz,<br> 
         Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and <br>
         Radio Science, Vol. 34, No 4, Page 1025, 1999
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   <li><b>H2O-air (H2O-ContMPM93)</b>:<br> 
         H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
         <i>Propagation modeling of moist air and suspended water/ice<br>
         particles at frequencies below 1000 GHz</i>,<br>
         AGARD 52nd Specialists Meeting of the Electromagnetic Wave<br>
         Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
         <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   <li><b>O2-air (O2-SelfContStandardType)</b>:<br> 
         P. W. Rosenkranz,<br> 
         Chapter 2, in M. A. Janssen, <br>
         <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
         John Wiley & Sons, Inc., 1993
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.<br>
         and also described in<br>
         H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
         <i>Propagation modeling of moist air and suspended water/ice<br>
         particles at frequencies below 1000 GHz</i>,<br>
         AGARD 52nd Specialists Meeting of the Electromagnetic Wave<br>
         Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
         <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   <li><b>N2-N2 (N2-SelfContStandardType)</b>:<br> 
         P. W. Rosenkranz,<br>
         Chapter 2, in M. A. Janssen, <br>
         <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
         John Wiley & Sons, Inc., 1993
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   <li><b>N2-N2 (N2-SelfContMPM93)</b>:<br> 
         H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
         <i>Propagation modeling of moist air and suspended water/ice<br>
         particles at frequencies below 1000 GHz</i>,<br>
         AGARD 52nd Specialists Meeting of the Electromagnetic Wave<br>
         Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
         <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   <li><b>CO2-CO2 (CO2-SelfContPWR93)</b>:<br> 
         P. W. Rosenkranz,<br>
         Chapter 2, in M. A. Janssen, <br>
         <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
         John Wiley & Sons, Inc., 1993 
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   <li><b>CO2-N2 (CO2-ForeignContPWR93)</b>: <br>
         P. W. Rosenkranz,<br> 
         Chapter 2, in M. A. Janssen, <br>
         <i>Atmospheric Remote Sensing by Microwave Radiometry</i>,
         John Wiley & Sons, Inc., 1993 
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   </ol>
   <br>
   <br>
   <br>
   <br>
   <H3>The following cloud absorption models are implemented:</H3><br>
   <ol>
   <li><b>Suspended water droplet (liquidcloud-MPM93)</b> 
      absorption parameterization from MPM93 model<br>
      H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
      <i>Propagation modeling of moist air and suspended water/ice
      particles at frequencies below 1000 GHz</i>,<br>
      AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21 
      <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   <li><b>Ice crystal absorption (icecloud-MPM93)</b> parameterization from MPM93 model<br>
      H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
      <i>Propagation modeling of moist air and suspended water/ice
      particles at frequencies below 1000 GHz</i>,<br>
      AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21 
      <a href="ftp://ftp.its.bldrdoc.gov/pub/mpm93/">(WWW access)</a>.
   </li>
   </ol>
<br>
<br>

   <b>The following unit conversions are used for the implemented models:</b><br>
   (SI units: meter, kilogram, second, ampere, Kelvin, candela)<br>

   \verbatim

   x g/cm³ = y kg/m³    <===>    y = x * 1.00e3
   x g/m³  = y kg/m³    <===>    y = x * 1.00e-3
   x GHz   = y Hz       <===>    y = x * 1.00e9
   x 1/GHz = y 1/Hz     <===>    y = x * 1.00e-9
   x hPa   = y Pa       <===>    y = x * 1.00e2
   x 1/hPa = y 1/Pa     <===>    y = x * 1.00e-2
   x 1/cm  = y 1/m      <===>    y = x * 1.0e2
   x 1/km  = y 1/m      <===>    y = x * 1.00e-3
   x dB    = y Np       <===>    y = x / [10.0 * log10(e)]
   x dB/km = y 1/m      <===>    y = x * 1.00e-3 / [10.0 * log10(e)]
   x Np/km = y 1/m      <===>    y = x * 1.00e-3

   and especially for the MPM model versions:

   (4 * pi / c) * 10 * log(e) = 0.1820 * 10^6  dB/km/GHz
                              = 0.1820 * 10^-6 dB/m/Hz
   \endverbatim

<br>

   \author Thomas Kuhn
   \date   2001-11-01
*/

#include <cmath>
#include "arts.h"
#include "matpackI.h"
#include "array.h"
#include "absorption.h"
#include "messages.h"
#include "continua.h"
#include "math_funcs.h"

// #################################################################################

// global constants as defined in constants.cc

extern const Numeric EULER_NUMBER;
extern const Numeric LOG10_EULER_NUMBER;
extern const Numeric NAT_LOG_TEN;
extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;

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
// ############################## WATER VAPOR MODELS ###############################
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            H2O (lines+continuum) according to MPM87 [1/m]
   \param    CCin           scaling factor for the H2O-continuum  [1]
   \param    CLin           scaling factor for the H2O-line strengths [1]
   \param    CWin           scaling factor for the H2O-line widths    [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, and CWin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for model 'user' the input parameters CCin, CLin, and CWin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM87', 'MPM87Lines', 'MPM87Continuum', and 'user'.
             See the user guide for detailed explanations.

   \remark   H. J. Liebe,<br> 
             <i>A contribution to modeling atmospheric millimeter-wave properties</i>,<br>
             Frequenz,  41, 1987, 31-36<br>
             and<br>
             H. J. Liebe and D. H. Layton,<br>
             <i>Millimeter-wave properties of the atmosphere:
             Laboratory studies and propagation modeling</i>,<br>
             U.S. Dept. of Commerce, National Telecommunications and Information
             Administration, Institute for Communication Sciences,<br>
             325 Broadway, Boulder, CO 80303-3328, report 87224.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM87H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,
                       ConstVectorView   f_mono,
                       ConstVectorView   p_abs,
                       ConstVectorView   t_abs,
                       ConstVectorView   vmr )
{
  //
  // Coefficients are from Liebe, Radio Science, 20(5), 1985, 1069
  //         0           1        2       3      
  //         f0          b1       b2      b3     
  //        [GHz]     [kHz/kPa]   [1]   [GHz/kPa]
  const Numeric mpm87[30][4] = { 
    {    22.235080,    0.1090,  2.143,   27.84e-3},
    {    67.813960,    0.0011,  8.730,   27.60e-3},
    {   119.995940,    0.0007,  8.347,   27.00e-3},
    {   183.310117,    2.3000,  0.653,   31.64e-3},
    {   321.225644,    0.0464,  6.156,   21.40e-3},
    {   325.152919,    1.5400,  1.515,   29.70e-3},
    {   336.187000,    0.0010,  9.802,   26.50e-3},
    {   380.197372,   11.9000,  1.018,   30.36e-3},
    {   390.134508,    0.0044,  7.318,   19.00e-3},
    {   437.346667,    0.0637,  5.015,   13.70e-3},
    {   439.150812,    0.9210,  3.561,   16.40e-3},
    {   443.018295,    0.1940,  5.015,   14.40e-3},
    {   448.001075,   10.6000,  1.370,   23.80e-3},
    {   470.888947,    0.3300,  3.561,   18.20e-3},
    {   474.689127,    1.2800,  2.342,   19.80e-3},
    {   488.491133,    0.2530,  2.814,   24.90e-3},
    {   503.568532,    0.0374,  6.693,   11.50e-3},
    {   504.482692,    0.0125,  6.693,   11.90e-3},
    {   556.936002,  510.0000,  0.114,   30.00e-3},
    {   620.700807,    5.0900,  2.150,   22.30e-3},
    {   658.006500,    0.2740,  7.767,   30.00e-3},
    {   752.033227,  250.0000,  0.336,   28.60e-3},
    {   841.073593,    0.0130,  8.113,   14.10e-3},
    {   859.865000,    0.1330,  7.989,   28.60e-3},
    {   899.407000,    0.0550,  7.845,   28.60e-3},
    {   902.555000,    0.0380,  8.360,   26.40e-3},
    {   906.205524,    0.1830,  5.039,   23.40e-3},
    {   916.171582,    8.5600,  1.369,   25.30e-3},
    {   970.315022,    9.1600,  1.842,   24.00e-3},
    {   987.926764,  138.0000,  0.178,   28.60e-3}};

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM87 model (Radio Science, 20(5), 1985, 1069):
  const Numeric CC_MPM87 = 1.00000;
  const Numeric CL_MPM87 = 1.00000;
  const Numeric CW_MPM87 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  if ( model == "MPM87" )
    {
      CC = CC_MPM87;
      CL = CL_MPM87;
      CW = CW_MPM87;
    }
  else if ( model == "MPM87Lines" )
    {
      CC = 0.000;
      CL = CL_MPM87;
      CW = CW_MPM87;
    }
  else if ( model == "MPM87Continuum" )
    {
      CC = CC_MPM87;
      CL = 0.000;
      CW = 0.000;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CL = CLin;
      CW = CWin;
    }
  else
    {
      ostringstream os;
      os << "H2O-MPM87: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM87', 'MPM87Lines', 'MPM87Continuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "H2O-MPM87: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n";
  
  
  // number of lines of liebe line catalog (30 lines)
  const Index i_first = 0;
  const Index i_last  = 29;
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in [hPa] therefore the factor 0.01)
  for ( Index i=0; i<n_p; ++i )
    {
      // here the total pressure is not multiplied by the H2O vmr for the 
      // P_H2O calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pwv_dummy = Pa_to_kPa * p_abs[i];
      // relative inverse temperature [1]
      Numeric theta = (300.0 / t_abs[i]);
      // H2O partial pressure [kPa]
      Numeric pwv   = Pa_to_kPa * p_abs[i] * vmr[i];
      // dry air partial pressure [kPa]
      Numeric pda   = (Pa_to_kPa * p_abs[i]) - pwv;
      // H2O continuum absorption [dB/km/GHz2] like in the original MPM87
      Numeric Nppc  = CC * pwv_dummy * pow(theta, 3.0) * 1.000e-5 *
                          ( (0.113 * pda) + (3.57 * pwv * pow(theta, 7.8)) );

      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
        {
          // input frequency in [GHz]
          Numeric ff   = f_mono[s] * Hz_to_GHz; 
          // H2O line contribution at position f
          Numeric Nppl = 0.000;
          
          // Loop over MPM89 H2O spectral lines
          for ( Index l = i_first; l <= i_last; ++l )
            {
              // line strength [kHz]
              Numeric strength = CL * pwv_dummy * mpm87[l][1] *
                                      pow(theta,3.5) * exp(mpm87[l][2]*(1.000-theta));
              // line broadening parameter [GHz]
              Numeric gam      = CW * mpm87[l][3] * 
                                     ( (4.80 * pwv * pow(theta, 1.1)) + 
                                       (       pda * pow(theta, 0.6)) );
              // effective line width with Doppler broadening [GHz]
              // gam              = sqrt(gam*gam + (2.14e-12 * mpm87[l][0] * mpm87[l][0] / theta));
              // H2O line absorption [dB/km/GHz] like in the original MPM87
              Nppl            += strength * MPMLineShapeFunction(gam, mpm87[l][0], ff); 
            }
          // xsec = abs/vmr [1/m] but MPM87 is in [dB/km] --> conversion necessary
          xsec(s,i)  += dB_km_to_1_m * 0.1820 * ff * ( Nppl + (Nppc * ff) );
        }
    }
  return;
}
//
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            H2O (lines+continuum) according to MPM89 [1/m]
   \param    CCin           scaling factor for the H2O-continuum  [1]
   \param    CLin           scaling factor for the line strengths [1]
   \param    CWin           scaling factor for the line widths    [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, and CWin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, and CWin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM89', 'MPM89Lines', 'MPM89Continuum', and 'user'.
             See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe, Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM89H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,     // model
                       ConstVectorView   f_mono,
                       ConstVectorView   p_abs,
                       ConstVectorView   t_abs,
                       ConstVectorView   vmr )
{
  //
  // Coefficients are from Liebe, Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631
  //         0           1        2       3        4      5      6
  //         f0          b1       b2      b3       b4     b5     b6
  //        [GHz]     [kHz/kPa]   [1]   [MHz/kPa]  [1]    [1]    [1]
  const Numeric mpm89[30][7] = { 
    {    22.235080,    0.1090,  2.143,   28.11,   0.69,  4.80,  1.00},
    {    67.813960,    0.0011,  8.735,   28.58,   0.69,  4.93,  0.82},
    {   119.995940,    0.0007,  8.356,   29.48,   0.70,  4.78,  0.79},
    {   183.310074,    2.3000,  0.668,   28.13,   0.64,  5.30,  0.85},
    {   321.225644,    0.0464,  6.181,   23.03,   0.67,  4.69,  0.54},
    {   325.152919,    1.5400,  1.540,   27.83,   0.68,  4.85,  0.74},
    {   336.187000,    0.0010,  9.829,   26.93,   0.69,  4.74,  0.61},
    {   380.197372,   11.9000,  1.048,   28.73,   0.69,  5.38,  0.84},
    {   390.134508,    0.0044,  7.350,   21.52,   0.63,  4.81,  0.55},
    {   437.346667,    0.0637,  5.050,   18.45,   0.60,  4.23,  0.48},
    {   439.150812,    0.9210,  3.596,   21.00,   0.63,  4.29,  0.52},
    {   443.018295,    0.1940,  5.050,   18.60,   0.60,  4.23,  0.50},
    {   448.001075,   10.6000,  1.405,   26.32,   0.66,  4.84,  0.67},
    {   470.888947,    0.3300,  3.599,   21.52,   0.66,  4.57,  0.65},
    {   474.689127,    1.2800,  2.381,   23.55,   0.65,  4.65,  0.64},
    {   488.491133,    0.2530,  2.853,   26.02,   0.69,  5.04,  0.72},
    {   503.568532,    0.0374,  6.733,   16.12,   0.61,  3.98,  0.43},
    {   504.482692,    0.0125,  6.733,   16.12,   0.61,  4.01,  0.45},
    {   556.936002,  510.0000,  0.159,   32.10,   0.69,  4.11,  1.00},
    {   620.700807,    5.0900,  2.200,   24.38,   0.71,  4.68,  0.68},
    {   658.006500,    0.2740,  7.820,   32.10,   0.69,  4.14,  1.00},
    {   752.033227,  250.0000,  0.396,   30.60,   0.68,  4.09,  0.84},
    {   841.073593,    0.0130,  8.180,   15.90,   0.33,  5.76,  0.45},
    {   859.865000,    0.1330,  7.989,   30.60,   0.68,  4.09,  0.84},
    {   899.407000,    0.0550,  7.917,   29.85,   0.68,  4.53,  0.90},
    {   902.555000,    0.0380,  8.432,   28.65,   0.70,  5.10,  0.95},
    {   906.205524,    0.1830,  5.111,   24.08,   0.70,  4.70,  0.53},
    {   916.171582,    8.5600,  1.442,   26.70,   0.70,  4.78,  0.78},
    {   970.315022,    9.1600,  1.920,   25.50,   0.64,  4.94,  0.67},
    {   987.926764,  138.0000,  0.258,   29.85,   0.68,  4.55,  0.90}};

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM89 model 
  // (Liebe, Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631):
  const Numeric CC_MPM89 = 1.00000;
  const Numeric CL_MPM89 = 1.00000;
  const Numeric CW_MPM89 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model goes for values!!):
  Numeric CC, CL, CW;
  if ( model == "MPM89" )
    {
      CC = CC_MPM89;
      CL = CL_MPM89;
      CW = CW_MPM89;
    }
  else if ( model == "MPM89Lines" )
    {
      CC = 0.000;
      CL = CL_MPM89;
      CW = CW_MPM89;
    }
  else if ( model == "MPM89Continuum" )
    {
      CC = CC_MPM89;
      CL = 0.000;
      CW = 0.000;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CL = CLin;
      CW = CWin;
    }
  else
    {
      ostringstream os;
      os << "H2O-MPM89: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM89', 'MPM89Lines', 'MPM89Continuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "H2O-MPM89: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n";
  
  
  // number of lines of Liebe line catalog (30 lines)
  const Index i_first = 0;
  const Index i_last  = 29;
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature (pressure in [hPa] therefore the factor 0.01)
  for ( Index i=0; i<n_p; ++i )
    {
      // here the total pressure is not multiplied by the H2O vmr for the 
      // P_H2O calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pwv_dummy = Pa_to_kPa * p_abs[i];
      // relative inverse temperature [1]
      Numeric theta     = (300.0 / t_abs[i]);
      // H2O partial pressure [kPa]
      Numeric pwv       = Pa_to_kPa * p_abs[i] * vmr[i];
      // dry air partial pressure [kPa]
      Numeric pda       = (Pa_to_kPa * p_abs[i]) - pwv;
      // H2O continuum absorption [dB/km/GHz^2] like in the original MPM89
      Numeric Nppc      = CC * pwv_dummy * pow(theta, 3.0) * 1.000e-5 *
                              ( (0.113 * pda) + (3.57 * pwv * pow(theta, 7.5)) );
      
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
        {
          // input frequency in [GHz]
          Numeric ff    = f_mono[s] * Hz_to_GHz; 
          // H2O line contribution at position f 
          Numeric Nppl  = 0.000;
          
          // Loop over MPM89 spectral lines:
          for ( Index l = i_first; l <= i_last; ++l )
            {
              // line strength [kHz]
              Numeric strength = CL * pwv_dummy * mpm89[l][1] * 
                                      pow(theta, 3.5) * exp(mpm89[l][2]*(1.000-theta));
              // line broadening parameter [GHz]
              Numeric gam      = CW * mpm89[l][3] * 0.001 * 
                                      ( mpm89[l][5] * pwv * pow(theta, mpm89[l][6]) +  
                                      pda * pow(theta, mpm89[l][4]) );
              // Doppler line width [GHz]
              // Numeric gamd     = 1.46e-6 * mpm89[l][0] / sqrt(theta);
              // effective line width [GHz]
              // gam              = 0.535 * gam + sqrt(0.217*gam*gam + gamd*gamd);  
              // H2O line absorption [dB/km/GHz] like in the original MPM89
              Nppl            += strength * MPMLineShapeFunction(gam, mpm89[l][0], ff); 
            }
          // xsec = abs/vmr [1/m] but MPM89 is in [dB/km] --> conversion necessary
          xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * ( Nppl + (Nppc * ff) );
        }
    }
  return;
}
//
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            H2O (lines+continuum) according to MPM93 [1/m]
   \param    CCin           scaling factor for the H2O-continuum  [1]
   \param    CLin           scaling factor for the line strengths [1]
   \param    CWin           scaling factor for the line widths    [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, and CWin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, and CWin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93', 'MPM93Lines', 'MPM93Continuum', 
             and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \attention The H2O lines at 547.676440 GHz and 552.020960 GHz are isotopic lines:<br>
              547 GHz is from the isotope 1-8-1 (HITRAN code 181, JPL code 20003) with an 
              isotopic ratio of 0.00199983 and <br>
              552 GHz is from the isotope 1-7-1  (HITRAN code 171, JPL code 19003) with an 
              isotopic ratio of 0.00037200.<br>
              The original source code of MPM93 has these isotopic ratios not included 
              in the line strength parameter b1, which is an error.<br> 
              In the arts implementation the line strength parameter b1 of these two lines 
              is multiplied with the appropriate isotopic ratio.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,
                       ConstVectorView   f_mono,
                       ConstVectorView   p_abs,
                       ConstVectorView   t_abs,
                       ConstVectorView   vmr )
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0           1        2       3        4      5      6
  //         f0          b1       b2      b3       b4     b5     b6
  //        [GHz]     [kHz/hPa]   [1]   [MHz/hPa]  [1]    [1]    [1]
  const Numeric mpm93[35][7] = { 
    {    22.235080,    0.01130,  2.143,   2.811,   4.80,  0.69,  1.00},
    {    67.803960,    0.00012,  8.735,   2.858,   4.93,  0.69,  0.82},
    {   119.995940,    0.00008,  8.356,   2.948,   4.78,  0.70,  0.79},
    {   183.310091,    0.24200,  0.668,   3.050,   5.30,  0.64,  0.85},
    {   321.225644,    0.00483,  6.181,   2.303,   4.69,  0.67,  0.54},
    {   325.152919,    0.14990,  1.540,   2.783,   4.85,  0.68,  0.74},
    {   336.222601,    0.00011,  9.829,   2.693,   4.74,  0.69,  0.61},
    {   380.197372,    1.15200,  1.048,   2.873,   5.38,  0.54,  0.89},
    {   390.134508,    0.00046,  7.350,   2.152,   4.81,  0.63,  0.55},
    {   437.346667,    0.00650,  5.050,   1.845,   4.23,  0.60,  0.48},
    {   439.150812,    0.09218,  3.596,   2.100,   4.29,  0.63,  0.52},
    {   443.018295,    0.01976,  5.050,   1.860,   4.23,  0.60,  0.50},
    {   448.001075,    1.03200,  1.405,   2.632,   4.84,  0.66,  0.67},
    {   470.888947,    0.03297,  3.599,   2.152,   4.57,  0.66,  0.65},
    {   474.689127,    0.12620,  2.381,   2.355,   4.65,  0.65,  0.64},
    {   488.491133,    0.02520,  2.853,   2.602,   5.04,  0.69,  0.72},
    {   503.568532,    0.00390,  6.733,   1.612,   3.98,  0.61,  0.43},
    {   504.482692,    0.00130,  6.733,   1.612,   4.01,  0.61,  0.45},
//    {   547.676440,    0.97010,  0.114,   2.600,   4.50,  0.70,  1.00},
//    {   552.020960,    1.47700,  0.114,   2.600,   4.50,  0.70,  1.00},
    {   547.676440,    0.97010*0.00199983,  0.114,   2.600,   4.50,  0.70,  1.00}, // isotopic ratio multiplied
    {   552.020960,    1.47700*0.00037200,  0.114,   2.600,   4.50,  0.70,  1.00}, // isotopic ratio multiplied
    {   556.936002,   48.74000,  0.159,   3.210,   4.11,  0.69,  1.00},
    {   620.700807,    0.50120,  2.200,   2.438,   4.68,  0.71,  0.68},
    {   645.866155,    0.00713,  8.580,   1.800,   4.00,  0.60,  0.50}, // ?? JPL tag 18003 (H2O) f_o = 645.7660100GHz
    {   658.005280,    0.03022,  7.820,   3.210,   4.14,  0.69,  1.00},
    {   752.033227,   23.96000,  0.396,   3.060,   4.09,  0.68,  0.84},
    {   841.053973,    0.00140,  8.180,   1.590,   5.76,  0.33,  0.45},
    {   859.962313,    0.01472,  7.989,   3.060,   4.09,  0.68,  0.84},
    {   899.306675,    0.00605,  7.917,   2.985,   4.53,  0.68,  0.90},
    {   902.616173,    0.00426,  8.432,   2.865,   5.10,  0.70,  0.95},
    {   906.207325,    0.01876,  5.111,   2.408,   4.70,  0.70,  0.53},
    {   916.171582,    0.83400,  1.442,   2.670,   4.78,  0.70,  0.78},
    {   923.118427,    0.00869, 10.220,   2.900,   5.00,  0.70,  0.80},
    {   970.315022,    0.89720,  1.920,   2.550,   4.94,  0.64,  0.67},
    {   987.926764,   13.21000,  0.258,   2.985,   4.55,  0.68,  0.90},
    {  1780.000000, 2230.00000,  0.952,  17.620,  30.50,  2.00,  5.00}};

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21)
  const Numeric CC_MPM93 = 1.00000;
  const Numeric CL_MPM93 = 1.00000;
  const Numeric CW_MPM93 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  // number of lines of Liebe line catalog (0-33 lines, 34 cont. pseudo line)
  Index i_first = 0;
  Index i_last  = 34;
  if ( model == "MPM93" )
    {
      CC      = CC_MPM93;
      CL      = CL_MPM93;
      CW      = CW_MPM93;
      i_first = 0;
      i_last  = 34;
    }
  else if ( model == "MPM93Lines" )
    {
      CC      = 0.000;
      CL      = CL_MPM93;
      CW      = CW_MPM93;
      i_first = 0;
      i_last  = 33;
    }
  else if ( model == "MPM93Continuum" )
    {
      CC      = CC_MPM93;
      CL      = 0.000;
      CW      = 0.000;
      i_first = 34;
      i_last  = 34;
    }
  else if ( model == "user" )
    {
      CC      = CCin;
      CL      = CLin;
      CW      = CWin;
      i_first = 0;
      i_last  = 34;

    }
  else
    {
      ostringstream os;
      os << "H2O-MPM93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM93', 'MPM93Lines', 'MPM93Continuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "H2O-MPM93: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n";
  
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in hPa therefore the factor 0.01)
  for ( Index i=0; i<n_p; ++i )
    {
      // here the total pressure is not multiplied by the H2O vmr for the 
      // P_H2O calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pwv_dummy = Pa_to_hPa * p_abs[i];
      // relative inverse temperature [1]
      Numeric theta    = (300.0 / t_abs[i]);
      // H2O partial pressure [hPa]
      Numeric pwv      = Pa_to_hPa * p_abs[i] * vmr[i];
      // dry air partial pressure [hPa]
      Numeric pda      = (Pa_to_hPa * p_abs[i]) - pwv;
      // Loop over MPM93 spectral lines:
      
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
        {
          // input frequency in [GHz]
          Numeric ff = f_mono[s] * Hz_to_GHz; 
          
          for ( Index l = i_first; l <= i_last; ++l )
            {
              // line strength [ppm]. The missing vmr of H2O will be multiplied 
              // at the stage of absorption calculation: abs / vmr * xsec.
              Numeric strength = 0.00;
              Numeric gam = 0.00;
              if ( (l >= 0) && (l <= 33) ) // ---- just the lines ------------------
                {
                  strength = CL * pwv_dummy * mpm93[l][1] * 
                                  pow(theta, 3.5)  * exp(mpm93[l][2]*(1.0-theta));
                  // line broadening parameter [GHz]
                  gam      = CW * mpm93[l][3] * 0.001 * 
                                  ( (mpm93[l][4] * pwv * pow(theta, mpm93[l][6])) +  
                                  (                pda * pow(theta, mpm93[l][5])) );
                }
              else if ( l == 34 ) // ----- just the continuum pseudo-line ----------
                {
                  strength = CC * pwv_dummy * mpm93[l][1] * 
                                  pow(theta, 3.5)  * exp(mpm93[l][2]*(1.0-theta));
                  // line broadening parameter [GHz]
                  gam      = mpm93[l][3] * 0.001 * 
                             ( (mpm93[l][4] * pwv * pow(theta, mpm93[l][6])) +  
                             (                pda * pow(theta, mpm93[l][5])) );
                }
              else // ----- if something strange happens ---------------------------
                {
                  ostringstream os;
                  os << "H2O-MPM93: wrong line number detected l=" << l << " (0-34)\n";
                  throw runtime_error(os.str());
                  return;
                } // ---------------------------------------------------------------
              // Doppler line width [GHz]
              // Numeric gamd     = 1.46e-6 * mpm93[l][0] / sqrt(theta);
              // effective line width [GHz]
              //gam              = 0.535 * gam + sqrt(0.217*gam*gam + gamd*gamd); 
              // absorption [dB/km] like in the original MPM93
              Numeric Npp = strength * MPMLineShapeFunction(gam, mpm93[l][0], ff); 
              // xsec = abs/vmr [1/m] but MPM89 is in [dB/km] --> conversion necessary
              xsec(s,i)   += dB_km_to_1_m * 0.1820 * ff * Npp;
            }
        }
    }
  return;
}
//
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            H2O (lines+continuum) according to P. W. Rosenkranz, 1998 [1/m]
   \param    CCin           scaling factor for the H2O-continuum  [1]
   \param    CLin           scaling factor for the line strengths [1]
   \param    CWin           scaling factor for the line widths    [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, and CWin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, and CWin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'RosenkranzLines', 'RosenkranzContinuum', 
             and 'user'. See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and 
             Radio Science, Vol. 34(4), 1025, 1999.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void PWR98H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
                       const Numeric     CLin,       // line strength scale factor
                       const Numeric     CWin,       // line broadening scale factor
                       const String&     model,
                       ConstVectorView   f_mono,
                       ConstVectorView   p_abs,
                       ConstVectorView   t_abs,
                       ConstVectorView   vmr )
{
  //   REFERENCES:
  //   LINE INTENSITIES FROM HITRAN92 (SELECTION THRESHOLD=
  //     HALF OF CONTINUUM ABSORPTION AT 1000 MB).
  //   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED:
  //     H.J.LIEBE AND T.A.DILLON, J.CHEM.PHYS. V.50, PP.727-732 (1969) &
  //     H.J.LIEBE ET AL., JQSRT V.9, PP. 31-47 (1969)  (22GHz);
  //     A.BAUER ET AL., JQSRT V.37, PP.531-539 (1987) & 
  //     ASA WORKSHOP (SEPT. 1989) (380GHz);
  //     AND A.BAUER ET AL., JQSRT V.41, PP.49-54 (1989) (OTHER LINES).
  //   AIR-BROADENED CONTINUUM BASED ON LIEBE & LAYTON, NTIA 
  //     REPORT 87-224 (1987); SELF-BROADENED CONTINUUM BASED ON 
  //     LIEBE ET AL, AGARD CONF. PROC. 542 (MAY 1993), 
  //     BUT READJUSTED FOR LINE SHAPE OF
  //     CLOUGH et al, ATMOS. RESEARCH V.23, PP.229-241 (1989).
  //
  // Coefficients are from P. W. Rosenkranz., Radio Science, 33(4), 919, 1998
  // line frequencies [GHz]
  const Numeric PWRfl[15] = {  22.2350800, 183.3101170, 321.2256400, 325.1529190, 380.1973720, 
                              439.1508120, 443.0182950, 448.0010750, 470.8889470,  474.6891270, 
                              488.4911330, 556.9360020,  620.7008070, 752.0332270, 916.1715820 };
  // line intensities at 300K [Hz * cm2] (see Janssen Appendix to Chap.2 for this)
  const Numeric PWRs1[15] = { 1.31e-14,  2.273e-12, 8.036e-14, 2.694e-12, 2.438e-11, 
                              2.179e-12, 4.624e-13, 2.562e-11, 8.369e-13, 3.263e-12, 
                              6.659e-13, 1.531e-9,  1.707e-11, 1.011e-9,  4.227e-11 };
  // T coeff. of intensities [1]
  const Numeric PWRb2[15] = { 2.144, 0.668, 6.179, 1.541, 1.048, 
                              3.595, 5.048, 1.405, 3.597, 2.379,
                              2.852, 0.159, 2.391, 0.396, 1.441 };
  // air-broadened width parameters at 300K [GHz/hPa]
  const Numeric PWRw3[15] = { 0.00281, 0.00281, 0.00230, 0.00278, 0.00287,
                              0.00210, 0.00186, 0.00263, 0.00215, 0.00236,
                              0.00260, 0.00321, 0.00244, 0.00306, 0.00267 };
  // T-exponent of air-broadening [1]
  const Numeric PWRx[15]  = { 0.69, 0.64, 0.67, 0.68, 0.54,
                              0.63, 0.60, 0.66, 0.66, 0.65,
                              0.69, 0.69, 0.71, 0.68, 0.70 };
  // self-broadened width parameters at 300K [GHz/hPa]
  const Numeric PWRws[15] = { 0.01349, 0.01491, 0.01080, 0.01350, 0.01541,
                              0.00900, 0.00788, 0.01275, 0.00983, 0.01095,
                              0.01313, 0.01320, 0.01140, 0.01253, 0.01275 };
  // T-exponent of self-broadening [1]
  const Numeric PWRxs[15] = { 0.61, 0.85, 0.54, 0.74, 0.89,
                              0.52, 0.50, 0.67, 0.65, 0.64,
                              0.72, 1.00, 0.68, 0.84, 0.78 };


  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM87 model (P. W. Rosenkranz., Radio Science, 33(4), 919, 1998):
  const Numeric CC_PWR98 = 1.00000;
  const Numeric CL_PWR98 = 1.00000;
  const Numeric CW_PWR98 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  if ( model == "Rosenkranz" )
    {
      CC = CC_PWR98;
      CL = CL_PWR98;
      CW = CW_PWR98;
    }
  else if ( model == "RosenkranzLines" )
    {
      CC = 0.000;
      CL = CL_PWR98;
      CW = CW_PWR98;
    }
  else if ( model == "RosenkranzContinuum" )
    {
      CC = CC_PWR98;
      CL = 0.000;
      CW = 0.000;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CL = CLin;
      CW = CWin;
    }
  else
    {
      ostringstream os;
      os << "H2O-PWR98: ERROR! Wrong model values given.\n"
         << "Valid models are: 'Rosenkranz', 'RosenkranzLines', 'RosenkranzContinuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "H2O-PWR98: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n";
  
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      // here the total pressure is not multiplied by the H2O vmr for the 
      // P_H2O calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pvap_dummy = Pa_to_hPa * p_abs[i];
      // water vapor partial pressure [hPa]
      Numeric pvap       = Pa_to_hPa * p_abs[i] * vmr[i];
      // dry air partial pressure [hPa]
      Numeric pda        = (Pa_to_hPa * p_abs[i]) - pvap;
      // Rosenkranz number density  (Rosenkranz H2O mass density in [g/m³])
      // [g/m³]    =  [g*K / Pa*m³]  *  [Pa/K]
      // rho       =   (M_H2O / R)   *  (P_H2O / T)
      // rho       =      2.1667     *  p_abs * vmr / t_abs
      // den       = 3.335e16 * rho

      //      Numeric den        = 3.335e16 * (2.1667 * p_abs[i] * vmr[i] / t_abs[i]);
      // FIXME: SAB+TKS: Commented this out, den_dummy is used currently.

      Numeric den_dummy  = 3.335e16 * (2.1667 * p_abs[i] / t_abs[i]);
      // inverse relative temperature [1]
      Numeric ti         = (300.0 / t_abs[i]);
      Numeric ti2        = pow(ti, 2.5);
      
      // continuum term [Np/km/GHz2]
      Numeric con = CC * pvap_dummy * pow(ti, 3.0) * 1.000e-9 * 
                         ( (0.543 * pda) + (17.96 * pvap * pow(ti, 4.5)) );
          
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
        {
          // input frequency in [GHz]
          Numeric ff  = f_mono[s] * Hz_to_GHz;
          // line contribution at position f
          Numeric sum = 0.000;
          
          // Loop over spectral lines
          for (Index l = 0; l < 15; l++) 
            {
              Numeric width    = CW * ( PWRw3[l] * pda  * pow(ti, PWRx[l]) + 
                                        PWRws[l] * pvap * pow(ti, PWRxs[l]) );
              Numeric wsq      = width * width;
              Numeric strength = CL * PWRs1[l] * ti2 * exp(PWRb2[l]*(1.0 - ti));
              // frequency differences
              Numeric df0      = ff - PWRfl[l];
              Numeric df1      = ff + PWRfl[l];
              // use Clough's definition of local line contribution
              Numeric base     = width / (wsq + 562500.000);
              // positive and negative resonances
              Numeric res      = 0.000;
              if (fabs(df0) < 750.0) res += width / (df0*df0 + wsq) - base;
              if (fabs(df1) < 750.0) res += width / (df1*df1 + wsq) - base;
              sum             += strength * res * pow( (ff/PWRfl[l]), 2.0 );
            }
          // line term [Np/km]
          Numeric absl = 0.3183e-4 * den_dummy * sum;
          // xsec = abs/vmr [1/m] (Rosenkranz model in [Np/km])
          // 4.1907e-5 = 0.230259 * 0.1820 * 1.0e-3    (1/(10*log(e)) = 0.230259)
          xsec(s,i)  += 1.000e-3 * ( absl + (con * ff * ff) );    
        }
    }
  return;
}
//
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            H2O (lines+continuum) according to Cruz-Pol 1998 [1/m]
   \param    CCin           scaling factor for the H2O-continuum  [1]
   \param    CLin           scaling factor for the line strengths [1]
   \param    CWin           scaling factor for the line widths    [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, and CWin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, and CWin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'CruzPol', 'CruzPolLines', 'CruzPolContinuum', 
             and 'user'. See the user guide for detailed explanations.

   \remark   Reference: S. L. Cruz-Pol et al., Radio Science, 33(5), 1319, 1998.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void CP98H2OAbsModel( MatrixView        xsec,
                      const Numeric     CCin,       // continuum scale factor 
                      const Numeric     CLin,       // line strength scale factor
                      const Numeric     CWin,       // line broadening scale factor
                      const String&     model,
                      ConstVectorView   f_mono,
                      ConstVectorView   p_abs,
                      ConstVectorView   t_abs,
                      ConstVectorView   vmr )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the CP98 model (S. L. Cruz-Pol et al., Radio Science, 33(5), 1319, 1998):
  const Numeric CC_CP98 = 1.2369; // +/- 0.155  !LARGE!
  const Numeric CL_CP98 = 1.0639; // +/- 0.016
  const Numeric CW_CP98 = 1.0658; // +/- 0.0096
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  if ( model == "CruzPol" )
    {
      CC = CC_CP98;
      CL = CL_CP98;
      CW = CW_CP98;
    }
  else if ( model == "CruzPolLine" )
    {
      CC = 0.000;
      CL = CL_CP98;
      CW = CW_CP98;
    }
  else if ( model == "CruzPolContinuum" )
    {
      CC = CC_CP98;
      CL = 0.000;
      CW = 0.000;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CL = CLin;
      CW = CWin;
    }
  else
    {
      ostringstream os;
      os << "H2O-CP98: ERROR! Wrong model values given.\n"
         << "Valid models are: 'CruzPol', 'CruzPolLine', 'CruzPolContinuum', and 'user'" << "\n";
      throw runtime_error(os.str());
    }
  out2  << "H2O-CP98: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n";
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in [hPa] therefore the factor 0.01)
  for ( Index i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(H2O) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
        {
          // relative inverse temperature [1]
          Numeric theta = (300.0 / t_abs[i]);
          // H2O partial pressure [hPa]
          Numeric pwv   = Pa_to_hPa * p_abs[i] * vmr[i];
          // dry air partial pressure [hPa]
          Numeric pda   = (Pa_to_hPa * p_abs[i]) - pwv;
          // line strength
          Numeric TL    = CL * 0.0109 * pwv * pow(theta,3.5) * exp(2.143*(1.0-theta));
          // line broadening parameter [GHz]
          Numeric gam   = CW * 0.002784 *  
                          ( (pda * pow(theta,0.6)) + (4.80 * pwv * pow(theta,1.1)) );
          // continuum term
          Numeric TC    = CC * pwv * pow(theta, 3.0) * 1.000e-7 * 
                          ( (0.113 * pda) + (3.57 * pwv * pow(theta,7.5)) );

          // Loop over input frequency
          for ( Index s=0; s<n_f; ++s )
            {
              // input frequency in [GHz]
              Numeric ff  = f_mono[s] * Hz_to_GHz; 
              Numeric TSf = MPMLineShapeFunction(gam, 22.235080, ff); 
              // xsec = abs/vmr [1/m] (Cruz-Pol model in [Np/km])
              xsec(s,i) += 4.1907e-5 * ff * ( (TL * TSf) + (ff * TC) ) / vmr[i];
            }
        }
    }
  return;
}
//
// #################################################################################
//
/**
   \param    xsec Output:          cross section (absorption/volume mixing ratio) of the 
                            H2O-H2O continuum [1/m]
   \param    C              constant absorption strength    [1/m / (Hz*Pa)²]
   \param    x              temperature exponent of (300/T) [1]
   \param    model          allows user defined input parameter set 
                            (C and x)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    t_abs          predefined temperature grid     [K] 
   \param    p_abs          predefined pressure grid       [Pa]
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters C and x 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'CruzPol', 'MPM89', 'MPM87', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and 
             Radio Science, Vol. 34(4), 1025, 1999.             
   \author Thomas Kuhn
   \date 2001-11-05
 */ 
void Standard_H2O_self_continuum( MatrixView        xsec,
                                  const Numeric     Cin,
                                  const Numeric     xin,
                                  const String&     model,
                                  ConstVectorView   f_mono,
                                  ConstVectorView   p_abs,
                                  ConstVectorView   t_abs,
                                  ConstVectorView   vmr  )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model (Radio Science, 33(4), 919, 1998):
  const Numeric Cs_PWR   = 1.796e-33;  //  [1/m / (Hz²*Pa²)]
  const Numeric xs_PWR   = 4.5;        //  [1]
  // standard values for the Cruz-Pol model (Radio Science, 33(5), 1319, 1998):
  const Numeric Cs_CP    = 1.851e-33;  //  [1/m / (Hz²*Pa²)]
  const Numeric xs_CP    = 7.5;        //  [1]
  // standard values for the MPM89 model (Int. J. Inf. and Millim. Waves, 10(6), 1989, 631):
  const Numeric Cs_MPM89 = 1.500e-33;  //  [1/m / (Hz²*Pa²)]
  const Numeric xs_MPM89 = 7.5;        //  [1]
  // standard values for the MPM87 model (Radio Science, 20(5), 1985, 1069):
  const Numeric Cs_MPM87 = 1.500e-33;  //  [1/m / (Hz²*Pa²)]
  const Numeric xs_MPM87 = 7.5;        //  [1]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model goes for values!!):
  Numeric C, x;
   if ( model == "Rosenkranz" )
     {
       C = Cs_PWR;
       x = xs_PWR;
     }
   else if ( model == "CruzPol" )
     {
       C = Cs_CP;
       x = xs_CP;
     }
   else if ( model == "MPM89" )
     {
       C = Cs_MPM89;
       x = xs_MPM89;
     }
   else if ( model == "MPM87" )
     {
       C = Cs_MPM87;
       x = xs_MPM87;
     }
   else if ( model == "user" )
     {
       C = Cin;
       x = xin;
     }
   else
     {
       ostringstream os;
       os << "H2O-SelfContStandardType: ERROR! Wrong model values given.\n"
          << "allowed models are: 'Rosenkranz', 'CruzPol', 'MPM89', 'MPM87', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out2  << "H2O-SelfContStandardType: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";



  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( Index i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The second vmr of H2O will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy =
        C * pow( 300./t_abs[i], x+3. ) * pow( p_abs[i], 2 ) * vmr[i];

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
        {
          xsec(s,i) += dummy * pow( f_mono[s], 2 );
          //      cout << "xsec(" << s << "," << i << "): " << xsec(s,i) << "\n";
        }
    }
}
//
// #################################################################################
//
/**
   \param   xsec Output:           cross section (absorption/volume mixing ratio) of the 
                            H2O-dry air continuum [1/m]
   \param    C              constant absorption strength [1/m / (Hz*Pa)²]
   \param    x              temperature exponent         [1] 
   \param    model          allows user defined input parameter set 
                            (C and x)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid    [Hz]
   \param    t_abs          predefined temperature grid  [K] 
   \param    p_abs          predefined pressure          [Pa]
   \param    vmr            H2O volume mixing ratio     [1] 

   \note     Except for  model 'user' the input parameters C and x 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'CruzPol', 'MPM89', 'MPM87', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz., Radio Science, 33(4), 919, 1998 and 
             Radio Science, Vol. 34(4), 1025, 1999.

   \author Thomas Kuhn
   \date 2001-08-03
 */ 
void Standard_H2O_foreign_continuum( MatrixView        xsec,
                                     const Numeric       Cin,
                                     const Numeric       xin,
                                     const String&     model,
                                     ConstVectorView   f_mono,
                                     ConstVectorView   p_abs,
                                     ConstVectorView   t_abs,
                                     ConstVectorView   vmr       )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model (Radio Science, 33(4), 919, 1998):
  const Numeric Cf_PWR   = 5.43e-35 ;  //  [1/m / (Hz²*Pa²)]
  const Numeric xf_PWR   = 0.0;        //  [1]
  // standard values for the Cruz-Pol model (Radio Science, 33(5), 1319, 1998):
  const Numeric Cf_CP    = 5.85e-35;  //  [1/m / (Hz²*Pa²)]
  const Numeric xf_CP    = 0.0;        //  [1]
  // standard values for the MPM89 model (Int. J. Inf. and Millim. Waves, 10(6), 1989, 631):
  const Numeric Cf_MPM89 = 4.74e-35;  //  [1/m / (Hz²*Pa²)]
  const Numeric xf_MPM89 = 0.0;        //  [1]
  // standard values for the MPM87 model (Radio Science, 20(5), 1985, 1069):
  const Numeric Cf_MPM87 = 4.74e-35;  //  [1/m / (Hz²*Pa²)]
  const Numeric xf_MPM87 = 0.0;        //  [1]
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model goes for values!!):
  Numeric C, x;
   if ( model == "Rosenkranz" )
     {
       C = Cf_PWR;
       x = xf_PWR;
     }
   else if ( model == "CruzPol" )
     {
       C = Cf_CP;
       x = xf_CP;
     }
   else if ( model == "MPM89" )
     {
       C = Cf_MPM89;
       x = xf_MPM89;
     }
   else if ( model == "MPM87" )
     {
       C = Cf_MPM87;
       x = xf_MPM87;
     }
   else if ( model == "user" )
     {
       C = Cin;
       x = xin;
     }
   else
     {
       ostringstream os;
       os << "H2O-ForeignContStandardType: ERROR! Wrong model values given.\n"
          << "allowed models are: 'Rosenkranz', 'CruzPol', 'MPM89', 'MPM87', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out2  << "H2O-ForeignContStandardType: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      // Dry air partial pressure: p_dry := p_tot - p_h2o.
      Numeric pdry  = p_abs[i] * (1.000e0-vmr[i]);
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The vmr of H2O will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy = C * pow( 300./t_abs[i], x+3. ) * p_abs[i] * pdry;

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
        {
          xsec(s,i) += dummy * pow( f_mono[s], 2 );
          //      cout << "xsec(" << s << "," << i << "): " << xsec(s,i) << "\n";
        }
    }
}
//
// #################################################################################
//
// MPM93 H2O pseudo continuum line parameters:
// see publication side of National Telecommunications and Information Administration
//   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
// and ftp side for downloading the MPM93 original source code:
//   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            H2O according to MPM87 [1/m]
   \param    fcenter        continuum pseudo-line center frequency [Hz]
   \param    b1             continuum pseudo-line line strength [Hz/Pa]
   \param    b2             continuum pseudo-line line strength temperature exponent [1]
   \param    b3             continuum pseudo-line line broadening parameter [Hz/Pa]
   \param    b4             continuum pseudo-line line broadening parameter [1]
   \param    b5             continuum pseudo-line line broadening parameter [1]
   \param    b6             continuum pseudo-line line broadening parameter [1]
   \param    model          allows user defined input parameter set 
                            (fcenter and b1 to b6)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters fcenter and b1 to b6
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93_H2O_continuum( MatrixView          xsec,
                          const Numeric       fcenter,
                          const Numeric       b1,
                          const Numeric       b2,
                          const Numeric       b3,
                          const Numeric       b4,
                          const Numeric       b5,
                          const Numeric       b6,
                          const String&       model,
                          ConstVectorView     f_mono,
                          ConstVectorView     p_abs,
                          ConstVectorView     t_abs,
                          ConstVectorView     vmr        )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 H2O continuum model 
  // (AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric MPM93fo_orig =  1780.000e9;  // [Hz]
  const Numeric MPM93b1_orig = 22300.000;    // [Hz/Pa]
  const Numeric MPM93b2_orig =     0.952;    // [1]
  const Numeric MPM93b3_orig =    17.600e4;  // [Hz/Pa]
  const Numeric MPM93b4_orig =    30.500;    // [1]
  const Numeric MPM93b5_orig =     2.000;    // [1]
  const Numeric MPM93b6_orig =     5.000;    // [1]
  // ---------------------------------------------------------------------------------------
  

  // select the parameter set (!!model goes for values!!):
  Numeric MPM93fopcl, MPM93b1pcl, MPM93b2pcl, 
          MPM93b3pcl, MPM93b4pcl, MPM93b5pcl,
          MPM93b6pcl;
  if ( model == "MPM93" )
    {
      MPM93fopcl =  MPM93fo_orig;
      MPM93b1pcl =  MPM93b1_orig;
      MPM93b2pcl =  MPM93b2_orig;
      MPM93b3pcl =  MPM93b3_orig;
      MPM93b4pcl =  MPM93b4_orig;
      MPM93b5pcl =  MPM93b5_orig;
      MPM93b6pcl =  MPM93b6_orig;
    }
  else if ( model == "user" )
     {
      MPM93fopcl =  fcenter;
      MPM93b1pcl =  b1;
      MPM93b2pcl =  b2;
      MPM93b3pcl =  b3;
      MPM93b4pcl =  b4;
      MPM93b5pcl =  b5;
      MPM93b6pcl =  b6;
    }
   else
     {
       ostringstream os;
       os << "H2O-ContMPM93: ERROR! Wrong model values given.\n"
          << "allowed models are: 'MPM93', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out2  << "H2O-ContMPM93: (model=" << model << ") parameter values in use:\n" 
         << " fo = " << MPM93fopcl << "\n"
         << " b1 = " << MPM93b1pcl << "\n"
         << " b2 = " << MPM93b2pcl << "\n"
         << " b3 = " << MPM93b3pcl << "\n"
         << " b4 = " << MPM93b4pcl << "\n"
         << " b5 = " << MPM93b5pcl << "\n"
         << " b6 = " << MPM93b6pcl << "\n";

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  

  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      Numeric th = 300.0 / t_abs[i];
      // the vmr of H2O will be multiplied at the stage of absorption calculation:
      // abs / vmr * xsec.
      Numeric strength =  MPM93b1pcl * p_abs[i] * pow( th, 3.5 ) * exp(MPM93b2pcl * (1 - th));
      Numeric gam      =  MPM93b3pcl * 0.001 * 
                          ( MPM93b4pcl * p_abs[i] * vmr[i]       * pow( th, MPM93b6pcl ) +  
                                         p_abs[i]*(1.000-vmr[i]) * pow( th, MPM93b5pcl ) );
      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
        {
          // xsec = abs/vmr [1/m] but MPM89 is in [dB/km] --> conversion necessary
          xsec(s,i) += dB_km_to_1_m * 0.1820 *
                       f_mono[s] * strength * 
                       MPMLineShapeFunction(gam, MPM93fopcl, f_mono[s]); 
        }
    }
  return;
}
//
// #################################################################################
// ################################# OXYGEN MODELS #################################
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            O2 according to MPM93 [1/m]
   \param    CCin           scaling factor for the O2-continuum   [1]
   \param    CLin           scaling factor for the O2-line strengths [1]
   \param    CWin           scaling factor for the O2-line widths    [1]
   \param    COin           scaling factor for the O2-line coupling  [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, CWin, and COin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid           [Hz]
   \param    p_abs          predefined pressure                 [Pa]
   \param    t_abs          predefined temperature grid         [K] 
   \param    h2o_abs        H2O volume mixing ratio profile    [1]
   \param    vmr            O2 volume mixing ratio profile     [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, CWin, and COin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93', 'MPM93Lines', 'MPM93Continuum', 'MPM93NoCoupling', 
             and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93O2AbsModel( MatrixView          xsec,
                      const Numeric       CCin,       // continuum scale factor 
                      const Numeric       CLin,       // line strength scale factor
                      const Numeric       CWin,       // line broadening scale factor
                      const Numeric       COin,       // line coupling scale factor
                      const String&       model,
                      ConstVectorView     f_mono,
                      ConstVectorView     p_abs,
                      ConstVectorView     t_abs,
                      ConstVectorView     h2o_abs,
                      ConstVectorView     vmr )
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0           1        2       3        4      5      6
  //         f0          a1       a2      a3       a4     a5     a6
  //        [GHz]     [kHz/hPa]   [1]   [MHz/hPa]  [1]    [10³/hPa]
  const Numeric mpm93[44][7] = { 
    {   50.474238,       0.094,      9.694,    0.890,     0.0,   0.240,    0.790},
    {   50.987749,       0.246,      8.694,    0.910,     0.0,   0.220,    0.780},
    {   51.503350,       0.608,      7.744,    0.940,     0.0,   0.197,    0.774},
    {   52.021410,       1.414,      6.844,    0.970,     0.0,   0.166,    0.764},
    {   52.542394,       3.102,      6.004,    0.990,     0.0,   0.136,    0.751},
    {   53.066907,       6.410,      5.224,    1.020,     0.0,   0.131,    0.714},
    {   53.595749,      12.470,      4.484,    1.050,     0.0,   0.230,    0.584},
    {   54.130000,      22.800,      3.814,    1.070,     0.0,   0.335,    0.431},
    {   54.671159,      39.180,      3.194,    1.100,     0.0,   0.374,    0.305},
    {   55.221367,      63.160,      2.624,    1.130,     0.0,   0.258,    0.339},
    {   55.783802,      95.350,      2.119,    1.170,     0.0,  -0.166,    0.705},
    {   56.264775,      54.890,      0.015,    1.730,     0.0,   0.390,   -0.113},
    {   56.363389,     134.400,      1.660,    1.200,     0.0,  -0.297,    0.753},
    {   56.968206,     176.300,      1.260,    1.240,     0.0,  -0.416,    0.742},
    {   57.612484,     214.100,      0.915,    1.280,     0.0,  -0.613,    0.697},
    {   58.323877,     238.600,      0.626,    1.330,     0.0,  -0.205,    0.051},
    {   58.446590,     145.700,      0.084,    1.520,     0.0,   0.748,   -0.146},
    {   59.164207,     240.400,      0.391,    1.390,     0.0,  -0.722,    0.266},
    {   59.590983,     211.200,      0.212,    1.430,     0.0,   0.765,   -0.090},
    {   60.306061,     212.400,      0.212,    1.450,     0.0,  -0.705,    0.081},
    {   60.434776,     246.100,      0.391,    1.360,     0.0,   0.697,   -0.324},
    {   61.150560,     250.400,      0.626,    1.310,     0.0,   0.104,   -0.067},
    {   61.800154,     229.800,      0.915,    1.270,     0.0,   0.570,   -0.761},
    {   62.411215,     193.300,      1.260,    1.230,     0.0,   0.360,   -0.777},
    {   62.486260,     151.700,      0.083,    1.540,     0.0,  -0.498,    0.097},
    {   62.997977,     150.300,      1.665,    1.200,     0.0,   0.239,   -0.768},
    {   63.568518,     108.700,      2.115,    1.170,     0.0,   0.108,   -0.706},
    {   64.127767,      73.350,      2.620,    1.130,     0.0,  -0.311,   -0.332},
    {   64.678903,      46.350,      3.195,    1.100,     0.0,  -0.421,   -0.298},
    {   65.224071,      27.480,      3.815,    1.070,     0.0,  -0.375,   -0.423},
    {   65.764772,      15.300,      4.485,    1.050,     0.0,  -0.267,   -0.575},
    {   66.302091,       8.009,      5.225,    1.020,     0.0,  -0.168,   -0.700},
    {   66.836830,       3.946,      6.005,    0.990,     0.0,  -0.169,   -0.735},
    {   67.369598,       1.832,      6.845,    0.970,     0.0,  -0.200,   -0.744},
    {   67.900867,       0.801,      7.745,    0.940,     0.0,  -0.228,   -0.753},
    {   68.431005,       0.330,      8.695,    0.920,     0.0,  -0.240,   -0.760},
    {   68.960311,       0.128,      9.695,    0.900,     0.0,  -0.250,   -0.765},
    {  118.750343,      94.500,      0.009,    1.630,     0.0,  -0.036,    0.009},
    {  368.498350,       6.790,      0.049,    1.920,     0.6,   0.000,    0.000},
    {  424.763124,      63.800,      0.044,    1.930,     0.6,   0.000,    0.000},
    {  487.249370,      23.500,      0.049,    1.920,     0.6,   0.000,    0.000},
    {  715.393150,       9.960,      0.145,    1.810,     0.6,   0.000,    0.000},
    {  773.839675,      67.100,      0.130,    1.820,     0.6,   0.000,    0.000},
    {  834.145330,      18.000,      0.147,    1.810 ,    0.6,   0.000,    0.000}};

  // number of lines of Liebe O2-line catalog (0-43 lines)
  const Index i_first = 0;
  const Index i_last  = 43;
  

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM87 model (Radio Science, 20(5), 1985, 1069):
  const Numeric CC_MPM93 = 1.00000;
  const Numeric CL_MPM93 = 1.00000;
  const Numeric CW_MPM93 = 1.00000;
  const Numeric CO_MPM93 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  if ( model == "MPM93" )
    {
      CC = CC_MPM93;
      CL = CL_MPM93;
      CW = CW_MPM93;
      CO = CO_MPM93;
    }
  else if ( model == "MPM93Lines" )
    {
      CC = 0.000;
      CL = CL_MPM93;
      CW = CW_MPM93;
      CO = CO_MPM93;
    }
  else if ( model == "MPM93Continuum" )
    {
      CC = CC_MPM93;
      CL = 0.000;
      CW = 0.000;
      CO = 0.000;
    }
  else if ( model == "MPM93NoCoupling" )
    {
      CC = CC_MPM93;
      CL = CL_MPM93;
      CW = CW_MPM93;
      CO = 0.000;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CL = CLin;
      CW = CWin;
      CO = COin;
    }
  else
    {
      ostringstream os;
      os << "O2-MPM93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM93' 'MPM93Lines' 'MPM93Continuum' 'MPM93NoCoupling' " 
         << "and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "O2-MPM93: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n"
        << " CO = " << CO << "\n";
  

  // O2 continuum parameters of MPM93:
  const Numeric S0 =  6.140e-5; // line strength                        [ppm]
  const Numeric G0 =  0.560e-3; // line width                           [GHz/hPa]
  const Numeric X0 =  0.800;    // temperature dependence of line width [1]

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in hPa therefore the factor 0.01)
  for ( Index i=0; i<n_p; ++i )
    {
      // check if O2-VMR will cause an underflow due to division by zero:
      if (vmr[i] < VMRCalcLimit)
        {
          ostringstream os;
          os << "ERROR: MPM93 O2 full absorption model has detected a O2 volume mixing ratio of " 
             << vmr[i] << " which is below the threshold of " << VMRCalcLimit << ".\n"
             << "Therefore no calculation is performed.\n";
          throw runtime_error(os.str());
          return;
        }

      // relative inverse temperature [1]
      Numeric theta     = (300.0 / t_abs[i]);
      // H2O partial pressure [hPa]
      Numeric pwv       = Pa_to_hPa * p_abs[i] * h2o_abs[i];
      // dry air partial pressure [hPa]
      Numeric pda       = (Pa_to_hPa * p_abs[i]) - pwv;
      // here the total pressure is devided by the O2 vmr for the 
      // P_dry calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pda_dummy = pda / vmr[i];
      // O2 continuum strength [ppm]
      Numeric strength_cont =  S0 * pda_dummy * pow( theta, 2 );
      // O2 continuum pseudo line broadening [GHz]
      Numeric gam_cont      =  G0 * (pwv+pda) *  pow( theta, X0 ); // GHz
      
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
        {
          // input frequency in [GHz]
          Numeric ff = f_mono[s] * Hz_to_GHz; 
          // O2 continuum absorption [1/m]
          // cross section: xsec = absorption / var
          // the vmr of O2 will be multiplied at the stage of absorption calculation:
          Numeric Nppc =  CC * strength_cont * ff * gam_cont /
            ( pow( ff, 2) + pow( gam_cont, 2) );
          
          // Loop over MPM93 O2 spectral lines:
          Numeric Nppl  = 0.0;
          for ( Index l = i_first; l <= i_last; ++l )
            {
              // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
              Numeric strength = CL * 1.000e-6  * pda_dummy * mpm93[l][1] / mpm93[l][0] * 
                                      pow(theta, 3) * exp(mpm93[l][2]*(1.0-theta));
              // line broadening parameter [GHz]
              Numeric gam      = CW * ( mpm93[l][3] * 0.001 * 
                                      ( (       pda * pow(theta, (0.8-mpm93[l][4]))) + 
                                        (1.10 * pwv * theta) ) );
              // line mixing parameter [1]
              //                  if (l < 11) CD = 1.1000;
              Numeric delta    = CO * ( (mpm93[l][5] + mpm93[l][6] * theta) * 
                                      (pda+pwv) * pow(theta, 0.8) * 0.001 );
              // absorption [dB/km] like in the original MPM93
              Nppl            += strength * MPMLineShapeO2Function(gam, mpm93[l][0], ff, delta); 
            }
          // in MPM93 there is a cutoff for O2 line absorption if abs_l < 0 
          if (Nppl < 0.0)  Nppl = 0.0; // absorption cannot be less than 0 according to MPM93.
          // O2 line absorption [1/m]
          // cross section: xsec = absorption / var
          // the vmr of O2 will be multiplied at the stage of absorption calculation:
          xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * (Nppl+Nppc);
        }
    }
  return;
}
//
// #################################################################################
// 
//   Oxygen complex at 60 GHz plus mm O2 lines plus O2 continuum
//
//   REFERENCE FOR EQUATIONS AND COEFFICIENTS:
//   P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
//   BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
//   AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
//   (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            O2 according to the P. W. Rosenkranz, 1993 [1/m]
   \param    CCin           O2-continuum scale factor  [1]
   \param    CLin           O2 line strength scale factor [1]
   \param    CWin           O2 line broadening scale factor [1]
   \param    COin           O2 line coupling scale factor [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, CWin, and COin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid        [Hz]
   \param    p_abs          predefined pressure              [Pa]
   \param    t_abs          predefined temperature grid      [K] 
   \param    vmrh2o         H2O volume mixing ratio profile [1]
   \param    vmr            O2 volume mixing ratio profile  [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, CWin, and COin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', 'RosenkranzLines', 
             'RosenkranzContinuum', 'RosenkranzNoCoupling', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference:  P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void PWR93O2AbsModel( MatrixView        xsec,
                      const Numeric     CCin,      // model parameter
                      const Numeric     CLin,      // model parameter
                      const Numeric     CWin,      // model parameter
                      const Numeric     COin,      // model parameter
                      const String&     model,
                      ConstVectorView   f_mono,
                      ConstVectorView   p_abs,
                      ConstVectorView   t_abs,
                      ConstVectorView   vmrh2o,
                      ConstVectorView   vmr )
{
  const Index n_lines = 40; // number of O2 lines in this model (range: 50-850 GHz)

  // lines are arranged 1-,1+,3-,3+,etc. in spin-rotation spectrum
  // line center frequency in [GHz]
  const Numeric F[n_lines] = { 118.7503,  56.2648,  62.4863,  58.4466, 
                           60.3061,  59.5910,  59.1642,  60.4348, 
                           58.3239,  61.1506,  57.6125,  61.8002,
                           56.9682,  62.4112,  56.3634,  62.9980, 
                           55.7838,  63.5685,  55.2214,  64.1278, 
                           54.6712,  64.6789,  54.1300,  65.2241,
                           53.5957,  65.7648,  53.0669,  66.3021, 
                           52.5424,  66.8368,  52.0214,  67.3696, 
                           51.5034,  67.9009, 368.4984, 424.7631,
                          487.2494, 715.3932, 773.8397, 834.1453};
  // line strength at T=300K in [cm² * Hz]
  const Numeric S300[n_lines] = { 0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14,
                             0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14,
                             0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,
                             0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14,
                             0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14,
                             0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,
                             0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15,
                             0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16,
                             0.8191E-17, 0.1423E-16, 0.6460E-15, 0.7047E-14, 
                             0.3011E-14, 0.1826E-14, 0.1152E-13, 0.3971E-14};
  // temperature exponent of the line strength in [1]
  const Numeric BE[n_lines] = {   0.009,   0.015,   0.083,   0.084, 
                             0.212,   0.212,   0.391,   0.391, 
                             0.626,   0.626,   0.915,   0.915, 
                             1.260,   1.260,   1.660,   1.665,   
                             2.119,   2.115,   2.624,   2.625, 
                             3.194,   3.194,   3.814,   3.814, 
                             4.484,   4.484,   5.224,   5.224, 
                             6.004,   6.004,   6.844,   6.844, 
                             7.744,   7.744,   0.048,   0.044, 
                             0.049,   0.145,   0.141,   0.145};

  // widths in MHz/mbar for the O2 continuum
  const Numeric WB300 = 0.56; // [MHz/mbar]
  const Numeric X     = 0.80; // [1]
  // line width parameter [GHz/bar]
  const Numeric W300[n_lines] = { 1.630, 1.646, 1.468, 1.449, 
                             1.382, 1.360, 1.319, 1.297, 
                             1.266, 1.248, 1.221, 1.207, 
                             1.181, 1.171, 1.144, 1.139, 
                             1.110, 1.108, 1.079, 1.078, 
                             1.050, 1.050, 1.020, 1.020, 
                             1.000, 1.000, 0.970, 0.970,
                             0.940, 0.940, 0.920, 0.920,
                             0.890, 0.890, 1.920, 1.920, 
                             1.920, 1.810, 1.810, 1.810};
  // y parameter for the calculation of Y [1/bar]
  const Numeric Y300[n_lines] = { -0.0233,  0.2408, -0.3486,  0.5227,
                             -0.5430,  0.5877, -0.3970,  0.3237, 
                             -0.1348,  0.0311,  0.0725, -0.1663, 
                              0.2832, -0.3629,  0.3970, -0.4599,
                              0.4695, -0.5199,  0.5187, -0.5597,
                              0.5903, -0.6246,  0.6656, -0.6942,
                              0.7086, -0.7325,  0.7348, -0.7546,
                              0.7702, -0.7864,  0.8083, -0.8210,
                              0.8439, -0.8529,  0.0000,  0.0000,
                              0.0000,  0.0000,  0.0000,  0.0000};

  // v parameter for the calculation of Y [1/bar]
  const Numeric V[n_lines] ={  0.0079, -0.0978,  0.0844, -0.1273,
                          0.0699, -0.0776,  0.2309, -0.2825, 
                          0.0436, -0.0584,  0.6056, -0.6619, 
                          0.6451, -0.6759,  0.6547, -0.6675,
                          0.6135, -0.6139,  0.2952, -0.2895, 
                          0.2654, -0.2590,  0.3750, -0.3680, 
                          0.5085, -0.5002,  0.6206, -0.6091,
                          0.6526, -0.6393,  0.6640, -0.6475,
                          0.6729, -0.6545,  0.0000,  0.0000,
                          0.0000,  0.0000,  0.0000,  0.0000};

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model 
  // (P. W. Rosenkranz, Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993):
  const Numeric CC_PWR93 = 1.00000;
  const Numeric CL_PWR93 = 1.00000;
  const Numeric CW_PWR93 = 1.00000;
  const Numeric CO_PWR93 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  if ( model == "Rosenkranz" )
    {
      CC = CC_PWR93;
      CL = CL_PWR93;
      CW = CW_PWR93;
      CO = CO_PWR93;
    }
  else if ( model == "RosenkranzLines" )
    {
      CC = 0.000;
      CL = CL_PWR93;
      CW = CW_PWR93;
      CO = CO_PWR93;
    }
  else if ( model == "RosenkranzContinuum" )
    {
      CC = CC_PWR93;
      CL = 0.000;
      CW = 0.000;
      CO = 0.000;
    }
  else if ( model == "RosenkranzNoCoupling" )
    {
      CC = CC_PWR93;
      CL = CL_PWR93;
      CW = CW_PWR93;
      CO = 0.000;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CL = CLin;
      CW = CWin;
      CO = COin;
    }
  else
    {
      ostringstream os;
      os << "O2-PWR93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'Rosenkranz', 'RosenkranzLines', RosenkranzContinuum, " 
         << "'RosenkranzNoCoupling', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "O2-PWR93: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CL = " << CL << "\n"
        << " CW = " << CW << "\n"
        << " CO = " << CO << "\n";
  
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      // check if O2-VMR will cause an underflow due to division by zero:
      if (vmr[i] < VMRCalcLimit)
        {
          ostringstream os;
          os << "ERROR: PWR93 O2 full absorption model has detected a O2 volume mixing ratio of " 
             << vmr[i] << " which is below the threshold of " << VMRCalcLimit << ".\n"
             << "Therefore no calculation is performed.\n";
          throw runtime_error(os.str());
          return;
        }
      // relative inverse temperature [1]
      Numeric TH  = 300.000 / t_abs[i];
      Numeric TH1 = TH-1.000;
      Numeric B   = pow(TH, X);
      // partial pressure of H2O and dry air [hPa]
      Numeric PRESWV = Pa_to_hPa * p_abs[i] * vmrh2o[i];
      Numeric PRESDA = Pa_to_hPa * p_abs[i] * (1.000 - vmrh2o[i]);
      Numeric DEN    = 0.001*(PRESDA*B + 1.1*PRESWV*TH);
      Numeric DFNR   = WB300*DEN;
      // initial O2 line absorption at position ff 
      Numeric O2ABS = 0.000;
      Numeric SUM   = 0.000;
      
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
        {
          // input frequency in [GHz]
          Numeric ff  = f_mono[s] * Hz_to_GHz; 
          // continuum absorption [Neper/km]
          SUM = CC * 1.6e-17 * ff * ff * DFNR / ( TH*(ff*ff + DFNR*DFNR) );
          
          // Loop over Rosnekranz '93 spectral line frequency:
          for ( Index l=0; l<n_lines; ++l )
            {
              Numeric DF   = CW * W300[l] * DEN;
              Numeric Y    = CO * 0.001 * 0.01 * p_abs[i] * B * ( Y300[l] + V[l]*TH1 );
              Numeric STR  = CL * S300[l] * exp(-BE[l] * TH1);
              Numeric SF1  = ( DF + (ff-F[l])*Y ) / ( (ff-F[l])*(ff-F[l]) + DF*DF );
              Numeric SF2  = ( DF - (ff+F[l])*Y ) / ( (ff+F[l])*(ff+F[l]) + DF*DF );
              SUM += STR * (SF1+SF2) * (ff/F[l]) * (ff/F[l]);
            }
          // O2 absorption [Neper/km] 
          O2ABS = 0.5034e12 * SUM * PRESDA * pow(TH, 3.0) / PI;
          // unit conversion x Nepers/km = y 1/m  --->  y = x * 1.000e-3 
          // xsec [1/m]  1.000e-3
          xsec(s,i) += 1.000e-3 * O2ABS / vmr[i];
        }
    }
  return;
}
//
// #################################################################################
//
// MPM93 O2 continuum:
// see publication side of National Telecommunications and Information Administration
//   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
// and ftp side for downloading the MPM93 original source code:
//   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            O2-continuum according to MPM93 [1/m]
   \param    S0in           O2-continuum strength [1/Pa]
   \param    G0in           O2-continuum width [Hz/Pa]
   \param    xS0in          O2-continuum strength temperature exponent [1]
   \param    xG0in          O2-continuum width temperature exponent    [1]
   \param    model          allows user defined input parameter set 
                            (S0in, G0in, xS0in, and xG0in)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid        [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    h2o_abs        H2O volume mixing ratio profile [1]
   \param    vmr            O2 volume mixing ratio profile  [1]

   \note     Except for  model 'user' the input parameters S0in, G0in, xS0in, and xG0in
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93_O2_continuum( MatrixView          xsec,
                         const Numeric       S0in,         // model parameter
                         const Numeric       G0in,         // model parameter
                         const Numeric       XS0in,        // model parameter
                         const Numeric       XG0in,        // model parameter
                         const String&       model,
                         ConstVectorView     f_mono,
                         ConstVectorView     p_abs,
                         ConstVectorView     t_abs,
                         ConstVectorView     h2o_abs,
                         ConstVectorView     vmr         )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric S0_MPM93  =  6.140e-13/0.20946; // line strength/VMR-O2 [1/Pa] 
  const Numeric G0_MPM93  =  0.560e4;           // line width [Hz/Pa]
  const Numeric XS0_MPM93 =  2.000;             // temperature dependence of line strength
  const Numeric XG0_MPM93 =  0.800;             // temperature dependence of line width
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates parameters!!):
  Numeric S0, G0, XS0, XG0;
  if ( model == "MPM93" )
    {
      S0  = S0_MPM93;
      G0  = G0_MPM93;
      XS0 = XS0_MPM93;
      XG0 = XG0_MPM93;
    }
  else if ( model == "user" )
    {
      S0  = S0in;
      G0  = G0in;
      XS0 = XS0in;
      XG0 = XG0in;
    }
  else
    {
      ostringstream os;
      os << "O2-SelfContMPM93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM93' and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "O2-SelfContMPM93: (model=" << model << ") parameter values in use:\n" 
        << " S0  = " << S0 <<  "\n"
        << " G0  = " << G0 <<  "\n"
        << " XS0 = " << XS0 << "\n"
        << " XG0 = " << XG0 << "\n";
  
  
  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  

  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      if (vmr[i] < VMRCalcLimit) // make sure that division by zero is excluded
        {
          ostringstream os;
          os << "ERROR: MPM93 O2 continuum absorption model has detected a O2 volume mixing ratio of " 
             << vmr[i] << " which is below the threshold of " << VMRCalcLimit << ".\n"
             << "Therefore no calculation is performed.\n";
          throw runtime_error(os.str());
          return;
        }
      Numeric th       = 300.0 / t_abs[i]; // Theta
      // S0 from the input has to be converted to unit 1/hPa --> * 1.0e-2
      Numeric strength =  S0 * p_abs[i] * (1.0000 - h2o_abs[i]) * pow( th, XS0 );
      // G0 from the input has to be converted to unit GHz/hPa --> * 1.0e-7
      Numeric gamma    =  G0 * p_abs[i] * pow( th, XG0 ); // Hz
          
      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
        {
          // the vmr of O2 will be multiplied at the stage of absorption calculation:
          // abs / vmr * xsec.
          xsec(s,i) +=  (4.0 * PI / SPEED_OF_LIGHT) *               // unit factor [1/(m*Hz)] 
                         strength *                                 // strength    [1]
                        ( pow( f_mono[s], 2) * gamma /              // line shape  [Hz]
                          ( pow( f_mono[s], 2) + pow( gamma, 2) ) );
        }
    }
  return;
}
//
// #################################################################################
//
//   3) O2-air : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
//               "Atmospheric Remote Sensing by Microwave Radiometry",
//               John Wiley & Sons, Inc., 1993. Also stated in 
//               Liebe et al. JQSRT, Vol 48, Nr 5/6, pp. 629-643, 1992.
//               Default continuum parameters are  C=1.6E-17*10E-9,  x=0.8
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            O2-continuum according to Rosenkranz 1993 [1/m]
   \param    S0in           line strength [K²/(Hz*Pa*m)]
   \param    G0in           line width [Hz/Pa]
   \param    XS0in          line strength temperature exponent  [1]
   \param    XG0in          line widths temperature exponent    [1]
   \param    model          allows user defined input parameter set 
                            (S0in, G0in, XS0in, and XG0in)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid        [Hz]
   \param    p_abs          predefined pressure grid         [Pa]
   \param    t_abs          predefined temperature grid      [K] 
   \param    h2o_abs        H2O volume mixing ratio profile  [1]
   \param    vmr            O2 volume mixing ratio profile   [1]

   \note     Except for  model 'user' the input parameters S0in, G0in, XS0in, and XG0in
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void Rosenkranz_O2_continuum( MatrixView        xsec,
                              const Numeric     S0in,          // model parameter
                              const Numeric     G0in,         // model parameter
                              const Numeric     XS0in,        // model parameter
                              const Numeric     XG0in,        // model parameter
                              const String&     model,
                              ConstVectorView   f_mono,
                              ConstVectorView   p_abs,
                              ConstVectorView   t_abs,
                              ConstVectorView   h2o_abs,
                              ConstVectorView   vmr      )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric S0_PWR93   =  1.108e-14; // [K²/(Hz*Pa*m)] line strength
  const Numeric G0_PWR93  =  5600.000;  // line width [Hz/Pa]
  const Numeric XS0_PWR93 =  2.000;    // temperature dependence of line strength
  const Numeric XG0_PWR93 =  0.800;    // temperature dependence of line width
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates values!!):
  Numeric S0, G0, XS0, XG0;
  if ( model == "Rosenkranz" )
    {
      S0  = S0_PWR93;
      G0  = G0_PWR93;
      XS0 = XS0_PWR93;
      XG0 = XG0_PWR93;
    }
  else if ( model == "user" )
    {
      S0  = S0in;
      G0  = G0in;
      XS0 = XS0in;
      XG0 = XG0in;
    }
  else
    {
      ostringstream os;
      os << "O2-SelfContPWR93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'Rosenkranz' and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "O2-SelfContPWR93: (model=" << model << ") parameter values in use:\n" 
        << " S0  = " << S0 <<  "\n"
        << " G0  = " << G0 <<  "\n"
        << " XS0 = " << XS0 << "\n"
        << " XG0 = " << XG0 << "\n";
  

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // loop over all pressure levels:
  for ( Index i=0; i<n_p; ++i )
    {
      Numeric TH = 300.00 / t_abs[i];         // relative temperature  [1]
      
      Numeric ph2o  = p_abs[i] * h2o_abs[i];  // water vapor partial pressure [Pa]
      Numeric pdry  = p_abs[i] - ph2o;        // dry air partial pressure     [Pa]


      // pseudo broadening term [Hz]
      Numeric gamma = G0 * (pdry * pow( TH, XG0 ) + 1.100 * ph2o * TH); 

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
        {
          // division by vmr of O2 is necessary because of the absorption calculation
          // abs = vmr * xsec.
          xsec(s,i) += S0 * p_abs[i] * gamma * pow( f_mono[s], 2 ) / 
                       ( pow( t_abs[i], XS0 ) * ( pow( f_mono[s], 2 ) + pow( gamma, 2 ) ) );
        }
    }
}
//
//
// #################################################################################
//
/**
   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            O2-continuum according to Rosenkranz 1993 [1/m]
   \param    S0in           line strength                             [1/(Hz*Pa*m)]
   \param    G0in           line width                                [Hz/Pa]
   \param    G0Ain          dry air broadening parameter              [1]
   \param    G0Bin          water vapor broadening parameter          [1]
   \param    XG0din         line strength temperature exponent        [1]
   \param    XG0win         line widths temperature exponent          [1]
   \param    model          allows user defined input parameter set 
                            (S0in, G0in, XS0in, and XG0in)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid                 [Hz]
   \param    p_abs          predefined pressure grid                  [Pa]
   \param    t_abs          predefined temperature grid               [K] 
   \param    h2o_abs        H2O volume mixing ratio profile           [1]
   \param    vmr            O2 volume mixing ratio profile            [1]

   \note     Except for  model 'user' the input parameters S0in, G0in, XS0in, and XG0in
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.<br>
             <br>
             Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void Standard_O2_continuum( MatrixView        xsec,         // cross section
                            const Numeric     Cin,          // model parameter
                            const Numeric     G0in,         // model parameter
                            const Numeric     G0Ain,        // model parameter
                            const Numeric     G0Bin,        // model parameter
                            const Numeric     XG0din,       // model parameter
                            const Numeric     XG0win,       // model parameter
                            const String&     model,  // model parameter
                            ConstVectorView   f_mono,       // frequency grid
                            ConstVectorView   p_abs,        // P_tot grid
                            ConstVectorView   t_abs,        // T grid
                            ConstVectorView   h2o_abs,      // VMR H2O profile
                            ConstVectorView   vmr       )   // VMR O2  profile
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry",
  // John Wiley & Sons, Inc., 1993. Also stated in 
  // Liebe et al. JQSRT, Vol 48, Nr 5/6, pp. 629-643, 1992.
  const Numeric C_PWR93    = 1.22e-19;  // line strength/VMR [1/m*1/Hz*1/Pa] 
  const Numeric G0_PWR93   = 5600.000;  // line width [Hz/Pa]
  const Numeric G0A_PWR93  =    1.000;  // line width [1]
  const Numeric G0B_PWR93  =    1.100;  // line width [1]
  const Numeric XG0d_PWR93 =    0.800;  // temperature dependence of line width [1]
  const Numeric XG0w_PWR93 =    1.000;  // temperature dependence of line width [1]
  //
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric C_MPM93    = 1.23e-19; // line strength/VMR [1/m*1/Hz*1/Pa] 
  const Numeric G0_MPM93   = 5600.000; // line width [Hz/Pa]
  const Numeric G0A_MPM93  =    1.000; // line width [1]
  const Numeric G0B_MPM93  =    1.000; // line width [1]
  const Numeric XG0d_MPM93 =    0.800; // temperature dependence of line strength [1]
  const Numeric XG0w_MPM93 =    0.800; // temperature dependence of line width    [1]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates values!!):
  Numeric C, G0, G0A, G0B, XG0d, XG0w;
  if ( model == "Rosenkranz" )
    {
      C    = C_PWR93;
      G0   = G0_PWR93;
      G0A  = G0A_PWR93;
      G0B  = G0B_PWR93;
      XG0d = XG0d_PWR93;
      XG0w = XG0w_PWR93;
    }
  else if ( model == "MPM93" )
    {
      C    = C_MPM93;
      G0   = G0_MPM93;
      G0A  = G0A_MPM93;
      G0B  = G0B_MPM93;
      XG0d = XG0d_MPM93;
      XG0w = XG0w_MPM93;
    }
  else if ( model == "user" )
    {
      C    = Cin;
      G0   = G0in;
      G0A  = G0Ain;
      G0B  = G0Bin;
      XG0d = XG0din;
      XG0w = XG0win;
    }
  else
    {
      ostringstream os;
      os << "O2-GenerealCont: ERROR! Wrong model values given.\n"
         << "Valid models are: 'Rosenkranz', 'MPM93' and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "O2-GeneralCont: (model=" << model << ") parameter values in use:\n" 
        << " C    = " << C <<  "\n"
        << " G0   = " << G0 <<  "\n"
        << " G0A  = " << G0A <<  "\n"
        << " G0B  = " << G0B <<  "\n"
        << " XG0d = " << XG0d << "\n"
        << " XG0w = " << XG0w << "\n";
  

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // loop over all pressure levels:
  for ( Index i=0; i<n_p; ++i )
    {
      Numeric TH = 300.00 / t_abs[i];         // relative temperature  [1]
      
      Numeric ph2o = p_abs[i] * h2o_abs[i];  // water vapor partial pressure [Pa]
      Numeric pdry = p_abs[i] - ph2o;        // dry air partial pressure     [Pa]


      // pseudo broadening term [Hz]
      Numeric gamma = G0 * (G0A * pdry * pow( TH, XG0d ) + G0B * ph2o * pow( TH, XG0w )); 

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
        {
          // division by vmr of O2 is necessary because of the absorption calculation
          // abs = vmr * xsec.
          xsec(s,i) += C * p_abs[i] * pow( TH, 2 ) * 
                       ( gamma * pow( f_mono[s], 2 ) / 
                         ( pow( f_mono[s], 2 ) + pow( gamma, 2 ) ) );
        }
    }
}
//
// #################################################################################
// ################################ NITROGEN MODELS ################################
// #################################################################################
//
// MPM93 N2 continuum:
// see publication side of National Telecommunications and Information Administration
//   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
// and ftp side for downloading the MPM93 original source code:
//   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            N2-continuum according to MPM93 [1/m]
   \param    Cin            continuum strength [ppm/GHz]
   \param    Gin            width parameter [Hz/Pa]
   \param    xTin           continuum strength temperature exponent [1]
   \param    xfin           continuum frequency exponent [1]
   \param    model          allows user defined input parameter set 
                            (Cin, Gin, xTin, and xfin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    h2o_abs        H2O volume mixing ratio profile      [1]
   \param    vmr            N2 volume mixing ratio profile       [1]

   \note     Except for  model 'user' the input parameters Cin, Gin, xTin, and xfin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93_N2_continuum( MatrixView          xsec,
                         const Numeric       Cin,
                         const Numeric       Gin,
                         const Numeric       xTin,
                         const Numeric       xfin,
                         const String&       model,
                         ConstVectorView     f_mono,
                         ConstVectorView     p_abs,
                         ConstVectorView     t_abs,
                         ConstVectorView     h2o_abs,
                         ConstVectorView     vmr         )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 H2O continuum model 
  // (AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric xT_MPM93  =  3.500;             // temperature exponent [1]
  const Numeric xf_MPM93  =  1.500;             // frequency exponent [1]
  const Numeric gxf_MPM93 =  9.000*xf_MPM93;    // needed for the unit conversion of G_MPM93
  const Numeric S_MPM93   =  2.296e-31;         // line strength  [1/Pa² * 1/Hz]
  const Numeric G_MPM93   =  1.930e-5*pow(10.000, -gxf_MPM93); // frequency factor [1/Hz^xf]
  // ---------------------------------------------------------------------------------------
  
  // select the parameter set (!!model dominates values!!):
  Numeric S0, G0, xf, xT, gxf;
  if ( model == "MPM93" )
    {
      S0  = S_MPM93;
      G0  = G_MPM93;
      xT  = xT_MPM93;
      xf  = xf_MPM93;
      gxf = gxf_MPM93;
    }
  else if ( model == "user" )
    {
      S0  = Cin;
      G0  = Gin;
      xT  = xTin;
      xf  = xfin;
      gxf = 9.000*xf;
    }
  else
    {
      ostringstream os;
      os << "N2-SelfContMPM93 : ERROR! Wrong model values given.\n"
         << "allowed models are: 'MPM93', 'user'" << '\n';
      throw runtime_error(os.str());
    }
  
  out2  << "N2-SelfContMPM93: (model=" << model << ") parameter values in use:\n" 
        << " S0 = " << S0 << "\n"
        << " G0 = " << G0 << "\n"
        << " xT = " << xT << "\n"
        << " xf = " << xf << "\n";
  
  // unit conversion internally:
  //const Numeric S0unitconv = 1.000e+13;  // x [1/(hPa²*GHz)] => y [1/(pa²*Hz)]
  //const Numeric G0unitconv = pow(10.000, gxf);

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  Numeric fac = 4.0 * PI / SPEED_OF_LIGHT;  //  = 4 * pi / c
  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      Numeric th = 300.0 / t_abs[i];
      Numeric strength =  S0 * 
                          pow( (p_abs[i] * (1.0000 - h2o_abs[i])), 2 ) * 
                          pow( th, xT );

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
        {

          //      Numeric f = f_mono[s] * Hz_to_GHz; // frequency in GHz
          // FIXME: SAB+TKS: Commented this out, f_mono is used directly

          // the vmr of N2 will be multiplied at the stage of absorption calculation:
          // abs / vmr * xsec.
          xsec(s,i) += fac * strength *                              // strength
                       pow(f_mono[s],2) /                      // frequency dependence
                       ( 1.000 + G0 * pow( f_mono[s], xf) ) *  
                        vmr[i];                                // N2 vmr
        }
    }
  return;
}
//
// #################################################################################
//
/** 
   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    xTin           continuum strength temperature exponent [1]
   \param    model          allows user defined input parameter set 
                            (Cin and xTin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid      [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid    [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters Cin and xTin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz, Chapter 2, pp 74, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void Rosenkranz_N2_self_continuum( MatrixView          xsec,
                                   const Numeric       Cin,
                                   const Numeric       xin,
                                   const String&       model,
                                   ConstVectorView     f_mono,
                                   ConstVectorView     p_abs,
                                   ConstVectorView     t_abs,
                                   ConstVectorView     vmr       )
{
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model (Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric C_PWR = 1.05e-38; // [1/(Pa²*Hz²*m)]
  const Numeric x_PWR = 3.55;     // [1]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates parameters!!):
  Numeric C, x;
   if ( model == "Rosenkranz" )
     {
       C = C_PWR;
       x = x_PWR;
     }
   else if ( model == "user" )
     {
       C = Cin;
       x = xin;
     }
   else
     {
       ostringstream os;
       os << "N2-SelfContPWR93: ERROR! Wrong model values given.\n"
          << "allowed models are: 'Rosenkranz', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out2  << "N2-SelfContPWR93: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";

   const Index n_p = p_abs.nelem();     // Number of pressure levels
   const Index n_f = f_mono.nelem();    // Number of frequencies

   // Check that dimensions of p_abs, t_abs, and vmr agree:
   assert ( n_p==t_abs.nelem() );
   assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( Index i=0; i<n_p; ++i )
    {
      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
        {
          // The second vmr of N2 will be multiplied at the stage of absorption 
          // calculation: abs = vmr * xsec.
          xsec(s,i) += C *                        // strength [1/(m*Hz²Pa²)] 
                       pow( f_mono[s], 2 ) *      // quadratic f dependence [Hz²]
                       pow( 300.0/t_abs[i], x ) * // free T dependence      [1]
                       pow( p_abs[i], 2 ) *       // quadratic p dependence [Pa²]
                       vmr[i];                    // second N2-VMR at the stage 
                                                  // of absorption calculation
        }
    }
}
//
// #################################################################################
//
// 4) N2-N2  : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
//    "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    xfin           continuum frequency exponent [1]
   \param    xtin           continuum strength temperature exponent [1]
   \param    xpin           continuum strength pressure exponent [1]
   \param    model          allows user defined input parameter set 
                            (Cin, xfin, xtin, and xpin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid      [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid    [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters Cin, xfin, xtin, and xpin
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: P. W. Rosenkranz, Chapter 2, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void Standard_N2_self_continuum( MatrixView          xsec,
                                 const Numeric       Cin,
                                 const Numeric       xfin,
                                 const Numeric       xtin,
                                 const Numeric       xpin,
                                 const String&       model,
                                 ConstVectorView     f_mono,
                                 ConstVectorView     p_abs,
                                 ConstVectorView     t_abs,
                                 ConstVectorView     vmr         )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model, Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric C_GM  = 1.05e-38; // [1/(Pa²*Hz²*m)]
  const Numeric xf_GM = 2.00;     // [1]
  const Numeric xt_GM = 3.55;     // [1]
  const Numeric xp_GM = 2.00;     // [1]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates over values!!):
  Numeric C, xt, xf, xp;
   if ( model == "Rosenkranz" )
     {
       C  = C_GM;
       xt = xt_GM;
       xf = xf_GM;
       xp = xp_GM;
     }
   else if ( model == "user" )
     {
       C  = Cin;
       xt = xtin;
       xf = xfin;
       xp = xpin;
     }
   else
     {
       ostringstream os;
       os << "N2-SelfContStandardType: ERROR! Wrong model values given.\n"
          << "allowed models are: 'Rosenkranz', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out2  << "N2-SelfContStandardType: (model=" << model << ") parameter values in use:\n" 
         << " C  = " << C  << "\n"
         << " xt = " << xt << "\n"
         << " xf = " << xf << "\n"
         << " xp = " << xp << "\n";


  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( Index i=0; i<n_p; ++i )
    {
      //cout << "vmr[" << i << "]= " << vmr[i] << "\n";
      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
        {
          // The second N2-VMR will be multiplied at the stage of absorption 
          // calculation: abs = vmr * xsec.
          xsec(s,i) += C                            * // strength [1/(m*Hz²Pa²)]  
                       pow( (300.00/t_abs[i]), xt ) * // T dependence        [1] 
                       pow( f_mono[s], xf )         * // f dependence    [Hz^xt]  
                       pow( p_abs[i], xp )          * // p dependence    [Pa^xp]
                       pow( vmr[i], (xp-1) );         // last N2-VMR at the stage 
                                                      // of absorption calculation
        }
    }
}
//
// #################################################################################
// ############################## CARBON DIOXIDE MODELS ############################
// #################################################################################
//
/** 
   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            CO2-CO2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    xin            continuum temperature exponent [1]
   \param    model          allows user defined input parameter set 
                            (Cin and xin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid      [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid    [K] 
   \param    vmr            CO2 volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters Cin and xin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference:  P. W. Rosenkranz, Chapter 2, pp 74, pp 74, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void Rosenkranz_CO2_self_continuum( MatrixView          xsec,
                                    const Numeric       Cin,
                                    const Numeric       xin,
                                    const String&       model,
                                    ConstVectorView     f_mono,
                                    ConstVectorView     p_abs,
                                    ConstVectorView     t_abs,
                                    ConstVectorView     vmr      )
{
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric C_PWR = 7.43e-37; // [ 1/(Pa²*Hz²*m) ]
  const Numeric x_PWR = 5.08;     // [ 1 ]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates values!!):
  Numeric C, x;
  if ( model == "Rosenkranz" )
    {
      C = C_PWR;
      x = x_PWR;
    }
  else if ( model == "user" )
    {
      C = Cin;
      x = xin;
    }
  else
    {
      ostringstream os;
      os << "CO2-SelfContPWR93 : ERROR! Wrong model values given.\n"
         << "allowed models are: 'Rosenkranz', 'user'" << "\n";
      throw runtime_error(os.str());
    }
  
  out2  << "CO2-SelfContPWR93: (model=" << model << ") parameter values in use:\n" 
        << " C = " << C << "\n"
        << " x = " << x << "\n";

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( Index i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The second vmr of CO2 will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy =
        C * pow( 300./t_abs[i], x ) * pow( p_abs[i], 2 ) * vmr[i];

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
        {
          xsec(s,i) += dummy * pow( f_mono[s], 2 );
        }
    }
}
//
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            CO2-N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    xin            continuum temperature exponent [1]
   \param    model          allows user defined input parameter set 
                            (Cin and xin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid        [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    n2_abs         N2 volume mixing ratio profile  [1]
   \param    vmr            CO2 volume mixing ratio profile [1]

   \note     Except for  model 'user' the input parameters Cin and xin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'Rosenkranz' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference:  P. W. Rosenkranz, Chapter 2, pp 74, in M. A. Janssen, <br>
             <I>Atmospheric Remote Sensing by Microwave Radiometry</i>,<br>
             John Wiley & Sons, Inc., 1993.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void Rosenkranz_CO2_foreign_continuum( MatrixView          xsec,
                                       const Numeric       Cin,
                                       const Numeric       xin,
                                       const String&       model,
                                       ConstVectorView     f_mono,
                                       ConstVectorView     p_abs,
                                       ConstVectorView     t_abs,
                                       ConstVectorView     n2_abs,
                                       ConstVectorView     vmr   )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric C_PWR = 2.71e-37; // default: 2.71*10^-37 1/(Pa²*Hz²*m)
  const Numeric x_PWR = 4.7;      // default: 4.7
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates values!!):
  Numeric C, x;
  if ( model == "Rosenkranz" )
    {
      C = C_PWR;
      x = x_PWR;
    }
  else if ( model == "user" )
    {
      C = Cin;
      x = xin;
    }
  else
    {
      ostringstream os;
      os << "CO2-ForeignContPWR93: ERROR! Wrong model values given.\n"
         << "allowed models are: 'Rosenkranz', 'user'" << "\n";
      throw runtime_error(os.str());
    }
  
  out2  << "CO2-ForeignContPWR93: (model=" << model << ") parameter values in use:\n" 
        << " C = " << C << "\n"
        << " x = " << x << "\n";

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The vmr of CO2 will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy = C * pow( 300./t_abs[i], x ) * p_abs[i] * p_abs[i] * n2_abs[i];

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
        {
          xsec(s,i) += dummy * pow( f_mono[s], 2 );
        }
    }
}
//
// #################################################################################
// ################################### CLOUD MODELS ################################
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            water clouds according to MPM93 [1/m]
   \param    CCin           scaling parameter of the calculated cross section [1]
   \param    CGin           scaling parameter of the first relaxation frequency 
                            (gamma_1, see page 3-6 in the reference) [1]
   \param    CEin           scaling parameter of the first permittivity component 
                            (epsilon_0, see page 3-6 in the reference) [1]
   \param    model          allows user defined input parameter 
                            (CCin, CGin, CEin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid        [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            suspended water droplet density profile (valid range: 0-0.001) [kg/m³]

   \note     Except for  model 'user' the input parameters CCin, CGin, and CEin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference:  H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93WaterDropletAbs( MatrixView         xsec,
                           const Numeric      CCin,   // input parameter
                           const Numeric      CGin,   // input parameter
                           const Numeric      CEin,   // input parameter
                           const String&      model, // model
                           ConstVectorView    f_mono, // frequency vector
                           ConstVectorView    p_abs,  // pressure vector
                           ConstVectorView    t_abs,  // temperature vector
                           ConstVectorView    vmr)    // suspended water droplet density vector
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21)
  const Numeric CC_MPM93 = 1.00000;
  const Numeric CG_MPM93 = 1.00000;
  const Numeric CE_MPM93 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CG, CE;
  if ( model == "MPM93" )
    {
      CC = CC_MPM93;
      CG = CG_MPM93;
      CE = CE_MPM93;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CG = CGin;
      CE = CEin;
    }
  else
    {
      ostringstream os;
      os << "liquidcloud-MPM93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM93' and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "liquidcloud-MPM93: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CG = " << CG << "\n"
        << " CE = " << CE << "\n";
  
  
  const Numeric m = 1.00e3; // specific weight of the droplet,  fixed value:  1.00e3 kg/m³
  const Numeric low_lim_den  =  0.000;   // lower limit of suspended droplet particle density vector [kg/m³]
  const Numeric high_lim_den = 10.00e-3; // lower limit of suspended droplet particle density vector [kg/m³]

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      // water vapor saturation pressure over liquid water [Pa]
      // Numeric es       = WVSatPressureLiquidWater(t_abs[i]);
      // water vapor partial pressure [Pa]
      // Numeric e        = p_abs[i] * vmr[i];      
      // relative humidity [1]
      // Numeric RH       = e / es;
  
      // Check limits of suspended water droplet density ("vmr") [kg/m³]
      if ( (vmr[i] > low_lim_den) && (vmr[i] < high_lim_den) ) 
        {
          // relative inverse temperature [1]
          Numeric theta    = 300.000 / t_abs[i];
          // relaxation frequencies [GHz]
          Numeric gamma1   = CG * 20.20 - 146.40*(theta-1.000) + 316.00*(theta-1.000)*(theta-1.000);
          // Numeric gamma1  = 20.1 * exp( 7.88 * theta ); // see Liebe et al. IJIMW, 1992, p667, Eq. (2b)
          Numeric gamma2   = 39.80 * gamma1; 
          // static and high-frequency permittivities
          Numeric epsilon0 = CE * 103.30 * (theta-1.000) + 77.66;
          Numeric epsilon1 = 0.0671 * epsilon0;
          Numeric epsilon2 = 3.52;
          
          // Loop frequency:
          for ( Index s=0; s<n_f; ++s )
            {
              // real part of the complex permittivity of water (double-debye model)
              Numeric Reepsilon  = epsilon0 - 
                pow((f_mono[s]*Hz_to_GHz),2) *
                ( ((epsilon0-epsilon1)/
                   (pow((f_mono[s]*Hz_to_GHz),2) + pow(gamma1,2))) + 
                  ((epsilon1-epsilon2)/
                   (pow((f_mono[s]*Hz_to_GHz),2) + pow(gamma2,2))) );
              // imaginary part of the complex permittivity of water (double-debye model)
              Numeric Imepsilon  = (f_mono[s]*Hz_to_GHz) *
                ( (gamma1*(epsilon0-epsilon1)/
                   (pow((f_mono[s]*Hz_to_GHz),2) + pow(gamma1,2))) + 
                  (gamma2*(epsilon1-epsilon2)/
                   (pow((f_mono[s]*Hz_to_GHz),2) + pow(gamma2,2))) );
              // the imaginary part of the complex refractivity of suspended liquid water particle [ppm]
              // In MPM93 w is in g/m³ and m is in g/cm³. Because of the units used in arts,
              // a factor of 1.000e6 must be multiplied with the ratio (w/m):
              // MPM93: (w/m)_MPM93  in   (g/m³)/(g/cm³)
              // arts:  (w/m)_arts   in  (kg/m³)/(kg/m³)
              // =====> (w/m)_MPM93   =   1.0e6 * (w/m)_arts
              // the factor of 1.0e6 is included below in xsec calculation.
              Numeric ImNw = 1.500 / m * 
                ( 3.000 * Imepsilon / ( pow((Reepsilon+2.000),2) + pow(Imepsilon,2) ) );
              // liquid water particle absorption cross section [1/m]
              // The vmr of H2O will be multiplied at the stage of absorption 
              // calculation: abs = vmr * xsec.
              // xsec = abs/vmr [1/m] but MPM93 is in [dB/km] --> conversion necessary
              xsec(s,i) += CC * 1.000e6 * dB_km_to_1_m * 0.1820 * (f_mono[s]*Hz_to_GHz) * ImNw;
            }
        } else
          {
            if ( (vmr[i] < low_lim_den) || (vmr[i] > high_lim_den) ) 
              {
                ostringstream os;
                os << "ERROR in MPM93WaterDropletAbs:\n"
                   << " suspended water droplet density (valid range 0.00-10.00e-3 kg/m3):" << vmr[i] << "\n"
                   << " ==> no calculation performed!\n";
                throw runtime_error(os.str());
              }
          }
    }

}
//
// #################################################################################
//
/** 

   \param    xsec Output:          cross section (absorption/volume mixing ratio) of 
                            ice clouds according to MPM93 [1/m]
   \param    CCin           scaling parameter of the calculated cross section [1]
   \param    CAin           scaling parameter of the permittivity function a
                            (see page 3-6 in the reference) [1]
   \param    CBin           scaling parameter of the permittivity function b
                            (see page 3-6 in the reference) [1]
   \param    model          allows user defined input parameter 
                            (CCin, CAin, CBin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid        [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            suspended water droplet density profile (valid range: 0-0.001) [kg/m³]

   \note     Except for  model 'user' the input parameters CCin, CAin, and CBin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference:  H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21.

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93IceCrystalAbs( MatrixView        xsec,
                         const Numeric     CCin,   // input parameter
                         const Numeric     CAin,   // input parameter
                         const Numeric     CBin,   // input parameter
                         const String&     model, // model
                         ConstVectorView   f_mono, // frequency vector
                         ConstVectorView   p_abs,  // pressure vector
                         ConstVectorView   t_abs,  // temperature vector
                         ConstVectorView   vmr   ) // suspended ice particle density vector, 
                                                   // valid range: 0-10.0e-3 kg/m³
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21)
  const Numeric CC_MPM93 = 1.00000;
  const Numeric CA_MPM93 = 1.00000;
  const Numeric CB_MPM93 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CA, CB;
  if ( model == "MPM93" )
    {
      CC = CC_MPM93;
      CA = CA_MPM93;
      CB = CB_MPM93;
    }
  else if ( model == "user" )
    {
      CC = CCin;
      CA = CAin;
      CB = CBin;
    }
  else
    {
      ostringstream os;
      os << "icecloud-MPM93: ERROR! Wrong model values given.\n"
         << "Valid models are: 'MPM93' and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out2  << "icecloud-MPM93: (model=" << model << ") parameter values in use:\n" 
        << " CC = " << CC << "\n"
        << " CA = " << CA << "\n"
        << " CB = " << CB << "\n";
  
  
  const Numeric m = 0.916e3;  // specific weight of ice particles,  fixed value:   0.916e3 kg/m³
  const Numeric low_lim_den  =  0.000;   // lower limit of suspended ice particle density vector [kg/m³]
  const Numeric high_lim_den = 10.00e-3; // lower limit of suspended ice particle density vector [kg/m³]

  const Index n_p = p_abs.nelem();      // Number of pressure levels
  const Index n_f = f_mono.nelem();     // Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  


  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      // water vapor saturation pressure over ice [Pa]
      // Numeric es = WVSatPressureIce(t_abs[i]);
      // water vapor partial pressure [Pa]
      // Numeric e  = p_abs[i] * vmr[i];
      // relative humidity [1]
      // Numeric RH = e / es;
  
      // Check limits of suspended water ice crystal density ("vmr") [kg/m³]
      if ( (vmr[i] > low_lim_den) && (vmr[i] < high_lim_den) ) 
        { 
          // relative inverse temperature [1]
          Numeric theta = 300.000 / t_abs[i];   
          // inverse frequency T-dependency function [Hz]
          Numeric ai = CA * (62.000 * theta - 11.600) * exp(-22.100 * (theta-1.000)) * 1.000e-4;
          // linear frequency T-dependency function [1/Hz]
          Numeric bi = CB * 0.542e-6 * 
                       ( -24.17 + (116.79/theta) + pow((theta/(theta-0.9927)),2) );
              
          // Loop frequency:
          for ( Index s=0; s<n_f; ++s )
            {
              // real part of the complex permittivity of ice
              Numeric Reepsilon  = 3.15;
              // imaginary part of the complex permittivity of water
              Numeric Imepsilon  = ( ( ai/(f_mono[s]*Hz_to_GHz) ) +
                                     ( bi*(f_mono[s]*Hz_to_GHz) ) );
              // the imaginary part of the complex refractivity of suspended ice particles.
              // In MPM93 w is in g/m³ and m is in g/cm³. Because of the units used in arts,
              // a factor of 1.000e6 must be multiplied with the ratio (w/m):
              // MPM93: (w/m)_MPM93  in   (g/m³)/(g/cm³)
              // arts:  (w/m)_arts   in  (kg/m³)/(kg/m³)
              // =====> (w/m)_MPM93   =   1.0e6 * (w/m)_arts
              // the factor of 1.0e6 is included below in xsec calculation.
              Numeric ImNw = 1.500 / m * 
                    ( 3.000 * Imepsilon / ( pow((Reepsilon+2.000),2) + pow(Imepsilon,2) ) );
              // ice particle absorption cross section [1/m]
              // The vmr of H2O will be multiplied at the stage of absorption 
              // calculation: abs = vmr * xsec.
              // xsec = abs/vmr [1/m] but MPM93 is in [dB/km] --> conversion necessary
              xsec(s,i) += CC * 1.000e6 * dB_km_to_1_m * 0.1820 * (f_mono[s]*Hz_to_GHz) * ImNw;
            }
        } else
          {
            if ( (vmr[i] < low_lim_den) || (vmr[i] > high_lim_den) ) 
              {
                ostringstream os;
                os << "ERROR in MPM93IceCrystalAbs:\n"
                   << " suspended ice particle density (valid range: 0-10.0e-3 kg/m3):" << vmr[i] << "\n"
                   << " ==> no calculation performed!\n";
                throw runtime_error(os.str());
              }
          }
    }
  return;
}
//
// #################################################################################
// ################################# HELP FUNCTIONS ################################
// #################################################################################
//
/** 

   \param   MPMLineShapeFunction Output:  H2O-line shape function value     [1/Hz]  
   \param    gamma                 H2O-line width                    [Hz]
   \param    fl                    H2O-line central frequency        [Hz]
   \param    f                     frequency position of calculation [Hz]

   \note     This function calculates the line shape function of Van Vleck and Weisskopf
             with the factor (f/fl)¹. for the MPM pseudo continuum line.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

Numeric MPMLineShapeFunction( const Numeric gamma, 
                              const Numeric fl, 
                              const Numeric f)
{
  /*
    this routine calculates the line shape function of Van Vleck and Weisskopf
    with the factor (f/f_o)¹. for the MPM pseudo continuum line.

    creation  TKS, 4.11.00 

    input:   gamma   [Hz]    line width of line L
             fl      [Hz]    central frequency of line L
             f       [Hz]    frequency position of calculation
             
    output:  value   [1/Hz]  line shape function value at f for the line parameters
                              of line L 
              
   */

  double f_minus, f_plus ;           /* internal variables */
  double value;                      /* return value       */

  // line at fl
  f_minus = 1.000 / ((f-fl)*(f-fl) + gamma*gamma);

  // mirror line at -fl
  f_plus  = 1.000 / ((f+fl)*(f+fl) + gamma*gamma);

  // VVW line shape function value
  value = fabs(f/fl) * gamma * (f_minus + f_plus);
  
  return value;
}
//
// #################################################################################
//
/** 

   \param   MPMLineShapeO2Function Output:  O2-line shape function value         [1]  
   \param    gamma                   O2-line width                        [Hz]
   \param    fl                      H2O-line central frequency of the    [Hz]
   \param    f                       frequency position of calculation    [Hz]
   \param    delta                   O2-line mixing parameter             [1]

   \note     This function calculates the line shape function of Van Vleck and Weisskopf
             for O2 with line mixing.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

Numeric MPMLineShapeO2Function( const Numeric gamma, 
                                const Numeric fl, 
                                const Numeric f,
                                const Numeric delta)
{
  /*
    this routine calculates the line shape function of Van Vleck and Weisskopf
    for O2 with line mixing.

    creation  TKS, 14.07.01 

    input:   gamma   [GHz]    line width of line L
             fl      [GHz]    central frequency of line L
             f       [GHz]    frequency position of calculation
             delta   [1]     line mixing parameter
             
    output:  value   [1]     line shape function value at f for the line parameters
                             of line L 
              
   */

  double f_minus, f_plus ;           /* internal variables */
  double value;                      /* return value       */

  // line at fl
  f_minus = (gamma - delta * (fl-f)) / ((fl-f)*(fl-f) + gamma*gamma);

  // mirror line at -fl
  f_plus  = (gamma - delta * (fl+f)) / ((fl+f)*(fl+f) + gamma*gamma);

  // VVW line shape function value
  value = f * (f_minus + f_plus);
  
  return value;
}
//
// #################################################################################
//
/** 

   \param   WVSatPressureLiquidWater Output:     water vapor saturation pressure over liquid water [Pa]  
   \param    t                            temperature                                       [K]

   \note     This function calculates the water vapor saturation pressure 
             over liquid water according to the 
             <a href="http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html">Goff-Gratch equation</a>.<br>
             Other links:<br>
             <a href="http://hydro.iis.u-tokyo.ac.jp/~taikan/Publication/HP95/HP95.html">Global Atmospheric Water Balance</a><br>
             <a href="http://dss.ucar.edu/datasets/ds540.1/software/supl_soft/profs">Schlatter (profsc.fsl.noaa.gov)</a><br>

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

Numeric WVSatPressureLiquidWater(const Numeric t)
{

  // check of temperature range
  if (t < 0.000)
    {
      ostringstream os;
      os << "In WVSatPressureLiquidWater:\n"
         << "temperature negative: T=" << t <<"K \n";
      throw runtime_error(os.str());
    }

  //  COMPUTES SATURATION H2O VAPOR PRESSURE (OVER LIQUID)
  //  USING LIEBE'S APPROXIMATION (CORRECTED)
  //  input  : T in Kelvin
  //  output : es in Pa
  //  PWR 4/8/92
  /*
  Numeric TH       = 300.0 / t;
  Numeric es_PWR98 = 100.00 * 35.3 * exp(22.64*(1.-TH)) * pow(TH,5.0);
  */

  // MPM93 calculation
  Numeric theta    = 373.16 / t;
  Numeric exponent = ( -7.90298 * (theta-1.000) +
                        5.02808 * log10(theta) -
                        1.3816e-7 * ( pow( 10.00, (11.344*(1.00-(1.00/theta))) ) - 1.000 ) +
                        8.1328e-3 * ( pow( 10.00, (-3.49149*(theta-1.00))) - 1.000) +
                        log10(1013.246) );
  Numeric es_MPM93 = 100.000 * pow(10.00,exponent);

  return es_MPM93; // [Pa]
}
//
// #################################################################################
//
/** 

   \param   WVSatPressureIce Output:     water vapor saturation pressure over liquid water [Pa]  
   \param    t                    temperature                                       [K]

   \note     This function calculates the water vapor saturation pressure 
             over ice water according to the 
             <a href="http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html">Goff-Gratch equation</a>.
             Other links:<br>
             <a href="http://hydro.iis.u-tokyo.ac.jp/~taikan/Publication/HP95/HP95.html">Global Atmospheric Water Balance</a><br>
             <a href="http://dss.ucar.edu/datasets/ds540.1/software/supl_soft/profs">Schlatter (profsc.fsl.noaa.gov)</a><br>

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

Numeric WVSatPressureIce(const Numeric t)
{

  // check of temperature range
  if (t < 0.000)
    {
      ostringstream os;
      os << "In WVSatPressureIce:\n"
         << "temperature negative: T=" << t <<"K \n";
      throw runtime_error(os.str());
    }

  // MPM93 calculation
  Numeric theta    = 273.16 / t;
  Numeric exponent = (-9.09718  * (theta-1.000) -
                       3.56654  * log10(theta)  +
                       0.876793 * (1.000-(1.000/theta)) +
                       log10(6.1071) );

  Numeric es_MPM93 = 100.000 * pow(10.00,exponent);

  return es_MPM93;
}
//
// #################################################################################
// #################### CONTROL OF ADDITIONAL ABSORPTION MODEL #####################
// #################################################################################
//
//
/** 
    Calculates model absorption for one continuum or full model tag. 
    Note, that only one tag can be taken at a time. 

    Calculated is the absorption cross section, that means you have to
    multiply this with the VMR in order to get the absorption
    coefficient in units of 1/m.

    \param xsec Output:       Cross section of one continuum tag,<br> 
                       xsec = alpha / VMR  [1/m * 1]

    \param  name       The name of the model to calculate (derived from the tag name)
    \param  parameters model parameters, as defined in method
                       cont_description_parameters.
    \param  model      model, related to model parameters
    \param  f_mono     Frequency grid [Hz]
    \param  p_abs      Pressure grid [Pa]
    \param  t_abs      Temperatures associated with the pressure grid, p_abs [K]
    \param  vmrh2o     Total volume mixing ratio profile of water vapor.<br> 
                       This will be needed only for the oxygen continuum <br> 
                       however one is forced to give this input [1]
    \param  vmrn2      Total volume mixing ratio profile of molecular nitrogen.<br> 
                       This will be needed only for the CO2 foreign continuum [1]<br> 
                       however one is forced to give this input [1]
    \param  vmr        Volume mixing ratio profile of the actual species [1]

   \author Stefan Bühler, Thomas Kuhn
   \date 2001-11-05
 */

void xsec_continuum_tag( MatrixView                 xsec,
                         const String&              name,
                         ConstVectorView            parameters,
                         const String&              model,
                         ConstVectorView            f_mono,
                         ConstVectorView            p_abs,
                         ConstVectorView            t_abs,
                         ConstVectorView            n2_abs,
                         ConstVectorView            h2o_abs,
                         ConstVectorView              vmr )
{
  //
  // In the following all the possible tags are listed here and 
  // after a first consistency check about the input parameters the 
  // appropriate internal function is called.
  //
  // ============= H2O continuum ========================================================
  if ( "H2O-SelfContStandardType"==name )
    {
      // 
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum coefficient (C_s)  [1/m / (Hz²*Pa²)]
      //     parameters[1] : temperature exponent  (x_s)  [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 2;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Standard_H2O_self_continuum( xsec,
                                       parameters[0],
                                       parameters[1],
                                       model,
                                       f_mono,
                                       p_abs,
                                       t_abs,
                                       vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Standard_H2O_self_continuum( xsec,
                                       0.00,
                                       0.00,
                                       model,
                                       f_mono,
                                       p_abs,
                                       t_abs,
                                       vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContStandardType"==name )
    {
      //      
      // specific continuum parameters units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1/m / (Hz²*Pa²)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 2;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Standard_H2O_foreign_continuum( xsec,
                                          parameters[0],
                                          parameters[1],
                                          model,
                                          f_mono,
                                          p_abs,
                                          t_abs,
                                          vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Standard_H2O_foreign_continuum( xsec,
                                          0.00,
                                          0.00,
                                          model,
                                          f_mono,
                                          p_abs,
                                          t_abs,
                                          vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ContMPM93"==name )
    {
      // self and foreign continuum term are simultaneously calculated
      // since the parameterization can not be divided up in these two 
      // terms because they are not additive terms.
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : pseudo continuum line frequency                      [Hz]
      //     parameters[1] : pseudo continuum line strength parameter             [Hz/Pa]
      //     parameters[2] : pseudo continuum line strength temperature parameter [1]
      //     parameters[3] : pseudo continuum line broadening parameter           [Hz/Pa]
      //     parameters[4] : pseudo continuum line broadening parameter           [1]
      //     parameters[5] : pseudo continuum line broadening parameter           [1]
      //     parameters[6] : pseudo continuum line broadening parameter           [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 7;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93_H2O_continuum( xsec,
                               parameters[0],
                               parameters[1],
                               parameters[2],
                               parameters[3],
                               parameters[4],
                               parameters[5],
                               parameters[6],
                               model,
                               f_mono,
                               p_abs,
                               t_abs,
                               vmr       );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM93_H2O_continuum( xsec,
                               0.00,
                               0.00,
                               0.00,
                               0.00,
                               0.00,
                               0.00,
                               0.00,
                               model,
                               f_mono,
                               p_abs,
                               t_abs,
                               vmr       );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters. " << "\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-SelfContCKD"==name )
    {
      // WWW resource: ftp.aer.com/aer_contnm_ckd
          ostringstream os;
          os << "CKD self continuum model not yet implemented"
             << ".\n";
          throw runtime_error(os.str());
          return;

    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContCKD"==name ) 
    {
          ostringstream os;
          os << "CKD foreign continuum model not yet implemented"
             << ".\n";
          throw runtime_error(os.str());
          return;

    }
  // ============= H2O full models ======================================================
  else if ( "H2O-CP98"==name )
    {
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum scale factor       (CC) [1]
      //     parameters[1] : line strength scale factor   (CL) [1]
      //     parameters[2] : line broadening scale factor (CW) [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          CP98H2OAbsModel( xsec,
                           parameters[0],
                           parameters[1],
                           parameters[2],
                           model,
                           f_mono,
                           p_abs,
                           t_abs,
                           vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          CP98H2OAbsModel( xsec,
                           0.00,
                           0.00,
                           0.00,
                           model,
                           f_mono,
                           p_abs,
                           t_abs,
                           vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-MPM87"==name )
    {
      //
      // specific continuum parameters and units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : continuum scale factor (CC)       [1]
      //     parameters[1] : line strength scale factor   (CL) [1]
      //     parameters[2] : line broadening scale factor (CW) [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM87H2OAbsModel( xsec,
                            parameters[0],
                            parameters[1],
                            parameters[2],
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM87H2OAbsModel( xsec,
                            0.00,
                            0.00,
                            0.00,
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-MPM89"==name )
    {
      //
      // specific continuum parameters and units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : continuum scale factor       (CC) [1]
      //     parameters[1] : line strength scale factor   (CL) [1]
      //     parameters[2] : line broadening scale factor (CW  [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM89H2OAbsModel( xsec,
                        parameters[0],
                        parameters[1],
                        parameters[2],
                        model,
                        f_mono,
                        p_abs,
                        t_abs,
                        vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM89H2OAbsModel( xsec,
                            0.00,
                            0.00,
                            0.00,
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-MPM93"==name )
    {
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum scale factor       (CC) [1]
      //     parameters[1] : line strength scale factor   (CL) [1]
      //     parameters[2] : line broadening scale factor (CW) [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93H2OAbsModel( xsec,
                            parameters[0],
                            parameters[1],
                            parameters[2],
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM93H2OAbsModel( xsec,
                            0.00,
                            0.00,
                            0.00,
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-PWR98"==name )
    {      
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum scale factor       (CC) [1]
      //     parameters[1] : line strength scale factor   (CL) [1]
      //     parameters[2] : line broadening scale factor (CW) [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          PWR98H2OAbsModel( xsec,
                            parameters[0],
                            parameters[1],
                            parameters[2],
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          PWR98H2OAbsModel( xsec,
                            0.00,
                            0.00,
                            0.00,
                            model,
                            f_mono,
                            p_abs,
                            t_abs,
                            vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ============= O2 continuum =========================================================
  else if ( "O2-SelfContStandardType"==name )
    {
      // MPM93, Rosenkranz 1993 O2 continuum:
      // see publication side of National Telecommunications and Information Administration
      //   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
      // and ftp side for downloading the MPM93 original source code:
      //   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
      //
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3
      // (see also JQSRT, Vol.48, No.5/6 pp.629-643, 1992)
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum coefficient (C)    [1/m*1/Hz*1/Pa] 
      //     parameters[1] : frequency coefficient (G0)   [Hz/Pa]
      //     parameters[3] : line width parameter  (G0A)  [1]
      //     parameters[3] : line width parameter  (G0B)  [1]
      //     parameters[2] : temperature exponent  (XG0d) [1]
      //     parameters[2] : temperature exponent  (x_s)  [1]
      //     parameters[5] : continuum coefficient (XG0w) [1]
      //     model   : model option ("MPM93", "Rosenkranz", or "user")
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs       : [1]
      //     vmr           : [1]
      //
      const int Nparam = 6;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Standard_O2_continuum( xsec,
                                 parameters[0],
                                 parameters[1], 
                                 parameters[2], 
                                 parameters[3], 
                                 parameters[4], 
                                 parameters[5], 
                                 model, 
                                 f_mono,
                                 p_abs,
                                 t_abs,
                                 h2o_abs,
                                 vmr );
        }  
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Standard_O2_continuum( xsec,
                                 0.00,
                                 0.00,
                                 0.00,
                                 0.00,
                                 0.00,
                                 0.00,
                                 model, 
                                 f_mono,
                                 p_abs,
                                 t_abs,
                                 h2o_abs,
                                 vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-SelfContMPM93"==name )
    {
      // MPM93 O2 continuum:
      // see publication side of National Telecommunications and Information Administration
      //   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
      // and ftp side for downloading the MPM93 original source code:
      //   ftp://ftp.its.bldrdoc.gov/pub/mpm93/

      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum coefficient (C)    [1/m / (Hz²*Pa²)]
      //     parameters[1] : temperature exponent  (x_s)  [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs       : [1]
      //     vmr           : [1]
      //
      const int Nparam = 4;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93_O2_continuum( xsec,
                              parameters[0],
                              parameters[1], 
                              parameters[2], 
                              parameters[3], 
                              model, 
                              f_mono,
                              p_abs,
                              t_abs,
                              h2o_abs,
                              vmr );
        }  
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM93_O2_continuum( xsec,
                              0.00,
                              0.00,
                              0.00,
                              0.00,
                              model, 
                              f_mono,
                              p_abs,
                              t_abs,
                              h2o_abs,
                              vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-SelfContPWR93"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3
      // (see also JQSRT, Vol.48, No.5/6 pp.629-643, 1992)
      // 
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum coefficient (C) [K²/(Hz*Pa*m)]
      //     parameters[1] : temperature exponent  (x) [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 4;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Rosenkranz_O2_continuum( xsec,
                                   parameters[0],
                                   parameters[1],
                                   parameters[2],
                                   parameters[3],
                                   model,
                                   f_mono,
                                   p_abs,
                                   t_abs,
                                   h2o_abs,
                                   vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Rosenkranz_O2_continuum( xsec,
                                   0.00,
                                   0.00,
                                   0.00,
                                   0.00,
                                   model,
                                   f_mono,
                                   p_abs,
                                   t_abs,
                                   h2o_abs,
                                   vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ============= O2 full model ========================================================
  else if ( "O2-PWR93"==name )
    {
      //  REFERENCE FOR EQUATIONS AND COEFFICIENTS:
      //  P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
      //  BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
      //  AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
      //  (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum term scale factor,   default CC = 1.000 [1]
      //     parameters[1] : line strength scale factor,    default CL = 1.000 [1]
      //     parameters[1] : line broadening scale factor,  default CW = 1.000 [1]
      //     parameters[1] : line coupling scale factor,    default CO = 1.000 [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs,      : [1]
      //     vmr           : [1]
      //
      const int Nparam = 4;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          PWR93O2AbsModel( xsec,
                           parameters[0], // continuum term scale factor
                           parameters[1], // line strength scale factor
                           parameters[2], // line broadening scale factor
                           parameters[3], // line coupling scale factor
                           model,
                           f_mono,
                           p_abs,
                           t_abs,
                           h2o_abs,
                           vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          PWR93O2AbsModel( xsec,
                           0.00,
                           0.00,
                           0.00,
                           0.00,
                           model,
                           f_mono,
                           p_abs,
                           t_abs,
                           h2o_abs,
                           vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-MPM93"==name )
    {
      //  H. J. Liebe and G. A. Hufford and M. G. Cotton,
      //  "Propagation modeling of moist air and suspended water/ice
      //   particles at frequencies below 1000 GHz",
      //  AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      //  Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21       
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum term scale factor,   default CC = 1.000 [1]
      //     parameters[1] : line strength scale factor,    default CL = 1.000 [1]
      //     parameters[2] : line broadening scale factor,  default CW = 1.000 [1]
      //     parameters[3] : line coupling scale factor,    default CO = 1.000 [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs,      : [1]
      //     vmr           : [1]
      //
      const int Nparam = 4;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Full model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93O2AbsModel( xsec,
                           parameters[0], // continuum term scale factor
                           parameters[1], // line strength scale factor
                           parameters[2], // line broadening scale factor
                           parameters[3], // line coupling scale factor
                           model,
                           f_mono,
                           p_abs,
                           t_abs,
                           h2o_abs,
                           vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Full model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM93O2AbsModel( xsec,
                           0.00,
                           0.00,
                           0.00,
                           0.00,
                           model,
                           f_mono,
                           p_abs,
                           t_abs,
                           h2o_abs,
                           vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Full model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ============= N2 continuum =========================================================
  else if ( "N2-SelfContMPM93"==name )
    {
      // MPM93 N2 continuum:
      // see publication side of National Telecommunications and Information Administration
      //   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
      // and ftp side for downloading the MPM93 original source code:
      //   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : strength parameter   [1/m * 1/(Hz²*Pa²)]
      //     parameters[1] : broadening parameter [1]
      //     parameters[2] : temperature exponent [1]
      //     parameters[3] : frequency exponent   [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs       : [1]
      //     vmr           : [1]
      //
      const int Nparam = 4;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93_N2_continuum( xsec,
                              parameters[0],
                              parameters[1],
                              parameters[2],                      
                              parameters[3],                      
                              model,
                              f_mono,
                              p_abs,
                              t_abs,
                              h2o_abs,
                              vmr );
        }  
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          MPM93_N2_continuum( xsec,
                              0.00,
                              0.00,
                              0.00,
                              0.00,
                              model,
                              f_mono,
                              p_abs,
                              t_abs,
                              h2o_abs,
                              vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "N2-SelfContPWR93"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3
      //
      // specific continuum parameters and units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : continuum strength coefficient  [1/m * 1/(Hz*Pa)²]
      //     parameters[1] : continuum temperature exponent  [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 2;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Rosenkranz_N2_self_continuum( xsec,
                                        parameters[0], // coefficient
                                        parameters[1], // temp. exponent
                                        model,
                                        f_mono,
                                        p_abs,
                                        t_abs,
                                        vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Rosenkranz_N2_self_continuum( xsec,
                                        0.00,
                                        0.00,
                                        model,
                                        f_mono,
                                        p_abs,
                                        t_abs,
                                        vmr );    
        } 
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "N2-SelfContStandardType"==name )
    {
      // data information about this continuum: 
      // A completely general expression for the N2 continuum
      //
      // specific continuum parameters and units:
      // OUTPUT 
      //     xsec          : [1/m],
      // INPUT
      //     parameters[0] : continuum coefficient (C)  [1/m * 1/(Hz*Pa)²]
      //     parameters[1] : frequency exponent    (xf) [1]
      //     parameters[2] : temperature exponent  (xt) [1]
      //     parameters[3] : pressure exponent     (xp) [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      const int Nparam = 4;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Standard_N2_self_continuum( xsec,
                                      parameters[0],
                                      parameters[1],
                                      parameters[2],
                                      parameters[3],
                                      model,
                                      f_mono,
                                      p_abs,
                                      t_abs,
                                      vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Standard_N2_self_continuum( xsec,
                                      0.000,
                                      0.000,
                                      0.000,
                                      0.000,
                                      model,
                                      f_mono,
                                      p_abs,
                                      t_abs,
                                      vmr );
        } 
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "N2-SelfContBorysow"==name )
    {
      // data information about this continuum: 
      // A. Borysow and L. Frommhold, The Astrophysical Journal,
      // Vol. 311, pp.1043-1057, 1986
      ostringstream os;
      os << "N2 continuum parameterization of A. Borysow and L. Frommhold\n"
         << "is not yet implemented  ==>  no calculation performed!\n";
      throw runtime_error(os.str());
      return;
      /*
      Borysow_Frommhold_N2_continuum( xsec,
                                      parameters[0],
                                      parameters[1],
                                      model,
                                      f_mono,
                                      p_abs,
                                      t_abs,
                                      vmr );
      */
    }  
  // ============= CO2 continuum ========================================================
  else if ( "CO2-SelfContPWR93"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT 
      //     parameters[0] : continuum strength coefficient [1/m * 1/(Hz*Pa)²]
      //     parameters[1] : continuum temperature exponent  [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 2;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Rosenkranz_CO2_self_continuum( xsec,
                                         parameters[0], // coefficient
                                         parameters[1], // temp. exponent
                                         model,
                                         f_mono,
                                         p_abs,
                                         t_abs,
                                         vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Rosenkranz_CO2_self_continuum( xsec,
                                         0.00,
                                         0.00,
                                         model,
                                         f_mono,
                                         p_abs,
                                         t_abs,
                                         vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters. " << "\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "CO2-ForeignContPWR93"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum strength coefficient [1/m * 1/(Hz*Pa)²]
      //     parameters[1] : continuum temperature exponent  [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     n2_abs        : [1]
      //     vmr           : [1]
      //
      const int Nparam = 2;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "Continuum model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          Rosenkranz_CO2_foreign_continuum( xsec,
                                            parameters[0],
                                            parameters[1],
                                            model,
                                            f_mono,
                                            p_abs,
                                            t_abs,
                                            n2_abs,
                                            vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "Continuum model " << name << " requires " << Nparam << " input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
          Rosenkranz_CO2_foreign_continuum( xsec,
                                            0.00,
                                            0.00,
                                            model,
                                            f_mono,
                                            p_abs,
                                            t_abs,
                                            n2_abs,
                                            vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: Continuum model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters. " << "\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ============= cloud and fog absorption from MPM93 ==================================
  else if ( "liquidcloud-MPM93"==name )
    {
      // Suspended water droplet absorption parameterization from MPM93 model
      // H. J. Liebe and G. A. Hufford and M. G. Cotton,
      // "Propagation modeling of moist air and suspended water/ice
      //  particles at frequencies below 1000 GHz",
      // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT 
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     model        : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      // liquid water droplet parameters:
      // suspended water droplet density   range: 0-10 g/m³
      // specific droplet weight           value:    1 g/cm³
      //
      // valid atmospheric condition:
      // temperature      : 233 to 323 K
      // relative humidity:   1 to 100 %
      // 
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "MPM93 liquid water cloud absorption model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93WaterDropletAbs(xsec,
                               parameters[0],     // scaling factror
                               parameters[1],     // scaling factror
                               parameters[2],     // scaling factror
                               model,       // model option
                               f_mono,
                               p_abs,
                               t_abs,
                               vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "MPM93 liquid water cloud absorption model " << name << " requires\n" 
             << Nparam << " input parameter for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "MPM93 liquid water cloud absorption model " << name << " running with \n" 
               << "the parameter for model " << model << ".\n";
          MPM93WaterDropletAbs(xsec,
                               0.000,        // scaling factror
                               0.000,        // scaling factror
                               0.000,        // scaling factror
                               model,  // model option
                               f_mono,
                               p_abs,
                               t_abs,
                               vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: MPM93 liquid water cloud absorption model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  // ============= ice particle absorption from MPM93 ===================================
  else if ( "icecloud-MPM93"==name )
    {
      // Ice particle absorption parameterization from MPM93 model
      // H. J. Liebe and G. A. Hufford and M. G. Cotton,
      // "Propagation modeling of moist air and suspended water/ice
      //  particles at frequencies below 1000 GHz",
      // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT 
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     model        : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      // ice crystal parameters:
      // suspended water droplet density   range: 0-10 g/m³
      // specific droplet weight           value:    1 g/cm³
      //
      // valid atmospheric condition:
      // temperature      : 233 to 323 K
      // relative humidity:   1 to 100 %
      // 
      const int Nparam = 3;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
        {
          out2 << "MPM93 ice water cloud absorption model " << name << " is running with \n"
               << "user defined parameters according to model " << model << ".\n";
          MPM93IceCrystalAbs(xsec,
                             parameters[0],     // scaling factror
                             parameters[1],     // scaling factror
                             parameters[2],     // scaling factror
                             model,       // model option
                             f_mono,
                             p_abs,
                             t_abs,
                             vmr );
        }
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
        {
          ostringstream os;
          os << "MPM93 ice water cloud absorption model " << name << " requires \n"
             << Nparam << " input parameter for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
          throw runtime_error(os.str());
          return;
        }
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
        {
          out2 << "MPM93 ice water cloud absorption model " << name << " running with \n" 
               << "the parameter for model " << model << ".\n";
          MPM93IceCrystalAbs(xsec,
                             0.000,       // scaling factror
                             0.000,       // scaling factror
                             0.000,       // scaling factror
                             model, // model option
                             f_mono,
                             p_abs,
                             t_abs,
                             vmr );
        }
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
        {
          ostringstream os;
          os << "ERROR: MPM93 ice water cloud absorption model " << name << " requires NO input\n"
             << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n"
             << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 2.\n";
          throw runtime_error(os.str());
          return;
        }
    }
  else // -----------------------------------------------------------------------
    {
      // none of the continuum or full model tags were selected -> error message.
      ostringstream os;
      os << "ERROR: Continuum/ full model tag `" << name << "' not yet implemented in arts!";
      throw runtime_error(os.str());
      return;
    }
 return;
}
//
// #################################################################################
//
/**
   An auxiliary functions that checks if a given continuum model is
   listed in species_data.cc. This is just in order to verify that this
   really represent a valid continuum model.

   The given name should be something like
   `ContStandardSelf'. The function simply checks if there is a
   species H2O with an isotope ContStandardSelf.

   For user-friendliness, the function also compiles a list of allowed
   continuum models and gives this as an error message if the model is
   not found. 

   The function has no return value, since, if the name does not match
   a valid model an error is thrown anyway.

   \param name The name of the continuum model to check.

   \throw runtime_error The model does not exist.

   \author Stefan Buehler
   \date   2001-03-12
*/
void check_continuum_model(const String& name)
{
  // The species lookup data:
  extern const Array<SpeciesRecord> species_data;

  // For the list of valid continuum models:
  ArrayOfString valid_models;

  bool found = false;

  // Loop all species:
  for ( Array<SpeciesRecord>::const_iterator i=species_data.begin();
        i<species_data.end();
        ++i )
    {
      String specnam = i->Name();

      // Loop all isotopes:
      for ( Array<IsotopeRecord>::const_iterator j=i->Isotope().begin();
            j<i->Isotope().end();
            ++j )
        {
          String isonam = j->Name();

          // The specified name consists of a species part and an
          // isotope part, e.g., H2O-ContStandardSelf. We need to
          // construct a similar String from the species lookup data
          // by concatenating species name and isotope name.

          String fullnam = specnam + "-" + isonam;
          //      cout << fullnam << "\n";

          // See if this is a continuum tag, so that we can add it to
          // the list:
          if ( 0 > j->Abundance() )
            {
              valid_models.push_back(fullnam);        
            }
          
          if ( name == fullnam )
            {
              found = true;
            }
        }
    }

  // ----------------------------------------------------------------------
  // Have we found it?
  if (!found)
    {
      ostringstream os;
      os << "The String `" << name << "' matches none of the known\n"
         << "continuum models. Known continuum models are:";
      for ( ArrayOfString::const_iterator i=valid_models.begin();
            i<valid_models.end();
            ++i )
        {
          os << "\n" << *i;
        }      
      throw runtime_error(os.str());
    }
}
//
// #################################################################################
