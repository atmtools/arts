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


   \retval xsec  Absorption cross section, defined such that the
                 absorption coefficient \f$\alpha\f$ (in units of 1/m) is:<br>
                 \f$\alpha\f$ = xsec * VMR.<br>
		 The absorption model functions adds absorption to xsec, 
                 rather than replacing the previous content. <br>
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
   <li><b>H2O-H2O (H2O-SelfContCKD222)</b>:<br> 
         CKDv2.2.2 H2O self continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>H2O-H2O (H2O-SelfContCKD242)</b>:<br> 
         CKDv2.4.2 H2O self continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>H2O-H2O (H2O-SelfContCKDMT100)</b>:<br> 
         CKD_MTv1.00 H2O self continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>H2O-air (H2O-ForeignContStandardType)</b>: <br>
         P. W. Rosenkranz,<br> 
         Radio Science, Vol. 33, No 4, Pages 919-928, 1998 and <br>
         Radio Science, Vol. 34, No 4, Page 1025, 1999
         <a href="ftp://mesa.mit.edu/phil/lbl_rt">(WWW access)</a>.
   </li>
   <li><b>H2O-air (H2O-foreignContCKD222)</b>:<br> 
         CKDv2.2.2 H2O foreign continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>H2O-air (H2O-foreignContCKD242)</b>:<br> 
         CKDv2.4.2 H2O foreign continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>H2O-air (H2O-foreignContCKDMT100)</b>:<br> 
         CKD_MTv1.00 H2O foreign continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>H2O-air (H2O-ForeignContMaTippingType)</b>: <br>
         Q. Ma and R. H. Tipping,<br> 
         J. Chem. Phys., 117(23), 10581, 2002.<br>
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
   <li><b>O2-air (O2-CIAfunCKDMT100)</b>:<br> 
         CKD_MT version 1.00 O2 CIA fundamental band from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>O2-air (O2-v0v0CKDMT100)</b>:<br> 
         CKD_MT version 1.00 O2 band absorption model for the \f$a^1\Delta_g\f$ 
	 \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
	 \f$X^3\Sigma^-_g\f$ band system 
	 (\f$\nu=0\f$
	 \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
	 \f$\nu=0\f$
	 transitions around 1.27 microns).<br>  
         Source code from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>O2-air (O2-v1v0CKDMT100)</b>:<br> 
                  CKD_MT version 1.00 O2 band absorption model for the \f$a^1\Delta_g\f$ 
	 \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
	 \f$X^3\Sigma^-_g\f$ band system 
	 (\f$\nu=1\f$
	 \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
	 \f$\nu=0\f$
	 transitions around 1.06 microns).<br>
         Source code from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
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
   <li><b>N2-N2 (N2-CIArotCKDMT100)</b>:<br> 
         CKD_MT version 1.00 N2-N2 CIA rotational band absorption model.<br>
         Source code from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   <li><b>N2-N2 (N2-CIAfunCKDMT100)</b>:<br> 
         CKD_MT version 1.00 N2-N2 CIA fundamental band absorption model.<br>
         Source code from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
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
   <li><b>CO2-air (CO2-CKD241)</b>:<br> 
         CKDv2.4.1 CO2 continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
   </li>
   <li><b>CO2-air (CO2-CKDMT100)</b>:<br> 
         CKD_MT version 1.00 CO2 continuum from the FORTRAN77 code written by<br>  
         <a href="http://www.rtweb.aer.com/continuum_frame.html">Atmospheric and 
         Environmental Research Inc. (AER),</a><br>
	 Radiation and Climate Group<br>
	 131 Hartwell Avenue<br>
	 Lexington, MA 02421, USA
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
   \date   2003-11-19
*/

#include <math.h>
#include "arts.h"
#include "matpackI.h"
#include "array.h"
#include "absorption.h"
#include "messages.h"
#include "continua.h"

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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
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
  out3  << "H2O-MPM87: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n";
  
  
  // number of lines of liebe line catalog (30 lines)
  const Index i_first = 0;
  const Index i_last  = 29;
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric Nppc  = CC * pwv_dummy * pow(theta, (Numeric)3.0) * 1.000e-5
        * ( (0.113 * pda) + (3.57 * pwv * pow(theta, (Numeric)7.8)) );

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
	      Numeric strength = CL * pwv_dummy * mpm87[l][1]
                * pow(theta,(Numeric)3.5) * exp(mpm87[l][2]*(1.000-theta));
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * mpm87[l][3] * 
                ( (4.80 * pwv * pow(theta, (Numeric)1.1)) + 
                  (       pda * pow(theta, (Numeric)0.6)) );
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
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
  out3  << "H2O-MPM89: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n";
  
  
  // number of lines of Liebe line catalog (30 lines)
  const Index i_first = 0;
  const Index i_last  = 29;
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric Nppc      = CC * pwv_dummy * pow(theta, (Numeric)3.0) * 1.000e-5
        * ( (0.113 * pda) + (3.57 * pwv * pow(theta, (Numeric)7.5)) );
      
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
	      Numeric strength = CL * pwv_dummy * mpm89[l][1]
                * pow(theta, (Numeric)3.5) * exp(mpm89[l][2]*(1.000-theta));
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * mpm89[l][3] * 0.001
                * ( mpm89[l][5] * pwv * pow(theta, mpm89[l][6]) +  
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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

   \attention Corrected version of MPM93 by TKS, iup, 2002
              The H2O lines at 547.676440 GHz and 552.020960 GHz are isotopic lines:<br>
              547 GHz is from the isotope 1-8-1 (HITRAN code 181, JPL code 20003) with an 
              isotopic ratio of 0.00199983 and <br>
              552 GHz is from the isotope 1-7-1  (HITRAN code 171, JPL code 19003) with an 
              isotopic ratio of 0.00037200.<br>
              The original source code of MPM93 has these isotopic ratios not included 
              in the line strength parameter b1, which is an error.<br> 
              In the arts implementation the line strength parameter b1 of these two lines 
              is multiplied with the appropriate isotopic ratio.

   \author Thomas Kuhn
   \date 2002-05-06
 */ 

void MPM02H2OAbsModel( MatrixView        xsec,
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,
		       ConstVectorView   f_mono,
		       ConstVectorView   p_abs,
		       ConstVectorView   t_abs,
		       ConstVectorView   vmr )
{
  //
  /*
CTKS  OTHER DATA USED IF NOT FROM THEORETICAL CALC. IN A. BAUER ET AL. 41(1989)49-54:
CTKS  --------------------------------------------------------------------------------------------------------------
CTKS           | T=300 K       | T=300 K      |   T=300 K        |
CTKS  F     ISO|GWVHZO   NWVHZO| GWVNZ    NWVNZ|   GWVAIR   NWVAIR| REFERENCE
CTKS  GHZ   1  |MHZ/TORR 1     | MHZ/TORR 1    |   MHZ/TORR 1     |
CTKS  --------------------------------------------------------------------------------------------------------------
CTKS   22.2 1   18.00(18) -      4.10     --       3.77     --       LIEBE ET AL., J.CHEM.PHYS., 50(1969)727
CTKS  183.3 1   19.88    0.85    4.07(7)  0.63(10) 3.75(6)  0.64(10) A. BAUER ET AL. JQSRT 41(1989)49-54
CTKS  183.3 1    -        -      4.19(17) 0.74(3)  3.89(14) 0.76(3)  T. M. GOYETTE ET AL. J. MOLEC. SPEC, 143(1990)346
CTKS  203.4 2   --       --      4.214    0.93     3.833    0.89     J.-M. COLMONT ET AL. J. MOLEC. SPEC. 193(1999)233-243
CTKS  225.9 4   --       --      4.21     0.70     3.798    0.75     J.-M. COLMONT ET AL. J. MOLEC. SPEC. 193(1999)233-243
CTKS  241.6 4   --       --      4.45     0.77     4.08     0.80     J.-M. COLMONT ET AL. J. MOLEC. SPEC. 193(1999)233-243
CTKS  241.9 4   --       --      3.47     0.67     3.07     0.70     J.-M. COLMONT ET AL. J. MOLEC. SPEC. 193(1999)233-243
CTKS  325.1 1   --       --      4.011    0.63     3.633    0.64     J.-M. COLMONT ET AL. J. MOLEC. SPEC. 193(1999)233-243
CTKS  380.2 1   20.61(7) 0.89(1) 4.24(7)  0.52(14) 3.83(6)  0.54(14) A. BAUER ET AL. JQSRT 41(1987) 531
CTKS  380.2 1    -        -      4.16(4)  0.70(3)  3.80     0.72     T. M. GOYETTE ET AL. JQSRT 41(1993)485
CTKS  439.2 1   12.95(25)0.62(9)  --      --       --       --       V. N. MARKOV, J. MOLEC. SPEC, 164(1994)233
CTKS  752.0 1                    4.16(18) --       3.75     --       S. S. D. GASSTER ET AL. JOSA, 5(1988)593       
CTKS  987.9 1                    4.42(23) --       4.01     --       S. S. D. GASSTER ET AL. JOSA, 5(1988)593  
*/
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0             1        2       3        4       5      6
  //         f0            b1       b2      b3       b4      b5     b6
  //        [MHz]       [kHz/kPa]   [1]       [MHz/hPa]     [1]    [1]
  //                                        air      self   air    self
  const Numeric mpm02[35][7] = { 
    {    22235.0800,    0.10947,  2.1678,   2.811,   4.80,  0.69,  0.61},
    {    67803.9600,    0.00111,  8.7518,   2.858,   4.93,  0.69,  0.82},
    {   119995.9400,    0.00072,  8.3688,   2.948,   4.78,  0.70,  0.79},
    {   183310.1170,    2.30351,  0.6794,   3.050,   5.30,  0.76,  0.85},
    {   321225.6400,    0.04646,  6.1792,   2.303,   4.69,  0.67,  0.54},
    {   325152.9190,    1.53869,  1.5408,   2.783,   4.85,  0.68,  0.74},
    {   336227.6200,    0.00099,  9.8233,   2.693,   4.74,  0.64,  0.61},
    {   380197.3720,    11.9079,  1.0439,   2.873,   5.38,  0.72,  0.89},
    {   390134.5080,    0.00437,  7.3408,   2.152,   4.81,  0.63,  0.55},
    {   437346.6670,    0.06378,  5.0384,   1.845,   4.23,  0.60,  0.48},
    {   439150.8120,    0.92144,  3.5853,   2.100,   4.29,  0.63,  0.62},
    {   443018.2950,    0.19384,  5.0384,   1.860,   4.23,  0.60,  0.50},
    {   448001.0750,    10.6190,  1.3952,   2.632,   4.84,  0.66,  0.67},
    {   470888.9470,    0.33005,  3.5853,   2.152,   4.57,  0.66,  0.65},
    {   474689.1270,    1.27660,  2.3674,   2.355,   4.65,  0.65,  0.64},
    {   488491.1330,    0.25312,  2.8391,   2.602,   5.04,  0.69,  0.72},
    {   503568.5320,    0.03746,  6.7158,   1.612,   3.98,  0.61,  0.43},
    {   504482.6920,    0.01250,  6.7158,   1.612,   4.01,  0.61,  0.45},
    {   547676.4400,    1.01467,  0.1427,   2.600,   4.50,  0.69,  1.00}, // *
    {   552020.9600,    0.18668,  0.1452,   2.600,   4.50,  0.69,  1.00}, // *
    {   556936.0020,  510.51086,  0.1405,   3.210,   4.11,  0.69,  1.00},
    {   620700.8070,    5.10539,  2.3673,   2.438,   4.68,  0.71,  0.68},
    {   645905.6200,    0.00667,  8.6065,   1.800,   4.00,  0.60,  0.43},
    {   658006.5500,    0.27451,  7.7889,   3.210,   4.14,  0.69,  1.00},
    {   752033.2270,  249.68466,  0.3625,   3.060,   4.09,  0.68,  0.84},
    {   841051.1620,    0.01308,  8.1347,   1.590,   5.76,  0.33,  0.45},
    {   859965.6490,    0.13326,  8.0114,   3.060,   4.09,  0.68,  0.84},
    {   899302.1710,    0.05492,  7.8676,   2.985,   4.53,  0.68,  0.90},
    {   902609.4360,    0.03854,  8.3823,   2.865,   5.10,  0.70,  0.95},
    {   906206.1180,    0.18323,  5.0628,   2.408,   4.70,  0.70,  0.53},
    {   916171.5820,    8.56444,  1.3943,   2.670,   4.78,  0.70,  0.78},
    {   923113.1900,    0.00784, 10.2441,   2.900,   5.00,  0.66,  0.67},
    {   970315.0220,    9.16280,  1.8673,   2.550,   4.94,  0.64,  0.67},
    {   987926.7640,  138.28461,  0.2045,   2.985,   4.55,  0.68,  0.90},
  //--------------------------------------------------------------------
    {  1780.000000, 2230.00000,  0.952,  17.620,  30.50,  2.00,  5.00}}; // pseudo continuum line

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


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW;
  // number of lines of Liebe line catalog (0-33 lines, 34 cont. pseudo line)
  Index i_first = 0;
  Index i_last  = 34;
  if ( model == "MPM02" )
    {
      CC      = CC_MPM02;
      CL      = CL_MPM02;
      CW      = CW_MPM02;
      i_first = 0;
      i_last  = 34;
    }
  else if ( model == "MPM02Lines" )
    {
      CC      = 0.000;
      CL      = CL_MPM02;
      CW      = CW_MPM02;
      i_first = 0;
      i_last  = 33;
    }
  else if ( model == "MPM02Continuum" )
    {
      CC      = CC_MPM02;
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
      os << "H2O-MPM02: ERROR! Wrong model values given.\n"
	 << "Valid models are: 'MPM02', 'MPM02Lines', 'MPM02Continuum', and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "H2O-MPM02: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n";
  
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      // Loop over MPM02 spectral lines:
      
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
		  strength = CL * pwv_dummy * mpm02[l][1] * 
			          pow(theta, (Numeric)3.5)  * exp(mpm02[l][2]*(1.0-theta));
		  // line broadening parameter [GHz]
		  gam      = CW * mpm02[l][3] * 0.001 * 
		                  ( (mpm02[l][4] * pwv * pow(theta, mpm02[l][6])) +  
		                  (                pda * pow(theta, mpm02[l][5])) );
		}
	      else if ( l == 34 ) // ----- just the continuum pseudo-line ----------
		{
		  strength = CC * pwv_dummy * mpm02[l][1] * 
		                  pow(theta, (Numeric)3.5)  * exp(mpm02[l][2]*(1.0-theta));
		  // line broadening parameter [GHz]
		  gam      = mpm02[l][3] * 0.001 * 
	                     ( (mpm02[l][4] * pwv * pow(theta, mpm02[l][6])) +  
	                     (                pda * pow(theta, mpm02[l][5])) );
		}
	      else // ----- if something strange happens ---------------------------
		{
		  ostringstream os;
		  os << "H2O-MPM02: wrong line number detected l=" << l << " (0-34)\n";
		  throw runtime_error(os.str());
		  return;
		} // ---------------------------------------------------------------
	      // Doppler line width [GHz]
	      // Numeric gamd     = 1.46e-6 * mpm02[l][0] / sqrt(theta);
	      // effective line width [GHz]
	      //gam              = 0.535 * gam + sqrt(0.217*gam*gam + gamd*gamd); 
	      // absorption [dB/km] like in the original MPM02
	      Numeric Npp = strength * MPMLineShapeFunction(gam, mpm02[l][0], ff); 
	      // xsec = abs/vmr [1/m] but MPM89 is in [dB/km] --> conversion necessary
	      xsec(s,i)   += dB_km_to_1_m * 0.1820 * ff * Npp;
	    }
	}
    }
  return;
}
//
//
// #################################################################################
//
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
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
  //--------------------------------------------------------------------
    {  1780.000000, 2230.00000,  0.952,  17.620,  30.50,  2.00,  5.00}}; // pseudo continuum line

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
  out3  << "H2O-MPM93: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n";
  
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
		  strength = CL * pwv_dummy * mpm93[l][1]
                    * pow(theta, (Numeric)3.5)  * exp(mpm93[l][2]*(1.0-theta));
		  // line broadening parameter [GHz]
		  gam      = CW * mpm93[l][3] * 0.001 * 
		                  ( (mpm93[l][4] * pwv * pow(theta, mpm93[l][6])) +  
		                  (                pda * pow(theta, mpm93[l][5])) );
		}
	      else if ( l == 34 ) // ----- just the continuum pseudo-line ----------
		{
		  strength = CC * pwv_dummy * mpm93[l][1]
                    * pow(theta, (Numeric)3.5)  * exp(mpm93[l][2]*(1.0-theta));
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
		       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
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
  out3  << "H2O-PWR98: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n";
  
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric den        = 3.335e16 * (2.1667 * p_abs[i] * vmr[i] / t_abs[i]);
      Numeric den_dummy  = 3.335e16 * (2.1667 * p_abs[i] / t_abs[i]);
      // inverse relative temperature [1]
      Numeric ti         = (300.0 / t_abs[i]);
      Numeric ti2        = pow(ti, (Numeric)2.5);
      
      // continuum term [Np/km/GHz2]
      Numeric con = CC * pvap_dummy * pow(ti, (Numeric)3.0) * 1.000e-9
        * ( (0.543 * pda) + (17.96 * pvap * pow(ti, (Numeric)4.5)) );
	  
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
	      Numeric width    = ( CW * PWRw3[l] * pda  * pow(ti, PWRx[l]) ) + 
				 (      PWRws[l] * pvap * pow(ti, PWRxs[l]));
	      //	      Numeric width    = CW * ( PWRw3[l] * pda  * pow(ti, PWRx[l]) + 
	      //					PWRws[l] * pvap * pow(ti, PWRxs[l]) );
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
	      sum             += strength * res * pow( (ff/PWRfl[l]),
                                                       (Numeric)2.0 );
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
                      const Numeric	CCin,       // continuum scale factor 
		      const Numeric     CLin,       // line strength scale factor
		      const Numeric	CWin,       // line broadening scale factor
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
  out3  << "H2O-CP98: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n";
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
	  Numeric TL    = CL * 0.0109 * pwv * pow(theta,(Numeric)3.5)
            * exp(2.143*(1.0-theta));
	  // line broadening parameter [GHz]
	  Numeric gam   = CW * 0.002784 *  
	                  ( (pda * pow(theta,(Numeric)0.6))
                            + (4.80 * pwv * pow(theta,(Numeric)1.1)) );
	  // continuum term
	  Numeric TC    = CC * pwv * pow(theta, (Numeric)3.0) * 1.000e-7
            * ( (0.113 * pda) + (3.57 * pwv * pow(theta,(Numeric)7.5)) );

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
   \retval   xsec           cross section (absorption/volume mixing ratio) of the 
                            H2O-H2O continuum [1/m]
   \param    Cin            constant absorption strength     [1/m / (Hz*Pa)²]
   \param    xin            temperature exponent of (300/T)  [1]
   \param    model          allows user defined input parameter set 
                            (C and x)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid        [Hz]
   \param    p_abs          predefined pressure grid         [Pa]
   \param    t_abs          predefined temperature grid      [K] 
   \param    vmr            H2O volume mixing ratio          [1]

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
				  ConstVectorView   vmr	 )
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
   out3  << "H2O-SelfContStandardType: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";



  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
	C * pow( (Numeric)300./t_abs[i], x+(Numeric)3. )
        * pow( p_abs[i], (Numeric)2. ) * vmr[i];

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
	{
	  xsec(s,i) += dummy * pow( f_mono[s], (Numeric)2. );
	  //	  cout << "xsec(" << s << "," << i << "): " << xsec(s,i) << "\n";
	}
    }
}
//
// #################################################################################
//
/**
   \retval   xsec           cross section (absorption/volume mixing ratio) of the 
                            H2O-dry air continuum [1/m]
   \param    Cin            constant absorption strength [1/m / (Hz*Pa)²]
   \param    xin            temperature exponent         [1] 
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
				     const Numeric     Cin,
				     const Numeric     xin,
				     const String&     model,
				     ConstVectorView   f_mono,
				     ConstVectorView   p_abs,
				     ConstVectorView   t_abs,
				     ConstVectorView   vmr	 )
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
   out3  << "H2O-ForeignContStandardType: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric dummy = C * pow( (Numeric)300./t_abs[i], x+(Numeric)3. )
        * p_abs[i] * pdry;

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  xsec(s,i) += dummy * pow( f_mono[s], (Numeric)2. );
	  //	  cout << "xsec(" << s << "," << i << "): " << xsec(s,i) << "\n";
	}
    }
}
//
//
// #################################################################################
//
/**
   \retval   xsec           cross section (absorption/volume mixing ratio) of the 
                            H2O-dry air continuum [1/m]
   \param    Cin            constant absorption strength [1/m / (Hz*Pa)²]
   \param    xin            temperature exponent         [1] 
   \param    model          allows user defined input parameter set 
                            (C and x)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid    [Hz]
   \param    p_abs          predefined pressure          [Pa]
   \param    t_abs          predefined temperature grid  [K] 
   \param    vmr            H2O volume mixing ratio     [1] 

   \note     Except for  model 'user' the input parameters C and x 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MaTipping', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: Q. Ma and R. H. Tipping, J. Chem. Phys., 117(23), 10581, 2002.

   \author Thomas Kuhn
   \date 2002-12-04
 */ 
void MaTipping_H2O_foreign_continuum( MatrixView        xsec,
				      const Numeric	 Cin,
				      const Numeric	 xin,
				      const String&     model,
				      ConstVectorView   f_mono,
				      ConstVectorView   p_abs,
				      ConstVectorView   t_abs,
				      ConstVectorView   vmr	 )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for Q. Ma and R. H. Tipping, J. Chem. Phys., 117(23), 10581, 2002:
  // the Cf value is originally given in dB/km/kPa^2/GHz^2.0389. the conversion factor is 
  // then 1.0283E-28 to get arts units. Additionally the Cf value is divided by 1.08 to 
  // get the Cf for air.
  const Numeric Cf_MaTipping = 1.8590e-35;   //  [1/m / (Hz²*Pa²)]
  const Numeric xf_MaTipping = 4.6019;       //  [1]
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model goes for values!!):
  Numeric C, x;
   if ( model == "MaTipping" )
     {
       C = Cf_MaTipping;
       x = xf_MaTipping;
     }
   else if ( model == "user" )
     {
       C = Cin;
       x = xin;
     }
   else
     {
       ostringstream os;
       os << "H2O-MaTipping_H2O_foreign_continuum: ERROR! Wrong model values given.\n"
	  << "allowed models are: 'MaTipping', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out3  << "H2O-MaTipping_H2O_foreign_continuum: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric dummy = C * pow( (Numeric)300./t_abs[i], x )
        * p_abs[i] * pdry;

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  xsec(s,i) += dummy * pow( f_mono[s], (Numeric)2.0389 );
	  //	  cout << "xsec(" << s << "," << i << "): " << xsec(s,i) << "\n";
	}
    }
}
//
// #################################################################################



// =================================================================================


Numeric XINT_FUN( const Numeric V1A,
                  const Numeric V2A,
                  const Numeric DVA,
		  const Numeric A[],
                  const Numeric VI)
{

// ----------------------------------------------------------------------
  //     THIS SUBROUTINE INTERPOLATES THE A ARRAY STORED               
  //     FROM V1A TO V2A IN INCREMENTS OF DVA INTO XINT
// ----------------------------------------------------------------------

  const Numeric ONEPL  = 1.001;     // original value given in F77 code
  const Numeric ONEMI  = 0.999;     // original value given in F77 code

  //const Numeric ONEPL  = 0.001;  // modified value for C/C++ code

  Numeric RECDVA = 1.00e0/DVA;
    
  int J      = (int) ((VI-V1A)*RECDVA + ONEPL) ; 
  Numeric VJ = V1A + DVA * (Numeric)(J-1);    
  Numeric P  = RECDVA * (VI-VJ);        
  Numeric C  = (3.00e0-2.00e0*P) * P * P;         
  Numeric B  = 0.500e0 * P * (1.00e0-P);          
  Numeric B1 = B * (1.00e0-P);              
  Numeric B2 = B * P;                   
  
  Numeric xint = -A[J-1] * B1             +
                  A[J]   * (1.00e0-C+B2)  + 
                  A[J+1] * (C+B1)         - 
                  A[J+2] * B2;

  /*
  cout << (J-1) << " <-> " << (J+2)  
       << ",  V=" << VI << ", VJ=" << VJ << "\n";
  cout << "xint=" << xint  << " " << A[J-1] << " " << A[J] << " " << A[J+1] << " " << A[J+2] << "\n";
  */

  return xint;
}

// =================================================================================

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
									   
  if (XKT > 0.0)
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

// CKD version 2.2.2 H2O self continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O self continuum according to CKD_2_2_2  [1/m]
   \param    Cin            strength scaling factor                    [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD version 2.2.2 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-31-10
*/ 
void CKD_222_self_h2o( MatrixView          xsec,
		       const Numeric       Cin,
		       const String&       model,
		       ConstVectorView     f_mono,
		       ConstVectorView     p_abs,
		       ConstVectorView     t_abs,
		       ConstVectorView     vmr,
		       ConstVectorView     n2_abs )
{


  // check the model name about consistency
  if ((model != "user") &&  (model != "CKD222"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKDv2.2.2 H2O self continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKD222\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the self H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt   = 2.686763e19; // [molecules/cm^3]
  const Numeric T1       =  273.0e0;
  const Numeric TO       =  296.0e0;
  const Numeric PO       = 1013.0e0;

  // CKD2.2.2 specific self continuum correction function parameters
  const Numeric ALPHA2 = 200.000 * 200.000;
  const Numeric ALPHS2 = 120.000 * 120.000;
  const Numeric BETAS  = 5.000e-06;
  const Numeric V0S    = 1310.000;
  const Numeric FACTRS = 0.150;
  
  // These are self-continuum modification factors from 700-1200 cm-1
  const Numeric XFAC[51] = {  
    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
    1.00000,1.00000,1.00000};

  // wavenumber range where CKD H2O self continuum is valid
  const Numeric VABS_min = SL260_ckd_0_v1; // [cm^-1]
  const Numeric VABS_max = SL260_ckd_0_v2; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD2.2.2 H2O self continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << SL296_ckd_0_v1 << "<->" << SL296_ckd_0_v2 << "cm^-1\n";
    }


  // ------------------- subroutine SL296/SL260 ----------------------------

  if (SL296_ckd_0_v1 != SL260_ckd_0_v1)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.2.2 H2O self continuum:\n"
	 << "parameter V1 not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (SL296_ckd_0_v2 != SL260_ckd_0_v2)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.2.2 H2O self continuum:\n"
	 << "parameter V2 not the same for different ref. temperatures.\n";
	throw runtime_error(os.str());
    }
  if (SL296_ckd_0_dv != SL260_ckd_0_dv)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.2.2 H2O self continuum:\n"
	 << "parameter DV not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (SL296_ckd_0_npt != SL260_ckd_0_npt)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.2.2 H2O self continuum:\n"
	 << "parameter NPT not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  
  // retrieve the appropriate array sequence of the self continuum
  // arrays of the CKD model.
  Numeric DVC = SL296_ckd_0_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-SL296_ckd_0_v1) / SL296_ckd_0_dv);
  if (V1C < SL296_ckd_0_v1) I1 = -1;
  V1C = SL296_ckd_0_v1 + (SL296_ckd_0_dv * (Numeric)I1);

  int I2 = (int) ((V2C-SL296_ckd_0_v1) / SL296_ckd_0_dv);

  int NPTC = I2-I1+3;
  if (NPTC > SL296_ckd_0_npt) NPTC = SL296_ckd_0_npt+1;

  V2C = V1C + SL296_ckd_0_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      ostringstream os;
      out3 << "WARNING:\n"  
	   << "  CKD2.2.2 H2O self continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  Numeric SH2OT0[NPTC+addF77fields]; // [cm^3/molecules]
  Numeric SH2OT1[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > SL296_ckd_0_npt) ) 
	{
	  SH2OT0[J] = 0.0e0;   // at T=296 K
	  SH2OT1[J] = 0.0e0;   // at T=260 K
	}
      else
	{
	  SH2OT0[J] = SL296_ckd_0[I];    // at T=296 K
	  SH2OT1[J] = SL260_ckd_0[I];    // at T=260 K 
	}
    }

  // ------------------- subroutine SL296/SL260 ----------------------------
  
  Numeric SFAC = 1.00e0;
  Numeric VS2  = 0.00e0;
  Numeric VS4  = 0.00e0;
  
  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {

      // atmospheric state parameters
      Numeric Tave   = t_abs[i];                                     // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                          // [hPa]
      Numeric Patm   = Pave/PO;                                      // [1]
      Numeric vmrh2o = vmr[i];                                       // [1]
      Numeric Ph2o   = Patm * vmrh2o;                                // [1]
      // second vmr in absCalc multiplied
      Numeric Rh2o   = Patm * (TO/Tave);                             // [1]
      Numeric Tfac   = (Tave-TO)/(260.0-TO);                         // [1]
      Numeric WTOT   = xLosmt * (Pave/1.013000e3) * (2.7300e2/Tave); // [molecules/cm^2]
      Numeric W1     = vmrh2o * WTOT;                                // [molecules/cm^2]
      Numeric XKT    = Tave / 1.4387752e0;                           // = (T*k_B)/(h*c)    
      
      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SH2O = 0.0e0;
	  if (SH2OT0[J] > 0.0e0)
	    {
	      SH2O = SH2OT0[J] * pow( (SH2OT1[J]/SH2OT0[J]), Tfac ); 
	      SFAC = 1.00e0;

	      if ( (VJ >= 700.0e0) && (VJ <= 1200.0e0) )
		{
		  int JFAC = (int)((VJ - 700.0e0)/10.0e0 + 0.00001e0);
		  if ( (JFAC >= 0) && (JFAC <= 50) )
		    SFAC = XFAC[JFAC];
		}

	      // ---------------------------------------------------------
	      // Correction to self continuum (1 SEPT 85); factor of    
	      // 0.78 at 1000 and  .......

	      VS2  = (VJ-V0S) * (VJ-V0S);

	      SFAC = SFAC * 
		( 1.000e0 + 0.3000e0 * (1.000e4 / ((VJ*VJ)                         + 1.000e4))  ) *   
		( 1.000e0 - 0.2333e0 * (ALPHA2  / ((VJ-1050.000e0)*(VJ-1050.000e0) + ALPHA2))   ) *   
		( 1.000e0 - FACTRS   * (ALPHS2  / (VS2+(BETAS*VS2*VS2)+ALPHS2))                 );

	      SH2O = SFAC * SH2O;
	    }

	  // CKD cross section with radiative field [1/cm]
	  // the VMRH2O will be multiplied in absCalc, hence Rh2o does not contain 
	  // VMRH2O as multiplicative term
	  k[J] = W1 * Rh2o * (SH2O*1.000e-20) * RADFN_FUN(VJ,XKT); // [1]

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V >= 0.000e0) && (V < SL296_ckd_0_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      // The factor 100 comes from the conversion from 1/cm to 1/m for
	      // the absorption coefficient
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}



// =================================================================================

// CKD version 2.2.2 H2O foreign continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O foreign continuum according to CKDv.2.2.2    [1/m]
   \param    Cin            strength scaling factor                          [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD version 2.2.2 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-28-08
*/ 
void CKD_222_foreign_h2o( MatrixView          xsec,
			  const Numeric       Cin,
			  const String&       model,
			  ConstVectorView     f_mono,
			  ConstVectorView     p_abs,
			  ConstVectorView     t_abs,
			  ConstVectorView     vmr,
			  ConstVectorView     n2_abs )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKD222"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKDv2.2.2 H2O foreign continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKD222\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the foreign H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // [molecules/cm^3]
  const Numeric T1     =  273.000e0;
  const Numeric TO     =  296.000e0;
  const Numeric PO     = 1013.000e0;

  // CKD2.2.2 foreign H2O continuum correction function parameters
  const Numeric HWSQF  = 330.000e0 * 330.000e0;
  const Numeric BETAF  = 8.000e-11;
  const Numeric V0F    = 1130.000e0;
  const Numeric FACTRF = 0.970e0;
  
  const Numeric V0F2   = 1900.000e0;
  const Numeric HWSQF2 = 150.000e0 * 150.000e0;
  const Numeric BETA2  = 3.000e-6;

  // wavenumber range where CKD H2O foreign continuum is valid
  const Numeric VABS_min = FH2O_ckd_0_v1; // [cm^-1]
  const Numeric VABS_max = FH2O_ckd_0_v2; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD2.2.2 H2O foreign continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << FH2O_ckd_0_v1 << "<->" << FH2O_ckd_0_v2 << "cm^-1\n";
    }


  // ---------------------- subroutine FRN296 ------------------------------

  // retrieve the appropriate array sequence of the foreign continuum
  // arrays of the CKD model.
  Numeric DVC = FH2O_ckd_0_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-FH2O_ckd_0_v1) / FH2O_ckd_0_dv);
  if (V1C < FH2O_ckd_0_v1) I1 = -1;
  V1C = FH2O_ckd_0_v1 + (FH2O_ckd_0_dv * (Numeric)I1);

  int I2 = (int) ((V2C-FH2O_ckd_0_v1) / FH2O_ckd_0_dv);

  int NPTC = I2-I1+3;
  if (NPTC > FH2O_ckd_0_npt) NPTC = FH2O_ckd_0_npt+1;

  V2C = V1C + FH2O_ckd_0_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n" 
	   << "  CKD2.2.2 H2O foreign continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  Numeric FH2OT0[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > FH2O_ckd_0_npt) ) 
	{
	  FH2OT0[J] = 0.0e0;
	}
      else
	{
	  FH2OT0[J] = FH2O_ckd_0[I];
	}
    }

  // ---------------------- subroutine FRN296 ------------------------------
  
  Numeric VF2   = 0.000e0;
  Numeric VF4   = 0.000e0;
  Numeric VF6   = 0.000e0;
  Numeric FSCAL = 0.000e0;
  Numeric FH2O  = 0.000e0;
  
  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {

      // atmospheric state parameters
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrh2o = vmr[i];                                 // [1]
      Numeric ph2o   = vmrh2o * Pave;                          // [hPa]
      Numeric PFRGN  = (Pave/PO) * (1.00000e0 - vmrh2o);       // dry air pressure [hPa]
      Numeric RFRGN  = PFRGN  * (TO/Tave);                     // [hPa]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric W1     = vmrh2o * WTOT;                          // [molecules/cm^2]
      Numeric XKT    = Tave   / 1.4387752;                     // = (T*k_B) / (h*c)
      
      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ = V1C + (DVC * (Numeric)(J-1)); 

	  // CORRECTION TO FOREIGN CONTINUUM
	  VF2   = (VJ-V0F)  * (VJ-V0F);
	  VF6   = VF2 * VF2 * VF2;
	  FSCAL = (1.000e0 - FACTRF*(HWSQF/(VF2+(BETAF*VF6)+HWSQF)));

	  VF2   = (VJ-V0F2) * (VJ-V0F2);
	  VF4   = VF2 * VF2;
	  FSCAL = FSCAL * (1.000e0 - 0.600e0*(HWSQF2/(VF2 + BETA2*VF4 + HWSQF2)));
     
	  FH2O  = FH2OT0[J] * FSCAL;

	  // CKD cross section with radiative field [1/cm]
	  // The VMRH2O will be multiplied in absCalc, hence WTOT and not W1
	  // as multiplicative term
	  k[J] = WTOT * RFRGN * (FH2O*1.000e-20) * RADFN_FUN(VJ,XKT);

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > 0.000e0) && (V < VABS_max) )
	    {
	      // arts CKD2.2.2 foreign H2O continuum cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      // The factor 100 comes from the conversion from (1/cm) to (1/m) 
	      // of the abs. coeff.
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}


// =================================================================================

// CKD version 2.4.2 H2O self continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O self continuum according to CKD_2_4_2  [1/m]
   \param    Cin            strength scaling factor                    [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD version 2.4.2 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-30-10
*/ 
void CKD_242_self_h2o( MatrixView          xsec,
		       const Numeric       Cin,
		       const String&       model,
		       ConstVectorView     f_mono,
		       ConstVectorView     p_abs,
		       ConstVectorView     t_abs,
		       ConstVectorView     vmr,
		       ConstVectorView     n2_abs )
{


  // check the model name about consistency
  if ((model != "user") &&  (model != "CKD242"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKDv2.4.2 H2O self continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKD242\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the self H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt   = 2.686763e19; // [molecules/cm^3]
  const Numeric T1       =  273.0e0;
  const Numeric TO       =  296.0e0;
  const Numeric PO       = 1013.0e0;

  // CKD2.4.2 specific correction functions
  const Numeric V0S1     = 0.000e+00;
  const Numeric HWSQ1    = (1.000e+02 * 1.000e+02);
  const Numeric BETAS1   = 1.000e-04;
  const Numeric FACTRS1  = 0.688e+00;
  
  const Numeric V0S2     = 1.050e+03;
  const Numeric HWSQ2    = (2.000e+02 * 2.000e+02);
  const Numeric FACTRS2  = -0.2333e+00;
  
  const Numeric V0S3     = 1.310e+03;
  const Numeric HWSQ3    = (1.200e+02 * 1.200e+02);
  const Numeric BETAS3   = 5.000e-06;
  const Numeric FACTRS3  = -0.150e+00;
  
  const Numeric XFAC[51] = {  
    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
    1.00000,1.00000,1.00000};

  // wavenumber range where CKD H2O self continuum is valid
  const Numeric VABS_min = SL260_ckd_0_v1; // [cm^-1]
  const Numeric VABS_max = SL260_ckd_0_v2; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD2.4.2 H2O self continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << SL296_ckd_0_v1 << "<->" << SL296_ckd_0_v2 << "cm^-1\n";
    }


  // ------------------- subroutine SL296/SL260 ----------------------------

  if (SL296_ckd_0_v1 != SL260_ckd_0_v1)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.4.2 H2O self continuum:\n"
	 << "parameter V1 not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (SL296_ckd_0_v2 != SL260_ckd_0_v2)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.4.2 H2O self continuum:\n"
	 << "parameter V2 not the same for different ref. temperatures.\n";
	throw runtime_error(os.str());
    }
  if (SL296_ckd_0_dv != SL260_ckd_0_dv)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.4.2 H2O self continuum:\n"
	 << "parameter DV not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (SL296_ckd_0_npt != SL260_ckd_0_npt)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD2.4.2 H2O self continuum:\n"
	 << "parameter NPT not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  
  // retrieve the appropriate array sequence of the self continuum
  // arrays of the CKD model.
  Numeric DVC = SL296_ckd_0_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-SL296_ckd_0_v1) / SL296_ckd_0_dv);
  if (V1C < SL296_ckd_0_v1) I1 = -1;
  V1C = SL296_ckd_0_v1 + (SL296_ckd_0_dv * (Numeric)I1);

  int I2 = (int) ((V2C-SL296_ckd_0_v1) / SL296_ckd_0_dv);

  int NPTC = I2-I1+3;
  if (NPTC > SL296_ckd_0_npt) NPTC = SL296_ckd_0_npt+1;

  V2C = V1C + SL296_ckd_0_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      ostringstream os;
      out3 << "WARNING:\n"  
	   << "  CKDv2.4.2 H2O self continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  Numeric SH2OT0[NPTC+addF77fields]; // [cm^3/molecules]
  Numeric SH2OT1[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > SL296_ckd_0_npt) ) 
	{
	  SH2OT0[J] = 0.0e0;   // at T=296 K
	  SH2OT1[J] = 0.0e0;   // at T=260 K
	}
      else
	{
	  SH2OT0[J] = SL296_ckd_0[I];    // at T=296 K
	  SH2OT1[J] = SL260_ckd_0[I];    // at T=260 K 
	}
    }

  // ------------------- subroutine SL296/SL260 ----------------------------
  
  Numeric SFAC = 1.00e0;
  Numeric VS2  = 0.00e0;
  Numeric VS4  = 0.00e0;
  
  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {

      // atmospheric state parameters
      Numeric Tave   = t_abs[i];                                     // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                          // [hPa]
      Numeric Patm   = Pave/PO;                                      // [1]
      Numeric vmrh2o = vmr[i];                                       // [1]
      Numeric Ph2o   = Patm * vmrh2o;                                // [1]
      // second vmr in absCalc multiplied
      Numeric Rh2o   = Patm * (TO/Tave);                             // [1]
      Numeric Tfac   = (Tave-TO)/(260.0-TO);                         // [1]
      Numeric WTOT   = xLosmt * (Pave/1.013000e3) * (2.7300e2/Tave); // [molecules/cm^2]
      Numeric W1     = vmrh2o * WTOT;                                // [molecules/cm^2]
      Numeric XKT    = Tave / 1.4387752e0;                           // = (T*k_B)/(h*c)   
      
      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SH2O = 0.0e0;
	  if (SH2OT0[J] > 0.0e0)
	    {
	      SH2O = SH2OT0[J] * pow( (SH2OT1[J]/SH2OT0[J]), Tfac ); 
	      SFAC = 1.00e0;

	      if ( (VJ >= 700.0e0) && (VJ <= 1200.0e0) )
		{
		  int JFAC = (int)((VJ - 700.0e0)/10.0e0 + 0.00001e0);
		  if ( (JFAC >= 0) && (JFAC <= 50) )
		    SFAC = XFAC[JFAC];
		}

	      // ---------------------------------------------------------
	      // Correction to self continuum (1 SEPT 85); factor of    
	      // 0.78 at 1000 and  .......

	      VS2  = (VJ-V0S1) * (VJ-V0S1);
	      VS4  = VS2*VS2;
	      SFAC = SFAC * 
		(1.000e0 + FACTRS1*(HWSQ1/((VJ*VJ)+(BETAS1*VS4)+HWSQ1))); 

	      VS2  = (VJ-V0S2) * (VJ-V0S2);
	      SFAC = SFAC *
		(1.000e0 + FACTRS2*(HWSQ2/(VS2+HWSQ2)));

	      VS2  = (VJ-V0S3) * (VJ-V0S3);
	      VS4  = VS2*VS2;
	      SFAC = SFAC *
		(1.000e0 + FACTRS3*(HWSQ3/(VS2+(BETAS3*VS4)+HWSQ3))); 
	      
	      SH2O = SFAC * SH2O;
	    }

	  // CKD cross section with radiative field [1/cm]
	  // The VMRH2O will be multiplied in absCalc, hence Rh2o does not contain 
	  // VMRH2O as multiplicative term
	  k[J] = W1 * Rh2o * (SH2O*1.000e-20) * RADFN_FUN(VJ,XKT);

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V >= 0.000e0) && (V < SL296_ckd_0_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      // The factor 100 comes from the conversion from 1/cm to 1/m for
	      // the absorption coefficient
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}


// =================================================================================

// CKD version 2.4.2 H2O foreign continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O foreign continuum according to CKDv.2.4.2    [1/m]
   \param    Cin            strength scaling factor                          [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD version 2.4.2 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-28-08
*/ 
void CKD_242_foreign_h2o( MatrixView          xsec,
			  const Numeric       Cin,
			  const String&       model,
			  ConstVectorView     f_mono,
			  ConstVectorView     p_abs,
			  ConstVectorView     t_abs,
			  ConstVectorView     vmr,
			  ConstVectorView     n2_abs )
{


  // check the model name about consistency
  if ((model != "user") &&  (model != "CKD242"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKDv2.4.2 H2O foreign continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKD242\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the foreign H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // [molecules/cm^3]
  const Numeric T1     =  273.0e0;
  const Numeric TO     =  296.0e0;
  const Numeric PO     = 1013.0e0;

  // CKD2.4.2 foreign H2O continuum correction function parameters
  const Numeric V0F1     = 350.000e0;
  const Numeric HWSQF1   = 200.000e0 * 200.000e0;
  const Numeric BETAF1   = 5.000e-9 ;
  const Numeric FACTRF1  = -0.700e0;
  
  const Numeric V0F1a    = 630.000e0;
  const Numeric HWSQF1a  = 65.000e0*65.000e0;
  const Numeric BETAF1a  = 2.000e-08 ;
  const Numeric FACTRF1a = 0.750e0;
  
  const Numeric V0F2     = 1130.000e0;
  const Numeric HWSQF2   = 330.000e0 * 330.000e0;
  const Numeric BETAF2   = 8.000e-11;
  const Numeric FACTRF2  = -0.970e0;
  
  const Numeric V0F3     = 1975.000e0;
  const Numeric HWSQF3   = 250.000e0 * 250.000e0;
  const Numeric BETAF3   = 5.000e-06;
  const Numeric FACTRF3  = -0.650e0;

  // wavenumber range where CKD H2O foreign continuum is valid
  const Numeric VABS_min = FH2O_ckd_0_v1; // [cm^-1]
  const Numeric VABS_max = FH2O_ckd_0_v2; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKDv2.4.2 H2O foreign continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << FH2O_ckd_0_v1 << "<->" << FH2O_ckd_0_v2 << "cm^-1\n";
    }


  // ---------------------- subroutine FRN296 ------------------------------

  // retrieve the appropriate array sequence of the foreign continuum
  // arrays of the CKD model.
  Numeric DVC = FH2O_ckd_0_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-FH2O_ckd_0_v1) / FH2O_ckd_0_dv);
  if (V1C < FH2O_ckd_0_v1) I1 = -1;
  V1C = FH2O_ckd_0_v1 + (FH2O_ckd_0_dv * (Numeric)I1);

  int I2 = (int) ((V2C-FH2O_ckd_0_v1) / FH2O_ckd_0_dv);

  int NPTC = I2-I1+3;
  if (NPTC > FH2O_ckd_0_npt) NPTC = FH2O_ckd_0_npt+1;

  V2C = V1C + FH2O_ckd_0_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n" 
	   << "  CKDv2.4.2 H2O foreign continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  Numeric FH2OT0[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > FH2O_ckd_0_npt) ) 
	{
	  FH2OT0[J] = 0.0e0;
	}
      else
	{
	  FH2OT0[J] = FH2O_ckd_0[I];
	}
    }

  // ---------------------- subroutine FRN296 ------------------------------
  
  Numeric VF2   = 0.000e0;
  Numeric VF4   = 0.000e0;
  Numeric VF6   = 0.000e0;
  Numeric FSCAL = 0.000e0;
  Numeric FH2O  = 0.000e0;
  
  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {

      // atmospheric state parameters
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrh2o = vmr[i];                                 // [1]
      Numeric ph2o   = vmrh2o * Pave;                          // [hPa]
      Numeric PFRGN  = (Pave/PO) * (1.00000e0 - vmrh2o);       // dry air pressure [hPa]
      Numeric RFRGN  = PFRGN  * (TO/Tave);                     // [hPa]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric W1     = vmrh2o * WTOT;                          // [molecules/cm^2]
      Numeric XKT    = Tave   / 1.4387752;                     // = (T*k_B) / (h*c)    
      
      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ = V1C + (DVC * (Numeric)(J-1)); 

	  // CORRECTION TO FOREIGN CONTINUUM
	  VF2   = (VJ-V0F1) * (VJ-V0F1);
	  VF6   = VF2 * VF2 * VF2;
	  FSCAL = (1.000e0 + FACTRF1*(HWSQF1/(VF2+(BETAF1*VF6)+HWSQF1)));

	  VF2   = (VJ-V0F1a) * (VJ-V0F1a);
	  VF6   = VF2 * VF2  * VF2;
	  FSCAL = FSCAL * 
	          (1.000e0 + FACTRF1a*(HWSQF1a/(VF2+(BETAF1a*VF6)+HWSQF1a)));

	  VF2   = (VJ-V0F2) * (VJ-V0F2);
	  VF6   = VF2 * VF2 * VF2;
	  FSCAL = FSCAL * 
	          (1.000e0 + FACTRF2*(HWSQF2/(VF2+(BETAF2*VF6)+HWSQF2)));

	  VF2   = (VJ-V0F3) * (VJ-V0F3);
	  VF4   = VF2 * VF2;
	  FSCAL = FSCAL * 
		  (1.000e0 + FACTRF3*(HWSQF3/(VF2+BETAF3*VF4+HWSQF3)));
     
	  FH2O  = FH2OT0[J] * FSCAL;

	  // CKD cross section without radiative field
	  // The VMRH2O will be multiplied in absCalc, hence WTOT and not W1
	  // as multiplicative term
	  k[J] = WTOT * RFRGN * (FH2O*1.000e-20) * RADFN_FUN(VJ,XKT);

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V >= 0.000e0) && (V < VABS_max) )
	    {
	      // arts CKD2.4.2 foreign H2O continuum cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}


// =================================================================================

// CKD version MT 1.00 H2O self continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O self continuum according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author   Thomas Kuhn
   \date     2002-28-08
*/ 
void CKD_mt_100_self_h2o( MatrixView          xsec,
			  const Numeric       Cin,
			  const String&       model,
			  ConstVectorView     f_mono,
			  ConstVectorView     p_abs,
			  ConstVectorView     t_abs,
			  ConstVectorView     vmr,
			  ConstVectorView     n2_abs )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 H2O self continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the self H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt =    2.68675e19; // [molecules/cm^3]
  const Numeric T1     =  273.000e0;    // [K]
  const Numeric TO     =  296.000e0;    // [K]
  const Numeric PO     = 1013.000e0;    // [hPa]

  const Numeric XFACREV[15] = 
    {1.003, 1.009, 1.015, 1.023, 1.029,1.033,
     1.037, 1.039, 1.040, 1.046, 1.036,1.027,
     1.01,  1.002, 1.00};

  // wavenumber range where CKD H2O self continuum is valid
  const Numeric VABS_min = -2.000e1; // [cm^-1]
  const Numeric VABS_max =  2.000e4; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
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

  Numeric SH2OT0[NPTC+addF77fields]; // [cm^3/molecules]
  Numeric SH2OT1[NPTC+addF77fields]; // [cm^3/molecules]

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
  
  Numeric SFAC = 1.00e0;  

  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {

      // atmospheric state parameters
      Numeric Tave   = t_abs[i];                                     // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                          // [hPa]
      Numeric Patm   = Pave/PO;                                      // [1]
      Numeric vmrh2o = vmr[i];                                       // [1]
      Numeric Ph2o   = Patm * vmrh2o;                                // [1]
      // second vmr in absCalc multiplied
      Numeric Rh2o   = Patm * (TO/Tave);                             // [1]
      Numeric Tfac   = (Tave-TO)/(260.0-TO);                         // [1]
      Numeric WTOT   = xLosmt * (Pave/1.013000e3) * (2.7300e2/Tave); // [molecules/cm^2]
      Numeric W1     = vmrh2o * WTOT;                                // [molecules/cm^2]
      Numeric XKT    = Tave / 1.4387752e0;                           // = (T*k_B)/(h*c)    

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SH2O = 0.0e0;
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

	  // CKD cross section with radiative field [1/cm]
	  // The VMRH2O will be multiplied in absCalc, hence Rh2o does not contain 
	  // VMRH2O as multiplicative term
	  k[J] = W1 * Rh2o * (SH2O*1.000e-20) * RADFN_FUN(VJ,XKT);

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > 0.000e0) && (V < SL296_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      // The factor 100 comes from the conversion from 1/cm to 1/m for
	      // the absorption coefficient
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

// =================================================================================

// CKD version MT 1.00 H2O foreign continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O foreign continuum according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                          [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-28-08
*/ 
void CKD_mt_100_foreign_h2o( MatrixView          xsec,
			     const Numeric       Cin,
			     const String&       model,
			     ConstVectorView     f_mono,
			     ConstVectorView     p_abs,
			     ConstVectorView     t_abs,
			     ConstVectorView     vmr,
			     ConstVectorView     n2_abs )
{


  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 H2O foreign continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the foreign H2O cont. absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.68675e19; // [molecules/cm^3]
  const Numeric T1     =  273.000e0;
  const Numeric TO     =  296.000e0;
  const Numeric PO     = 1013.000e0;

  // wavenumber range where CKD H2O self continuum is valid
  const Numeric VABS_min = -2.000e1; // [cm^-1]
  const Numeric VABS_max =  2.000e4; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
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

  Numeric FH2OT0[NPTC+addF77fields]; // [cm^3/molecules]

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
	}
    }

  // ---------------------- subroutine FRN296 ------------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      // atmospheric state parameters
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrh2o = vmr[i];                                 // [1]
      Numeric ph2o   = vmrh2o * Pave;                          // [hPa]
      Numeric PFRGN  = (Pave/PO) * (1.00000e0 - vmrh2o);       // dry air pressure [hPa]
      Numeric RFRGN  = PFRGN  * (TO/Tave);                     // [hPa]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric W1     = vmrh2o * WTOT;                          // [molecules/cm^2]
      Numeric XKT    = Tave   / 1.4387752;                     // = (T*k_B) / (h*c)    
      
      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric FH2O = FH2OT0[J];

	  // CKD cross section with radiative field [1/cm]
	  // The VMRH2O will be multiplied in absCalc, hence WTOT and not W1
	  // as multiplicative term
	  k[J] = WTOT * RFRGN * (FH2O*1.000e-20) * RADFN_FUN(VJ,XKT);

	}
      
      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V >= 0.000e0) && (V < VABS_max) )
	    {
	      // arts CKD_MT.100 cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      // The factor 100 comes from the conversion from (1/cm) to (1/m) 
	      // of the abs. coeff.
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

//

// =================================================================================

// CKD version 2.4.1 CO2 continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            CO2 continuum according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                          [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            CO2 volume mixing ratio profile      [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD version 2.4.1 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html

   \author Thomas Kuhn
   \date 2002-28-08
 */ 
void CKD_241_co2( MatrixView         xsec,
		 const Numeric       Cin,
		 const String&       model,
		 ConstVectorView     f_mono,
		 ConstVectorView     p_abs,
		 ConstVectorView     t_abs,
		 ConstVectorView     vmr )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKD241"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKDv2.4.1 CO2 continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKD241\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the CO2 absorption
  Numeric  ScalingFac = 0.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }
  else
    {
      ScalingFac = 1.0000e0;
    }

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // [molecules/cm^3]
  const Numeric T1     =  273.0e0;
  const Numeric TO     =  296.0e0;
  const Numeric PO     = 1013.0e0;

  // wavenumber range where CKD CO2 continuum is valid
  const Numeric VABS_min = -2.000e1; // [cm^-1]
  const Numeric VABS_max =  1.000e4; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKDv2.4.1 CO2 continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << FCO2_ckd_mt_100_v1 << "<->" << FCO2_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ---------------------- subroutine FRNCO2 ------------------------------

  // retrieve the appropriate array sequence of the CO2 continuum
  // arrays of the CKD model.
  Numeric DVC = FCO2_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-FCO2_ckd_mt_100_v1) / FCO2_ckd_mt_100_dv);
  if (V1C < FCO2_ckd_mt_100_v1) I1 = -1;
  V1C = FCO2_ckd_mt_100_v1 + (FCO2_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-FCO2_ckd_mt_100_v1) / FCO2_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > FCO2_ckd_mt_100_npt) NPTC = FCO2_ckd_mt_100_npt+1;

  V2C = V1C + FCO2_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n"  
	   << "  CKDv2.4.1 CO2 continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  Numeric FCO2T0[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > FCO2_ckd_mt_100_npt) ) 
	{
	  FCO2T0[J] = 0.0e0;
	}
      else
	{
	  FCO2T0[J] = FCO2_ckd_mt_100[I];
	}
    }

  // ---------------------- subroutine FRNCO2 ------------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrco2 = vmr[i];                                 // [1]
      Numeric Rhoave = (Pave/PO) * (TO/Tave);                  // [hPa]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric XKT    = Tave / 1.4387752;                       // = (T*k_B) / (h*c)    
      

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric FCO2 = FCO2T0[J];

	  // CKD cross section times number density with radiative field [1]
	  // the VMRCO2 will be multiplied in absCalc
	  k[J] = ((WTOT * Rhoave) * (FCO2*1.000e-20) * RADFN_FUN(VJ,XKT));

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > 0.000e0) && (V < FCO2_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}


// =================================================================================


// CKD version MT 1.00 CO2 continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            CO2 continuum according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                          [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            CO2 volume mixing ratio profile      [1]

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html

   \author   Thomas Kuhn
   \date     2002-28-08
 */ 
void CKD_mt_co2( MatrixView          xsec,
		 const Numeric       Cin,
		 const String&       model,
		 ConstVectorView     f_mono,
		 ConstVectorView     p_abs,
		 ConstVectorView     t_abs,
		 ConstVectorView     vmr)
{


  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT.1.00 CO2 continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the CO2 absorption
  Numeric  ScalingFac = 0.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }
  else
    {
      ScalingFac = 1.0000e0;
    }

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // [molecules/cm^3]
  const Numeric T1     =  273.0e0;
  const Numeric TO     =  296.0e0;
  const Numeric PO     = 1013.0e0;

  // wavenumber range where CKD CO2 continuum is valid
  const Numeric VABS_min = FCO2_ckd_mt_100_v1; // [cm^-1]
  const Numeric VABS_max = FCO2_ckd_mt_100_v2; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD_MT 1.00 CO2 continuum:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << FCO2_ckd_mt_100_v1 << "<->" << FCO2_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ---------------------- subroutine FRNCO2 ------------------------------

  // retrieve the appropriate array sequence of the CO2 continuum
  // arrays of the CKD model.
  Numeric DVC = FCO2_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-FCO2_ckd_mt_100_v1) / FCO2_ckd_mt_100_dv);
  if (V1C < FCO2_ckd_mt_100_v1) I1 = -1;
  V1C = FCO2_ckd_mt_100_v1 + (FCO2_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-FCO2_ckd_mt_100_v1) / FCO2_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > FCO2_ckd_mt_100_npt) NPTC = FCO2_ckd_mt_100_npt+1;

  V2C = V1C + FCO2_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n"  
	   << "  CKD_MT 1.00 CO2 continuum:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.";
      return;
    }

  Numeric FCO2T0[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > FCO2_ckd_mt_100_npt) ) 
	{
	  FCO2T0[J] = 0.0e0;
	}
      else
	{
	  FCO2T0[J] = FCO2_ckd_mt_100[I];
	}
    }

  // ---------------------- subroutine FRNCO2 ------------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrco2 = vmr[i];                                 // [1]
      Numeric Rhoave = (Pave/PO) * (TO/Tave);                  // [hPa]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric XKT    = Tave / 1.4387752;                       // = (T*k_B) / (h*c)    
      

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ   = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric FCO2 = FCO2T0[J];

	  // continuum has been increased in the nu2 band by a factor of 7
	  if ( (VJ > 500.0e0) && (VJ < 900.0e0) )
	    {
	      FCO2 = 7.000e0 * FCO2; 
	    }

	  // CKD cross section times number density with radiative field [1]
	  // the VMRCO2 will be multiplied in absCalc
	  k[J] = ((WTOT * Rhoave) * (FCO2*1.000e-20) * RADFN_FUN(VJ,XKT));

	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > 0.000e0) && (V < FCO2_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}


// =================================================================================
// CKD version MT 1.00 N2-N2 collision induced absorption (rotational band)
// Model reference:
//  Borysow, A, and L. Frommhold, 
//  "Collision-induced rototranslational absorption spectra of N2-N2
//  pairs for temperatures from 50 to 300 K", The
//  Astrophysical Journal, 311, 1043-1057, 1986.
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            N2-N2 CIA rot. band according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            N2 volume mixing ratio profile       [1]

   \remark   Borysow, A, and L. Frommhold,<br> 
             Collision-induced rototranslational absorption spectra of N2-N2
             pairs for temperatures from 50 to 300 K,<br> 
             The Astrophysical Journal, 311, 1043-1057, 1986.

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-28-08
 */ 
void CKD_mt_CIArot_n2( MatrixView         xsec,
		      const Numeric       Cin,
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     vmr )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 N2 CIA rotational band:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the N2-N2 CIA rot. band absorption
  Numeric  ScalingFac = 0.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }
  else
    {
      ScalingFac = 1.0000e0;
    }

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // Loschmidt Number [molecules/cm^3]
  const Numeric T1     =  273.0e0;
  const Numeric TO     =  296.0e0;
  const Numeric PO     = 1013.0e0;


  // wavenumber range where CKD H2O self continuum is valid
  const Numeric VABS_min = -1.000e1; // [cm^-1]
  const Numeric VABS_max =  3.500e2; // [cm^-1]


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < VABS_min) || (V1ABS > VABS_max) ||
       (V2ABS < VABS_min) || (V2ABS > VABS_max) )
    {
      out3  << "WARNING:\n"
            << "  CKD_MT 1.00 N2-N2 CIA rotational band:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << N2N2_CT296_ckd_mt_100_v1 << "<->" << N2N2_CT220_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ------------------- subroutine N2R296/N2R220 ----------------------------

  if (N2N2_CT296_ckd_mt_100_v1 != N2N2_CT220_ckd_mt_100_v1)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 N2-N2 CIA rotational band:\n"
	 << "parameter V1 not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (N2N2_CT296_ckd_mt_100_v2 != N2N2_CT220_ckd_mt_100_v2)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 N2-N2 CIA rotational band:\n"
	 << "parameter V2 not the same for different ref. temperatures.\n";
	throw runtime_error(os.str());
    }
  if (N2N2_CT296_ckd_mt_100_dv != N2N2_CT220_ckd_mt_100_dv)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 N2-N2 CIA rotational band:\n"
	 << "parameter DV not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  if (N2N2_CT296_ckd_mt_100_npt != N2N2_CT220_ckd_mt_100_npt)
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT 1.00 N2-N2 CIA rotational band:\n"
	 << "parameter NPT not the same for different ref. temperatures.\n";
      throw runtime_error(os.str());
    }
  
  // retrieve the appropriate array sequence of the self continuum
  // arrays of the CKD model.
  Numeric DVC = N2N2_CT296_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-N2N2_CT296_ckd_mt_100_v1) / N2N2_CT296_ckd_mt_100_dv);
  if (V1C < N2N2_CT296_ckd_mt_100_v1) I1 = -1;
  V1C    = N2N2_CT296_ckd_mt_100_v1 + (N2N2_CT296_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-N2N2_CT296_ckd_mt_100_v1) / N2N2_CT296_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > N2N2_CT296_ckd_mt_100_npt) NPTC = N2N2_CT296_ckd_mt_100_npt+1;

  V2C    = V1C + N2N2_CT296_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n" 
	   << "  CKD_MT 1.00 N2-N2 CIA rotational band:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.\n";
      return;
    }

  Numeric C0[NPTC+addF77fields]; // [cm^3/molecules]
  Numeric C1[NPTC+addF77fields]; // [cm^3/molecules]

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      Index I = I1+J;
      if ( (I < 1) || (I > N2N2_CT296_ckd_mt_100_npt) ) 
	{
	  C0[J] = 0.0e0;   // at T=296 K
	  C1[J] = 0.0e0;   // at T=260 K
	}
      else
	{
	  C0[J] = N2N2_CT296_ckd_mt_100[I];    // at T=296 K
	  C1[J] = N2N2_CT220_ckd_mt_100[I];    // at T=260 K 
	}
    }

  // ------------------- subroutine N2R296/N2R220 ----------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrn2  = vmr[i];                                 // [1]
      Numeric facfac = vmrn2 * (Pave/PO) * (Pave/PO) * 
                               (T1/Tave) * (T1/Tave);

      Numeric XKT    = Tave / 1.4387752;                       // = (T*k_B) / (h*c)    
      Numeric Tfac   = (Tave - TO) / (220.0e0 - TO);

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  k[J] = 0.000e0;
	  Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SN2 = 0.000e0;
	  if ( (C0[J] > 0.000e0) && (C1[J] > 0.000e0) )
	    {
	      SN2 = facfac* C0[J] * pow( (C1[J]/C0[J]), Tfac ); 
	    }

	  // CKD cross section with radiative field
	  k[J] = SN2 * RADFN_FUN(VJ,XKT); // [1]
	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > 0.000e0) && (V < N2N2_CT220_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

// =================================================================================

// CKD version MT 1.00 N2-N2 collision induced absorption (fundamental band)
// Model reference:
//  version_1 of the Nitrogen Collision Induced Fundamental
//  Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and 
//  J._M. Hartmann, Infrared collision-induced absorption by 
//  N2 near 4.3 microns for atmospheric applications: 
//  Measurements and emprirical modeling, Appl. Optics, 35, 
//  5911-5917, (1996).
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            N2-N2 CIA fundamental band according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            N2 volume mixing ratio profile       [1]

   \remark   Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and 
             J._M. Hartmann,<br> 
             Infrared collision-induced absorption by 
             N2 near 4.3 microns for atmospheric applications: 
             Measurements and emprirical modeling, <br> 
             Appl. Optics, 35, 5911-5917, (1996)

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             

   \author Thomas Kuhn
   \date 2002-28-08
 */ 
void CKD_mt_CIAfun_n2( MatrixView         xsec,
		      const Numeric       Cin,
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     vmr )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 N2 CIA fundamental band:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the N2-N2 CIA fundamental band absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // Loschmidt Number [molecules/cm^3]
  const Numeric T1     =  273.0e0;
  const Numeric TO     =  296.0e0;
  const Numeric PO     = 1013.0e0;
  const Numeric a1     = 0.8387e0;
  const Numeric a2     = 0.0754e0;


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < N2N2_N2F_ckd_mt_100_v1) || (V1ABS > N2N2_N2F_ckd_mt_100_v2) ||
       (V2ABS < N2N2_N2F_ckd_mt_100_v1) || (V2ABS > N2N2_N2F_ckd_mt_100_v2) )
    {
      out3  << "WARNING:\n"
            << "   CKD_MT 1.00 N2-N2 CIA fundamental band:\n"
	    << "   input frequency vector exceeds range of model validity\n"
	    << "  " << N2N2_N2F_ckd_mt_100_v1 << "<->" << N2N2_N2F_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ------------------- subroutine N2_VER_1 ----------------------------
  
  // retrieve the appropriate array sequence of the self continuum
  // arrays of the CKD model.
  Numeric DVC = N2N2_N2F_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-N2N2_N2F_ckd_mt_100_v1) / N2N2_N2F_ckd_mt_100_dv);
  if (V1C < N2N2_N2F_ckd_mt_100_v1) I1 = -1;
  V1C    = N2N2_N2F_ckd_mt_100_v1 + (N2N2_N2F_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-N2N2_N2F_ckd_mt_100_v1) / N2N2_N2F_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > N2N2_N2F_ckd_mt_100_npt) NPTC = N2N2_N2F_ckd_mt_100_npt+1;

  V2C    = V1C + N2N2_N2F_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n"  
	   << "  CKD_MT 1.00 N2-N2 CIA fundamental band:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n";
      return;
    }

  Numeric xn2[NPTC+addF77fields];
  Numeric xn2t[NPTC+addF77fields];

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      xn2[J]  = 0.000e0;
      xn2t[J] = 0.000e0;
      Index I = I1+J;
      if ( (I > 0) && (I <= N2N2_N2F_ckd_mt_100_npt) ) 
	{
	  xn2[J]  = N2N2_N2F_ckd_mt_100[I];
	  xn2t[J] = N2N2_N2Ft_ckd_mt_100[I];
	}
    }

  // ------------------- subroutine N2_VER_1 ----------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmrn2  = vmr[i];                                 // [1]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric tau_fac= WTOT   * (Pave/PO) * (T1/Tave); 

      Numeric XKT    = Tave / 1.4387752e0;                     // = (T*k_B) / (h*c)    

      Numeric Tfac   = (Tave - TO) / (220.0e0 - TO);           // [1]
      Numeric xktfac = (1.000e0/TO) - (1.000e0/Tave);          // [1/K]
      Numeric factor = 0.000e0;
      if (vmrn2 > VMRCalcLimit)
	{
	  factor = (1.000e0 / xLosmt) * (1.000e0/vmrn2) * (a1 - a2*(Tave/TO));
	}

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.000e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  k[J] = 0.000e0;
	  Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SN2 = 0.000e0;
	  if (xn2[J] > 0.000e0)
	    {
	      Numeric C0 = factor * xn2[J] * exp(xn2t[J]*xktfac) / VJ;
	      SN2 = tau_fac * C0;
	    }

	  // CKD cross section with radiative field
	  k[J] = SN2 * RADFN_FUN(VJ,XKT); // [1/cm]
	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > N2N2_N2F_ckd_mt_100_v1) && (V < N2N2_N2F_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

// =================================================================================

// CKD version MT 1.00 O2-O2 collision induced absorption (fundamental band)
// Model reference:
// F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, 
// J.-M. Hartmann, Ch. Boulet, 
// "Infrared collision-induced absorption by O2 near 6.4 microns for
// atmospheric applications: measurements and emprirical modeling", 
// Appl. Optics, 35, 5911-5917, (1996).
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            O2-O2 CIA fundamental band according to CKD_MT 1.00   [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            O2 volume mixing ratio profile       [1]

   \remark   F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, 
             J.-M. Hartmann, Ch. Boulet,<br>
             Infrared collision-induced absorption by O2 near 6.4 microns for
             atmospheric applications: measurements and emprirical modeling,<br>
             Appl. Optics, 35, 5911-5917, (1996).

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html
             
   \author Thomas Kuhn
   \date 2002-28-08
 */ 
void CKD_mt_CIAfun_o2( MatrixView         xsec,
		      const Numeric       Cin,
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     vmr )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 O2 CIA fundamental band:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the O2-O2 CIA fundamental band absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt = 2.686763e19; // Loschmidt Number [molecules/cm^3]
  const Numeric T1     =  273.0e0;
  const Numeric TO     =  296.0e0;
  const Numeric PO     = 1013.0e0;


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < O2O2_O2F_ckd_mt_100_v1) || (V1ABS > O2O2_O2F_ckd_mt_100_v2) ||
       (V2ABS < O2O2_O2F_ckd_mt_100_v1) || (V2ABS > O2O2_O2F_ckd_mt_100_v2) )
    {
      out3  << "WARNING:\n"
            << "  CKD_MT 1.00 O2-O2 CIA fundamental band:\n"
	    << "  input frequency vector exceeds range of model validity\n"
	    << "  " << O2O2_O2F_ckd_mt_100_v1 << "<->" << O2O2_O2F_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ------------------- subroutine O2_VER_1 ----------------------------
  
  // retrieve the appropriate array sequence of the CKD model array.
  Numeric DVC = O2O2_O2F_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-O2O2_O2F_ckd_mt_100_v1) / O2O2_O2F_ckd_mt_100_dv);
  if (V1C < O2O2_O2F_ckd_mt_100_v1) I1 = -1;
  V1C    = O2O2_O2F_ckd_mt_100_v1 + (O2O2_O2F_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-O2O2_O2F_ckd_mt_100_v1) / O2O2_O2F_ckd_mt_100_dv);

  int NPTC = I2-I1+3;
  if (NPTC > O2O2_O2F_ckd_mt_100_npt) NPTC = O2O2_O2F_ckd_mt_100_npt+1;

  V2C    = V1C + O2O2_O2F_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n"
	   << "  CKD_MT 1.00 O2 CIA fundamental band:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.\n";
      return;
    }

  Numeric xo2[NPTC+addF77fields];
  Numeric xo2t[NPTC+addF77fields];

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      xo2[J]  = 0.000e0;
      xo2t[J] = 0.000e0;
      Index I = I1+J;
      if ( (I > 0) && (I <= O2O2_O2F_ckd_mt_100_npt) ) 
	{
	  xo2[J]  = O2O2_O2Fo_ckd_mt_100[I];
	  xo2t[J] = O2O2_O2Ft_ckd_mt_100[I];
	}
    }

  // ------------------- subroutine O2_VER_1 ----------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmro2  = vmr[i];                                 // [1]
      Numeric WTOT   = xLosmt * (Pave/PO) * (T1/Tave);         // [molecules/cm^2]
      Numeric tau_fac= WTOT * (Pave/PO) * (T1/Tave);

      Numeric XKT    = Tave / 1.4387752;                       // = (T*k_B) / (h*c)    

      Numeric xktfac = (1.000e0/TO) - (1.000e0/Tave);          // [1/K]
      Numeric factor = (1.000e0 / xLosmt);

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SO2 = 0.0e0;
	  if (xo2[J] > 0.0e0)
	    {
	      Numeric C0 = factor * xo2[J] * exp(xo2t[J]*xktfac) / VJ;
	      SO2 = tau_fac * C0;
	    }

	  // CKD cross section without radiative field
	  k[J] = SO2 * RADFN_FUN(VJ,XKT); // [1]
	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > O2O2_O2F_ckd_mt_100_v1) && (V < O2O2_O2F_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

// =================================================================================

// CKD version MT 1.00 O2 v0<-v0 band absorption
// Model reference:
// CKD_MT 1.00 implementation of oxygen collision induced fundamental model of 
// O2 continuum formulated by 
//   Mate et al. over the spectral region 7550-8486 cm-1: 
//   B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,
//   "Absolute Intensities for the O2 1.27 micron
//   continuum absorption",  
//   J. Geophys. Res., 104, 30,585-30,590, 1999. 
//
// The units of these continua coefficients are  1 / (amagat_O2*amagat_air)
//
// Also, refer to the paper "Observed  Atmospheric
// Collision Induced Absorption in Near Infrared Oxygen Bands",
// Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
// Journal of Geophysical Research (1997).
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            O2 v0<-v0 band according to CKD_MT 1.00  [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            O2 volume mixing ratio profile       [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \remark   B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,<br> 
             Absolute Intensities for the O2 1.27 micron continuum absorption,<br>   
             J. Geophys. Res., 104, 30,585-30,590, 1999. 

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
             http://www.rtweb.aer.com/continuum_frame.html<br>
             <br> 
	     Oxygen band absorption model for the \f$a^1\Delta_g\f$ 
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
             \f$X^3\Sigma^-_g\f$ band system considering the 
             \f$\nu=0\f$
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
             \f$\nu=0\f$
             transitions.
             
   \author   Thomas Kuhn
   \date     2002-28-08
 */ 
void CKD_mt_v0v0_o2( MatrixView          xsec,
		     const Numeric       Cin,
		     const String&       model,
		     ConstVectorView     f_mono,
		     ConstVectorView     p_abs,
		     ConstVectorView     t_abs,
		     ConstVectorView     vmr,
		     ConstVectorView     n2_abs)
{


  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 O2 band at 1.27 micrometer:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the O2 v0<-v0 band absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    };

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt    = 2.686763e19; // Loschmidt Number [molecules/cm^3]
  const Numeric T1        =  273.0e0;
  const Numeric TO        =  296.0e0;
  const Numeric PO        = 1013.0e0;

  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < O2_00_ckd_mt_100_v1) || (V1ABS > O2_00_ckd_mt_100_v2) ||
       (V2ABS < O2_00_ckd_mt_100_v1) || (V2ABS > O2_00_ckd_mt_100_v2) )
    {
      out3  << "WARNING:\n"
            << "   CKD_MT 1.00 O2 v0<-v0 band:\n"
	    << "   input frequency vector exceeds range of model validity\n"
	    << "  " << O2_00_ckd_mt_100_v1 << "<->" << O2_00_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ------------------- subroutine O2INF1 ----------------------------
  
  // retrieve the appropriate array sequence of the CKD model array.
  Numeric DVC = O2_00_ckd_mt_100_dv;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int I1 = (int) ((V1C-O2_00_ckd_mt_100_v1) / O2_00_ckd_mt_100_dv);
  if (V1C < O2_00_ckd_mt_100_v1) I1 = I1-1;
  V1C    = O2_00_ckd_mt_100_v1 + (O2_00_ckd_mt_100_dv * (Numeric)I1);

  int I2 = (int) ((V2C-O2_00_ckd_mt_100_v1) / O2_00_ckd_mt_100_dv);

  int NPTC = I2-I1+3;

  V2C    = V1C + O2_00_ckd_mt_100_dv * (Numeric)(NPTC-1);  
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n"  
	   << "  CKD_MT 1.00 O2 v0<-v0 band:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.\n";
      return;
    }

  Numeric CO[(int)(NPTC+addF77fields)];


  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      CO[J]  = 0.000e0;
      Index I = I1+J;
      if ( (I > 0) && (I <= O2_00_ckd_mt_100_npt) ) 
	{
	  Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 
	  CO[J]  = O2_00_ckd_mt_100[I] / VJ;
	}
    }

  // ------------------- subroutine O2INF1 ----------------------------
  
  


  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmro2  = vmr[i];                                 // [1]
      Numeric vmrn2  = n2_abs[i];                              // [1]
      Numeric ADJWO2 = (vmro2 + 0.300e0*vmrn2) / 0.446e0 * 
                       (Pave/PO) * (Pave/PO) * (T1/Tave) * (T1/Tave);
      Numeric XKT    = Tave / 1.4387752e0;                     // = (T*k_B) / (h*c)    

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid. The abs. coeff. is then the 
      // cross section times the number density.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SO2 = 0.0e0;
	  if (CO[J] > 0.0e0)
	    {
	      SO2 = ADJWO2 * CO[J];
	    }

	  // CKD (cross section * number density) with radiative field
	  k[J] = SO2 * RADFN_FUN(VJ,XKT); // [1/cm]
	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > O2_00_ckd_mt_100_v1) && (V < O2_00_ckd_mt_100_v2) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

// =================================================================================

// CKD version MT 1.00 O2 v1<-v0 band absorption
// Model reference:
// CKD_MT 1.00 implementation of oxygen v1<-v0 band model of 
// Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,
// "Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands",
// Journal of Geophysical Research, vol 103, no. D4, pp. 3859-3863, 1998.
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            O2 v1<-v0 band according to CKD_MT 1.00  [1/m]
   \param    Cin            strength scaling factor                  [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> 
                            or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            O2 volume mixing ratio profile       [1]

   \remark   Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,<br>
             Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands,<br>
             J. Geophys. Res., 103, D4, 3859-3863, 1998. 

   \note     This absorption model is taken from the FORTRAN77 code of 
             CKD_MT version 1.00 written by<br>  
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421, USA<br> 
	     http://www.rtweb.aer.com/continuum_frame.html<br>
             <br> 
	     Oxygen band absorption model for the \f$a^1\Delta_g\f$ 
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
             \f$X^3\Sigma^-_g\f$ band system considering the 
             \f$\nu=0\f$
             \htmlonly&larr;\endhtmlonly \latexonly$\leftarrow$\endlatexonly 
             \f$\nu=1\f$
             transitions.
             
   \author Thomas Kuhn
   \date 2002-28-08
 */ 
void CKD_mt_v1v0_o2( MatrixView          xsec,
		     const Numeric       Cin,
		     const String&       model,
		     ConstVectorView     f_mono,
		     ConstVectorView     p_abs,
		     ConstVectorView     t_abs,
		     ConstVectorView     vmr )
{

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKDMT100"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKD_MT1.00 O2 band at 1.06 micrometer:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKDMT100\n";
      throw runtime_error(os.str());
    }


  // scaling factor of the O2 v1<-v0 band absorption
  Numeric  ScalingFac = 1.0000e0;
  if ( model == "user" )
    {
      ScalingFac = Cin; // input scaling factor of calculated absorption
    };

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );


  // ************************** CKD stuff ************************************

  const Numeric xLosmt    = 2.686763e19; // Loschmidt Number [molecules/cm^3]
  const Numeric T1        =  273.0e0;
  const Numeric TO        =  296.0e0;
  const Numeric PO        = 1013.0e0;
  const Numeric vmr_argon = 9.000e-3;    // VMR of argon is assumed to be const.


  // CKD_MT 1.00 implementation of oxygen v1<-v0 band model of
  // Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,
  // "Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands",
  // Journal of Geophysical Research, vol 103, no. D4, pp. 3859-3863, 1998.
  const Numeric V1S    = O2_10_ckd_mt_100_v1;
  const Numeric V2S    = O2_10_ckd_mt_100_v2;
  const Numeric DVS    = O2_10_ckd_mt_100_dv;
  const Numeric V1_osc =  9375.000e0;
  const Numeric HW1    =    58.960e0;
  const Numeric V2_osc =  9439.000e0;
  const Numeric HW2    =    45.040e0;
  const Numeric S1     =     1.166e-4;
  const Numeric S2     =     3.086e-5;


  // It is assumed here that f_mono is monotonically increasing with index!
  // In future change this return into a change of the loop over
  // the frequency f_mono. n_f_new < n_f
  Numeric V1ABS = f_mono[0]     / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  Numeric V2ABS = f_mono[n_f-1] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
  if ( (V1ABS < O2_10_ckd_mt_100_v1) || (V1ABS > O2_10_ckd_mt_100_v2) ||
       (V2ABS < O2_10_ckd_mt_100_v1) || (V2ABS > O2_10_ckd_mt_100_v2) )
    {
      out3  << "WARNING:\n"
            << "   CKD_MT 1.00 O2 v1<-v0 band:\n"
	    << "   input frequency vector exceeds range of model validity\n"
	    << "  " << O2_10_ckd_mt_100_v1 << "<->" << O2_10_ckd_mt_100_v2 << "cm^-1\n";
    }


  // ------------------- subroutine O2INF2 ----------------------------
  
  // retrieve the appropriate array sequence of the CKD model array.
  Numeric DVC = DVS;
  Numeric V1C = V1ABS - DVC;
  Numeric V2C = V2ABS + DVC;
  
  int NPTC = (int)( ((V2C-V1C)/DVC) + 3 );

  V2C = V1C + ( DVC * (Numeric)(NPTC-1) );
  
  if (NPTC < 1)
    {
      out3 << "WARNING:\n" 
	   << "  CKD_MT 1.00 O2 v1<-v0 band:\n"
	   << "  no elements of internal continuum coefficients could be found for the\n"
	   << "  input frequency range.\n"
	   << "  Leave the function without calculating the absorption.\n";
      return;
    }

  Numeric C[NPTC+addF77fields];
  C[0] = 0.000e0; // not used field of array

  for (Index J = 1 ; J <= NPTC ; ++J)
    {
      C[J]  = 0.000e0;
      Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 

      if ( (VJ > V1S) && (VJ < V2S) ) 
      	{
	  Numeric DV1 = VJ - V1_osc;
	  Numeric DV2 = VJ - V2_osc;

	  Numeric DAMP1 = 1.00e0;
	  Numeric DAMP2 = 1.00e0;

	  if ( DV1 < 0.00e0 ) 
	    {
	      DAMP1 = exp(DV1 / 176.100e0);
	    }

	  if ( DV2 < 0.00e0 ) 
	    {
	      DAMP2 = exp(DV2 / 176.100e0);
	    }
	  
	  Numeric O2INF = 0.31831e0 * 
                 ( ((S1 * DAMP1 / HW1)/(1.000e0 + pow((DV1/HW1),(Numeric)2.0e0) )) + 
		   ((S2 * DAMP2 / HW2)/(1.000e0 + pow((DV2/HW2),(Numeric)2.0e0) )) ) * 1.054e0;
	  C[J] = O2INF / VJ;
	}
    }


  // ------------------- subroutine O2INF2 ----------------------------
  
  
  // Loop pressure/temperature:
  for ( Index i = 0 ; i < n_p ; ++i )
    {
      Numeric Tave   = t_abs[i];                               // [K]
      Numeric Pave   = (p_abs[i]*1.000e-2);                    // [hPa]
      Numeric vmro2  = vmr[i];                                 // [1]
      Numeric WTOT   = 1.000e-20 * xLosmt * (Pave/PO) * (T1/Tave); // [molecules/cm^2]
      Numeric ADJWO2 = (vmro2 / 0.209e0) * WTOT * (Pave/PO) * (TO/Tave);
      Numeric XKT    = Tave / 1.4387752;                       // = (T*k_B) / (h*c)    

      // Molecular cross section calculated by CKD.
      // The cross sectionis calculated on the predefined 
      // CKD wavenumber grid.
      Numeric k[NPTC+addF77fields]; // [1/cm]
      k[0] = 0.00e0; // not used array field
      for (Index J = 1 ; J <= NPTC ; ++J)
	{
	  Numeric VJ  = V1C + (DVC * (Numeric)(J-1)); 
	  Numeric SO2 = 0.0e0;
	  if (C[J] > 0.0e0)
	    {
	      SO2 = ADJWO2 * C[J];
	    }

	  // CKD cross section without radiative field
	  k[J] = SO2 * RADFN_FUN(VJ,XKT); // [1]
	}
      

      // Loop input frequency array. The previously calculated cross section 
      // has therefore to be interpolated on the input frequencies.
      for ( Index s = 0 ; s < n_f ; ++s )
	{
	  // calculate the associated wave number (= 1/wavelength)
	  Numeric V = f_mono[s] / (SPEED_OF_LIGHT * 1.00e2); // [cm^-1]
	  if ( (V > V1S) && (V < V2S) )
	    {
	      // arts cross section [1/m]
	      // interpolate the k vector on the f_mono grid
	      xsec(s,i) +=  ScalingFac * 1.000e2 * XINT_FUN(V1C,V2C,DVC,k,V);
	    }
	}
    }
  
}

// #################################################################################

// CKD version 2.4 H2O continuum absorption model
/**

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            H2O continuum according to CKD2.4    [1/m]
   \param    isf            =0 self continuum, =1 foreign continuum
   \param    Cin            strength scaling factor              [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            H2O volume mixing ratio profile      [1]
   \param    n2_abs         N2 volume mixing ratio profile       [1]

   \note     this "crude" version of the CKD2.4 model is a f2c 
             conversion of the F77 code taken out of MonoRTM RT-model written by<br> 
             S. BOUKABARA, S.A. CLOUGH, and R. HOFFMAN<br>
             Atmospheric and Environmental Research Inc. (AER),<br> 
             Radiation and Climate Group<br>
             131 Hartwell Avenue<br>
             Lexington, MA 02421<br>
             USA<br>
             E-mail: sboukaba@aer.com, clough@aer.com

//   \remark   Reference: A. Borysow and L. Frommhold, 
//           The Astrophysical Journal, vol.311, pp.1043-1057, 1986
//           see <a href="http://adsabs.harvard.edu/article_service.html">for a scanned 
//           version of the paper</a>.

   \author Thomas Kuhn
   \date 2002-03-06
 */ 

void CKD24_H20( MatrixView          xsec,
		int                 isf,
		const Numeric       Cin,
		const String&       model,
		ConstVectorView     f_mono,
		ConstVectorView     p_abs,
		ConstVectorView     t_abs,
		ConstVectorView     vmr,
                ConstVectorView     n2_abs )
{
  //
  //
  // external function to call (original F77 code translated with f2c)
  /* INPUT PARAMETERS:                           */
  /*  P          [hPa]  TOTAL PRESSURE           */
  /*  T          [K]    TEMPERATURE              */
  /*  VMRH2O     [1]    H2O VOLUME MIXING RATIO  */
  /*  VMRN2      [1]    N2  VOLUME MIXING RATIO  */
  /*  VMRO2      [1]    O2  VOLUME MIXING RATIO  */
  /*  FREQ       [Hz]   FREQUENCY OF CALCULATION */
  extern double artsckd_(double p, double t, 
			 double vmrh2o, double vmrn2, double vmro2, 
			 double freq, int ivc);
  //
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  Numeric  XFAC  =  1.0000;             // scaling factor
  // ---------------------------------------------------------------------------------------
  

  // check the model name about consistency
  if ((model != "user") &&  (model != "CKD24"))
    {
      ostringstream os;
      os << "!!ERROR!!\n"  
	 << "CKDv2.4.2 H2O self/foreign continuum:\n"
	 << "INPUT model name is: " << model << ".\n"
	 << "VALID model names are user and CKD24\n";
      throw runtime_error(os.str());
    }


  // select the parameter set (!!model dominates values!!):
  if ( model == "CKD24" )
    {
      XFAC =  1.0000;
    }
  else if ( model == "user" )
    {
      XFAC =  Cin;
    }
  else
    {
      if (isf == 0) {
	ostringstream os;
	os << "H2O-SelfContCKD24: ERROR! Wrong model values given.\n"
	   << "allowed models are: 'CKD24', 'user'" << '\n';
	throw runtime_error(os.str());
      }
      if (isf == 1) {
	ostringstream os;
	os << "H2O-ForeignContCKD: ERROR! Wrong model values given.\n"
	   << "allowed models are: 'CKD24', 'user'" << '\n';
	throw runtime_error(os.str());
      }
    }
  
  if (isf == 0) {
    out3  << "H2O-SelfContCKD24: (model=" << model << ") parameter values in use:\n" 
	  << " XFAC = " << XFAC << "\n";
  }
  if (isf == 1) {
    out3  << "H2O-ForeignContCKD: (model=" << model << ") parameter values in use:\n" 
	  << " XFAC = " << XFAC << "\n";
  }


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  //  ivc = 1     : N2-N2 CKD version of Borysow-Fromhold model
  //  ivc = 21    : H2O CKD2.4  self cont part
  //  ivc = 22    : H2O CKD2.4  foreign cont part 
  //  ivc = 31    : MPMf87/s93  self cont part
  //  ivc = 32    : MPMf87/s93  foreign cont part 
  int ivc = 55;
  if (isf == 0) {
    ivc = 21; // CKD2.4  self continuum
    // ivc = 31; // MPMf87/s93  self continuum
  }
  if (isf == 1) {
    ivc = 22; // CKD2.4  foreign continuum
    //ivc = 32; // MPMf87/s93  foreign continuum
  }
  if ((ivc != 1) && (ivc != 21) && (ivc != 22) && (ivc != 31) && (ivc != 32)) {
    ostringstream os;
    os << "!!ERROR: CKD24 H2O model: wrong input parameter isf (=0,1) given!\n"
       << "retrun without calculation!" << "\n"
       << "actual value of isf is " << isf << "\n";
    throw runtime_error(os.str());
    return;
  }
  // ivc = 1;

  // Loop pressure/temperature:
  for ( Index i=0; i<n_p; ++i )
    {
      double T      = (double) t_abs[i];            // [K]
      double p      = (double) (p_abs[i]*1.000e-2); // [hPa]
      double vmrh2o = (double) vmr[i];              // [1]
      double vmrn2  = (double) n2_abs[i];           // [1]
      double vmro2  = 0.0e0;                        // [1]

      //cout << "------------------------------------------------\n";
      //cout << "CKD2.4 H2O: ivc   =" << ivc << "\n";
      //cout << "CKD2.4 H2O: T     =" << T << " K\n";
      //cout << "CKD2.4 H2O: p     =" << p << " hPa\n";
      //cout << "CKD2.4 H2O: vmrh2o=" << vmrh2o << "\n";
      //cout << "CKD2.4 H2O: vmrn2 =" << vmrn2 << "\n";
      //cout << "CKD2.4 H2O: vmro2 =" << vmro2 << "\n";
      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  // the second vmr of N2 will be multiplied at the stage of
	  // absorption calculation: abs =  vmr * xsec.
	  double f = (double) f_mono[s];            // [Hz]
	  if (ivc == 1) { // ---------- N2 -----------------
	    if (n2_abs[i] > 0.0e0) {
	      //cout << "CKD2.4 N2: f   =" << f << " Hz\n";
	      double cont = artsckd_(p, T, vmrh2o, vmrn2, vmro2, f, ivc);
	      xsec(s,i) +=  (Numeric) (cont / vmr[i]);
	      //cout << "CKD2.4 N2: abs =" << cont << " 1/m\n";
	    }
	  } else { // ---------------- H2O -----------------
	    if (vmr[i] > 0.0e0) {
	      //cout << "CKD2.4 H2O: f   =" << f << " Hz\n";
	      double cont = artsckd_(p, T, vmrh2o, vmrn2, vmro2, f, ivc);
	      xsec(s,i) +=  (Numeric) (cont / vmr[i]);
	      //cout << "CKD2.4 H2O: abs =" << cont << " 1/m\n";
	    }
	  }
	}
    }
  return;
}
//
// #################################################################################
//
/** 
   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    model          allows user defined input parameter set 
                            (Cin and xTin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid      [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid    [K] 
   \param    vmr            H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters Cin and xTin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'ATM', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: Pardo et al., IEEE, Trans. Ant. Prop., <br>
             Vol 49, No 12, pp. 1683-1694, 2001.

   \author Thomas Kuhn
   \date 2001-04-10
 */ 

void Pardo_ATM_H2O_ForeignContinuum( MatrixView          xsec,
				     const Numeric       Cin,
				     const String&       model,
				     ConstVectorView     f_mono,
				     ConstVectorView     p_abs,
				     ConstVectorView     t_abs,
				     ConstVectorView     vmr)
{
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Pardo et al. model (IEEE, Trans. Ant. Prop., 
  // Vol 49, No 12, pp. 1683-1694, 2001)
  const Numeric	C_ATM = 0.0315; // [1/m]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates parameters!!):
  Numeric C;
   if ( model == "ATM" )
     {
       C = C_ATM;
     }
   else if ( model == "user" )
     {
       C = Cin;
     }
   else
     {
       ostringstream os;
       os << "H2O-ForeignContATM01: ERROR! Wrong model values given.\n"
	  << "allowed models are: 'ATM', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out3  << "H2O-ForeignContATM01: (model=" << model << ") parameter values in use:\n" 
         << " C_f = " << C << "\n";

   const Index n_p = p_abs.nelem();	// Number of pressure levels
   const Index n_f = f_mono.nelem();	// Number of frequencies

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
      // since this is an effective "dry air" continuum, it is not really
      // it is not specifically attributed to N2, so we need the total 
      // dry air part in total which is equal to the total minus the 
      // water vapor pressure:
      Numeric  pd = p_abs[i] * ( 1.00000e0 - vmr[i] ); // [Pa]
      // since the H2O VMR will be multiplied in absCalc, we omit it here
      Numeric  pwdummy = p_abs[i]                    ; // [Pa]
      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
	{
	  // Becaue this is an effective "dry air" continuum, it is not really
	  // specific N2 but mainly caused by N2. Therefore the N2 vmr must be 
	  // canceled out here which is later in absCalc multiplied 
	  // (calculation: abs = vmr * xsec):
	  xsec(s,i) += C *                  // strength [1/(m*Hz²Pa²)] 
	    pow( (f_mono[s]/(Numeric)2.25e11), (Numeric)2. ) * // quadratic f dependence [1]
	    pow( ((Numeric)300.0/t_abs[i]), (Numeric)3. ) * // free T dependence      [1]
	    (pd/1.01300e5)                * // p_dry dependence       [1]
	    (pwdummy/1.01300e5);            // p_H2O dependence       [1]
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
			  ConstVectorView     vmr	 )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 H2O continuum model 
  // (AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric	MPM93fo_orig =  1780.000e9;  // [Hz]
  const Numeric	MPM93b1_orig = 22300.000;    // [Hz/Pa]
  const Numeric	MPM93b2_orig =     0.952;    // [1]
  const Numeric	MPM93b3_orig =    17.600e4;  // [Hz/Pa]
  const Numeric	MPM93b4_orig =    30.500;    // [1]
  const Numeric	MPM93b5_orig =     2.000;    // [1]
  const Numeric	MPM93b6_orig =     5.000;    // [1]
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
   out3  << "H2O-ContMPM93: (model=" << model << ") parameter values in use:\n" 
         << " fo = " << MPM93fopcl << "\n"
         << " b1 = " << MPM93b1pcl << "\n"
         << " b2 = " << MPM93b2pcl << "\n"
         << " b3 = " << MPM93b3pcl << "\n"
         << " b4 = " << MPM93b4pcl << "\n"
         << " b5 = " << MPM93b5pcl << "\n"
         << " b6 = " << MPM93b6pcl << "\n";

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric strength =  MPM93b1pcl * p_abs[i] * pow( th, (Numeric)3.5 )
        * exp(MPM93b2pcl * (1 - th));
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
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
                            O2 according to MPM89 [1/m]
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
             Allowed models: 'MPM85', 'MPM85Lines', 'MPM85Continuum', 'MPM85NoCoupling', 
             'MPM85NoCutoff', and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe,<br>
             <i>An updated model for millimeter wave propagation in moist air,</i>,<br>
	     Radio Science, vol. 20, pp. 1069-1089, 1985
   \author Thomas Kuhn
   \date 2002-04-05
 */ 

void MPM85O2AbsModel( MatrixView          xsec,
		      const Numeric	  CCin,       // continuum scale factor 
		      const Numeric	  CLin,       // line strength scale factor
		      const Numeric	  CWin,       // line broadening scale factor
		      const Numeric	  COin,       // line coupling scale factor
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     h2o_abs,
		      ConstVectorView     vmr )
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //      0             1           2         3          4      5         6
  //      f0           a1          a2        a3         a4     a5        a6
  //    [GHz]      [kHz/hPa]      [1]     [MHz/hPa]    [1]   [1/kPa]     [1]
  const Numeric mpm85[48][7] = { 
   { 49.452379 ,       0.12 ,    11.830 ,    8.40 ,    0.0 ,  5.600 ,   1.700 },
   { 49.962257 ,       0.34 ,    10.720 ,    8.50 ,    0.0 ,  5.600 ,   1.700 },
   { 50.474238 ,       0.94 ,     9.690 ,    8.60 ,    0.0 ,  5.600 ,   1.700 },
   { 50.987748 ,       2.46 ,     8.690 ,    8.70 ,    0.0 ,  5.500 ,   1.700 },
   { 51.503350 ,       6.08 ,     7.740 ,    8.90 ,    0.0 ,  5.600 ,   1.800 },
   { 52.021409 ,      14.14 ,     6.840 ,    9.20 ,    0.0 ,  5.500 ,   1.800 },
   { 52.542393 ,      31.02 ,     6.000 ,    9.40 ,    0.0 ,  5.700 ,   1.800 },
   { 53.066906 ,      64.10 ,     5.220 ,    9.70 ,    0.0 ,  5.300 ,   1.900 },
   { 53.595748 ,     124.70 ,     4.480 ,   10.00 ,    0.0 ,  5.400 ,   1.800 },
   { 54.129999 ,     228.00 ,     3.810 ,   10.20 ,    0.0 ,  4.800 ,   2.000 },
   { 54.671157 ,     391.80 ,     3.190 ,   10.50 ,    0.0 ,  4.800 ,   1.900 },
   { 55.221365 ,     631.60 ,     2.620 ,   10.79 ,    0.0 ,  4.170 ,   2.100 },
   { 55.783800 ,     953.50 ,     2.115 ,   11.10 ,    0.0 ,  3.750 ,   2.100 },
   { 56.264777 ,     548.90 ,     0.010 ,   16.46 ,    0.0 ,  7.740 ,   0.900 },
   { 56.363387 ,    1344.00 ,     1.655 ,   11.44 ,    0.0 ,  2.970 ,   2.300 },
   { 56.968180 ,    1763.00 ,     1.255 ,   11.81 ,    0.0 ,  2.120 ,   2.500 },
   { 57.612481 ,    2141.00 ,     0.910 ,   12.21 ,    0.0 ,  0.940 ,   3.700 },
   { 58.323874 ,    2386.00 ,     0.621 ,   12.66 ,    0.0 , -0.550 ,  -3.100 },
   { 58.446589 ,    1457.00 ,     0.079 ,   14.49 ,    0.0 ,  5.970 ,   0.800 },
   { 59.164204 ,    2404.00 ,     0.386 ,   13.19 ,    0.0 , -2.440 ,   0.100 },
   { 59.590982 ,    2112.00 ,     0.207 ,   13.60 ,    0.0 ,  3.440 ,   0.500 },
   { 60.306057 ,    2124.00 ,     0.207 ,   13.82 ,    0.0 , -4.130 ,   0.700 },
   { 60.434775 ,    2461.00 ,     0.386 ,   12.97 ,    0.0 ,  1.320 ,  -1.000 },
   { 61.150558 ,    2504.00 ,     0.621 ,   12.48 ,    0.0 , -0.360 ,   5.800 },
   { 61.800152 ,    2298.00 ,     0.910 ,   12.07 ,    0.0 , -1.590 ,   2.900 },
   { 62.411212 ,    1933.00 ,     1.255 ,   11.71 ,    0.0 , -2.660 ,   2.300 },
   { 62.486253 ,    1517.00 ,     0.078 ,   14.68 ,    0.0 , -4.770 ,   0.900 },
   { 62.997974 ,    1503.00 ,     1.660 ,   11.39 ,    0.0 , -3.340 ,   2.200 },
   { 63.568515 ,    1087.00 ,     2.110 ,   11.08 ,    0.0 , -4.170 ,   2.000 },
   { 64.127764 ,     733.50 ,     2.620 ,   10.78 ,    0.0 , -4.480 ,   2.000 },
   { 64.678900 ,     463.50 ,     3.190 ,   10.50 ,    0.0 , -5.100 ,   1.800 },
   { 65.224067 ,     274.80 ,     3.810 ,   10.20 ,    0.0 , -5.100 ,   1.900 },
   { 65.764769 ,     153.00 ,     4.480 ,   10.00 ,    0.0 , -5.700 ,   1.800 },
   { 66.302088 ,      80.09 ,     5.220 ,    9.70 ,    0.0 , -5.500 ,   1.800 },
   { 66.836827 ,      39.46 ,     6.000 ,    9.40 ,    0.0 , -5.900 ,   1.700 },
   { 67.369595 ,      18.32 ,     6.840 ,    9.20 ,    0.0 , -5.600 ,   1.800 },
   { 67.900862 ,       8.01 ,     7.740 ,    8.90 ,    0.0 , -5.800 ,   1.700 },
   { 68.431001 ,       3.30 ,     8.690 ,    8.70 ,    0.0 , -5.700 ,   1.700 },
   { 68.960306 ,       1.28 ,     9.690 ,    8.60 ,    0.0 , -5.600 ,   1.700 },
   { 69.489021 ,       0.47 ,    10.720 ,    8.50 ,    0.0 , -5.600 ,   1.700 },
   { 70.017342 ,       0.16 ,    11.830 ,    8.40 ,    0.0 , -5.600 ,   1.700 },
  { 118.750341 ,     945.00 ,     0.000 ,   15.92 ,    0.0 , -0.440 ,   0.900 },
  { 368.498350 ,      67.90 ,     0.020 ,   19.20 ,    0.6 ,  0.000 ,   0.000 },
  { 424.763120 ,     638.00 ,     0.011 ,   19.16 ,    0.6 ,  0.000 ,   0.000 },
  { 487.249370 ,     235.00 ,     0.011 ,   19.20 ,    0.6 ,  0.000 ,   0.000 },
  { 715.393150 ,      99.60 ,     0.089 ,   18.10 ,    0.6 ,  0.000 ,   0.000 },
  { 773.838730 ,     671.00 ,     0.079 ,   18.10 ,    0.6 ,  0.000 ,   0.000 },
  { 834.145330 ,     180.00 ,     0.079 ,   18.10 ,    0.6 ,  0.000 ,   0.000 },
};

  // number of lines of Liebe O2-line catalog (0-47 lines)
  const Index i_first = 0;
  const Index i_last  = 47; // all the spec. lines up to 1THz
  // const Index i_last  = 40; // only the 60GHz complex + 118GHz line
  

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM85 model (Liebe, Radio Science, 20, 1069-1089, 1985):
  const Numeric CC_MPM85 = 1.00000;
  const Numeric CL_MPM85 = 1.00000;
  const Numeric CW_MPM85 = 1.00000;
  const Numeric CO_MPM85 = 1.00000;
  int   AppCutoff = 0;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  if ( model == "MPM85" )
    {
      CC = CC_MPM85;
      CL = CL_MPM85;
      CW = CW_MPM85;
      CO = CO_MPM85;
    }
  else if ( model == "MPM85Lines" )
    {
      CC = 0.000;
      CL = CL_MPM85;
      CW = CW_MPM85;
      CO = CO_MPM85;
    }
  else if ( model == "MPM85Continuum" )
    {
      CC = CC_MPM85;
      CL = 0.000;
      CW = 0.000;
      CO = 0.000;
    }
  else if ( model == "MPM85NoCoupling" )
    {
      CC = CC_MPM85;
      CL = CL_MPM85;
      CW = CW_MPM85;
      CO = 0.000;
    }
  else if ( model == "MPM85NoCutoff" )
    {
      CC = CC_MPM85;
      CL = CL_MPM85;
      CW = CW_MPM85;
      CO = CO_MPM85;
      AppCutoff = 1;
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
      os << "O2-MPM85: ERROR! Wrong model values given.\n"
	 << "Valid models are: 'MPM85' 'MPM85Lines' 'MPM85Continuum' 'MPM85NoCoupling' 'MPM85NoCutoff'" 
         << "and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "O2-MPM85: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n"
	<< " CO = " << CO << "\n";
  

  // O2 continuum parameters of MPM92:
  const Numeric	S0 =  6.140e-4; // line strength                        [ppm]
  const Numeric G0 =  5.600e-3;  // line width                           [GHz/kPa]
  const Numeric	X0 =  0.800;    // temperature dependence of line width [1]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;

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
	  os << "ERROR: MPM87 O2 full absorption model has detected a O2 volume mixing ratio of " 
	     << vmr[i] << " which is below the threshold of " << VMRCalcLimit << ".\n"
	     << "Therefore no calculation is performed.\n";
	  throw runtime_error(os.str());
	  return;
	}

      // relative inverse temperature [1]
      Numeric theta     = (300.0 / t_abs[i]);
      // H2O partial pressure [kPa]
      Numeric pwv       = Pa_to_kPa * p_abs[i] * h2o_abs[i];
      // dry air partial pressure [kPa]
      Numeric pda       = (Pa_to_kPa * p_abs[i]) - pwv;
      // here the total pressure is devided by the O2 vmr for the 
      // P_dry calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pda_dummy = pda;
      // O2 continuum strength [ppm]
      Numeric strength_cont =  S0 * pda_dummy * pow( theta, (Numeric)2. );
      // O2 continuum pseudo line broadening [GHz]
      Numeric gam_cont      =  G0 * ( pda + 1.10*pwv ) *  pow( theta, X0 ); // GHz
      
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
	{
	  // input frequency in [GHz]
	  Numeric ff = f_mono[s] * Hz_to_GHz; 
	  // O2 continuum absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  // here the rolloff parameter FAC is implemented!
	  // Numeric FAC =  1.000 / ( pow( ff, 2) + pow( 60.000, 2) );
	  // if we let the non-proofen rollofff away as in further version: 
	  Numeric FAC =  1.000 ;
	  Numeric Nppc =  CC * strength_cont * FAC * ff * gam_cont /
	                  ( pow( ff, (Numeric)2.)
                            + pow( gam_cont, (Numeric)2.) );
	  
	  // Loop over MPM85 O2 spectral lines:
	  Numeric Nppl  = 0.0;
	  for ( Index l = i_first; l <= i_last; ++l )
	    {
	      // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
	      Numeric strength = CL * mpm85[l][1] * 1.000e-6  * pda_dummy * 
		                      pow(theta, (Numeric)3.) * exp(mpm85[l][2]*(1.000-theta)) /
		                      mpm85[l][0];
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * ( mpm85[l][3] * 1.000e-3 * 
	                              ( (       pda * pow(theta, ((Numeric)0.80-mpm85[l][4]))) + 
                                        (1.10 * pwv * theta) ) );
	      // line mixing parameter [1]
	      Numeric delta    = CO * mpm85[l][5] * 1.000e-3 * 
			              pda * pow(theta, mpm85[l][6]);
	      // absorption [dB/km] like in the original MPM92
	      Nppl            += strength * MPMLineShapeO2Function(gam, mpm85[l][0], ff, delta); 
	    }
	  // in MPM85 there is a cutoff for O2 line absorption if abs_l < 0 
	  // absorption cannot be less than 0 according to MPM87 philosophy.
	  // since this cutoff is only 'detectable' in the source code and not in the
          // publications we assume this cutoff also for MPM85 since it is also 
          // implemented in MPM87. 
	  if (AppCutoff == 0) 
	    {
	      if (Nppl < 0.000)  Nppl = 0.0000;
	    }
	  //
	  // O2 line absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * (Nppl+Nppc) / VMRISO;
	}
    }
  return;
}
//
// #################################################################################
// 
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
                            O2 according to MPM89 [1/m]
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
             Allowed models: 'MPM87', 'MPM87Lines', 'MPM87Continuum', 'MPM87NoCoupling', 
             'MPM87NoCutoff', and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and D. H. Layton,<br>
             <i>Millimeter-wave properties of the atmosphere: 
                Laboratory studies and propagation modelling</i>,<br>
	     U.S. Dept. of Commerce, National Telecommunications and Information
	     Administration, Institute for Communication Sciences, rep. 87-224,<br>
             325 Broadway, Boulder, CO 80303-3328

   \author Thomas Kuhn
   \date 2002-04-05
 */ 

void MPM87O2AbsModel( MatrixView          xsec,
		      const Numeric	  CCin,       // continuum scale factor 
		      const Numeric	  CLin,       // line strength scale factor
		      const Numeric	  CWin,       // line broadening scale factor
		      const Numeric	  COin,       // line coupling scale factor
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     h2o_abs,
		      ConstVectorView     vmr )
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0            1           2        3          4        5        6
  //         f0           a1          a2       a3         a4      a5        a6
  //        [GHz]      [kHz/hPa]     [1]    [MHz/hPa]    [1]         [1/kPa]
  const Numeric mpm87[48][7] = { 
    { 49.452379 ,       0.12 ,    11.830 ,    8.40 ,    0.0 ,  6.600 ,   1.700}, // 0
    { 49.962257 ,       0.34 ,    10.720 ,    8.50 ,    0.0 ,  6.600 ,   1.700}, // 1
    { 50.474238 ,       0.94 ,     9.690 ,    8.60 ,    0.0 ,  6.600 ,   1.700}, // 2 
    { 50.987748 ,       2.46 ,     8.690 ,    8.70 ,    0.0 ,  6.500 ,   1.700}, // 3
    { 51.503350 ,       6.08 ,     7.740 ,    8.90 ,    0.0 ,  6.627 ,   1.800}, // 4
    { 52.021409 ,      14.14 ,     6.840 ,    9.20 ,    0.0 ,  6.347 ,   1.800}, // 5
    { 52.542393 ,      31.02 ,     6.000 ,    9.40 ,    0.0 ,  6.046 ,   1.800}, // 6
    { 53.066906 ,      64.10 ,     5.220 ,    9.70 ,    0.0 ,  5.719 ,   1.900}, // 7
    { 53.595748 ,     124.70 ,     4.480 ,   10.00 ,    0.0 ,  5.400 ,   1.800}, // 8
    { 54.129999 ,     228.00 ,     3.810 ,   10.20 ,    0.0 ,  5.157 ,   2.000}, // 9
    { 54.671157 ,     391.80 ,     3.190 ,   10.50 ,    0.0 ,  4.783 ,   1.900}, // 10
    { 55.221365 ,     631.60 ,     2.620 ,   10.79 ,    0.0 ,  4.339 ,   2.100}, // 11
    { 55.783800 ,     953.50 ,     2.115 ,   11.10 ,    0.0 ,  4.011 ,   2.100}, // 12
    { 56.264777 ,     548.90 ,     0.010 ,   16.46 ,    0.0 ,  2.772 ,   0.900}, // 13
    { 56.363387 ,    1344.00 ,     1.655 ,   11.44 ,    0.0 ,  3.922 ,   2.300}, // 14
    { 56.968180 ,    1763.00 ,     1.255 ,   11.81 ,    0.0 ,  3.398 ,   2.500}, // 15
    { 57.612481 ,    2141.00 ,     0.910 ,   12.21 ,    0.0 ,  1.145 ,   3.200}, // 16
    { 58.323874 ,    2386.00 ,     0.621 ,   12.66 ,    0.0 , -0.317 ,  -2.500}, // 17
    { 58.446589 ,    1457.00 ,     0.079 ,   14.49 ,    0.0 ,  6.270 ,   0.800}, // 18
    { 59.164204 ,    2404.00 ,     0.386 ,   13.19 ,    0.0 , -4.119 ,   0.100}, // 19
    { 59.590982 ,    2112.00 ,     0.207 ,   13.60 ,    0.0 ,  6.766 ,   0.500}, // 20
    { 60.306057 ,    2124.00 ,     0.207 ,   13.82 ,    0.0 , -6.183 ,   0.700}, // 21
    { 60.434775 ,    2461.00 ,     0.386 ,   12.97 ,    0.0 ,  3.290 ,  -0.400}, // 22
    { 61.150558 ,    2504.00 ,     0.621 ,   12.48 ,    0.0 , -1.591 ,   3.500}, // 23
    { 61.800152 ,    2298.00 ,     0.910 ,   12.07 ,    0.0 , -2.068 ,   2.900}, // 24
    { 62.411212 ,    1933.00 ,     1.255 ,   11.71 ,    0.0 , -4.158 ,   2.300}, // 25
    { 62.486253 ,    1517.00 ,     0.078 ,   14.68 ,    0.0 , -4.068 ,   0.900}, // 26
    { 62.997974 ,    1503.00 ,     1.660 ,   11.39 ,    0.0 , -4.482 ,   2.200}, // 27
    { 63.568515 ,    1087.00 ,     2.110 ,   11.08 ,    0.0 , -4.442 ,   2.000}, // 28
    { 64.127764 ,     733.50 ,     2.620 ,   10.78 ,    0.0 , -4.687 ,   2.000}, // 29
    { 64.678900 ,     463.50 ,     3.190 ,   10.50 ,    0.0 , -5.074 ,   1.800}, // 30
    { 65.224067 ,     274.80 ,     3.810 ,   10.20 ,    0.0 , -5.403 ,   1.900}, // 31
    { 65.764769 ,     153.00 ,     4.480 ,   10.00 ,    0.0 , -5.610 ,   1.800}, // 32
    { 66.302088 ,      80.09 ,     5.220 ,    9.70 ,    0.0 , -5.896 ,   1.800}, // 33
    { 66.836827 ,      39.46 ,     6.000 ,    9.40 ,    0.0 , -6.194 ,   1.700}, // 34
    { 67.369595 ,      18.32 ,     6.840 ,    9.20 ,    0.0 , -6.468 ,   1.800}, // 35
    { 67.900862 ,       8.01 ,     7.740 ,    8.90 ,    0.0 , -6.718 ,   1.700}, // 36
    { 68.431001 ,       3.30 ,     8.690 ,    8.70 ,    0.0 , -6.700 ,   1.700}, // 37
    { 68.960306 ,       1.28 ,     9.690 ,    8.60 ,    0.0 , -6.600 ,   1.700}, // 38
    { 69.489021 ,       0.47 ,    10.720 ,    8.50 ,    0.0 , -6.600 ,   1.700}, // 39
    { 70.017342 ,       0.16 ,    11.830 ,    8.40 ,    0.0 , -6.600 ,   1.700}, // 40
   { 118.750341 ,     945.00 ,     0.000 ,   16.30 ,    0.0 , -0.134 ,   0.800}, // 41
  {  368.498350 ,      67.90 ,     0.020 ,   19.20 ,    0.6 ,  0.000 ,   0.000}, // 42
  {  424.763120 ,     638.00 ,     0.011 ,   19.16 ,    0.6 ,  0.000 ,   0.000}, // 43
  {  487.249370 ,     235.00 ,     0.011 ,   19.20 ,    0.6 ,  0.000 ,   0.000}, // 44
  {  715.393150 ,      99.60 ,     0.089 ,   18.10 ,    0.6 ,  0.000 ,   0.000}, // 45
  {  773.838730 ,     671.00 ,     0.079 ,   18.10 ,    0.6 ,  0.000 ,   0.000}, // 46
  {  834.145330 ,     180.00 ,     0.079 ,   18.10 ,    0.6 ,  0.000 ,   0.000}  // 47
  };

  // number of lines of Liebe O2-line catalog (0-47 lines)
  const Index i_first = 0;
  const Index i_last  = 47; // all the spec. lines up to 1THz
  // const Index i_last  = 40; // only the 60GHz complex + 118GHz line
  

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM87 model (NITA Report 87-224):
  const Numeric CC_MPM87 = 1.00000;
  const Numeric CL_MPM87 = 1.00000;
  const Numeric CW_MPM87 = 1.00000;
  const Numeric CO_MPM87 = 1.00000;
  int   AppCutoff = 0;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  if ( model == "MPM87" )
    {
      CC = CC_MPM87;
      CL = CL_MPM87;
      CW = CW_MPM87;
      CO = CO_MPM87;
    }
  else if ( model == "MPM87Lines" )
    {
      CC = 0.000;
      CL = CL_MPM87;
      CW = CW_MPM87;
      CO = CO_MPM87;
    }
  else if ( model == "MPM87Continuum" )
    {
      CC = CC_MPM87;
      CL = 0.000;
      CW = 0.000;
      CO = 0.000;
    }
  else if ( model == "MPM87NoCoupling" )
    {
      CC = CC_MPM87;
      CL = CL_MPM87;
      CW = CW_MPM87;
      CO = 0.000;
    }
  else if ( model == "MPM87NoCutoff" )
    {
      // !!ATTENTION!!
      // In the window regions the total absorption can get negative values.
      // So be carefull with this selection!
      CC = CC_MPM87;
      CL = CL_MPM87;
      CW = CW_MPM87;
      CO = CO_MPM87;
      AppCutoff = 1;
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
      os << "O2-MPM87: ERROR! Wrong model values given.\n"
	 << "Valid models are: 'MPM87' 'MPM87Lines' 'MPM87Continuum' 'MPM87NoCoupling' 'MPM87NoCutoff'" 
         << "and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "O2-MPM87: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n"
	<< " CO = " << CO << "\n";
  

  // O2 continuum parameters of MPM92:
  const Numeric	S0 =  6.140e-4; // line strength                        [ppm]
  const Numeric G0 =  4.800e-3;  // line width [GHz/kPa] !! 14% lower than in all the other versions !!
  const Numeric	X0 =  0.800;    // temperature dependence of line width [1]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;

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
	  os << "ERROR: MPM87 O2 full absorption model has detected a O2 volume mixing ratio of " 
	     << vmr[i] << " which is below the threshold of " << VMRCalcLimit << ".\n"
	     << "Therefore no calculation is performed.\n";
	  throw runtime_error(os.str());
	  return;
	}

      // relative inverse temperature [1]
      Numeric theta     = (300.0 / t_abs[i]);
      // H2O partial pressure [kPa]
      Numeric pwv       = Pa_to_kPa * p_abs[i] * h2o_abs[i];
      // dry air partial pressure [kPa]
      Numeric pda       = (Pa_to_kPa * p_abs[i]) - pwv;
      // here the total pressure is devided by the O2 vmr for the 
      // P_dry calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pda_dummy = pda;
      // O2 continuum strength [ppm]
      Numeric strength_cont =  S0 * pda_dummy * pow( theta, (Numeric)2. );
      // O2 continuum pseudo line broadening [GHz]
      Numeric gam_cont      =  G0 * ( pda + 1.10*pwv ) *  pow( theta, X0 ); // GHz
      
      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
	{
	  // input frequency in [GHz]
	  Numeric ff = f_mono[s] * Hz_to_GHz; 
	  // O2 continuum absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  Numeric Nppc =  CC * strength_cont * ff * gam_cont /
	                  ( pow( ff, (Numeric)2.) + pow( gam_cont, (Numeric)2.) );
	  
	  // Loop over MPM87 O2 spectral lines:
	  Numeric Nppl  = 0.0;
	  for ( Index l = i_first; l <= i_last; ++l )
	    {
	      // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
	      Numeric strength = CL * mpm87[l][1] * 1.000e-6  * pda_dummy * 
		                      pow(theta, (Numeric)3.) * exp(mpm87[l][2]*(1.000-theta)) /
		                      mpm87[l][0];
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * ( mpm87[l][3] * 1.000e-3 * 
	                              ( (       pda * pow(theta, ((Numeric)0.80-mpm87[l][4]))) + 
                                        (1.10 * pwv * theta) ) );
	      // line mixing parameter [1]
	      Numeric delta    = CO * mpm87[l][5] * 1.000e-3 * 
			              pda * pow(theta, mpm87[l][6]);
	      // absorption [dB/km] like in the original MPM92
	      Nppl            += strength * MPMLineShapeO2Function(gam, mpm87[l][0], ff, delta); 
	    }
	  // in MPM87 there is a cutoff for O2 line absorption if abs_l < 0 
	  // absorption cannot be less than 0 according to MPM87 source code.
	  if (AppCutoff == 0) 
	    {
	      if (Nppl < 0.000)  Nppl = 0.0000;
	    }
	  //
	  // O2 line absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * (Nppl+Nppc) / VMRISO;
	}
    }
  return;
}
//
// #################################################################################
// 
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
                            O2 according to MPM89 [1/m]
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
             Allowed models: 'MPM89', 'MPM89Lines', 'MPM89Continuum', 'MPM89NoCoupling', 
             'MPM89NoCutoff', and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe,<br>
             <i>MPM - an atmospheric millimeter-wave propagation model</i>,<br>
             Int. J. Infrared and Mill. Waves, Vol 10, pp. 631-650, 1989.

   \author Thomas Kuhn
   \date 2002-04-05
 */ 

void MPM89O2AbsModel( MatrixView          xsec,
		      const Numeric	  CCin,       // continuum scale factor 
		      const Numeric	  CLin,       // line strength scale factor
		      const Numeric	  CWin,       // line broadening scale factor
		      const Numeric	  COin,       // line coupling scale factor
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     h2o_abs,
		      ConstVectorView     vmr )
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0            1           2        3          4        5        6
  //         f0           a1          a2       a3         a4      a5        a6
  //        [GHz]      [kHz/hPa]     [1]    [MHz/hPa]    [1]         [1/kPa]
  const Numeric mpm89[44][7] = { 
    {   50.474238,        0.94 ,    9.694 ,    8.60 ,    0.0 ,  1.600 ,   5.520 }, // 0
    {   50.987749,        2.46 ,    8.694 ,    8.70 ,    0.0 ,  1.400 ,   5.520 }, // 1
    {   51.503350,        6.08 ,    7.744 ,    8.90 ,    0.0 ,  1.165 ,   5.520 }, // 2
    {   52.021410,       14.14 ,    6.844 ,    9.20 ,    0.0 ,  0.883 ,   5.520 }, // 3
    {   52.542394,       31.02 ,    6.004 ,    9.40 ,    0.0 ,  0.579 ,   5.520 }, // 4
    {   53.066907,       64.10 ,    5.224 ,    9.70 ,    0.0 ,  0.252 ,   5.520 }, // 5
    {   53.595749,      124.70 ,    4.484 ,   10.00 ,    0.0 , -0.066 ,   5.520 }, // 6
    {   54.130000,      228.00 ,    3.814 ,   10.20 ,    0.0 , -0.314 ,   5.520 }, // 7
    {   54.671159,      391.80 ,    3.194 ,   10.50 ,    0.0 , -0.706 ,   5.520 }, // 8
    {   55.221367,      631.60 ,    2.624 ,   10.79 ,    0.0 , -1.151 ,   5.514 }, // 9
    {   55.783802,      953.50 ,    2.119 ,   11.10 ,    0.0 , -0.920 ,   5.025 }, // 10
    {   56.264775,      548.90 ,    0.015 ,   16.46 ,    0.0 ,  2.881 ,  -0.069 }, // 11
    {   56.363389,     1344.00 ,    1.660 ,   11.44 ,    0.0 , -0.596 ,   4.750 }, // 12
    {   56.968206,     1763.00 ,    1.260 ,   11.81 ,    0.0 , -0.556 ,   4.104 }, // 13
    {   57.612484,     2141.00 ,    0.915 ,   12.21 ,    0.0 , -2.414 ,   3.536 }, // 14
    {   58.323877,     2386.00 ,    0.626 ,   12.66 ,    0.0 , -2.635 ,   2.686 }, // 15
    {   58.446590,     1457.00 ,    0.084 ,   14.49 ,    0.0 ,  6.848 ,  -0.647 }, // 16
    {   59.164207,     2404.00 ,    0.391 ,   13.19 ,    0.0 , -6.032 ,   1.858 }, // 17
    {   59.590983,     2112.00 ,    0.212 ,   13.60 ,    0.0 ,  8.266 ,  -1.413 }, // 18
    {   60.306061,     2124.00 ,    0.212 ,   13.82 ,    0.0 , -7.170 ,   0.916 }, // 19
    {   60.434776,     2461.00 ,    0.391 ,   12.97 ,    0.0 ,  5.664 ,  -2.323 }, // 20
    {   61.150560,     2504.00 ,    0.626 ,   12.48 ,    0.0 ,  1.731 ,  -3.039 }, // 21
    {   61.800154,     2298.00 ,    0.915 ,   12.07 ,    0.0 ,  1.738 ,  -3.797 }, // 22
    {   62.411215,     1933.00 ,    1.260 ,   11.71 ,    0.0 , -0.048 ,  -4.277 }, // 23
    {   62.486260,     1517.00 ,    0.083 ,   14.68 ,    0.0 , -4.290 ,   0.238 }, // 24
    {   62.997977,     1503.00 ,    1.665 ,   11.39 ,    0.0 ,  0.134 ,  -4.860 }, // 25
    {   63.568518,     1087.00 ,    2.115 ,   11.08 ,    0.0 ,  0.541 ,  -5.079 }, // 26
    {   64.127767,      733.50 ,    2.620 ,   10.78 ,    0.0 ,  0.814 ,  -5.525 }, // 27
    {   64.678903,      463.50 ,    3.195 ,   10.50 ,    0.0 ,  0.415 ,  -5.520 }, // 28
    {   65.224071,      274.80 ,    3.815 ,   10.20 ,    0.0 ,  0.069 ,  -5.520 }, // 29
    {   65.764772,      153.00 ,    4.485 ,   10.00 ,    0.0 , -0.143 ,  -5.520 }, // 30
    {   66.302091,       80.09 ,    5.225 ,    9.70 ,    0.0 , -0.428 ,  -5.520 }, // 31
    {   66.836830,       39.46 ,    6.005 ,    9.40 ,    0.0 , -0.726 ,  -5.520 }, // 32
    {   67.369598,       18.32 ,    6.845 ,    9.20 ,    0.0 , -1.002 ,  -5.520 }, // 33
    {   67.900867,        8.01 ,    7.745 ,    8.90 ,    0.0 , -1.255 ,  -5.520 }, // 34
    {   68.431005,        3.30 ,    8.695 ,    8.70 ,    0.0 , -1.500 ,  -5.520 }, // 35
    {   68.960311,        1.28 ,    9.695 ,    8.60 ,    0.0 , -1.700 ,  -5.520 }, // 36
    {  118.750343,      945.00 ,    0.009 ,   16.30 ,    0.0 , -0.247 ,   0.003 }, // 37
    {  368.498350,       67.90 ,    0.049 ,   19.20 ,    0.6 ,  0.000 ,   0.000 }, // 38
    {  424.763124,      638.00 ,    0.044 ,   19.16 ,    0.6 ,  0.000 ,   0.000 }, // 39
    {  487.249370,      235.00 ,    0.049 ,   19.20 ,    0.6 ,  0.000 ,   0.000 }, // 40
    {  715.393150,       99.60 ,    0.145 ,   18.10 ,    0.6 ,  0.000 ,   0.000 }, // 41
    {  773.839675,      671.00 ,    0.130 ,   18.10 ,    0.6 ,  0.000 ,   0.000 }, // 42
    {  834.145330,      180.00 ,    0.147 ,   18.10 ,    0.6 ,  0.000 ,   0.000 }  // 43
  };

  // number of lines of Liebe O2-line catalog (0-43 lines)
  const Index i_first = 0;
  const Index i_last  = 43; // all the spec. lines up to 1THz
  // const Index i_last  = 37; // only the 60GHz complex + 118GHz line
  

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM89 model (IJIMW, Vol 10, pp. 631-650, 1989):
  const Numeric CC_MPM89 = 1.00000;
  const Numeric CL_MPM89 = 1.00000;
  const Numeric CW_MPM89 = 1.00000;
  const Numeric CO_MPM89 = 1.00000;
  int   AppCutoff = 0; 
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  if ( model == "MPM89" )
    {
      CC = CC_MPM89;
      CL = CL_MPM89;
      CW = CW_MPM89;
      CO = CO_MPM89;
    }
  else if ( model == "MPM89Lines" )
    {
      CC = 0.000;
      CL = CL_MPM89;
      CW = CW_MPM89;
      CO = CO_MPM89;
    }
  else if ( model == "MPM89Continuum" )
    {
      CC = CC_MPM89;
      CL = 0.000;
      CW = 0.000;
      CO = 0.000;
    }
  else if ( model == "MPM89NoCoupling" )
    {
      CC = CC_MPM89;
      CL = CL_MPM89;
      CW = CW_MPM89;
      CO = 0.000;
    }
  else if ( model == "MPM89NoCutoff" )
    {
      CC = CC_MPM89;
      CL = CL_MPM89;
      CW = CW_MPM89;
      CO = CO_MPM89;
      AppCutoff = 1; 
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
      os << "O2-MPM89: ERROR! Wrong model values given.\n"
	 << "Valid models are: 'MPM89' 'MPM89Lines' 'MPM89Continuum' 'MPM89NoCoupling' 'MPM89NoCutoff'" 
         << "and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "O2-MPM89: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n"
	<< " CO = " << CO << "\n";
  

  // O2 continuum parameters of MPM92:
  const Numeric	S0 =  6.140e-4; // line strength                        [ppm]
  const Numeric G0 =  5.60e-3;  // line width                           [GHz/kPa]
  const Numeric	X0 =  0.800;    // temperature dependence of line width [1]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;

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
	  os << "ERROR: MPM89 O2 full absorption model has detected a O2 volume mixing ratio of " 
	     << vmr[i] << " which is below the threshold of " << VMRCalcLimit << ".\n"
	     << "Therefore no calculation is performed.\n";
	  throw runtime_error(os.str());
	  return;
	}

      // relative inverse temperature [1]
      Numeric theta     = (300.0 / t_abs[i]);
      // H2O partial pressure [kPa]
      Numeric pwv       = Pa_to_kPa * p_abs[i] * h2o_abs[i];
      // dry air partial pressure [kPa]
      Numeric pda       = (Pa_to_kPa * p_abs[i]) - pwv;
      // here the total pressure is devided by the O2 vmr for the 
      // P_dry calculation because we calculate xsec and not abs: abs = vmr * xsec
      Numeric pda_dummy = pda;
      // O2 continuum strength [ppm]
      Numeric strength_cont =  S0 * pda_dummy * pow( theta, (Numeric)2. );
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
	                  ( pow( ff, (Numeric)2.) + pow( gam_cont, (Numeric)2.) );
	  
	  // Loop over MPM89 O2 spectral lines:
	  Numeric Nppl  = 0.0;
	  for ( Index l = i_first; l <= i_last; ++l )
	    {
	      // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
	      Numeric strength = CL * mpm89[l][1] * 1.000e-6  * pda_dummy * 
		                      pow(theta, (Numeric)3.) * exp(mpm89[l][2]*(1.000-theta)) /
		                      mpm89[l][0];
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * ( mpm89[l][3] * 1.000e-3 * 
	                              ( (       pda * pow(theta, ((Numeric)0.80-mpm89[l][4]))) + 
                                        (1.10 * pwv * theta) ) );
	      // line mixing parameter [1]
	      Numeric delta    = CO * ( (mpm89[l][5] + mpm89[l][6] * theta) * 1.000e-3 * 
			                pda * pow(theta, (Numeric)0.8) );
	      // absorption [dB/km] like in the original MPM92
	      Nppl            += strength * MPMLineShapeO2Function(gam, mpm89[l][0], ff, delta); 
	    }
	  // in MPM89 we adopt the cutoff for O2 line absorption if abs_l < 0 
	  // absorption cannot be less than 0 according to MPM87 source code.
	  if (AppCutoff == 0) 
	    {
	      if (Nppl < 0.000)  Nppl = 0.0000;
	    }
	  //
	  // O2 line absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * (Nppl+Nppc) / VMRISO;
	}
    }
  return;
}
//
// #################################################################################
// 
//
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
             Allowed models: 'MPM92', 'MPM92Lines', 'MPM92Continuum', 'MPM92NoCoupling', 
             'MPM92NoCutoff', and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe, P. W. Rosenkranz and G. A. Hufford,<br>
             <i>Atmospheric 60-GHz Oxygen Spectrum: New Laboratory Measurements 
             and Line Parameters</i>,<br>
             JQSRT, Vol 48, pp. 629-643, 1992

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM92O2AbsModel( MatrixView          xsec,
		      const Numeric	  CCin,       // continuum scale factor 
		      const Numeric	  CLin,       // line strength scale factor
		      const Numeric	  CWin,       // line broadening scale factor
		      const Numeric	  COin,       // line coupling scale factor
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     h2o_abs,
		      ConstVectorView     vmr )
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0           1            2       3        4      5      6
  //         f0          a1           a2      a3       a4     a5     a6
  //        [GHz]     [kHz/hPa]      [1]   [MHz/hPa]  [1]    [10³/hPa]
  const Numeric mpm92[44][7] = { 
    {   50.474238,       0.094,      9.694,    0.850,     0.0,   0.210,    0.685}, // 0
    {   50.987749,       0.246,      8.694,    0.870,     0.0,   0.190,    0.680}, // 1
    {   51.503350,       0.608,      7.744,    0.890,     0.0,   0.171,    0.673}, // 2 
    {   52.021410,       1.414,      6.844,    0.920,     0.0,   0.144,    0.664}, // 3
    {   52.542394,       3.102,      6.004,    0.940,     0.0,   0.118,    0.653}, // 4
    {   53.066907,       6.410,      5.224,    0.970,     0.0,   0.114,    0.621}, // 5
    {   53.595749,      12.470,      4.484,    1.000,     0.0,   0.200,    0.508}, // 6
    {   54.130000,      22.800,      3.814,    1.020,     0.0,   0.291,    0.375}, // 7
    {   54.671159,      39.180,      3.194,    1.050,     0.0,   0.325,    0.265}, // 8
    {   55.221367,      63.160,      2.624,    1.080,     0.0,   0.224,    0.295}, // 9
    {   55.783802,      95.350,      2.119,    1.110,     0.0,  -0.144,    0.613}, // 0
    {   56.264775,      54.890,      0.015,    1.646,     0.0,   0.339,   -0.098}, // 11
    {   56.363389,     134.400,      1.660,    1.144,     0.0,  -0.258,    0.655}, // 12
    {   56.968206,     176.300,      1.260,    1.181,     0.0,  -0.362,    0.645}, // 13
    {   57.612484,     214.100,      0.915,    1.221,     0.0,  -0.533,    0.606}, // 14
    {   58.323877,     238.600,      0.626,    1.266,     0.0,  -0.178,    0.044}, // 15
    {   58.446590,     145.700,      0.084,    1.449,     0.0,   0.650,   -0.127}, // 16
    {   59.164207,     240.400,      0.391,    1.319,     0.0,  -0.628,    0.231}, // 17
    {   59.590983,     211.200,      0.212,    1.360,     0.0,   0.665,   -0.078}, // 18
    {   60.306061,     212.400,      0.212,    1.382,     0.0,  -0.613,    0.070}, // 19
    {   60.434776,     246.100,      0.391,    1.297,     0.0,   0.606,   -0.282}, // 20
    {   61.150560,     250.400,      0.626,    1.248,     0.0,   0.090,   -0.058}, // 21
    {   61.800154,     229.800,      0.915,    1.207,     0.0,   0.496,   -0.662}, // 22
    {   62.411215,     193.300,      1.260,    1.171,     0.0,   0.313,   -0.676}, // 23
    {   62.486260,     151.700,      0.083,    1.468,     0.0,  -0.433,    0.084}, // 24
    {   62.997977,     150.300,      1.665,    1.139,     0.0,   0.208,   -0.668}, // 25
    {   63.568518,     108.700,      2.115,    1.110,     0.0,   0.094,   -0.614}, // 26
    {   64.127767,      73.350,      2.620,    1.080,     0.0,  -0.270,   -0.289}, // 27
    {   64.678903,      46.350,      3.195,    1.050,     0.0,  -0.366,   -0.259}, // 28
    {   65.224071,      27.480,      3.815,    1.020,     0.0,  -0.326,   -0.368}, // 29
    {   65.764772,      15.300,      4.485,    1.000,     0.0,  -0.232,   -0.500}, // 30
    {   66.302091,       8.009,      5.225,    0.970,     0.0,  -0.146,   -0.609}, // 31
    {   66.836830,       3.946,      6.005,    0.940,     0.0,  -0.147,   -0.639}, // 32
    {   67.369598,       1.832,      6.845,    0.920,     0.0,  -0.174,   -0.647}, // 33
    {   67.900867,       0.801,      7.745,    0.890,     0.0,  -0.198,   -0.655}, // 34
    {   68.431005,       0.330,      8.695,    0.870,     0.0,  -0.210,   -0.660}, // 35
    {   68.960311,       0.128,      9.695,    0.850,     0.0,  -0.220,   -0.665}, // 36
    {  118.750343,      94.500,      0.009,    1.630,     0.0,  -0.031,    0.008}, // 37
    {  368.498350,       6.790,      0.049,    1.920,     0.6,   0.000,    0.000}, // 38
    {  424.763124,      63.800,      0.044,    1.926,     0.6,   0.000,    0.000}, // 39
    {  487.249370,      23.500,      0.049,    1.920,     0.6,   0.000,    0.000}, // 40
    {  715.393150,       9.960,      0.145,    1.810,     0.6,   0.000,    0.000}, // 41
    {  773.839675,      67.100,      0.130,    1.810,     0.6,   0.000,    0.000}, // 42
    {  834.145330,      18.000,      0.147,    1.810,     0.6,   0.000,    0.000}}; // 43

  // number of lines of Liebe O2-line catalog (0-43 lines)
  const Index i_first = 0;
  const Index i_last  = 43; // all the spec. lines up to 1THz
  // const Index i_last  = 37; // only the 60GHz complex + 118GHz line
  

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM92 model (JQSRT, Vol 48, pp. 629-643, 1992):
  const Numeric CC_MPM92 = 1.00000;
  const Numeric CL_MPM92 = 1.00000;
  const Numeric CW_MPM92 = 1.00000;
  const Numeric CO_MPM92 = 1.00000;
  int AppCutoff = 0;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CC, CL, CW, CO;
  if ( model == "MPM92" )
    {
      CC = CC_MPM92;
      CL = CL_MPM92;
      CW = CW_MPM92;
      CO = CO_MPM92;
    }
  else if ( model == "MPM92Lines" )
    {
      CC = 0.000;
      CL = CL_MPM92;
      CW = CW_MPM92;
      CO = CO_MPM92;
    }
  else if ( model == "MPM92Continuum" )
    {
      CC = CC_MPM92;
      CL = 0.000;
      CW = 0.000;
      CO = 0.000;
    }
  else if ( model == "MPM92NoCoupling" )
    {
      CC = CC_MPM92;
      CL = CL_MPM92;
      CW = CW_MPM92;
      CO = 0.000;
    }
  else if ( model == "MPM92NoCutoff" )
    {
      CC = CC_MPM92;
      CL = CL_MPM92;
      CW = CW_MPM92;
      CO = CO_MPM92;
      AppCutoff = 1;
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
      os << "O2-MPM92: ERROR! Wrong model values given.\n"
	 << "Valid models are: 'MPM92' 'MPM92Lines' 'MPM92Continuum' 'MPM92NoCoupling' 'MPM92NoCutoff'" 
         << "and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "O2-MPM92: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n"
	<< " CO = " << CO << "\n";
  

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;

  // O2 continuum parameters of MPM92:
  const Numeric	S0 =  6.140e-5; // line strength                        [ppm]
  const Numeric G0 =  0.560e-3; // line width                           [GHz/hPa]
  const Numeric	X0 =  0.800;    // temperature dependence of line width [1]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
	  os << "ERROR: MPM92 O2 full absorption model has detected a O2 volume mixing ratio of " 
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
      Numeric pda_dummy = pda;
      // O2 continuum strength [ppm]
      Numeric strength_cont =  S0 * pda_dummy * pow( theta, (Numeric)2. );
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
	                  ( pow( ff, (Numeric)2.) + pow( gam_cont, (Numeric)2.) );
	  
	  // Loop over MPM92 O2 spectral lines:
	  Numeric Nppl  = 0.0;
	  for ( Index l = i_first; l <= i_last; ++l )
	    {
	      // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
	      Numeric strength = CL * 1.000e-6  * pda_dummy * mpm92[l][1] / mpm92[l][0] * 
		                      pow(theta, (Numeric)3.) * exp(mpm92[l][2]*(1.0-theta));
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * ( mpm92[l][3] * 0.001 * 
	                              ( (       pda * pow(theta, ((Numeric)0.8-mpm92[l][4]))) + 
                                        (1.10 * pwv * theta) ) );
	      // line mixing parameter [1]
	      //		  if (l < 11) CD = 1.1000;
	      Numeric delta    = CO * ( (mpm92[l][5] + mpm92[l][6] * theta) * 
			              (pda+pwv) * 0.001 * pow(theta, (Numeric)0.8) );
	      // absorption [dB/km] like in the original MPM92
	      Nppl            += strength * MPMLineShapeO2Function(gam, mpm92[l][0], ff, delta); 
	    }
	  // in MPM92 we adopt the cutoff for O2 line absorption if abs_l < 0 
	  // absorption cannot be less than 0 according to MPM87 and MPM93 source code.
	  if (AppCutoff == 0) 
	    {
	      if (Nppl < 0.000)  Nppl = 0.0000;
	    }
	  //
	  // O2 line absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * (Nppl+Nppc) / VMRISO;
	}
    }
  return;
}
//
// #################################################################################
// 
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
             'MPM93NoCutoff', and 'user'. See the user guide for detailed explanations.

   \remark   Reference: H. J. Liebe and G. A. Hufford and M. G. Cotton,<br>
             <i>Propagation modeling of moist air and suspended water/ice
             particles at frequencies below 1000 GHz</i>,<br>
             AGARD 52nd Specialists Meeting of the Electromagnetic Wave
             Propagation Panel,<br> Palma de Mallorca, Spain, 1993, May 17-21

   \author Thomas Kuhn
   \date 2001-11-05
 */ 

void MPM93O2AbsModel( MatrixView          xsec,
		      const Numeric	  CCin,       // continuum scale factor 
		      const Numeric	  CLin,       // line strength scale factor
		      const Numeric	  CWin,       // line broadening scale factor
		      const Numeric	  COin,       // line coupling scale factor
		      const String&       model,
		      ConstVectorView     f_mono,
		      ConstVectorView     p_abs,
		      ConstVectorView     t_abs,
		      ConstVectorView     h2o_abs,    // VMR 0f H2O
		      ConstVectorView     vmr )       // VMR of O2
{
  //
  // Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
  //         0             1           2         3         4      5        6
  //         f0           a1           a2       a3        a4      a5       a6
  //        [GHz]      [kHz/hPa]      [1]    [MHz/hPa]    [1]       [10³/hPa]
  const Numeric mpm93[44][7] = { 
    {   50.474238,       0.094,      9.694,    0.890,     0.0,   0.240,    0.790}, // 0
    {   50.987749,       0.246,      8.694,    0.910,     0.0,   0.220,    0.780}, // 1
    {   51.503350,       0.608,      7.744,    0.940,     0.0,   0.197,    0.774}, // 2 
    {   52.021410,       1.414,      6.844,    0.970,     0.0,   0.166,    0.764}, // 3
    {   52.542394,       3.102,      6.004,    0.990,     0.0,   0.136,    0.751}, // 4
    {   53.066907,       6.410,      5.224,    1.020,     0.0,   0.131,    0.714}, // 5
    {   53.595749,      12.470,      4.484,    1.050,     0.0,   0.230,    0.584}, // 6
    {   54.130000,      22.800,      3.814,    1.070,     0.0,   0.335,    0.431}, // 7
    {   54.671159,      39.180,      3.194,    1.100,     0.0,   0.374,    0.305}, // 8
    {   55.221367,      63.160,      2.624,    1.130,     0.0,   0.258,    0.339}, // 9
    {   55.783802,      95.350,      2.119,    1.170,     0.0,  -0.166,    0.705}, // 10
    {   56.264775,      54.890,      0.015,    1.730,     0.0,   0.390,   -0.113}, // 11
    {   56.363389,     134.400,      1.660,    1.200,     0.0,  -0.297,    0.753}, // 12
    {   56.968206,     176.300,      1.260,    1.240,     0.0,  -0.416,    0.742}, // 13
    {   57.612484,     214.100,      0.915,    1.280,     0.0,  -0.613,    0.697}, // 14
    {   58.323877,     238.600,      0.626,    1.330,     0.0,  -0.205,    0.051}, // 15
    {   58.446590,     145.700,      0.084,    1.520,     0.0,   0.748,   -0.146}, // 16
    {   59.164207,     240.400,      0.391,    1.390,     0.0,  -0.722,    0.266}, // 17
    {   59.590983,     211.200,      0.212,    1.430,     0.0,   0.765,   -0.090}, // 18
    {   60.306061,     212.400,      0.212,    1.450,     0.0,  -0.705,    0.081}, // 19
    {   60.434776,     246.100,      0.391,    1.360,     0.0,   0.697,   -0.324}, // 20
    {   61.150560,     250.400,      0.626,    1.310,     0.0,   0.104,   -0.067}, // 21
    {   61.800154,     229.800,      0.915,    1.270,     0.0,   0.570,   -0.761}, // 22
    {   62.411215,     193.300,      1.260,    1.230,     0.0,   0.360,   -0.777}, // 23
    {   62.486260,     151.700,      0.083,    1.540,     0.0,  -0.498,    0.097}, // 24
    {   62.997977,     150.300,      1.665,    1.200,     0.0,   0.239,   -0.768}, // 25
    {   63.568518,     108.700,      2.115,    1.170,     0.0,   0.108,   -0.706}, // 26
    {   64.127767,      73.350,      2.620,    1.130,     0.0,  -0.311,   -0.332}, // 27
    {   64.678903,      46.350,      3.195,    1.100,     0.0,  -0.421,   -0.298}, // 28
    {   65.224071,      27.480,      3.815,    1.070,     0.0,  -0.375,   -0.423}, // 29
    {   65.764772,      15.300,      4.485,    1.050,     0.0,  -0.267,   -0.575}, // 30
    {   66.302091,       8.009,      5.225,    1.020,     0.0,  -0.168,   -0.700}, // 31
    {   66.836830,       3.946,      6.005,    0.990,     0.0,  -0.169,   -0.735}, // 32
    {   67.369598,       1.832,      6.845,    0.970,     0.0,  -0.200,   -0.744}, // 33
    {   67.900867,       0.801,      7.745,    0.940,     0.0,  -0.228,   -0.753}, // 34
    {   68.431005,       0.330,      8.695,    0.920,     0.0,  -0.240,   -0.760}, // 35
    {   68.960311,       0.128,      9.695,    0.900,     0.0,  -0.250,   -0.765}, // 36
    {  118.750343,      94.500,      0.009,    1.630,     0.0,  -0.036,    0.009}, // 37
    {  368.498350,       6.790,      0.049,    1.920,     0.6,   0.000,    0.000}, // 38
    {  424.763124,      63.800,      0.044,    1.930,     0.6,   0.000,    0.000}, // 39
    {  487.249370,      23.500,      0.049,    1.920,     0.6,   0.000,    0.000}, // 40
    {  715.393150,       9.960,      0.145,    1.810,     0.6,   0.000,    0.000}, // 41
    {  773.839675,      67.100,      0.130,    1.820,     0.6,   0.000,    0.000}, // 42
    {  834.145330,      18.000,      0.147,    1.810 ,    0.6,   0.000,    0.000}}; // 43
  // number of lines of Liebe O2-line catalog (0-43 lines)
  const Index i_first = 0;
  const Index i_last  = 43; // all the spec. lines up to 1THz
  // const Index i_last  = 37; // only the 60GHz complex + 118GHz line
  

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM87 model (Radio Science, 20(5), 1985, 1069):
  const Numeric CC_MPM93 = 1.00000;
  const Numeric CL_MPM93 = 1.00000;
  const Numeric CW_MPM93 = 1.00000;
  const Numeric CO_MPM93 = 1.00000;
  int   AppCutoff = 0;
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
  else if ( model == "MPM93NoCutoff" )
    {
      // !!ATTENTION!!
      // In the window regions the total absorption can get negative values.
      // So be carefull with this selection!
      CC = CC_MPM93;
      CL = CL_MPM93;
      CW = CW_MPM93;
      CO = CO_MPM93;
      AppCutoff = 1;
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
	 << "Valid models are: 'MPM93' 'MPM93Lines' 'MPM93Continuum' 'MPM93NoCoupling' 'MPM93NoCutoff'" 
         << "and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "O2-MPM93: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n"
	<< " CO = " << CO << "\n";
  

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;

  // O2 continuum parameters of MPM93:
  const Numeric	S0 =  6.140e-5; // line strength                        [ppm]
  const Numeric G0 =  0.560e-3; // line width                           [GHz/hPa]
  const Numeric	X0 =  0.800;    // temperature dependence of line width [1]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      // old version without VMRISO: Numeric pda_dummy = pda / vmr[i];
      Numeric pda_dummy = pda;
      // O2 continuum strength [ppm]
      Numeric strength_cont =  S0 * pda_dummy * pow( theta, (Numeric)2. );
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
	                  ( pow( ff, (Numeric)2.)
                            + pow( gam_cont, (Numeric)2.) );
	  
	  // Loop over MPM93 O2 spectral lines:
	  Numeric Nppl  = 0.0;
	  for ( Index l = i_first; l <= i_last; ++l )
	    {
	      // line strength [ppm]   S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
	      Numeric strength = CL * 1.000e-6  * pda_dummy * 
		                      mpm93[l][1] / mpm93[l][0] * 
		                      pow(theta, (Numeric)3.) * exp(mpm93[l][2]*(1.0-theta));
	      // line broadening parameter [GHz]
	      Numeric gam      = CW * ( mpm93[l][3] * 0.001 * 
                                      ( (       pda * pow(theta, ((Numeric)0.8-mpm93[l][4]))) + 
                                        (1.10 * pwv * theta) ) );
	      // line mixing parameter [1]
	      //		  if (l < 11) CD = 1.1000;
	      Numeric delta    = CO * ( (mpm93[l][5] + mpm93[l][6] * theta) * 
			              (pda+pwv) * pow(theta, (Numeric)0.8)
                                        * (Numeric)0.001 );
	      // absorption [dB/km] like in the original MPM93
	      Nppl            += strength * MPMLineShapeO2Function(gam, mpm93[l][0], ff, delta); 
	    }
	  // in MPM93 there is a cutoff for O2 line absorption if abs_l < 0 
	  // absorption cannot be less than 0 according to MPM93 philosophy.
	  if (AppCutoff == 0) 
	    {
	      if (Nppl < 0.000)  Nppl = 0.0000;// <---!!IMPORTANT FEATURE!!
	    }
	  //
	  // O2 line absorption [1/m]
	  // cross section: xsec = absorption / var
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  xsec(s,i) += dB_km_to_1_m * 0.1820 * ff * (Nppl+Nppc) / VMRISO;
	}
    }
  return;
}
//
// #################################################################################
// 
//   Oxygen complex at 60 GHz plus mm O2 lines plus O2 continuum
//
//    REFERENCES FOR EQUATIONS AND COEFFICIENTS:
//    P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
//     BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
//    H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
//    M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
//    SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
//    This version differs from Liebe's MPM92 in two significant respects:
//    1. It uses the modification of the 1- line width temperature dependence
//    recommended by Schwartz: (1/T).
//    2. It uses the same temperature dependence (X) for submillimeter 
//    line widths as in the 60 GHz band: (1/T)**0.8 
//
//   history:
//   05-01-95  P. Rosenkranz 
//   11-05-97  P. Rosenkranz - 1- line modification.
//   12-16-98  pwr - updated submm freq's and intensities from HITRAN96
//
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
                            O2 according to the P. W. Rosenkranz, 1993 [1/m]
   \param    CCin           O2-continuum scale factor  [1]
   \param    CLin           O2 line strength scale factor [1]
   \param    CWin           O2 line broadening scale factor [1]
   \param    COin           O2 line coupling scale factor [1]
   \param    model          allows user defined input parameter set 
                            (CCin, CLin, CWin, and COin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    version        determines model version: 1988, 1993, 1998    			    
   \param    f_mono         predefined frequency grid        [Hz]
   \param    p_abs          predefined pressure              [Pa]
   \param    t_abs          predefined temperature grid      [K] 
   \param    vmrh2o         H2O volume mixing ratio profile  [1]
   \param    vmr            O2 volume mixing ratio profile   [1]

   \note     Except for  model 'user' the input parameters CCin, CLin, CWin, and COin 
             are neglected (model dominates over parameters).<br>
             Allowed models:<br> 
             'Rosenkranz', 'RosenkranzLines', 'RosenkranzContinuum',
             'RosenkranzNoCoupling', and 'user'. <br> 
	     For the parameter  version the following three string values are allowed:
	     'PWR88', 'PWR93', 'PWR98'.<br> 
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
		      const String&     model,     // model selection string
		      const String&     version,   // PWR98, PWR93 or PWR88
		      ConstVectorView   f_mono,
		      ConstVectorView   p_abs,
		      ConstVectorView   t_abs,
		      ConstVectorView   vmrh2o,
                      ConstVectorView   vmr )
{
  const Index n_lines = 40; // all O2 lines in this model (range: 50-850 GHz)
  //
  // lines are arranged 1-,1+,3-,3+,etc. in spin-rotation spectrum
  // line center frequency in [GHz]
  const Numeric F93[n_lines] = { 118.7503,  56.2648,  62.4863,  58.4466,   // 00-03
			          60.3061,  59.5910,  59.1642,  60.4348,   // 04-07
	  		          58.3239,  61.1506,  57.6125,  61.8002,   // 08-11
			          56.9682,  62.4112,  56.3634,  62.9980,   // 12-15
			          55.7838,  63.5685,  55.2214,  64.1278,   // 16-19
			          54.6712,  64.6789,  54.1300,  65.2241,   // 20-23
			          53.5957,  65.7648,  53.0669,  66.3021,   // 24-27
			          52.5424,  66.8368,  52.0214,  67.3696,   // 28-31
			          51.5034,  67.9009, 368.4984, 424.7631,   // 32-35
			         487.2494, 715.3932, 773.8397, 834.1453};  // 36-39

  // intensities in the submm range are updated according to HITRAN96
  const Numeric F98[n_lines] = { 118.7503,  56.2648,  62.4863,  58.4466,  60.3061,  59.5910,
				  59.1642,  60.4348,  58.3239,  61.1506,  57.6125,  61.8002,
				  56.9682,  62.4112,  56.3634,  62.9980,  55.7838,  63.5685,
			          55.2214,  64.1278,  54.6712,  64.6789,  54.1300,  65.2241,
			          53.5957,  65.7648,  53.0669,  66.3021,  52.5424,  66.8368,
			          52.0214,  67.3696,  51.5034,  67.9009, 368.4984, 424.7632,
				 487.2494, 715.3931, 773.8397, 834.1458};


  // line strength at T=300K in [cm² * Hz] 
  const Numeric S93[n_lines] = { 0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14,
		                 0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14,
		        	 0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,
		                 0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14,
		                 0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14,
		                 0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,
		                 0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15,
		                 0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16,
		                 0.8191E-17, 0.1423E-16, 0.6460E-15, 0.7047E-14, 
		                 0.3011E-14, 0.1826E-14, 0.1152E-13, 0.3971E-14};

  // intensities in the submm range are updated according to HITRAN96
  const Numeric S98[n_lines] = { 0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14,
				 0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14,
			         0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,
			         0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14,
			         0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14,
			         0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,
			         0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15,
			         0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16,
			         0.8191E-17, 0.1423E-16, 0.6494E-15, 0.7083E-14, 
                                 0.3025E-14, 0.1835E-14, 0.1158E-13, 0.3993E-14};

  // temperature exponent of the line strength in [1]
  const Numeric BE[n_lines] = { 0.009,   0.015,   0.083,   0.084, 
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
  const Numeric WB300 = 0.56; // [MHz/mbar]=[MHz/hPa]
  const Numeric X     = 0.80; // [1]

  // line width parameter [GHz/bar]
  const Numeric W300[n_lines] = {   1.630, 1.646, 1.468, 1.449, 
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
  const Numeric Y93[n_lines] = { -0.0233,  0.2408, -0.3486,  0.5227,
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
  const Numeric Y88[n_lines] = { -0.0244,  0.2772, -0.4068,  0.6270,
                                 -0.6183,  0.6766, -0.4119,  0.3290, 
                                  0.0317, -0.1591,  0.1145, -0.2068, 
                                  0.3398, -0.4158,  0.3922, -0.4482,
                                  0.4011, -0.4442,  0.4339, -0.4687,
                                  0.4783, -0.5074,  0.5157, -0.5403,
                                  0.5400, -0.5610,  0.5719, -0.5896,
                                  0.6046, -0.6194,  0.6347, -0.6468, 
                                  0.6627, -0.6718,  0.0000,  0.0000,
                                  0.0000,  0.0000,  0.0000,  0.0000};

  // v parameter for the calculation of Y [1/bar]
  const Numeric V[n_lines] ={    0.0079, -0.0978,  0.0844, -0.1273,
				 0.0699, -0.0776,  0.2309, -0.2825, 
				 0.0436, -0.0584,  0.6056, -0.6619, 
				 0.6451, -0.6759,  0.6547, -0.6675,
				 0.6135, -0.6139,  0.2952, -0.2895, 
				 0.2654, -0.2590,  0.3750, -0.3680, 
				 0.5085, -0.5002,  0.6206, -0.6091,
				 0.6526, -0.6393,  0.6640, -0.6475,
				 0.6729, -0.6545,  0.0000,  0.0000,
				 0.0000,  0.0000,  0.0000,  0.0000};
  // range of lines to take into account for the line absorption part
  const Index first_line = 0;  // first line for calculation
  const Index last_line  = 39; // last line for calculation

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
  Numeric CC, CL, CW, CO, Y300[n_lines], S300[n_lines], F[n_lines];
  int oldnewflag = 0;

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
  out3  << "O2-PWR93: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CL = " << CL << "\n"
	<< " CW = " << CW << "\n"
	<< " CO = " << CO << "\n";
  

  // determin if version Rosenkranz 1993 or Rosenkranz 1988 is selected
  if ( (version != "PWR98") && (version != "PWR93") && (version != "PWR88") )
    {
      ostringstream os;
      os << "O2-PWR93/PWR88: ERROR! Wrong version is selected.\n"
	 << "Valid versions are:\n" 
	 << "  'PWR98'  updates of F and S to HISTRAN96 and M.J.Schwartz, MIT, 1997\n"
         << "           suggestions implemented.\n"
         << "  'PWR93'  for the oxygen absorption model described in  \n" 
         << "           P. W. Rosenkranz, Chapter 2, in M. A. Janssen,\n"
         << "           Atmospheric Remote Sensing by Microwave Radiometry,\n"
         << "           John Wiley & Sons, Inc., 1993.\n"
         << "  'PWR88'  for the oxygen absorption model described in \n" 
         << "           P. W. Rosenkranz, Interference coefficients for the \n" 
	 << "           overlapping oxygen lines in air, \n"
         << "           JQSRT, 1988, Volume 39, 287-297.\n";
      throw runtime_error(os.str());
    }


  // select version dependent parameters
  if ( version == "PWR88" ) {
    for ( Index i=0; i<n_lines; ++i ) 
      {
	F[i]    = F93[i];
	S300[i] = S93[i];
        Y300[i] = Y88[i];
      };
  }
  if ( version == "PWR93" ) {
    for ( Index i=0; i<n_lines; ++i )
      {
	F[i]    = F93[i];
        S300[i] = S93[i];
        Y300[i] = Y93[i];
      };
  }
  if ( version == "PWR98" ) {
    for ( Index i=0; i<n_lines; ++i )
      {
	F[i]    = F98[i];
	S300[i] = S98[i];
        Y300[i] = Y93[i];
      };
  }

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric TH     = 3.0000e2 / t_abs[i];
      Numeric TH1    = (TH-1.000e0);
      Numeric B      = pow(TH, X);
      // partial pressure of H2O and dry air [hPa]
      Numeric PRESWV = Pa_to_hPa * (p_abs[i] * vmrh2o[i]);
      Numeric PRESDA = Pa_to_hPa * (p_abs[i] * (1.000e0 - vmrh2o[i]));
      Numeric DEN    = 0.001*(PRESDA*B + 1.1*PRESWV*TH); // [hPa]
      Numeric DENS   = 0.001*(PRESDA   + 1.1*PRESWV)*TH; // [hPa]
      Numeric DFNR   = WB300*DEN; // [GHz]
      
      // continuum absorption [1/m/GHz]
      Numeric CCONT  = CC * 1.23e-10 * pow( TH, (Numeric)2. ) * p_abs[i];

      // Loop over input frequency
      for ( Index s=0; s<n_f; ++s )
	{
	  // initial O2 line absorption at position ff 
	  // Numeric O2ABS  = 0.000e0;cd safff

	  // input frequency in [GHz]
	  Numeric ff   = Hz_to_GHz * f_mono[s]; 

	  // continuum absorption [Neper/km]
	  Numeric CONT = CCONT * (ff * ff * DFNR / (ff*ff + DFNR*DFNR));

	  // Loop over Rosnekranz '93 spectral line frequency:
	  Numeric SUM  = 0.000e0;
	  for ( Index l=first_line; l<=last_line; ++l )
	    {
	      Numeric DF = CW * W300[l] * DEN; // [hPa]
	      // 118 line update according to M. J. Schwartz, MIT, 1997
	      if ( (version == "PWR98") && (fabs((F[l]-118.75)) < 0.10) )
		{
		  DF = CW * W300[l] * DENS; // [hPa]
		}
	      Numeric Y    = CO * 0.001 * 0.01 * p_abs[i] * B * ( Y300[l] + V[l]*TH1 );
	      Numeric STR  = CL * S300[l] * exp(-BE[l] * TH1);
	      Numeric SF1  = ( DF + (ff-F[l])*Y ) / ( (ff-F[l])*(ff-F[l]) + DF*DF );
	      Numeric SF2  = ( DF - (ff+F[l])*Y ) / ( (ff+F[l])*(ff+F[l]) + DF*DF );
	      SUM         += STR * (SF1+SF2) * (ff/F[l]) * (ff/F[l]);
	    }
	  // O2 absorption [Neper/km]
	  // Rosenkranz uses the factor 0.5034e12 in the calculation of the abs coeff.
	  // This factor is the product of several terms:
          // 0.5034e12 = ISORATIO *   VMR   * (Hz/GHz) * (k_B*300K)^-1 
          //           = 0.995262 * 0.20946 *   10^-9  * 2.414322e21(hPa*cm^2*km)^-1
	  //             |---- 0.2085 ----|   |---- 2.414322e12(hPa*cm^2*km)^-1 ---|
	  //             |---- 0.2085 ----|   |---- 2.414322e10( Pa*cm^2*km)^-1 ---|
	  // O2ABS = 2.4143e12 * SUM * PRESDA * pow(TH, 3.0) / PI;
	  // O2ABS = CONT + (2.414322e10 * SUM * p_abs[i] * pow(TH, 3.0) / PI);
	  // unit conversion x Nepers/km = y 1/m  --->  y = x * 1.000e-3 
	  // therefore 2.414322e10 --> 2.414322e7
	  // xsec [1/m] 
	  xsec(s,i) += CONT + (2.414322e7 * SUM * p_abs[i] * pow(TH, (Numeric)3.) / PI);
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

   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            O2-continuum according to MPM93 [1/m]
   \param    S0in           O2-continuum strength [1/Pa]
   \param    G0in           O2-continuum width [Hz/Pa]
   \param    XS0in          O2-continuum strength temperature exponent [1]
   \param    XG0in          O2-continuum width temperature exponent    [1]
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
			 ConstVectorView     vmr	 )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  // const Numeric	S0_MPM93  =  6.140e-13/0.20946; // line strength/VMR-O2 [1/Pa] 
  const Numeric	S0_MPM93  =  6.140e-13;         // line strength [1/Pa] 
  const Numeric G0_MPM93  =  0.560e4;           // line width [Hz/Pa]
  const Numeric	XS0_MPM93 =  2.000;             // temperature dependence of line strength
  const Numeric	XG0_MPM93 =  0.800;             // temperature dependence of line width
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
  out3  << "O2-SelfContMPM93: (model=" << model << ") parameter values in use:\n" 
	<< " S0  = " << S0 <<  "\n"
	<< " G0  = " << G0 <<  "\n"
	<< " XS0 = " << XS0 << "\n"
	<< " XG0 = " << XG0 << "\n";
  
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;
  

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
      // continuum strength
      Numeric strength =  S0 * p_abs[i] * (1.0000 - h2o_abs[i]) * pow( th, XS0 );
      // G0 from the input has to be converted to unit GHz/hPa --> * 1.0e-7
      Numeric gamma    =  G0 * p_abs[i] * pow( th, XG0 ); // Hz
	  
      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  // the vmr of O2 will be multiplied at the stage of absorption calculation:
	  // abs / vmr * xsec.
	  xsec(s,i) +=  (4.0 * PI / SPEED_OF_LIGHT)  *              // unit factor [1/(m*Hz)] 
	                (strength / VMRISO)          *              // strength    [1]
	                ( pow( f_mono[s], (Numeric)2.) * gamma /              // line shape  [Hz]
	                  ( pow( f_mono[s], (Numeric)2.) + pow( gamma, (Numeric)2.) ) );
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
			      const Numeric     S0in,         // model parameter
			      const Numeric     G0in,         // model parameter
			      const Numeric     XS0in,        // model parameter
			      const Numeric     XG0in,        // model parameter
			      const String&     model,
			      ConstVectorView  	f_mono,
			      ConstVectorView  	p_abs,        // total pressure [Pa]
			      ConstVectorView  	t_abs,
			      ConstVectorView   h2o_abs,      // H2O VMR
			      ConstVectorView   vmr	 )    // O2 VMR
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // P. W. Rosenkranz, Chapter 2, in M. A. Janssen,
  // Atmospheric Remote Sensing by Microwave Radiometry, John Wiley & Sons, Inc., 1993
  // ftp://mesa.mit.edu/phil/lbl_rt
  const Numeric	S0_PWR93   =  1.11e-14; // [K²/(Hz*Pa*m)] line strength
  const Numeric G0_PWR93  =  5600.000;  // line width [Hz/Pa]
  const Numeric	XS0_PWR93 =  2.000;    // temperature dependence of line strength
  const Numeric	XG0_PWR93 =  0.800;    // temperature dependence of line width
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
  out3  << "O2-SelfContPWR93: (model=" << model << ") parameter values in use:\n" 
	<< " S0  = " << S0 <<  "\n"
	<< " G0  = " << G0 <<  "\n"
	<< " XS0 = " << XS0 << "\n"
	<< " XG0 = " << XG0 << "\n";
  

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
	  xsec(s,i) += S0 * p_abs[i] / pow( t_abs[i], XS0 ) * 
                        ( pow( f_mono[s], (Numeric)2. )
                          * gamma / ( pow( f_mono[s], 2 )
                                      + pow( gamma, (Numeric)2. ) ) ) ;
	}
    }
}
//
//
// #################################################################################
//
/**
   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            O2-continuum according to Rosenkranz 1993 [1/m]
   \param    Cin            O2-continuum coefficient                  [1/(Hz*Pa*m)]
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
			    const String&     model,        // model parameter
			    ConstVectorView   f_mono,       // frequency grid
			    ConstVectorView   p_abs,        // P_tot grid
			    ConstVectorView   t_abs,        // T grid
			    ConstVectorView   h2o_abs,      // VMR H2O profile
			    ConstVectorView   vmr	)   // VMR O2  profile
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // P. W. Rosenkranz, Chapter 2, in M. A. Janssen,
  // Atmospheric Remote Sensing by Microwave Radiometry, John Wiley & Sons, Inc., 1993
  // ftp://mesa.mit.edu/phil/lbl_rt
  const Numeric	C_PWR93    = (1.108e-14/pow((Numeric)3.0e2,(Numeric)2.)); // [1/(Hz*Pa*m)] line strength
  const Numeric G0_PWR93   = 5600.000;  // line width [Hz/Pa]
  const Numeric G0A_PWR93  =    1.000;  // line width [1]
  const Numeric G0B_PWR93  =    1.100;  // line width [1]
  const Numeric	XG0d_PWR93 =    0.800;  // temperature dependence of line width [1]
  const Numeric	XG0w_PWR93 =    1.000;  // temperature dependence of line width [1]
  //
  // standard values for the MPM93 model (J. Liebe and G. A. Hufford and M. G. Cotton,
  // "Propagation modeling of moist air and suspended water/ice
  // particles at frequencies below 1000 GHz",
  // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  // const Numeric	C_MPM93    = 1.23e-19; // line strength/VMR [1/m*1/Hz*1/Pa] 
  const Numeric	C_MPM93    = 6.14e-13*(4.0*PI/SPEED_OF_LIGHT)/0.2085; // line strength [1/m*1/Hz*1/Pa]
  // 0.2085 = VMR * ISORATIO = 0.20946 * 0.99519
  const Numeric G0_MPM93   = 5600.000; // line width [Hz/Pa]
  const Numeric G0A_MPM93  =    1.000; // line width [1]
  const Numeric G0B_MPM93  =    1.000; // line width [1]
  const Numeric	XG0d_MPM93 =    0.800; // temperature dependence of line strength [1]
  const Numeric	XG0w_MPM93 =    0.800; // temperature dependence of line width    [1]
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
  out3  << "O2-GeneralCont: (model=" << model << ") parameter values in use:\n" 
	<< " C    = " << C <<  "\n"
	<< " G0   = " << G0 <<  "\n"
	<< " G0A  = " << G0A <<  "\n"
	<< " G0B  = " << G0B <<  "\n"
	<< " XG0d = " << XG0d << "\n"
	<< " XG0w = " << XG0w << "\n";
  

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.nelem() );
  assert ( n_p==vmr.nelem()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // const = VMR * ISORATIO = 0.20946 * 0.99519
  // this constant is already incorporated into the line strength, so we 
  // have top devide the line strength by this value since arts multiplies xsec
  // by these variables later in absCalc.
  const Numeric	VMRISO = 0.2085;

  // loop over all pressure levels:
  for ( Index i=0; i<n_p; ++i )
    {
      Numeric TH = 3.0e2 / t_abs[i];         // relative temperature  [1]
      
      Numeric ph2o = p_abs[i] * h2o_abs[i];  // water vapor partial pressure [Pa]
      Numeric pdry = p_abs[i] - ph2o;        // dry air partial pressure     [Pa]


      // pseudo broadening term [Hz]
      Numeric gamma = G0 * (G0A * pdry * pow( TH, XG0d ) + G0B * ph2o * pow( TH, XG0w )); 

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
	{
	  // division by vmr of O2 is necessary because of the absorption calculation
          // abs = vmr * xsec.
	  xsec(s,i) += C * p_abs[i] * pow( TH, (Numeric)2. ) * 
                       ( gamma * pow( f_mono[s], (Numeric)2. ) / 
                         ( pow( f_mono[s], 2 ) + pow( gamma, (Numeric)2. ) ) );
	}
    }
}
//
// #################################################################################
// ################################ NITROGEN MODELS ################################
// #################################################################################
//
// Borysow-Frommhold 1986 N2-N2 CIA absorption model;
// see publication A. Borysow and L. Frommhold, 
//                 The Astrophysical Journal, vol.311, pp.1043-1057, 1986
//                 see http://adsabs.harvard.edu/article_service.html for a scanned 
//                 version of the paper
/**

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
                            N2-CIA according to BF-86 model      [1/m]
   \param    Cin            strength scaling factor              [1]
   \param    model          allows user defined input parameter set 
                            (Cin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid            [Hz]
   \param    p_abs          predefined pressure grid             [Pa]
   \param    t_abs          predefined temperature grid          [K] 
   \param    vmr            N2 volume mixing ratio profile       [1]

   \note     this "crude" version of the N2-N2 model is a f2c 
             conversion of the N2-N2 F77 code of Prof. A. Borysow.
             The original code can be downloaded at 
             <a href="http://www.astro.ku.dk/~aborysow/">F77 code</a>.

   \remark   Reference: A. Borysow and L. Frommhold, 
//           The Astrophysical Journal, vol.311, pp.1043-1057, 1986
//           see <a href="http://adsabs.harvard.edu/article_service.html">for a scanned 
//           version of the paper</a>.

   \author Thomas Kuhn
   \date 2002-03-05
 */ 

void BF86_CIA_N2( MatrixView          xsec,
		  const Numeric       Cin,
		  const String&       model,
		  ConstVectorView     f_mono,
		  ConstVectorView     p_abs,
		  ConstVectorView     t_abs,
		  ConstVectorView     vmr   )
{
  //
  //
  // external function to call (original F77 code translated with f2c)
  extern Numeric n2n2tks_(double t, double f);
  //
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 H2O continuum model 
  // (AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  Numeric  XFAC  =  1.0000;             // scaling factor
  // ---------------------------------------------------------------------------------------
  
  // select the parameter set (!!model dominates values!!):
  if ( model == "BF86" )
    {
      XFAC =  1.0000;
    }
  else if ( model == "user" )
    {
      XFAC =  Cin;
    }
  else
    {
      ostringstream os;
      os << "N2-SelfContBorysow: ERROR! Wrong model values given.\n"
	 << "allowed models are: 'BF86', 'user'" << '\n';
      throw runtime_error(os.str());
    }
  
  out3  << "N2-SelfContBorysow: (model=" << model << ") parameter values in use:\n" 
	<< " XFAC = " << XFAC << "\n";
  
  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies


  const Numeric AMAG2DEN = 44.53807; // inverse of N2 mol volume at std p/T
  const Numeric RIDGAS =  8.314510;  // ideal gas constant

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
      //cout << "------------------------------------------------\n";
      double T = (double) t_abs[i];
      //cout << "N2-N2 BF86: T     =" << T << " K\n";
      //cout << "N2-N2 BF86: p     =" << p_abs[i] << " Pa\n";
      //cout << "N2-N2 BF86: VMR   =" << vmr[i] << "\n";
      Numeric XAMA  = (p_abs[i]) / ( AMAG2DEN * RIDGAS * t_abs[i] );
      Numeric XAMA2 = pow(XAMA,(Numeric)2.);
      //cout << "N2-N2 BF86: XAMA  =" << XAMA << "\n";

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  // the second vmr of N2 will be multiplied at the stage of
	  // absorption calculation: abs =  vmr * xsec.
	  double f = (double) f_mono[s];
	  //cout << "N2-N2 BF86: f     =" << f << " Hz\n";
	  double cont = n2n2tks_(T, f);
	  xsec(s,i) += (Numeric) (cont * 1.000e2 * vmr[i] * XAMA2);
	  //cout << "N2-N2 BF86: cont  =" << cont << " cm-1 * amagat-2\n";
	  //cout << "N2-N2 BF86: abs   =" << (vmr[i] * xsec(s,i)) << " m-1\n";
	}
    }
  return;
}
//
// #################################################################################
//
// MPM93 N2 continuum:
// see publication side of National Telecommunications and Information Administration
//   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
// and ftp side for downloading the MPM93 original source code:
//   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
			 ConstVectorView     vmr	 )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 H2O continuum model 
  // (AGARD 52nd Specialists Meeting of the Electromagnetic Wave
  // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21):
  const Numeric	xT_MPM93  =  3.500;             // temperature exponent [1]
  const Numeric	xf_MPM93  =  1.500;             // frequency exponent [1]
  const Numeric	gxf_MPM93 =  9.000*xf_MPM93;    // needed for the unit conversion of G_MPM93
  const Numeric	S_MPM93   =  2.296e-31;         // line strength  [1/Pa² * 1/Hz]
  const Numeric G_MPM93   =  1.930e-5*pow((Numeric)10.000, -gxf_MPM93); // frequency factor [1/Hz^xf]
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
  else if ( model == "MPM93Scale" )
    {
      S0  = Cin * S_MPM93;
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
	 << "allowed models are: 'MPM93', 'MPM93Scale' or 'user'" << '\n';
      throw runtime_error(os.str());
    }
  
  out3  << "N2-SelfContMPM93: (model=" << model << ") parameter values in use:\n" 
	<< " S0 = " << S0 << "\n"
	<< " G0 = " << G0 << "\n"
	<< " xT = " << xT << "\n"
	<< " xf = " << xf << "\n";
  
  // unit conversion internally:
  //const Numeric S0unitconv = 1.000e+13;  // x [1/(hPa²*GHz)] => y [1/(pa²*Hz)]
  //const Numeric G0unitconv = pow(10.000, gxf);

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
                          pow( (p_abs[i] * ((Numeric)1.0000 - h2o_abs[i])),
                               (Numeric)2. )
                          * pow( th, xT );

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  Numeric f = f_mono[s] * Hz_to_GHz; // frequency in GHz
	  // the vmr of N2 will be multiplied at the stage of absorption calculation:
	  // abs / vmr * xsec.
	  xsec(s,i) += fac * strength *                              // strength
                       pow(f_mono[s], (Numeric)2.) /                      // frequency dependence
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
   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    model          allows user defined input parameter set 
                            (Cin and xTin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid      [Hz]
   \param    p_abs          predefined pressure grid       [Pa]
   \param    t_abs          predefined temperature grid    [K] 
   \param    vmr            N2  volume mixing ratio        [1]
   \param    h2ovmr         H2O volume mixing ratio        [1]

   \note     Except for  model 'user' the input parameters Cin and xTin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'ATM', and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference: Pardo et al., IEEE, Trans. Ant. Prop., <br>
             Vol 49, No 12, pp. 1683-1694, 2001.

   \author Thomas Kuhn
   \date 2001-04-10
 */ 

void Pardo_ATM_N2_dry_continuum( MatrixView          xsec,
				 const Numeric       Cin,
				 const String&       model,
				 ConstVectorView     f_mono,
				 ConstVectorView     p_abs,
				 ConstVectorView     t_abs,
				 ConstVectorView     vmr,
				 ConstVectorView     h2ovmr	 )
{
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Pardo et al. model (IEEE, Trans. Ant. Prop., 
  // Vol 49, No 12, pp. 1683-1694, 2001)
  const Numeric	C_ATM = 2.612e-6; // [1/m]
  // ---------------------------------------------------------------------------------------

  // select the parameter set (!!model dominates parameters!!):
  Numeric C;
   if ( model == "ATM" )
     {
       C = C_ATM;
     }
   else if ( model == "user" )
     {
       C = Cin;
     }
   else
     {
       ostringstream os;
       os << "N2-DryContATM01: ERROR! Wrong model values given.\n"
	  << "allowed models are: 'ATM', 'user'" << '\n';
       throw runtime_error(os.str());
     }
   out3  << "N2-DryContATM01: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n";

   const Index n_p = p_abs.nelem();	// Number of pressure levels
   const Index n_f = f_mono.nelem();	// Number of frequencies

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
      // since this is an effective "dry air" continuum, it is not really
      // it is not specifically attributed to N2, so we need the total 
      // dry air part in total which is equal to the total minus the 
      // water vapor pressure:
      Numeric  pd = p_abs[i] * ( 1.00000e0 - h2ovmr[i] ); // [Pa]
      // Loop over frequency grid:
      if (vmr[i] > VMRCalcLimit )
	{
	  for ( Index s=0; s<n_f; ++s )
	    {
	      // Becaue this is an effective "dry air" continuum, it is not really
	      // specific N2 but mainly caused by N2. Therefore the N2 vmr must be 
	      // canceled out here which is later in absCalc multiplied 
	      // (calculation: abs = vmr * xsec):
	      xsec(s,i) += C *                    // strength [1/(m*Hz²Pa²)] 
		pow( (f_mono[s]/(Numeric)2.25e11), (Numeric)2. ) * // quadratic f dependence [Hz²]
		pow( ((Numeric)300.0/t_abs[i]), (Numeric)3.5 )   * // free T dependence      [1]
		pow( (pd/(Numeric)1.01300e5), (Numeric)2. )      / // quadratic p dependence [Pa²]
		vmr[i];                                            // cancel the vmr dependency
	    }
	}
    }
}
//
// #################################################################################
//
/** 
   \retval   xsec           cross section (absorption/volume mixing ratio) of 
                            N2-continuum according to Rosenkranz, 1993 [1/m]
   \param    Cin            continuum strength [1/m * 1/(Hz*Pa)²]
   \param    xin            temperature exponent of N2-continuum [1]
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
				   ConstVectorView     vmr	 )
{
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model (Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric	C_PWR = 1.05e-38; // [1/(Pa²*Hz²*m)]
  const Numeric	x_PWR = 3.55;     // [1]
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
   out3  << "N2-SelfContPWR93: (model=" << model << ") parameter values in use:\n" 
         << " C_s = " << C << "\n"
         << " x_s = " << x << "\n";

   const Index n_p = p_abs.nelem();	// Number of pressure levels
   const Index n_f = f_mono.nelem();	// Number of frequencies

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
	               pow( f_mono[s], (Numeric)2. ) *      // quadratic f dependence [Hz²]
	               pow( (Numeric)300.0/t_abs[i], x ) * // free T dependence      [1]
                       pow( p_abs[i], (Numeric)2. ) *       // quadratic p dependence [Pa²]
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
				 ConstVectorView     vmr	 )
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the Rosenkranz model, Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric	C_GM  = 1.05e-38; // [1/(Pa²*Hz²*m)]
  const Numeric	xf_GM = 2.00;     // [1]
  const Numeric	xt_GM = 3.55;     // [1]
  const Numeric	xp_GM = 2.00;     // [1]
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
   out3  << "N2-SelfContStandardType: (model=" << model << ") parameter values in use:\n" 
         << " C  = " << C  << "\n"
         << " xt = " << xt << "\n"
         << " xf = " << xf << "\n"
         << " xp = " << xp << "\n";


  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
                       pow( ((Numeric)300.00/t_abs[i]), xt ) * // T dependence        [1] 
                       pow( f_mono[s], xf )         * // f dependence    [Hz^xt]  
                       pow( p_abs[i], xp )          * // p dependence    [Pa^xp]
                       pow( vmr[i], (xp-(Numeric)1.) );         // last N2-VMR at the stage 
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
   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
				    ConstVectorView     vmr	 )
{
  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
  // "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
  const Numeric	C_PWR = 7.43e-37; // [ 1/(Pa²*Hz²*m) ]
  const Numeric	x_PWR = 5.08;     // [ 1 ]
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
  
  out3  << "CO2-SelfContPWR93: (model=" << model << ") parameter values in use:\n" 
	<< " C = " << C << "\n"
	<< " x = " << x << "\n";

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
	C * pow( (Numeric)300./t_abs[i], x ) * pow( p_abs[i], (Numeric)2. ) * vmr[i];

      // Loop over frequency grid:
      for ( Index s=0; s<n_f; ++s )
	{
	  xsec(s,i) += dummy * pow( f_mono[s], (Numeric)2. );
	}
    }
}
//
// #################################################################################
//
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
				       ConstVectorView     vmr	 )
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
  
  out3  << "CO2-ForeignContPWR93: (model=" << model << ") parameter values in use:\n" 
	<< " C = " << C << "\n"
	<< " x = " << x << "\n";

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      Numeric dummy = C * pow( (Numeric)300./t_abs[i], x ) * p_abs[i] * p_abs[i] * n2_abs[i];

      // Loop frequency:
      for ( Index s=0; s<n_f; ++s )
	{
	  xsec(s,i) += dummy * pow( f_mono[s], (Numeric)2. );
	}
    }
}
//
// #################################################################################
// ################################### CLOUD AND RAIN MODELS #######################
// #################################################################################
//
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
  out3  << "liquidcloud-MPM93: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CG = " << CG << "\n"
	<< " CE = " << CE << "\n";
  
  
  const Numeric m = 1.00e3; // specific weight of the droplet,  fixed value:  1.00e3 kg/m³
  const Numeric low_lim_den  =  0.000;   // lower limit of suspended droplet particle density vector [kg/m³]
  const Numeric high_lim_den = 10.00e-3; // lower limit of suspended droplet particle density vector [kg/m³]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
		pow((f_mono[s]*Hz_to_GHz),(Numeric)2.) *
		( ((epsilon0-epsilon1)/
		   (pow((f_mono[s]*Hz_to_GHz),(Numeric)2.)
                    + pow(gamma1,(Numeric)2.))) + 
		  ((epsilon1-epsilon2)/
		   (pow((f_mono[s]*Hz_to_GHz),(Numeric)2.)
                    + pow(gamma2,(Numeric)2.))) );
	      // imaginary part of the complex permittivity of water (double-debye model)
	      Numeric Imepsilon  = (f_mono[s]*Hz_to_GHz) *
		( (gamma1*(epsilon0-epsilon1)/
		   (pow((f_mono[s]*Hz_to_GHz),(Numeric)2.)
                    + pow(gamma1,(Numeric)2.))) + 
		  (gamma2*(epsilon1-epsilon2)/
		   (pow((f_mono[s]*Hz_to_GHz),(Numeric)2.)
                    + pow(gamma2,(Numeric)2.))) );
	      // the imaginary part of the complex refractivity of suspended liquid water particle [ppm]
	      // In MPM93 w is in g/m³ and m is in g/cm³. Because of the units used in arts,
	      // a factor of 1.000e6 must be multiplied with the ratio (w/m):
	      // MPM93: (w/m)_MPM93  in   (g/m³)/(g/cm³)
	      // arts:  (w/m)_arts   in  (kg/m³)/(kg/m³)
	      // =====> (w/m)_MPM93   =   1.0e6 * (w/m)_arts
	      // the factor of 1.0e6 is included below in xsec calculation.
	      Numeric ImNw = 1.500 / m * 
		( 3.000 * Imepsilon
                  / ( pow((Reepsilon+(Numeric)2.000),(Numeric)2.)
                      + pow(Imepsilon,(Numeric)2.) ) );
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

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
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
			 ConstVectorView   vmr	 ) // suspended ice particle density vector, 
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
  out3  << "icecloud-MPM93: (model=" << model << ") parameter values in use:\n" 
	<< " CC = " << CC << "\n"
	<< " CA = " << CA << "\n"
	<< " CB = " << CB << "\n";
  
  
  const Numeric m = 0.916e3;  // specific weight of ice particles,  fixed value:   0.916e3 kg/m³
  const Numeric low_lim_den  =  0.000;   // lower limit of suspended ice particle density vector [kg/m³]
  const Numeric high_lim_den = 10.00e-3; // lower limit of suspended ice particle density vector [kg/m³]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
		       ( -24.17 + (116.79/theta)
                         + pow((theta/(theta-(Numeric)0.9927)),(Numeric)2.) );
	      
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
		    ( 3.000 * Imepsilon
                      / ( pow((Reepsilon+(Numeric)2.000),(Numeric)2.)
                          + pow(Imepsilon,(Numeric)2.) ) );
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
//
/** 

   \retval    xsec          cross section (absorption/volume mixing ratio) of 
                            water clouds according to MPM93 [1/m]
   \param    CEin           scaling parameter of the calculated cross section [1]
   \param    CAin           scaling parameter of the factor a_rain [1]
   \param    CBin           scaling parameter of the exponent b_rain [1]
   \param    model          allows user defined input parameter 
                            (CEin, CAin, CBin)<br> or choice of 
                            pre-defined parameters of specific models (see note below).
   \param    f_mono         predefined frequency grid       [Hz]
   \param    p_abs          predefined pressure grid        [Pa]
   \param    t_abs          predefined temperature grid     [K] 
   \param    vmr            rain rate vector (i.e. vertical profile),
                            (valid range: 0-150) [mm/h]

   \note     Except for  model 'user' the input parameters CEin, CAin, and CBin 
             are neglected (model dominates over parameters).<br>
             Allowed models: 'MPM93' and 'user'. 
             See the user guide for detailed explanations.

   \remark   Reference:  R. L. Olsen and D.V. Rogers and D. B. Hodge,<br>
             <i> The aR^b relation in the calculation of rain attenuation</i>,<br>
	     IEEE Trans. Antennas Propagat., vol. AP-26, pp. 318-329, 1978.

   \author Christian Melsheimer
   \date 2003-22-05
 */ 

void MPM93RainExt( MatrixView         xsec,
		   const Numeric      CEin,   // input parameter
		   const Numeric      CAin,   // input parameter
		   const Numeric      CBin,   // input parameter
		   const String&      model, // model
		   ConstVectorView    f_mono, // frequency vector
		   ConstVectorView    p_abs,  // pressure vector
		   ConstVectorView    t_abs,  // temperature vector
		   ConstVectorView    vmr)    // rain rate profile [mm/h]
{

  // --------- STANDARD MODEL PARAMETERS ---------------------------------------------------
  // standard values for the MPM93 model based on Olsen, R.L.,
  // D.V. Rogers, and D. B. Hodge, "The aR^b relation in the
  // calculation of rain attenuation", IEEE Trans. Antennas Propagat.,
  // vol. AP-26, pp. 318-329, 1978,
  const Numeric CE_MPM93 = 1.00000;
  const Numeric CA_MPM93 = 1.00000;
  const Numeric CB_MPM93 = 1.00000;
  // ---------------------------------------------------------------------------------------


  // select the parameter set (!!model dominates values!!):
  Numeric CE, CA, CB;
  if ( model == "MPM93" )
    {
      CE = CE_MPM93;
      CA = CA_MPM93;
      CB = CB_MPM93;
    }
  else if ( model == "user" )
    {
      CE = CEin;
      CA = CAin;
      CB = CBin;
    }
  else
    {
      ostringstream os;
      os << "rain-MPM93: ERROR! Wrong model values given.\n"
	 << "Valid models are: 'MPM93' and 'user'" << '\n';
      throw runtime_error(os.str());
    }
  out3  << "rain-MPM93: (model=" << model << ") parameter values in use:\n" 
	<< " CE = " << CE << "\n"
	<< " CA = " << CA << "\n"
	<< " CB = " << CB << "\n";
  
  
  const Numeric low_lim_rr  =  0.000;   // lower limit of allowed rain rate  [mm/h]
  const Numeric high_lim_rr = 150.000;  // upper limit of allowed rain rate  [mm/h]

  const Index n_p = p_abs.nelem();	// Number of pressure levels
  const Index n_f = f_mono.nelem();	// Number of frequencies

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
      // Extinction by rain is parameterized as:
      //  ext_rain = a_rain * rr ^ b_rain
      // a_rain and b_rain each depend on frequency by  power laws:
      //  a_rain = Ga * freq ^ Ea
      //  b_rain = Gb * freq ^ Eb
      Numeric Ga;
      Numeric Ea;
      Numeric Gb;
      Numeric Eb;
      Numeric a_rain;
      Numeric b_rain;
      Numeric ext_rain;

      // Check limits of rain rate ("vmr") [mm/h]
      if ( (vmr[i] >= low_lim_rr) && (vmr[i] < high_lim_rr) ) 
	{
	  // Loop frequency:
	  for ( Index s=0; s<n_f; ++s )
	    {
	      // for rain rate < 25 mm/h, take parameters from Olsen et al.'s
	      // own power law fit to their Laws-Parsons-Low data;
	      // for rain rate > 25 mm/h, take C. Melsheimer's power law fit
	      // to Olsen et al.'s Laws-Parson-High data
	      if ( vmr[i] <= 25 ) 
		{	      
		  // power law coeff. Ga and exponent Ea for a, piecewise:
		  if ( f_mono[s] <= 2.9e9 ) 
		    {
		      Ga = 6.39e-5;
		      Ea = 2.03;
		    }
		  else if ( f_mono[s] <= 54.0e9 )
		    {
		      Ga = 4.21e-5;
		      Ea = 2.42;
		    }
		  else if ( f_mono[s] <= 180e9 )
		    {
		      Ga = 4.09e-2;
		      Ea = 0.699;
		    }
		  else if ( f_mono[s] <= 1000e9 )
		    {
		      Ga = 3.38;
		      Ea = -0.151;
		    }
		  else 
		    {
		      ostringstream os;
		      os << "ERROR in MPM93RainExt:\n"
			 << " frequency (valid range 0-1000 GHz):" << f_mono[s]*Hz_to_GHz << "\n"
			 << " ==> no calculation performed!\n";
		      throw runtime_error(os.str());
		    }
		  // power law coeff. Gb and exponent Eb for b, piecewise:
		  if ( f_mono[s] <= 8.5e9 ) 
		    {
		      Gb = 0.851;
		      Eb = 0.158;
		    }
		  else if ( f_mono[s] <= 25.0e9 )
		    {
		      Gb = 1.41;	   
		      Eb = -0.0779;
		    }
		  else if ( f_mono[s] <= 164.0e9 )
		    {
		      Gb = 2.63;	  
		      Eb = -0.272;
		    }
		  else if ( f_mono[s] <= 1000e9 )
		    {
		      Gb = 0.616;
		      Eb = 0.0126;
		    }
		  else 
		    {
		      ostringstream os;
		      os << "ERROR in MPM93RainExt:\n"
			 << " frequency (valid range 0-1000 GHz):" << f_mono[s]*Hz_to_GHz << "\n"
			 << " ==> no calculation performed!\n";
		      throw runtime_error(os.str());
		    }
		  
		}
	      else if (vmr[i] > 25)
		{
		  // power law coeff. Ga and exponent Ea for a, piecewise:
		  if ( f_mono[s] <= 4.9e9 ) 
		    {
		      Ga = 5.30e-5;
		      Ea = 1.87;
		    }
		  else if ( f_mono[s] <= 10.7e9 )
		    {
		      Ga = 5.03e-6;
		      Ea = 3.35;
		    }
		  else if ( f_mono[s] <= 40.1e9 )
		    {
		      Ga = 2.53e-5;
		      Ea = 2.67;
		    }
		  else if ( f_mono[s] <= 59.1e9 )
		    {
		      Ga = 3.58e-3;
		      Ea = 1.33;
		    }
		  else if ( f_mono[s] <= 100e9 )
		    {
		      Ga = 0.143;
		      Ea = 0.422;
		    }
		  else 
		    {
		      ostringstream os;
		      os << "ERROR in MPM93RainExt:\n"
			 << " frequency (valid range for rain rate > 25mm/h: 0-100 GHz):" << f_mono[s]*Hz_to_GHz << "\n"
			 << " ==> no calculation performed!\n";
		      throw runtime_error(os.str());
		    }
		  // power law coeff. Gb and exponent Eb for b, piecewise:
		  if ( f_mono[s] <= 6.2e9 ) 
		    {
		      Gb = 0.911;
		      Eb = 0.190;
		    }
		  else if ( f_mono[s] <= 23.8e9 )
		    {
		      Gb = 1.71;	   
		      Eb = -0.156;
		    }
		  else if ( f_mono[s] <= 48.4e9 )
		    {
		      Gb = 3.08;	  
		      Eb = -0.342;
		    }
		  else if ( f_mono[s] <= 68.2e9 )
		    {
		      Gb = 1.28;
		      Eb = -0.116;
		    }
		  else if ( f_mono[s] <= 100e9 )
		    {
		      Gb =  0.932;
		      Eb =  -0.0408;
		    }
		  else 
		    {
		      ostringstream os;
		      os << "ERROR in MPM93RainExt:\n"
			 << " frequency (valid range for rain rate > 25mm/h: 0-100 GHz):" << f_mono[s]*Hz_to_GHz << "\n"
			 << " ==> no calculation performed!\n";
		      throw runtime_error(os.str());
		    }
		}
	      //Factor a_rain
	      Numeric a_rain = Ga * pow((f_mono[s]*Hz_to_GHz),Ea);
	      //Factor b_rain
	      Numeric b_rain = Gb * pow((f_mono[s]*Hz_to_GHz),Eb);
	      // Extinction coefficient [dB/km], with scaling
	      // parameters CA and CB
	      Numeric ext_rain = CA * a_rain * pow(vmr[i],(CB*b_rain));
	      // rain extinction cross section [1/m]
	      // The vmr will be multiplied at the stage of extinction
	      // calculation: ext = vmr * xsec.
	      // xsec = ext/vmr [1/m] but MPM93 is in [dB/km] --> conversion necessary
	      xsec(s,i) += CE * dB_km_to_1_m * ext_rain / vmr[i];
	    }
	} else
	  {
	    if ( (vmr[i] < low_lim_rr) || (vmr[i] > high_lim_rr) ) 
	      {
		ostringstream os;
		os << "ERROR in MPM93RainExt:\n"
		   << " rain rate (valid range 0.00-150.00 mm/h):" << vmr[i] << "\n"
		   << " ==> no calculation performed!\n";
		throw runtime_error(os.str());
	      }
	  }
    }
  
}
//
// #################################################################################
// ################################# HELP FUNCTIONS ################################
// #################################################################################
//
/** 

   \retval   MPMLineShapeFunction  H2O-line shape function value     [1/Hz]  
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

   \retval   MPMLineShapeO2Function  O2-line shape function value         [1]  
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

   \retval   WVSatPressureLiquidWater     water vapor saturation pressure over liquid water [Pa]  
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
		        1.3816e-7 * ( pow( (Numeric)10.00,
                                           ((Numeric)11.344*
                                            ((Numeric)1.00-((Numeric)1.00
                                                             /theta))) )
                                      - (Numeric)1.000 ) +
		        8.1328e-3 * ( pow( (Numeric)10.00,
                                           ((Numeric)-3.49149
                                            *(theta-(Numeric)1.00)))
                                      - 1.000) +
		        log10(1013.246) );
  Numeric es_MPM93 = 100.000 * pow((Numeric)10.00,exponent);

  return es_MPM93; // [Pa]
}
//
// #################################################################################
//
/** 

   \retval   WVSatPressureIce     water vapor saturation pressure over liquid water [Pa]  
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

  Numeric es_MPM93 = 100.000 * pow((Numeric)10.00,exponent);

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

    \retval xsec       Cross section of one continuum tag,<br> 
                       xsec = alpha / VMR  [1/m * 1]

    \param  name       The name of the model to calculate (derived from the tag name)
    \param  parameters model parameters, as defined in method
                       cont_description_parameters.
    \param  model      model, related to model parameters
    \param  f_mono     Frequency grid [Hz]
    \param  p_abs      Pressure grid [Pa]
    \param  t_abs      Temperatures associated with the pressure grid, p_abs [K]
    \param  n2_abs     Total volume mixing ratio profile of molecular nitrogen.<br> 
                       This will be needed only for the CO2 foreign continuum [1]<br> 
                       however one is forced to give this input [1]
    \param  h2o_abs    Total volume mixing ratio profile of water vapor.<br> 
                       This will be needed only for the oxygen continuum <br> 
                       however one is forced to give this input [1]
    \param  vmr        Volume mixing ratio profile of the actual species [1]

   \author Stefan Bühler, Thomas Kuhn
   \date 2001-11-05
 */

void xsec_continuum_tag( MatrixView                 xsec,
			 const String&              name,
			 ConstVectorView            parameters,
                         const String&              model,
			 ConstVectorView  	    f_mono,
			 ConstVectorView  	    p_abs,
			 ConstVectorView  	    t_abs,
			 ConstVectorView  	    n2_abs,
			 ConstVectorView  	    h2o_abs,
			 ConstVectorView            vmr )
{
  //
  /* In the following all the possible tags are listed here and 
   after a first consistency check about the input parameters the 
   appropriate internal function is called,

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ATTENTION PLEASE UPDATE THIS COMMENT IF ANY CHANGES ARE MADE CONCERNING 
  THE ASSOCIATED MODELS TO EACH TAG
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ----------------------------------------------------------------------------------------------------
                    TAG                      VALID MODELS  
  ----------------------------------------------------------------------------------------------------
  *CONTAGMODINFO*   H2O-SelfContStandardType: Rosenkranz, user
  *CONTAGMODINFO*   H2O-ForeignContStandardType: Rosenkranz, user
  *CONTAGMODINFO*   H2O-ForeignContMaTippingType: MaTipping, user
  *CONTAGMODINFO*   H2O-MPM87:               MPM87, MPM87Lines, MPM87Continuum, user
  *CONTAGMODINFO*   H2O-MPM89:               MPM89, MPM89Lines, MPM89Continuum, user
  *CONTAGMODINFO*   H2O-MPM93:               MPM93, MPM93Lines, MPM93Continuum, user
  *CONTAGMODINFO*   H2O-PWR98:               Rosenkranz, RosenkranzLines, RosenkranzContinuum, user
  *CONTAGMODINFO*   H2O-CP98:                CruzPol, CruzPolLine, CruzPolContinuum, user
  *CONTAGMODINFO*   H2O-CKD24:               CKD24, user
  *CONTAGMODINFO*   O2-MPM85:                MPM85, MPM85Lines, MPM85Continuum, MPM85NoCoupling, MPM85NoCutoff, user
  *CONTAGMODINFO*   O2-MPM87:                MPM87, MPM87Lines, MPM87Continuum, MPM87NoCoupling, MPM87NoCutoff, user
  *CONTAGMODINFO*   O2-MPM89:                MPM89, MPM89Lines, MPM89Continuum, MPM89NoCoupling, MPM89NoCutoff, user
  *CONTAGMODINFO*   O2-MPM92:                MPM92, MPM92Lines, MPM92Continuum, MPM92NoCoupling, MPM92NoCutoff, user
  *CONTAGMODINFO*   O2-MPM93:                MPM93, MPM93Lines, MPM93Continuum, MPM93NoCoupling, MPM92NoCutoff, user
  *CONTAGMODINFO*   O2-PWR93:                Rosenkranz, RosenkranzLines, RosenkranzContinuum, user
  *CONTAGMODINFO*   O2-PWR88:                Rosenkranz, RosenkranzLines, RosenkranzContinuum, user
  *CONTAGMODINFO*   O2-SelfContMPM93:        MPM93, user
  *CONTAGMODINFO*   O2-SelfContPWR93:        Rosenkranz, user
  *CONTAGMODINFO*   O2-GenerealCont:         Rosenkranz, MPM93, user
  *CONTAGMODINFO*   N2-BFCIA86:              BF86, user
  *CONTAGMODINFO*   N2-SelfContMPM93:        MPM93, user, MPM93Scale
  *CONTAGMODINFO*   N2-SelfContPWR93:        Rosenkranz, user
  *CONTAGMODINFO*   N2-SelfContStandardType: Rosenkranz, user
  *CONTAGMODINFO*   CO2-SelfContPWR93:       Rosenkranz, user
  *CONTAGMODINFO*   CO2-ForeignContPWR93:    Rosenkranz, user
  *CONTAGMODINFO*   liquidcloud-MPM93:       MPM93, user
  *CONTAGMODINFO*   icecloud-MPM93:          MPM93, user
  *CONTAGMODINFO*   rain-MPM93:              MPM93, user
  ----------------------------------------------------------------------------------------------------
  */
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContMaTippingType"==name )
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
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  MaTipping_H2O_foreign_continuum( xsec,
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  MaTipping_H2O_foreign_continuum( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
			       vmr	 );
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
	  out3 << "Continuum model " << name << " running with \n" 
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
			       vmr	 );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: Continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters. " << "\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContATM01"==name )
    {
      // Foreign wet continuum term.
      //
      // Pardo et al., IEEE, Trans. Ant. Prop., 
      // Vol 49, No 12, pp. 1683-1694, 2001.
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : pseudo continuum line frequency                      [Hz]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  Pardo_ATM_H2O_ForeignContinuum( xsec,
					  parameters[0],
					  model,
					  f_mono,
					  p_abs,
					  t_abs,
					  vmr	 );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  Pardo_ATM_H2O_ForeignContinuum( xsec,
					  0.000,
					  model,
					  f_mono,
					  p_abs,
					  t_abs,
					  vmr	 );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: Continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters. " << "\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-SelfContCKD222"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O self continuum according to CKD2.2.2    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_222_self_h2o( xsec,
			    parameters[0],
			    model,
			    f_mono,
			    p_abs,
			    t_abs,
			    vmr,
			    n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_222_self_h2o( xsec,
			    0.000,
			    model,
			    f_mono,
			    p_abs,
			    t_abs,
			    vmr,
			    n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
       {
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContCKD222"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O foreign continuum according to CKD2.2.2    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_222_foreign_h2o( xsec,
			       parameters[0],
			       model,
			       f_mono,
			       p_abs,
			       t_abs,
			       vmr,
			       n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_222_foreign_h2o( xsec,
			       0.000,
			       model,
			       f_mono,
			       p_abs,
			       t_abs,
			       vmr,
			       n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
       {
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-SelfContCKD242"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O self continuum according to CKD2.4.2    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_242_self_h2o( xsec,
			    parameters[0],
			    model,
			    f_mono,
			    p_abs,
			    t_abs,
			    vmr,
			    n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_242_self_h2o( xsec,
			    0.000,
			    model,
			    f_mono,
			    p_abs,
			    t_abs,
			    vmr,
			    n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
       {
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContCKD242"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O foreign continuum according to CKD2.4.2    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_242_foreign_h2o( xsec,
			       parameters[0],
			       model,
			       f_mono,
			       p_abs,
			       t_abs,
			       vmr,
			       n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_242_foreign_h2o( xsec,
			       0.000,
			       model,
			       f_mono,
			       p_abs,
			       t_abs,
			       vmr,
			       n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-SelfContCKDMT100"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O self continuum according to CKD MT 1.00    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_100_self_h2o( xsec,
			       parameters[0],
			       model,
			       f_mono,
			       p_abs,
			       t_abs,
			       vmr,
			       n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_100_self_h2o( xsec,
			       0.000,
			       model,
			       f_mono,
			       p_abs,
			       t_abs,
			       vmr,
			       n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
       {
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContCKDMT100"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O foreign continuum according to CKD MT 1.00    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_100_foreign_h2o( xsec,
				  parameters[0],
				  model,
				  f_mono,
				  p_abs,
				  t_abs,
				  vmr,
				  n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_100_foreign_h2o( xsec,
				  0.000,
				  model,
				  f_mono,
				  p_abs,
				  t_abs,
				  vmr,
				  n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-SelfContCKD24"==name )
    {
      // OUTPUT:
      //   xsec           cross section (absorption/volume mixing ratio) of 
      //                  H2O continuum according to CKD2.4    [1/m]
      // INPUT:
      //   parameters[0]  strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD24_H20( xsec,
		     0,
		     parameters[0],
		     model,
		     f_mono,
		     p_abs,
		     t_abs,
		     vmr,
		     n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD24_H20( xsec,
		     0,
		     0.000,
		     model,
		     f_mono,
		     p_abs,
		     t_abs,
		     vmr,
		     n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
       {
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
       }
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "H2O-ForeignContCKD24"==name ) 
    {
      // OUTPUT:
      //   xsec             cross section (absorption/volume mixing ratio) of 
      //                    H2O continuum according to CKD2.4    [1/m]
      // INPUT:
      //   Cin            strength scaling factor              [1]
      //   model          allows user defined input parameter set 
      //                  (Cin) or choice of 
      //                  pre-defined parameters of specific models (see note below).
      //   f_mono         predefined frequency grid            [Hz]
      //   p_abs          predefined pressure grid             [Pa]
      //   t_abs          predefined temperature grid          [K] 
      //   vmr            H2O volume mixing ratio profile      [1]
      //   n2_abs         N2 volume mixing ratio profile       [1]
      //
      // WWW resource: ftp.aer.com/aer_contnm_ckd
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD24_H20( xsec,
		     1,
		     parameters[0],
		     model,
		     f_mono,
		     p_abs,
		     t_abs,
		     vmr,
		     n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD24_H20( xsec,
		     1,
		     0,
		     model,
		     f_mono,
		     p_abs,
		     t_abs,
		     vmr,
		     n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	}
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
	  out3 << "Full model " << name << " is running with \n"
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
	  out3 << "Full model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Full model " << name << " is running with \n"
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
	  out3 << "Full model " << name << " running with \n" 
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
	     << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Full model " << name << " is running with \n"
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
	  out3 << "Full model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Full model " << name << " is running with \n"
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
	  out3 << "Full model " << name << " running with \n" 
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
	     << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Full model " << name << " is running with \n"
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
	  out3 << "Full model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ============= O2 continuum =========================================================
  else if ( "O2-CIAfunCKDMT100"==name )
    {
      // Model reference:
      // F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, 
      // J.-M. Hartmann, Ch. Boulet, 
      // "Infrared collision-induced absorption by O2 near 6.4 microns for
      // atmospheric applications: measurements and emprirical modeling", 
      // Appl. Optics, 35, 5911-5917, (1996).
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum scaling
      //     model         : model option ("CKD" or "user")
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs       : [1]
      //     vmr           : [1]
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_CIAfun_o2( xsec,
			    parameters[0],
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_CIAfun_o2( xsec,
			    0.00e0,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-v0v0CKDMT100"==name )
    {
      // Model reference:
      //   B. Mate, C. Lugez, G.T. Fraser, W.J. Lafferty,
      //   "Absolute Intensities for the O2 1.27 micron
      //   continuum absorption",  
      //   J. Geophys. Res., 104, 30,585-30,590, 1999. 
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum scaling
      //     model         : model option ("CKD" or "user")
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //     n2_abs        : [1]
      //     h2o_abs       : [1]
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_v0v0_o2( xsec,
			  parameters[0],
			  model, 
			  f_mono,
			  p_abs,
			  t_abs,
			  vmr,
			  n2_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_v0v0_o2( xsec,
			  0.0e0,
			  model, 
			  f_mono,
			  p_abs,
			  t_abs,
			  vmr,
			  n2_abs );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: Continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-v1v0CKDMT100"==name )
    {
      // Model reference:
      //   Mlawer, Clough, Brown, Stephen, Landry, Goldman, Murcray,
      //   "Observed  Atmospheric Collision Induced Absorption in Near Infrared Oxygen Bands",
      //   Journal of Geophysical Research, vol 103, no. D4, pp. 3859-3863, 1998.
      //
      //  specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT
      //     parameters[0] : continuum scaling
      //     model         : model option ("CKD" or "user")
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //     n2_abs        : [1]
      //     h2o_abs       : [1]
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_v1v0_o2( xsec,
			  parameters[0],
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_v1v0_o2( xsec,
			  0.0e0,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ============= O2 full model ========================================================
  else if ( "O2-PWR88"==name )
    {
      //  REFERENCE FOR EQUATIONS AND COEFFICIENTS:
      //  P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
      //  BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
      //  AND 
      //  H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
      //  (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
      //  AND
      //  P. W. ROSENKRANZ, INTERFERENCE COEFFICIENTS FOR THE 
      //  OVERLAPPING OXYGEN LINES IN AIR, JQSRT, 1988, VOLUME 39, 287-297.
      //
      //  the only difference to the 1993 version is the line mixing 
      //  parameter Y, which is taken from the above reference JQSRT, 1988.
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
      const char *version="PWR88";
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  PWR93O2AbsModel( xsec,
			   parameters[0], // continuum term scale factor
			   parameters[1], // line strength scale factor
			   parameters[2], // line broadening scale factor
			   parameters[3], // line coupling scale factor
			   model,
			   version, 
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  PWR93O2AbsModel( xsec,
			   0.00,
			   0.00,
			   0.00,
			   0.00,
			   model,
			   version,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      const char *version="PWR93";
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  PWR93O2AbsModel( xsec,
			   parameters[0], // continuum term scale factor
			   parameters[1], // line strength scale factor
			   parameters[2], // line broadening scale factor
			   parameters[3], // line coupling scale factor
			   model,
			   version,
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  PWR93O2AbsModel( xsec,
			   0.00,
			   0.00,
			   0.00,
			   0.00,
			   model,
			   version,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-PWR98"==name )
    {
      //  REFERENCES FOR EQUATIONS AND COEFFICIENTS:
      //    P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
      //     BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
      //    H.J. Liebe et al, JQSRT V.48, PP.629-643 (1992).
      //    M.J. Schwartz, Ph.D. thesis, M.I.T. (1997).
      //    SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
      //    This version differs from Liebe's MPM92 in two significant respects:
      //    1. It uses the modification of the 1- line width temperature dependence
      //    recommended by Schwartz: (1/T).
      //    2. It uses the same temperature dependence (X) for submillimeter 
      //    line widths as in the 60 GHz band: (1/T)**0.8 
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
      const char *version="PWR98";
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  PWR93O2AbsModel( xsec,
			   parameters[0], // continuum term scale factor
			   parameters[1], // line strength scale factor
			   parameters[2], // line broadening scale factor
			   parameters[3], // line coupling scale factor
			   model,
			   version, 
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  PWR93O2AbsModel( xsec,
			   0.00,
			   0.00,
			   0.00,
			   0.00,
			   model,
			   version,
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Full model " << name << " is running with \n"
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
	  out3 << "Full model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-MPM92"==name )
    {
      //   H. J. Liebe, P. W. Rosenkranz and G. A. Hufford,
      //   Atmospheric 60-GHz Oxygen Spectrum: New Laboratory 
      //   Measurements and Line Parameters
      //   JQSRT, Vol 48, pp. 629-643, 1992
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
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  MPM92O2AbsModel( xsec,
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  MPM92O2AbsModel( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-MPM89"==name )
    {
      //   H. J. Liebe,
      //   MPM - an atmospheric millimeter-wave propagation model,
      //   Int. J. Infrared and Mill. Waves, Vol 10, pp. 631-650, 1989.
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
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  MPM89O2AbsModel( xsec,
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  MPM89O2AbsModel( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-MPM87"==name )
    {
      //   H. J. Liebe and D. H. Layton,
      //   Millimeter-wave properties of the atmosphere: 
      //   Laboratory studies and propagation modelling,
      //   NITA Report 87-224, 
      //   U.S. Dept. of Commerce, National Telecommunications and Information
      //   Administration, Institute for Communication Sciences, rep. 87-224,
      //   325 Broadway, Boulder, CO 80303-3328
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
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  MPM87O2AbsModel( xsec,
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  MPM87O2AbsModel( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "O2-MPM85"==name )
    {
      //   H. J. Liebe and D. H. Layton,
      //   An updated model for millimeter wave propagation in moist air
      //   Radio Science, vol. 20, pp. 1069-1089, 1985
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
	  out3 << "Full model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  MPM85O2AbsModel( xsec,
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
	  out3 << "Full model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  MPM85O2AbsModel( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
      else if ( (model == "MPM93Scale") && (parameters.nelem() == 1) ) // --------------------
	{
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  MPM93_N2_continuum( xsec,
			      parameters[0],
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
      else if ( (model == "MPM93Scale") && (parameters.nelem() != 1) ) // --------------------
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires 1 scaling input\n"
	     << "parameters for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
	  throw runtime_error(os.str());
	  return;
	}
      else if ( (model != "user") && (model != "MPM93Scale") && (parameters.nelem() == 0) ) // --
	{
	  out3 << "Continuum model " << name << " running with \n" 
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
      /* --------------------------------------------------------------------------
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: Continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
      ----------------------------------------------------------------------*/
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "N2-DryContATM01"==name )
    {
      // data information about this continuum: 
      // Pardo et al. model model (IEEE, Trans. Ant. Prop., 
      // Vol 49, No 12, pp. 1683-1694, 2001)
      //
      // specific continuum parameters and units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : continuum strength coefficient  [1/m * 1/(Hz*Pa)²]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]   for  N2
      //     h2ovmr        : [1]   for  H2O
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  Pardo_ATM_N2_dry_continuum( xsec,
				      parameters[0], // coefficient
				      model,
				      f_mono,
				      p_abs,
				      t_abs,
				      vmr,
				      h2o_abs );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  Pardo_ATM_N2_dry_continuum( xsec,
				      0.000, // coefficient
				      model,
				      f_mono,
				      p_abs,
				      t_abs,
				      vmr,
				      h2o_abs );
	} 
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: Continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  BF86_CIA_N2( xsec,
		       parameters[0], // scaling factor
		       model,
		       f_mono,
		       p_abs,
		       t_abs,
		       vmr   );
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  BF86_CIA_N2( xsec,
		       0.0,
		       model,
		       f_mono,
		       p_abs,
		       t_abs,
		       vmr   );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: Continuum model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters. " << "\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "N2-CIArotCKDMT100"==name )
    {
      // data information about this continuum: 
      // A. Borysow and L. Frommhold, The Astrophysical Journal,
      // Vol. 311, pp.1043-1057, 1986
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_CIArot_n2( xsec,
			    parameters[0], // scaling factor
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_CIArot_n2( xsec,
			    0.0,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }  
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "N2-CIAfunCKDMT100"==name )
    {
      // data information about this continuum: 
      // Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M. Hartmann,
      // Infrared collision-induced absorption by 
      // N2 near 4.3 microns for atmospheric applications: 
      // Measurements and emprirical modeling,
      // Appl. Optics, 35, 5911-5917, (1996)

      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_CIAfun_n2( xsec,
			    parameters[0], // scaling factor
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_CIAfun_n2( xsec,
			    0.0,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }  
  // ============= CO2 continuum ========================================================
  else if ( "CO2-CKD241"==name )
    {
      // data information about this continuum: 
      // CKDv2.4.1 model at http://www.rtweb.aer.com/continuum_frame.html
      // This continuum accounts for the far wings of the many COS lines/bands since
      // the line is used with a cutoff in the line shape with +/- 25 cm^-1.
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
      //     h2o_abs       : [1]
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_241_co2( xsec,
		      parameters[0], // abs. scaling
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_241_co2( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  else if ( "CO2-CKDMT100"==name )
    {
      // data information about this continuum: 
      // CKD model at http://www.rtweb.aer.com/continuum_frame.html
      // This continuum accounts for the far wings of the many COS lines/bands since
      // the line is used with a cutoff in the line shape with +/- 25 cm^-1.
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
      //     h2o_abs       : [1]
      //
      const int Nparam = 1;
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "Continuum model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  CKD_mt_co2( xsec,
		      parameters[0], // abs. scaling
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
	  out3 << "Continuum model " << name << " running with \n" 
               << "the parameters for model " << model << ".\n";
	  CKD_mt_co2( xsec,
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
             << "Please see the arts user guide chapter 3.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "Continuum model " << name << " is running with \n"
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
	  out3 << "Continuum model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 3.\n";
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
	  out3 << "MPM93 liquid water cloud absorption model " << name << " is running with \n"
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
	  out3 << "MPM93 liquid water cloud absorption model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 4.\n";
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
	  out3 << "MPM93 ice water cloud absorption model " << name << " is running with \n"
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
	  out3 << "MPM93 ice water cloud absorption model " << name << " running with \n" 
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
             << "Please see the arts user guide chapter 4.\n";
	  throw runtime_error(os.str());
	  return;
	}
    }
  // ============= rain extinction from MPM93 ===========================================
  else if ( "rain-MPM93"==name )
    {
      // Rain extinction parameterization from MPM93 model, described in
      //  H. J. Liebe, 
      //  "MPM - An Atmospheric Millimeter-Wave Propagation Model",
      //  Int. J. Infrared and Millimeter Waves, vol. 10(6),
      //  pp. 631-650, 1989
      // and based on 
      //  Olsen, R.L., D.V. Rogers, and D. B. Hodge, "The aR^b relation in the
      //  calculation of rain attenuation", IEEE Trans. Antennas Propagat.,
      // vol. AP-26, pp. 318-329, 1978.
      //
      // specific continuum parameters and units:
      //  OUTPUT 
      //     xsec          : [1/m],
      //  INPUT	
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     model         : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [mm/h]
      //
      // rain parameters:
      // rain rate                         range: 0-150 mm/h
      //
      // valid atmospheric condition:
      // temperature      : (preferably above 273 K...)
      // 
      const int Nparam = 3; 
      if ( (model == "user") && (parameters.nelem() == Nparam) ) // -------------------------
	{
	  out3 << "MPM93 rain extinction model " << name << " is running with \n"
	       << "user defined parameters according to model " << model << ".\n";
	  MPM93RainExt(xsec,
		       parameters[0],     // scaling factror
		       parameters[1],     // scaling factror
		       parameters[2],     // scaling factror
		       model,             // model option
		       f_mono,
		       p_abs,
		       t_abs,
		       vmr );
	}
      else if ( (model == "user") && (parameters.nelem() != Nparam) ) // --------------------
	{
	  ostringstream os;
	  os << "MPM93 rain extinction model  " << name << " requires \n"
             << Nparam << " input parameter for the model " << model << ",\n"
             << "but you specified " << parameters.nelem() << " parameters.\n";
	  throw runtime_error(os.str());
	  return;
	}
      else if ( (model != "user") && (parameters.nelem() == 0) ) // --------------------
	{
	  out3 << "MPM93 rain extinction model " << name << " running with \n" 
               << "the parameter for model " << model << ".\n";
	 MPM93RainExt(xsec,
		      0.000,       // scaling factror
		      0.000,       // scaling factror
		      0.000,       // scaling factror
		      model,       // model option
		      f_mono,
		      p_abs,
		      t_abs,
		      vmr );
	}
      else if ( (model != "user") && (parameters.nelem() != 0) ) // --------------------
	{
	  ostringstream os;
	  os << "ERROR: MPM93 rain extinction model " << name << " requires NO input\n"
	     << "parameters for the model " << model << ",\n"
	     << "but you specified " << parameters.nelem() << " parameters.\n"
	     << "This ambiguity can not be solved by arts.\n"
             << "Please see the arts user guide chapter 4.\n";
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
	  //	  cout << fullnam << "\n";

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
//
// #################################################################################
// ############################# f2c code implementation ###########################
// #################################################################################
//
//
// ------------------- begin of f2c.h file --------------------------------
//
/* f2c.h  --  Standard Fortran to C header file */
#ifndef F2C_INCLUDE
#define F2C_INCLUDE

typedef long int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
#ifdef INTEGER_STAR_8	/* Adjust for integer*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif

#define TRUE_ (1)
#define FALSE_ (0)

/* Extern is for use with -E */
#ifndef Extern
#define Extern extern
#endif

/* I/O stuff */

#ifdef f2c_i2
/* for -i2 */
typedef short flag;
typedef short ftnlen;
typedef short ftnint;
#else
typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};

typedef union Multitype Multitype;

/*typedef long int Long;*/	/* No longer used; formerly in Namelist */

struct Vardesc {	/* for Namelist */
	char *name;
	char *addr;
	ftnlen *dims;
	int  type;
	};
typedef struct Vardesc Vardesc;

struct Namelist {
	char *name;
	Vardesc **vars;
	int nvars;
	};
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

/* procedure parameter types for -A and -C++ */

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int /* Unknown procedure type */ (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef /* Complex */ VOID (*C_fp)(...);
typedef /* Double Complex */ VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef /* Character */ VOID (*H_fp)(...);
typedef /* Subroutine */ int (*S_fp)(...);
#else
typedef int /* Unknown procedure type */ (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef /* Complex */ VOID (*C_fp)();
typedef /* Double Complex */ VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef /* Character */ VOID (*H_fp)();
typedef /* Subroutine */ int (*S_fp)();
#endif
/* E_fp is for real functions when -R is not specified */
typedef VOID C_f;	/* complex function */
typedef VOID H_f;	/* character function */
typedef VOID Z_f;	/* double complex function */
typedef doublereal E_f;	/* real function with -R not specified */

/* undef any lower-case symbols that your C compiler predefines, e.g.: */

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif
#endif

// ------------------- end of f2c.h file --------------------------------


// ------------------ begin of Borysow N2N2 F77 code --------------------


/* n2n2tks.f -- translated by f2c (version 20010821).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* Common Block Declarations */

struct {
    double temp, fnumin, fnumax, dnu;
} blockin_;

#define blockin_1 blockin_

struct {
    double slit, dx, wnrmax3;
} app3a_;

#define app3a_1 app3a_

struct {
    int nsri, ns, nsriup;
} app3b_;

#define app3b_1 app3b_

struct {
    double rsilo[201];
} rsilo_;

#define rsilo_1 rsilo_

struct {
    int initb;
} bou43_;

#define bou43_1 bou43_

union {
    struct {
	double omeg[201], rsi[201], rsigg[201], alfa;
    } m_1;
    struct {
	double omeg[201], rsi[201], rsigg[201], beta;
    } m_2;
} bba_;

#define bba_1 (bba_.m_1)
#define bba_2 (bba_.m_2)

struct {
    int nsol;
} bbc_;

#define bbc_1 bbc_

struct {
    double g0bf, delbf, om0;
} bf_;

#define bf_1 bf_

struct like_1_ {
    int like;
    char lgas[5];
};

#define like_1 (*(struct like_1_ *) &like_)

struct {
    int ik1k0;
} k1k0_;

#define k1k0_1 k1k0_

struct {
    int ibound;
} bbb_;

#define bbb_1 bbb_

struct energ_1_ {
    double eb[246]	/* was [41][6] */;
    int niv[6];
};

#define energ_1 (*(struct energ_1_ *) &energ_)

struct {
    int nlines;
} dimer_;

#define dimer_1 dimer_

struct n2part_1_ {
    double q1, wn2[2], b01, d01;
    int jrange2;
};
struct n2part_2_ {
    double q, wn2[2], b0, d0;
    int jrange1;
};

#define n2part_1 (*(struct n2part_1_ *) &n2part_)
#define n2part_2 (*(struct n2part_2_ *) &n2part_)

union {
    struct {
	double rsi[401];
    } m_1;
    struct {
	double rsibb[401];
    } m_2;
} bl3_;

#define bl3_1 (bl3_.m_1)
#define bl3_2 (bl3_.m_2)

union {
    struct {
	int idelv, iv, ivp, idell, il, ilp;
    } m_1;
    struct {
	int ldelvi, ivi, ivip, ldelel, ll, llp;
    } m_2;
} bbbb_;

#define bbbb_1 (bbbb_.m_1)
#define bbbb_2 (bbbb_.m_2)

/* Initialized data */

struct {
    double e_1[246];
    int e_2[6];
    } energ_ = { -54.99996, -54.86228, -54.58697, -54.17413, -53.62391, 
	    -52.93648, -52.11211, -51.15108, -50.05374, -48.82049, -47.45179, 
	    -45.94815, -44.31014, -42.53841, -40.63365, -38.59665, -36.42824, 
	    -34.12937, -31.70105, -29.14439, -26.46061, -23.65103, -20.71709, 
	    -17.66041, -14.48271, -11.18593, -7.77221, -4.24393, -.60374, 
	    3.14531, 6.99978, 10.95566, 15.00818, 19.15136, 23.37787, 
	    27.67681, 32.03237, 36.42278, 40.83668, 45.29436, 49.79246, 
	    -31.89437, -31.77215, -31.52779, -31.16143, -30.67334, -30.06382, 
	    -29.33328, -28.48222, -27.51123, -26.42099, -25.21229, -23.88603, 
	    -22.44322, -20.88502, -19.21272, -17.42777, -15.53182, -13.52669, 
	    -11.41446, -9.1975, -6.87848, -4.46049, -1.94714, .65736, 3.34788,
	     6.11816, 8.95978, 11.8613, 14.80383, 17.75924, 20.71774, 
	    23.71589, 0., 0., 0., 0., 0., 0., 0., 0., 0., -16.05019, -15.9464,
	     -15.73896, -15.42815, -15.0144, -14.4983, -13.88057, -13.16213, 
	    -12.34407, -11.42771, -10.41455, -9.30639, -8.10531, -6.81376, 
	    -5.43459, -3.97121, -2.42768, -.80899, .87859, 2.62689, 4.42334, 
	    6.24733, 8.06983, 9.90464, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
	     0., 0., 0., 0., 0., 0., 0., -6.49343, -6.41131, -6.24732, 
	    -6.00202, -5.67623, -5.27111, -4.78813, -4.22919, -3.59665, 
	    -2.89345, -2.12325, -1.29074, -.40202, .5345, 1.50455, 2.48212, 
	    3.46665, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1.76583, -1.70887, 
	    -1.59552, -1.427, -1.20523, -.93302, -.61434, -.25504, .13641, 0.,
	     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    -.17133, -.14341, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 41, 32, 24, 17, 9, 2 }
	    ;

struct {
    double fill_1[1];
    double e_2[4];
    int fill_3[1];
    } n2part_ = { {0}, 2., 1., 1.98957, 5.8e-6 };

struct {
    int fill_1[1];
    char e_2[5];
    } like_ = { {0}, "N2N2" };


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static int cs__1 = 1;
static int cs__0 = 0;
static double c_b24 = 2.9723;
static double c_b25 = -.99569;
static double c_b26 = .09464;
static double c_b27 = 1.2962e-12;
static double c_b28 = -.13048;
static double c_b29 = -.03128;
static double c_b30 = 3.7969e-14;
static double c_b31 = 1.03681;
static double c_b32 = -.14336;
static int cs__2 = 2;
static int cs__3 = 3;
static double c_b43 = .180926;
static double c_b44 = -1.69153;
static double c_b45 = .18605;
static double c_b46 = .3;
static double c_b47 = 0.;
static double c_b49 = 6.6017e-16;
static double c_b50 = 2.59982;
static double c_b51 = -.31831;
static double c_b52 = 1.2481e-12;
static double c_b53 = -.57028;
static double c_b54 = .05983;
static double c_b55 = 5.2681e-13;
static double c_b56 = -.24719;
static double c_b57 = .00519;
static double c_b58 = 2.7518e15;
static double c_b59 = -25.38969;
static double c_b60 = 2.46542;
static int cs__4 = 4;
static int cs__5 = 5;
static integer c__2 = 2;
static double c_b78 = .0825299;
static double c_b79 = -1.25562;
static double c_b80 = .12981;
static double c_b84 = 3.6611e-15;
static double c_b85 = 1.47688;
static double c_b86 = -.16537;
static double c_b87 = 6.1264e-10;
static double c_b88 = -2.25011;
static double c_b89 = .15289;
static double c_b90 = 7.982e-10;
static double c_b91 = -2.76152;
static double c_b92 = .21847;
static double c_b93 = 5.2868e-22;
static double c_b94 = 7.66253;
static double c_b95 = -.77527;
static double c_b112 = 119.261;
static double c_b113 = -3.78587;
static double c_b114 = .34024;
static double c_b115 = 9.3777e-12;
static double c_b116 = -.66548;
static double c_b117 = .0033;
static double c_b118 = 3.0395e-13;
static double c_b119 = .24728;
static double c_b120 = -.06607;
static double c_b183 = 1e-6;
static double c_b186 = 1.5;

#define temp (blockin_1.temp)
#define fnumin (blockin_1.fnumin)
#define fnumax (blockin_1.fnumax)
#define dnu (blockin_1.dnu)
#define slit (app3a_1.slit)
#define dx (app3a_1.dx)
#define rsilo (rsilo_1.rsilo)
#define omeg (bba_1.omeg)
#define rsi (bba_1.rsi)
#define rsigg (bba_1.rsigg)
#define nsol (bbc_1.nsol)
#define like (like_1.like)
#define ik1k0 (k1k0_1.ik1k0)
#define ibound (bbb_1.ibound)

/* TKS ****** SUBROUTINE N2N2TKS(T, F) ***************************************/
Numeric n2n2tks_(double t, double f)
{
    /* System generated locals */
    int s__1;
    double ret_val;

    /* Local variables */
    double hexa[10], quad[10], freq[10], e;
    int i__;
    double s, x, t1, t2, t3, t4;
    int ij, nf, jj;
    double rslow1, si;
    int nr;
    double ss[1], tt[2];
    extern /* Subroutine */ int bound32_(double *, double *, int 
	    *), bound54_(double *, double *, int *);
    double tksabs[5];
    extern /* Subroutine */ int spline_(int *, int *, int *, 
	    double *, double *, double *, double *, 
	    double *, double *, int *, double *);
    double dtrans[10], abscoef[10];
    extern /* Subroutine */ int addspec_(double *, double *, 
	    double *, double *, double *, double *, 
	    double *, int *, double *, double *, int *, 
	    int *, int *, int *, int *, int *);
    double eps, alfatot[10];
    extern /* Subroutine */ int partsum_(double *);

/*     ========================================= */
/*     Copyright (C) Aleksandra Borysow, 1987) */
/*     ==================================================================== */
/*     PROGRAM PREPARED BY ALEKSANDRA BORYSOW (APRIL'1987) */
/*     (UNIVERSITY OF TEXAS AT AUSTIN, PHYSICS DEPARTMENT) */
/*     ORIGINAL VERSION: WRITTEN ON CYBER */

/*     PROGRAM GENERATES N2-N2 COLLISION-INDUCED SPECTRA AT */
/*     TEMPERATURES BETWEEN 50 TO 300 K. */
/*     CIA SPECTRA MODELED AFTER PAPER (*) */
/*     ALEKSANDRA BORYSOW AND LOTHAR FROMMHOLD, */
/*     ASTROPHYSICAL JOURNAL, VOL. 311, PAGES 1043-1057, (1986) */

/*     REVISED BY GLENN ORTON (1989) - TO WORK ON SUN WORKSTATIONS */
/*     AND ON THE VAX MACHINES (FORTRAN-77) */
/*     PASSES STANDARD TEST ON SUN, AT 200K (JULY 1992) */
/*     ==================================================================== */

/*     ALSO IN REVISION: DOUBLE PRECISION FOR ALL F.P. VARIABLES */

/*     ==================================================================== */

/*     HISTORY: */

/*     2001-02-28 THOMAS KUHN: */
/*     CHANGE OF LINES */
/*       RSILO(I)=DLOG(RSI(I)*1.E80) */
/*     TO */
/*       RSILO(I)=(DLOG(RSI(I))+80.0D0*DLOG(10.0D0)) */
/*     BECAUSE OF OVERFLOW PROBLEMS. */
/*     COSMETICS FOR THE CODE TO BE FASTER READABLE. */
/*     CHANGE OF OUTPUT FILE NAME. */
/*     CHANGE OF OUTPUT FILE CONTENT */

/*     ==================================================================== */

/* TKS*      IMPLICIT REAL*8 (A-H,O-Z) */

/* TKS  INPUT/OUTPUT VARIABLES */
/*      REAL T, F */

    ret_val = 0.;

/*     TEMP   = TEMPERATURE IN KELVIN, SHOULD BE BETWEEN 50. AND 300. */
/*     FNUMIN = LOWEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU) */
/*     FNUMAX = HIGHEST FREQUENCY IN CM-1, FOR LISTING OF ALPHA(FNU) */
/*     LINE SHAPE MODELLING WILL BE MOST ACCURATE WITHIN RANGE OF */
/*     R-T SPECTRAL INTENSITIES AS 1:100. */
/*     DNU    = FREQUENCY INCREMENT IN CM-1. DNU SHOULD BE CHOSEN SO */
/*              THAT NOT MORE THAN 10 STEPS ARE NEEDED TO GO FROM */
/*              FNUMIN TO FNUMAX (ELSE ARRAY DIMENSIONS OF FREQ,ABSCOEF */
/*              MUST BE ADJUSTED IN ADDEM). */


/*      	USER: */
/*       ----- */
/* 	EDIT ONLY HERE: TEMP (K), MIN. FREQ. (CM^-1)= FNUMIN, */
/* 	MAX. FREQ. =  FNUMAX, STEP = DNU, SLITWIDTH (CM^-1)=SLIT */
/* 	(SLIT=4.3 IS EQUIVALENT TO THAT OF VOYAGER SPECTRA, ONLY BOUND BOUND */
/* 	SPECTRA ARE CONVOLUTED WITH THIS SLITWIDTH, THE FREE FREE SPECTRA */
/* 	ARE FAR TOO BROAD FOR THE SLITWIDTH FUNCTION TO MATTER, */
/* 	LEAVE LIKE = 1 (FOR LIKE PAIRS,  AS N2-N2) */
/* 	THE PROGRAM WILL ASSUME EQUILIBRIUM N2, */
/* 	ALLOWED TEMPERATURE RANGE: 50-300K (DO NOT  EXTEND IT BEYOND THESE LIMITS!) */
/* 	IF QUESTIONS: CONTACT ABORYSOW@NBI.DK */
/* TKS*      NF=INT((FNUMAX-FNUMIN)/DNU+0.5)+1 */
/* TKS*      IF (NF.GT.10) NF=10 */
/* TKS*      FNUMAX=FNUMIN+FLOAT(NF-1)*DNU */

/* TKS  INPUT TEMPERATURE (K) CHECK OF RANGE */
    if (t < 50. || t > 300.) {
      ostringstream os;
      os  << "out of T range ( 50<T<300)! return without calc.!" <<"\n";
      throw runtime_error(os.str());
      goto L999;
    }
    temp = t;

/*     *********************** INPUT DATA FROM USER *********************** */
/*     FNUMIN = MINIMUM FREQENCY [CM^-1] */
    fnumin = f / 29979245800.;
/*     FNUMAX = MAXIMUM FREQENCY [CM^-1] */
    fnumax = fnumin;
/*     ONLY ONE FREQUENCY PER CALL */
    nf = 1;
/*     DEFAULT VALUE OF FREQUENCY STEP [CM^-1] */
    dnu = 10.;
/*     DEFAULT VALUE */
    like = 1;
/*     SLIT = SLITWIDTH [CM^-1] */
/*     SLIT=4.3 IS EQUIVALENT TO THAT OF VOYAGER SPECTRA, ONLY BOUND BOUND */
/*     SPECTRA ARE CONVOLUTED WITH THIS SLITWIDTH, THE FREE FREE SPECTRA */
/*     ARE FAR TOO BROAD FOR THE SLITWIDTH FUNCTION TO MATTER. */
    slit = 4.3;
/*     ******************************************************************** */

/* TKS*      WRITE (6,14) LGAS,TEMP,FNUMIN,FNUMAX,DNU,NF-1 */
/* TKS*14    FORMAT(' ABSORPTION SPECTRA OF ',A5,' AT',F8.1,' K'/ */
/* TKS*     $   1X,43(1H=),/ */
/* TKS*     1' MIN.FREQ.=',F8.1,' CM-1',10X,'MAX.FREQ.=',F8.1,' CM-1',10X, */
/* TKS*     2'FREQ.INCREMENT=',F8.2,' CM-1',5X,'IN',I5,' STEPS'//) */


    partsum_(&temp);


/*     THE N2-N2 SPECTRA   FOR 50-300K */
/*     ================================ */

    x = log(temp);
    s__1 = nf;
    for (i__ = 1; i__ <= s__1; ++i__) {
/*         FREQ(I)=FNUMIN+FLOAT(I-1)*DNU */
	freq[i__ - 1] = fnumin;
	alfatot[i__ - 1] = 0.;
/* L10: */
	abscoef[i__ - 1] = 0.;
    }


/*    ==================================================================== */


    jj = 1;
L442:
/* L1023: */
    ++jj;
    if (jj == 42) {
	goto L444;
    }
    goto L442;
/*     EB(JJ,IV) JJ-ROTATIONAL LEVEL "L", IV- VIBRATIONAL LEVEL "V"; */
L444:



/*     ==================================================================== */

/*     QUADRUPOLAR INDUCTION: (50-300K) LAMBDA1,LAMBDA2,LAMBDA,L=2023&0223 */
/*     ------------------------------------------------------------------- */

    eps = 1e-5;
    tt[0] = 10.;
    bound32_(&temp, rsi, &nsol);
    ij = 0;
    rslow1 = 0.;
    s__1 = nsol;
    for (i__ = 1; i__ <= s__1; ++i__) {
	++ij;
/*        MOD CAN BE ONLY 0 OR 1 OR 2 */
	if (ij % 3 == 0) {
	    rslow1 = 1.5e-60;
	}
	if (ij % 3 == 1) {
	    rslow1 = 1.7e-60;
	}
	if (ij % 3 == 2) {
	    rslow1 = 1.6e-60;
	}
	if (rsi[i__ - 1] < 1e-60) {
	    rsi[i__ - 1] = rslow1;
	}
/* TKS*      RSILO(I)=DLOG(RSI(I)*1.D80) */
	rsilo[i__ - 1] = log(rsi[i__ - 1]) + log(10.) * 80.;
	omeg[i__ - 1] = (double) (i__ - 1) * dx;
/* L88: */
    }

/* L9991: */

    spline_(&nsol, &cs__1, &cs__0, &eps, omeg, rsilo, tt, ss, &si, &nr, rsigg)
	    ;

    ik1k0 = 1;
    ibound = 1;
/*       B-C LINESHAPE HERE */
/*       THESE VALUES (S,T1,T2) REPLACE VALUES GIVEN IN PAPER (*): */
/*       PUBLISHED IN AN ERRATUM, ASTROPHYSICAL JOURNAL, VOL.320, P.437 */
/*       (1987) */
    s = c_b24 * exp((c_b26 * x + c_b25) * x);
    t1 = c_b27 * exp((c_b29 * x + c_b28) * x);
    t2 = c_b30 * exp((c_b32 * x + c_b31) * x);
    e = 0.;
    t3 = 0.;
    t4 = 0.;

    addspec_(&s, &e, &t1, &t2, &t3, &t4, &temp, &nf, freq, abscoef, &cs__0, &
	    like, &cs__2, &cs__0, &cs__2, &cs__3);
    s__1 = nf;
    for (i__ = 1; i__ <= s__1; ++i__) {
	quad[i__ - 1] = abscoef[i__ - 1];
/* L20: */
	alfatot[i__ - 1] = abscoef[i__ - 1] + alfatot[i__ - 1];
    }



/*     ==================================================================== */

/*     HEXADECAPOLE COMPONENTS: LAMBDA1,LAMBDA2,LAMBDA,L=4045&0445 */
/*     ----------------------------------------------------------- */

    bound54_(&temp, rsi, &nsol);
    ij = 0;
    s__1 = nsol;
    for (i__ = 1; i__ <= s__1; ++i__) {
	++ij;
/*     MOD CAN BE ONLY 0 OR 1 OR 2 */
	if (ij % 3 == 0) {
	    rslow1 = 1.5e-60;
	}
	if (ij % 3 == 1) {
	    rslow1 = 1.7e-60;
	}
	if (ij % 3 == 2) {
	    rslow1 = 1.6e-60;
	}
	if (rsi[i__ - 1] < 1e-60) {
	    rsi[i__ - 1] = rslow1;
	}
/* TKS         RSILO(I)=DLOG(RSI(I)*1.E80) */
	rsilo[i__ - 1] = log(rsi[i__ - 1]) + log(10.) * 80.;
/* L111: */
	omeg[i__ - 1] = (double) (i__ - 1) * dx;
    }
    spline_(&nsol, &cs__1, &cs__0, &eps, omeg, rsilo, tt, ss, &si, &nr, rsigg)
	    ;

/*     --------------------------- */
/*       TEMPERATURES 50-140K */
/*     --------------------------- */

    if (temp >= 140.) {
	goto L333;
    }

    s = c_b43 * exp((c_b45 * x + c_b44) * x);
    e = c_b46 * exp((c_b47 * x + c_b47) * x);
    t1 = c_b49 * exp((c_b51 * x + c_b50) * x);
    t2 = c_b52 * exp((c_b54 * x + c_b53) * x);
    t3 = c_b55 * exp((c_b57 * x + c_b56) * x);
    t4 = c_b58 * exp((c_b60 * x + c_b59) * x);

    ik1k0 = 0;
    ibound = 1;
    addspec_(&s, &e, &t1, &t2, &t3, &t4, &temp, &nf, freq, abscoef, &cs__0, &
	    like, &cs__4, &cs__0, &cs__4, &cs__5);
    s__1 = nf;
    for (i__ = 1; i__ <= s__1; ++i__) {
	hexa[i__ - 1] = abscoef[i__ - 1];
	/*
	  s_wsle(&io___25);
	  do_lio(&c__9, &c__1, " T=50-140K: HEXA(", (ftnlen)17);
	  do_lio(&c__2, &c__1, (char *)&i__, (ftnlen)sizeof(int));
	  do_lio(&c__9, &c__1, ") =", (ftnlen)3);
	  do_lio(&c__5, &c__1, (char *)&abscoef[i__ - 1], (ftnlen)sizeof(
	  double));
	  e_wsle();
	*/
/* L50: */
	alfatot[i__ - 1] += abscoef[i__ - 1];
    }
    goto L334;

/*     --------------------------- */
/*       TEMPERATURES 140-300K */
/*     --------------------------- */

L333:
    ik1k0 = 0;
    ibound = 1;
    s = c_b78 * exp((c_b80 * x + c_b79) * x);
    e = c_b46 * exp((c_b47 * x + c_b47) * x);
    t1 = c_b84 * exp((c_b86 * x + c_b85) * x);
    t2 = c_b87 * exp((c_b89 * x + c_b88) * x);
    t3 = c_b90 * exp((c_b92 * x + c_b91) * x);
    t4 = c_b93 * exp((c_b95 * x + c_b94) * x);

    addspec_(&s, &e, &t1, &t2, &t3, &t4, &temp, &nf, freq, abscoef, &cs__0, &
	    like, &cs__4, &cs__0, &cs__4, &cs__5);
    s__1 = nf;
    for (i__ = 1; i__ <= s__1; ++i__) {
	hexa[i__ - 1] = abscoef[i__ - 1];
	/*
	  s_wsle(&io___26);
	  do_lio(&c__9, &c__1, " T=140-300K: HEXA(", (ftnlen)18);
	  do_lio(&c__2, &c__1, (char *)&i__, (ftnlen)sizeof(int));
	  do_lio(&c__9, &c__1, ") =", (ftnlen)3);
	  do_lio(&c__5, &c__1, (char *)&abscoef[i__ - 1], (ftnlen)sizeof(
	  double));
	  e_wsle();
	*/
/* L550: */
	alfatot[i__ - 1] += abscoef[i__ - 1];
    }

/*     ==================================================================== */

/*     DOUBLE TRANSITIONS: LAMBDA1,LAMBDA2,LAMBDA,L=2,2,3,3 */
/*     ---------------------------------------------------- */

/*     --------------------------- */
/*       TEMPERATURES 50-300K */
/*     --------------------------- */

L334:
    ik1k0 = 1;
    ibound = 0;
/* X        S=Y(X,1.19261D-58, -3.78587,0.34024) */
    s = c_b112 * exp((c_b114 * x + c_b113) * x);
    t1 = c_b115 * exp((c_b117 * x + c_b116) * x);
    t2 = c_b118 * exp((c_b120 * x + c_b119) * x);
    t3 = 0.;
    t4 = 0.;
    addspec_(&s, &e, &t1, &t2, &t3, &t4, &temp, &nf, freq, abscoef, &cs__0, &
	    like, &cs__2, &cs__2, &cs__3, &cs__3);
    s__1 = nf;
    for (i__ = 1; i__ <= s__1; ++i__) {
	dtrans[i__ - 1] = abscoef[i__ - 1];
/* L650: */
	alfatot[i__ - 1] += abscoef[i__ - 1];
    }

/*     ==================================================================== */

/*     ANISOTROPIC OVERLAP NEGLECTED (LAMBDA1,LAMBDA2,LAMBDA,L=0212) */
/*     SINCE THIS TERM IS EXTREMELY SMALL */

/*     ==================================================================== */

/* TKS*      PRINT 154, FREQ(1),FREQ(NF),FREQ(2)-FREQ(1),TEMP */
/* TKS*      PRINT 156, (ALFATOT(I),I=1,NF) */
/* TKS*154   FORMAT(///' ABSORPTION COEFFICIENT ALPHA(FNU), FROM ',F5.1, */
/* TKS*     $' CM-1 TO',F7.1,' CM-1, AT',F6.2,' CM-1 INCREMENTS, AT T=', */
/* TKS*     $F7.2,' K, IN UNITS OF CM-1 AMAGAT-2'/) */
/* TKS*156   FORMAT(' ',10D13.5) */
/* TKS*	OPEN(UNIT=10, FILE='OUT', STATUS='UNKNOWN') */
/* TKS*	WRITE(10, 2929) (FREQ(I), ALFATOT(I), I=1, NF) */
/* TKS*C	WRITE(10, 2929) (FREQ(I), QUAD(I), HEXA(I) ALFATOT(I), I=1, NF) */
/* TKS*2929	FORMAT(F10.3, E12.4) */

/* TKS*      STOP */



/* TKS FILL OUTPUT VARIABLE */
    tksabs[0] = quad[0];
    tksabs[1] = hexa[0];
    tksabs[2] = dtrans[0];
    tksabs[3] = alfatot[0];
    ret_val = alfatot[0];
/* TKS      print*,'QUAD(1),HEXA(1),DTRANS(1)=',QUAD(1),HEXA(1),DTRANS(1) */
L999:
    return ret_val;
} /* n2n2tks_ */

#undef temp
#undef fnumin
#undef fnumax
#undef dnu
#undef slit
#undef dx
#undef rsilo
#undef omeg
#undef rsi
#undef rsigg
#undef nsol
#undef like
#undef ik1k0
#undef ibound

#define wnrmax3 (app3a_1.wnrmax3)
#define rsilo (rsilo_1.rsilo)
#define omeg (bba_2.omeg)
#define rsigg (bba_2.rsigg)
#define beta (bba_2.beta)
#define nsol (bbc_1.nsol)
#define ibound (bbb_1.ibound)
#define q1 (n2part_1.q1)
#define wn2 (n2part_1.wn2)
#define b01 (n2part_1.b01)
#define d01 (n2part_1.d01)
#define jrange2 (n2part_1.jrange2)



/* ########################################################################## */


/* Subroutine */ int addspec_(double *g0, double *ep, double *
	tau1, double *tau2, double *tau5, double *tau6, 
	double *temp, int *nf, double *freq, double *abscoef,
	 int *mp, int *like, int *lambda1, int *lambda2, 
	int *lambda, int *lvalue)
{
    /* Initialized data */

    static double closchm = 2.68675484e19;
    static double boltzwn = .6950304;
    static double hbar = 1.054588757e-27;
    static double pi = 3.1415926535898;
    static double clight = 2.997925e10;

    /* Format strings */
    /*
    static char fmt_20[] = "(/\002 LAMBDA1,LAMBDA2, LAMBDA,LVALUE=\002,2i3,2"
	    "x,2i3,\002 COMPONENT.\002/15x,\002LINE SHAPE PARAMETERS:\002,6d1"
	    "2.3,5x,\002G(0)=\002,d12.3/)";
    */

    /* System generated locals */
    int s__1, s__2, s__3, s__4, s__5, s__6;
    double d__1, d__2;

    /* Builtin functions */
    /*
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    */

    /* Local variables */
    int list, jsum;
    extern double bgama_(double *, double *, double *, 
	    double *, double *, double *, double *);
    int i__, j;
    double calib, p;
    int i1, j1, i2, j2;
    double p1, p2, omega1, omega2;
    int ip, jp, iq;
    double twopic;
    int jplusl, ip1, jp1, ip2, jp2;
    double fac, cgs, xbg, wkf, frq, wki;
    extern double specfct_(double *, double *, double *, 
	    double *, int *, double *), clebsqr_(int *, 
	    int *, int *);
    double cg1s, cg2s;

    /* Fortran I/O blocks */
    /*
    static cilist io___38 = { 0, 6, 0, fmt_20, 0 };
    static cilist io___40 = { 0, 6, 0, 0, 0 };
    */


/*     THIS PROGRAM GENERATES LISTING OF R-T CIA  ALFA(OMEGA) */
/*     IF EITHER LAMBDA1 OR LAMBDA2 EQUAL TO ZERO - SINGLE TRANSITIONS; */
/*     DOUBLE TRANSITIONS ARE ASSUMED OTHERWISE. */
/*     LIKE=1 FOR LIKE SYSTEMS (AS H2-H2), SET LIKE=0 ELSE. */

/*     COMMON/BB/OMEG(201),RSI(201),RSIGG(201),NSOL,BETA */
/*     COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP */
/*     DIMENSION ABSCOEF(NF),FREQ(NF) */
    /* Parameter adjustments */
    --abscoef;
    --freq;

    /* Function Body */

    twopic = 2. * pi * clight;
    if (*like != 1) {
	*like = 0;
    }
/*        TAKE CARE OF FACTOR OF 1.E-60 HERE. */
/* Computing 2nd power */
    d__1 = pi;
/* Computing 2nd power */
    d__2 = closchm * 1e-30;
    calib = twopic * (d__1 * d__1 * 4. / (hbar * 3. * clight)) * (d__2 * d__2)
	    ;
    calib /= (double) (*like + 1);
    beta = 1. / (boltzwn * *temp);
    list = *nf;
    s__1 = list;
    for (i__ = 1; i__ <= s__1; ++i__) {
/* L88: */
	abscoef[i__] = 0.;
    }

/*     ROTATIONAL SPECTRUM FOR THE DETAILED LISTING   ******************* */
    /*
      s_wsfe(&io___38);
      do_fio(&c__1, (char *)&(*lambda1), (ftnlen)sizeof(int));
      do_fio(&c__1, (char *)&(*lambda2), (ftnlen)sizeof(int));
      do_fio(&c__1, (char *)&(*lambda), (ftnlen)sizeof(int));
      do_fio(&c__1, (char *)&(*lvalue), (ftnlen)sizeof(int));
      do_fio(&c__1, (char *)&(*g0), (ftnlen)sizeof(double));
      do_fio(&c__1, (char *)&(*ep), (ftnlen)sizeof(double));
      do_fio(&c__1, (char *)&(*tau1), (ftnlen)sizeof(double));
      do_fio(&c__1, (char *)&(*tau2), (ftnlen)sizeof(double));
      do_fio(&c__1, (char *)&(*tau5), (ftnlen)sizeof(double));
      do_fio(&c__1, (char *)&(*tau6), (ftnlen)sizeof(double));
      d__1 = *g0 * bgama_(&c_b47, tau1, tau2, ep, tau5, tau6, temp);
      do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(double));
      e_wsfe();
    */
    if (*lambda1 == 0 || *lambda2 == 0) {
	goto L152;
    }
    jplusl = jrange2 + max(*lambda1,*lambda2);

    /*
      s_wsle(&io___40);
      do_lio(&c__9, &c__1, "LAMBDA1,LAMBDA2,ABSCOEF(1)=", (ftnlen)27);
      do_lio(&c__2, &c__1, (char *)&(*lambda1), (ftnlen)sizeof(int));
      do_lio(&c__2, &c__1, (char *)&(*lambda2), (ftnlen)sizeof(int));
      do_lio(&c__5, &c__1, (char *)&abscoef[1], (ftnlen)sizeof(double));
      e_wsle();
    */
    jsum = 0;
    s__1 = jrange2;
    for (i1 = 1; i1 <= s__1; ++i1) {
	j1 = i1 - 1;
	s__2 = jplusl;
	for (ip1 = 1; ip1 <= s__2; ++ip1) {
	    jp1 = ip1 - 1;
	    cg1s = clebsqr_(&j1, lambda1, &jp1);
	    if (cg1s <= 0.) {
		goto L150;
	    } else {
		goto L130;
	    }
L130:
	    s__3 = j1 * (j1 + 1);
	    p1 = (double) (2 * j1 + 1) * wn2[1 + j1 % 2 - 1] * exp(
		    -1.4387859 / *temp * ((b01 - (double) s__3 * d01) * (
		    double) s__3)) / q1;
	    ++jsum;
	    s__3 = jp1 * ip1;
	    s__4 = j1 * i1;
	    omega1 = (b01 - (double) s__3 * d01) * (double) s__3 - (
		    b01 - (double) s__4 * d01) * (double) s__4;
	    s__3 = jrange2;
	    for (i2 = 1; i2 <= s__3; ++i2) {
		j2 = i2 - 1;
		s__4 = jplusl;
		for (ip2 = 1; ip2 <= s__4; ++ip2) {
		    jp2 = ip2 - 1;
		    cg2s = clebsqr_(&j2, lambda2, &jp2);
		    if (cg2s <= 0.) {
			goto L148;
		    } else {
			goto L132;
		    }
L132:
		    s__5 = j2 * (j2 + 1);
		    p2 = (double) (2 * j2 + 1) * wn2[1 + j2 % 2 - 1] * 
			    exp(-1.4387859 / *temp * ((b01 - (double) 
			    s__5 * d01) * (double) s__5)) / q1;
		    s__5 = jp2 * ip2;
		    s__6 = j2 * i2;
		    omega2 = (b01 - (double) s__5 * d01) * (double) 
			    s__5 - (b01 - (double) s__6 * d01) * (
			    double) s__6;
		    fac = calib * p1 * p2 * cg1s * cg2s;
		    s__5 = list;
		    for (i__ = 1; i__ <= s__5; ++i__) {
			frq = freq[i__] - omega1 - omega2;
			wki = freq[i__] * (1. - exp(-beta * freq[i__]));
			wkf = wki * fac;
			xbg = *g0 * bgama_(&frq, tau1, tau2, ep, tau5, tau6, 
				temp);
			if (ibound == 0) {
			    goto L555;
			}
			if (abs(frq) <= wnrmax3) {
			    xbg += specfct_(&frq, omeg, rsilo, rsigg, &nsol, &
				    beta);
			}
L555:
			abscoef[i__] += xbg * wkf;
/* L146: */
		    }
L148:
		    ;
		}
	    }
L150:
	    ;
	}
    }
    goto L2222;
/*     SINGLE TRANSITIONS AT NITROGEN'S ROTATIONAL FREQUENCIES */
/*     ======================================================= */
L152:
    jplusl = jrange2 + *lambda;
    s__2 = jrange2;
    for (i__ = 1; i__ <= s__2; ++i__) {
	j = i__ - 1;
	s__1 = jplusl;
	for (ip = 1; ip <= s__1; ++ip) {
	    jp = ip - 1;
	    cgs = clebsqr_(&j, lambda, &jp);
	    if (cgs <= 0.) {
		goto L200;
	    } else {
		goto L210;
	    }
L210:
	    s__4 = j * (j + 1);
	    p = (double) (2 * j + 1) * wn2[1 + j % 2 - 1] * exp(
		    -1.4387859 / *temp * ((b01 - (double) s__4 * d01) * (
		    double) s__4)) / q1;
	    ++jsum;
	    s__4 = jp * ip;
	    s__3 = j * i__;
	    omega1 = (b01 - (double) s__4 * d01) * (double) s__4 - (
		    b01 - (double) s__3 * d01) * (double) s__3;
	    fac = calib * p * cgs;
	    s__4 = list;
	    for (iq = 1; iq <= s__4; ++iq) {
		frq = freq[iq] - omega1;
/*               XWKI=FREQ(IQ)*(1.-EXP(-BETA*FREQ(IQ))) */
		wki = freq[iq] * (1. - exp(-beta * freq[iq]));
		wkf = wki * fac;
		xbg = *g0 * bgama_(&frq, tau1, tau2, ep, tau5, tau6, temp);
		if (ibound == 0) {
		    goto L444;
		}
		if (abs(frq) <= wnrmax3) {
		    xbg += specfct_(&frq, omeg, rsilo, rsigg, &nsol, &beta);
		}
L444:
		abscoef[iq] += xbg * wkf;
/* L199: */
	    }
L200:
	    ;
	}
    }

L2222:
/* TKS 2222   PRINT 44,(ABSCOEF(I),I=1,LIST) */
/* L44: */
    return 0;
} /* addspec_ */

#undef wnrmax3
#undef rsilo
#undef omeg
#undef rsigg
#undef beta
#undef nsol
#undef ibound
#undef q1
#undef wn2
#undef b01
#undef d01
#undef jrange2

#define q (n2part_2.q)
#define wn2 (n2part_2.wn2)
#define b0 (n2part_2.b0)
#define d0 (n2part_2.d0)
#define jrange1 (n2part_2.jrange1)

/* Subroutine */ int partsum_(double *temp)
{
    /* System generated locals */
    int s__1;

    /* Local variables */
    int j;
    double dq;

/*     N2 ROTATIONAL PARTITION SUM Q = Q(T). */



/*     Q,B0,D0,WN2 - PARTITION FCT., ROT.CONSTANTS, WEIGHTS FOR N2 */
    q = 0.;
    j = 0;
L50:
    s__1 = j * (j + 1);
    dq = (double) (2 * j + 1) * wn2[1 + j % 2 - 1] * exp(-1.4387859 * ((
	    b0 - (double) s__1 * d0) * (double) s__1) / *temp);
    q += dq;
    ++j;
    if (dq > q / 900.) {
	goto L50;
    }
    jrange1 = j;
/* TKS      PRINT 30, Q, JRANGE1 */

/* L30: */
    return 0;
} /* partsum_ */

#undef q
#undef wn2
#undef b0
#undef d0
#undef jrange1

#define slit (app3a_1.slit)
#define dx (app3a_1.dx)
#define wnrmax3 (app3a_1.wnrmax3)
#define nsri (app3b_1.nsri)
#define ns (app3b_1.ns)
#define nsriup (app3b_1.nsriup)
#define rsi (bl3_1.rsi)

/* Subroutine */ int profile_(double *x, double *y)
{
    /* System generated locals */
    int s__1;

    /* Local variables */
    int i__;
    double slope;
    int n1;
    double x0;
    int nc;
    double dr;
    int no;
    double xi;
    int nu;

/*     A TRIANGULAR SLIT FUNCTION IS USED. */


/*     COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP */

    if (*y < 0.) {
	goto L105;
    } else if (*y == 0) {
	goto L106;
    } else {
	goto L1;
    }
L1:
    x0 = nsri + 1. + *x / dx;
    nc = (int) x0;
    n1 = nc + 1;
    slope = *y / slit;
    nu = (int) (x0 - ns);
    if (nu < 1) {
	nu = 1;
    }
    if (nu > nsriup) {
	return 0;
    }
    no = (int) (x0 + ns);
    if (no > nsriup) {
	no = nsriup;
    }
    if (no < 1) {
	return 0;
    }
    if (nc > nsriup) {
	nc = nsriup;
    }
    if (nc <= 1) {
	goto L101;
    }
    s__1 = nc;
    for (i__ = nu; i__ <= s__1; ++i__) {
	xi = (i__ - 1.) * dx - wnrmax3;
	dr = slope * (xi - (*x - slit));
	if (dr <= 0.) {
	    goto L100;
	}
	rsi[i__ - 1] += dr;
L100:
	;
    }
L101:

    if (nc >= nsriup) {
	return 0;
    }
    if (n1 < 1) {
	n1 = 1;
    }
    s__1 = no;
    for (i__ = n1; i__ <= s__1; ++i__) {
	xi = (i__ - 1.) * dx - wnrmax3;
	dr = *y - slope * (xi - *x);
	if (dr <= 0.) {
	    goto L102;
	}
	rsi[i__ - 1] += dr;
L102:
	;
    }
    return 0;
L105:
/* TKS  105 PRINT 10,SLIT */
/* TKS   10 FORMAT(/' A TRIANGULAR SLIT FUNCTION OF',F6.3,' CM-1 HALFWIDTH IS */
/* TKS     ' USED'/) */
L106:
    return 0;
} /* profile_ */

#undef slit
#undef dx
#undef wnrmax3
#undef nsri
#undef ns
#undef nsriup
#undef rsi


/* X    FUNCTION SPECFCT(FREQ,OMEGA,PHI,PHI2,N,RTEMP) */

double specfct_(double *freq, double *omega, double *phi, 
	double *phi2, int *n, double *rtemp)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double tfac, f, gp, si;
    int nr;
    extern /* Subroutine */ int ixpolat_(int *, int *, int *, 
	    double *, double *, double *, double *, 
	    double *, double *, int *, double *);


/*     THIS INTERPOLATES THE SPECTRAL FUNCTION PHI(FREQ) DEFINED AT */
/*     OMEGA(N) AS PHI(N). PHI2 IS THE SECOND DERIVATIVE AT OMEGA */
/*     WHICH MUST BE OBTAINED FIRST (USE SPLINE FOR THAT PURPOSE). */
/*     RTEMP IS THE RECIPROCAL TEMPERATURE IN CM-1 UNITS. */
/*     NOTE THAT WE INTERPOLATE 1.E80 TIMES THE LOGARITHM OF PHI(OMEGA) */
/*     NOTE THAT IN GSO'S REVISION, THIS FACTOR IS REMOVED. */
/*     (REVISION MODIFIED) */



/*      PRINT*,'FREQ =',FREQ */
/*      PRINT*,'OMEGA=',OMEGA */
/*      PRINT*,'RTEMP=',RTEMP */
/*      PRINT*,'PHI  =',PHI */
/*      PRINT*,'PHI2 =',PHI2 */
    /* Parameter adjustments */
    --phi2;
    --phi;
    --omega;

    /* Function Body */
    tfac = 0.;
    f = *freq;
    if (f >= 0.) {
	goto L20;
    } else {
	goto L10;
    }
L10:
    f = abs(f);
    tfac = -(*rtemp) * f;
L20:
    if (f <= omega[*n]) {
	goto L30;
    }
    ret_val = exp(-(phi[*n - 1] - phi[*n]) * (f - omega[*n]) / (omega[*n] - 
	    omega[*n - 1]) + phi[*n] + tfac) * 1e-80;
/*      print*,' (A) SPECFCT=',SPECFCT */
/* X    SPECFCT=DEXP(-(PHI(N-1)-PHI(N))*(F-OMEGA(N))/ */
/* X   $(OMEGA(N)-OMEGA(N-1))+PHI(N)+TFAC) */
    return ret_val;
/* 30   PRINT*,'SI,NR,PHI2=', */
/*     &        SI,NR,PHI2 */
/*      CALL IXPOLAT(N,1,0,1.D-6,OMEGA,PHI,F,GP,SI,NR,PHI2) */
L30:
    ixpolat_(n, &cs__1, &cs__0, &c_b183, &omega[1], &phi[1], &f, &gp, &si, &
	    nr, &phi2[1]);
    ret_val = exp(tfac + gp) * 1e-80;
/* X    SPECFCT=DEXP(TFAC+GP) */
/*      print*,' (B) GP,SPECFCT=',GP,SPECFCT */

    return ret_val;
} /* specfct_ */
#define slit (app3a_1.slit)
#define dx (app3a_1.dx)
#define wnrmax3 (app3a_1.wnrmax3)
#define nsri (app3b_1.nsri)
#define ns (app3b_1.ns)
#define nsriup (app3b_1.nsriup)
#define eb (energ_1.eb)
#define niv (energ_1.niv)
#define nlines (dimer_1.nlines)
#define rsibb (bl3_2.rsibb)
#define ldelvi (bbbb_2.ldelvi)
#define ivi (bbbb_2.ivi)
#define ivip (bbbb_2.ivip)
#define ldelel (bbbb_2.ldelel)
#define ll (bbbb_2.ll)
#define llp (bbbb_2.llp)

/* Subroutine */ int bound32_(double *temp, double *rsi, int *
	nsol)
{
    /* Initialized data */

    static int ldelvis[63] = { 0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,
	    1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,
	    3,4,4,4,4,4,4,4,4 };
    static int ivis[63] = { 0,0,1,1,2,2,3,3,4,4,0,0,0,0,1,1,1,1,2,2,2,2,
	    3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,1,2,0,1,2,0,1,0,1,0,1,2,
	    0,1,0,1,0,1,0,1 };
    static int ivips[63] = { 0,0,1,1,2,2,3,3,4,4,1,1,1,1,2,2,2,2,3,3,3,3,
	    4,4,4,4,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,3,4,5,3,4,5,3,4,3,4,3,4,5,
	    4,5,4,5,4,5,4,5 };
    static int ldelels[63] = { 1,3,1,3,1,3,1,3,1,3,-3,-1,1,3,-3,-1,1,3,
	    -3,-1,1,3,3,1,-1,-3,-3,-1,1,3,-3,-1,1,3,-3,-1,1,3,-1,-3,1,3,1,1,1,
	    -1,-1,-1,-3,-3,3,3,-3,-3,-3,-3,-3,3,3,1,1,-1,-1 };
    static double as[63] = { 4.4844e-40,4.4356e-40,2.9345e-40,2.885e-40,
	    1.6441e-40,1.5899e-40,7.2882e-41,6.7748e-41,1.0378e-41,1.3041e-42,
	    1.5006e-41,1.537e-41,1.6139e-41,1.7143e-41,1.9985e-41,2.0169e-41,
	    2.0994e-41,2.2094e-41,1.636e-41,1.6281e-41,1.6714e-41,1.7326e-41,
	    8.0425e-42,8.0862e-42,8.0093e-42,8.1366e-42,2.4471e-42,2.5406e-42,
	    2.6629e-42,2.8064e-42,4.6227e-42,4.715e-42,4.8513e-42,5.0133e-42,
	    3.9968e-42,3.984e-42,3.981e-42,3.9687e-42,1.1806e-42,1.3458e-42,
	    3.8746e-42,3.9219e-42,7.3334e-43,1.339e-42,1.3041e-42,7.1401e-43,
	    1.3461e-42,6.5776e-43,6.9847e-43,1.3517e-42,7.5545e-43,1.3268e-42,
	    6.9847e-43,1.3517e-42,7.464e-43,2.1322e-43,2.6037e-43,2.0823e-43,
	    2.0632e-43,2.1067e-43,2.0531e-43,2.1218e-43,2.3006e-43 };
    static double bs[63] = { 4.3e-4,4.6e-4,8.3e-4,8.9e-4,.0017,.00186,
	    .0041,.00457,0.,0.,9.99e-4,5.23e-4,1.49e-4,-1.68e-4,.001837,
	    .001153,6.6e-4,2.54e-4,.003603,.002677,.002101,.001738,.00595,
	    .006843,0.,.007035,.001025,6.42e-4,2.54e-4,-1.64e-4,.002342,
	    .001975,.00164,.001328,.004943,.004999,.005461,.006839,0.,.010993,
	    0.,0.,.001367,.005262,0.,.001601,.00451,0.,.001828,.004175,.04816,
	    .007033,.001828,.004175,.009338,.003733,.008508,.006979,0.,
	    .005035,0.,.004169,0. };
    static double twopic = 1.88365183e11;
    static double pi = 3.141592654;

    /* System generated locals */
    int s__1;
    double d__1;

    /* Local variables */
    double alfa;
    int nnii, nsol2;
    double a, b;
    int i__, l, n;
    double stoke, stoki, am, pf;
    int lp;
    double rm;
    int nr, iv;
    double stokip;
    int ivp;
    extern double clebsqr_(int *, int *, int *);
    extern /* Subroutine */ int profile_(double *, double *);


#define eb_ref(a_1,a_2) eb[(a_2)*41 + a_1 - 42]



/*     COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP */


/*          STORED VALUES */
    /* Parameter adjustments */
    --rsi;

    /* Function Body */


/*     EB(I,K) - BOUND ENERGIES */
/*     AM =  (MTX.EL. (L,BETA,L') )**2 */
/*     M(L,L',V,V') OF PAPER (*) ARE TO BE CORRECTED: */
/*     M(L,L',V,V') = AM * (2L+1) * C(L,3,L')**2 */
/*     NSRI - HOW MANY POINTS FOR B-B SPECTRAL FUNCTION TO BE GIVEN */
/*     WNRMAX3 - THE FREQUENCY RANGE OF B-B CONTRIBUTION */
/*     SLIT- THE HALFWIDTH OF THE SPECTRAL PROFILE CONVOLUTED WITH */
/*     B-B SPECTRUM, IN [CM-1]. */
/*       A,B, COEFFICIENTS , EQ. 7, A.BORYSOW, L.FROMMHOLD, */
/*       AP. J. VOL.311, 1043-1057, (1986) */

    nsri = 190;
    wnrmax3 = 45.;
    nsriup = (nsri << 1) + 1;
    dx = wnrmax3 / (double) nsri;
    ns = (int) (slit / dx);

    for (i__ = 1; i__ <= 401; ++i__) {
/* L300: */
	rsibb[i__ - 1] = 0.;
    }

    alfa = 1. / (*temp * .69519);
    rm = 2.32498211e-23;
/*     RM - REDUCED MASS FOR N2-N2 */

    d__1 = rm * 1.380662e-16 * *temp * 2. * pi / 4.3906208382975998e-53;
    pf = pow(d__1, c_b186);
/*       SCALE DOWN BY 1.D60 */
    pf *= 1e-60;

    nr = 0;
L555:
    ++nr;
    ldelvi = ldelvis[nr - 1];
    ivi = ivis[nr - 1];
    ivip = ivips[nr - 1];
    ldelel = ldelels[nr - 1];
    a = as[nr - 1];
    b = bs[nr - 1];
/*     WRITE (6,334) LDELVI,IVI,IVIP, LDELEL,A,B */
/* L334: */

/*       LDELVI=DELTA(V)=V'-V */
/*       IVI = V; IVIP = V' */
/*       LDELEL = DELTA(L) = L'-L */

    iv = ivi + 1;
    ivp = ivip + 1;
/*       THESE ARE ENERGY COLUMNS (V, V') */

    nnii = niv[iv - 1];

    s__1 = nnii;
    for (l = 1; l <= s__1; ++l) {
/*       LOOP OVER INITIAL L-VALUES... */

	am = a * exp(-b * (double) ((l - 1) * l));
/*       L - NUMBER OF ROW, (L-1) - ROTATIONAL LEVEL */
	lp = l + ldelel;
	ll = l - 1;
	llp = lp - 1;
/*       LL,LLP ARE INITIAL L AND FINAL L' ANGULAR MOMENTUM QUANTUM NUMBERS */

	if (lp > niv[ivp - 1] || lp < 1) {
	    goto L20;
	}
	if (eb_ref(lp, ivp) == 0.) {
	    goto L20;
	}
	if (eb_ref(l, iv) == 0.) {
	    goto L20;
	}
	stoke = eb_ref(lp, ivp) - eb_ref(l, iv);

	stoki = am * exp(-alfa * eb_ref(l, iv)) / pf * (double) ((ll << 1)
		 + 1) * clebsqr_(&ll, &cs__3, &llp);

	profile_(&stoke, &stoki);
	if (stoki > 0.) {
	    ++nlines;
	}

	stokip = am * exp(-alfa * eb_ref(lp, ivp)) / pf * (double) ((llp 
		<< 1) + 1) * clebsqr_(&llp, &cs__3, &ll);

	d__1 = -stoke;
	profile_(&d__1, &stokip);
	if (stokip > 0.) {
	    ++nlines;
	}

L20:
	;
    }
    if (nr == 63) {
	goto L56;
    }
    goto L555;
L56:

/*     32 ENTRIES FOR (3220)+(3202) IN TABLE 6, 63 IN ALL: */
/*     DATA EXPANDED AND INCLUDE NOW ALL POSSIBLE B-B TRANSITIONS */

    s__1 = nsriup;
    for (n = 1; n <= s__1; ++n) {
/* L90: */
	rsibb[n - 1] = rsibb[n - 1] / twopic / slit;
    }

    *nsol = nsri + 1;
    nsol2 = *nsol + 1;
/*     RSI - CONTRIBUTION FOR POSITIVE FREQUENCY SHIFTS */

    s__1 = *nsol;
    for (i__ = 1; i__ <= s__1; ++i__) {
/* L22: */
	rsi[i__] = rsibb[*nsol - 1 + i__ - 1];
    }

/*        PRINT 999, (RSI(I),I=1,NSOL) */
/* L999: */

    return 0;
} /* bound32_ */

#undef slit
#undef dx
#undef wnrmax3
#undef nsri
#undef ns
#undef nsriup
#undef eb
#undef niv
#undef nlines
#undef rsibb
#undef ldelvi
#undef ivi
#undef ivip
#undef ldelel
#undef ll
#undef llp


#undef eb_ref

#define slit (app3a_1.slit)
#define dx (app3a_1.dx)
#define wnrmax3 (app3a_1.wnrmax3)
#define nsri (app3b_1.nsri)
#define ns (app3b_1.ns)
#define nsriup (app3b_1.nsriup)
#define eb (energ_1.eb)
#define niv (energ_1.niv)
#define nlines (dimer_1.nlines)
#define rsibb (bl3_2.rsibb)

/* Subroutine */ int bound54_(double *temp, double *rsi, int *
	nsol)
{
    /* Initialized data */

    static int ldelvis[54] = { 0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,
	    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 
	    };
    static int ivis[54] = { 0,0,0,1,1,1,2,2,2,3,3,3,0,0,0,0,0,0,1,1,1,1,
	    1,1,2,2,2,2,2,2,3,3,3,3,3,3,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2 };
    static int ivips[54] = { 0,0,0,1,1,1,2,2,2,3,3,3,1,1,1,1,1,1,2,2,2,2,
	    2,2,3,3,3,3,3,3,4,4,4,4,4,4,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4 };
    static int ldelels[54] = { 1,3,5,1,3,5,1,3,5,1,3,5,-5,-3,-1,1,3,5,-5,
	    -3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,1,3,5,-5,-3,-1,
	    1,3,5,-5,-3,-1,1,3,5 };
    static double as[54] = { 7.9332e-42,7.8207e-42,7.7235e-42,4.5815e-42,
	    4.4834e-42,4.4059e-42,2.173e-42,2.0824e-42,2.025e-42,7.7222e-43,
	    7.0351e-43,6.6815e-43,4.9611e-43,5.2232e-43,5.2979e-43,5.4652e-43,
	    5.6827e-43,5.9277e-43,5.733e-43,6.062e-43,6.0862e-43,6.2104e-43,
	    6.3809e-43,6.5698e-43,3.9501e-43,4.1599e-43,4.1033e-43,4.1097e-43,
	    4.1339e-43,4.153e-43,1.5858e-43,1.5976e-43,1.5478e-43,1.5066e-43,
	    1.4554e-43,1.3848e-43,9.9241e-44,1.0109e-43,1.0396e-43,1.0758e-43,
	    1.1176e-43,1.1636e-43,1.646e-43,1.647e-43,1.6617e-43,1.6837e-43,
	    1.7085e-43,1.7327e-43,1.1797e-43,1.1593e-43,1.1405e-43,1.1174e-43,
	    1.0853e-43,1.0401e-43 };
    static double bs[54] = { 6.12e-4,6.35e-4,6.77e-4,.001137,.001201,
	    .001341,.00229,.002449,.00287,.005426,.005876,.00745,.001,8.83e-4,
	    6.09e-4,3.92e-4,2.07e-4,3.7e-5,.001625,.001624,.001305,.001084,
	    9.27e-4,8.21e-4,.002978,.003273,.002994,.002954,.003153,.003668,
	    .005799,.006423,.006733,.00796,.010937,.019179,.001229,9.93e-4,
	    7.67e-4,5.43e-4,3.09e-4,5.1e-5,.002456,.0023,.00221,.002193,
	    .002273,.002506,.004556,.004825,.005454,.006725,.009431,.016672 };
    static double twopic = 1.88365183e11;
    static double pi = 3.141592654;

    /* System generated locals */
    int s__1;
    double d__1;

    /* Local variables */
    double alfa;
    int nnii, ivip, nsol2;
    double a, b;
    int i__, l, n;
    double stoke, stoki, am, pf;
    int ll, lp;
    double rm;
    int nr, ldelel, iv, ldelvi;
    double stokip;
    int ivi, llp, ivp;
    extern double clebsqr_(int *, int *, int *);
    extern /* Subroutine */ int profile_(double *, double *);


#define eb_ref(a_1,a_2) eb[(a_2)*41 + a_1 - 42]



/*     COMMON/APP3/SLIT,DX,NSRI,WNRMAX3,NS,NSRIUP */
/*          STORED VALUES */

/* TKS     THIS DATA STRUCTURE HAS 55 ENTRIES AND NOT 54!! */
/* TKS      DATA LDELVIS /0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1, */
/* TKS     &              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2, */
/* TKS     &              2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/ */
/* TKS  WE HAVE SSKIPPED THE LAST "2" IN THIS DATA STATEMENT */
    /* Parameter adjustments */
    --rsi;

    /* Function Body */








    nsri = 190;
    wnrmax3 = 47.;
    nsriup = (nsri << 1) + 1;
    dx = wnrmax3 / (double) nsri;
    ns = (int) (slit / dx);

    for (i__ = 1; i__ <= 401; ++i__) {
/* L300: */
	rsibb[i__ - 1] = 0.;
    }
    alfa = 1. / (*temp * .69519);
    rm = 2.32498211e-23;
/*     RM - REDUCED MASS FOR N2-N2 */

    d__1 = rm * 1.380662e-16 * *temp * 2. * pi / 4.3906208382975998e-53;
    pf = pow(d__1, c_b186);
/*       SCALE DOWNWARD BY 1.D60 */
    pf *= 1e-60;

    nr = 0;
L555:
    ++nr;
    ldelvi = ldelvis[nr - 1];
    ivi = ivis[nr - 1];
    ivip = ivips[nr - 1];
    ldelel = ldelels[nr - 1];
    a = as[nr - 1];
    b = bs[nr - 1];
/*       WRITE (6,334) LDELVI,IVI,IVIP, LDELEL,A,B */
/* L334: */
/*       LDELVI=DELTA(V)=V'-V */
/*       IVI = V */
/*       LDELEL = DELTA(L) = L'-L */

    iv = ivi + 1;
    ivp = ivip + 1;
    nnii = niv[iv - 1];
    s__1 = nnii;
    for (l = 1; l <= s__1; ++l) {
	am = a * exp(-b * (double) (l * (l + 1)));
	lp = l + ldelel;
	ll = l - 1;
	llp = lp - 1;
	if (lp > niv[ivp - 1] || lp < 1) {
	    goto L20;
	}
	if (eb_ref(lp, ivp) == 0.) {
	    goto L20;
	}
	if (eb_ref(l, iv) == 0.) {
	    goto L20;
	}

	stoke = eb_ref(lp, ivp) - eb_ref(l, iv);
	stoki = am * exp(-alfa * eb_ref(l, iv)) / pf * (double) ((ll << 1)
		 + 1) * clebsqr_(&ll, &cs__5, &llp);
	profile_(&stoke, &stoki);
	if (stoki > 0.) {
	    ++nlines;
	}
	stokip = am * exp(-alfa * eb_ref(lp, ivp)) / pf * (double) ((llp 
		<< 1) + 1) * clebsqr_(&llp, &cs__5, &ll);
	d__1 = -stoke;
	profile_(&d__1, &stokip);
	if (stokip > 0.) {
	    ++nlines;
	}
L20:
	;
    }
    if (nr == 54) {
	goto L56;
    }

    goto L555;
L56:
/*       54 ENTRIES FOR (5440)=(5404) IN TABLE 6 */

    s__1 = nsriup;
    for (n = 1; n <= s__1; ++n) {
/* L90: */
	rsibb[n - 1] = rsibb[n - 1] / twopic / slit;
    }

    *nsol = nsri + 1;
    nsol2 = *nsol + 1;

    s__1 = *nsol;
    for (i__ = 1; i__ <= s__1; ++i__) {
/* L22: */
	rsi[i__] = rsibb[*nsol - 1 + i__ - 1];
    }

/*     PRINT 999, (RSI(I),I=1,NSOL) */
/* L999: */

    return 0;
} /* bound54_ */

#undef slit
#undef dx
#undef wnrmax3
#undef nsri
#undef ns
#undef nsriup
#undef eb
#undef niv
#undef nlines
#undef rsibb


#undef eb_ref


double clebsqr_0_(int n__, int *l, int *lambda, int *lp)
{
    /* System generated locals */
    int s__1, s__2;
    double ret_val, d__1;

    /* Local variables */
    extern double fctl_(int *);
    double f;
    int i__;
    double p;
    int i0, i1;
    double fc;

/*     SQUARE OF CLEBSCH-GORDAN COEFFICIENT (L,LAMBDA,0,0;LP,0) */
/*     FOR INTEGER ARGUMENTS ONLY */
/*     NOTE THAT LAMBDA SHOULD BE SMALL, MAYBE @10 OR SO. */


    switch(n__) {
	case 1: goto L_threej2;
	}

    fc = (double) ((*lp << 1) + 1);
    goto L2;


L_threej2:
/*     THIS ENTRY RETURNS THE SQUARED 3-J SYMBOL   L LAMBDA LP */
/*                                                 0    0    0 */
/*     INSTEAD OF THE CLEBSCH-GORDAN COEFFICIENT */
/*     (LIMITATION TO INTEGER ARGUMENTS ONLY) */

/*     NOTE THAT THE THREE-J SYMBOLS ARE COMPLETELY SYMMETRIC IN THE */
/*     ARGUMENTS. IT WOULD BE ADVANTAGEOUS TO REORDER THE INPUT ARGUMENT */
/*     LIST SO THAT LAMBDA BECOMES THE SMALLEST OF THE 3 ARGUMENTS. */
    fc = 1.;
L2:
    ret_val = 0.;
    if (*l + *lambda < *lp || *lambda + *lp < *l || *l + *lp < *lambda) {
	return ret_val;
    }
    if ((*l + *lp + *lambda) % 2 != 0) {
	return ret_val;
    }
    if (*l < 0 || *lp < 0 || *lambda < 0) {
	return ret_val;
    }
    f = 1. / (double) (*l + *lp + 1 - *lambda);
    if (*lambda == 0) {
	goto L22;
    }
    i1 = (*l + *lp + *lambda) / 2;
    i0 = (*l + *lp - *lambda) / 2 + 1;
    s__1 = i1;
    for (i__ = i0; i__ <= s__1; ++i__) {
/* L20: */
	f = f * (double) i__ / (double) ((i__ << 1) + 1 << 1);
    }
L22:
    s__1 = *lambda + *l - *lp;
    s__2 = *lambda + *lp - *l;
    p = fc * f * fctl_(&s__1) * fctl_(&s__2);
    s__1 = (*lambda + *l - *lp) / 2;
    s__2 = (*lambda + *lp - *l) / 2;
/* Computing 2nd power */
    d__1 = fctl_(&s__1) * fctl_(&s__2);
    ret_val = p / (d__1 * d__1);
    return ret_val;
} /* clebsqr_ */

double clebsqr_(int *l, int *lambda, int *lp)
{
    return clebsqr_0_(0, l, lambda, lp);
    }

double threej2_(void)
{
    return clebsqr_0_(1, (int *)0, (int *)0, (int *)0);
    }

double fctl_(int *n)
{
    /* System generated locals */
    int s__1;
    double ret_val, d__1;

    /* Local variables */
    int i__, j;
    double z__;



    ret_val = 1.;
    if (*n <= 1) {
	return ret_val;
    }
    if (*n > 15) {
	goto L20;
    }
    j = 1;
    s__1 = *n;
    for (i__ = 2; i__ <= s__1; ++i__) {
/* L10: */
	j *= i__;
    }
    ret_val = (double) j;
    return ret_val;
L20:
    z__ = (double) (*n + 1);
    d__1 = z__ - .5;
    ret_val = exp(-z__) * pow(z__, d__1) * ((((-2.294720936e-4 / z__ - 
	    .00268132716) / z__ + .003472222222) / z__ + .08333333333) / z__ 
	    + 1.) * 2.506628274631;
    return ret_val;
} /* fctl_ */
#define ik1k0 (k1k0_1.ik1k0)

double bgama_(double *fnu, double *t1, double *t2, double 
	*eps, double *t3, double *t4, double *temp)
{
    /* Initialized data */

    static double pi = 3.1415926535898;
    static double clight = 29979245800.;
    static double hbar = 1.0545887e-27;
    static double boltz = 1.380662e-16;

    /* System generated locals */
    double ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    double omega, z__, k0, t0, bgambc, zp, xk1;

/*     SPECTRAL FUNCTION "EBC", FOR REFERENCE: */
/*     SEE "PHENOMENA INDUCED BY INTERMOLECULAR INTERACTIONS", */
/*     ED. G. BIRNBAUM; J. BORYSOW AND L. FROMMHOLD, P.67, (1985) */
/*       ============================================ */


/*       IF IK1K0=1 ONLY B-C; EBC OTHERWISE */

    omega = 2. * pi * clight * *fnu;
    t0 = hbar / (boltz * 2. * *temp);
/* Computing 2nd power */
    d__1 = omega * *t1;
    z__ = sqrt((d__1 * d__1 + 1.) * (*t2 * *t2 + t0 * t0)) / *t1;
    if (z__ - 2. <= 0.) {
	goto L10;
    } else {
	goto L12;
    }
L10:
/* Computing 2nd power */
    d__2 = z__ / 3.75;
    d__1 = d__2 * d__2;
/* Computing 2nd power */
    d__4 = z__ / 2.;
    d__3 = d__4 * d__4;
    xk1 = z__ * z__ * log(z__ / 2.) * ((((((3.2411e-4 * d__1 + .00301532) * 
	    d__1 + .02658733) * d__1 + .15084934) * d__1 + .51498869) * d__1 
	    + .87890594) * d__1 + .5) + ((((((-4.686e-5 * d__3 - .00110404) * 
	    d__3 - .01919402) * d__3 - .18156897) * d__3 - .67278579) * d__3 
	    + .15443144) * d__3 + 1.);
    goto L20;
L12:
    d__1 = 2. / z__;
    xk1 = sqrt(z__) * exp(-z__) * ((((((-6.8245e-4 * d__1 + .00325614) * d__1 
	    - .00780353) * d__1 + .01504268) * d__1 - .0365562) * d__1 + 
	    .23498619) * d__1 + 1.25331414);
L20:
/* Computing 2nd power */
    d__1 = *t1 * omega;
    bgambc = *t1 / pi * exp(*t2 / *t1 + t0 * omega) * xk1 / (d__1 * d__1 + 1.)
	    ;
    if (ik1k0 == 1) {
	goto L55;
    }
/* Computing 2nd power */
    d__1 = omega * *t4;
    zp = sqrt((d__1 * d__1 + 1.) * (*t3 * *t3 + t0 * t0)) / *t4;
    if (zp - 2. <= 0.) {
	goto L22;
    } else {
	goto L24;
    }
L22:
/* Computing 2nd power */
    d__2 = zp / 3.75;
    d__1 = d__2 * d__2;
/* Computing 2nd power */
    d__4 = zp / 2.;
    d__3 = d__4 * d__4;
    k0 = -log(zp / 2.) * ((((((.0045813 * d__1 + .0360768) * d__1 + .2659732) 
	    * d__1 + 1.2067492) * d__1 + 3.0899424) * d__1 + 3.5156229) * 
	    d__1 + 1.) + ((((((7.4e-6 * d__3 + 1.075e-4) * d__3 + .00262698) *
	     d__3 + .0348859) * d__3 + .23069756) * d__3 + .4227842) * d__3 - 
	    .57721566);
    goto L30;
L24:
    d__1 = 2. / zp;
    k0 = exp(-zp) * ((((((5.3208e-4 * d__1 - .0025154) * d__1 + .00587872) * 
	    d__1 - .01062446) * d__1 + .02189568) * d__1 - .07832358) * d__1 
	    + 1.25331414) / sqrt(zp);
L30:
    ret_val = (bgambc + *eps * (*t3 / pi) * exp(*t3 / *t4 + t0 * omega) * k0) 
	    / (*eps + 1.);
    goto L66;
L55:
    ret_val = bgambc;
L66:
    return ret_val;
} /* bgama_ */

#undef ik1k0


/* Subroutine */ int spline_0_(int n__, int *l, int *m, int *k,
	 double *eps, double *x, double *y, double *t, 
	double *ss, double *si, int *nr, double *s2)
{
    /* System generated locals */
    int s__1, s__2;
    double d__1, d__2;

    /* Local variables */
    double epsi, prod, h__;
    int i__, j, n;
    double w, omega;
    int n1;
    double s3;
    int ic;
    double sm, delsqs, ht1, ht2, ss2, yp1, eta, ypn;



/*     SPLINE INTERPOLATION AND QUADRATURE, THIRD ORDER AFTER GREVILLE. */
/*     INPUT ARGUMENTS L...Y, OUTPUT SS...NR. */
/*     L DATA POINTS X(1), Y(1) ... X(L),Y(L) */
/*     EPS=ERROR CRITERION, TYPICALLY EPS=1.D-5 FOR 5 DECI. PLACES ACCURA */
/*     M ARGUMENTS T(1)..T(M) FOR WHICH FUNCTION VALUES SS(1)..SS(M), FOR */
/*     K=0; OR FIRST OR SECOND DERIVATIVE FOR K=1 OR -1, RESPECTIVELY. */
/*     NOTE THAT M HAS TO BE AT LEAST EQUAL TO 1. */
/*     SI=INTEGRAL (OVER WHOLE INTERVAL) FOR K=2 ONLY. */
/*     FOR 'NATURAL' SPLINE FUNCTIONS, S2(1)=S2(L)=0. MUST BE INPUT*NOTE* */
/*     N0 INDICATES THE NUMBER OF OUT-OF-RANGE CALLS. X(1)<T(I)<X(L) */
/*     EXTRAPOLATE WITH CAUTION. (ASSUMPTION D2Y/DX2 = 0.) */
/*     S2(I) IS THE 2ND DERIVATIVE AT X=X(I) AND IS COMPUTED INTERNALLY. */
/*     DIMENSION X(L),Y(L),T(M),SS(M),S2(L) */

    /* Parameter adjustments */
    --x;
    --y;
    --t;
    --ss;
    --s2;

    /* Function Body */
    switch(n__) {
	case 1: goto L_ixpolat;
	}

    n = *l;
    n1 = n - 1;
    *nr = 0;
/* L4: */
    s__1 = n1;
    for (i__ = 2; i__ <= s__1; ++i__) {
/* L52: */
	s__2 = i__ - 1;
	s2[i__] = 3. * ((y[i__ + 1] - y[i__]) / (x[i__ + 1] - x[i__]) - (y[
		s__2 + 1] - y[s__2]) / (x[s__2 + 1] - x[s__2])) / (x[i__ + 1] 
		- x[i__ - 1]) / 1.5;
    }
    omega = 1.0717968;
    ic = 0;
/*     'NATURAL' SPLINE FUNCTIONS OF THIRD ORDER. */
    s2[1] = 0.;
    s2[n] = 0.;
L5:
    eta = 0.;
    ++ic;
    sm = abs(s2[1]);
    s__2 = n;
    for (i__ = 2; i__ <= s__2; ++i__) {
	if ((d__1 = s2[i__], abs(d__1)) > sm) {
	    sm = (d__2 = s2[i__], abs(d__2));
	}
/* L25: */
    }
    epsi = *eps * sm;
/* L6: */
    s__2 = n1;
    for (i__ = 2; i__ <= s__2; ++i__) {
/* L7: */
	s__1 = i__ - 1;
	w = (3. * ((y[i__ + 1] - y[i__]) / (x[i__ + 1] - x[i__]) - (y[s__1 + 
		1] - y[s__1]) / (x[s__1 + 1] - x[s__1])) / (x[i__ + 1] - x[
		i__ - 1]) - (x[i__] - x[i__ - 1]) * .5 / (x[i__ + 1] - x[i__ 
		- 1]) * s2[i__ - 1] - (.5 - (x[i__] - x[i__ - 1]) * .5 / (x[
		i__ + 1] - x[i__ - 1])) * s2[i__ + 1] - s2[i__]) * omega;
/* L8: */
	if (abs(w) - eta <= 0.) {
	    goto L10;
	} else {
	    goto L9;
	}
L9:
	eta = abs(w);
L10:
	s2[i__] += w;
    }
/* L13: */
    if (eta - epsi >= 0.) {
	goto L5;
    } else {
	goto L14;
    }
/*     ENTRY IXPOLAT */

L_ixpolat:
/*     THIS ENTRY USEFUL WHEN ITERATION PREVIOUSLY COMPLETED */

    n = *l;
    n1 = n - 1;
    *nr = 0;
    ic = -1;
L14:
    if (*k - 2 != 0) {
	goto L15;
    } else {
	goto L20;
    }
L15:
    s__2 = *m;
    for (j = 1; j <= s__2; ++j) {
/* L16: */
	i__ = 1;
/* L54: */
	if ((d__1 = t[j] - x[1]) < 0.) {
	    goto L58;
	} else if (d__1 == 0) {
	    goto L17;
	} else {
	    goto L55;
	}
L55:
	if ((d__1 = t[j] - x[n]) < 0.) {
	    goto L57;
	} else if (d__1 == 0) {
	    goto L59;
	} else {
	    goto L158;
	}
L56:
	if ((d__1 = t[j] - x[i__]) < 0.) {
	    goto L60;
	} else if (d__1 == 0) {
	    goto L17;
	} else {
	    goto L57;
	}
L57:
	++i__;
	goto L56;

L58:
	++(*nr);
	ht1 = t[j] - x[1];
	ht2 = t[j] - x[2];
	yp1 = (y[cs__1 + 1] - y[cs__1]) / (x[cs__1 + 1] - x[cs__1]) + (x[1] - 
		x[2]) * (s2[1] * 2. + s2[2]) / 6.;
	if (*k < 0) {
	    goto L72;
	} else if (*k == 0) {
	    goto L70;
	} else {
	    goto L71;
	}
L71:
	ss[j] = yp1 + ht1 * s2[1];
	goto L61;
L70:
	ss[j] = y[1] + yp1 * ht1 + s2[1] * ht1 * ht1 / 2.;
	goto L61;
L72:
	ss[j] = s2[i__];
	goto L61;
L158:
	ht2 = t[j] - x[n];
	ht1 = t[j] - x[n1];
	++(*nr);
	ypn = (y[n1 + 1] - y[n1]) / (x[n1 + 1] - x[n1]) + (x[n] - x[n1]) * (
		s2[n1] + s2[n] * 2.) / 6.;
	if (*k < 0) {
	    goto L82;
	} else if (*k == 0) {
	    goto L80;
	} else {
	    goto L81;
	}
L81:
	ss[j] = ypn + ht2 * s2[n];
	goto L61;
L80:
	ss[j] = y[n] + ypn * ht2 + s2[n] * ht2 * ht2 / 2.;
	goto L61;
L82:
	ss[j] = s2[n];
	goto L61;

L59:
	i__ = n;
L60:
	--i__;
L17:
	ht1 = t[j] - x[i__];
	ht2 = t[j] - x[i__ + 1];
	prod = ht1 * ht2;
	s3 = (s2[i__ + 1] - s2[i__]) / (x[i__ + 1] - x[i__]);
	ss2 = s2[i__] + ht1 * s3;
	delsqs = (s2[i__] + s2[i__ + 1] + ss2) / 6.;

	if (*k < 0) {
	    goto L43;
	} else if (*k == 0) {
	    goto L41;
	} else {
	    goto L42;
	}
L41:
	ss[j] = y[i__] + ht1 * ((y[i__ + 1] - y[i__]) / (x[i__ + 1] - x[i__]))
		 + prod * delsqs;
	goto L61;
L42:
	ss[j] = (y[i__ + 1] - y[i__]) / (x[i__ + 1] - x[i__]) + (ht1 + ht2) * 
		delsqs + prod * s3 / 6.;
	goto L61;
L43:
	ss[j] = ss2;
L61:
	;
    }
L20:
    *si = 0.;

    s__2 = n1;
    for (i__ = 1; i__ <= s__2; ++i__) {
	h__ = x[i__ + 1] - x[i__];
/* L62: */
/* Computing 3rd power */
	d__1 = h__;
	*si = *si + h__ * .5 * (y[i__] + y[i__ + 1]) - d__1 * (d__1 * d__1) * 
		(s2[i__] + s2[i__ + 1]) / 24.;
    }

    if (*k == 2) {
	*nr = ic;
    }

    return 0;
} /* spline_ */

/* Subroutine */ int spline_(int *l, int *m, int *k, 
	double *eps, double *x, double *y, double *t, 
	double *ss, double *si, int *nr, double *s2)
{
    return spline_0_(0, l, m, k, eps, x, y, t, ss, si, nr, s2);
    }

/* Subroutine */ int ixpolat_(int *l, int *m, int *k, 
	double *eps, double *x, double *y, double *t, 
	double *ss, double *si, int *nr, double *s2)
{
    return spline_0_(1, l, m, k, eps, x, y, t, ss, si, nr, s2);
    }


// ---------------------- end of Borysow N2N2 F77 code -------------------------


// ---------------------- begin of monortm CKD F77 code -------------------------


/* Common Block Declarations */

struct fh2oa_1_ {
    double fh2o[2003];
};

#define fh2oa_1 (*(struct fh2oa_1_ *) &fh2oa_)

struct fh2ob_1_ {
    double v1, v2, dv;
    int nptfh2o;
};
struct fh2ob_2_ {
    double v1, v2, dv;
    int npts;
};

#define fh2ob_1 (*(struct fh2ob_1_ *) &fh2ob_)
#define fh2ob_2 (*(struct fh2ob_2_ *) &fh2ob_)

struct sh2oa_1_ {
    double swv296[2003];
};

#define sh2oa_1 (*(struct sh2oa_1_ *) &sh2oa_)

struct sh2ob_1_ {
    double v1, v2, dv;
    int nptslfwv;
};
struct sh2ob_2_ {
    double v1, v2, dv;
    int npts;
};

#define sh2ob_1 (*(struct sh2ob_1_ *) &sh2ob_)
#define sh2ob_2 (*(struct sh2ob_2_ *) &sh2ob_)

struct s260a_1_ {
    double swv260[2003];
};

#define s260a_1 (*(struct s260a_1_ *) &s260a_)

struct s260b_1_ {
    double v1___, v2___, dv___;
    int nptslfwv___;
};
struct s260b_2_ {
    double v1, v2, dv;
    int npts;
};

#define s260b_1 (*(struct s260b_1_ *) &s260b_)
#define s260b_2 (*(struct s260b_2_ *) &s260b_)

struct consts_1_ {
    double pi, planck, boltz, clight, avogad, alosmt, gascon, radcn1, 
	    radcn2;
};

#define consts_1 (*(struct consts_1_ *) &consts_)

/* Initialized data */

struct {
    double e_1[2003];
    } fh2oa_ = { .012859, .011715, .011038, .011715, .012859, .015326, 
	    .016999, .018321, .019402, .01957, .019432, .017572, .01676, 
	    .01548, .013984, .012266, .010467, .0094526, .0080485, .0069484, 
	    .0061416, .0050941, .0044836, .0038133, .0034608, .0031487, 
	    .0024555, .0020977, .0017266, .001492, .0012709, 9.8081e-4, 
	    8.5063e-4, 6.8822e-4, 5.3809e-4, 4.4679e-4, 3.3774e-4, 2.7979e-4, 
	    2.1047e-4, 1.6511e-4, 1.2993e-4, 9.3033e-5, 7.436e-5, 5.6428e-5, 
	    4.5442e-5, 3.4575e-5, 2.7903e-5, 2.1374e-5, 1.6075e-5, 1.3022e-5, 
	    1.0962e-5, 8.5959e-6, 6.9125e-6, 5.3808e-6, 4.3586e-6, 3.6394e-6, 
	    2.9552e-6, 2.3547e-6, 1.8463e-6, 1.6036e-6, 1.3483e-6, 1.1968e-6, 
	    1.0333e-6, 8.4484e-7, 6.7195e-7, 5.0947e-7, 4.2343e-7, 3.4453e-7, 
	    2.783e-7, 2.3063e-7, 1.9951e-7, 1.7087e-7, 1.4393e-7, 1.2575e-7, 
	    1.075e-7, 8.2325e-8, 5.7524e-8, 4.4482e-8, 3.8106e-8, 3.4315e-8, 
	    2.9422e-8, 2.5069e-8, 2.2402e-8, 1.9349e-8, 1.6152e-8, 1.2208e-8, 
	    8.966e-9, 7.1322e-9, 6.1028e-9, 5.2938e-9, 4.535e-9, 3.4977e-9, 
	    2.9511e-9, 2.4734e-9, 2.0508e-9, 1.8507e-9, 1.6373e-9, 1.5171e-9, 
	    1.3071e-9, 1.2462e-9, 1.2148e-9, 1.259e-9, 1.3153e-9, 1.3301e-9, 
	    1.4483e-9, 1.6944e-9, 2.0559e-9, 2.2954e-9, 2.6221e-9, 3.2606e-9, 
	    4.2392e-9, 5.2171e-9, 6.2553e-9, 8.2548e-9, 9.5842e-9, 1.128e-8, 
	    1.3628e-8, 1.7635e-8, 2.1576e-8, 2.4835e-8, 3.0014e-8, 3.8485e-8, 
	    4.744e-8, 5.5202e-8, 7.0897e-8, 9.6578e-8, 1.3976e-7, 1.8391e-7, 
	    2.3207e-7, 2.996e-7, 4.0408e-7, 5.926e-7, 7.8487e-7, 1.0947e-6, 
	    1.4676e-6, 1.9325e-6, 2.6587e-6, 3.4534e-6, 4.4376e-6, 5.8061e-6, 
	    7.0141e-6, 8.4937e-6, 1.0186e-5, 1.2034e-5, 1.3837e-5, 1.6595e-5, 
	    1.9259e-5, 2.162e-5, 2.3681e-5, 2.7064e-5, 3.251e-5, 3.546e-5, 
	    3.9109e-5, 4.2891e-5, 4.7757e-5, 5.0981e-5, 5.0527e-5, 4.8618e-5, 
	    4.4001e-5, 3.7982e-5, 3.2667e-5, 2.7794e-5, 2.491e-5, 2.4375e-5, 
	    2.7316e-5, 3.2579e-5, 3.5499e-5, 3.801e-5, 4.1353e-5, 4.3323e-5, 
	    4.3004e-5, 3.979e-5, 3.7718e-5, 3.636e-5, 3.2386e-5, 2.7409e-5, 
	    2.3626e-5, 2.0631e-5, 1.8371e-5, 1.5445e-5, 1.2989e-5, 1.1098e-5, 
	    9.6552e-6, 8.0649e-6, 7.2365e-6, 5.9137e-6, 5.2759e-6, 4.886e-6, 
	    4.1321e-6, 3.5918e-6, 2.764e-6, 2.4892e-6, 2.1018e-6, 1.7848e-6, 
	    1.5855e-6, 1.3569e-6, 1.1986e-6, 9.4693e-7, 7.4097e-7, 6.3443e-7, 
	    4.8131e-7, 4.0942e-7, 3.3316e-7, 2.8488e-7, 2.3461e-7, 1.7397e-7, 
	    1.4684e-7, 1.0953e-7, 8.5396e-8, 6.9261e-8, 5.4001e-8, 4.543e-8, 
	    3.2791e-8, 2.5995e-8, 2.0225e-8, 1.571e-8, 1.3027e-8, 1.0229e-8, 
	    8.5277e-9, 6.5249e-9, 5.0117e-9, 3.9906e-9, 3.2332e-9, 2.7847e-9, 
	    2.457e-9, 2.3359e-9, 2.0599e-9, 1.8436e-9, 1.6559e-9, 1.491e-9, 
	    1.2794e-9, 9.8229e-10, 8.0054e-10, 6.0769e-10, 4.5646e-10, 
	    3.3111e-10, 2.4428e-10, 1.8007e-10, 1.3291e-10, 9.7974e-11, 
	    7.8271e-11, 6.3833e-11, 5.4425e-11, 4.6471e-11, 4.0209e-11, 
	    3.5227e-11, 3.1212e-11, 2.884e-11, 2.7762e-11, 2.7935e-11, 
	    3.2012e-11, 3.9525e-11, 5.0303e-11, 6.8027e-11, 9.3954e-11, 
	    1.2986e-10, 1.8478e-10, 2.5331e-10, 3.4827e-10, 4.6968e-10, 
	    6.238e-10, 7.9106e-10, 1.0026e-9, 1.2102e-9, 1.4146e-9, 1.6154e-9,
	     1.751e-9, 1.8575e-9, 1.8742e-9, 1.87e-9, 1.8582e-9, 1.9657e-9, 
	    2.1204e-9, 2.0381e-9, 2.0122e-9, 2.0436e-9, 2.1213e-9, 2.0742e-9, 
	    1.987e-9, 2.0465e-9, 2.1556e-9, 2.2222e-9, 2.1977e-9, 2.1047e-9, 
	    1.9334e-9, 1.7357e-9, 1.5754e-9, 1.4398e-9, 1.4018e-9, 1.5459e-9, 
	    1.7576e-9, 2.1645e-9, 2.948e-9, 4.4439e-9, 5.8341e-9, 8.0757e-9, 
	    1.1658e-8, 1.6793e-8, 2.2694e-8, 2.9468e-8, 3.9278e-8, 5.2145e-8, 
	    6.4378e-8, 7.7947e-8, 8.5321e-8, 9.7848e-8, 1.0999e-7, 1.1489e-7, 
	    1.2082e-7, 1.2822e-7, 1.4053e-7, 1.5238e-7, 1.5454e-7, 1.5018e-7, 
	    1.4048e-7, 1.2359e-7, 1.0858e-7, 9.3486e-8, 8.1638e-8, 7.769e-8, 
	    8.4625e-8, 1.0114e-7, 1.143e-7, 1.2263e-7, 1.3084e-7, 1.338e-7, 
	    1.3573e-7, 1.3441e-7, 1.2962e-7, 1.2638e-7, 1.1934e-7, 1.1371e-7, 
	    1.0871e-7, 9.8843e-8, 9.1877e-8, 9.105e-8, 9.3213e-8, 9.2929e-8, 
	    1.0155e-7, 1.1263e-7, 1.237e-7, 1.3636e-7, 1.54e-7, 1.7656e-7, 
	    2.1329e-7, 2.3045e-7, 2.5811e-7, 2.9261e-7, 3.4259e-7, 4.077e-7, 
	    4.8771e-7, 5.8081e-7, 7.2895e-7, 8.7482e-7, 1.0795e-6, 1.3384e-6, 
	    1.7208e-6, 2.0677e-6, 2.5294e-6, 3.1123e-6, 3.79e-6, 4.7752e-6, 
	    5.6891e-6, 6.6261e-6, 7.6246e-6, 8.773e-6, 9.6672e-6, 1.098e-5, 
	    1.1287e-5, 1.167e-5, 1.1635e-5, 1.1768e-5, 1.2039e-5, 1.2253e-5, 
	    1.3294e-5, 1.4005e-5, 1.3854e-5, 1.342e-5, 1.3003e-5, 1.2645e-5, 
	    1.1715e-5, 1.1258e-5, 1.1516e-5, 1.2494e-5, 1.3655e-5, 1.4931e-5, 
	    1.4649e-5, 1.3857e-5, 1.312e-5, 1.1791e-5, 1.0637e-5, 8.276e-6, 
	    6.5821e-6, 5.1959e-6, 4.0158e-6, 3.0131e-6, 2.0462e-6, 1.4853e-6, 
	    1.0365e-6, 7.3938e-7, 4.9752e-7, 3.4148e-7, 2.4992e-7, 1.8363e-7, 
	    1.4591e-7, 1.138e-7, 9.0588e-8, 7.3697e-8, 6.0252e-8, 5.1868e-8, 
	    4.266e-8, 3.6163e-8, 3.2512e-8, 2.9258e-8, 2.4238e-8, 2.1209e-8, 
	    1.6362e-8, 1.3871e-8, 1.2355e-8, 9.694e-9, 7.7735e-9, 6.2278e-9, 
	    5.2282e-9, 4.3799e-9, 3.5545e-9, 2.7527e-9, 2.095e-9, 1.6344e-9, 
	    1.2689e-9, 1.0403e-9, 8.488e-10, 6.3461e-10, 4.7657e-10, 
	    3.522e-10, 2.7879e-10, 2.3021e-10, 1.6167e-10, 1.1732e-10, 
	    8.9206e-11, 7.0596e-11, 5.831e-11, 4.4084e-11, 3.1534e-11, 
	    2.5068e-11, 2.2088e-11, 2.2579e-11, 2.2637e-11, 2.5705e-11, 
	    3.2415e-11, 4.6116e-11, 6.5346e-11, 9.4842e-11, 1.2809e-10, 
	    1.8211e-10, 2.4052e-10, 3.027e-10, 3.5531e-10, 4.2402e-10, 
	    4.673e-10, 4.7942e-10, 4.6813e-10, 4.5997e-10, 4.5788e-10, 
	    4.0311e-10, 3.7367e-10, 3.3149e-10, 2.9281e-10, 2.5231e-10, 
	    2.1152e-10, 1.9799e-10, 1.8636e-10, 1.9085e-10, 2.0786e-10, 
	    2.2464e-10, 2.3785e-10, 2.5684e-10, 2.7499e-10, 2.6962e-10, 
	    2.6378e-10, 2.6297e-10, 2.6903e-10, 2.7035e-10, 2.5394e-10, 
	    2.5655e-10, 2.7184e-10, 2.9013e-10, 3.0585e-10, 3.0791e-10, 
	    3.1667e-10, 3.4343e-10, 3.7365e-10, 4.0269e-10, 4.726e-10, 
	    5.6584e-10, 6.9791e-10, 8.6569e-10, 1.0393e-9, 1.2067e-9, 
	    1.5047e-9, 1.8583e-9, 2.2357e-9, 2.6498e-9, 3.2483e-9, 3.9927e-9, 
	    4.6618e-9, 5.5555e-9, 6.6609e-9, 8.2139e-9, 1.0285e-8, 1.3919e-8, 
	    1.8786e-8, 2.515e-8, 3.313e-8, 4.5442e-8, 6.337e-8, 9.0628e-8, 
	    1.2118e-7, 1.5927e-7, 2.1358e-7, 2.7825e-7, 3.7671e-7, 4.4894e-7, 
	    5.4442e-7, 6.224e-7, 7.3004e-7, 8.3384e-7, 8.7933e-7, 8.808e-7, 
	    8.6939e-7, 8.6541e-7, 8.2055e-7, 7.7278e-7, 7.5989e-7, 8.6909e-7, 
	    9.7945e-7, 1.0394e-6, 1.0646e-6, 1.1509e-6, 1.2017e-6, 1.1915e-6, 
	    1.1259e-6, 1.1549e-6, 1.1938e-6, 1.2356e-6, 1.2404e-6, 1.1716e-6, 
	    1.1149e-6, 1.0073e-6, 8.9845e-7, 7.6639e-7, 6.1517e-7, 5.0887e-7, 
	    4.1269e-7, 3.2474e-7, 2.5698e-7, 1.8893e-7, 1.4009e-7, 1.034e-7, 
	    7.7724e-8, 5.7302e-8, 4.2178e-8, 2.9603e-8, 2.1945e-8, 1.6301e-8, 
	    1.2806e-8, 1.0048e-8, 7.897e-9, 6.1133e-9, 4.9054e-9, 4.1985e-9, 
	    3.6944e-9, 3.2586e-9, 2.7362e-9, 2.3647e-9, 2.1249e-9, 1.8172e-9, 
	    1.6224e-9, 1.5158e-9, 1.2361e-9, 1.0682e-9, 9.2312e-10, 7.922e-10,
	     6.8174e-10, 5.6147e-10, 4.8268e-10, 4.1534e-10, 3.3106e-10, 
	    2.8275e-10, 2.4584e-10, 2.0742e-10, 1.784e-10, 1.4664e-10, 
	    1.239e-10, 1.0497e-10, 8.5038e-11, 6.7008e-11, 5.6355e-11, 
	    4.3323e-11, 3.6914e-11, 3.2262e-11, 3.0749e-11, 3.0318e-11, 
	    2.9447e-11, 2.9918e-11, 3.0668e-11, 3.1315e-11, 3.0329e-11, 
	    2.8259e-11, 2.6065e-11, 2.3578e-11, 2.0469e-11, 1.6908e-11, 
	    1.4912e-11, 1.1867e-11, 9.973e-12, 8.1014e-12, 6.7528e-12, 
	    6.3133e-12, 5.8599e-12, 6.0145e-12, 6.5105e-12, 7.0537e-12, 
	    7.4973e-12, 7.8519e-12, 8.5039e-12, 9.1995e-12, 1.0694e-11, 
	    1.1659e-11, 1.2685e-11, 1.3087e-11, 1.3222e-11, 1.2634e-11, 
	    1.1077e-11, 9.6259e-12, 8.3202e-12, 7.4857e-12, 6.8069e-12, 
	    6.7496e-12, 7.3116e-12, 8.0171e-12, 8.6394e-12, 9.2659e-12, 
	    1.0048e-11, 1.0941e-11, 1.2226e-11, 1.3058e-11, 1.5193e-11, 
	    1.8923e-11, 2.3334e-11, 2.8787e-11, 3.6693e-11, 4.8295e-11, 
	    6.426e-11, 8.8269e-11, 1.1865e-10, 1.5961e-10, 2.0605e-10, 
	    2.7349e-10, 3.7193e-10, 4.8216e-10, 6.1966e-10, 7.715e-10, 
	    1.0195e-9, 1.2859e-9, 1.6535e-9, 2.0316e-9, 2.3913e-9, 3.0114e-9, 
	    3.7495e-9, 4.6504e-9, 5.9145e-9, 7.684e-9, 1.0304e-8, 1.301e-8, 
	    1.6441e-8, 2.1475e-8, 2.5892e-8, 2.9788e-8, 3.382e-8, 4.0007e-8, 
	    4.4888e-8, 4.5765e-8, 4.6131e-8, 4.6239e-8, 4.4849e-8, 4.0729e-8, 
	    3.6856e-8, 3.6164e-8, 3.7606e-8, 4.1457e-8, 4.375e-8, 5.115e-8, 
	    5.6054e-8, 6.1586e-8, 6.4521e-8, 6.6494e-8, 6.9024e-8, 6.8893e-8, 
	    7.0901e-8, 6.976e-8, 7.1485e-8, 7.074e-8, 7.3764e-8, 7.6618e-8, 
	    8.4182e-8, 9.3838e-8, 1.0761e-7, 1.2851e-7, 1.4748e-7, 1.8407e-7, 
	    2.2109e-7, 2.6392e-7, 2.9887e-7, 3.4493e-7, 4.0336e-7, 4.3551e-7, 
	    4.9231e-7, 5.0728e-7, 5.3781e-7, 5.3285e-7, 5.4496e-7, 5.5707e-7, 
	    5.6944e-7, 6.1123e-7, 6.4317e-7, 6.4581e-7, 6.1999e-7, 6.0191e-7, 
	    5.7762e-7, 5.7241e-7, 5.7013e-7, 6.016e-7, 6.6905e-7, 7.4095e-7, 
	    8.2121e-7, 8.0947e-7, 7.6145e-7, 7.2193e-7, 6.3722e-7, 5.4316e-7, 
	    4.2186e-7, 3.2528e-7, 2.5207e-7, 1.8213e-7, 1.2658e-7, 8.6746e-8, 
	    6.0216e-8, 4.1122e-8, 2.8899e-8, 2.174e-8, 1.799e-8, 1.5593e-8, 
	    1.397e-8, 1.2238e-8, 1.0539e-8, 9.2386e-9, 7.8481e-9, 6.8704e-9, 
	    5.7615e-9, 5.0434e-9, 4.6886e-9, 4.377e-9, 3.9768e-9, 3.5202e-9, 
	    3.1854e-9, 2.9009e-9, 2.5763e-9, 2.2135e-9, 1.9455e-9, 1.6248e-9, 
	    1.3368e-9, 1.0842e-9, 8.4254e-10, 6.7414e-10, 5.4667e-10, 
	    4.5005e-10, 3.4932e-10, 2.6745e-10, 2.2053e-10, 1.8162e-10, 
	    1.4935e-10, 1.1618e-10, 9.1888e-11, 8.0672e-11, 6.8746e-11, 
	    6.2668e-11, 5.5715e-11, 4.5074e-11, 3.7669e-11, 3.2082e-11, 
	    2.8085e-11, 2.4838e-11, 1.9791e-11, 1.6964e-11, 1.3887e-11, 
	    1.1179e-11, 9.7499e-12, 7.8255e-12, 6.3698e-12, 5.3265e-12, 
	    4.6588e-12, 4.4498e-12, 3.9984e-12, 3.7513e-12, 3.7176e-12, 
	    3.9148e-12, 4.2702e-12, 5.009e-12, 6.5801e-12, 8.7787e-12, 
	    1.2718e-11, 1.8375e-11, 2.5304e-11, 3.5403e-11, 4.8842e-11, 
	    6.484e-11, 8.0911e-11, 1.0136e-10, 1.2311e-10, 1.4203e-10, 
	    1.5869e-10, 1.8093e-10, 2.137e-10, 2.5228e-10, 2.8816e-10, 
	    3.4556e-10, 3.986e-10, 4.435e-10, 4.776e-10, 5.2357e-10, 
	    6.0827e-10, 6.3635e-10, 6.5886e-10, 6.8753e-10, 7.2349e-10, 
	    7.2789e-10, 6.8232e-10, 6.6081e-10, 6.4232e-10, 6.3485e-10, 
	    6.4311e-10, 7.2235e-10, 7.7263e-10, 8.1668e-10, 9.0324e-10, 
	    9.7643e-10, 1.0535e-9, 1.0195e-9, 1.0194e-9, 1.0156e-9, 
	    9.6792e-10, 9.2725e-10, 8.7347e-10, 8.4484e-10, 8.2647e-10, 
	    8.4363e-10, 9.1261e-10, 1.0051e-9, 1.1511e-9, 1.4037e-9, 
	    1.8066e-9, 2.4483e-9, 3.2739e-9, 4.3194e-9, 5.6902e-9, 7.7924e-9, 
	    9.7376e-9, 1.2055e-8, 1.4303e-8, 1.6956e-8, 1.9542e-8, 2.2233e-8, 
	    2.5186e-8, 2.7777e-8, 2.8943e-8, 2.8873e-8, 2.9417e-8, 2.7954e-8, 
	    2.7524e-8, 2.704e-8, 3.1254e-8, 3.6843e-8, 3.7797e-8, 3.8713e-8, 
	    4.0135e-8, 4.2824e-8, 4.3004e-8, 4.0279e-8, 4.2781e-8, 4.522e-8, 
	    4.8948e-8, 5.0172e-8, 4.8499e-8, 4.7182e-8, 4.2204e-8, 3.7701e-8, 
	    3.0972e-8, 2.4654e-8, 1.9543e-8, 1.4609e-8, 1.1171e-8, 8.3367e-9, 
	    6.3791e-9, 5.079e-9, 4.0655e-9, 3.3658e-9, 2.7882e-9, 2.4749e-9, 
	    2.2287e-9, 2.0217e-9, 1.8191e-9, 1.5897e-9, 1.4191e-9, 1.2448e-9, 
	    1.0884e-9, 9.3585e-10, 7.9429e-10, 7.3214e-10, 6.5008e-10, 
	    5.7549e-10, 5.43e-10, 4.7251e-10, 4.3451e-10, 3.8446e-10, 
	    3.5589e-10, 3.4432e-10, 2.8209e-10, 2.462e-10, 2.1278e-10, 
	    1.8406e-10, 1.6314e-10, 1.3261e-10, 1.1696e-10, 9.6865e-11, 
	    7.6814e-11, 6.6411e-11, 5.0903e-11, 4.0827e-11, 3.0476e-11, 
	    2.323e-11, 1.7707e-11, 1.3548e-11, 1.0719e-11, 9.3026e-12, 
	    8.7967e-12, 8.3136e-12, 7.3918e-12, 6.5293e-12, 5.9243e-12, 
	    5.3595e-12, 3.5266e-12, 2.2571e-12, 1.615e-12, 1.1413e-12, 
	    8.4998e-13, 7.0803e-13, 5.1747e-13, 4.0694e-13, 3.6528e-13, 
	    3.367e-13, 3.1341e-13, 2.939e-13, 2.868e-13, 3.1283e-13, 
	    3.7294e-13, 5.0194e-13, 6.7919e-13, 1.0455e-12, 1.523e-12, 
	    2.3932e-12, 3.4231e-12, 5.0515e-12, 7.3193e-12, 9.9406e-12, 
	    1.2193e-11, 1.4742e-11, 1.9269e-11, 2.1816e-11, 2.275e-11, 
	    2.2902e-11, 2.3888e-11, 2.4902e-11, 2.216e-11, 2.0381e-11, 
	    1.9903e-11, 2.0086e-11, 1.9304e-11, 2.0023e-11, 2.2244e-11, 
	    2.545e-11, 3.1228e-11, 3.456e-11, 3.6923e-11, 3.7486e-11, 
	    3.8124e-11, 3.8317e-11, 3.4737e-11, 3.3037e-11, 3.1724e-11, 
	    2.984e-11, 2.8301e-11, 2.5857e-11, 2.3708e-11, 1.9452e-11, 
	    1.6232e-11, 1.5174e-11, 1.4206e-11, 1.4408e-11, 1.5483e-11, 
	    1.8642e-11, 2.3664e-11, 3.0181e-11, 4.016e-11, 5.2287e-11, 
	    7.2754e-11, 1.0511e-10, 1.4531e-10, 2.0998e-10, 2.6883e-10, 
	    3.3082e-10, 4.2638e-10, 5.3132e-10, 6.3617e-10, 7.1413e-10, 
	    8.5953e-10, 9.9715e-10, 1.0796e-9, 1.0978e-9, 1.1052e-9, 
	    1.1095e-9, 1.0641e-9, 9.7881e-10, 9.659e-10, 1.0332e-9, 1.1974e-9,
	     1.3612e-9, 1.5829e-9, 1.8655e-9, 2.1465e-9, 2.4779e-9, 2.737e-9, 
	    2.9915e-9, 3.3037e-9, 3.6347e-9, 3.9587e-9, 4.4701e-9, 5.0122e-9, 
	    5.8044e-9, 6.1916e-9, 6.9613e-9, 7.7863e-9, 8.282e-9, 9.4359e-9, 
	    9.7387e-9, 1.0656e-8, 1.0746e-8, 1.121e-8, 1.1905e-8, 1.2194e-8, 
	    1.3145e-8, 1.3738e-8, 1.3634e-8, 1.3011e-8, 1.2511e-8, 1.1805e-8, 
	    1.2159e-8, 1.239e-8, 1.3625e-8, 1.5678e-8, 1.7886e-8, 1.9933e-8, 
	    1.9865e-8, 1.9e-8, 1.7812e-8, 1.5521e-8, 1.2593e-8, 9.5635e-9, 
	    7.2987e-9, 5.2489e-9, 3.5673e-9, 2.4206e-9, 1.6977e-9, 1.2456e-9, 
	    9.3744e-10, 7.8379e-10, 6.996e-10, 6.6451e-10, 6.8521e-10, 
	    7.4234e-10, 8.6658e-10, 9.4972e-10, 1.0791e-9, 1.2359e-9, 
	    1.3363e-9, 1.5025e-9, 1.5368e-9, 1.6152e-9, 1.6184e-9, 1.6557e-9, 
	    1.7035e-9, 1.6916e-9, 1.7237e-9, 1.7175e-9, 1.6475e-9, 1.5335e-9, 
	    1.4272e-9, 1.3282e-9, 1.3459e-9, 1.4028e-9, 1.5192e-9, 1.7068e-9, 
	    1.9085e-9, 2.1318e-9, 2.102e-9, 1.9942e-9, 1.8654e-9, 1.6391e-9, 
	    1.3552e-9, 1.0186e-9, 7.854e-10, 5.7022e-10, 3.9247e-10, 
	    2.5441e-10, 1.6699e-10, 1.1132e-10, 6.8989e-11, 4.5255e-11, 
	    3.1106e-11, 2.3161e-11, 1.7618e-11, 1.438e-11, 1.1601e-11, 
	    9.7148e-12, 8.4519e-12, 6.5392e-12, 5.4113e-12, 4.7624e-12, 
	    4.0617e-12, 3.6173e-12, 2.8608e-12, 2.2724e-12, 1.7436e-12, 
	    1.3424e-12, 1.0358e-12, 7.3064e-13, 5.45e-13, 4.0551e-13, 
	    2.8642e-13, 2.1831e-13, 1.686e-13, 1.2086e-13, 1.015e-13, 
	    9.355e-14, 8.4105e-14, 7.3051e-14, 6.9796e-14, 7.9949e-14, 
	    1.0742e-13, 1.5639e-13, 2.1308e-13, 3.1226e-13, 4.6853e-13, 
	    6.6917e-13, 1.0088e-12, 1.4824e-12, 2.2763e-12, 3.3917e-12, 
	    4.4585e-12, 6.3187e-12, 8.4189e-12, 1.1302e-11, 1.3431e-11, 
	    1.5679e-11, 1.9044e-11, 2.2463e-11, 2.3605e-11, 2.3619e-11, 
	    2.3505e-11, 2.3805e-11, 2.2549e-11, 1.9304e-11, 1.8382e-11, 
	    1.7795e-11, 1.8439e-11, 1.9146e-11, 2.1966e-11, 2.6109e-11, 
	    3.1883e-11, 3.7872e-11, 4.3966e-11, 4.8789e-11, 5.3264e-11, 
	    5.9705e-11, 6.3744e-11, 7.0163e-11, 7.9114e-11, 8.8287e-11, 
	    9.9726e-11, 1.1498e-10, 1.37e-10, 1.6145e-10, 1.9913e-10, 
	    2.2778e-10, 2.6216e-10, 2.977e-10, 3.3405e-10, 3.7821e-10, 
	    3.9552e-10, 4.1322e-10, 4.0293e-10, 4.0259e-10, 3.8853e-10, 
	    3.7842e-10, 3.8551e-10, 4.4618e-10, 5.0527e-10, 5.0695e-10, 
	    5.1216e-10, 5.193e-10, 5.5794e-10, 5.332e-10, 5.2008e-10, 
	    5.6888e-10, 6.1883e-10, 6.9006e-10, 6.9505e-10, 6.6768e-10, 
	    6.329e-10, 5.6753e-10, 5.0327e-10, 3.983e-10, 3.1147e-10, 
	    2.4416e-10, 1.886e-10, 1.3908e-10, 9.9156e-11, 7.3779e-11, 
	    5.6048e-11, 4.2457e-11, 3.4505e-11, 2.9881e-11, 2.7865e-11, 
	    2.8471e-11, 3.1065e-11, 3.4204e-11, 3.914e-11, 4.3606e-11, 
	    4.9075e-11, 5.3069e-11, 5.5236e-11, 5.5309e-11, 5.3832e-11, 
	    5.3183e-11, 5.1783e-11, 5.2042e-11, 5.4422e-11, 5.5656e-11, 
	    5.4409e-11, 5.2659e-11, 5.1696e-11, 5.1726e-11, 4.9003e-11, 
	    4.905e-11, 5.17e-11, 5.6818e-11, 6.3129e-11, 6.6542e-11, 
	    6.4367e-11, 5.9908e-11, 5.447e-11, 4.7903e-11, 3.9669e-11, 
	    2.9651e-11, 2.2286e-11, 1.6742e-11, 1.1827e-11, 7.7739e-12, 
	    4.8805e-12, 3.1747e-12, 2.0057e-12, 1.255e-12, 8.7434e-13, 
	    6.2755e-13, 4.9752e-13, 4.0047e-13, 3.5602e-13, 3.093e-13, 
	    2.4903e-13, 1.9316e-13, 1.4995e-13, 1.2059e-13, 8.7242e-14, 
	    6.4511e-14, 5.33e-14, 4.3741e-14, 3.4916e-14, 2.656e-14, 
	    1.6923e-14, 1.1816e-14, 6.7071e-15, 3.6474e-15, 2.0686e-15, 
	    1.1925e-15, 6.8948e-16, 3.9661e-16, 2.2576e-16, 1.2669e-16, 
	    6.9908e-17, 3.7896e-17, 2.028e-17, 1.1016e-17, 6.7816e-18, 
	    6.0958e-18, 8.9913e-18, 1.7201e-17, 3.4964e-17, 7.0722e-17, 
	    1.402e-16, 2.7167e-16, 5.1478e-16, 9.55e-16, 1.7376e-15, 
	    3.1074e-15, 5.4789e-15, 9.564e-15, 1.6635e-14, 2.9145e-14, 
	    5.2179e-14, 8.8554e-14, 1.4764e-13, 2.3331e-13, 3.5996e-13, 
	    5.2132e-13, 6.3519e-13, 7.3174e-13, 8.3752e-13, 9.8916e-13, 
	    1.1515e-12, 1.4034e-12, 1.6594e-12, 2.1021e-12, 2.7416e-12, 
	    3.4135e-12, 4.5517e-12, 5.5832e-12, 7.2303e-12, 9.9484e-12, 
	    1.2724e-11, 1.6478e-11, 2.0588e-11, 2.5543e-11, 3.3625e-11, 
	    4.1788e-11, 5.0081e-11, 6.0144e-11, 6.9599e-11, 8.4408e-11, 
	    9.7143e-11, 1.0805e-10, 1.1713e-10, 1.2711e-10, 1.3727e-10, 
	    1.4539e-10, 1.6049e-10, 1.768e-10, 2.0557e-10, 2.4967e-10, 
	    3.0096e-10, 3.5816e-10, 4.0851e-10, 4.6111e-10, 5.2197e-10, 
	    5.5043e-10, 6.0324e-10, 6.4983e-10, 6.7498e-10, 7.0545e-10, 
	    7.068e-10, 7.5218e-10, 7.5723e-10, 7.784e-10, 8.0081e-10, 
	    8.0223e-10, 7.7271e-10, 7.1676e-10, 6.7819e-10, 6.4753e-10, 
	    6.5844e-10, 7.0163e-10, 7.7503e-10, 8.8152e-10, 9.9022e-10, 
	    1.0229e-9, 9.9296e-10, 8.9911e-10, 7.7813e-10, 6.3785e-10, 
	    4.7491e-10, 3.528e-10, 2.4349e-10, 1.6502e-10, 1.1622e-10, 
	    8.6715e-11, 6.736e-11, 5.391e-11, 4.5554e-11, 4.13e-11, 
	    3.9728e-11, 3.9e-11, 3.9803e-11, 4.1514e-11, 4.3374e-11, 
	    4.6831e-11, 4.8921e-11, 5.1995e-11, 5.7242e-11, 6.2759e-11, 
	    7.0801e-11, 7.4555e-11, 7.9754e-11, 8.7616e-11, 9.1171e-11, 
	    1.0349e-10, 1.1047e-10, 1.2024e-10, 1.299e-10, 1.3725e-10, 
	    1.5005e-10, 1.5268e-10, 1.5535e-10, 1.5623e-10, 1.5009e-10, 
	    1.4034e-10, 1.3002e-10, 1.2225e-10, 1.1989e-10, 1.2411e-10, 
	    1.3612e-10, 1.5225e-10, 1.7202e-10, 1.9471e-10, 1.9931e-10, 
	    1.9079e-10, 1.7478e-10, 1.5259e-10, 1.2625e-10, 9.3332e-11, 
	    6.8796e-11, 4.6466e-11, 2.9723e-11, 1.8508e-11, 1.2106e-11, 
	    8.0142e-12, 5.4066e-12, 3.9329e-12, 3.1665e-12, 2.742e-12, 
	    2.3996e-12, 2.3804e-12, 2.3242e-12, 2.4476e-12, 2.5331e-12, 
	    2.3595e-12, 2.2575e-12, 2.1298e-12, 2.0088e-12, 1.8263e-12, 
	    1.6114e-12, 1.4422e-12, 1.2946e-12, 1.0837e-12, 9.1282e-13, 
	    7.2359e-13, 5.3307e-13, 3.8837e-13, 2.6678e-13, 1.6769e-13, 
	    1.0826e-13, 7.2364e-14, 4.5201e-14, 3.0808e-14, 2.2377e-14, 
	    1.704e-14, 9.2181e-15, 5.2934e-15, 3.5774e-15, 3.1431e-15, 
	    3.7647e-15, 5.6428e-15, 9.5139e-15, 1.7322e-14, 2.8829e-14, 
	    4.7708e-14, 6.9789e-14, 9.7267e-14, 1.4662e-13, 1.9429e-13, 
	    2.5998e-13, 3.6636e-13, 4.796e-13, 6.5129e-13, 7.7638e-13, 
	    9.3774e-13, 1.1467e-12, 1.3547e-12, 1.5686e-12, 1.6893e-12, 
	    1.9069e-12, 2.1352e-12, 2.3071e-12, 2.4759e-12, 2.8247e-12, 
	    3.4365e-12, 4.3181e-12, 5.6107e-12, 7.0017e-12, 8.6408e-12, 
	    1.0974e-11, 1.3742e-11, 1.6337e-11, 2.0157e-11, 2.3441e-11, 
	    2.6733e-11, 3.0247e-11, 3.3737e-11, 3.8618e-11, 4.1343e-11, 
	    4.387e-11, 4.4685e-11, 4.4881e-11, 4.5526e-11, 4.3628e-11, 
	    4.4268e-11, 4.6865e-11, 5.3426e-11, 5.402e-11, 5.3218e-11, 
	    5.4587e-11, 5.636e-11, 5.774e-11, 5.6426e-11, 6.0399e-11, 
	    6.6981e-11, 7.4319e-11, 7.7977e-11, 7.5539e-11, 7.161e-11, 
	    6.4606e-11, 5.5498e-11, 4.3944e-11, 3.3769e-11, 2.5771e-11, 
	    1.9162e-11, 1.3698e-11, 1.0173e-11, 7.8925e-12, 6.1938e-12, 
	    4.7962e-12, 4.0811e-12, 3.3912e-12, 2.8625e-12, 2.4504e-12, 
	    2.2188e-12, 2.2139e-12, 2.2499e-12, 2.2766e-12, 2.3985e-12, 
	    2.5459e-12, 2.9295e-12, 3.4196e-12, 3.6155e-12, 4.0733e-12, 
	    4.461e-12, 4.9372e-12, 5.4372e-12, 5.7304e-12, 6.164e-12, 
	    6.1278e-12, 6.294e-12, 6.4947e-12, 6.8174e-12, 7.519e-12, 
	    8.2608e-12, 8.4971e-12, 8.3484e-12, 8.1888e-12, 7.8552e-12, 
	    7.8468e-12, 7.5943e-12, 7.9096e-12, 8.6869e-12, 9.1303e-12, 
	    9.2547e-12, 8.9322e-12, 8.2177e-12, 7.3408e-12, 5.7956e-12, 
	    4.447e-12, 3.5881e-12, 2.6748e-12, 1.7074e-12, 9.67e-13, 
	    5.2645e-13, 2.9943e-13, 1.7316e-13, 1.0039e-13, 5.7859e-14, 
	    3.2968e-14, 1.8499e-14, 1.0192e-14, 5.5015e-15, 2.904e-15, 
	    1.4968e-15, 7.5244e-16, 3.6852e-16, 1.7568e-16, 8.1464e-17, 
	    3.6717e-17, 1.6076e-17, 6.8341e-18, 2.8195e-18, 1.1286e-18, 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 1.407e-18, 3.0405e-18, 6.4059e-18,
	     1.3169e-17, 2.6443e-17, 5.1917e-17, 9.9785e-17, 1.8802e-16, 
	    3.4788e-16, 6.3328e-16, 1.137e-15, 2.0198e-15, 3.5665e-15, 
	    6.3053e-15, 1.1309e-14, 2.1206e-14, 3.2858e-14, 5.5165e-14, 
	    8.6231e-14, 1.2776e-13, 1.778e-13, 2.5266e-13, 3.6254e-13, 
	    5.1398e-13, 6.8289e-13, 8.7481e-13, 1.1914e-12, 1.6086e-12, 
	    2.0469e-12, 2.5761e-12, 3.4964e-12, 4.498e-12, 5.5356e-12, 
	    6.7963e-12, 8.572e-12, 1.07e-11, 1.2983e-11, 1.627e-11, 
	    1.9609e-11, 2.2668e-11, 2.5963e-11, 3.0918e-11, 3.493e-11, 
	    3.933e-11, 4.4208e-11, 4.6431e-11, 5.1141e-11, 5.4108e-11, 
	    5.8077e-11, 6.505e-11, 7.2126e-11, 8.1064e-11, 8.1973e-11, 
	    8.1694e-11, 8.3081e-11, 8.024e-11, 7.9225e-11, 7.6256e-11, 
	    7.8468e-11, 8.0041e-11, 8.1585e-11, 8.3485e-11, 8.3774e-11, 
	    8.587e-11, 8.6104e-11, 8.8516e-11, 9.0814e-11, 9.2522e-11, 
	    8.8913e-11, 7.8381e-11, 6.8568e-11, 5.6797e-11, 4.4163e-11, 
	    3.2369e-11, 2.3259e-11, 1.6835e-11, 1.1733e-11, 8.5273e-12, 
	    6.3805e-12, 4.8983e-12, 3.8831e-12, 3.261e-12, 2.8577e-12, 
	    2.521e-12, 2.2913e-12, 2.0341e-12, 1.8167e-12, 1.6395e-12, 
	    1.489e-12, 1.3516e-12, 1.2542e-12, 1.291e-12, 1.3471e-12, 
	    1.4689e-12, 1.5889e-12, 1.6989e-12, 1.8843e-12, 2.0902e-12, 
	    2.3874e-12, 2.7294e-12, 3.3353e-12, 4.0186e-12, 4.5868e-12, 
	    5.2212e-12, 5.8856e-12, 6.5991e-12, 7.2505e-12, 7.6637e-12, 
	    8.5113e-12, 9.4832e-12, 9.9678e-12, 1.0723e-11, 1.0749e-11, 
	    1.138e-11, 1.1774e-11, 1.1743e-11, 1.2493e-11, 1.2559e-11, 
	    1.2332e-11, 1.1782e-11, 1.1086e-11, 1.0945e-11, 1.1178e-11, 
	    1.2083e-11, 1.3037e-11, 1.473e-11, 1.645e-11, 1.7403e-11, 
	    1.7004e-11, 1.5117e-11, 1.3339e-11, 1.0844e-11, 8.0915e-12, 
	    5.6615e-12, 3.7196e-12, 2.5194e-12, 1.6569e-12, 1.1201e-12, 
	    8.2335e-13, 6.027e-13, 4.8205e-13, 4.1313e-13, 3.6243e-13, 
	    3.2575e-13, 2.773e-13, 2.5292e-13, 2.3062e-13, 2.1126e-13, 
	    2.1556e-13, 2.1213e-13, 2.2103e-13, 2.1927e-13, 2.0794e-13, 
	    1.9533e-13, 1.6592e-13, 1.4521e-13, 1.1393e-13, 8.3772e-14, 
	    6.2077e-14, 4.3337e-14, 2.7165e-14, 1.6821e-14, 9.5407e-15, 
	    5.3093e-15, 3.032e-15, 1.7429e-15, 9.9828e-16, 5.6622e-16, 
	    3.1672e-16, 1.7419e-16, 9.3985e-17, 4.9656e-17, 2.5652e-17, 
	    1.2942e-17, 6.3695e-18, 3.0554e-18, 1.4273e-18, -0., -0., -0., 
	    -0., -0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
	    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };

struct {
    double e_1[3];
    int e_2;
    } fh2ob_ = { -20., 2e4, 10., 2003 };

struct {
    double e_1[2003];
    } sh2oa_ = { .11109, .10573, .10162, .10573, .11109, .12574, .13499, 
	    .14327, .15065, .15164, .15022, .13677, .13115, .12253, .11271, 
	    .1007, .087495, .080118, .06994, .062034, .056051, .047663, 
	    .04245, .03669, .033441, .030711, .025205, .022113, .01888, 
	    .016653, .014626, .012065, .010709, .0091783, .0077274, .0067302, 
	    .0056164, .0049089, .0041497, .0035823, .0031124, .0026414, 
	    .0023167, .0020156, .0017829, .0015666, .0013928, .0012338, 
	    .0010932, 9.7939e-4, 8.8241e-4, 7.9173e-4, 7.1296e-4, 6.4179e-4, 
	    5.8031e-4, 5.2647e-4, 4.7762e-4, 4.3349e-4, 3.9355e-4, 3.5887e-4, 
	    3.2723e-4, 2.9919e-4, 2.7363e-4, 2.5013e-4, 2.2876e-4, 2.0924e-4, 
	    1.9193e-4, 1.7618e-4, 1.6188e-4, 1.4891e-4, 1.3717e-4, 1.2647e-4, 
	    1.1671e-4, 1.0786e-4, 9.9785e-5, 9.235e-5, 8.5539e-5, 7.9377e-5, 
	    7.3781e-5, 6.8677e-5, 6.3993e-5, 5.9705e-5, 5.5788e-5, 5.2196e-5, 
	    4.8899e-5, 4.5865e-5, 4.3079e-5, 4.0526e-5, 3.8182e-5, 3.6025e-5, 
	    3.4038e-5, 3.2203e-5, 3.0511e-5, 2.8949e-5, 2.7505e-5, 2.617e-5, 
	    2.4933e-5, 2.3786e-5, 2.2722e-5, 2.1736e-5, 2.0819e-5, 1.9968e-5, 
	    1.9178e-5, 1.8442e-5, 1.776e-5, 1.7127e-5, 1.6541e-5, 1.5997e-5, 
	    1.5495e-5, 1.5034e-5, 1.4614e-5, 1.423e-5, 1.3883e-5, 1.3578e-5, 
	    1.3304e-5, 1.3069e-5, 1.2876e-5, 1.2732e-5, 1.2626e-5, 1.2556e-5, 
	    1.2544e-5, 1.2604e-5, 1.2719e-5, 1.2883e-5, 1.3164e-5, 1.3581e-5, 
	    1.4187e-5, 1.4866e-5, 1.5669e-5, 1.6717e-5, 1.8148e-5, 2.0268e-5, 
	    2.2456e-5, 2.5582e-5, 2.9183e-5, 3.3612e-5, 3.9996e-5, 4.6829e-5, 
	    5.5055e-5, 6.5897e-5, 7.536e-5, 8.7213e-5, 1.0046e-4, 1.1496e-4, 
	    1.2943e-4, 1.5049e-4, 1.6973e-4, 1.8711e-4, 2.0286e-4, 2.2823e-4, 
	    2.678e-4, 2.8766e-4, 3.1164e-4, 3.364e-4, 3.6884e-4, 3.9159e-4, 
	    3.8712e-4, 3.7433e-4, 3.4503e-4, 3.1003e-4, 2.8027e-4, 2.5253e-4, 
	    2.3408e-4, 2.2836e-4, 2.4442e-4, 2.7521e-4, 2.9048e-4, 3.0489e-4, 
	    3.2646e-4, 3.388e-4, 3.3492e-4, 3.0987e-4, 2.9482e-4, 2.8711e-4, 
	    2.6068e-4, 2.2683e-4, 1.9996e-4, 1.7788e-4, 1.6101e-4, 1.3911e-4, 
	    1.2013e-4, 1.0544e-4, 9.4224e-5, 8.1256e-5, 7.3667e-5, 6.2233e-5, 
	    5.5906e-5, 5.1619e-5, 4.514e-5, 4.0273e-5, 3.3268e-5, 3.0258e-5, 
	    2.644e-5, 2.3103e-5, 2.0749e-5, 1.8258e-5, 1.6459e-5, 1.4097e-5, 
	    1.2052e-5, 1.0759e-5, 9.14e-6, 8.1432e-6, 7.146e-6, 6.4006e-6, 
	    5.6995e-6, 4.9372e-6, 4.4455e-6, 3.9033e-6, 3.474e-6, 3.1269e-6, 
	    2.8059e-6, 2.5558e-6, 2.2919e-6, 2.0846e-6, 1.8983e-6, 1.7329e-6, 
	    1.5929e-6, 1.4631e-6, 1.3513e-6, 1.2461e-6, 1.1519e-6, 1.0682e-6, 
	    9.9256e-7, 9.2505e-7, 8.6367e-7, 8.0857e-7, 7.5674e-7, 7.0934e-7, 
	    6.658e-7, 6.258e-7, 5.8853e-7, 5.5333e-7, 5.2143e-7, 4.9169e-7, 
	    4.6431e-7, 4.3898e-7, 4.1564e-7, 3.9405e-7, 3.7403e-7, 3.5544e-7, 
	    3.3819e-7, 3.2212e-7, 3.0714e-7, 2.9313e-7, 2.8003e-7, 2.6777e-7, 
	    2.5628e-7, 2.4551e-7, 2.354e-7, 2.2591e-7, 2.1701e-7, 2.0866e-7, 
	    2.0082e-7, 1.9349e-7, 1.8665e-7, 1.8027e-7, 1.7439e-7, 1.6894e-7, 
	    1.64e-7, 1.5953e-7, 1.5557e-7, 1.5195e-7, 1.4888e-7, 1.4603e-7, 
	    1.4337e-7, 1.4093e-7, 1.3828e-7, 1.3569e-7, 1.327e-7, 1.2984e-7, 
	    1.2714e-7, 1.2541e-7, 1.2399e-7, 1.2102e-7, 1.1878e-7, 1.1728e-7, 
	    1.1644e-7, 1.1491e-7, 1.1305e-7, 1.1235e-7, 1.1228e-7, 1.1224e-7, 
	    1.1191e-7, 1.1151e-7, 1.1098e-7, 1.1068e-7, 1.1109e-7, 1.1213e-7, 
	    1.1431e-7, 1.1826e-7, 1.2322e-7, 1.3025e-7, 1.4066e-7, 1.5657e-7, 
	    1.7214e-7, 1.9449e-7, 2.2662e-7, 2.6953e-7, 3.1723e-7, 3.7028e-7, 
	    4.4482e-7, 5.3852e-7, 6.2639e-7, 7.2175e-7, 7.7626e-7, 8.7248e-7, 
	    9.6759e-7, 1.0102e-6, 1.062e-6, 1.1201e-6, 1.2107e-6, 1.2998e-6, 
	    1.313e-6, 1.2856e-6, 1.235e-6, 1.1489e-6, 1.0819e-6, 1.012e-6, 
	    9.4795e-7, 9.2858e-7, 9.806e-7, 1.0999e-6, 1.1967e-6, 1.2672e-6, 
	    1.3418e-6, 1.3864e-6, 1.433e-6, 1.4592e-6, 1.4598e-6, 1.4774e-6, 
	    1.4726e-6, 1.482e-6, 1.505e-6, 1.4984e-6, 1.5181e-6, 1.5888e-6, 
	    1.685e-6, 1.769e-6, 1.9277e-6, 2.1107e-6, 2.3068e-6, 2.5347e-6, 
	    2.8069e-6, 3.1345e-6, 3.5822e-6, 3.9051e-6, 4.3422e-6, 4.8704e-6, 
	    5.5351e-6, 6.3454e-6, 7.269e-6, 8.2974e-6, 9.7609e-6, 1.1237e-5, 
	    1.3187e-5, 1.5548e-5, 1.8784e-5, 2.1694e-5, 2.5487e-5, 3.0092e-5, 
	    3.5385e-5, 4.2764e-5, 4.9313e-5, 5.58e-5, 6.2968e-5, 7.106e-5, 
	    7.7699e-5, 8.7216e-5, 8.9335e-5, 9.2151e-5, 9.2779e-5, 9.4643e-5, 
	    9.7978e-5, 1.0008e-4, 1.0702e-4, 1.1026e-4, 1.0828e-4, 1.055e-4, 
	    1.0432e-4, 1.0428e-4, 9.898e-5, 9.4992e-5, 9.5159e-5, 1.0058e-4, 
	    1.0738e-4, 1.155e-4, 1.1229e-4, 1.0596e-4, 1.0062e-4, 9.1742e-5, 
	    8.4492e-5, 6.8099e-5, 5.6295e-5, 4.6502e-5, 3.8071e-5, 3.0721e-5, 
	    2.3297e-5, 1.8688e-5, 1.483e-5, 1.2049e-5, 9.6754e-6, 7.9192e-6, 
	    6.6673e-6, 5.6468e-6, 4.8904e-6, 4.2289e-6, 3.688e-6, 3.2396e-6, 
	    2.8525e-6, 2.5363e-6, 2.2431e-6, 1.9949e-6, 1.7931e-6, 1.6164e-6, 
	    1.4431e-6, 1.2997e-6, 1.1559e-6, 1.0404e-6, 9.43e-7, 8.4597e-7, 
	    7.6133e-7, 6.8623e-7, 6.2137e-7, 5.6345e-7, 5.1076e-7, 4.6246e-7, 
	    4.1906e-7, 3.8063e-7, 3.461e-7, 3.1554e-7, 2.8795e-7, 2.6252e-7, 
	    2.3967e-7, 2.1901e-7, 2.0052e-7, 1.8384e-7, 1.6847e-7, 1.5459e-7, 
	    1.4204e-7, 1.3068e-7, 1.2036e-7, 1.1095e-7, 1.0237e-7, 9.4592e-8, 
	    8.753e-8, 8.1121e-8, 7.5282e-8, 6.9985e-8, 6.5189e-8, 6.0874e-8, 
	    5.6989e-8, 5.353e-8, 5.0418e-8, 4.7745e-8, 4.5367e-8, 4.3253e-8, 
	    4.1309e-8, 3.9695e-8, 3.8094e-8, 3.6482e-8, 3.4897e-8, 3.35e-8, 
	    3.2302e-8, 3.0854e-8, 2.9698e-8, 2.8567e-8, 2.76e-8, 2.6746e-8, 
	    2.5982e-8, 2.551e-8, 2.5121e-8, 2.4922e-8, 2.4909e-8, 2.5013e-8, 
	    2.5216e-8, 2.5589e-8, 2.6049e-8, 2.6451e-8, 2.6978e-8, 2.7687e-8, 
	    2.86e-8, 2.9643e-8, 3.0701e-8, 3.2058e-8, 3.3695e-8, 3.5558e-8, 
	    3.7634e-8, 3.9875e-8, 4.2458e-8, 4.548e-8, 4.8858e-8, 5.2599e-8, 
	    5.703e-8, 6.2067e-8, 6.7911e-8, 7.4579e-8, 8.1902e-8, 8.9978e-8, 
	    9.987e-8, 1.1102e-7, 1.2343e-7, 1.3732e-7, 1.5394e-7, 1.7318e-7, 
	    1.9383e-7, 2.1819e-7, 2.4666e-7, 2.8109e-7, 3.2236e-7, 3.776e-7, 
	    4.4417e-7, 5.2422e-7, 6.1941e-7, 7.4897e-7, 9.2041e-7, 1.1574e-6, 
	    1.4126e-6, 1.7197e-6, 2.1399e-6, 2.6266e-6, 3.3424e-6, 3.8418e-6, 
	    4.514e-6, 5.0653e-6, 5.8485e-6, 6.5856e-6, 6.8937e-6, 6.9121e-6, 
	    6.9005e-6, 6.9861e-6, 6.82e-6, 6.6089e-6, 6.5809e-6, 7.3496e-6, 
	    8.0311e-6, 8.3186e-6, 8.426e-6, 9.0644e-6, 9.4965e-6, 9.4909e-6, 
	    9.016e-6, 9.1494e-6, 9.3629e-6, 9.5944e-6, 9.5459e-6, 8.9919e-6, 
	    8.604e-6, 7.8613e-6, 7.1567e-6, 6.2677e-6, 5.1899e-6, 4.4188e-6, 
	    3.7167e-6, 3.0636e-6, 2.5573e-6, 2.0317e-6, 1.6371e-6, 1.3257e-6, 
	    1.0928e-6, 8.9986e-7, 7.4653e-7, 6.1111e-7, 5.1395e-7, 4.35e-7, 
	    3.7584e-7, 3.2633e-7, 2.8413e-7, 2.4723e-7, 2.1709e-7, 1.9294e-7, 
	    1.7258e-7, 1.5492e-7, 1.382e-7, 1.2389e-7, 1.1189e-7, 1.0046e-7, 
	    9.0832e-8, 8.2764e-8, 7.4191e-8, 6.7085e-8, 6.0708e-8, 5.4963e-8, 
	    4.9851e-8, 4.5044e-8, 4.0916e-8, 3.722e-8, 3.3678e-8, 3.0663e-8, 
	    2.7979e-8, 2.5495e-8, 2.3286e-8, 2.1233e-8, 1.9409e-8, 1.777e-8, 
	    1.626e-8, 1.4885e-8, 1.3674e-8, 1.2543e-8, 1.1551e-8, 1.0655e-8, 
	    9.8585e-9, 9.1398e-9, 8.4806e-9, 7.8899e-9, 7.3547e-9, 6.867e-9, 
	    6.4131e-9, 5.993e-9, 5.6096e-9, 5.2592e-9, 4.9352e-9, 4.6354e-9, 
	    4.3722e-9, 4.125e-9, 3.9081e-9, 3.7118e-9, 3.5372e-9, 3.3862e-9, 
	    3.2499e-9, 3.1324e-9, 3.0313e-9, 2.9438e-9, 2.8686e-9, 2.805e-9, 
	    2.7545e-9, 2.7149e-9, 2.6907e-9, 2.6724e-9, 2.6649e-9, 2.6642e-9, 
	    2.6725e-9, 2.6871e-9, 2.7056e-9, 2.7357e-9, 2.7781e-9, 2.8358e-9, 
	    2.9067e-9, 2.9952e-9, 3.102e-9, 3.2253e-9, 3.3647e-9, 3.5232e-9, 
	    3.7037e-9, 3.9076e-9, 4.1385e-9, 4.3927e-9, 4.6861e-9, 5.0238e-9, 
	    5.4027e-9, 5.8303e-9, 6.3208e-9, 6.8878e-9, 7.5419e-9, 8.313e-9, 
	    9.1952e-9, 1.0228e-8, 1.1386e-8, 1.2792e-8, 1.4521e-8, 1.6437e-8, 
	    1.8674e-8, 2.116e-8, 2.4506e-8, 2.8113e-8, 3.2636e-8, 3.7355e-8, 
	    4.2234e-8, 4.9282e-8, 5.7358e-8, 6.6743e-8, 7.8821e-8, 9.4264e-8, 
	    1.1542e-7, 1.3684e-7, 1.6337e-7, 2.0056e-7, 2.3252e-7, 2.6127e-7, 
	    2.9211e-7, 3.3804e-7, 3.7397e-7, 3.8205e-7, 3.881e-7, 3.9499e-7, 
	    3.9508e-7, 3.7652e-7, 3.5859e-7, 3.6198e-7, 3.7871e-7, 4.0925e-7, 
	    4.2717e-7, 4.8241e-7, 5.2008e-7, 5.653e-7, 5.9531e-7, 6.1994e-7, 
	    6.508e-7, 6.6355e-7, 6.9193e-7, 6.993e-7, 7.3058e-7, 7.4678e-7, 
	    7.9193e-7, 8.3627e-7, 9.1267e-7, 1.0021e-6, 1.1218e-6, 1.2899e-6, 
	    1.4447e-6, 1.7268e-6, 2.0025e-6, 2.3139e-6, 2.5599e-6, 2.892e-6, 
	    3.3059e-6, 3.5425e-6, 3.9522e-6, 4.0551e-6, 4.2818e-6, 4.2892e-6, 
	    4.421e-6, 4.5614e-6, 4.6739e-6, 4.9482e-6, 5.1118e-6, 5.0986e-6, 
	    4.9417e-6, 4.9022e-6, 4.8449e-6, 4.8694e-6, 4.8111e-6, 4.9378e-6, 
	    5.3231e-6, 5.7362e-6, 6.235e-6, 6.0951e-6, 5.7281e-6, 5.4585e-6, 
	    4.9032e-6, 4.3009e-6, 3.4776e-6, 2.8108e-6, 2.2993e-6, 1.7999e-6, 
	    1.387e-6, 1.075e-6, 8.5191e-7, 6.7951e-7, 5.5336e-7, 4.6439e-7, 
	    4.0243e-7, 3.5368e-7, 3.1427e-7, 2.7775e-7, 2.4486e-7, 2.1788e-7, 
	    1.9249e-7, 1.7162e-7, 1.5115e-7, 1.3478e-7, 1.2236e-7, 1.1139e-7, 
	    1.0092e-7, 9.0795e-8, 8.2214e-8, 7.4691e-8, 6.7486e-8, 6.0414e-8, 
	    5.4584e-8, 4.8754e-8, 4.3501e-8, 3.8767e-8, 3.4363e-8, 3.0703e-8, 
	    2.7562e-8, 2.4831e-8, 2.2241e-8, 1.9939e-8, 1.8049e-8, 1.6368e-8, 
	    1.4863e-8, 1.346e-8, 1.2212e-8, 1.1155e-8, 1.0185e-8, 9.3417e-9, 
	    8.5671e-9, 7.8292e-9, 7.1749e-9, 6.5856e-9, 6.0588e-9, 5.5835e-9, 
	    5.135e-9, 4.7395e-9, 4.3771e-9, 4.0476e-9, 3.756e-9, 3.4861e-9, 
	    3.2427e-9, 3.024e-9, 2.8278e-9, 2.6531e-9, 2.4937e-9, 2.3511e-9, 
	    2.2245e-9, 2.1133e-9, 2.0159e-9, 1.933e-9, 1.8669e-9, 1.8152e-9, 
	    1.7852e-9, 1.7752e-9, 1.7823e-9, 1.8194e-9, 1.8866e-9, 1.9759e-9, 
	    2.0736e-9, 2.2083e-9, 2.3587e-9, 2.4984e-9, 2.6333e-9, 2.816e-9, 
	    3.0759e-9, 3.372e-9, 3.6457e-9, 4.0668e-9, 4.4541e-9, 4.7976e-9, 
	    5.0908e-9, 5.4811e-9, 6.1394e-9, 6.3669e-9, 6.5714e-9, 6.8384e-9, 
	    7.1918e-9, 7.3741e-9, 7.2079e-9, 7.2172e-9, 7.2572e-9, 7.3912e-9, 
	    7.6188e-9, 8.3291e-9, 8.7885e-9, 9.2412e-9, 1.0021e-8, 1.0752e-8, 
	    1.1546e-8, 1.1607e-8, 1.1949e-8, 1.2346e-8, 1.2516e-8, 1.2826e-8, 
	    1.3053e-8, 1.3556e-8, 1.4221e-8, 1.5201e-8, 1.6661e-8, 1.8385e-8, 
	    2.0585e-8, 2.3674e-8, 2.7928e-8, 3.3901e-8, 4.1017e-8, 4.9595e-8, 
	    6.0432e-8, 7.6304e-8, 9.0764e-8, 1.0798e-7, 1.2442e-7, 1.4404e-7, 
	    1.6331e-7, 1.8339e-7, 2.0445e-7, 2.2288e-7, 2.3083e-7, 2.3196e-7, 
	    2.3919e-7, 2.3339e-7, 2.3502e-7, 2.3444e-7, 2.6395e-7, 2.9928e-7, 
	    3.0025e-7, 3.0496e-7, 3.1777e-7, 3.4198e-7, 3.4739e-7, 3.2696e-7, 
	    3.41e-7, 3.5405e-7, 3.7774e-7, 3.8285e-7, 3.6797e-7, 3.58e-7, 
	    3.2283e-7, 2.9361e-7, 2.4881e-7, 2.0599e-7, 1.7121e-7, 1.3641e-7, 
	    1.1111e-7, 8.9413e-8, 7.3455e-8, 6.2078e-8, 5.2538e-8, 4.5325e-8, 
	    3.9005e-8, 3.4772e-8, 3.1203e-8, 2.8132e-8, 2.525e-8, 2.2371e-8, 
	    2.0131e-8, 1.7992e-8, 1.6076e-8, 1.4222e-8, 1.249e-8, 1.1401e-8, 
	    1.0249e-8, 9.2279e-9, 8.5654e-9, 7.6227e-9, 6.9648e-9, 6.2466e-9, 
	    5.7252e-9, 5.38e-9, 4.696e-9, 4.2194e-9, 3.7746e-9, 3.3813e-9, 
	    3.0656e-9, 2.6885e-9, 2.4311e-9, 2.1572e-9, 1.8892e-9, 1.7038e-9, 
	    1.4914e-9, 1.3277e-9, 1.1694e-9, 1.0391e-9, 9.2779e-10, 
	    8.3123e-10, 7.4968e-10, 6.8385e-10, 6.2915e-10, 5.7784e-10, 
	    5.2838e-10, 4.8382e-10, 4.4543e-10, 4.1155e-10, 3.7158e-10, 
	    3.3731e-10, 3.0969e-10, 2.8535e-10, 2.6416e-10, 2.4583e-10, 
	    2.2878e-10, 2.1379e-10, 2.0073e-10, 1.8907e-10, 1.7866e-10, 
	    1.6936e-10, 1.6119e-10, 1.5424e-10, 1.4847e-10, 1.4401e-10, 
	    1.4068e-10, 1.3937e-10, 1.3943e-10, 1.4281e-10, 1.4766e-10, 
	    1.5701e-10, 1.7079e-10, 1.8691e-10, 2.0081e-10, 2.174e-10, 
	    2.4847e-10, 2.6463e-10, 2.7087e-10, 2.7313e-10, 2.8352e-10, 
	    2.9511e-10, 2.8058e-10, 2.7227e-10, 2.7356e-10, 2.8012e-10, 
	    2.8034e-10, 2.9031e-10, 3.103e-10, 3.3745e-10, 3.8152e-10, 
	    4.0622e-10, 4.2673e-10, 4.3879e-10, 4.5488e-10, 4.7179e-10, 
	    4.614e-10, 4.6339e-10, 4.6716e-10, 4.7024e-10, 4.7931e-10, 
	    4.8503e-10, 4.9589e-10, 4.9499e-10, 5.0363e-10, 5.3184e-10, 
	    5.6451e-10, 6.0932e-10, 6.6469e-10, 7.4076e-10, 8.3605e-10, 
	    9.4898e-10, 1.0935e-9, 1.2593e-9, 1.4913e-9, 1.8099e-9, 2.1842e-9,
	     2.7284e-9, 3.2159e-9, 3.7426e-9, 4.5226e-9, 5.3512e-9, 6.1787e-9,
	     6.8237e-9, 7.9421e-9, 9.0002e-9, 9.6841e-9, 9.9558e-9, 1.0232e-8,
	     1.0591e-8, 1.0657e-8, 1.0441e-8, 1.0719e-8, 1.1526e-8, 1.2962e-8,
	     1.4336e-8, 1.615e-8, 1.8417e-8, 2.0725e-8, 2.3426e-8, 2.5619e-8, 
	    2.7828e-8, 3.0563e-8, 3.3438e-8, 3.6317e-8, 4.04e-8, 4.4556e-8, 
	    5.0397e-8, 5.3315e-8, 5.9185e-8, 6.5311e-8, 6.9188e-8, 7.7728e-8, 
	    7.9789e-8, 8.6598e-8, 8.7768e-8, 9.1773e-8, 9.7533e-8, 1.0007e-7, 
	    1.065e-7, 1.0992e-7, 1.0864e-7, 1.0494e-7, 1.0303e-7, 1.0031e-7, 
	    1.0436e-7, 1.0537e-7, 1.1184e-7, 1.2364e-7, 1.3651e-7, 1.4881e-7, 
	    1.4723e-7, 1.4118e-7, 1.3371e-7, 1.1902e-7, 1.0007e-7, 7.9628e-8, 
	    6.4362e-8, 5.0243e-8, 3.8133e-8, 2.94e-8, 2.3443e-8, 1.9319e-8, 
	    1.6196e-8, 1.4221e-8, 1.2817e-8, 1.1863e-8, 1.1383e-8, 1.1221e-8, 
	    1.1574e-8, 1.1661e-8, 1.2157e-8, 1.2883e-8, 1.3295e-8, 1.4243e-8, 
	    1.424e-8, 1.4614e-8, 1.4529e-8, 1.4685e-8, 1.4974e-8, 1.479e-8, 
	    1.489e-8, 1.4704e-8, 1.4142e-8, 1.3374e-8, 1.2746e-8, 1.2172e-8, 
	    1.2336e-8, 1.2546e-8, 1.3065e-8, 1.409e-8, 1.5215e-8, 1.654e-8, 
	    1.6144e-8, 1.5282e-8, 1.4358e-8, 1.2849e-8, 1.0998e-8, 8.6956e-9, 
	    7.0881e-9, 5.5767e-9, 4.2792e-9, 3.2233e-9, 2.502e-9, 1.9985e-9, 
	    1.5834e-9, 1.3015e-9, 1.0948e-9, 9.4141e-10, 8.1465e-10, 
	    7.1517e-10, 6.2906e-10, 5.5756e-10, 4.9805e-10, 4.3961e-10, 
	    3.9181e-10, 3.5227e-10, 3.167e-10, 2.8667e-10, 2.5745e-10, 
	    2.3212e-10, 2.0948e-10, 1.897e-10, 1.7239e-10, 1.5659e-10, 
	    1.4301e-10, 1.3104e-10, 1.2031e-10, 1.1095e-10, 1.0262e-10, 
	    9.513e-11, 8.8595e-11, 8.2842e-11, 7.7727e-11, 7.3199e-11, 
	    6.9286e-11, 6.5994e-11, 6.3316e-11, 6.1244e-11, 5.9669e-11, 
	    5.8843e-11, 5.8832e-11, 5.9547e-11, 6.1635e-11, 6.4926e-11, 
	    7.0745e-11, 7.8802e-11, 8.6724e-11, 1.0052e-10, 1.1575e-10, 
	    1.3626e-10, 1.5126e-10, 1.6751e-10, 1.9239e-10, 2.1748e-10, 
	    2.2654e-10, 2.2902e-10, 2.324e-10, 2.4081e-10, 2.393e-10, 
	    2.2378e-10, 2.2476e-10, 2.2791e-10, 2.4047e-10, 2.5305e-10, 
	    2.8073e-10, 3.1741e-10, 3.6592e-10, 4.1495e-10, 4.6565e-10, 
	    5.099e-10, 5.5607e-10, 6.1928e-10, 6.6779e-10, 7.335e-10, 
	    8.1434e-10, 8.9635e-10, 9.9678e-10, 1.1256e-9, 1.2999e-9, 
	    1.4888e-9, 1.7642e-9, 1.9606e-9, 2.2066e-9, 2.4601e-9, 2.7218e-9, 
	    3.0375e-9, 3.1591e-9, 3.2852e-9, 3.2464e-9, 3.3046e-9, 3.271e-9, 
	    3.2601e-9, 3.3398e-9, 3.7446e-9, 4.0795e-9, 4.0284e-9, 4.0584e-9, 
	    4.1677e-9, 4.5358e-9, 4.4097e-9, 4.2744e-9, 4.5449e-9, 4.8147e-9, 
	    5.2656e-9, 5.2476e-9, 5.0275e-9, 4.7968e-9, 4.3654e-9, 3.953e-9, 
	    3.2447e-9, 2.6489e-9, 2.1795e-9, 1.788e-9, 1.4309e-9, 1.1256e-9, 
	    9.1903e-10, 7.6533e-10, 6.3989e-10, 5.5496e-10, 4.9581e-10, 
	    4.5722e-10, 4.3898e-10, 4.3505e-10, 4.3671e-10, 4.5329e-10, 
	    4.6827e-10, 4.9394e-10, 5.1122e-10, 5.1649e-10, 5.0965e-10, 
	    4.9551e-10, 4.8928e-10, 4.7947e-10, 4.7989e-10, 4.9071e-10, 
	    4.8867e-10, 4.726e-10, 4.5756e-10, 4.54e-10, 4.5993e-10, 
	    4.4042e-10, 4.3309e-10, 4.4182e-10, 4.6735e-10, 5.0378e-10, 
	    5.2204e-10, 5.0166e-10, 4.6799e-10, 4.3119e-10, 3.8803e-10, 
	    3.3291e-10, 2.6289e-10, 2.1029e-10, 1.7011e-10, 1.3345e-10, 
	    1.0224e-10, 7.8207e-11, 6.2451e-11, 5.0481e-11, 4.1507e-11, 
	    3.5419e-11, 3.0582e-11, 2.69e-11, 2.3778e-11, 2.1343e-11, 
	    1.9182e-11, 1.7162e-11, 1.5391e-11, 1.3877e-11, 1.2619e-11, 
	    1.145e-11, 1.0461e-11, 9.6578e-12, 8.9579e-12, 8.3463e-12, 
	    7.8127e-12, 7.3322e-12, 6.9414e-12, 6.6037e-12, 6.3285e-12, 
	    6.1095e-12, 5.9387e-12, 5.8118e-12, 5.726e-12, 5.6794e-12, 
	    5.6711e-12, 5.7003e-12, 5.767e-12, 5.8717e-12, 6.0151e-12, 
	    6.1984e-12, 6.4232e-12, 6.6918e-12, 7.0065e-12, 7.3705e-12, 
	    7.7873e-12, 8.2612e-12, 8.7972e-12, 9.4009e-12, 1.0079e-11, 
	    1.084e-11, 1.1692e-11, 1.2648e-11, 1.3723e-11, 1.4935e-11, 
	    1.6313e-11, 1.7905e-11, 1.974e-11, 2.1898e-11, 2.4419e-11, 
	    2.7426e-11, 3.0869e-11, 3.4235e-11, 3.7841e-11, 4.1929e-11, 
	    4.6776e-11, 5.2123e-11, 5.8497e-11, 6.5294e-11, 7.4038e-11, 
	    8.4793e-11, 9.6453e-11, 1.1223e-10, 1.2786e-10, 1.4882e-10, 
	    1.7799e-10, 2.0766e-10, 2.4523e-10, 2.8591e-10, 3.3386e-10, 
	    4.0531e-10, 4.7663e-10, 5.4858e-10, 6.3377e-10, 7.1688e-10, 
	    8.4184e-10, 9.5144e-10, 1.0481e-9, 1.1356e-9, 1.2339e-9, 
	    1.3396e-9, 1.4375e-9, 1.5831e-9, 1.7323e-9, 1.9671e-9, 2.2976e-9, 
	    2.6679e-9, 3.0777e-9, 3.4321e-9, 3.8192e-9, 4.2711e-9, 4.4903e-9, 
	    4.8931e-9, 5.2253e-9, 5.404e-9, 5.6387e-9, 5.6704e-9, 6.0345e-9, 
	    6.1079e-9, 6.2576e-9, 6.4039e-9, 6.3776e-9, 6.1878e-9, 5.8616e-9, 
	    5.7036e-9, 5.584e-9, 5.6905e-9, 5.8931e-9, 6.2478e-9, 6.8291e-9, 
	    7.4528e-9, 7.6078e-9, 7.3898e-9, 6.7573e-9, 5.9827e-9, 5.0927e-9, 
	    4.0099e-9, 3.1933e-9, 2.4296e-9, 1.8485e-9, 1.4595e-9, 1.2017e-9, 
	    1.0164e-9, 8.7433e-10, 7.7108e-10, 7.0049e-10, 6.5291e-10, 
	    6.1477e-10, 5.9254e-10, 5.815e-10, 5.7591e-10, 5.849e-10, 
	    5.8587e-10, 5.9636e-10, 6.2408e-10, 6.5479e-10, 7.048e-10, 
	    7.2313e-10, 7.5524e-10, 8.0863e-10, 8.3386e-10, 9.2342e-10, 
	    9.6754e-10, 1.0293e-9, 1.0895e-9, 1.133e-9, 1.221e-9, 1.2413e-9, 
	    1.2613e-9, 1.2671e-9, 1.2225e-9, 1.1609e-9, 1.0991e-9, 1.06e-9, 
	    1.057e-9, 1.0818e-9, 1.1421e-9, 1.227e-9, 1.337e-9, 1.4742e-9, 
	    1.4946e-9, 1.4322e-9, 1.321e-9, 1.1749e-9, 1.0051e-9, 7.8387e-10, 
	    6.1844e-10, 4.6288e-10, 3.4164e-10, 2.5412e-10, 1.9857e-10, 
	    1.5876e-10, 1.2966e-10, 1.092e-10, 9.4811e-11, 8.3733e-11, 
	    7.3906e-11, 6.7259e-11, 6.1146e-11, 5.7119e-11, 5.3546e-11, 
	    4.8625e-11, 4.4749e-11, 4.1089e-11, 3.7825e-11, 3.4465e-11, 
	    3.1018e-11, 2.8109e-11, 2.561e-11, 2.2859e-11, 2.049e-11, 
	    1.8133e-11, 1.5835e-11, 1.3949e-11, 1.2295e-11, 1.0799e-11, 
	    9.6544e-12, 8.7597e-12, 7.999e-12, 7.3973e-12, 6.9035e-12, 
	    6.4935e-12, 6.1195e-12, 5.8235e-12, 5.5928e-12, 5.4191e-12, 
	    5.2993e-12, 5.2338e-12, 5.2272e-12, 5.2923e-12, 5.4252e-12, 
	    5.6523e-12, 5.9433e-12, 6.3197e-12, 6.9016e-12, 7.5016e-12, 
	    8.2885e-12, 9.405e-12, 1.0605e-11, 1.2257e-11, 1.3622e-11, 
	    1.5353e-11, 1.7543e-11, 1.9809e-11, 2.2197e-11, 2.4065e-11, 
	    2.6777e-11, 2.9751e-11, 3.2543e-11, 3.5536e-11, 3.9942e-11, 
	    4.6283e-11, 5.4556e-11, 6.549e-11, 7.6803e-11, 9.0053e-11, 
	    1.0852e-10, 1.2946e-10, 1.4916e-10, 1.7748e-10, 2.0073e-10, 
	    2.2485e-10, 2.5114e-10, 2.7715e-10, 3.1319e-10, 3.3305e-10, 
	    3.5059e-10, 3.5746e-10, 3.6311e-10, 3.7344e-10, 3.6574e-10, 
	    3.7539e-10, 3.9434e-10, 4.351e-10, 4.334e-10, 4.2588e-10, 
	    4.3977e-10, 4.6062e-10, 4.7687e-10, 4.6457e-10, 4.8578e-10, 
	    5.2344e-10, 5.6752e-10, 5.8702e-10, 5.6603e-10, 5.3784e-10, 
	    4.9181e-10, 4.3272e-10, 3.5681e-10, 2.8814e-10, 2.332e-10, 
	    1.8631e-10, 1.4587e-10, 1.1782e-10, 9.8132e-11, 8.2528e-11, 
	    6.9174e-11, 6.1056e-11, 5.3459e-11, 4.7116e-11, 4.1878e-11, 
	    3.8125e-11, 3.6347e-11, 3.5071e-11, 3.3897e-11, 3.3541e-11, 
	    3.3563e-11, 3.5469e-11, 3.8111e-11, 3.8675e-11, 4.1333e-11, 
	    4.3475e-11, 4.6476e-11, 4.9761e-11, 5.138e-11, 5.4135e-11, 
	    5.3802e-11, 5.5158e-11, 5.6864e-11, 5.9311e-11, 6.3827e-11, 
	    6.7893e-11, 6.823e-11, 6.6694e-11, 6.6018e-11, 6.4863e-11, 
	    6.5893e-11, 6.3813e-11, 6.4741e-11, 6.863e-11, 7.0255e-11, 
	    7.0667e-11, 6.881e-11, 6.4104e-11, 5.8136e-11, 4.7242e-11, 
	    3.7625e-11, 3.1742e-11, 2.5581e-11, 1.8824e-11, 1.3303e-11, 
	    9.6919e-12, 7.5353e-12, 6.0986e-12, 5.0742e-12, 4.3094e-12, 
	    3.719e-12, 3.252e-12, 2.8756e-12, 2.568e-12, 2.3139e-12, 
	    2.1025e-12, 1.9257e-12, 1.7777e-12, 1.6539e-12, 1.5508e-12, 
	    1.4657e-12, 1.3966e-12, 1.3417e-12, 1.2998e-12, 1.27e-12, 
	    1.2514e-12, 1.2437e-12, 1.2463e-12, 1.2592e-12, 1.2823e-12, 
	    1.3157e-12, 1.3596e-12, 1.4144e-12, 1.4806e-12, 1.5588e-12, 
	    1.6497e-12, 1.7544e-12, 1.8738e-12, 2.0094e-12, 2.1626e-12, 
	    2.3354e-12, 2.5297e-12, 2.7483e-12, 2.9941e-12, 3.2708e-12, 
	    3.5833e-12, 3.9374e-12, 4.3415e-12, 4.8079e-12, 5.3602e-12, 
	    5.9816e-12, 6.7436e-12, 7.6368e-12, 8.6812e-12, 9.8747e-12, 
	    1.135e-11, 1.3181e-11, 1.5406e-11, 1.7868e-11, 2.0651e-11, 
	    2.4504e-11, 2.9184e-11, 3.4159e-11, 3.9979e-11, 4.8704e-11, 
	    5.7856e-11, 6.7576e-11, 7.9103e-11, 9.437e-11, 1.1224e-10, 
	    1.3112e-10, 1.5674e-10, 1.8206e-10, 2.0576e-10, 2.3187e-10, 
	    2.7005e-10, 3.0055e-10, 3.3423e-10, 3.6956e-10, 3.8737e-10, 
	    4.263e-10, 4.5154e-10, 4.8383e-10, 5.3582e-10, 5.8109e-10, 
	    6.3741e-10, 6.3874e-10, 6.387e-10, 6.5818e-10, 6.5056e-10, 
	    6.5291e-10, 6.3159e-10, 6.3984e-10, 6.4549e-10, 6.5444e-10, 
	    6.7035e-10, 6.7665e-10, 6.9124e-10, 6.8451e-10, 6.9255e-10, 
	    6.9923e-10, 7.0396e-10, 6.7715e-10, 6.0371e-10, 5.3774e-10, 
	    4.6043e-10, 3.7635e-10, 2.9484e-10, 2.2968e-10, 1.8185e-10, 
	    1.4191e-10, 1.1471e-10, 9.479e-11, 7.9613e-11, 6.7989e-11, 
	    5.9391e-11, 5.281e-11, 4.7136e-11, 4.2618e-11, 3.8313e-11, 
	    3.4686e-11, 3.1669e-11, 2.911e-11, 2.6871e-11, 2.5074e-11, 
	    2.4368e-11, 2.3925e-11, 2.4067e-11, 2.4336e-11, 2.4704e-11, 
	    2.5823e-11, 2.7177e-11, 2.9227e-11, 3.1593e-11, 3.573e-11, 
	    4.0221e-11, 4.3994e-11, 4.8448e-11, 5.3191e-11, 5.8552e-11, 
	    6.3458e-11, 6.6335e-11, 7.2457e-11, 7.9091e-11, 8.2234e-11, 
	    8.7668e-11, 8.7951e-11, 9.2952e-11, 9.6157e-11, 9.5926e-11, 
	    1.012e-10, 1.0115e-10, 9.9577e-11, 9.6633e-11, 9.2891e-11, 
	    9.3315e-11, 9.5584e-11, 1.0064e-10, 1.0509e-10, 1.1455e-10, 
	    1.2443e-10, 1.2963e-10, 1.2632e-10, 1.1308e-10, 1.0186e-10, 
	    8.588e-11, 6.7863e-11, 5.1521e-11, 3.778e-11, 2.8842e-11, 
	    2.2052e-11, 1.7402e-11, 1.4406e-11, 1.1934e-11, 1.0223e-11, 
	    8.9544e-12, 7.9088e-12, 7.0675e-12, 6.2222e-12, 5.6051e-12, 
	    5.0502e-12, 4.5578e-12, 4.2636e-12, 3.9461e-12, 3.7599e-12, 
	    3.5215e-12, 3.2467e-12, 3.0018e-12, 2.6558e-12, 2.3928e-12, 
	    2.0707e-12, 1.7575e-12, 1.5114e-12, 1.2941e-12, 1.1004e-12, 
	    9.5175e-13, 8.2894e-13, 7.3253e-13, 6.5551e-13, 5.9098e-13, 
	    5.3548e-13, 4.8697e-13, 4.4413e-13, 4.06e-13, 3.7188e-13, 
	    3.4121e-13, 3.1356e-13, 2.8856e-13, 2.659e-13, 2.4533e-13, 
	    2.2663e-13, 2.096e-13, 1.9407e-13, 1.799e-13, 1.6695e-13, 
	    1.5512e-13, 1.4429e-13, 1.3437e-13, 1.2527e-13, 1.1693e-13, 
	    1.0927e-13, 1.0224e-13, 9.5767e-14, 8.9816e-14, 8.4335e-14, 
	    7.9285e-14, 7.4626e-14, 7.0325e-14, 6.6352e-14, 6.2676e-14, 
	    5.9274e-14, 5.6121e-14, 5.3195e-14, 5.0479e-14, 4.7953e-14, 
	    4.5602e-14, 4.3411e-14, 4.1367e-14, 3.9456e-14, 3.767e-14, 
	    3.5996e-14, 3.4427e-14, 3.2952e-14, 3.1566e-14, 3.0261e-14, 
	    2.903e-14, 2.7868e-14, 2.677e-14, 2.573e-14, 2.4745e-14, 
	    2.3809e-14, 2.2921e-14, 2.2076e-14, 2.1271e-14, 2.0504e-14, 
	    1.9772e-14, 1.9073e-14, 1.8404e-14, 1.7764e-14, 1.7151e-14, 
	    1.6564e-14, 1.6e-14, 1.5459e-14, 1.4939e-14, 1.4439e-14, 
	    1.3958e-14, 1.3495e-14, 1.3049e-14, 1.262e-14, 1.2206e-14, 
	    1.1807e-14, 1.1422e-14, 1.105e-14, 1.0691e-14, 1.0345e-14, 
	    1.001e-14, 9.687e-15, 9.3747e-15, 9.0727e-15, 8.7808e-15, 
	    8.4986e-15, 8.2257e-15, 7.9617e-15, 7.7064e-15, 7.4594e-15, 
	    7.2204e-15, 6.9891e-15, 6.7653e-15, 6.5488e-15, 6.3392e-15, 
	    6.1363e-15, 5.9399e-15, 5.7499e-15, 5.5659e-15, 5.3878e-15, 
	    5.2153e-15, 5.0484e-15, 4.8868e-15, 4.7303e-15, 4.5788e-15, 
	    4.4322e-15, 4.2902e-15, 4.1527e-15, 4.0196e-15, 3.8907e-15, 
	    3.7659e-15, 3.6451e-15, 3.5281e-15, 3.4149e-15, 3.3052e-15, 
	    3.1991e-15, 3.0963e-15, 2.9967e-15, 2.9004e-15, 2.8071e-15, 
	    2.7167e-15, 2.6293e-15, 2.5446e-15, 2.4626e-15, 2.3833e-15, 
	    2.3064e-15, 2.232e-15, 2.16e-15, 2.0903e-15, 2.0228e-15, 
	    1.9574e-15, 1.8942e-15, 1.8329e-15, 1.7736e-15, 1.7163e-15, 
	    1.6607e-15, 1.6069e-15, 1.5548e-15, 1.5044e-15, 1.4557e-15, 
	    1.4084e-15, 1.3627e-15, 1.3185e-15, 1.2757e-15, 1.2342e-15, 
	    1.1941e-15, 1.1552e-15, 1.1177e-15, 1.0813e-15, 1.0461e-15, 
	    1.012e-15, 9.79e-16, 9.4707e-16, 9.1618e-16, 8.8628e-16, 
	    8.5734e-16, 8.2933e-16, 8.0223e-16, 7.76e-16, 7.5062e-16, 
	    7.2606e-16, 7.0229e-16, 6.7929e-16, 6.5703e-16, 6.355e-16, 
	    6.1466e-16, 5.9449e-16, 5.7498e-16, 5.561e-16, 5.3783e-16, 
	    5.2015e-16, 5.0305e-16, 4.865e-16, 4.7049e-16, 4.55e-16, 
	    4.4002e-16, 4.2552e-16, 4.1149e-16, 3.9792e-16, 3.8479e-16, 
	    3.7209e-16, 3.5981e-16, 3.4792e-16, 3.3642e-16, 3.253e-16, 
	    3.1454e-16, 3.0413e-16, 2.9406e-16, 2.8432e-16, 2.749e-16, 
	    2.6579e-16, 2.5697e-16, 2.4845e-16, 2.402e-16, 2.3223e-16, 
	    2.2451e-16, 2.1705e-16, 2.0984e-16, 2.0286e-16, 1.9611e-16, 
	    1.8958e-16, 1.8327e-16, 1.7716e-16, 1.7126e-16, 1.6555e-16, 
	    1.6003e-16, 1.5469e-16, 1.4952e-16, 1.4453e-16, 1.397e-16, 
	    1.3503e-16 };

struct {
    double e_1[3];
    int e_2;
    } sh2ob_ = { -20., 2e4, 10., 2003 };

struct {
    double e_1[2003];
    } s260a_ = { .1775, .17045, .16457, .17045, .1775, .20036, .21347, .22454,
	     .23428, .23399, .23022, .20724, .19712, .18317, .16724, .1478, 
	    .12757, .11626, .10098, .089033, .07977, .067416, .059588, 
	    .051117, .046218, .042179, .034372, .029863, .025252, .022075, 
	    .019209, .015816, .013932, .011943, .010079, .0087667, .0074094, 
	    .0064967, .0055711, .0048444, .0042552, .0036953, .0032824, 
	    .0029124, .0026102, .002337, .00211, .0019008, .0017145, .0015573,
	     .0014206, .0012931, .0011803, .0010774, 9.8616e-4, 9.0496e-4, 
	    8.3071e-4, 7.6319e-4, 7.0149e-4, 6.4637e-4, 5.9566e-4, 5.4987e-4, 
	    5.0768e-4, 4.688e-4, 4.3317e-4, 4.0037e-4, 3.7064e-4, 3.4325e-4, 
	    3.1809e-4, 2.9501e-4, 2.7382e-4, 2.543e-4, 2.363e-4, 2.1977e-4, 
	    2.0452e-4, 1.9042e-4, 1.774e-4, 1.6544e-4, 1.5442e-4, 1.4425e-4, 
	    1.3486e-4, 1.2618e-4, 1.1817e-4, 1.1076e-4, 1.0391e-4, 9.7563e-5, 
	    9.1696e-5, 8.6272e-5, 8.1253e-5, 7.6607e-5, 7.2302e-5, 6.8311e-5, 
	    6.4613e-5, 6.1183e-5, 5.8001e-5, 5.5048e-5, 5.2307e-5, 4.9761e-5, 
	    4.7395e-5, 4.5197e-5, 4.3155e-5, 4.1256e-5, 3.9491e-5, 3.7849e-5, 
	    3.6324e-5, 3.4908e-5, 3.3594e-5, 3.2374e-5, 3.1244e-5, 3.0201e-5, 
	    2.924e-5, 2.8356e-5, 2.7547e-5, 2.6814e-5, 2.6147e-5, 2.5551e-5, 
	    2.5029e-5, 2.4582e-5, 2.4203e-5, 2.3891e-5, 2.3663e-5, 2.3531e-5, 
	    2.3483e-5, 2.3516e-5, 2.3694e-5, 2.4032e-5, 2.4579e-5, 2.5234e-5, 
	    2.6032e-5, 2.7119e-5, 2.8631e-5, 3.0848e-5, 3.3262e-5, 3.6635e-5, 
	    4.0732e-5, 4.5923e-5, 5.3373e-5, 6.1875e-5, 7.2031e-5, 8.598e-5, 
	    9.8642e-5, 1.1469e-4, 1.3327e-4, 1.539e-4, 1.7513e-4, 2.0665e-4, 
	    2.3609e-4, 2.622e-4, 2.8677e-4, 3.259e-4, 3.8624e-4, 4.157e-4, 
	    4.5207e-4, 4.9336e-4, 5.45e-4, 5.8258e-4, 5.8086e-4, 5.6977e-4, 
	    5.3085e-4, 4.802e-4, 4.3915e-4, 4.0343e-4, 3.7853e-4, 3.7025e-4, 
	    3.9637e-4, 4.4675e-4, 4.7072e-4, 4.9022e-4, 5.2076e-4, 5.3676e-4, 
	    5.2755e-4, 4.8244e-4, 4.5473e-4, 4.3952e-4, 3.9614e-4, 3.4086e-4, 
	    2.9733e-4, 2.6367e-4, 2.3767e-4, 2.0427e-4, 1.7595e-4, 1.5493e-4, 
	    1.3851e-4, 1.1874e-4, 1.0735e-4, 9.049e-5, 8.1149e-5, 7.4788e-5, 
	    6.5438e-5, 5.8248e-5, 4.8076e-5, 4.3488e-5, 3.7856e-5, 3.3034e-5, 
	    2.9592e-5, 2.6088e-5, 2.3497e-5, 2.0279e-5, 1.7526e-5, 1.5714e-5, 
	    1.3553e-5, 1.2145e-5, 1.0802e-5, 9.7681e-6, 8.8196e-6, 7.8291e-6, 
	    7.1335e-6, 6.4234e-6, 5.8391e-6, 5.3532e-6, 4.9079e-6, 4.5378e-6, 
	    4.1716e-6, 3.8649e-6, 3.5893e-6, 3.3406e-6, 3.1199e-6, 2.9172e-6, 
	    2.7348e-6, 2.5644e-6, 2.4086e-6, 2.2664e-6, 2.1359e-6, 2.0159e-6, 
	    1.9051e-6, 1.8031e-6, 1.7074e-6, 1.6185e-6, 1.5356e-6, 1.4584e-6, 
	    1.3861e-6, 1.3179e-6, 1.2545e-6, 1.1951e-6, 1.1395e-6, 1.0873e-6, 
	    1.0384e-6, 9.925e-7, 9.4935e-7, 9.0873e-7, 8.705e-7, 8.3446e-7, 
	    8.0046e-7, 7.6834e-7, 7.38e-7, 7.0931e-7, 6.8217e-7, 6.5648e-7, 
	    6.3214e-7, 6.0909e-7, 5.8725e-7, 5.6655e-7, 5.4693e-7, 5.2835e-7, 
	    5.1077e-7, 4.9416e-7, 4.7853e-7, 4.6381e-7, 4.5007e-7, 4.3728e-7, 
	    4.255e-7, 4.145e-7, 4.0459e-7, 3.9532e-7, 3.8662e-7, 3.7855e-7, 
	    3.7041e-7, 3.6254e-7, 3.542e-7, 3.4617e-7, 3.3838e-7, 3.3212e-7, 
	    3.2655e-7, 3.1865e-7, 3.1203e-7, 3.067e-7, 3.0252e-7, 2.9749e-7, 
	    2.9184e-7, 2.8795e-7, 2.8501e-7, 2.8202e-7, 2.7856e-7, 2.7509e-7, 
	    2.7152e-7, 2.6844e-7, 2.6642e-7, 2.6548e-7, 2.6617e-7, 2.6916e-7, 
	    2.7372e-7, 2.8094e-7, 2.9236e-7, 3.1035e-7, 3.2854e-7, 3.5481e-7, 
	    3.9377e-7, 4.4692e-7, 5.0761e-7, 5.7715e-7, 6.7725e-7, 8.0668e-7, 
	    9.3716e-7, 1.0797e-6, 1.1689e-6, 1.3217e-6, 1.4814e-6, 1.5627e-6, 
	    1.6519e-6, 1.7601e-6, 1.906e-6, 2.0474e-6, 2.0716e-6, 2.0433e-6, 
	    1.9752e-6, 1.8466e-6, 1.7526e-6, 1.6657e-6, 1.587e-6, 1.5633e-6, 
	    1.652e-6, 1.8471e-6, 1.9953e-6, 2.0975e-6, 2.2016e-6, 2.2542e-6, 
	    2.3081e-6, 2.3209e-6, 2.2998e-6, 2.3056e-6, 2.2757e-6, 2.2685e-6, 
	    2.2779e-6, 2.2348e-6, 2.2445e-6, 2.3174e-6, 2.4284e-6, 2.529e-6, 
	    2.734e-6, 2.972e-6, 3.2332e-6, 3.5392e-6, 3.9013e-6, 4.3334e-6, 
	    4.9088e-6, 5.3428e-6, 5.9142e-6, 6.6106e-6, 7.4709e-6, 8.5019e-6, 
	    9.6835e-6, 1.0984e-5, 1.2831e-5, 1.4664e-5, 1.708e-5, 2.0103e-5, 
	    2.4148e-5, 2.7948e-5, 3.2855e-5, 3.9046e-5, 4.6429e-5, 5.6633e-5, 
	    6.6305e-5, 7.6048e-5, 8.7398e-5, 1.0034e-4, 1.1169e-4, 1.2813e-4, 
	    1.3354e-4, 1.3952e-4, 1.4204e-4, 1.4615e-4, 1.5144e-4, 1.5475e-4, 
	    1.6561e-4, 1.7135e-4, 1.6831e-4, 1.6429e-4, 1.6353e-4, 1.6543e-4, 
	    1.5944e-4, 1.5404e-4, 1.5458e-4, 1.6287e-4, 1.7277e-4, 1.8387e-4, 
	    1.7622e-4, 1.636e-4, 1.5273e-4, 1.3667e-4, 1.2364e-4, 9.7576e-5, 
	    7.914e-5, 6.4241e-5, 5.1826e-5, 4.1415e-5, 3.1347e-5, 2.5125e-5, 
	    2.0027e-5, 1.6362e-5, 1.3364e-5, 1.1117e-5, 9.4992e-6, 8.1581e-6, 
	    7.1512e-6, 6.2692e-6, 5.5285e-6, 4.9e-6, 4.3447e-6, 3.8906e-6, 
	    3.4679e-6, 3.1089e-6, 2.8115e-6, 2.5496e-6, 2.2982e-6, 2.0861e-6, 
	    1.8763e-6, 1.7035e-6, 1.5548e-6, 1.4107e-6, 1.2839e-6, 1.1706e-6, 
	    1.0709e-6, 9.8099e-7, 8.9901e-7, 8.2394e-7, 7.5567e-7, 6.9434e-7, 
	    6.3867e-7, 5.8845e-7, 5.4263e-7, 5.0033e-7, 4.6181e-7, 4.2652e-7, 
	    3.9437e-7, 3.6497e-7, 3.3781e-7, 3.1292e-7, 2.9011e-7, 2.6915e-7, 
	    2.4989e-7, 2.3215e-7, 2.1582e-7, 2.0081e-7, 1.87e-7, 1.7432e-7, 
	    1.6264e-7, 1.5191e-7, 1.4207e-7, 1.3306e-7, 1.2484e-7, 1.1737e-7, 
	    1.1056e-7, 1.0451e-7, 9.906e-8, 9.4135e-8, 8.9608e-8, 8.5697e-8, 
	    8.1945e-8, 7.8308e-8, 7.4808e-8, 7.1686e-8, 6.8923e-8, 6.5869e-8, 
	    6.3308e-8, 6.084e-8, 5.8676e-8, 5.6744e-8, 5.5016e-8, 5.3813e-8, 
	    5.2792e-8, 5.2097e-8, 5.1737e-8, 5.1603e-8, 5.1656e-8, 5.1989e-8, 
	    5.2467e-8, 5.2918e-8, 5.3589e-8, 5.456e-8, 5.5869e-8, 5.7403e-8, 
	    5.8968e-8, 6.0973e-8, 6.3432e-8, 6.6245e-8, 6.9353e-8, 7.2686e-8, 
	    7.6541e-8, 8.0991e-8, 8.595e-8, 9.1429e-8, 9.7851e-8, 1.0516e-7, 
	    1.1349e-7, 1.2295e-7, 1.3335e-7, 1.4488e-7, 1.5864e-7, 1.7412e-7, 
	    1.914e-7, 2.1078e-7, 2.3369e-7, 2.5996e-7, 2.8848e-7, 3.2169e-7, 
	    3.5991e-7, 4.0566e-7, 4.5969e-7, 5.3094e-7, 6.1458e-7, 7.1155e-7, 
	    8.3045e-7, 9.9021e-7, 1.2042e-6, 1.4914e-6, 1.8145e-6, 2.221e-6, 
	    2.7831e-6, 3.4533e-6, 4.4446e-6, 5.1989e-6, 6.2289e-6, 7.1167e-6, 
	    8.3949e-6, 9.6417e-6, 1.0313e-5, 1.0485e-5, 1.0641e-5, 1.0898e-5, 
	    1.0763e-5, 1.0506e-5, 1.0497e-5, 1.1696e-5, 1.2654e-5, 1.3029e-5, 
	    1.3175e-5, 1.4264e-5, 1.4985e-5, 1.4999e-5, 1.4317e-5, 1.4616e-5, 
	    1.4963e-5, 1.5208e-5, 1.4942e-5, 1.3879e-5, 1.3087e-5, 1.1727e-5, 
	    1.0515e-5, 9.0073e-6, 7.3133e-6, 6.1181e-6, 5.0623e-6, 4.1105e-6, 
	    3.3915e-6, 2.6711e-6, 2.1464e-6, 1.7335e-6, 1.4302e-6, 1.1847e-6, 
	    9.9434e-7, 8.2689e-7, 7.0589e-7, 6.075e-7, 5.3176e-7, 4.6936e-7, 
	    4.1541e-7, 3.6625e-7, 3.2509e-7, 2.9156e-7, 2.6308e-7, 2.3819e-7, 
	    2.1421e-7, 1.9366e-7, 1.7626e-7, 1.5982e-7, 1.4567e-7, 1.3354e-7, 
	    1.2097e-7, 1.1029e-7, 1.0063e-7, 9.2003e-8, 8.4245e-8, 7.7004e-8, 
	    7.0636e-8, 6.4923e-8, 5.9503e-8, 5.4742e-8, 5.045e-8, 4.647e-8, 
	    4.2881e-8, 3.955e-8, 3.6541e-8, 3.3803e-8, 3.1279e-8, 2.8955e-8, 
	    2.6858e-8, 2.4905e-8, 2.3146e-8, 2.1539e-8, 2.0079e-8, 1.8746e-8, 
	    1.7517e-8, 1.6396e-8, 1.5369e-8, 1.4426e-8, 1.3543e-8, 1.2724e-8, 
	    1.1965e-8, 1.1267e-8, 1.0617e-8, 1.001e-8, 9.4662e-9, 8.9553e-9, 
	    8.4988e-9, 8.0807e-9, 7.7043e-9, 7.3721e-9, 7.0707e-9, 6.8047e-9, 
	    6.5702e-9, 6.3634e-9, 6.1817e-9, 6.0239e-9, 5.8922e-9, 5.7824e-9, 
	    5.7019e-9, 5.6368e-9, 5.594e-9, 5.5669e-9, 5.5583e-9, 5.5653e-9, 
	    5.5837e-9, 5.6243e-9, 5.6883e-9, 5.78e-9, 5.8964e-9, 6.0429e-9, 
	    6.2211e-9, 6.4282e-9, 6.6634e-9, 6.9306e-9, 7.2336e-9, 7.5739e-9, 
	    7.9562e-9, 8.3779e-9, 8.8575e-9, 9.3992e-9, 1.0004e-8, 1.0684e-8, 
	    1.145e-8, 1.232e-8, 1.3311e-8, 1.4455e-8, 1.5758e-8, 1.7254e-8, 
	    1.8927e-8, 2.093e-8, 2.3348e-8, 2.6074e-8, 2.9221e-8, 3.277e-8, 
	    3.7485e-8, 4.2569e-8, 4.8981e-8, 5.5606e-8, 6.2393e-8, 7.1901e-8, 
	    8.2921e-8, 9.5513e-8, 1.1111e-7, 1.3143e-7, 1.5971e-7, 1.8927e-7, 
	    2.2643e-7, 2.786e-7, 3.2591e-7, 3.7024e-7, 4.2059e-7, 4.9432e-7, 
	    5.5543e-7, 5.7498e-7, 5.921e-7, 6.1005e-7, 6.1577e-7, 5.9193e-7, 
	    5.6602e-7, 5.7403e-7, 6.005e-7, 6.4723e-7, 6.7073e-7, 7.5415e-7, 
	    8.0982e-7, 8.7658e-7, 9.143e-7, 9.4459e-7, 9.8347e-7, 9.8768e-7, 
	    1.0153e-6, 1.0066e-6, 1.0353e-6, 1.0353e-6, 1.0722e-6, 1.1138e-6, 
	    1.1923e-6, 1.2947e-6, 1.4431e-6, 1.6537e-6, 1.8662e-6, 2.2473e-6, 
	    2.6464e-6, 3.1041e-6, 3.4858e-6, 4.0167e-6, 4.6675e-6, 5.0983e-6, 
	    5.7997e-6, 6.0503e-6, 6.4687e-6, 6.5396e-6, 6.7986e-6, 7.0244e-6, 
	    7.2305e-6, 7.6732e-6, 7.9783e-6, 7.9846e-6, 7.7617e-6, 7.7657e-6, 
	    7.7411e-6, 7.8816e-6, 7.8136e-6, 8.0051e-6, 8.5799e-6, 9.1659e-6, 
	    9.8646e-6, 9.492e-6, 8.767e-6, 8.2034e-6, 7.2297e-6, 6.2324e-6, 
	    4.9315e-6, 3.9128e-6, 3.1517e-6, 2.4469e-6, 1.8815e-6, 1.4627e-6, 
	    1.1698e-6, 9.4686e-7, 7.8486e-7, 6.697e-7, 5.8811e-7, 5.2198e-7, 
	    4.6809e-7, 4.1671e-7, 3.7006e-7, 3.3066e-7, 2.9387e-7, 2.6415e-7, 
	    2.3409e-7, 2.0991e-7, 1.9132e-7, 1.7519e-7, 1.5939e-7, 1.4368e-7, 
	    1.305e-7, 1.1883e-7, 1.0772e-7, 9.6884e-8, 8.7888e-8, 7.8956e-8, 
	    7.1024e-8, 6.3824e-8, 5.7256e-8, 5.1769e-8, 4.7037e-8, 4.2901e-8, 
	    3.897e-8, 3.5467e-8, 3.2502e-8, 2.9827e-8, 2.7389e-8, 2.5111e-8, 
	    2.3056e-8, 2.1267e-8, 1.961e-8, 1.8133e-8, 1.6775e-8, 1.5491e-8, 
	    1.4329e-8, 1.3265e-8, 1.23e-8, 1.142e-8, 1.0593e-8, 9.8475e-9, 
	    9.1585e-9, 8.5256e-9, 7.9525e-9, 7.4226e-9, 6.9379e-9, 6.495e-9, 
	    6.0911e-9, 5.7242e-9, 5.3877e-9, 5.0821e-9, 4.8051e-9, 4.5554e-9, 
	    4.3315e-9, 4.1336e-9, 3.9632e-9, 3.8185e-9, 3.708e-9, 3.6296e-9, 
	    3.5804e-9, 3.5776e-9, 3.6253e-9, 3.7115e-9, 3.8151e-9, 3.9804e-9, 
	    4.1742e-9, 4.3581e-9, 4.5306e-9, 4.7736e-9, 5.1297e-9, 5.5291e-9, 
	    5.9125e-9, 6.4956e-9, 7.0362e-9, 7.5318e-9, 7.9947e-9, 8.6438e-9, 
	    9.7227e-9, 1.013e-8, 1.0549e-8, 1.1064e-8, 1.1702e-8, 1.2043e-8, 
	    1.1781e-8, 1.1838e-8, 1.1917e-8, 1.2131e-8, 1.2476e-8, 1.3611e-8, 
	    1.436e-8, 1.5057e-8, 1.6247e-8, 1.7284e-8, 1.842e-8, 1.8352e-8, 
	    1.8722e-8, 1.9112e-8, 1.9092e-8, 1.9311e-8, 1.9411e-8, 1.9884e-8, 
	    2.0508e-8, 2.151e-8, 2.3143e-8, 2.505e-8, 2.7596e-8, 3.1231e-8, 
	    3.626e-8, 4.341e-8, 5.224e-8, 6.3236e-8, 7.7522e-8, 9.8688e-8, 
	    1.1859e-7, 1.4341e-7, 1.6798e-7, 1.9825e-7, 2.2898e-7, 2.6257e-7, 
	    2.9884e-7, 3.3247e-7, 3.4936e-7, 3.5583e-7, 3.715e-7, 3.658e-7, 
	    3.7124e-7, 3.703e-7, 4.1536e-7, 4.6656e-7, 4.6677e-7, 4.7507e-7, 
	    4.9653e-7, 5.3795e-7, 5.4957e-7, 5.2238e-7, 5.469e-7, 5.6569e-7, 
	    5.9844e-7, 5.9835e-7, 5.6522e-7, 5.4123e-7, 4.7904e-7, 4.2851e-7, 
	    3.5603e-7, 2.8932e-7, 2.3655e-7, 1.8592e-7, 1.4943e-7, 1.1971e-7, 
	    9.8482e-8, 8.3675e-8, 7.127e-8, 6.2496e-8, 5.4999e-8, 4.9821e-8, 
	    4.5387e-8, 4.134e-8, 3.7453e-8, 3.3298e-8, 3.012e-8, 2.7032e-8, 
	    2.4236e-8, 2.15e-8, 1.8988e-8, 1.7414e-8, 1.5706e-8, 1.4192e-8, 
	    1.3204e-8, 1.1759e-8, 1.0737e-8, 9.6309e-9, 8.8179e-9, 8.2619e-9, 
	    7.2264e-9, 6.4856e-9, 5.8037e-9, 5.2093e-9, 4.7205e-9, 4.1749e-9, 
	    3.7852e-9, 3.3915e-9, 3.0089e-9, 2.7335e-9, 2.4398e-9, 2.2031e-9, 
	    1.9786e-9, 1.789e-9, 1.6266e-9, 1.483e-9, 1.3576e-9, 1.2518e-9, 
	    1.1587e-9, 1.0726e-9, 9.9106e-10, 9.1673e-10, 8.5084e-10, 
	    7.9147e-10, 7.2882e-10, 6.7342e-10, 6.2593e-10, 5.8294e-10, 
	    5.4435e-10, 5.0997e-10, 4.7806e-10, 4.4931e-10, 4.2357e-10, 
	    4.0023e-10, 3.7909e-10, 3.5999e-10, 3.4285e-10, 3.2776e-10, 
	    3.1468e-10, 3.0377e-10, 2.9479e-10, 2.8877e-10, 2.8512e-10, 
	    2.8617e-10, 2.8976e-10, 3.0001e-10, 3.1718e-10, 3.3898e-10, 
	    3.5857e-10, 3.8358e-10, 4.3131e-10, 4.5741e-10, 4.6948e-10, 
	    4.7594e-10, 4.9529e-10, 5.1563e-10, 4.9475e-10, 4.8369e-10, 
	    4.8829e-10, 5.0047e-10, 5.0203e-10, 5.1954e-10, 5.5352e-10, 
	    5.9928e-10, 6.7148e-10, 7.1121e-10, 7.4317e-10, 7.6039e-10, 
	    7.8313e-10, 8.0684e-10, 7.8553e-10, 7.8312e-10, 7.8537e-10, 
	    7.8872e-10, 8.0185e-10, 8.1004e-10, 8.2608e-10, 8.2525e-10, 
	    8.3857e-10, 8.792e-10, 9.2451e-10, 9.8661e-10, 1.0629e-9, 
	    1.1659e-9, 1.2922e-9, 1.4387e-9, 1.6254e-9, 1.8425e-9, 2.1428e-9, 
	    2.5477e-9, 3.0379e-9, 3.757e-9, 4.4354e-9, 5.1802e-9, 6.2769e-9, 
	    7.4894e-9, 8.7474e-9, 9.8037e-9, 1.1582e-8, 1.3293e-8, 1.4471e-8, 
	    1.5025e-8, 1.558e-8, 1.6228e-8, 1.6413e-8, 1.602e-8, 1.6393e-8, 
	    1.7545e-8, 1.959e-8, 2.1449e-8, 2.3856e-8, 2.705e-8, 3.0214e-8, 
	    3.3733e-8, 3.6487e-8, 3.9353e-8, 4.266e-8, 4.6385e-8, 4.9955e-8, 
	    5.5313e-8, 6.0923e-8, 6.8948e-8, 7.3649e-8, 8.2602e-8, 9.2212e-8, 
	    9.908e-8, 1.1319e-7, 1.179e-7, 1.2941e-7, 1.3199e-7, 1.3914e-7, 
	    1.4843e-7, 1.53e-7, 1.6419e-7, 1.7095e-7, 1.6988e-7, 1.6494e-7, 
	    1.6327e-7, 1.6067e-7, 1.6909e-7, 1.7118e-7, 1.8106e-7, 1.9857e-7, 
	    2.1696e-7, 2.3385e-7, 2.2776e-7, 2.1402e-7, 1.9882e-7, 1.7362e-7, 
	    1.4308e-7, 1.1158e-7, 8.8781e-8, 6.8689e-8, 5.2062e-8, 4.0427e-8, 
	    3.2669e-8, 2.7354e-8, 2.32e-8, 2.058e-8, 1.8676e-8, 1.7329e-8, 
	    1.6621e-8, 1.6433e-8, 1.6953e-8, 1.7134e-8, 1.7948e-8, 1.9107e-8, 
	    1.9875e-8, 2.1416e-8, 2.1556e-8, 2.2265e-8, 2.2171e-8, 2.2534e-8, 
	    2.3029e-8, 2.2828e-8, 2.3143e-8, 2.2965e-8, 2.2223e-8, 2.1108e-8, 
	    2.0265e-8, 1.9516e-8, 1.9941e-8, 2.0312e-8, 2.108e-8, 2.2611e-8, 
	    2.421e-8, 2.6069e-8, 2.5097e-8, 2.3318e-8, 2.1543e-8, 1.8942e-8, 
	    1.596e-8, 1.2386e-8, 9.934e-9, 7.7502e-9, 5.9462e-9, 4.5113e-9, 
	    3.5523e-9, 2.8844e-9, 2.3394e-9, 1.9584e-9, 1.6749e-9, 1.4624e-9, 
	    1.2809e-9, 1.1359e-9, 1.0087e-9, 9.0166e-10, 8.1079e-10, 
	    7.2219e-10, 6.4922e-10, 5.8803e-10, 5.329e-10, 4.859e-10, 
	    4.4111e-10, 4.0184e-10, 3.6644e-10, 3.3529e-10, 3.0789e-10, 
	    2.8286e-10, 2.6089e-10, 2.4125e-10, 2.2355e-10, 2.0783e-10, 
	    1.937e-10, 1.8088e-10, 1.6948e-10, 1.5929e-10, 1.5013e-10, 
	    1.4193e-10, 1.347e-10, 1.2841e-10, 1.2307e-10, 1.1865e-10, 
	    1.1502e-10, 1.1243e-10, 1.1099e-10, 1.1066e-10, 1.1216e-10, 
	    1.1529e-10, 1.2171e-10, 1.3128e-10, 1.4153e-10, 1.5962e-10, 
	    1.8048e-10, 2.0936e-10, 2.3165e-10, 2.5746e-10, 2.96e-10, 
	    3.3707e-10, 3.5267e-10, 3.5953e-10, 3.6822e-10, 3.8363e-10, 
	    3.8286e-10, 3.5883e-10, 3.6154e-10, 3.6653e-10, 3.8507e-10, 
	    4.025e-10, 4.4435e-10, 4.9889e-10, 5.6932e-10, 6.3599e-10, 
	    7.0281e-10, 7.5777e-10, 8.1279e-10, 8.891e-10, 9.34e-10, 
	    1.0076e-9, 1.0945e-9, 1.1898e-9, 1.3108e-9, 1.4725e-9, 1.7028e-9, 
	    1.9619e-9, 2.3527e-9, 2.6488e-9, 3.0327e-9, 3.4396e-9, 3.8797e-9, 
	    4.4115e-9, 4.6853e-9, 4.9553e-9, 4.9551e-9, 5.1062e-9, 5.0996e-9, 
	    5.1119e-9, 5.2283e-9, 5.8297e-9, 6.3439e-9, 6.2675e-9, 6.3296e-9, 
	    6.5173e-9, 7.1685e-9, 7.0528e-9, 6.8856e-9, 7.3182e-9, 7.699e-9, 
	    8.3461e-9, 8.1946e-9, 7.7153e-9, 7.2411e-9, 6.4511e-9, 5.7336e-9, 
	    4.6105e-9, 3.6962e-9, 2.9944e-9, 2.4317e-9, 1.9399e-9, 1.5331e-9, 
	    1.2633e-9, 1.0613e-9, 9.0136e-10, 7.9313e-10, 7.1543e-10, 
	    6.6485e-10, 6.4225e-10, 6.398e-10, 6.4598e-10, 6.7428e-10, 
	    7.027e-10, 7.4694e-10, 7.7946e-10, 7.9395e-10, 7.8716e-10, 
	    7.6933e-10, 7.622e-10, 7.4825e-10, 7.4805e-10, 7.6511e-10, 
	    7.6492e-10, 7.4103e-10, 7.1979e-10, 7.1686e-10, 7.3403e-10, 
	    7.1142e-10, 7.0212e-10, 7.1548e-10, 7.5253e-10, 8.0444e-10, 
	    8.2378e-10, 7.8004e-10, 7.1712e-10, 6.4978e-10, 5.7573e-10, 
	    4.8675e-10, 3.7945e-10, 3.0118e-10, 2.4241e-10, 1.91e-10, 
	    1.4816e-10, 1.1567e-10, 9.4183e-11, 7.766e-11, 6.527e-11, 
	    5.6616e-11, 4.9576e-11, 4.4137e-11, 3.9459e-11, 3.5759e-11, 
	    3.2478e-11, 2.9419e-11, 2.6703e-11, 2.4365e-11, 2.2412e-11, 
	    2.0606e-11, 1.9067e-11, 1.78e-11, 1.6695e-11, 1.5729e-11, 
	    1.4887e-11, 1.4135e-11, 1.3519e-11, 1.2992e-11, 1.2563e-11, 
	    1.2223e-11, 1.1962e-11, 1.1775e-11, 1.1657e-11, 1.1605e-11, 
	    1.1619e-11, 1.1697e-11, 1.1839e-11, 1.2046e-11, 1.2319e-11, 
	    1.2659e-11, 1.307e-11, 1.3553e-11, 1.4113e-11, 1.4754e-11, 
	    1.548e-11, 1.6298e-11, 1.7214e-11, 1.8236e-11, 1.9372e-11, 
	    2.0635e-11, 2.2036e-11, 2.359e-11, 2.5317e-11, 2.7242e-11, 
	    2.94e-11, 3.1849e-11, 3.4654e-11, 3.7923e-11, 4.1695e-11, 
	    4.6055e-11, 5.094e-11, 5.5624e-11, 6.0667e-11, 6.6261e-11, 
	    7.2692e-11, 7.9711e-11, 8.7976e-11, 9.6884e-11, 1.0775e-10, 
	    1.2093e-10, 1.3531e-10, 1.5404e-10, 1.7315e-10, 1.9862e-10, 
	    2.3341e-10, 2.7014e-10, 3.1716e-10, 3.6957e-10, 4.3233e-10, 
	    5.2566e-10, 6.2251e-10, 7.2149e-10, 8.3958e-10, 9.5931e-10, 
	    1.1388e-9, 1.2973e-9, 1.4442e-9, 1.5638e-9, 1.6974e-9, 1.8489e-9, 
	    1.983e-9, 2.172e-9, 2.3662e-9, 2.6987e-9, 3.1697e-9, 3.6907e-9, 
	    4.2625e-9, 4.7946e-9, 5.3848e-9, 6.0897e-9, 6.473e-9, 7.1483e-9, 
	    7.7432e-9, 8.0851e-9, 8.5013e-9, 8.5909e-9, 9.189e-9, 9.3124e-9, 
	    9.5936e-9, 9.8787e-9, 9.9036e-9, 9.6712e-9, 9.2036e-9, 9.0466e-9, 
	    8.938e-9, 9.1815e-9, 9.5092e-9, 1.0027e-8, 1.0876e-8, 1.1744e-8, 
	    1.1853e-8, 1.1296e-8, 1.0134e-8, 8.8245e-9, 7.393e-9, 5.715e-9, 
	    4.4884e-9, 3.4027e-9, 2.6054e-9, 2.079e-9, 1.7267e-9, 1.4724e-9, 
	    1.2722e-9, 1.1234e-9, 1.0186e-9, 9.468e-10, 8.8854e-10, 
	    8.5127e-10, 8.3157e-10, 8.2226e-10, 8.3395e-10, 8.3294e-10, 
	    8.4725e-10, 8.8814e-10, 9.3697e-10, 1.0112e-9, 1.0412e-9, 
	    1.0948e-9, 1.181e-9, 1.2267e-9, 1.369e-9, 1.4512e-9, 1.5568e-9, 
	    1.6552e-9, 1.7321e-9, 1.8797e-9, 1.921e-9, 1.9686e-9, 1.9917e-9, 
	    1.9357e-9, 1.8486e-9, 1.7575e-9, 1.7113e-9, 1.7163e-9, 1.7623e-9, 
	    1.8536e-9, 1.9765e-9, 2.1334e-9, 2.3237e-9, 2.3259e-9, 2.1833e-9, 
	    1.9785e-9, 1.7308e-9, 1.4596e-9, 1.1198e-9, 8.7375e-10, 
	    6.5381e-10, 4.8677e-10, 3.6756e-10, 2.9155e-10, 2.3735e-10, 
	    1.959e-10, 1.6638e-10, 1.4549e-10, 1.2947e-10, 1.1511e-10, 
	    1.0548e-10, 9.6511e-11, 9.0469e-11, 8.517e-11, 7.7804e-11, 
	    7.1971e-11, 6.6213e-11, 6.1063e-11, 5.5881e-11, 5.0508e-11, 
	    4.5932e-11, 4.1997e-11, 3.7672e-11, 3.3972e-11, 3.0318e-11, 
	    2.6769e-11, 2.3874e-11, 2.1336e-11, 1.9073e-11, 1.7313e-11, 
	    1.5904e-11, 1.4684e-11, 1.3698e-11, 1.2873e-11, 1.2175e-11, 
	    1.1542e-11, 1.1024e-11, 1.0602e-11, 1.0267e-11, 1.0012e-11, 
	    9.8379e-12, 9.7482e-12, 9.7564e-12, 9.8613e-12, 1.0092e-11, 
	    1.0418e-11, 1.0868e-11, 1.1585e-11, 1.2351e-11, 1.3372e-11, 
	    1.4841e-11, 1.6457e-11, 1.8681e-11, 2.055e-11, 2.2912e-11, 
	    2.5958e-11, 2.9137e-11, 3.2368e-11, 3.4848e-11, 3.8462e-11, 
	    4.219e-11, 4.5629e-11, 4.9022e-11, 5.4232e-11, 6.19e-11, 
	    7.1953e-11, 8.5368e-11, 9.9699e-11, 1.1734e-10, 1.4185e-10, 
	    1.7017e-10, 1.9813e-10, 2.3859e-10, 2.7304e-10, 3.0971e-10, 
	    3.5129e-10, 3.9405e-10, 4.5194e-10, 4.8932e-10, 5.2436e-10, 
	    5.4098e-10, 5.5542e-10, 5.7794e-10, 5.6992e-10, 5.879e-10, 
	    6.1526e-10, 6.8034e-10, 6.7956e-10, 6.6864e-10, 6.9329e-10, 
	    7.2971e-10, 7.6546e-10, 7.5078e-10, 7.8406e-10, 8.3896e-10, 
	    9.0111e-10, 9.1994e-10, 8.7189e-10, 8.1426e-10, 7.3097e-10, 
	    6.3357e-10, 5.1371e-10, 4.0936e-10, 3.2918e-10, 2.6255e-10, 
	    2.0724e-10, 1.6879e-10, 1.4165e-10, 1.1989e-10, 1.0125e-10, 
	    8.9629e-11, 7.8458e-11, 6.8826e-11, 6.0935e-11, 5.5208e-11, 
	    5.2262e-11, 5.026e-11, 4.8457e-11, 4.7888e-11, 4.8032e-11, 
	    5.0838e-11, 5.4668e-11, 5.579e-11, 6.0056e-11, 6.3811e-11, 
	    6.8848e-11, 7.459e-11, 7.8249e-11, 8.3371e-11, 8.3641e-11, 
	    8.6591e-11, 8.9599e-11, 9.3487e-11, 1.0066e-10, 1.0765e-10, 
	    1.0851e-10, 1.0619e-10, 1.0557e-10, 1.046e-10, 1.0796e-10, 
	    1.0523e-10, 1.0674e-10, 1.1261e-10, 1.1431e-10, 1.1408e-10, 
	    1.0901e-10, 9.9105e-11, 8.8077e-11, 6.9928e-11, 5.4595e-11, 
	    4.5401e-11, 3.6313e-11, 2.6986e-11, 1.9463e-11, 1.4577e-11, 
	    1.1583e-11, 9.5492e-12, 8.077e-12, 6.9642e-12, 6.0966e-12, 
	    5.4046e-12, 4.8431e-12, 4.3815e-12, 3.9987e-12, 3.679e-12, 
	    3.4113e-12, 3.1868e-12, 2.9992e-12, 2.8434e-12, 2.7153e-12, 
	    2.612e-12, 2.5311e-12, 2.4705e-12, 2.429e-12, 2.4053e-12, 
	    2.3988e-12, 2.4087e-12, 2.4349e-12, 2.4771e-12, 2.5355e-12, 
	    2.6103e-12, 2.7019e-12, 2.811e-12, 2.9383e-12, 3.0848e-12, 
	    3.2518e-12, 3.4405e-12, 3.6527e-12, 3.8902e-12, 4.1555e-12, 
	    4.451e-12, 4.7801e-12, 5.1462e-12, 5.5539e-12, 6.0086e-12, 
	    6.5171e-12, 7.0884e-12, 7.7357e-12, 8.4831e-12, 9.3096e-12, 
	    1.0282e-11, 1.1407e-11, 1.269e-11, 1.4148e-11, 1.5888e-11, 
	    1.7992e-11, 2.0523e-11, 2.3342e-11, 2.6578e-11, 3.0909e-11, 
	    3.6228e-11, 4.2053e-11, 4.9059e-11, 5.9273e-11, 7.0166e-11, 
	    8.2298e-11, 9.7071e-11, 1.1673e-10, 1.401e-10, 1.6621e-10, 
	    2.0127e-10, 2.3586e-10, 2.705e-10, 3.095e-10, 3.6584e-10, 
	    4.1278e-10, 4.6591e-10, 5.222e-10, 5.5246e-10, 6.15e-10, 
	    6.5878e-10, 7.1167e-10, 7.9372e-10, 8.6975e-10, 9.6459e-10, 
	    9.7368e-10, 9.8142e-10, 1.0202e-9, 1.02e-9, 1.0356e-9, 1.0092e-9, 
	    1.0269e-9, 1.0366e-9, 1.049e-9, 1.0717e-9, 1.0792e-9, 1.1016e-9, 
	    1.0849e-9, 1.0929e-9, 1.0971e-9, 1.0969e-9, 1.046e-9, 9.2026e-10, 
	    8.1113e-10, 6.8635e-10, 5.5369e-10, 4.2908e-10, 3.3384e-10, 
	    2.648e-10, 2.081e-10, 1.6915e-10, 1.4051e-10, 1.1867e-10, 
	    1.0158e-10, 8.899e-11, 7.9175e-11, 7.044e-11, 6.3453e-11, 
	    5.7009e-11, 5.1662e-11, 4.7219e-11, 4.3454e-11, 4.0229e-11, 
	    3.7689e-11, 3.6567e-11, 3.5865e-11, 3.5955e-11, 3.5928e-11, 
	    3.6298e-11, 3.7629e-11, 3.93e-11, 4.1829e-11, 4.4806e-11, 
	    5.0534e-11, 5.6672e-11, 6.2138e-11, 6.8678e-11, 7.6111e-11, 
	    8.4591e-11, 9.2634e-11, 9.8085e-11, 1.083e-10, 1.1949e-10, 
	    1.2511e-10, 1.3394e-10, 1.3505e-10, 1.4342e-10, 1.4874e-10, 
	    1.492e-10, 1.5872e-10, 1.5972e-10, 1.5821e-10, 1.5425e-10, 
	    1.4937e-10, 1.5089e-10, 1.5521e-10, 1.6325e-10, 1.6924e-10, 
	    1.8265e-10, 1.9612e-10, 2.0176e-10, 1.9359e-10, 1.7085e-10, 
	    1.5197e-10, 1.2646e-10, 9.8552e-11, 7.453e-11, 5.5052e-11, 
	    4.2315e-11, 3.2736e-11, 2.6171e-11, 2.1909e-11, 1.8286e-11, 
	    1.5752e-11, 1.3859e-11, 1.2288e-11, 1.1002e-11, 9.7534e-12, 
	    8.8412e-12, 8.0169e-12, 7.2855e-12, 6.8734e-12, 6.4121e-12, 
	    6.1471e-12, 5.778e-12, 5.3478e-12, 4.9652e-12, 4.4043e-12, 
	    3.9862e-12, 3.4684e-12, 2.9681e-12, 2.5791e-12, 2.2339e-12, 
	    1.9247e-12, 1.6849e-12, 1.4863e-12, 1.3291e-12, 1.2021e-12, 
	    1.0947e-12, 1.0015e-12, 9.1935e-13, 8.4612e-13, 7.8036e-13, 
	    7.21e-13, 6.6718e-13, 6.1821e-13, 5.7353e-13, 5.3269e-13, 
	    4.9526e-13, 4.6093e-13, 4.2937e-13, 4.0034e-13, 3.7361e-13, 
	    3.4895e-13, 3.2621e-13, 3.052e-13, 2.8578e-13, 2.6782e-13, 
	    2.512e-13, 2.3581e-13, 2.2154e-13, 2.0832e-13, 1.9605e-13, 
	    1.8466e-13, 1.7408e-13, 1.6425e-13, 1.5511e-13, 1.4661e-13, 
	    1.3869e-13, 1.3131e-13, 1.2444e-13, 1.1803e-13, 1.1205e-13, 
	    1.0646e-13, 1.0124e-13, 9.6358e-14, 9.1789e-14, 8.7509e-14, 
	    8.3498e-14, 7.9735e-14, 7.6202e-14, 7.2882e-14, 6.976e-14, 
	    6.6822e-14, 6.4053e-14, 6.1442e-14, 5.8978e-14, 5.665e-14, 
	    5.4448e-14, 5.2364e-14, 5.0389e-14, 4.8516e-14, 4.6738e-14, 
	    4.5048e-14, 4.3441e-14, 4.1911e-14, 4.0453e-14, 3.9063e-14, 
	    3.7735e-14, 3.6467e-14, 3.5254e-14, 3.4093e-14, 3.298e-14, 
	    3.1914e-14, 3.0891e-14, 2.9909e-14, 2.8965e-14, 2.8058e-14, 
	    2.7185e-14, 2.6344e-14, 2.5535e-14, 2.4755e-14, 2.4002e-14, 
	    2.3276e-14, 2.2576e-14, 2.1899e-14, 2.1245e-14, 2.0613e-14, 
	    2.0002e-14, 1.9411e-14, 1.8839e-14, 1.8285e-14, 1.7749e-14, 
	    1.723e-14, 1.6727e-14, 1.624e-14, 1.5768e-14, 1.531e-14, 
	    1.4867e-14, 1.4436e-14, 1.4019e-14, 1.3614e-14, 1.3221e-14, 
	    1.284e-14, 1.2471e-14, 1.2112e-14, 1.1764e-14, 1.1425e-14, 
	    1.1097e-14, 1.0779e-14, 1.0469e-14, 1.0169e-14, 9.8775e-15, 
	    9.5943e-15, 9.3193e-15, 9.0522e-15, 8.7928e-15, 8.5409e-15, 
	    8.2962e-15, 8.0586e-15, 7.8278e-15, 7.6036e-15, 7.3858e-15, 
	    7.1742e-15, 6.9687e-15, 6.7691e-15, 6.5752e-15, 6.3868e-15, 
	    6.2038e-15, 6.026e-15, 5.8533e-15, 5.6856e-15, 5.5226e-15, 
	    5.3642e-15, 5.2104e-15, 5.061e-15, 4.9158e-15, 4.7748e-15, 
	    4.6378e-15, 4.5047e-15, 4.3753e-15, 4.2497e-15, 4.1277e-15, 
	    4.0091e-15, 3.8939e-15, 3.782e-15, 3.6733e-15, 3.5677e-15, 
	    3.4651e-15, 3.3655e-15, 3.2686e-15, 3.1746e-15, 3.0832e-15, 
	    2.9944e-15, 2.9082e-15, 2.8244e-15, 2.7431e-15, 2.664e-15, 
	    2.5872e-15, 2.5126e-15, 2.4401e-15, 2.3697e-15, 2.3014e-15, 
	    2.2349e-15, 2.1704e-15, 2.1077e-15, 2.0468e-15, 1.9877e-15, 
	    1.9302e-15, 1.8744e-15, 1.8202e-15, 1.7675e-15, 1.7164e-15, 
	    1.6667e-15, 1.6184e-15, 1.5716e-15, 1.526e-15, 1.4818e-15, 
	    1.4389e-15, 1.3971e-15, 1.3566e-15, 1.3172e-15, 1.279e-15, 
	    1.2419e-15, 1.2058e-15, 1.1708e-15, 1.1368e-15, 1.1037e-15, 
	    1.0716e-15, 1.0405e-15, 1.0102e-15, 9.8079e-16, 9.5224e-16, 
	    9.2451e-16, 8.9758e-16, 8.7142e-16, 8.4602e-16, 8.2136e-16, 
	    7.974e-16, 7.7414e-16, 7.5154e-16, 7.2961e-16, 7.083e-16, 
	    6.8761e-16, 6.6752e-16, 6.4801e-16, 6.2906e-16, 6.1066e-16, 
	    5.928e-16, 5.7545e-16, 5.586e-16, 5.4224e-16, 5.2636e-16, 
	    5.1094e-16, 4.9596e-16 };

struct {
    double e_1[3];
    int e_2;
    } s260b_ = { -20., 2e4, 10., 2003 };

struct {
    double e_1[9];
    } consts_ = { 3.1415927410125732, 6.62606876e-27, 1.3806503e-16, 
	    29979245800., 6.02214199e23, 2.6867775e19, 83144720., 
	    1.191042722e-12, 1.4387752 };


/* Table of constant values */

/*
static integer c__9 = 9;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__5 = 5;
static int cs__0 = 0;
*/
static double c_b125 = 0.;

/* ############################################################################ */
/*     path:		$Source: /srv/svn/cvs/cvsroot/arts/src/continua.cc,v $ */
/*     author:		$Author $ */
/*     revision:	        $Revision: 1.26.2.17 $ */
/*     created:	        $Date: 2003/11/21 11:28:04 $ */
/* ############################################################################ */

/* CKD2.4 TEST */
/* TKS, 2002-02-28 */
/* CALL : g77 -c testckd.f ; g77 testckd.o -o testckd */

/* ############################################################################ */
/* ----------------------------------------------------------------------------- */

/* INPUT PARAMETERS: */
/*  P          [hPa]  TOTAL PRESSURE */
/*  T          [K]    TEMPERATURE */
/*  VMRH2O     [1]    H2O VOLUME MIXING RATIO */
/*  VMRN2      [1]    N2  VOLUME MIXING RATIO */
/*  VMRO2      [1]    O2  VOLUME MIXING RATIO */
/*  FREQ       [Hz]   FREQUENCY OF ABSORPTION CALCULATION */


/* OUTPUT PARAMETER: */
/*  artsckd_     [1/m]  ABSORPTION COEFFICIENT */

/* ----------------------------------------------------------------------------- */

double artsckd_(double p, double t, double vmrh2o, 
	double vmrn2, double vmro2, double freq, int ivc)
{
    /* Initialized data */

    static double xslf = 1.;
    static double xfrg = 1.;
    static double xcn2 = 1.;

    /* System generated locals */
    double ret_val=0.0e0;
    double d__1, d__2, d__3, d__4;

    /* Local variables */
    int iosa;
    double w_wv__, oc_n2, radct;
    double w_other__, w_n2__, w_o2__;
    double of_wv, os_wv, p0, xn_wv__, t0, rhofac, wn, xn, xn0, tksvpt, rft;

    extern int initi_(double, double , double *, 
		      double *, double *, double *, double *, 
		      double *, double *, double *, double *, 
		      double *, double *);
    extern double fwv_(int , double , double *, double *, 
		       double *, double *, double *, double *);
    extern double swv_(int , double , double , double *, double *
		       , double *, double *, double *, double *, 
		       double *);
    extern double conti_n2__(double , double , double *, 
	    double *, double *, double *, double *);


/*     PROGRAM:  MODM */
/*     ------- */

/*     AUTHOR: Sid-Ahmed Boukabara */
/*     ------ */

/*     AFFILIATION: ATMOSPHERIC AND ENVIRONMENTAL RESEARCH INC. */
/*     ----------- */

/*     DATE OF CREATION : October 1998 */
/*     ---------------- */

/*     AIM: This program is aimed at the calculation of the */
/*     ---  atmospheric optical depths. The spectral validity depends */
/*          only on the region covered by the file:"spectral_lines.dat" */
/*          The components treated here are the water vapor, the */
/*          oxygen, the ozone, the nitrogen and nitrogen dioxide. */

/*     - IVC      : Flag of contin. vers.: CKD2.4(if=2)  MPMf87/s93 (if=3) */
/*     - ICP      : Flag to take(if =1) or not (if=0) the line coupling */
/*     - NWN      : Number of wavenumbers to be treated */
/*     - WN       : Vector of NWN wavenumbers [in cm-1], one should note that */
/*                  this input could be a scalar (associated with NWN=1) */
/*     - NLAY     : Number of layers to be treated. */
/*     - P        : Vector of NLAY pressures (in mbar), one should note that */
/*                  this input could be a scalar (associated with NLAY=1) */
/*     - T        : Vector of NLAY temperatures [in Kelvin] */
/*     - W_WV     : Vector of NLAY water vapor column amounts [in molecules/cm2] */
/*     - W_O2     : Vector of NLAY oxygen column amounts [in molecules/cm2] */
/*     - W_N2     : Vector of NLAY nitrogen column amounts [in molecules/cm2] */
/*     - W_OTHER  : Vector of NLAY of other species column amounts [in molecules/cm2] */
/*     - CLW      : Vector of NLAY Cloud Liquid Water amounts [in kg/m2 or mm] */
/*                  When Cloud is present, the frequency must be consistent */
/*                  with Rayleigh absorption (no scattering performed in */
/*                  monortm). */
/*     - XSLF     : Scaling factor of the self WV continuum (usually XSLF=1) */
/*     - XFRG     : Scaling factor of the foreign WV continuum (usually XFRG=1) */
/*     - XCN2     : Scaling factor of the N2 continuum (usually XCN2=1) */
/*     - O        : An array of NWNxNLAY elts containing the total optical depths */
/*                  due to all the active species [in nepers] */
/*     - OS_WV    : An array of NWNxNLAY elts containing the water vapor optical */
/*                  depth (due to self continuum), [in Nepers] */
/*     - OF_WV    : An array of NWNxNLAY elts containing the water vapor optical */
/*                  depth (due to foreign continuum), [in Nepers] */
/*     - OC_N2    : An array of NWNxNLAY elts containing the nitrogen optical */
/*                  depth (due to continuum), [in Nepers] */
/*     - O_CLW    : An array of NWNxNLAY elts containing the CLW optical */
/*                  depth , [in Nepers] */

/*     History of the modifications: */
/*     ***************************** */
/*     - written in 1999 by Sid Ahmed Boukabara, Ross Hoffman */
/* 	and Tony Clough. */
/*     - validated against ARM sondes in the */
/* 	microwave spectrum (23.8 and 31.4 GHz). SAB, 2000. */
/*     - extended to more species by Sid Ahmed Boukabara in 03/2000. */
/*     - cleaned up and commented in 2001 for first public release. */
/* 	Also put under CVS configuration management. SAB. */
/*     - Extended O2 lines to submillimeter. Extensive validation */
/* 	by comparison to Rosenkranz model and MWR data. */
/* 	Update of the LBLATM module (accepts inputs at pressure */
/* 	grid, along with altitude grid). */
/* 	Fixed the handling of N2 amount coming from LBLATM (which */
/* 	depends on the number of molecules NMOL). */
/* 	Adopted accurate constants values. */
/* 	Sid Ahmed Boukabara. Dec 14th 2001. */
/*     - Updated on January 7th 2002. ARM option (INP=2) updated and */
/*       made more efficient after Jim's comments. (INP=3) option optimized. */
/*       WV line intensities modified in the microwave (see Tony's email). */

/*     Comments should be forwarded to Sid Ahmed Boukabara (sboukaba@aer.com) */
/*     or Tony Clough (clough@aer.com). */

/* ============================================================================ */


/* TKS functions: */

/*     scaling factor (SLF cont) */
/*     scaling factor (FRG cont) */
/*     scaling factor (N2 cont) */


    w_wv__ = 0.0e0;
    w_o2__ = 0.0e0;
    w_n2__ = 0.0e0;
    w_other__ = 0.0e0;
    ret_val = 0.0e0;
    rft = 0.0e0;
    os_wv = 0.0e0;
    of_wv = 0.0e0;
    oc_n2 = 0.0e0;

/*      ---INPUTS & GENERAL CONTROL PARAMETERS */

    /*     set H2O, O2 and N2 number density to column amount [molec/cm2] */
    /*     TKSVPT = P[Pa] / T[K] */
    tksvpt    = (p * 100.0) / t;
    /*     7.242923e16 = k_B [J/K] * 1.0e-6 [m^3/cm^3] */
    w_wv__    = vmrh2o * 7.242923e16 * tksvpt;
    w_o2__    = vmro2  * 7.242923e16 * tksvpt;
    w_n2__    = vmrn2  * 7.242923e16 * tksvpt;
    w_other__ = (1.0000E0-vmrh2o-vmro2-vmrn2) * 7.242923e16 * tksvpt; 

      /*     frequency [Hz] to wave number [cm-1] */
    wn = freq / 29979245800.0;
    //cout << "CKD2.4 H2O column amounts    [molec/cm2] =" << w_wv__ << "\n";
    //cout << "CKD2.4 O2 column amounts     [molec/cm2] =" << w_o2__ << "\n";
    //cout << "CKD2.4 H2O column amounts    [molec/cm2] =" << w_n2__ << "\n";
    //cout << "CKD2.4 others column amounts [molec/cm2] =" << w_other__ << "\n";
    //cout << "freq=" << freq << " Hz,   wave num=" << wn << " cm-1\n";

/* ---------------------------------------------------------------------------- */

/*     --- INITIALIZATION ----------------------------------------- */
    initi_(p, t, &radct, &t0, &p0, &w_wv__, &w_o2__, &w_n2__, &w_other__, &
	    xn0, &xn, &xn_wv__, &rhofac);
    //cout << "CKD2.4 t0=" << t0 << "  p0=" << p0 << "\n";
    //cout << "radct   =" << radct << "\n";
    //cout << "xn0     =" << xn0 << "\n";
    //cout << "xn      =" << xn << "\n";
    //cout << "xn_wv__ =" << xn_wv__ << "\n";
    //cout << "rhofac  =" << rhofac << "\n";

    /*     --- RAD_FIELD_TERM ----------------------------------------- */
    rft = wn * tanh(radct * wn / (t * 2));
    //cout << "rft =" << rft << "\n";

    /*     --- H2O CONTINUUM TERM ------------------------------------- */

    if (ivc == 21) {
      /* CKD2.4  CONT_SELF_WV [Np/m] */
	os_wv = 1.0000e2 * swv_(2, wn, t, &t0, &w_wv__, &rft, &xn, &xn_wv__, &xn0, &xslf);
	//cout << "CKD2.4 ivc=21, H2O self cont  [in Np/m]   =" << os_wv << "\n";
	return os_wv;
    }
    if (ivc == 31) {
      /* MPMf87/s93  CONT_SELF_WV [Np/m] */
	os_wv = 1.0000e2 * swv_(3, wn, t, &t0, &w_wv__, &rft, &xn, &xn_wv__, &xn0, &xslf);
	//cout << "CKD2.4 ivc=31, H2O self cont  [in Np/m]   =" << os_wv << "\n";
	return os_wv;
    }
    if (ivc == 22) {
      /* CKD2.4  CONT_FRGN_WV [Np/m] */
	of_wv = 1.0000e2 * fwv_(2, wn, &w_wv__, &rft, &xn, &xn_wv__, &xn0, &xfrg);
	//cout << "CKD2.4 ivc=22, H2O foreign cont [in Np/m] =" << of_wv << "\n";
	return of_wv;
      }
    if (ivc == 32) {
      /* MPMf87/s93  CONT_FRGN_WV [Np/m] */
	of_wv = 1.0000e2 * fwv_(3, wn, &w_wv__, &rft, &xn, &xn_wv__, &xn0, &xfrg);
	//cout << "CKD2.4 ivc=32, H2O foreign cont [in Np/m] =" << of_wv << "\n";
	return of_wv ;
    }

    /* --- N2 CONTINUUM TERM [Np/m] ----------------------------------- */
    if (ivc == 1) {
	oc_n2 = 1.0000e2 * conti_n2__(wn, t, &t0, &w_n2__, &rft, &rhofac, &xcn2);
	//cout << "CKD2.4 ivc=1, N2 cont           [in Np/m] =" << oc_n2 << "\n";
	return oc_n2;
    }

    /* --- TOTAL ABSORPTION IN [in Np/m] --------------------------- */
    //    cout << "CKD2.4 H2O s+f cont         [in Np/m]         =" << ((os_wv+of_wv) * 1.0000e2) << "\n";
    //ret_val = ((os_wv + of_wv + oc_n2) * 1.0000e2);

L999:

    return ret_val;  // [Np/m] 
} /* artsckd_ */


/* ############################################################################ */
/*     foreign continuum functions -------------------------------------------- */
double fwv_(int ivc, double wn, double *w_wv__, double *rft, 
	    double *xn, double *xn_wv__, double *xn0, double *xfrg)
{
    /* System generated locals */
    double ret_val = 0.0e0;

    /* Local variables */
    extern double fwv24_(double , double *, double *, 
	    double *, double *, double *, double *), 
	    fwv_mpmf87s93__(double , double *, double *, 
	    double *, double *, double *, double *);

    ret_val = 0.0e0;

/*     --- CKD2.4 CONTINUUM ------------------------------------- */
    if (ivc == 2 && *w_wv__ > 0.) {
	ret_val = fwv24_(wn, w_wv__, rft, xn, xn_wv__, xn0, xfrg);
    }


/*     --- MPMf87s93 CONTINUUM ---------------------------------- */
    if (ivc == 3 && *w_wv__ > 0.) {
	ret_val = fwv_mpmf87s93__(wn, w_wv__, rft, xn, xn_wv__, xn0, xfrg);
    }

    return ret_val;
} /* fwv_ */

double fwv_mpmf87s93__(double wn, double *w_wv__, double *rft, 
		       double *xn, double *xn_wv__, double *xn0, double *xfrg)
{
    /* System generated locals */
    double ret_val=0.0e0;

    /* Local variables */
    extern double xlgr_(double *, double *);
    int i__, j;
    double x[4], fscal, xf;

    ret_val = 0.0e0;

    j = (int) ((wn - fh2ob_1.v1) / fh2ob_1.dv) + 1;

    for (i__ = 1; i__ <= 4; ++i__) {
	x[i__ - 1] = fh2oa_1.fh2o[j + i__ - 3];
    }

    xf = (wn - (fh2ob_1.v1 + fh2ob_1.dv * (double) (j - 1))) / 
	    fh2ob_1.dv;
    fscal = .8;
    ret_val = xlgr_(&xf, x) * 1e-20 * (*w_wv__ * *rft * ((*xn - *xn_wv__) / *
	    xn0)) * fscal * *xfrg;

/* L999: */
    return ret_val;
} /* fwv_mpmf87s93__ */

double fwv24_(double wn, double *w_wv__, double *rft, 
	      double *xn, double *xn_wv__, double *xn0, double *
	      xfrg)
{
    /* Initialized data */

  static double v0f1 = 350.;
  static double hwsqf1 = 4e4;
  static double betaf1 = 5e-9;
  static double factrf1 = -.7;
  static double v0f1a = 630.;
  static double hwsqf1a = 4225.;
  static double betaf1a = 2e-8;
  static double factrf1a = .75;
  static double v0f2 = 1130.;
  static double hwsqf2 = 108900.;
  static double betaf2 = 8e-11;
  static double factrf2 = -.97;
  static double v0f3 = 1975.;
  static double hwsqf3 = 62500.;
  static double betaf3 = 5e-6;
  static double factrf3 = -.65;
  
  /* System generated locals */
  double ret_val=0.0e0;
  double d__1;
  
  /* Local variables */
  extern double xlgr_(double *, double *);
  int i__, j;
  double x[4], fscal, xf, vf2, vf4, vf6;

  ret_val = 0.0e0;
  
  j = (int) ((wn - fh2ob_1.v1) / fh2ob_1.dv) + 1;
  for (i__ = 1; i__ <= 4; ++i__) {
    x[i__ - 1] = fh2oa_1.fh2o[j + i__ - 3];
  }

  xf = (wn - (fh2ob_1.v1 + fh2ob_1.dv * (double) (j - 1))) / 
        fh2ob_1.dv;

/*     ---added correction to the forgn continuum */
/* Computing 2nd power */
    d__1  = wn - v0f1;
    vf2   = d__1 * d__1;
    vf6   = vf2 * vf2 * vf2;
    fscal = factrf1 * (hwsqf1 / (vf2 + betaf1 * vf6 + hwsqf1)) + 1.;
/* Computing 2nd power */
    d__1   = wn - v0f1a;
    vf2    = d__1 * d__1;
    vf6    = vf2 * vf2 * vf2;
    fscal *= factrf1a * (hwsqf1a / (vf2 + betaf1a * vf6 + hwsqf1a)) + 1.;
/* Computing 2nd power */
    d__1 = wn - v0f2;
    vf2 = d__1 * d__1;
    vf6 = vf2 * vf2 * vf2;
    fscal *= factrf2 * (hwsqf2 / (vf2 + betaf2 * vf6 + hwsqf2)) + 1.;
/* Computing 2nd power */
    d__1 = wn - v0f3;
    vf2 = d__1 * d__1;
    vf4 = vf2 * vf2;
    fscal *= factrf3 * (hwsqf3 / (vf2 + betaf3 * vf4 + hwsqf3)) + 1.;
    ret_val = xlgr_(&xf, x) * 1e-20 * (*w_wv__ * *rft * ((*xn - *xn_wv__) / *
	    xn0)) * fscal * *xfrg;

/* L999: */
    return ret_val;
} /* fwv24_ */



/*     self continuum function ------------------------------------------------ */
double swv_(int ivc, double wn, double t, double *t0, 
	    double *w_wv__, double *rft, double *xn, double *
	    xn_wv__, double *xn0, double *xslf)
{
  /* System generated locals */
  double ret_val;

    /* Local variables */
    extern double swv_mpmf87s93__(double , double , double *
				  , double *, double *, double *, double *, 
				  double *, double *); 
    extern double swv24_(double , double , 
			 double *, double *, double *, double *, 
			 double *, double *, double *);
    
    ret_val = 0.;

/*     CKD2.4 CONTINUUM */
    if (ivc == 2 && *w_wv__ > 0.) {
/*     CNT_SLF_WV CKD2.4 */
	ret_val = swv24_(wn, t, t0, w_wv__, rft, xn, xn_wv__, xn0, xslf);
    }

    if (ivc == 3 && *w_wv__ > 0.) {
/*     MPMf87s93 CKD2.4 CONT. */
/*     CNT_SLF_WV */
	ret_val = swv_mpmf87s93__(wn, t, t0, w_wv__, rft, xn, xn_wv__, xn0, 
		xslf);
    }

    return ret_val;
} /* swv_ */



double swv24_(double wn, double t, double *t0, double *
	      w_wv__, double *rft, double *xn, double *xn_wv__, 
	      double *xn0, double *xslf)
{
    /* Initialized data */

  static double v0s1 = 0.;
  static double hwsq1 = 1e4;
  static double betas1 = 1e-4;
  static double factrs1 = .688;
  static double v0s2 = 1050.;
  static double hwsq2 = 4e4;
  static double factrs2 = -.2333;
  static double v0s3 = 1310.;
  static double hwsq3 = 14400.;
  static double betas3 = 5e-6;
  static double factrs3 = -.15;

    /* System generated locals */
    double ret_val, d__1, d__2;

    /* Local variables */
    double sfac;
    extern double xlgr_(double *, double *);
    int j;
    double x[4], xf, vs2, vs4;

/*     ---UNITS(CM**3/MOL)*1.E-20 */

    ret_val = 0.;

    j = (int) ((wn - sh2ob_1.v1) / sh2ob_1.dv) + 1;
    d__1 = s260a_1.swv260[j - 2] / sh2oa_1.swv296[j - 2];
    d__2 = (t - *t0) / (260. - *t0);
    x[0] = sh2oa_1.swv296[j - 2] * pow(d__1, d__2);
    d__1 = s260a_1.swv260[j - 1] / sh2oa_1.swv296[j - 1];
    d__2 = (t - *t0) / (260. - *t0);
    x[1] = sh2oa_1.swv296[j - 1] * pow(d__1, d__2);
    d__1 = s260a_1.swv260[j] / sh2oa_1.swv296[j];
    d__2 = (t - *t0) / (260. - *t0);
    x[2] = sh2oa_1.swv296[j] * pow(d__1, d__2);
    d__1 = s260a_1.swv260[j + 1] / sh2oa_1.swv296[j + 1];
    d__2 = (t - *t0) / (260. - *t0);
    x[3] = sh2oa_1.swv296[j + 1] * pow(d__1, d__2);
    xf = (wn - (sh2ob_1.v1 + sh2ob_1.dv * (double) (j - 1))) / 
	    sh2ob_1.dv;
    sfac = 1.;
/* Computing 2nd power */
    d__1 = wn - v0s1;
    vs2 = d__1 * d__1;
    vs4 = vs2 * vs2;
/* Computing 2nd power */
    d__1 = wn;
    sfac *= factrs1 * (hwsq1 / (d__1 * d__1 + betas1 * vs4 + hwsq1)) + 1.;
/* Computing 2nd power */
    d__1 = wn - v0s2;
    vs2 = d__1 * d__1;
    sfac *= factrs2 * (hwsq2 / (vs2 + hwsq2)) + 1.;
/* Computing 2nd power */
    d__1 = wn - v0s3;
    vs2 = d__1 * d__1;
    vs4 = vs2 * vs2;
    sfac *= factrs3 * (hwsq3 / (vs2 + betas3 * vs4 + hwsq3)) + 1.;
    ret_val = *w_wv__ * *rft * (*xn_wv__ / *xn0) * xlgr_(&xf, x) * 1e-20 * 
	    sfac * *xslf;

    return ret_val;
} /* swv24_ */

double swv_mpmf87s93__(double wn, double t, double *t0, 
	double *w_wv__, double *rft, double *xn, double *
	xn_wv__, double *xn0, double *xslf)
{
    /* System generated locals */
    double ret_val, d__1, d__2;

    /* Local variables */
    double sfac;
    extern double xlgr_(double *, double *);
    int j;
    double x[4], xf;

/*     ---UNITS(CM**3/MOL)*1.E-20 */

    ret_val = 0.;

    j = (int) ((wn - sh2ob_1.v1) / sh2ob_1.dv) + 1;
    d__1 = s260a_1.swv260[j - 2] / sh2oa_1.swv296[j - 2];
    d__2 = (t - *t0) / (260. - *t0);
    x[0] = sh2oa_1.swv296[j - 2] * pow(d__1, d__2);
    d__1 = s260a_1.swv260[j - 1] / sh2oa_1.swv296[j - 1];
    d__2 = (t - *t0) / (260. - *t0);
    x[1] = sh2oa_1.swv296[j - 1] * pow(d__1, d__2);
    d__1 = s260a_1.swv260[j] / sh2oa_1.swv296[j];
    d__2 = (t - *t0) / (260. - *t0);
    x[2] = sh2oa_1.swv296[j] * pow(d__1, d__2);
    d__1 = s260a_1.swv260[j + 1] / sh2oa_1.swv296[j + 1];
    d__2 = (t - *t0) / (260. - *t0);
    x[3] = sh2oa_1.swv296[j + 1] * pow(d__1, d__2);
    xf = (wn - (sh2ob_1.v1 + sh2ob_1.dv * (double) (j - 1))) / 
	    sh2ob_1.dv;
    sfac = 3.;
    ret_val = *w_wv__ * *rft * (*xn_wv__ / *xn0) * xlgr_(&xf, x) * 1e-20 * 
	    sfac * *xslf;

/* L999: */
    return ret_val;
} /* swv_mpmf87s93__ */

/*     --- N2 continuum ------------------------------------------------------- */
double conti_n2__(double wn, double t, double *t0, 
	double *w_n2__, double *rft, double *rhofac, double *
	xcn2)
{
    /* Initialized data */

    static double v1 = -10.;
    static double dv = 5.;
    static double ct296[73] = { 4.303e-7,4.85e-7,4.979e-7,4.85e-7,
	    4.303e-7,3.715e-7,3.292e-7,3.086e-7,2.92e-7,2.813e-7,2.804e-7,
	    2.738e-7,2.726e-7,2.724e-7,2.635e-7,2.621e-7,2.547e-7,2.428e-7,
	    2.371e-7,2.228e-7,2.1e-7,1.991e-7,1.822e-7,1.697e-7,1.555e-7,
	    1.398e-7,1.281e-7,1.138e-7,1.012e-7,9.078e-8,7.879e-8,6.944e-8,
	    6.084e-8,5.207e-8,4.54e-8,3.897e-8,3.313e-8,2.852e-8,2.413e-8,
	    2.045e-8,1.737e-8,1.458e-8,1.231e-8,1.031e-8,8.586e-9,7.162e-9,
	    5.963e-9,4.999e-9,4.226e-9,3.607e-9,3.09e-9,2.669e-9,2.325e-9,
	    2.024e-9,1.783e-9,1.574e-9,1.387e-9,1.236e-9,1.098e-9,9.777e-10,
	    8.765e-10,7.833e-10,7.022e-10,6.317e-10,5.65e-10,5.1e-10,
	    4.572e-10,4.115e-10,3.721e-10,3.339e-10,3.005e-10,2.715e-10,
	    2.428e-10 };
    static double ct220[73] = { 4.946e-7,5.756e-7,5.964e-7,5.756e-7,
	    4.946e-7,4.145e-7,3.641e-7,3.482e-7,3.34e-7,3.252e-7,3.299e-7,
	    3.206e-7,3.184e-7,3.167e-7,2.994e-7,2.943e-7,2.794e-7,2.582e-7,
	    2.468e-7,2.237e-7,2.038e-7,1.873e-7,1.641e-7,1.474e-7,1.297e-7,
	    1.114e-7,9.813e-8,8.309e-8,7.059e-8,6.068e-8,5.008e-8,4.221e-8,
	    3.537e-8,2.885e-8,2.407e-8,1.977e-8,1.605e-8,1.313e-8,1.057e-8,
	    8.482e-9,6.844e-9,5.595e-9,4.616e-9,3.854e-9,3.257e-9,2.757e-9,
	    2.372e-9,2.039e-9,1.767e-9,1.548e-9,1.346e-9,1.181e-9,1.043e-9,
	    9.11e-10,8.103e-10,7.189e-10,6.314e-10,5.635e-10,4.976e-10,
	    4.401e-10,3.926e-10,3.477e-10,3.085e-10,2.745e-10,2.416e-10,
	    2.155e-10,1.895e-10,1.678e-10,1.493e-10,1.31e-10,1.154e-10,
	    1.019e-10,8.855e-11 };

    /* System generated locals */
    double ret_val, d__1, d__2;

    /* Local variables */
    extern double xlgr_(double *, double *);
    int j;
    double x[4], xf;

/* TKS      INTEGER NPTCONTN2 */
/* TKS      INTEGER NPTCONTN2B */
/* TKS     &                 V1B,V2B,DVB */
/* TKS      DATA V2        / 350.0 / */
/* TKS      DATA NPTCONTN2 / 73    / */
/* TKS      DATA V1B, V2B, DVB, NPTCONTN2B / -10., 350., 5.0, 73 / */



    ret_val = 0.;
/* TKS  -- begin implementation of TKS */
    if (wn <= 0.) {
	ret_val = 0.;
	return ret_val;
    }
    if (wn > 350.) {
	ret_val = 0.;
	return ret_val;
    }
/* TKS  -- end implementation of TKS */

    if (*w_n2__ == 0.) {
	ret_val = 0.;
	return ret_val;
    }

    j = (int) ((wn - v1) / dv) + 1;
    d__1 = ct296[j - 2] / ct220[j - 2];
    d__2 = (t - *t0) / (220. - *t0);
    x[0] = ct296[j - 2] * pow(d__1, d__2);
    d__1 = ct296[j - 1] / ct220[j - 1];
    d__2 = (t - *t0) / (220. - *t0);
    x[1] = ct296[j - 1] * pow(d__1, d__2);
    d__1 = ct296[j] / ct220[j];
    d__2 = (t - *t0) / (220. - *t0);
    x[2] = ct296[j] * pow(d__1, d__2);
    d__1 = ct296[j + 1] / ct220[j + 1];
    d__2 = (t - *t0) / (220. - *t0);
    x[3] = ct296[j + 1] * pow(d__1, d__2);
    xf = (wn - (v1 + dv * (double) (j - 1))) / dv;
    ret_val = xlgr_(&xf, x) * 1e-20 * (*w_n2__ / .26867775 * *rft * *rhofac) *
	     *xcn2;

    return ret_val;
} /* conti_n2__ */

/*     --- 4 points Lagrange interpolation ----------------------------------- */
double xlgr_(double *xf, double *x)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    double a[4], b;




/*     with continous derivatives */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    b = *xf * .5 * (1. - *xf);
    a[0] = -b * (1. - *xf);
    a[1] = 1. - (3. - *xf * 2.) * *xf * *xf + b * *xf;
    a[2] = (3. - *xf * 2.) * *xf * *xf + b * (1. - *xf);
    a[3] = -(b * *xf);
    ret_val = a[0] * x[1] + a[1] * x[2] + a[2] * x[3] + a[3] * x[4];

/* L999: */
    return ret_val;
} /* xlgr_ */

/*     --- initializations ---------------------------------------------------- */
int initi_(double p, double t, double *radct, 
	double *t0, double *p0, double *w_wv__, double *
	w_o2__, double *w_n2__, double *w_other__, double *xn0, 
	double *xn, double *xn_wv__, double *rhofac)
{
    /* Initialized data */

    static double wvmolmass = 18.016;
    static double drymolmass = 28.97;

    double wdry, ratiomix, wvpress;

/*     [K] */
    *t0 = 296.;
/*     [hPa] */
    *p0 = 1013.25;

/*     [K/cm-1] */
    *radct   = consts_1.planck * consts_1.clight / consts_1.boltz;
    *xn0     = *p0 / (consts_1.boltz * *t0) * 1e3;
    *xn      = p / (consts_1.boltz * t) * 1e3;
    wdry     = *w_o2__ + *w_n2__ + *w_other__;
    ratiomix = *w_wv__ * wvmolmass / (wdry * drymolmass);
    wvpress  = ratiomix / (ratiomix + wvmolmass / drymolmass) * p;
    *xn_wv__ = wvpress / (consts_1.boltz * t) * 1e3;
    *rhofac  = *w_n2__ / (wdry + *w_wv__) * (p / *p0) * (273.15 / t);

/* L999: */
    return 0;
} /* initi_ */

/* ############################################################################ */
/* ---Block data to be consistent with LBLRTM/LBLATM */
/* Subroutine */ int phys_consts__(void)
{
    return 0;
} /* phys_consts__ */


/*     Pi was obtained from PI = 2.*ASIN(1.) */
/* --------------------------------------------- */
/*       Constants from NIST 01/11/2002 */
/* --------------------------------------------- */
/* --------------------------------------------- */
/*       units are generally cgs */
/*       The first and second radiation constants are taken from NIST. */
/*       They were previously obtained from the relations: */
/*       RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07 */
/*       RADCN2 = PLANCK*CLIGHT/BOLTZ */
/* --------------------------------------------- */

/* ############################################################################ */
/* Subroutine */ int bsa296_(void)
{
    return 0;
} /* bsa296_ */







/* Subroutine */ int bsb296_(void)
{
    return 0;
} /* bsb296_ */







/* ---------------------------------------------------------------------------- */


/* Subroutine */ int bs260a_(void)
{
    return 0;
} /* bs260a_ */








/* Subroutine */ int bs260b_(void)
{
    return 0;
} /* bs260b_ */







/* ---------------------------------------------------------------------------- */


/* Subroutine */ int bfh2oa_(void)
{
    return 0;
} /* bfh2oa_ */








/* Subroutine */ int bfh2ob_(void)
{
    return 0;
} /* bfh2ob_ */


// ---------------------- end of monortm CKD F77  code --------------------------
