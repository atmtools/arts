/**
   \file   continua.cc

   The following continua parameterizations are implemented:
   1) H2O-H2O: P. W. Rosenkranz, Radio Science, Vol. 33, No 4, Pages 919-928, 1998.
   2) H2O-air: P. W. Rosenkranz, Radio Science, Vol. 33, No 4, Pages 919-928, 1998.
                            and  Radio Science, Vol. 34, No 4, Page  1025,    1999.
   3) O2-air : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
                                "Atmospheric Remote Sensing by Microwave Radiometry",
                                John Wiley & Sons, Inc., 1993
   4) N2-N2  : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
                                "Atmospheric Remote Sensing by Microwave Radiometry",
                                John Wiley & Sons, Inc., 1993
   5) CO2-CO2: P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
                                "Atmospheric Remote Sensing by Microwave Radiometry",
                                John Wiley & Sons, Inc., 1993
   6) CO2-N2 : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
                                "Atmospheric Remote Sensing by Microwave Radiometry",
                                John Wiley & Sons, Inc., 1993

   A) Suspended water droplet absorption parameterization from MPM93 model
      H. J. Liebe and G. A. Hufford and M. G. Cotton,
      "Propagation modeling of moist air and suspended water/ice
      particles at frequencies below 1000 GHz",
      AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
      
   B) Ice crystal absorption parameterization from MPM93 model
      H. J. Liebe and G. A. Hufford and M. G. Cotton,
      "Propagation modeling of moist air and suspended water/ice
      particles at frequencies below 1000 GHz",
      AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 

   i) MPM87 H2O absorption model (line and continuum) according to 
      Liebe, Radio Science, 20(5), 1985, 1069.

  ii) MPM89 H2O absorption model (line and continuum) according to 
      Liebe, Int. J. Infrared and Millimeter Waves, 10(6), 1989, 631

 iii) MPM93 H2O absorption model (line and continuum) according to 
      . J. Liebe and G. A. Hufford and M. G. Cotton,
      "Propagation modeling of moist air and suspended water/ice
      particles at frequencies below 1000 GHz",
      AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 

  iv) Cruz-Pol H2O absorption model (line and continuum) according to 
      S. L. Cruz-Pol et al., Radio Science, 33(5), 1319, 1998

   v) Rosenkranz H2O absorption model (line and continuum) according to 
      P. W. Rosenkranz., Radio Science, 33(4), 919, 1998

   ------------------------------------------------------------------

   The following unit conversions are used:
   (SI units: meter, kilogram, second, ampere, Kelvin, candela) 
   x g/cm3 = y kg/m3    <===>      y = x * 1.00e3
   x g/m3  = y kg/m3    <===>      y = x * 1.00e-3
   x GHz   = y Hz       <===>      y = x * 1.00e9
   x 1/GHz = y 1/Hz     <===>      y = x * 1.00e-9
   x hPa   = y Pa       <===>      y = x * 1.00e2
   x 1/hPa = y 1/Pa     <===>      y = x * 1.00e-2
   x 1/cm  = y 1/m      <===>      y = x * 1.0e2
   x 1/km  = y 1/m      <===>      y = x * 1.00e-3
   x dB    = y Np       <===>      y = x / (10.0 * log10(e))
   x dB/km = y 1/m      <===>      y = x * 1.00e-3 / (10.0 * log10(e))
   x Np/km = y 1/m      <===>      y = x * 1.00e-3

   \author Thomas Kuhn
   \date   2001-03-27
*/

#include "vecmat.h"
#include "absorption.h"
#include "messages.h"
#include "continua.h"

/**
   Calculates the continuum according to the simple
   empirical function as formulated by Rosenkranz98:

   self broadened continuum:
   xsec += C * (300/t_abs)^(x+3) * f_mono^2 * (300/t_abs)^3 * p_abs^2 * vmr

   foreign broadened continuum:
   pdry = p_abs * (1-vmr)
   xsec += C * (300/t_abs)^(x+3) * f_mono^2 * (300/t_abs)^3 * p_abs * pdry

   See the list of parameters below for the meaning of these
   variables. The equation is slightly modified from Rosenkranz's
   original, because we want to return only the cross section, not the
   absorption coefficient. (Only vmr, not vmr^2 at the end of the
   equation.) 

   This xsec is understood to parameterize the difference between
   observed absorption cross section and the one calculated from the
   pure line spectrum. The continuum is added to the previous content
   of xsec!

   \retval xsec  Absorption cross section, defined such that the
                 absorption coefficient alpha is:<br>
                 alpaha [1/m] = xsec * VMR.<br>
		 The functions adds to xsec, rather than replacing the
		 previous content. 

   \param  C       Continuum coefficient.
   \param  x       Temperature exponenet.
   \param  f_mono  Frequency grid.
   \param  p_abs   Pressure grid.
   \param  t_abs   Temperatures associated with p_abs.
   \param  vmr     Volume mixing ratio of the calculated species (H2O).

   References:
   
   P.W. Rosenkranz, `Water vapor microwave continuum absorption: A
   comparison of measurements and models', Radio Science, Vol. 33, No
   4, Pages 919-928, July-August 1998.


   \author Stefan Buehler
   \date   2001-01-17
*/
//
// #################################################################################
//
// 1) H2O-H2O: P. W. Rosenkranz, Radio Science, Vol. 33, No 4, Pages 919-928, 1998.
void Rosenkranz_h2o_self_continuum( MATRIX&           xsec,
				    Numeric	      C,
				    Numeric	      x,
				    const VECTOR&     f_mono,
				    const VECTOR&     p_abs,
				    const VECTOR&     t_abs,
				    const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( size_t i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The second vmr of H2O will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy =
	C * pow( 300./t_abs[i], x+3. ) * pow( p_abs[i], 2 ) * vmr[i];

      // Loop over frequency grid:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += dummy * pow( f_mono[s], 2 );
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
//   2) H2O-air: P. W. Rosenkranz, Radio Science, Vol. 33, No 4, Pages 919-928, 1998.
//                            and  Radio Science, Vol. 34, No 4, Page  1025,    1999.
void Rosenkranz_h2o_foreign_continuum( MATRIX&           xsec,
				       Numeric	         C,
				       Numeric	         x,
				       const VECTOR&     f_mono,
				       const VECTOR&     p_abs,
				       const VECTOR&     t_abs,
				       const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // Dry air partial pressure: p_dry := p_tot - p_h2o.
      Numeric pdry  = p_abs[i] * (1.000e0-vmr[i]);
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The vmr of H2O will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy = C * pow( 300./t_abs[i], x+3. ) * p_abs[i] * pdry;

      // Loop frequency:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += dummy * pow( f_mono[s], 2 );
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
// 2b) H2O-air: Coefficients are from Liebe et al., AGARD CP-May93, Paper 3/1-10
//
Numeric MPMLineShapeFunction( Numeric gamma, 
				Numeric fl, 
				Numeric f)
{
  /*
    this routine calculates the line shape function of Van Vleck and Weisskopf
    with the factor (f/f_o)^1. for the MPM pseudo continuum line.

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
// MPM93 H2O pseudo continuum line parameters:
// see publication side of National Telecommunications and Information Administration
//   http://www.its.bldrdoc.gov/pub/all_pubs/all_pubs.html
// and ftp side for downloading the MPM93 original source code:
//   ftp://ftp.its.bldrdoc.gov/pub/mpm93/
void MPM93_h2o_continuum( MATRIX&           xsec,
			  Numeric	    MPM93fopcl, // default: 1780.0*10^9 Hz
			  Numeric	    MPM93b1pcl, // default: 22300.0 Hz/Pa
			  Numeric	    MPM93b2pcl, // default: 0.952
			  Numeric	    MPM93b3pcl, // default: 17.6*10^4 Hz/Pa
			  Numeric	    MPM93b4pcl, // default: 30.5
			  Numeric	    MPM93b5pcl, // default: 2
			  Numeric	    MPM93b6pcl, // default: 5
			  const VECTOR&     f_mono,
			  const VECTOR&     p_abs,
			  const VECTOR&     t_abs,
			  const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  

  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      Numeric th = 300.0 / t_abs[i];
      // the vmr of H2O will be multiplied at the stage of absorption calculation:
      // abs / vmr * xsec.
      Numeric strength =  MPM93b1pcl * p_abs[i] * pow( th, 3.5 ) * exp(MPM93b2pcl * (1 - th));
      Numeric gam =  MPM93b3pcl * 0.001 * 
	           ( MPM93b4pcl * p_abs[i] * vmr[i]       * pow( th, MPM93b6pcl ) +  
	                          p_abs[i]*(1.000-vmr[i]) * pow( th, MPM93b5pcl ) );
      // Loop frequency:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += 0.182 * 0.001 / (10.000*log10(2.718281828))  * 
	                f_mono[s] * strength * 
	                MPMLineShapeFunction(gam, MPM93fopcl, f_mono[s]); 
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
//   3) O2-air : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
//               "Atmospheric Remote Sensing by Microwave Radiometry",
//               John Wiley & Sons, Inc., 1993. Also stated in 
//               Liebe et al. JQSRT, Vol 48, Nr 5/6, pp. 629-643, 1992.
//               Default continuum parameters are  C=1.6E-17*10E-9,  x=0.8
void Rosenkranz_o2_continuum( MATRIX&           xsec,
			      Numeric		C, // default: 1.108*10^-14 K^2/(Hz*Pa*m)
			      Numeric		x, // default: 0.8
			      const VECTOR&  	f_mono,
			      const VECTOR&  	p_abs,
			      const VECTOR&  	t_abs,
			      const VECTOR&     h2o_abs,
			      const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // loop over all pressure levels:
  for ( size_t i=0; i<n_p; ++i )
    {
      Numeric TH = 300.00 / t_abs[i];        // relative temperature  [1]
      
      Numeric ph2o  = p_abs[i] * h2o_abs[i];  // water vapor partial pressure [Pa]
      Numeric pdry  = p_abs[i] - ph2o;       // dry air partial pressure     [Pa]
      // pseudo broadening term [Hz]
      Numeric gamma = 5600.000 * (pdry * pow( TH, x ) + 1.100 * ph2o * TH); 

      // Loop over frequency grid:
      for ( size_t s=0; s<n_f; ++s )
	{
	  // division by vmr of O2 is necessary because of the absorption calculation
          // abs = vmr * xsec.
	  xsec[s][i] += C * p_abs[i] * gamma * pow( f_mono[s], 2 ) / 
                          ( pow( t_abs[i], 2 ) * ( pow( f_mono[s], 2 ) + pow( gamma, 2 ) ) );
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
// 4) N2-N2  : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
//    "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
//
void Rosenkranz_n2_self_continuum( MATRIX&           xsec,
				   Numeric	     C, // default: 1.05*10^-38 1/(Pa^2*Hz^2*m)
				   Numeric	     x, // default: 3.55
				   const VECTOR&     f_mono,
				   const VECTOR&     p_abs,
				   const VECTOR&     t_abs,
				   const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( size_t i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The second vmr of N2 will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy =
	C * pow( 300./t_abs[i], x ) * pow( p_abs[i], 2 ) * vmr[i];

      // Loop over frequency grid:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += dummy * pow( f_mono[s], 2 );
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
// 5) CO2-CO2: P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
// "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
//
void Rosenkranz_co2_self_continuum( MATRIX&           xsec,
				    Numeric	      C, // default: 7.43*10^-37 1/(Pa^2*Hz^2*m)
				    Numeric	      x, // default: 5.08
				    const VECTOR&     f_mono,
				    const VECTOR&     p_abs,
				    const VECTOR&     t_abs,
				    const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop over pressure/temperature grid:
  for ( size_t i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The second vmr of CO2 will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy =
	C * pow( 300./t_abs[i], x ) * pow( p_abs[i], 2 ) * vmr[i];

      // Loop over frequency grid:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += dummy * pow( f_mono[s], 2 );
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
// 6) CO2-N2 : P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
//    "Atmospheric Remote Sensing by Microwave Radiometry", John Wiley & Sons, Inc., 1993
//
void Rosenkranz_co2_foreign_continuum( MATRIX&           xsec,
				       Numeric	         C, // default: 2.71*10^-37 1/(Pa^2*Hz^2*m)
				       Numeric	         x, // default: 4.7
				       const VECTOR&     f_mono,
				       const VECTOR&     p_abs,
				       const VECTOR&     t_abs,
				       const VECTOR&     n2_abs,
				       const VECTOR&     vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // Dummy scalar holds everything except the quadratic frequency dependence.
      // The vmr of CO2 will be multiplied at the stage of absorption 
      // calculation: abs = vmr * xsec.
      Numeric dummy = C * pow( 300./t_abs[i], x ) * p_abs[i] * p_abs[i] * n2_abs[i];

      // Loop frequency:
      for ( size_t s=0; s<n_f; ++s )
	{
	  xsec[s][i] += dummy * pow( f_mono[s], 2 );
	  //	  cout << "xsec[" << s << "][" << i << "]: " << xsec[s][i] << "\n";
	}
    }
}
//
// #################################################################################
//
// saturation water vapor pressure over liquid water,
// calculated according to Goff and Gratch formula.
// The saturation water vapor pressure is in units of Pa
Numeric WVSatPressureLiquidWater(Numeric t)
{
  Numeric theta  = 373.16 / (t + 273.16);
  Numeric es     = 100.000 * 
                   ( -7.90298 * (theta-1.000) +
		     5.02808 * log10(theta) -
		     1.3816e-7 * ( pow( 10.00, (11.344*(1.00-(1.00/theta))) ) - 1.000 ) +
		     8.1328e-3 * ( pow( 10.00, (-3.49149*(theta-1.00))) - 1.000) +
		     log10(1013.246) );
  return es;
}
//
// #################################################################################
//
// saturation water vapor pressure over ice,
// calculated according to Goff and Gratch formula.
// The saturation water vapor pressure is in units of Pa
Numeric WVSatPressureIce(Numeric t)
{
  Numeric theta  = 273.16 / (t + 273.16);
  Numeric es     = 100.000 * 
                  (-9.09718  * (theta-1.000) -
		   3.56654  * log10(theta)  +
		   0.876793 * (1.000-(1.000/theta)) +
		   log10(6.1071) );
  return es;
}
// #################################################################################
//
//   A) cloud and fog absorption parameterization from MPM93 model
//      input parameters:
//      w :  suspended water droplet density, valid range: 0-10.00e-3 kg/m3
//      m :  specific weight of the droplet,  fixed value:     1.00e3 kg/m3
//      The internal numerical values (and units) are the same as in MPM93
//
void MPM93WaterDropletAbs( MATRIX&           xsec,
			   Numeric	        w, // suspended water droplet density
			   Numeric	        m, // specific droplet weight 
			   const VECTOR&   f_mono,
			   const VECTOR&    p_abs,
			   const VECTOR&    t_abs,
			   const VECTOR&      vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Check that suspended water droplet density and specific weight of the droplet
  // are in the correct limits
  if ((w < 0.00) || (w > 10.00e-3) || (fabs(m-1.000e3)> 0.100)) {
    ostringstream os;
    os << "MPM93WaterDropletAbs: \n"
       << "suspended liquid water particle density,       valid range: 0-10.00e-3 kg/m3,  w=" << w << "\n"
       << "specific weight of the liquid water particle,  fixed value:     1.00e3 kg/m3,  m=" << m << "\n"
       << ".";
    throw runtime_error(os.str());
    return;
  }

  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // water vapor saturation pressure over liquid water [Pa]
      // Numeric es       = WVSatPressureLiquidWater(t_abs[i]);
      // water vapor partial pressure [Pa]
      // Numeric e        = p_abs[i] * vmr[i];      
      // relative humidity [1]
      // Numeric RH       = e / es;
  
      // relative inverse temperature [1]
      Numeric theta    = 300.000 / t_abs[i];
      // relaxation frequencies [GHz]
      Numeric gamma1   = 20.20 - 146.40*(theta-1.000) + 316.00*(theta-1.000)*(theta-1.000);
      Numeric gamma2   = 39.80 * gamma1; 
      // static and high-frequency permittivities
      Numeric epsilon0 = 103.30 * (theta-1.000) + 77.66;
      Numeric epsilon1 = 0.0671 * epsilon0;
      Numeric epsilon2 = 3.52;

      // Loop frequency:
      for ( size_t s=0; s<n_f; ++s )
	{
	  // real part of the complex permittivity of water (double-debye model)
	  Numeric Reepsilon  = epsilon0 - 
	                       pow((f_mono[s]*1.000e-9),2) *
	                       ( ((epsilon0-epsilon1)/
				  (pow((f_mono[s]*1.000e-9),2) + pow(gamma1,2))) + 
                                 ((epsilon1-epsilon2)/
				  (pow((f_mono[s]*1.000e-9),2) + pow(gamma2,2))) );
	  // imaginary part of the complex permittivity of water (double-debye model)
	  Numeric Imepsilon  = (f_mono[s]*1.000e-9) *
	                       ( (gamma1*(epsilon0-epsilon1)/
				  (pow((f_mono[s]*1.000e-9),2) + pow(gamma1,2))) + 
                                 (gamma2*(epsilon1-epsilon2)/
				 (pow((f_mono[s]*1.000e-9),2) + pow(gamma2,2))) );
	  // the imaginary part of the complex refractivity of suspended liquid water particle [ppm]
	  // In MPM93 w is in g/m3 and m is in g/cm3. Because of the units used in arts,
          // a factor of 1.000e6 must be multiplied with the ratio (w/m):
          // MPM93: (w/m)_MPM93  in   (g/m3)/(g/cm3)
          // arts:  (w/m)_arts   in  (kg/m3)/(kg/m3)
          // ===> (w/m)_MPM93 = 1.0e6 * (w/m)_arts
          // the factor of 1.0e6 is included below in the constant 41.90705.
	  Numeric ImNw = 1.500 * (w/m) * 
 	                ( 3.000 * Imepsilon / ( pow((Reepsilon+2.000),2) + pow(Imepsilon,2) ) );
	  // liquid water particle absorption cross section [1/m]
	  // The vmr of H2O will be multiplied at the stage of absorption 
	  // calculation: abs = vmr * xsec.
          // 41.90705 = (0.182 * 0.001 / (10.000*log10(2.718281828))) * 1.000e6
	  xsec[s][i] += 41.90705 * (f_mono[s]*1.000e-9) * ImNw / vmr[i];
	}
    }

}
//
// #################################################################################
//
//   A) ice particle absorption parameterization from MPM93 model
//      input parameters:
//      w :  suspended ice particle density,    valid range: 0-10.0e-3 kg/m3
//      m :  specific weight of ice particles,  fixed value:   0.916e3 kg/m3
//      The internal numerical values (and units) are the same as in MPM93
//
void MPM93IceCrystalAbs( MATRIX&           xsec,
			 Numeric	      w, // suspended ice particle density
			 Numeric	      m, // specific ice particles weight 
			 const VECTOR&   f_mono,
			 const VECTOR&    p_abs,
			 const VECTOR&    t_abs,
			 const VECTOR&      vmr	 )
{
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Check that suspended water droplet density and specific weight of the droplet
  // are in the correct limits
  if ((w < 0.00) || (w > 10.00e-3) || (fabs(m-0.916e3)> 0.100)) {
    ostringstream os;
    os << "MPM93IceCrystalAbs: \n"
       << "suspended ice particle density,       valid range: 0-10.0e-3 kg/m3,   w=" << w << "\n"
       << "specific weight of the ice particle,  fixed value:   0.916e3 kg/m3,   m=" << m << "\n"
       << ".";
    throw runtime_error(os.str());
    return;
  }


  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // water vapor saturation pressure over ice [Pa]
      // Numeric es = WVSatPressureIce(t_abs[i]);
      // water vapor partial pressure [Pa]
      // Numeric e  = p_abs[i] * vmr[i];
      // relative humidity [1]
      // Numeric RH = e / es;
  


      // if the temperature is below 0 degree C (273 K) then 
      // absorption due to ice particles is added:
      if (t_abs[i] < 273.15)
	{
	  // relative inverse temperature [1]
	  Numeric theta = 300.000 / t_abs[i];	
	  // inverse frequency T-dependency function [Hz]
	  Numeric ai = (62.000 * theta - 11.600) * exp(-22.100 * (theta-1.000)) * 1.000e-4;
	  // linear frequency T-dependency function [1/Hz]
	  Numeric bi = 0.542e-6 * 
	               ( -24.17 + (116.79/theta) + pow((theta/(theta-0.9927)),2) );

	  // Loop frequency:
	  for ( size_t s=0; s<n_f; ++s )
	    {
	      // real part of the complex permittivity of ice
	      Numeric Reepsilon  = 3.15;
	      // imaginary part of the complex permittivity of water
	      Numeric Imepsilon  = ( ( ai/(f_mono[s]*1.000e-9) ) +
				     ( bi*(f_mono[s]*1.000e-9) ) );
	      // the imaginary part of the complex refractivity of suspended ice particles.
	      // In MPM93 w is in g/m3 and m is in g/cm3. Because of the units used in arts,
	      // a factor of 1.000e6 must be multiplied with the ratio (w/m):
	      // MPM93: (w/m)_MPM93  in   (g/m3)/(g/cm3)
	      // arts:  (w/m)_arts   in  (kg/m3)/(kg/m3)
	      // ===> (w/m)_MPM93 = 1.0e6 * (w/m)_arts
	      // the factor of 1.0e6 is included below in the constant 41.90705.
	      Numeric ImNw = 1.500 * (w/m) * 
 	                     ( 3.000 * Imepsilon / ( pow((Reepsilon+2.000),2) + pow(Imepsilon,2) ) );
	      // ice particle absorption cross section [1/m]
	      // The vmr of H2O will be multiplied at the stage of absorption 
	      // calculation: abs = vmr * xsec.
	      // 41.90705 = (0.182 * 0.001 / (10.000*log10(2.718281828))) * 1.000e6
	      xsec[s][i] += 41.90705 * (f_mono[s]*1.000e-9) * ImNw / vmr[i];
	    }
	}
    }
  return;
}
//
// #################################################################################
// ####################### additional line absorption models #######################
// implemented models are: 
//   1.) MPM87:    H. J. Liebe and D. H. Layton, NITA Report 87-224, 1987
//   2.) MPM89:    H. J. Liebe, Int. J. Infrared and Mill. Waves, 10(6), 631, 1989
//   3.) MPM93:    H. J. Liebe and G. A. Hufford and M. G. Cotton,
//                 AGARD 52nd Specialists Meeting of the Electromagnetic Wave
//                 Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 
//   4.) Cruz-Pol: S. L. Cruz-Pol et al. Radio Science 33(5), 1319, 1998
// 
// #################################################################################
//
void MPM87H2OAbsModel( MATRIX&           xsec,
		       Numeric             CC, // continuum term scale factor
		       Numeric             CL, // line term scale factor
		       Numeric             CW, // line broadening scale factor
		       const VECTOR&   f_mono,
		       const VECTOR&    p_abs,
		       const VECTOR&    t_abs,
		       const VECTOR&      vmr )
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
    {   183.310074,    2.3000,  0.653,   31.64^e-3},
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

  // number of lines of liebe line catalogue (30 lines)
  const size_t i_first = 0;
  const size_t i_last  = 29;
  
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in [hPa] therefore the factor 0.01)
  for ( size_t i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(H2O) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
	{
	  // relative inverse temperature [1]
	  Numeric theta = (300.0 / t_abs[i]);
	  // H2O partial pressure [kPa]
	  Numeric pwv   = 0.001 * p_abs[i] * vmr[i];
	  // dry air partial pressure [kPa]
	  Numeric pda   = (0.001 * p_abs[i]) - pwv;
	  // H2O continuum absorption [dB/km/GHz2] like in the original MPM87
	  Numeric absc  = pwv * pow(theta, 3.0) * 1.000e-5 *
	                  ( (0.113 * pda) + (3.57 * pwv * pow(theta, 7.5)) );

	  // Loop over input frequency
	  for ( size_t s=0; s<n_f; ++s )
	    {
	      // input frequency in [GHz]
	      Numeric ff   = f_mono[s] * 1.00e-9; 
	      // H2O line contribution at position f
	      Numeric absl = 0.000;
	      
	      // Loop over MPM89 H2O spectral lines
	      for ( size_t l = i_first; l <= i_last; ++l )
		{
		  // line strength [kHz]
		  Numeric strength = CL * mpm87[l][1] * pwv * 
	                             pow(theta,3.5) * exp(mpm87[l][2]*(1.000-theta));
		  // line broadening parameter [GHz]
		  Numeric gam      = CW * mpm87[l][3] * 
 		                     ( (4.80 * pwv * pow(theta, 1.1)) + 
				       (       pda * pow(theta, 0.6)) );
		  // effective line width with Doppler broadening [GHz]
		  // gam              = sqrt(gam*gam + (2.14e-12 * mpm87[l][0] * mpm87[l][0] / theta));
		  // H2O line absorption [dB/km/GHz] like in the original MPM87
		  absl            += strength * MPMLineShapeFunction(gam, mpm87[l][0], ff); 
		}
	      // xsec = abs/vmr [1/m]
	      // 4.1907e-5 = 0.230259 * 1.0e-3 * 0.1820     (1/(10*log(e)) = 0.230259)
	      xsec[s][i]  += 4.1907e-5 * ff * ( absl + (CC * absc * ff) ) / vmr[i];
	    }
	}
    }
  return;
}
//
// #################################################################################
//
void MPM89H2OAbsModel( MATRIX&           xsec,
		       Numeric             CC, // continuum term scale factor
		       Numeric             CL, // line term scale factor
		       Numeric             CW, // line broadening scale factor
		       const VECTOR&   f_mono,
		       const VECTOR&    p_abs,
		       const VECTOR&    t_abs,
		       const VECTOR&      vmr )
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

  // number of lines of liebe line catalogue (30 lines)
  const size_t i_first = 0;
  const size_t i_last  = 29;
  
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature (pressure in [hPa] therefore the factor 0.01)
  for ( size_t i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(H2O) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
	{
	  // relative inverse temperature [1]
	  Numeric theta = (300.0 / t_abs[i]);
	  // H2O partial pressure [kPa]
	  Numeric pwv   = 0.001 * p_abs[i] * vmr[i];
	  // dry air partial pressure [kPa]
	  Numeric pda   = (0.001 * p_abs[i]) - pwv;
	  // H2O continuum absorption [dB/km/GHz2] like in the original MPM89
	  Numeric absc  = pwv * pow(theta, 3.0) * 1.000e-5 *
	                  ( (0.113 * pda) + (3.57 * pwv * pow(theta, 7.5)) );

	  // Loop over input frequency
	  for ( size_t s=0; s<n_f; ++s )
	    {
	      // input frequency in [GHz]
	      Numeric ff   = f_mono[s] * 1.00e-9; 
	      // H2O line contribution at position f 
	      Numeric absl = 0.000;
	      
	      // Loop over MPM89 spectral lines:
	      for ( size_t l = i_first; l <= i_last; ++l )
		{
		  // line strength [kHz]
		  Numeric strength = CL * mpm89[l][1] * pwv * 
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
		  absl            += strength * MPMLineShapeFunction(gam, mpm89[l][0], ff); 
		}
	      // xsec = abs/vmr [1/m]
	      // 4.1907e-5 = 0.230259 * 1.0e-3 * 0.1820    (1/(10*log(e)) = 0.230259)
	      xsec[s][i] += 4.1907e-5 * ff * ( absl + (CC * absc * ff) ) / vmr[i];
	    }
	}
    }
  return;
}
//
// #################################################################################
//
void MPM93H2OAbsModel( MATRIX&           xsec,
		       Numeric             CC, // continuum term scale factor
		       Numeric             CL, // line term scale factor
		       Numeric             CW, // line broadening scale factor
		       const VECTOR&   f_mono,
		       const VECTOR&    p_abs,
		       const VECTOR&    t_abs,
		       const VECTOR&      vmr )
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
    {   547.676440,    0.97010,  0.114,   2.600,   4.50,  0.70,  1.00},
    {   552.020960,    1.47700,  0.114,   2.600,   4.50,  0.70,  1.00},
    {   556.936002,   48.74000,  0.159,   3.210,   4.11,  0.69,  1.00},
    {   620.700807,    0.50120,  2.200,   2.438,   4.68,  0.71,  0.68},
    {   645.866155,    0.00713,  8.580,   1.800,   4.00,  0.60,  0.50},
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

  
  // number of lines of liebe line catalogue (0-33 lines, 34 cont. pseudo line)
  const size_t i_first = 0;
  const size_t i_last  = 34;
  
  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in hPa therefore the factor 0.01)
  for ( size_t i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(H2O) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
	{
	  // relative inverse temperature [1]
	  Numeric theta    = (300.0 / t_abs[i]);
	  // H2O partial pressure [hPa]
	  Numeric pwv      = 0.01 * p_abs[i] * vmr[i];
	  // dry air partial pressure [hPa]
	  Numeric pda      = (0.01 * p_abs[i]) - pwv;
	  // Loop over MPM93 spectral lines:
	  
	  // Loop over input frequency
	  for ( size_t s=0; s<n_f; ++s )
	    {
	      // input frequency in [GHz]
	      Numeric ff = f_mono[s] * 1.00e-9; 

	      for ( size_t l = i_first; l <= i_last; ++l )
		{
		  // line strength [ppm]. The missing vmr of H2O will be multiplied 
		  // at the stage of absorption calculation: abs / vmr * xsec.
		  Numeric strength = CL * mpm93[l][1] * pwv * 
	                             pow(theta, 3.5)  * exp(mpm93[l][2]*(1.0-theta));
		  // line broadening parameter [GHz]
		  Numeric gam      = CW * mpm93[l][3] * 0.001 * 
	                             ( (mpm93[l][4] * pwv * pow(theta, mpm93[l][6])) +  
	                               (              pda * pow(theta, mpm93[l][5])) );
		  // Doppler line width [GHz]
		  // Numeric gamd     = 1.46e-6 * mpm93[l][0] / sqrt(theta);
		  // effective line width [GHz]
		  //gam              = 0.535 * gam + sqrt(0.217*gam*gam + gamd*gamd); 
		  // absorption [dB/km] like in the original MPM93
		  Numeric abs = strength * MPMLineShapeFunction(gam, mpm93[l][0], ff); 
		  if (l == 34)
		    {
		      // H2O pseudo continuum line
		      // 4.1907e-5 = 0.230259 * 1.0e-3 * 0.1820    (1/(10*log(e)) = 0.230259)
		      // xsec = abs/vmr [1/m]
		      xsec[s][i]   += 4.1907e-5 * CC * ff * abs / vmr[i];
		    } else 
		      {
			// H2O spectral line
			// xsec = abs/vmr [1/m]
			xsec[s][i]   += 4.1907e-5 * ff * abs / vmr[i];
		      }
		}
	    }
	}
    }
  return;
}
//
// #################################################################################
//
void CP98H2OAbsModel( MATRIX&           xsec,
		      Numeric             CC, // continuum scale factor
		      Numeric             CL, // line strength scale factor
		      Numeric             CW, // line broadening scale factor
		      const VECTOR&   f_mono,
		      const VECTOR&    p_abs,
		      const VECTOR&    t_abs,
		      const VECTOR&      vmr )
{
  //
  // Coefficients are from S. L. Cruz-Pol et al., Radio Science, 33(5), 1319, 1998
  // nominal values for the scale factors:
  //  CC = 1.2369 +/- 0.155
  //  CL = 1.0639 +/- 0.016
  //  CW = 1.0658 +/- 0.0096

  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature (pressure in [hPa] therefore the factor 0.01)
  for ( size_t i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(H2O) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
	{
	  // relative inverse temperature [1]
	  Numeric theta = (300.0 / t_abs[i]);
	  // H2O partial pressure [hPa]
	  Numeric pwv   = 0.01 * p_abs[i] * vmr[i];
	  // dry air partial pressure [hPa]
	  Numeric pda   = (0.01 * p_abs[i]) - pwv;
	  // line strength
	  Numeric TL    = CL * 0.0109 * pwv * pow(theta,3.5) * exp(2.143*(1.0-theta));
	  // line broadening parameter [GHz]
	  Numeric gam   = CW * 0.002784 *  
	                  ( (pda * pow(theta,0.6)) + (4.80 * pwv * pow(theta,1.1)) );
	  // continuum term
	  Numeric TC    = CC * pwv * pow(theta, 3.0) * 1.000e-7 * 
	                  ( (0.113 * pda) + (3.57 * pwv * pow(theta,7.5)) );
	  // Loop over input frequency
	  for ( size_t s=0; s<n_f; ++s )
	    {
	      // input frequency in [GHz]
	      Numeric ff  = f_mono[s] * 1.00e-9; 
	      Numeric TSf = MPMLineShapeFunction(gam, 22.235080, ff); 
	      // xsec = abs/vmr [1/m] (Cruz-Pol model in [Np/km])
	      xsec[s][i] += 4.1907e-5 * ff * ( (TL * TSf) + (ff * TC) ) / vmr[i];
	    }
	}
    }
  return;
}
//
// #################################################################################
//
void PWR98H2OAbsModel( MATRIX&           xsec,
		       Numeric             CC, // continuum term scale factor
		       Numeric             CL, // line term scale factor
		       Numeric             CW, // line broadening scale factor
		       const VECTOR&   f_mono,
		       const VECTOR&    p_abs,
		       const VECTOR&    t_abs,
		       const VECTOR&      vmr )
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
  const Numeric PWRfl[15] = {  22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 
                              439.1508, 443.0183, 448.0011, 470.889,  474.6891, 
                              488.4911, 556.936,  620.7008, 752.0332, 916.1712 };
  // line intensities at 300K
  const Numeric PWRs1[15] = { 1.31e-14,  2.273e-12, 8.036e-14, 2.694e-12, 2.438e-11, 
                              2.179e-12, 4.624e-13, 2.562e-11, 8.369e-13, 3.263e-12, 
                              6.659e-13, 1.531e-9,  1.707e-11, 1.011e-9,  4.227e-11 };
  // T coeff. of intensities
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

  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );

  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(H2O) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
	{
	  // water vapor partial pressure [hPa]
	  Numeric pvap = 0.01 * p_abs[i] * vmr[i];
	  // dry air partial pressure [hPa]
	  Numeric pda  = (0.01 *p_abs[i]) - pvap;
	  // Rosenkranz number density  (Rosenkranz H2O mass density in [g/m3])
	  // [g/m3]    =  [g*K / Pa*m3]  *  [Pa/K]
	  // rho       =   (M_H2O / R)   *  (P_H2O / T)
	  // rho       =      2.1667     *  p_abs * vmr / t_abs
	  // den       = 3.335e16 * rho
	  Numeric den  = 3.335e16 * (2.1667 * p_abs[i] * vmr[i] / t_abs[i]);
	  // inverse relative temperature [1]
	  Numeric ti   = (300.0 / t_abs[i]);
	  Numeric ti2  = pow(ti, 2.5);
	  
	  // continuum term [Np/km/GHz2]
	  Numeric con = pvap * pow(ti, 3.0) * 1.000e-9 * 
	                ( (0.543 * pda) + (17.96 * pvap * pow(ti, 4.5)) );
	  
	  // Loop over input frequency
	  for ( size_t s=0; s<n_f; ++s )
	    {
	      // input frequency in [GHz]
	      Numeric ff  = f_mono[s] * 1.00e-9;
	      // line contribution at position f
	      Numeric sum = 0.000;

	      // Loop over spectral lines
	      for (size_t l = 0; l < 15; l++) 
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
	      Numeric absl = 0.3183e-4 * den * sum;
	      // xsec = abs/vmr [1/m] (Rosenkranz model in [Np/km])
	      // 4.1907e-5 = 0.230259 * 0.1820 * 1.0e-3    (1/(10*log(e)) = 0.230259)
	      xsec[s][i]  += 1.000e-3 * ( absl + (CC * con * ff * ff) ) / vmr[i];	  
	    }
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
void PWR93O2AbsModel( MATRIX&           xsec,
		      Numeric             CC, // continuum term scale factor
		      Numeric             CL, // line strength scale factor
		      Numeric             CW, // line broadening scale factor
		      Numeric             CO, // line coupling scale factor
		      const VECTOR&   f_mono,
		      const VECTOR&    p_abs,
		      const VECTOR&    t_abs,
		      const VECTOR&   vmrh2o,
                      const VECTOR&      vmr )
{
  const size_t n_lines = 40; // number of O2 lines in this model (range: 50-850 GHz)

  // lines are arranged 1-,1+,3-,3+,etc. in spin-rotation spectrum
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

  // widths in MHz/mbar
  const Numeric WB300 = 0.56;
  const Numeric X     = 0.80;

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

  const size_t n_p = p_abs.size();	// Number of pressure levels
  const size_t n_f = f_mono.size();	// Number of frequencies

  // Check that dimensions of p_abs, t_abs, and vmr agree:
  assert ( n_p==t_abs.size() );
  assert ( n_p==vmr.size()   );

  // Check that dimensions of xsec are consistent with n_f
  // and n_p. It should be [n_f,n_p]:
  assert ( n_f==xsec.nrows() );
  assert ( n_p==xsec.ncols() );
  
  // Loop pressure/temperature:
  for ( size_t i=0; i<n_p; ++i )
    {
      // calculate xsec only if VMR(O2) > VMRCalcLimit
      if (vmr[i] > VMRCalcLimit)
	{
	  // relative inverse temperature [1]
	  Numeric TH  = 300.000 / t_abs[i];
	  Numeric TH1 = TH-1.000;
	  Numeric B   = pow(TH, X);
	  //cout << "TH,TH1,B= " << TH << ", " << TH1 << ", " << B << "\n";
	  // partial pressure of H2O and dry air [hPa]
	  Numeric PRESWV = 0.01 * p_abs[i] * vmrh2o[i];
	  Numeric PRESDA = 0.01 * p_abs[i] * (1.000 - vmrh2o[i]);
	  //cout << "PRESWV,PRESDA=" << PRESWV << ", " << PRESDA << "\n";
	  Numeric DEN    = .001*(PRESDA*B + 1.1*PRESWV*TH);
	  Numeric DFNR   = WB300*DEN;
	  // initial O2 line absorption at position ff 
	  Numeric O2ABS = 0.000;

	    // Loop over input frequency
	    for ( size_t s=0; s<n_f; ++s )
	      {
		// input frequency in [GHz]
		Numeric ff  = f_mono[s] * 1.00e-9; 
		// continuum absorption [Neper/km]
		Numeric SUM = CC * 1.6e-17 * ff * ff * DFNR / ( TH*(ff*ff + DFNR*DFNR) );
		//cout << "DEN,DFNR,SUM=" << DEN  << ", " << DFNR  << ", " << SUM << "\n";

		//cout << "ff =" << ff << ",  O2ABS=" << O2ABS << ",  SUM=" << SUM << "\n";
		// Loop over Rosnekranz '93 spectral line frequency:
		for ( size_t l=0; l<n_lines; ++l )
		  {
		    Numeric DF   = CO * W300[l] * DEN;
		    Numeric Y    = CW * 0.001 * 0.01 * p_abs[i] * B * ( Y300[l] + V[l]*TH1 );
		    Numeric STR  = CL * S300[l] * exp(-BE[l] * TH1);
		    Numeric SF1  = ( DF + (ff-F[l])*Y ) / ( (ff-F[l])*(ff-F[l]) + DF*DF );
		    Numeric SF2  = ( DF - (ff+F[l])*Y ) / ( (ff+F[l])*(ff+F[l]) + DF*DF );
		            SUM += STR * (SF1+SF2) * pow((ff/F[l]), 2.0);
		    //cout << "k,DF,Y,STR,SF1,SF2=" << l << ", " << DF << ", " << Y << ", " << STR
		    //<< ", " << SF1  << ", " << SF2 << "\n";
		  }
		// O2 absorption [Neper/km] 
		O2ABS += 0.5034e12 * SUM * PRESDA * pow(TH, 3.0) / 3.14159;
		//cout << "O2ABS=" << O2ABS << "\n";
		// unit conversion x Nepers/km = y 1/m  --->  y = x * 1.000e-3 
		// xsec [1/m]  1.000e-3
		xsec[s][i] += 1.000e-3 * O2ABS / vmr[i];
	      }
	}
    }
  return;
}
//
// #################################################################################
// #################### CONTROL OF ADDITIONAL ABSORPTION MODEL #####################
// #################################################################################
//
/** Calculates continuum absorption for one continuum tag. Note, that
    only one tag can be taken at a time. That means for water vapor
    you will have to call this function two times, once with the
    self-continuum tag and once with the foreign continuum tag.

    Calculated is the absorption cross section, that means you have to
    multiply this with the VMR in order to get the absorption
    coefficient:

    alpaha [1/m] = xsec * VMR

    \retval xsec       Cross section of one continuum tag.
    \param  name       The name of the continuum to calculate.
    \param  parameters Continuum model parameters, as defined in
                       cont_description_parameters. 
    \param  f_mono     Frequency grid.
    \param  p_abs      Pressure grid.
    \param  t_abs      Temperatures associated with p_abs.
    \param  vmrh2o     Total volume mixing ratio of water vapor. This
                       will be needed only for the oxygen continuum
    \param  vmrn2      Total volume mixing ratio of molecular nitrogen. This
                       will be needed only for the CO2 foreign continuum
    \param  vmr        Volume mixing ratio of the calculated species.

    \date   2001-01-16
    \author Stefan Buehler */
void xsec_continuum_tag( MATRIX&                    xsec,
			 const string&              name,
			 const VECTOR&              parameters,
			 const VECTOR&  	    f_mono,
			 const VECTOR&  	    p_abs,
			 const VECTOR&  	    t_abs,
			 const VECTOR&  	    n2_abs,
			 const VECTOR&  	    h2o_abs,
			 const VECTOR&              vmr )
{
  //
  // out3 << "  Continuum paramters are: " << parameters << "\n";
  //
  // A simple switch statement does not work here, because the
  // switching condition is not a simple value. So we need to use a
  // chain of if-else statements.
  //
  // ============= H2O continuum ================================================
  if ( "H2O-ContStandardSelf"==name )
    {
      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum coefficient (C_s)
      // parameters[1] : temperature exponent  (x_s)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1/m / (Hz^2*Pa^2)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      Rosenkranz_h2o_self_continuum( xsec,
				     parameters[0],
				     parameters[1],
				     f_mono,
				     p_abs,
				     t_abs,
				     vmr );
    }
  else if ( "H2O-ContStandardForeign"==name ) // ------------------------------
    {
      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum coefficient (C_f)
      // parameters[1] : temperature exponent  (x_f)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1/m / (Hz^2*Pa^2)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      Rosenkranz_h2o_foreign_continuum( xsec,
					parameters[0],
					parameters[1],
					f_mono,
					p_abs,
					t_abs,
					vmr );
    }
  else if ( "H2O-ContMPM93"==name ) // --------------------------------------
    {
      // only self and foreign continuum term is only simultaneously to calculated
      // since the parameterization is not devided up in these two terms.

      // Check if the right number of paramters has been specified:
      if ( 7 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : pseudo continuum line frequency                      default: 1780.0*10^9 Hz
      // parameters[1] : pseudo continuum line strength parameter             default: 22300.0 Hz/Pa
      // parameters[2] : pseudo continuum line strength temperature parameter default: 0.952
      // parameters[3] : pseudo continuum line broadening parameter           default: 17.6*10^4 Hz/Pa
      // parameters[4] : pseudo continuum line broadening parameter           default: 30.5
      // parameters[5] : pseudo continuum line broadening parameter           default: 2
      // parameters[6] : pseudo continuum line broadening parameter           default: 5
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [Hz]
      //     parameters[1] : [Hz/Pa]
      //     parameters[2] : [1]
      //     parameters[3] : [Hz/Pa]
      //     parameters[4] : [1]
      //     parameters[5] : [1]
      //     parameters[6] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      MPM93_h2o_continuum( xsec,
			   parameters[0],
			   parameters[1],
			   parameters[2],
			   parameters[3],
			   parameters[4],
			   parameters[5],
			   parameters[6],
			   f_mono,
			   p_abs,
			   t_abs,
			   vmr	 );
    }
  else if ( "H2O-ContCKDSelf"==name ) // ----------------------------------------
    {
	  ostringstream os;
	  os << "CKD self continuum model not yet implemented"
	     << ".";
	  throw runtime_error(os.str());
	  return;

    }
  else if ( "H2O-ContCKDForeign"==name ) // -------------------------------------
    {
	  ostringstream os;
	  os << "CKD foreign continuum model not yet implemented"
	     << ".";
	  throw runtime_error(os.str());
	  return;

    }
  else if ( "H2O-CP98Model"==name ) // ------------------------------
    {
      // Check if the right number of paramters has been specified:
      if ( 3 != parameters.size() )
	{
	  ostringstream os;
	  os << "Cruz-Pol absorption model " << name << " requires three input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum scale factor       CC (default = 1.2369 +/- 0.155)
      // parameters[1] : line strength scale factor   CL (default = 1.0639 +/- 0.016)
      // parameters[2] : line broadening scale factor CW (default = 1.0658 +/- 0.0096)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      CP98H2OAbsModel( xsec,
		       parameters[0],
		       parameters[1],
		       parameters[2],
		       f_mono,
		       p_abs,
		       t_abs,
		       vmr );
    }
  else if ( "H2O-MPM87Model"==name ) // ------------------------------
    {
      // Check if the right number of paramters has been specified:
      if ( 3 != parameters.size() )
	{
	  ostringstream os;
	  os << "MPM87 absorption model " << name << " requires three input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum scale factor       CC (default = 1.000)
      // parameters[1] : line strength scale factor   CL (default = 1.000)
      // parameters[2] : line broadening scale factor CW (default = 1.000)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      MPM87H2OAbsModel( xsec,
			parameters[0],
			parameters[1],
			parameters[2],
			f_mono,
			p_abs,
			t_abs,
			vmr );
    }
  else if ( "H2O-MPM89Model"==name ) // ------------------------------
    {
      // Check if the right number of paramters has been specified:
      if ( 3 != parameters.size() )
	{
	  ostringstream os;
	  os << "MPM89 absorption model " << name << " requires three input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum scale factor       CC (default = 1.000)
      // parameters[1] : line strength scale factor   CL (default = 1.000)
      // parameters[2] : line broadening scale factor CW (default = 1.000)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      MPM89H2OAbsModel( xsec,
			parameters[0],
			parameters[1],
			parameters[2],
			f_mono,
			p_abs,
			t_abs,
			vmr );
    }
  else if ( "H2O-MPM93Model"==name ) // ------------------------------
    {
      // Check if the right number of paramters has been specified:
      if ( 3 != parameters.size() )
	{
	  ostringstream os;
	  os << "MPM93 absorption model " << name << " requires three input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum scale factor       CC (default = 1.000)
      // parameters[1] : line strength scale factor   CL (default = 1.000)
      // parameters[2] : line broadening scale factor CW (default = 1.000)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      MPM93H2OAbsModel( xsec,
			parameters[0],
			parameters[1],
			parameters[2],
			f_mono,
			p_abs,
			t_abs,
			vmr );
    }
  else if ( "H2O-PWR98Model"==name ) // ------------------------------
    {
      // Check if the right number of paramters has been specified:
      if ( 3 != parameters.size() )
	{
	  ostringstream os;
	  os << "Rosenkranz98 absorption model " << name << " requires three input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum scale factor       CC (default = 1.000)
      // parameters[1] : line strength scale factor   CL (default = 1.000)
      // parameters[2] : line broadening scale factor CW (default = 1.000)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[2] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      PWR98H2OAbsModel( xsec,
			parameters[0],
			parameters[1],
			parameters[2],
			f_mono,
			p_abs,
			t_abs,
			vmr );
    }
  // ============= O2 continuum =================================================
  else if ( "O2-ContRosenkranz"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3
      // (see also JQSRT, Vol.48, No.5/6 pp.629-643, 1992)

      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum coefficient (C)
      // parameters[1] : temperature exponent  (x)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [K^2/(Hz*Pa*m)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      Rosenkranz_o2_continuum( xsec,
			       parameters[0], // coefficient
			       parameters[1], // temp. exponent
			       f_mono,
			       p_abs,
			       t_abs,
			       h2o_abs,
			       vmr );
    }
  else if ( "O2-PWR93O2Model"==name )
    {
      //  REFERENCE FOR EQUATIONS AND COEFFICIENTS:
      //  P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
      //  BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
      //  AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
      //  (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)

      // Check if the right number of paramters has been specified:
      if ( 4 != parameters.size() )
	{
	  ostringstream os;
	  os << "Rosenkranz O2 abs. model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum term scale factor,   default CC = 1.000
      // parameters[1] : line strength scale factor,    default CL = 1.000
      // parameters[2] : line broadening scale factor,  default CW = 1.000
      // parameters[3] : line coupling scale factor,    default CO = 1.000
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1]
      //     parameters[1] : [1]
      //     parameters[1] : [1]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     h2o_abs,      : [1]
      //     vmr           : [1]
      //
      PWR93O2AbsModel( xsec,
		       parameters[0], // continuum term scale factor
		       parameters[1], // line strength scale factor
		       parameters[2], // line broadening scale factor
		       parameters[3], // line coupling scale factor
		       f_mono,
		       p_abs,
		       t_abs,
		       h2o_abs,
		       vmr );
    }
  // ============= N2 continuum =================================================
  else if ( "N2-ContRosenkranzSelf"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3

      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum coefficient (C)
      // parameters[1] : temperature exponent  (x)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1/(Hz^2*Pa^2*m)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      Rosenkranz_n2_self_continuum( xsec,
				    parameters[0], // coefficient
				    parameters[1], // temp. exponent
				    f_mono,
				    p_abs,
				    t_abs,
				    vmr );
    }  
  else if ( "N2-ContBorysowSelf"==name ) // -------------------------------------
    {
      // data information about this continuum: 
      // A. Borysow and L. Frommhold, The Astrophysical Journal,
      // Vol. 311, pp.1043-1057, 1986
      ostringstream os;
      os << "N2 continuum parameterization of A. Borysow and L. Frommhold\n"
         << "is not yet implemented  ==>  no calculation performed!\n";
      /*
      Borysow_Frommhold_n2_continuum( xsec,
                                      parameters[0],
				      parameters[1],
				      f_mono,
				      p_abs,
				      t_abs,
				      vmr );
      */
      return;
    }  
  // ============= CO2 continuum ================================================
  else if ( "CO2-ContRosenkranzSelf"==name )
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3

      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum coefficient (C)
      // parameters[1] : temperature exponent  (x)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1/(Hz^2*Pa^2*m)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     vmr           : [1]
      //
      Rosenkranz_co2_self_continuum( xsec,
				     parameters[0], // coefficient
				     parameters[1], // temp. exponent
				     f_mono,
				     p_abs,
				     t_abs,
				     vmr );
    }
  else if ( "CO2-ContRosenkranzForeign"==name ) // ------------------------------
    {
      // data information about this continuum: 
      // P. W. Rosenkranz Chapter 2, pp 74, in M. A. Janssen, 
      // "Atmospheric Remote Sensing by Microwave Radiometry",
      // John Wiley & Sons, Inc., 1993, ISBN 0-471-62891-3

      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "Continuum model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // specific continuum parameters:
      // parameters[0] : continuum coefficient (C)
      // parameters[1] : temperature exponent  (x)
      //
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [1/(Hz^2*Pa^2*m)]
      //     parameters[1] : [1]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     n2_abs        : [1]
      //     vmr           : [1]
      //
      Rosenkranz_co2_foreign_continuum( xsec,
					parameters[0], // coefficient
					parameters[1], // temp. exponent
					f_mono,
					p_abs,
					t_abs,
					n2_abs,
					vmr );
    }
  // ============= cloud and fog absorption from MPM93 ==========================
  else if ( "H2O-MPM93droplet"==name )
    {
      // Suspended water droplet absorption parameterization from MPM93 model
      // H. J. Liebe and G. A. Hufford and M. G. Cotton,
      // "Propagation modeling of moist air and suspended water/ice
      //  particles at frequencies below 1000 GHz",
      // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 

      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "MPM93 liquid water particle absorption model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // liquid water droplet parameters:
      // parameters[0] : suspended water droplet density   range: 0-10 g/m3
      // parameters[1] : specific droplet weight           value:    1 g/cm3
      //
      // valid atmospheric condition:
      // temperature      : 233 to 323 K
      // relative humidity:   1 to 100 %
      // 
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [g/m3]
      //     parameters[1] : [g/m3]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     n2_abs        : [1]
      //     vmr           : [1]
      //
      MPM93WaterDropletAbs(xsec,
			   parameters[0],     // suspended water droplet density
			   parameters[1],     // specific droplet weight 
			   f_mono,
			   p_abs,
			   t_abs,
			   vmr );
    }
  // ============= ice particle absorption from MPM93 ============================
  else if ( "H2O-MPM93ice"==name )
    {
      // Ice particle absorption parameterization from MPM93 model
      // H. J. Liebe and G. A. Hufford and M. G. Cotton,
      // "Propagation modeling of moist air and suspended water/ice
      //  particles at frequencies below 1000 GHz",
      // AGARD 52nd Specialists Meeting of the Electromagnetic Wave
      // Propagation Panel, Palma de Mallorca, Spain, 1993, May 17-21 

      // Check if the right number of paramters has been specified:
      if ( 2 != parameters.size() )
	{
	  ostringstream os;
	  os << "MPM93 ice particle absorption model " << name << " requires two input\n"
	     << "parameters, but you specified " << parameters.size()
	     << ".";
	  throw runtime_error(os.str());
	  return;
	}
      
      // ice crystal parameters:
      // parameters[0] : suspended ice particle density   range:  0-10 g/m3
      // parameters[1] : specific ice particle weight     value: 0.916 g/cm3
      //
      // valid atmospheric condition:
      // temperature      : 233 to 323 K
      // relative humidity:   1 to 100 %
      // 
      // units:
      //  a) output 
      //     xsec          : [1/m],
      //  b) input
      //     parameters[0] : [g/m3]
      //     parameters[1] : [g/cm3]
      //     f_mono        : [Hz]
      //     p_abs         : [Pa]
      //     t_abs         : [K]
      //     n2_abs        : [1]
      //     vmr           : [1]
      //
      MPM93IceCrystalAbs(xsec,
			 parameters[0],     // suspended water droplet density
			 parameters[1],     // specific droplet weight 
			 f_mono,
			 p_abs,
			 t_abs,
			 vmr );
    }
  else // -----------------------------------------------------------------------
    {
      ostringstream os;
      os << "Continuum tag `" << name << "' not implemented.";
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
void check_continuum_model(const string& name)
{
  // The species lookup data:
  extern const ARRAY<SpeciesRecord> species_data;

  // For the list of valid continuum models:
  ARRAY<string> valid_models;

  bool found = false;

  // Loop all species:
  for ( ARRAY<SpeciesRecord>::const_iterator i=species_data.begin();
	i<species_data.end();
	++i )
    {
      string specnam = i->Name();

      // Loop all isotopes:
      for ( ARRAY<IsotopeRecord>::const_iterator j=i->Isotope().begin();
	    j<i->Isotope().end();
	    ++j )
	{
	  string isonam = j->Name();

	  // The specified name consists of a species part and an
	  // isotope part, e.g., H2O-ContStandardSelf. We need to
	  // construct a similar string from the species lookup data
	  // by concatenating species name and isotope name.

	  string fullnam = specnam + "-" + isonam;
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
      os << "The string `" << name << "' matches none of the known\n"
	 << "continuum models. Known continuum models are:";
      for ( ARRAY<string>::const_iterator i=valid_models.begin();
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
