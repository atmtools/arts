/* Copyright (C) 2000-2012
   Axel von Engeln <engeln@uni-bremen.de>
   Stefan Buehler  <sbuehler@ltu.se>

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
  \file   lineshapes.cc
  \brief  Stuff related to lineshape functions.

  This file contains both the lineshape functions themselves and the
  function define_lineshape_data which sets the lineshape lookup
  data. 

  This is the file from arts-1-0, back-ported to arts-1-1.

  \author Stefan Buehler
  \date   2000-08-21
*/

#include <cmath>
#include "arts.h"
#include "matpackI.h"
#include "array.h"
#include "absorption.h"
#include "Faddeeva.hh"

/*! The dummy line shape. This lineshape does nothing. It only exists,
    because formally you have to specify a lineshape also for
    continuum tags. It has to have the same arguments as all the other
    lineshapes, though...

    \retval ls_attenuation              The shape function.
    \retval ls_phase                    The shape function.
    \retval X                           Auxillary parameter, only used in Voigt fct.
    \param  f0                          Line center frequency.
    \param  gamma                       The pressure broadening parameter.
    \param  sigma                       The Doppler broadening parameter. (Not used.)
    \param  f_grid                      The frequency grid.

    \throw  runtime_error This exception is always thrown when the
                          function is called.  

    \date   2001-01-16 
    \author Stefan Buehler 
*/
void lineshape_no_shape(  Vector&,
                          Vector&,
                          Vector&,
                          Vector&,
                          Vector&,
                          Vector&,
                          Vector&,
                          const Numeric,
                          const Numeric,
                          const Numeric,
                          const Numeric,
                          const Numeric,
                          const Numeric,
                          const Numeric,
                          const Numeric,
                          ConstVectorView,
                          const bool,
                          const bool)
{
  // This function should never be called so throw an error here: 
  throw std::runtime_error("The no_shape lineshape is only a placeholder, but you tried\n"
                      "to use it like a real lineshape.");
}


/*! The Lorentz line shape. This is a quick and dirty implementation.

    \retval ls_attenuation              The shape function.
    \retval ls_phase                    The shape function.
    \retval X                           Auxillary parameter, only used in Voigt fct.
    \param  f0                          Line center frequency.
    \param  gamma                       The pressure broadening parameter.
    \param  sigma                       The Doppler broadening parameter. (Not used.)
    \param  f_grid                      The frequency grid.

    \author Stefan Buehler 
    \date 2000-06-16 */
void lineshape_lorentz(Vector&         ls_attenuation,
                       Vector&         ls_phase,
                       Vector&,        //ls_dattenuation_dfrequency_term
                       Vector&,        //ls_dphase_dfrequency_term
                       Vector&,        //ls_dattenuation_dpressure_term
                       Vector&,        //ls_dphase_dpressure_term
                       Vector&,        //X
                       const Numeric   f0,
                       const Numeric   gamma,
                       const Numeric,
                       const Numeric,
                       const Numeric,
                       const Numeric,
                       const Numeric,  //sigma
                       const Numeric,
                       ConstVectorView f_grid,
                       const bool      do_phase,
                       const bool)
{
  const Index nf = f_grid.nelem();

  // PI:
  extern const Numeric PI;

  //  assert( ls.nelem() == nf );

  Numeric gamma2 = gamma * gamma;
  Numeric fac = gamma/PI;
  
  if(do_phase)
  {
    for ( Index i=0; i<nf; ++i )
    {
        ls_attenuation[i] =  fac / ( (f_grid[i]-f0) * (f_grid[i]-f0) + gamma2 );
        ls_phase[i] = (f_grid[i]-f0) / PI / ( (f_grid[i]-f0) * (f_grid[i]-f0) + gamma2 );
    }
  }
  else
  {
      for ( Index i=0; i<nf; ++i )
      {
          ls_attenuation[i] =  fac / ( (f_grid[i]-f0) * (f_grid[i]-f0) + gamma2 );
      }
  }     
}

/*! The Mirrored Lorentz line shape. This is a quick and dirty implementation.
 * 
 *   \retval ls_attenuation              The shape function.
 *   \retval ls_phase                    The shape function.
 *   \param  f0                          Line center frequency.
 *   \param  gamma                       The pressure broadening parameter.
 *   \param  f_grid                      The frequency grid.
 * 
 *   \author Richard Larsson
 *   \date 2015-08-13 */
void lineshape_mirrored_lorentz(Vector&         ls_attenuation,
                                Vector&         ls_phase,
                                Vector&,        //ls_dattenuation_dfrequency_term
                                Vector&,        //ls_dphase_dfrequency_term
                                Vector&,        //ls_dattenuation_dpressure_term
                                Vector&,        //ls_dphase_dpressure_term
                                Vector&,
                                const Numeric   f0,
                                const Numeric   gamma,
                                const Numeric,
                                const Numeric,
                                const Numeric,
                                const Numeric,
                                const Numeric,
                                const Numeric,
                                ConstVectorView f_grid,
                                const bool      do_phase,
                                const bool)
{
    const Index nf = f_grid.nelem();
    
    // PI:
    extern const Numeric PI;
    
    //  assert( ls.nelem() == nf );
    
    Numeric gamma2 = gamma * gamma;
    Numeric fac = gamma/PI;
    
    if(do_phase)
    {
        for ( Index i=0; i<nf; ++i )
        {
            ls_attenuation[i] =  fac / ( (f_grid[i]-f0) * (f_grid[i]-f0) + gamma2 ) + fac / ( (f_grid[i]+f0) * (f_grid[i]+f0) + gamma2 );
            ls_phase[i] = (f_grid[i]-f0) / PI / ( (f_grid[i]-f0) * (f_grid[i]-f0) + gamma2 ) + (f_grid[i]+f0) / PI / ( (f_grid[i]+f0) * (f_grid[i]+f0) + gamma2 );//uncertain on this part
        }
    }
    else 
    {
        for ( Index i=0; i<nf; ++i )
        {
            ls_attenuation[i] =  fac / ( (f_grid[i]-f0) * (f_grid[i]-f0) + gamma2 ) + fac / ( (f_grid[i]+f0) * (f_grid[i]+f0) + gamma2 );
        }
    }
}

/*! The Doppler line shape.

    \retval ls_attenuation              The shape function.
    \retval ls_phase                    The shape function.
    \retval x                           Auxillary parameter, only used in Voigt fct.
    \param  f0                          Line center frequency.
    \param  gamma                       The pressure broadening parameter. (Not used.)
    \param  sigma                       The Doppler broadening parameter.
    \param  f_grid                      The frequency grid.

    \author Axel von Engeln
    \date 2000-12-06 */
void lineshape_doppler(Vector&         ls_attenuation,
                       Vector&,        //ls_phase
                       Vector&,        //ls_dattenuation_dfrequency_term
                       Vector&,        //ls_dphase_dfrequency_term
                       Vector&,        //ls_dattenuation_dpressure_term
                       Vector&,        //ls_dphase_dpressure_term
                       Vector&,        //x
                       const Numeric   f0,
                       const Numeric,  //gamma
                       const Numeric,
                       const Numeric,
                       const Numeric,
                       const Numeric,
                       const Numeric   sigma,
                       const Numeric,
                       ConstVectorView f_grid,
                       const bool,
                       const bool)
{
  const Index nf = f_grid.nelem();

  // SQRT(PI):
  extern const Numeric PI;
  const Numeric sqrtPI = sqrt(PI);

  //  assert( ls_attenuation.nelem() == nf );

  Numeric sigma2 = sigma * sigma;
  Numeric fac = 1.0 / (sqrtPI * sigma);
  
  for ( Index i=0; i<nf ; ++i )
    {
      ls_attenuation[i] = fac * exp( - pow( f_grid[i]-f0, 2) / sigma2 );
    }
}


/*! The O2 non-resonant line shape.  Should be VVW/2 so do not use this call...
 * 
 *   \retval ls_attenuation              The shape function.
 *   \retval ls_phase                    The shape function.
 *   \retval X                           Auxillary parameter, only used in Voigt fct.
 *   \param  f0                          Line center frequency.
 *   \param  gamma                       The pressure broadening parameter.
 *   \param  sigma                       The Doppler broadening parameter. (Not used.)
 *   \param  f_grid                      The frequency grid.
 * 
 *   \author Richard Larsson
 *   \date 2014-10-02 */
void deprecated( Vector&, Vector&, Vector&, Vector&, Vector&, Vector&, Vector&,
                 const Numeric, const Numeric, const Numeric, const Numeric,
                 const Numeric, const Numeric, const Numeric, const Numeric,
                 ConstVectorView, const bool, const bool)
{
    throw std::runtime_error("You are running a deprecated lineshape.\n"
        "The deprecated line shapes are:\n"
        "\t Voigt_Kuntz3 --- please replace with Voigt_Kuntz6\n"
        "\t Voigt_Kuntz4 --- please replace with Voigt_Kuntz6\n"
    );
}
//------------------------------------------------------------------------

// help function for lineshape_voigt_kuntz1
long bfun6_(Numeric y, Numeric x)
{
  /* System generated locals */
  long int ret_val;

  /* Local variables */
  Numeric s = 0;

  /* -------------------------------------------------------------------- */
  s = fabs(x) + y;
  if (s >= 15.f) {
    ret_val = 1;
  } else if (s >= 5.5f) {
    ret_val = 2;
  } else if (y >= fabs(x) * .195f - .176f) {
    ret_val = 3;
  } else {
    ret_val = 4;
  }
  return ret_val;
} /* bfun6_ */



/*! The Voigt line shape. Kuntz approximation of the Voigt line
  shape. 

    \retval ls_attenuation            The shape function.
    \retval x             Auxillary parameter to store frequency grid.
    \param  f0            Line center frequency.
    \param  gamma         The pressure broadening parameter.
    \param  sigma         The Doppler broadening parameter.
    \param  f_grid        The frequency grid.

    Original c function call and documention:

    int voigt
    (
    long nx,
    float *x,
    float y,
    float *prb,
    float fak
    )
    
    Calculates the Voigt-Function times the user-definied value 
    fac with a relative accuracy better than 2*10-6.
    
    If this subroutine is called several times with the same 
    parameter y the numerically expensive coefficents a1..t8 
    are only calculated once thus further accelerating the 
    algorithm

    \verbatim
    --------------------------------------------------------------------
    x(nx)   (in)    :Distance from line center in units of Doppler
                    :halfwidths
    y       (in)    :Ratio of the Doppler halfwidth to the Lorentz
                    :halfwidth  
    prb     (out)   :voigt-function times fak
    fak     (in)    :factor to be specified by the user
    --------------------------------------------------------------------

    author: M. Kuntz, 
            Institut fuer Meteorologie und Klimaforschung, 
            Forschungszentrum Karlsruhe, 
            Postfach 3640, 
            76021 Karlsruhe, Germany. 
            email: kuntz@imk.fzk.de 

    \endverbatim


    About 'voigt' : The program was originally written by M. Kuntz in
    Fortran77 but has been translated into C by Dietrich Feist (f2c)
    and into C++ by Oliver Lemke and Axel von Engeln. fak is removed
    from program code. Replaced nx by nf. Replaced prb by ls.
    Multiplied ls with the factor fac.
    

    \author Oliver Lemke and Axel von Engeln
    \date 2000-09-27
*/ 
void lineshape_voigt_kuntz6(Vector&         ls_attenuation,
                            Vector&,        //ls_phase
                            Vector&         ls_dattenuation_dfrequency_term,
                            Vector&,        //ls_dphase_dfrequency_term
                            Vector&         ls_dattenuation_dpressure_term,
                            Vector&,        //ls_dphase_dpressure_term
                            Vector&         x,
                            const Numeric   f0,
                            const Numeric   gamma,
                            const Numeric,
                            const Numeric,
                            const Numeric,
                            const Numeric,
                            const Numeric   sigma,
                            const Numeric,
                            ConstVectorView f_grid,
                            const bool,
                            const bool do_partials)

{
  const Index nf = f_grid.nelem();

  // seems not necessary for Doppler correction
  //    extern const Numeric SQRT_NAT_LOG_2;

  // PI
  extern const Numeric PI;

  // constant sqrt(1/pi)
  const Numeric sqrt_invPI =  sqrt(1/PI);

  // constant normalization factor for voigt
  Numeric fac = 1.0 / sigma * sqrt_invPI;


  /* Initialized data */

  Numeric yps1 = -1.0;
  Numeric yps2 = -1.0;
  Numeric yps3 = -1.0;
  Numeric yps4 = -1.0;

  /* System generated locals */
  long int i__1, i__2;
  Numeric r__1;

  /* Local variables */
  long int bmin = 0, lauf[16] = {0}        /* was [4][4] */, bmax;
  long int imin = 0, imax = 0, stack[80] = {0} /* was [20][4] */;
  Numeric a1, a2, a3, a4, a5, a6, a8, b8, c8, d8, e8, f8, g8, h8, a7, 
    b7, c7, d7, e7, f7, o8, p8, q8, r8, s8, t8, g7, h7, o7, p7, q7, 
    r7, s7, t7, b6, c6, d6, e6, b5, c5, d5, e5, b4, c4, d4, b3, c3, 
    d3, b1, y2;
  a1 = a2 = a3 = a4 = a5 = a6 = a8 = b8 = c8 = d8 = e8 = f8 = g8 = h8 = a7
    = b7 = c7 = d7 = e7 = f7 = o8 = p8 = q8 = r8 = s8 = t8 = g7 = h7 = o7 = p7
    = q7 = r7 = s7 = t7 = b6 = c6 = d6 = e6 = b5 = c5 = d5 = e5 = b4 = c4 = d4
    = b3 = c3 = d3 = b1 = y2 = 0;
  long int i2 = 0, i1 = 0;
  Numeric x2 = 0, b2 = 0, c1 = 0;
  long int stackp = 0, imitte = 0;
  Numeric ym2 = 0;


  // variables needed in original c routine:

  // Ratio of the Lorentz halfwidth to the Doppler halfwidth
  Numeric y = gamma / sigma;

  // Helper bools
  bool do_x=false, do_y=false;
  
  // Forcing 0.0001 here, or 1/10000th of relevant broadening.
  // This is from Wfuns agreeing on less than 1e-6 with this
  // level of perturbation (faddeeva_algorithm_916 and Voigt_Kuntz6 were compared)
  const Numeric dx=0.0001, dy=0.0001;
  
  // do_partials is put here because later on x will change...
  if(do_partials)
  {
      // pass variable for self-calling
      Vector empty;
      
      // x = (f-f0)/sigma, so x+dx from f0+df0 means dx = - df0/sigma
      const Numeric df0 = - dx*sigma; 
      
      // Use this dx to see if reasonable to calculate dFa/dx.  
      // Too small dx and we will have constant
      // Jacobian that are unreasonably large...
      do_x = (f0>(-df0*10));
      
      if(do_x)
      {
          lineshape_voigt_kuntz6(ls_dattenuation_dfrequency_term,
                                 empty,empty,empty,empty,empty,
                                 x,
                                 f0+df0,
                                 gamma,
                                 0.0,0.0,0.0,0.0,
                                 sigma,
                                 0.0,
                                 f_grid,
                                 false,false);
      }
      
      // y = gamma/sigma so y+dy from gamma+dgamma means dy=dgamma/sigma
      const Numeric dgamma = dy*sigma; 
      
      // Use this dy to see if reasonable to calculate dFa/dy.  
      // Too small dy and we will have constant
      // Jacobian that are unreasonably large
      do_y = (y>(dy*10));
      
      if(do_y)
      {
          lineshape_voigt_kuntz6(ls_dattenuation_dpressure_term,
                                 empty,empty,empty,empty,empty,
                                 x,
                                 f0,
                                 gamma+dgamma,
                                 0.0,0.0,0.0,0.0,
                                 sigma,
                                 0.0,
                                 f_grid,
                                 false,false);
      }
  }
  
  
  // frequency in units of Doppler 
  for (i1=0; i1< (int) nf; i1++)
    {
      x[i1] = (f_grid[i1] - f0) / sigma;
    }

  /* Parameter adjustments */
  // this does not work for variables type Vector, corrected by
  // adjusting the index to array ls and x
  //--ls;
  //--x;

  /* Function Body */
  y2 = y * y;
  if (y >= 15.0 || x[0] >= 15.0 || x[nf-1] <= -15.0) {
    lauf[0] = 1;
    lauf[4] = nf;
    lauf[8] = nf;
    lauf[12] = 0;
    goto L7;
  }
  for (i2 = 1; i2 <= 4; ++i2) {
    for (i1 = 1; i1 <= 4; ++i1) {
      lauf[i1 + (i2 << 2) - 5] = i2 % 2 * (nf + 1);
      /* L1: */
    }
  }
  stackp = 1;
  stack[stackp - 1] = 1;
  stack[stackp + 19] = nf;
  stack[stackp + 39] = bfun6_(y, x[0]);
  stack[stackp + 59] = bfun6_(y, x[nf-1]);
 L2:
  imin = stack[stackp - 1];
  imax = stack[stackp + 19];
  bmin = stack[stackp + 39];
  bmax = stack[stackp + 59];
  if (bmin == bmax) {
    if (x[imax-1] < 0.f) {
      /* Computing MIN */
      i__1 = imin, i__2 = lauf[bmin - 1];
      lauf[bmin - 1] = min(i__1,i__2);
      /* Computing MAX */
      i__1 = imax, i__2 = lauf[bmax + 3];
      lauf[bmax + 3] = max(i__1,i__2);
      --stackp;
      goto L3;
    } else if (x[imin-1] >= 0.f) {
      /* Computing MIN */
      i__1 = imin, i__2 = lauf[bmin + 7];
      lauf[bmin + 7] = min(i__1,i__2);
      /* Computing MAX */
      i__1 = imax, i__2 = lauf[bmax + 11];
      lauf[bmax + 11] = max(i__1,i__2);
      --stackp;
      goto L3;
    }
  }
  imitte = (imax + imin) / 2;
  stack[stackp - 1] = imitte + 1;
  stack[stackp + 19] = imax;
  stack[stackp + 39] = bfun6_(y, x[imitte]);
  stack[stackp + 59] = bmax;
  ++stackp;
  stack[stackp - 1] = imin;
  stack[stackp + 19] = imitte;
  stack[stackp + 39] = bmin;
  stack[stackp + 59] = bfun6_(y, x[imitte-1]);
 L3:
  if (stackp > 0) {
    goto L2;
  }
  /* ---- Region 4 */
  /* -------------------------------------------------------------------- */
  if (lauf[7] >= lauf[3] || lauf[15] >= lauf[11]) {
    if ((r__1 = y - yps4, fabs(r__1)) > 1e-8f) {
      //yps4 = y;
      a7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * (y2 * (y2 * (2.35944f - y2 * .56419f) - 72.9359f) + 
                    571.687f) - 5860.68f) + 40649.2f) - 320772.f) + 1684100.f)
                     - 9694630.f) + 40816800.f) - 1.53575e8f) + 4.56662e8f) - 
                    9.86604e8f) + 1.16028e9f);
      b7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * (y2 * (23.0312f - y2 * 7.33447f) - 234.143f) - 
                    2269.19f) + 45251.3f) - 234417.f) + 3599150.f) - 
                    7723590.f) + 86482900.f) - 2.91876e8f) + 8.06985e8f) - 
                    9.85386e8f) - 5.60505e8f);
      c7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * (97.6203f - y2 * 44.0068f) + 1097.77f) - 25338.3f) + 
                    98079.1f) + 576054.f) - 2.3818e7f) + 22930200.f) - 
                    2.04467e8f) + 2.94262e8f) + 2.47157e8f) - 6.51523e8f);
      d7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    228.563f - y2 * 161.358f) + 8381.97f) - 66431.2f) - 
                    303569.f) + 2240400.f) + 38311200.f) - 41501300.f) - 
                    99622400.f) + 2.70167e8f) - 2.63894e8f);
      e7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    296.38f - y2 * 403.396f) + 23507.6f) - 66212.1f) - 
                    1.003e6f) + 468142.f) + 24620100.f) + 5569650.f) + 
                    1.40677e8f) - 63177100.f);
      f7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (125.591f - 
                    y2 * 726.113f) + 37544.8f) + 8820.94f) - 934717.f) - 
                    1931140.f) - 33289600.f) + 4073820.f) - 16984600.f);
      g7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (-260.198f - y2 * 
                    968.15f) + 37371.9f) + 79902.5f) - 186682.f) - 900010.f) 
                    + 7528830.f) - 1231650.f);
      h7 = y * (y2 * (y2 * (y2 * (y2 * (y2 * (-571.645f - y2 * 968.15f)
                     + 23137.1f) + 72520.9f) + 153468.f) + 86407.6f) - 
                    610622.f);
      o7 = y * (y2 * (y2 * (y2 * (y2 * (-575.164f - y2 * 726.113f) + 
                    8073.15f) + 26538.5f) + 49883.8f) - 23586.5f);
      p7 = y * (y2 * (y2 * (y2 * (-352.467f - y2 * 403.396f) + 
                    953.655f) + 2198.86f) - 8009.1f);
      q7 = y * (y2 * (y2 * (-134.792f - y2 * 161.358f) - 271.202f) - 
                    622.056f);
      r7 = y * (y2 * (-29.7896f - y2 * 44.0068f) - 77.0535f);
      s7 = y * (-2.92264f - y2 * 7.33447f);
      t7 = y * -.56419f;
      a8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * (y2 * (y2 * (y2 - 3.68288f) + 126.532f) - 955.194f) 
                    + 9504.65f) - 70946.1f) + 483737.f) - 2857210.f) + 
                    14464700.f) - 61114800.f) + 2.11107e8f) - 5.79099e8f) + 
                    1.17022e9f) - 1.5599e9f) + 1.02827e9f;
      b8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * (y2 * (y2 * 14.f - 40.5117f) + 533.254f) + 3058.26f) 
                    - 55600.f) + 498334.f) - 2849540.f) + 13946500.f) - 
                    70135800.f) + 2.89676e8f) - 7.53828e8f) + 1.66421e9f) - 
                    2.28855e9f) + 1.5599e9f;
      c8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * (y2 * 91 - 198.876f) - 1500.17f) + 48153.3f) - 
                    217801.f) - 1063520.f) + 1.4841e7f) - 46039600.f) + 
                    63349600.f) - 6.60078e8f) + 1.06002e9f) - 1.66421e9f) + 
                    1.17022e9f;
      d8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (
                    y2 * 364 - 567.164f) - 16493.7f) + 161461.f) + 280428.f) 
                    - 6890020.f) - 6876560.f) + 1.99846e8f) + 54036700.f) + 
                    6.60078e8f) - 7.53828e8f) + 5.79099e8f;
      e8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * 
                    1001 - 1012.79f) - 55582.f) + 240373.f) + 1954700.f) - 
                    5257220.f) - 50101700.f) - 1.99846e8f) + 63349600.f) - 
                    2.89676e8f) + 2.11107e8f;
      f8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * 2002 - 
                    1093.82f) - 106663.f) + 123052.f) + 3043160.f) + 
                    5257220.f) - 6876560.f) + 46039600.f) - 70135800.f) + 
                    61114800.f;
      g8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * 3003 - 
                    486.14f) - 131337.f) - 123052.f) + 1954700.f) + 6890020.f)
                     + 1.4841e7f) - 13946500.f) + 14464700.f;
      h8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * (y2 * 3432 + 486.14f) - 
                    106663.f) - 240373.f) + 280428.f) + 1063520.f) - 
                    2849540.f) + 2857210.f;
      o8 = y2 * (y2 * (y2 * (y2 * (y2 * (y2 * 3003 + 1093.82f) - 
                    55582.f) - 161461.f) - 217801.f) - 498334.f) + 483737.f;
      p8 = y2 * (y2 * (y2 * (y2 * (y2 * 2002 + 1012.79f) - 16493.7f) - 
                    48153.3f) - 55600.f) + 70946.1f;
      q8 = y2 * (y2 * (y2 * (y2 * 1001.f + 567.164f) - 1500.17f) - 
                    3058.26f) + 9504.65f;
      r8 = y2 * (y2 * (y2 * 364 + 198.876f) + 533.254f) + 955.194f;
      s8 = y2 * (y2 * 91.f + 40.5117f) + 126.532f;
      t8 = y2 * 14.f + 3.68288f;
    }
    ym2 = y * 2;
    for (i2 = 1; i2 <= 3; i2 += 2) {
      i__1 = lauf[((i2 + 1) << 2) - 1];
      for (i1 = lauf[(i2 << 2) - 1]; i1 <= i__1; ++i1) {
        x2 = x[i1-1] * x[i1-1];
        ls_attenuation[i1-1] = fac * (exp(y2 - x2) * cos(x[i1-1] * ym2) - (a7 + x2 *
                         (b7 + x2 * (c7 + x2 * (d7 + x2 * (e7 + x2 * (f7 + x2 
                        * (g7 + x2 * (h7 + x2 * (o7 + x2 * (p7 + x2 * (q7 + 
                        x2 * (r7 + x2 * (s7 + x2 * t7))))))))))))) / (a8 + x2 
                        * (b8 + x2 * (c8 + x2 * (d8 + x2 * (e8 + x2 * (f8 + 
                        x2 * (g8 + x2 * (h8 + x2 * (o8 + x2 * (p8 + x2 * (q8 
                        + x2 * (r8 + x2 * (s8 + x2 * (t8 + x2)))))))))))))));
        /* L4: */
      }
    }
  }
  /* ---- Region 3 */
  /* -------------------------------------------------------------------- */
  if (lauf[6] >= lauf[2] || lauf[14] >= lauf[10]) {
    if ((r__1 = y - yps3, fabs(r__1)) > 1e-8f) {
      //yps3 = y;
      a5 = y * (y * (y * (y * (y * (y * (y * (y * (y * 
                    .564224f + 7.55895f) + 49.5213f) + 204.501f) + 581.746f) 
                    + 1174.8f) + 1678.33f) + 1629.76f) + 973.778f) + 272.102f;
      b5 = y * (y * (y * (y * (y * (y * (y * 2.25689f + 22.6778f)
                     + 100.705f) + 247.198f) + 336.364f) + 220.843f) - 
                    2.34403f) - 60.5644f;
      c5 = y * (y * (y * (y * (y * 3.38534f + 22.6798f) + 52.8454f)
                     + 42.5683f) + 18.546f) + 4.58029f;
      d5 = y * (y * (y * 2.25689f + 7.56186f) + 1.66203f) - .128922f;
      e5 = y * .564224f + 9.71457e-4f;
      a6 = y * (y * (y * (y * (y * (y * (y * (y * (y * (y + 
                    13.3988f) + 88.2674f) + 369.199f) + 1074.41f) + 2256.98f) 
                    + 3447.63f) + 3764.97f) + 2802.87f) + 1280.83f) + 
                    272.102f;
      b6 = y * (y * (y * (y * (y * (y * (y * (y * 5.f + 
                    53.5952f) + 266.299f) + 793.427f) + 1549.68f) + 2037.31f) 
                    + 1758.34f) + 902.306f) + 211.678f;
      c6 = y * (y * (y * (y * (y * (y * 10.f + 80.3928f) + 
                    269.292f) + 479.258f) + 497.302f) + 308.186f) + 78.866f;
      d6 = y * (y * (y * (y * 10.f + 53.5952f) + 92.7568f) + 
                    55.0293f) + 22.0353f;
      e6 = y * (y * 5.f + 13.3988f) + 1.49645f;
    }
    for (i2 = 1; i2 <= 3; i2 += 2) {
      i__1 = lauf[((i2 + 1) << 2) - 2];
      for (i1 = lauf[(i2 << 2) - 2]; i1 <= i__1; ++i1) {
        x2 = x[i1-1] * x[i1-1];
        ls_attenuation[i1-1] = fac * (a5 + x2 * (b5 + x2 * (c5 + x2 * (d5 + x2 * 
                        e5)))) / (a6 + x2 * (b6 + x2 * (c6 + x2 * (d6 + x2 * (
                        e6 + x2)))));
        /* L5: */
      }
    }
  }
  /* ---- Region 2 */
  /* -------------------------------------------------------------------- */
  if (lauf[5] >= lauf[1] || lauf[13] >= lauf[9]) {
    if ((r__1 = y - yps2, fabs(r__1)) > 1e-8f) {
      //yps2 = y;
      a3 = y * (y2 * (y2 * (y2 * .56419f + 3.10304f) + 4.65456f) + 
                    1.05786f);
      b3 = y * (y2 * (y2 * 1.69257f + .56419f) + 2.962f);
      c3 = y * (y2 * 1.69257f - 2.53885f);
      d3 = y * .56419f;
      a4 = y2 * (y2 * (y2 * (y2 + 6.f) + 10.5f) + 4.5f) + .5625f;
      b4 = y2 * (y2 * (y2 * 4.f + 6.f) + 9.f) - 4.5f;
      c4 = y2 * (y2 * 6.f - 6.f) + 10.5f;
      d4 = y2 * 4.f - 6.f;
    }
    for (i2 = 1; i2 <= 3; i2 += 2) {
      i__1 = lauf[((i2 + 1) << 2) - 3];
      for (i1 = lauf[(i2 << 2) - 3]; i1 <= i__1; ++i1) {
        x2 = x[i1-1] * x[i1-1];
        ls_attenuation[i1-1] = fac * (a3 + x2 * (b3 + x2 * (c3 + x2 * d3))) / (a4 
                        + x2 * (b4 + x2 * (c4 + x2 * (d4 + x2))));
        /* L6: */
      }
    }
  }
  /* ---- Region 1 */
  /* -------------------------------------------------------------------- */
 L7:
  if (lauf[4] >= lauf[0] || lauf[12] >= lauf[8]) {
    if ((r__1 = y - yps1, fabs(r__1)) > 1e-8f) {
      //yps1 = y;
      a1 = y * .5641896f;
      b1 = y2 + .5f;
      a2 = y2 * 4;
    }

    c1 = fac * a1;
    for (i2 = 1; i2 <= 3; i2 += 2) {
      i__1 = lauf[((i2 + 1) << 2) - 4];
      for (i1 = lauf[(i2 << 2) - 4]; i1 <= i__1; ++i1) {
        x2 = x[i1-1] * x[i1-1];
        b2 = b1 - x2;
        ls_attenuation[i1-1] = c1 * (b1 + x2) / (b2 * b2 + a2 * x2);
        /* L8: */
      }
    }
  }
  
  if(do_partials)
  {
    // Set things right depending on if possible
    if(do_y&&do_x)
    {
        for(Index iv=0; iv<nf; iv++)
        {
            ls_dattenuation_dfrequency_term[iv] -= ls_attenuation[iv];
            ls_dattenuation_dfrequency_term[iv] /= dx;
            ls_dattenuation_dpressure_term[iv]  -=  ls_attenuation[iv];
            ls_dattenuation_dpressure_term[iv]  /= dy;
        }
    }
    else if(do_x)
    {
        for(Index iv=0; iv<nf; iv++)
        {
            ls_dattenuation_dfrequency_term[iv] -= ls_attenuation[iv];
            ls_dattenuation_dfrequency_term[iv] /= dx;
            ls_dattenuation_dpressure_term[iv] = 0.0;
        }
    }
    else if(do_y)
    {
        for(Index iv=0; iv<nf; iv++)
        {
            ls_dattenuation_dfrequency_term[iv] = 0.0;
            ls_dattenuation_dpressure_term[iv]  -=  ls_attenuation[iv];
            ls_dattenuation_dpressure_term[iv]  /= dy;
        }
    }
    else
    {
        for(Index iv=0; iv<nf; iv++)
        {
            ls_dattenuation_dfrequency_term[iv] = 0.0;
            ls_dattenuation_dpressure_term[iv] = 0.0;
        }
    }
  }
}


//---------------------------------------------------------------
// help function for lineshape_voigt_kuntz3

long int bfun3_(Numeric y, Numeric x)
{
    /* System generated locals */
    long int ret_val;

    /* Local variables */
    Numeric x2 = 0, y2 = 0;

/* -------------------------------------------------------------------- */
    x2 = x * x;
    y2 = y * y;
    if (x2 * .4081676f + y2 > 21.159543f) {
        if (x2 * .7019639f + y2 > 1123.14221f) {
            ret_val = 0;
        } else {
            ret_val = 1;
        }
    } else {
        if (x2 * .20753051f + y2 > 4.20249292f) {
            ret_val = 2;
        } else if (y >= fabs(x) * .08f - .12f) {
            ret_val = 3;
        } else {
            ret_val = 4;
        }
    }
    return ret_val;
} /* bfun3_ */


/*! The Voigt line shape. Drayson approximation of the Voigt line
  shape.

    \retval ls_attenuation            The shape function.
    \retval x             Auxillary parameter to store frequency grid.
    \param  f0            Line center frequency.
    \param  gamma         The pressure broadening parameter.
    \param  sigma         The Doppler broadening parameter.
    \param  f_grid        The frequency grid.

    Original ife function call and documention:

    int voigt_vector(int nx, Numeric *X, Numeric Y, Numeric *Z, Numeric fac)
    
    direct translation of the FORTRAN algorithm given in
    Drayson, S. R., Rapid Computation of the Voigt Profile,
    J. Quant. Spectrosc. Radiat. Transfer, Vol. 16, pp. 611-614, 1976
    by Bjoern-Martin Sinnhuber, 13.Mar.96 in Ny-Aalesund, Spitsbergen.

    Modified for vector calculation of a frequency array:

    \verbatim
    --------------------------------------------------------------------
            int    nx        (in)           number of frequencies
            double *X        (in)           distance from line center in 
                                            units of (Doppler halfwidths 
                                            times sqrt(ln 2) )
            double  Y        (in)           Ratio of the collision halfwidth 
                                            to the ( Doppler halfwidth times sqrt(ln2) )
            double *Z        (out)          return array of voigt 
            double fac       (in)           no function, included to be
                                            consistent with the other voigt
                                            profile calculations
    --------------------------------------------------------------------
    \endverbatim

    23.02.98 AvE

    Replaced nx by nf, Z by ls, X by x, and multiplied ls with the factor fac.
    \author Axel von Engeln
    \date 2000-12-06 */ 

/***  ROUTINE COMPUTES THE VOIGT FUNCTION Y/PI*INTEGRAL FROM ***/
/***   - TO + INFINITY OF EXP(-T*T)/(Y*Y+(X-T)*(X-T)) DT     ***/
void lineshape_voigt_drayson(Vector&         ls_attenuation,
                             Vector&,        //ls_phase
                             Vector&,        //ls_dattenuation_dfrequency_term
                             Vector&,        //ls_dphase_dfrequency_term
                             Vector&,        //ls_dattenuation_dpressure_term
                             Vector&,        //ls_dphase_dpressure_term
                             Vector&         x,
                             const Numeric   f0,
                             const Numeric   gamma,
                             const Numeric,
                             const Numeric,
                             const Numeric,
                             const Numeric,
                             const Numeric   sigma,
                             const Numeric,
                             ConstVectorView f_grid,
                             const bool,
                             const bool)

{
  const Index nf = f_grid.nelem();

  // seems not necessary for Doppler correction
  //    extern const Numeric SQRT_NAT_LOG_2;

  // PI
  extern const Numeric PI;

  // constant sqrt(1/pi)
  const Numeric sqrt_invPI =  sqrt(1/PI);

  // constant normalization factor for voigt
  Numeric fac = 1.0 / sigma * sqrt_invPI;

      Numeric B[22+1] = {0.,0.,.7093602e-7};
      Numeric RI[15+1] = {0};
      const  Numeric XN[15+1] = {0.,10.,9.,8.,8.,7.,6.,5.,4.,3.,3.,3.,3.,3.,3.,3.};
      const  Numeric YN[15+1] = {0.,.6,.6,.6,.5,.4,.4,.3,.3,.3,.3,1.,.9,.8,.7,.7};
      Numeric D0[25+1] = {0}, D1[25+1] = {0}, D2[25+1] = {0}, D3[25+1] = {0}, D4[25+1] = {0};
      Numeric HN[25+1] = {0};
      Numeric H = .201;
      const  Numeric XX[3+1] = {0.,.5246476,1.65068,.7071068};
      const  Numeric HH[3+1] = {0.,.2562121,.2588268e-1,.2820948};
      const  Numeric NBY2[19+1] = {0.,9.5,9.,8.5,8.,7.5,7.,6.5
          ,6.,5.5,5.,4.5,4.,3.5,3.,2.5,2.,1.5,1.,.5};
      const  Numeric C[21+1] = {0.,.7093602e-7,-.2518434e-6,.856687e-6,
          -.2787638e-5,.866074e-5,-.2565551e-4,.7228775e-4,
          -.1933631e-3,.4899520e-3,-.1173267e-2,.2648762e-2,
          -.5623190e-2,.1119601e-1,-.2084976e-1,.3621573e-1,
          -.5851412e-1,.8770816e-1,-.121664,.15584,-.184,.2};
      Numeric CO = 0;
      Numeric U, V, UU, VV, Y2, DX;
      int I, J, K, MAX, MIN, N, i1;
      //int TRU = 0;


      // variables needed in original c routine:

      // Ratio of the Lorentz halfwidth to the Doppler halfwidth
      Numeric Y = gamma / sigma;

      // frequency in units of Doppler, Note: Drayson accepts only
      // positive values
      for (i1=0; i1< (int) nf; i1++)
        {
          x[i1] = fabs( (f_grid[i1] - f0) )/ sigma;
        }


      //if (TRU) goto L104;
/***  REGION I. COMPUTE DAWSON'S FUNCTION AT MESH POINTS ***/
      //TRU = 1;
      for (I=1; I<=15; I++)
        RI[I] = -I/2.;
      for (I=1; I<=25; I++)
      {
        HN[I] = H*(I-.5);
        CO = 4.*HN[I]*HN[I]/25.-2.;
        for (J=2; J<=21; J++)
          B[J+1] = CO*B[J]-B[J-1]+C[J];
        D0[I] = HN[I]*(B[22]-B[21])/5.;
        D1[I] = 1.-2.*HN[I]*D0[I];
        D2[I] = (HN[I]*D1[I]+D0[I])/RI[2];
        D3[I] = (HN[I]*D2[I]+D1[I])/RI[3];
        D4[I] = (HN[I]*D3[I]+D2[I])/RI[4];
      }
//L104:
    for (K=0; K<(int) nf; K++)
  {
    if ((x[K]-5.) < 0.) goto L105; else goto L112;
  L105: if ((Y-1.) <= 0.) goto L110; else goto L106;
  L106: if (x[K] > 1.85*(3.6-Y)) goto L112;
      /***  REGION II: CONTINUED FRACTION. COMPUTE NUMBER OF TERMS NEEDED. ***/
      if (Y < 1.45) goto L107;
      I = (int) (Y+Y);
      goto L108;
  L107: I = (int) (11.*Y);
  L108: J = (int) (x[K]+x[K]+1.85);
      MAX = (int) (XN[J]*YN[I]+.46);
      MIN = (int) ( (16 < 21-2*MAX) ? 16 : 21-2*MAX );
      /***  EVALUATE CONTINUED FRACTION ***/
      UU = Y;
      VV = x[K];
      for (J=MIN; J <= 19; J++)
        {
          U = NBY2[J]/(UU*UU+VV*VV);
          UU = Y+U*UU;
          VV = x[K]-U*VV;
        }
      ls_attenuation[K] = UU/(UU*UU+VV*VV)/1.772454*fac;
      continue;
  L110: Y2 = Y*Y;
      if (x[K]+Y >= 5.) goto L113;
      /***  REGION I. COMPUTE DAWSON'S FUNCTION AT X FROM TAYLOR SERIES. ***/
      N = (int) (x[K]/H);
      DX = x[K]-HN[N+1];
      U = (((D4[N+1]*DX+D3[N+1])*DX+D2[N+1])*DX+D1[N+1])*DX+D0[N+1];
      V = 1.-2.*x[K]*U;
      /***  TAYLOR SERIES EXPANSION ABOUT Y=0.0 ***/
      VV = exp(Y2-x[K]*x[K])*cos(2.*x[K]*Y)/1.128379-Y*V;
      UU = -Y;
      MAX = (int) (5.+(12.5-x[K])*.8*Y);
      for (I=2; I<=MAX; I+=2)
        {
          U = (x[K]*V+U)/RI[I];
          V = (x[K]*U+V)/RI[I+1];
          UU = -UU*Y2;
          VV = VV+V*UU;
        }
      ls_attenuation[K] = 1.128379*VV*fac;
      continue;
  L112: Y2 = Y*Y;
      if (Y < 11.-.6875*x[K]) goto L113;
      /***  REGION IIIB: 2-POINT GAUSS-HERMITE QUADRATURE. ***/
      U = x[K]-XX[3];
      V = x[K]+XX[3];
      ls_attenuation[K] = Y*(HH[3]/(Y2+U*U)+HH[3]/(Y2+V*V))*fac;
      continue;
      /***  REGION IIIA: 4-POINT GAUSS-HERMITE QUADRATURE. ***/
  L113: U = x[K]-XX[1];
      V = x[K]+XX[1];
      UU = x[K]-XX[2];
      VV = x[K]+XX[2];
      ls_attenuation[K] = Y*(HH[1]/(Y2+U*U)+HH[1]/(Y2+V*V)+HH[2]/(Y2+UU*UU)+HH[2]/(Y2+
                                                                      VV*VV))*fac;
      continue;
  }
}


/*! Chi factors according to Cousin

    The CO2-CO2 self-broadening is neglected. Broadening factors for both
    O2 and N2 are considered, assuming 79% N2 and 21% O2.

    \param  chi  Out: The chi factor
    \param  df   (f-f0) in Hz.

    \author Patrick Eriksson 
    \date 2000-09-07 
*/
void chi_cousin(
            Numeric&   chi,
      const Numeric&   df )
{
  // Conversion factor from Hz to cm-1
  extern const Numeric HZ2CM;

  const Numeric n2 = 0.79;
  const Numeric o2 = 0.21;

  const Numeric   df_cm     = df * HZ2CM;
  const Numeric   df_cm_abs = fabs( df_cm );

  chi = 0;  

  // N2 term
  if( df_cm_abs <= 5 )
    { chi += n2; }
  else if( df_cm_abs <= 22 )
    { chi += n2 * 1.968 * exp( -0.1354 * df_cm_abs ); }
  else if( df_cm_abs <= 50 )
    { chi += n2 * 0.160 * exp( -0.0214 * df_cm_abs ); }
  else
    { chi += n2 * 0.162 * exp( -0.0216 * df_cm_abs ); }

  // O2 term
  if( df_cm_abs <= 5 )
    { chi += o2; }
  else if( df_cm_abs <= 22 )
    { chi += o2 * 1.968 * exp( -0.1354 * df_cm_abs ); }
  else if( df_cm_abs <= 50 )
    { chi += o2 * 0.160 * exp( -0.0214 * df_cm_abs ); }
  else
    { chi += o2 * 0.162 * exp( -0.0216 * df_cm_abs ); }
}



/*! A CO2 IR line shape.

    \retval ls_attenuation            The shape function.
    \param  f0            Line center frequency.
    \param  gamma         The pressure broadening parameter.
    \param  f_grid        The frequency grid.

    \author Patrick Eriksson 
    \date 2000-09-04 */
void lineshape_CO2_lorentz(Vector&         ls_attenuation,
                           Vector&,        //ls_phase
                           Vector&,        //ls_dattenuation_dfrequency_term
                           Vector&,        //ls_dattenuation_dfrequency_term
                           Vector&,        //ls_dphase_dpressure_term
                           Vector&,        //ls_dphase_dpressure_term
                           Vector&,        //X
                           const Numeric   f0,
                           const Numeric   gamma,
                           const Numeric,
                           const Numeric,
                           const Numeric,
                           const Numeric,
                           const Numeric,  //sigma
                           const Numeric,
                           ConstVectorView f_grid,
                           const bool,
                           const bool)
{
  const Index nf = f_grid.nelem();

  // PI:
  extern const Numeric PI;

  // Some constant variables
  const Numeric gamma2 = gamma * gamma;
  const Numeric fac = gamma/PI;

  for ( Index i=0; i<nf; ++i )
    {
      // f-f0
      const Numeric df = f_grid[i] - f0;

      // The chi factor
      Numeric chi;
      chi_cousin( chi, df );      

      // chi * Lorentz
      ls_attenuation[i] =  chi * fac / ( df * df + gamma2 );
    }
}



/*! A CO2 IR line shape.

    \retval ls_attenuation            The shape function.
    \retval X             Auxillary parameter, only used in Voigt fct.
    \param  f0            Line center frequency.
    \param  gamma         The pressure broadening parameter.
    \param  sigma         The Doppler broadening parameter. (Not used.)
    \param  f_grid        The frequency grid.

    \author Patrick Eriksson 
    \date 2000-09-04 */
void lineshape_CO2_drayson(Vector&         ls_attenuation,
                           Vector&,        //ls_phase,
                           Vector&,        //ls_dattenuation_dfrequency_term
                           Vector&,        //ls_dphase_dfrequency_term
                           Vector&,        //ls_dattenuation_dpressure_term
                           Vector&,        //ls_dphase_dpressure_term
                           Vector&         X,
                           const Numeric   f0,
                           const Numeric   gamma,
                           const Numeric,
                           const Numeric,
                           const Numeric,
                           const Numeric,
                           const Numeric   sigma,
                           const Numeric,
                           ConstVectorView f_grid,
                           const bool,
                           const bool)
{
    Vector tmp(0);
    lineshape_voigt_drayson(  ls_attenuation, tmp,tmp,tmp,tmp,tmp, X, f0, gamma,0.,1.,0.,0., sigma, 0., f_grid ,0,0);

  const Index nf = f_grid.nelem();
  for ( Index i=0; i<nf; ++i )
    {
      // f-f0
      const Numeric df = f_grid[i] - f0;

      // The chi factor
      Numeric chi;
      chi_cousin( chi, df );      

      // chi * Lorentz
      ls_attenuation[i] *=  chi;
    }
}


/*! The Voigt and Faraday-Voigt line shape. Based on rewritten Faddeeva 
  function algorithm 916. For more information 
  read:
  MOFREH R. ZAGHLOUL and AHMED N. ALI
  Algorithm 916: Computing the Faddeyeva and Voigt Functions
  ACM Transactions on Mathematical Software, Vol. 38, No. 2, Article 15, Publication date: December 2011.

  The main bulk of code is in Faddeeva.cc and written by Steven G. Johnson of MIT.
  Keep Faddeeva.{cc,hh} up to date in speed and accuracy by watching:
  http://ab-initio.mit.edu/Faddeeva sometimes.

  \retval ls_attenuation        The shape function.
  \retval ls_phase              The phase shape function.
  \param  f0                    Line center frequency.
  \param  gamma                 The pressure broadening parameter.
  \param  sigma                 The Doppler broadening parameter.
  \param  f_grid                The frequency grid.

  \author Richard Larsson 2013-01-17

 */ 
void faddeeva_algorithm_916(    Vector&         ls_attenuation,
                                Vector&         ls_phase,
                                Vector&         ls_dattenuation_dfrequency_term,
                                Vector&         ls_dphase_dfrequency_term,
                                Vector&         ls_dattenuation_dpressure_term,
                                Vector&         ls_dphase_dpressure_term,
                                Vector&,
                                const Numeric   f0,
                                const Numeric   gamma,
                                const Numeric,
                                const Numeric,
                                const Numeric,
                                const Numeric,
                                const Numeric   sigma,
                                const Numeric,
                                ConstVectorView f_grid,
                                const bool      do_phase,
                                const bool      do_partials)

{
    const Index nf = f_grid.nelem();
    
    // PI
    extern const Numeric PI;
    
    // constant sqrt(1/pi)
    const Numeric sqrt_invPI =  sqrt(1/PI);
    
    // constant normalization factor for voigt
    const Numeric fac = sqrt_invPI / sigma;
    
    // Ratio of the Lorentz halfwidth to the Doppler halfwidth
    const Numeric y = gamma / (sigma);
    
    // frequency in units of Doppler
    for (Index ii=0; ii<nf; ii++)
    {
        const Numeric x = (f_grid[ii] - f0) / sigma;
        const std::complex<Numeric> w =
        Faddeeva::w(std::complex<Numeric>(x, y));
        
        ls_attenuation[ii] = fac * w.real();
        if(do_phase||do_partials)
            ls_phase[ii]   = fac * w.imag();
        if(do_partials)
        {
            // Derivatives are from a paper that gives w(y-ix) = Fa + iFb but our formalism use
            // w(x+iy) = Fa + iFb.  Thus signs on Fb-derivatives are difficult...
            
            // dFa/dx
            ls_dattenuation_dfrequency_term[ii] = 2.*(y*ls_phase[ii]-x*ls_attenuation[ii]);
            
            // dFa/dy
            ls_dattenuation_dpressure_term[ii]  = 2.*( y*ls_attenuation[ii] + x*ls_phase[ii] - fac*sqrt_invPI);
            
            if(do_phase)
            {
                // dFb/dx
                ls_dphase_dfrequency_term[ii] = - ls_dattenuation_dpressure_term[ii];
                
                // dFb/dy
                ls_dphase_dpressure_term[ii]  =   ls_dattenuation_dfrequency_term[ii];
            }
        }
    }
}


// Helpers connected to internal partials.
void w_x_plus_iy_dT(Vector& dx_dT,
                    Numeric& dy_dT,
                    Numeric& dFu_dT,
                    ConstVectorView f,     
                    const Numeric& f0,    // Note that this is NOT the line center, but the shifted center
                    const Numeric& sigma, 
                    const Numeric& dPF_dT,    // pressure shift temperature derivative
                    const Numeric& dDF_dT,    // line mixing shift temperature derivative
                    const Numeric& dsigma_dT, // derivative of sigma with temperature
                    const Numeric& gamma,
                    const Numeric& dgamma_dT) // derivative of gamma with temperature
{
    // Function is w(x+iy), and dw/dx and dw/dy are both known already.
    // This function serves to find dx/dT
    
    // x = (f-f0) / sigma;
    
    // f0 contains many different terms:
    // 1) Magnetic field splitting (not T-dependent)
    // 2) Pressure frequency shift (T-dependent)
    // 3) Pressure line mixing shift (T-dependent)
    // f contains the frequency grid influenced by wind
    // sigma is the halfwidth Doppler broadening
    
    // Since both denominator and nominator depends on temperature, we need to split the derivative into two terms
    // d(f-f0)/dT * 1/sigma  and  (f-f0) * d(1/sigma)/dT
    
    for(Index iv=0; iv<f.nelem(); iv++)
        dx_dT[iv] = ( - dPF_dT - dDF_dT  -  (f[iv] - f0) / sigma * dsigma_dT ) / sigma;
    
    // y = gamma / sigma
    
    // Since both denominator and nominator depends on temperature, we need to split the derivative into two terms
    // d(gamma)/dT * 1/sigma  and  (gamma) * d(1/sigma)/dT
    
    dy_dT = ( dgamma_dT  -  gamma / sigma * dsigma_dT )  / sigma;
    
    // This line shape is normalization depends on temperature and is
    dFu_dT = - 1.0 / sigma * dsigma_dT; 
    // NOTE: The factor is 1/sqrt(pi)/sigma, its derivative is then - 1 / sqrt(pi) / sigma^2 * dsigma_dT, but
    //       the factor is already in the lineshape that this will be multiplied to, so
    //       derivative divided by the factor and only the above remains!  Which should only be -1/2T for this case.
}
void w_x_plus_iy_dF(Numeric& dx_dF,
                    const Numeric& sigma)
{
    // Function is w(x+iy), and dw/dx and dw/dy are both known already.
    // This function serves to find dx/dT
    
    // x = (f-f0) / sigma;
    
    dx_dF = 1.0/sigma;
    
    // No other terms depend on f
}
void w_x_plus_iy_dF0(Numeric& dx_dF0,
                     const Numeric& sigma)
{
    // Function is w(x+iy), and dw/dx and dw/dy are both known already.
    // This function serves to find dx/dT
    
    // x = (f-f0) / sigma;
    
    dx_dF0 = -1.0/sigma;
    
    // No other terms depend on f
}
void w_x_plus_iy_dgamma(Numeric& dy_dgamma,
                        const Numeric& sigma)
{
    // Function is w(x+iy), and dw/dx and dw/dy are both known already.
    // This function serves to find dy/dgamma
    
    // y = gamma / sigma;
    
    dy_dgamma = 1.0/sigma;
    
    // No other terms depend on gamma
}
void w_x_plus_iy_dH(Numeric& dx_dH, 
                    const Numeric& sigma,
                    const Numeric& df0_dH)
{
    // Function is w(x+iy), and dw/dx and dw/dy are both known already.
    // This function serves to find dx/dmag
    
    // x = (f-f0) / sigma;
    
    // f0 contains many different terms:
    // 1) Magnetic field splitting (mag-dependent)
    // 2) Pressure frequency shift (not mag-dependent)
    // 3) Pressure line mixing shift (not mag-dependent)
    // f contains the frequency grid influenced by wind
    // sigma is the halfwidth Doppler broadening
    
    dx_dH = - df0_dH/sigma;
    
    // No other terms depend on mag
}
void w_x_plus_iy_dDF(Numeric& dx_dDF,       
                     const Numeric& sigma)
{
    // Function is w(x+iy), and dw/dx and dw/dy are both known already.
    // This function serves to find dx/dDF
    
    // x = (f-f0) / sigma;
    
    // f0 contains many different terms:
    // 1) Magnetic field splitting (not DF-dependent)
    // 2) Pressure frequency shift (not DF-dependent)
    // 3) Pressure line mixing shift (DF-dependent)
    // f contains the frequency grid influenced by wind
    // sigma is the halfwidth Doppler broadening
    
    dx_dDF = - 1.0/sigma;
    
    // No other terms depend on DF
}


/*! The Voigt and Faraday-Voigt line shape. Based on Complex Error Function of
  X+i*Y.
  
  Y must be positive.
  
  Reference for equations:
  Hui, Armstrong and Wray, JQSRT V.19, P.509 (1978).
  
  ----------------------
  CERROR(z)=CERROR(x+iy)=CERROR(x,y)=w(x,y)
  g(nu,nu_0) = (1/DopHW)* sqrt(ln2/pi)*(Re(w(x,y) + Py*Im(w(x,y)) 
  y=(P*LorHW /DopHW)*sqrt(ln2)
  x=(nu-nu_0)/DopHW)*sqrt(ln2)                      
  NB Re(w)=Voigt
  ----------------------
  
  \retval ls_attenuation        The shape function.
  \retval ls_phase              The phase shape function.
  \retval xvector               Auxillary parameter to store frequency grid.
  \param  f0                    Line center frequency.
  \param  gamma                 The pressure broadening parameter.
  \param  sigma                 The Doppler broadening parameter.
  \param  f_grid                The frequency grid.

  \author P. Rosenkranz 1988-02-19 
  
  (Added by Richard Larsson 2014-03-13 with slight moderation)

 */ 
void hui_etal_1978_lineshape( Vector&         ls_attenuation,
                              Vector&         ls_phase,
                              Vector&,        //ls_dattenuation_dfrequency_term
                              Vector&,        //ls_dphase_dfrequency_term
                              Vector&,        //ls_dattenuation_dpressure_term
                              Vector&,        //ls_dphase_dpressure_term
                              Vector&         xvector,
                              const Numeric   f0,
                              const Numeric   gamma,
                              const Numeric,
                              const Numeric,
                              const Numeric,
                              const Numeric,
                              const Numeric   sigma,
                              const Numeric,
                              ConstVectorView f_grid,
                              const bool do_phase,
                              const bool)

{
    const Index nf = f_grid.nelem();
    
    // PI
    extern const Numeric PI;
    
    // constant sqrt(1/pi)
    const Numeric sqrt_invPI =  sqrt(1/PI);
    
    // constant normalization factor for voigt
    const Numeric fac = sqrt_invPI / sigma;
    
    // Ratio of the Lorentz halfwidth to the Doppler halfwidth
    const Numeric y = gamma / (sigma);
    
    // frequency in units of Doppler
    for (Index ii=0; ii<nf; ii++)
    {
        xvector[ii] = (f_grid[ii] - f0) / sigma;
        
        // Note that this works but I don't know why there is a difference
        // between the theory described above and this practical implementation.
        const std::complex<Numeric> z(y , - xvector[ii]);
        
        const std::complex<Numeric> A = (((((.5641896*z+5.912626)*z+30.18014)*z+
              93.15558)*z+181.9285)*z+214.3824)*z+122.6079;
        const std::complex<Numeric> B = ((((((z+10.47986)*z+53.99291)*z+170.3540)*z+
              348.7039)*z+457.3345)*z+352.7306)*z+122.6079;
        const std::complex<Numeric> C = A / B;
        
        ls_attenuation[ii] = fac * C.real();
        if(do_phase)
            ls_phase[ii]       = fac * C.imag();
    }
}


/*! Hartmann-Tran line shape. Based on
    [1] Ngo NH, Lisak D, Tran H, Hartmann J-M. An isolated line-shape model
    to go beyond the Voigt profile in spectroscopic databases and radiative 
    transfer codes. J Quant Spec & Radiat Transfer 2013;129:89-100.                
    [2] Tran H, Ngo NH, Hartmann J-M. Efficient computation of some speed-dependent 
    isolated line profiles. J Quant Spec & Radiat Transfer 2013;129:199-203.
      
    N.B.  Input is not handled properly yet but assumed standard where nothing is known.
    This is indicated by all numerics set to constants at the start of the code.
    
  \retval ls_attenuation        The shape function.
  \retval ls_phase              The phase shape function.
  \param  f0                    Line center frequency.
  \param  gamma_0               The pressure broadening parameter.
  \param  gamma_2               The pressure broadening speed dependent parameter.
  \param  eta                   A correlation value.
  \param  df_0                  The pressure frequency shift parameter.
  \param  df_2                  The pressure frequency shift speed-dependent parameter.
  \param  gamma_D               The Doppler broadening parameter.
  \param  f_VC                  Collisional frequency parameter
  \param  f_grid                The frequency grid.

  \author Richard Larsson 2015-07-02

 */ 
void hartmann_tran_lineshape(   Vector&         ls_attenuation,
                                Vector&         ls_phase,
                                Vector&         ls_dattenuation_dfrequency_term,
                                Vector&         ls_dphase_dfrequency_term,
                                Vector&         ls_dattenuation_dpressure_term,
                                Vector&         ls_dphase_dpressure_term,
                                Vector&         xvector, 
                                const Numeric   f0,      
                                const Numeric   gamma_0, 
                                const Numeric   gamma_2, //untested
                                const Numeric   eta,     //untested
                                const Numeric   df_0,    //untested
                                const Numeric   df_2,    //untested
                                const Numeric   gamma_D, 
                                const Numeric   f_VC,    //untested
                                ConstVectorView f_grid,  
                                const bool      do_phase,
                                const bool      do_partials)

{
    // NOTE: f_VC is the only non-scaled variable in the reference document.  
    // If given as something other than cm-1 it might behave strange.  Still not tested,
    // just a future "fixme".  The others are probably OK, but I am very unsure on this 
    // as well.
    
    // NOTE: need more outputs for the partial derivatives?  This requires redesign so that
    // there is no longer such abstract derivatives.  That would slow the other methods down.
    // So a better solution is to just allow complex vector outputs, and then save:
    // w(iZ-) and w(iZ+).  Their partials can be found from knowing these two --- Y and X has to
    // be re-derived in code.  The Forthomme etal 2015 method of finding the outputs with 
    // regards to the other parameters of the equations can be done.  See their supplementary 
    // material for how to calculate the many partial derivatives (or redo the work yourself).
    // This does mean that, internally, this function will output derivatives that have 
    // a very different meaning from all the others.
    // This might not work out, since Line mixing makes mixing of these parameters excruciating.
    
    const bool test1 = gamma_2==df_2&&gamma_2==f_VC&&gamma_2==eta&&gamma_2==0;
    
    if( gamma_0==df_0 && gamma_0==gamma_2 && test1 )
    {
        throw std::runtime_error("Doppler line shape not supported in ARTS using HTP.\n");
    }
    else if( test1 )  // Standard Voigt function so use standard Voigt function call.  NOTE:  Partial derivations must know this.
    {
        faddeeva_algorithm_916(   ls_attenuation, ls_phase, ls_dattenuation_dfrequency_term,
                                  ls_dphase_dfrequency_term, ls_dattenuation_dpressure_term,
                                  ls_dphase_dpressure_term, xvector, f0, gamma_0,
                                  gamma_2, eta, df_0, df_2, gamma_D, f_VC, f_grid,
                                  do_phase, do_partials);
        return;
    }
    else //This is mostly untested and there are likely errors.  Like for c_2 == 0 and for eta == 1
    {
        if(do_partials)
            throw std::runtime_error("It is not yet supported to do partial derivation with\n"
            "HTP line shape parameters.");
        
        // Constants
        extern const Numeric PI;
        extern const Numeric SPEED_OF_LIGHT;
        const std::complex<Numeric> i(0.,1.);
        
        // Direct pressure term
        const std::complex<Numeric> c_0(gamma_0, df_0);
        
        // Indirect pressure term
        const std::complex<Numeric> c_2(gamma_2, df_2);
        
        // Constant c_0 helper
        const std::complex<Numeric> c_0_tilde = (1.-eta)*(c_0-1.5*c_2)+f_VC;
        
        // Constant c_2 helper
        const std::complex<Numeric> c_2_tilde = (1.-eta)*c_2;
        
        // Average speed of molecule
        const Numeric v_a0 = SPEED_OF_LIGHT * gamma_D / f0;
        
        // Constant sqrt(pi)
        const Numeric sqrtPI =  sqrt(PI);
        
        // Some calculation constants
        const std::complex<Numeric> c1 = sqrtPI/gamma_D;
        const std::complex<Numeric> v_a02 = v_a0*v_a0;
        const std::complex<Numeric> c4 = f_VC-eta*(c_0-1.5*c_2);
        const std::complex<Numeric> c5 = eta*c_2/v_a02;
        
        // Special case when no speed dependency or no correlation
        if(c_2_tilde == std::complex<Numeric>(0., 0.))
        {
            for(Index nf = 0; nf<f_grid.nelem(); nf++)
            {
                const std::complex<Numeric> Z1 = 
                std::complex<Numeric>(0.,(f0-f_grid[nf])/gamma_D) + c_0_tilde;
                
                // NB.  Z1 → infinity is treated specially in Tran etal.  This I think already happens in
                // the Faddeeeva function as well.  So I will do nothing special
                
                const std::complex<Numeric> w1 = Faddeeva::w(i*Z1);
                
                const std::complex<Numeric> A = c1*w1;
                
                const std::complex<Numeric> B = c1*v_a02* ( (1.-Z1*Z1)*w1 + Z1/sqrtPI );
                
                // Hartmann-Tran line shape
                const std::complex<Numeric> HTP = 1./PI * A / ( 1. - c4*A + c5*B );
                
                // Output
                ls_attenuation[nf] = HTP.real();
                if(do_phase||do_partials)
                    ls_phase[nf] = HTP.imag();
            }
            return;
        }
        
        // Y variable otherwise
        const std::complex<Numeric> sqrtY = gamma_D/(2. * c_2_tilde);
        const std::complex<Numeric> Y = sqrtY*sqrtY;
        const Numeric absY = abs(Y);
        const std::complex<Numeric> sqrtY_t2 = 2.*sqrtY;
        const std::complex<Numeric> c2 = v_a02/c_2_tilde/c_2_tilde;
        const std::complex<Numeric> c3 = sqrtPI/sqrtY_t2;
        
        for(Index nf = 0; nf<f_grid.nelem(); nf++)
        {
            const Numeric df = f0-f_grid[nf];
            
            // X variable
            const std::complex<Numeric> X = ( std::complex<Numeric>(0.,df) + c_0_tilde )/c_2_tilde;
            
            // Tran etal says this must be between 3e-8 and 10e15 for numerical stability in their routine.
            // This is not the case for us since we have a different code, but is arbitrarily used anyways.
            const Numeric numerical_test = abs(X)/absY;
            
            // Same variables in all cases
            std::complex<Numeric> A, B;
            if(numerical_test < 3e-8)
            {
                // Z_plus 
                // Note that sqrt(X+Y)+sqrt(Y) → 2*sqrt(Y) in this case
                // This simplification is ignored below but is faster. Should it be used?
                const std::complex<Numeric> Z_plus  = sqrt(X+Y)+sqrtY;
                
                // Z_minus
                const std::complex<Numeric> Z_minus = 
                std::complex<Numeric>(0.,df/gamma_D)+c_0_tilde;
                
                // w(Z_minus)
                const std::complex<Numeric> w_plus  =  Faddeeva::w(i*Z_plus);
                
                // w(Z_minus)
                const std::complex<Numeric> w_minus =  Faddeeva::w(i*Z_minus);
                
                // A
                A = c1 * (w_minus - w_plus);
                
                // B
                B = c2 *(-1.+ c3*((1.-Z_minus*Z_minus)*w_minus - (1.-Z_plus*Z_plus)*w_plus));
            }
            else if(numerical_test > 1e15)
            {
                // Note that the Z:s are very similar but 
                // Tran etal still makes a difference between them
                const std::complex<Numeric> sqrtX = sqrt(X);
                const std::complex<Numeric>  Z2   = sqrt(X+Y);
                
                const std::complex<Numeric> w1 = Faddeeva::w(i*sqrtX);
                const std::complex<Numeric> w2 = Faddeeva::w(i*Z2);
                
                // Note that there is a warning for these equations as X → infinity that is ignored for now.
                
                std::complex<Numeric> Aconst = (1./sqrtPI - sqrtX*w1);
                
                // A
                A = 2.*c1 * Aconst;
                
                // B
                B = c2 * (-1. + 2.*sqrtPI*(1.-X-2.*Y)*Aconst + 2.*sqrtPI*Z2*w2);
                
            }
            else
            {
                // Z_plus
                const std::complex<Numeric> Z_plus  = sqrt(X+Y)+sqrtY;
                
                // Z_minus
                const std::complex<Numeric> Z_minus = Z_plus-sqrtY_t2;
                
                // w(Z_minus)
                const std::complex<Numeric> w_plus  =  Faddeeva::w(i*Z_plus);
                
                // w(Z_minus)
                const std::complex<Numeric> w_minus =  Faddeeva::w(i*Z_minus);
                
                // A
                A = c1 * (w_minus - w_plus);
                
                // B
                B = c2 *(-1.+ c3*((1.-Z_minus*Z_minus)*w_minus - (1.-Z_plus*Z_plus)*w_plus));
            }
            
            // Hartmann-Tran line shape
            const std::complex<Numeric> HTP = 1./PI * A / ( 1. - c4*A + c5*B);
            
            // Output
            ls_attenuation[nf] = HTP.real();
            if(do_phase||do_partials)
                ls_phase[nf] = HTP.imag();
        }
    }
}


/*! The O2 non-resonant line shape.  Should be VVW/2 so do not use this call...

    \retval ls_attenuation              The shape function.
    \retval ls_phase                    The shape function.
    \retval X                           Auxillary parameter, only used in Voigt fct.
    \param  f0                          Line center frequency.
    \param  gamma                       The pressure broadening parameter.
    \param  sigma                       The Doppler broadening parameter. (Not used.)
    \param  f_grid                      The frequency grid.

    \author Richard Larsson
    \date 2014-10-02 */
void lineshape_o2nonresonant( Vector&         ls_attenuation,
                              Vector&,        //ls_phase _U_,
                              Vector&,        //ls_dattenuation_dfrequency_term
                              Vector&,        //ls_dphase_dfrequency_term
                              Vector&,        //ls_dattenuation_dpressure_term
                              Vector&,        //ls_dphase_dpressure_term
                              Vector&,        //X,
                              const Numeric   f0,
                              const Numeric   gamma,
                              const Numeric,
                              const Numeric,
                              const Numeric,
                              const Numeric,
                              const Numeric,   //sigma
                              const Numeric,
                              ConstVectorView f_grid,
                              const bool,
                              const bool)
{
  const Index nf = f_grid.nelem();

  // PI:
  extern const Numeric PI;

  //  assert( ls.nelem() == nf );
  
  Numeric fac = gamma/PI;

  for ( Index i=0; i<nf; ++i )
    {
      ls_attenuation[i] =  fac / ( (f_grid[i]-f0) * (f_grid[i]-f0) );
    }
}


//------------------------------------------------------------------------
// Normalization Functions 
//------------------------------------------------------------------------


/*!  No normalization of the lineshape function.

    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_grid The frequency grid.
    \param  T      Temperature (unused here)

    \author Axel von Engeln 30.11.2000 */
void lineshape_norm_no_norm(Vector&         fac,
                            const Numeric,  //f0
                            ConstVectorView,//f_grid
                            const Numeric)  //T
{
  fac = 1.0;
}
void lineshape_norm_no_norm_dT(Vector&         fac,
                               const Numeric,  //f0
                               ConstVectorView,//f_grid
                               const Numeric)  //T
{
    fac = 0.0;
}
void lineshape_norm_no_norm_dF(Vector&         fac,
                               const Numeric,  //f0
                               ConstVectorView,//f_grid
                               const Numeric)  //T
{
    fac = 0.0;
}
void lineshape_norm_no_norm_dF0(Vector&         fac,
                                const Numeric,  //f0
                                ConstVectorView,//f_grid
                                const Numeric)  //T
{
    fac = 0.0;
}


/*!  Linear normalization factor of the lineshape function with f/f0.

    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_grid The frequency grid.
    \param  T      (unused here)

    \author Axel von Engeln 30.11.2000 */
void lineshape_norm_linear(Vector&         fac,
                           const Numeric   f0,
                           ConstVectorView f_grid,
                           const Numeric)  //T
{
  const Index nf = f_grid.nelem();

  // Abs(f0) is constant in the loop:
  const Numeric abs_f0 = abs(f0);

  for ( Index i=0; i<nf; ++i )
    {
      fac[i] = f_grid[i] / abs_f0;
    }
}
void lineshape_norm_linear_dT(Vector&         fac,
                             const Numeric, // f0
                             ConstVectorView, //f_grid
                             const Numeric)  //T
{
    fac = 0.0;
}
void lineshape_norm_linear_dF(Vector&         fac,
                           const Numeric   f0,
                           ConstVectorView f_grid,
                           const Numeric)  //T
{
    const Index nf = f_grid.nelem();
    
    // Abs(f0) is constant in the loop:
    const Numeric abs_f0 = abs(f0);
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = 1.0 / abs_f0;
    }
}
void lineshape_norm_linear_dF0(Vector&         fac,
                               const Numeric   f0,
                               ConstVectorView f_grid,
                               const Numeric)  //T
{
    const Index nf = f_grid.nelem();
    
    // Abs(f0) is constant in the loop:
    const Numeric dabs_f0 = -f0/abs(f0)/abs(f0)/abs(f0); //-sign(f0)/abs(f0)^2 Why is this abs?
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = f_grid[i] * dabs_f0;
    }
}

/*!  Quadratic normalization factor of the lineshape function with (f/f0)^2*hf0/2kT/sinh(hf0/2kT).
     in case of the quadratic normalization factor use the
     so called 'microwave approximation' of the line intensity
     given by
     P. W. Rosenkranz, Chapter 2, Eq.(2.16), in M. A. Janssen,
     Atmospheric Remote Sensing by Microwave Radiometry,
     John Wiley & Sons, Inc., 1993
 
    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_grid The frequency grid.
    \param  T      Temperature (unused here)

    \author Axel von Engeln 30.11.2000 */
void lineshape_norm_quadratic_Rosenkranz(Vector&         fac,
                                         const Numeric   f0,
                                         ConstVectorView f_grid,
                                         const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    const Numeric mafac = (PLANCK_CONST * f0) / (2.000e0 * BOLTZMAN_CONST * T) /
                                sinh((PLANCK_CONST * f0) / (2.000e0 * BOLTZMAN_CONST * T));
    const Index nf = f_grid.nelem();
    // don't do this for the whole loop
    const Numeric f0_2 = f0 * f0;
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = mafac * (f_grid[i] * f_grid[i]) / f0_2;
    }
}
void lineshape_norm_quadratic_Rosenkranz_dT(Vector&        fac,
                                           const Numeric   f0,
                                           ConstVectorView f_grid,
                                           const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    const Numeric hf0 = PLANCK_CONST*f0, kT =T*BOLTZMAN_CONST;
    const Numeric sinh_term=sinh((hf0)/(2*kT));
    const Numeric cosh_term = cosh((hf0)/(2.0*kT));
    const Numeric dmafac_dT = -(hf0*(2.0*kT*sinh_term - hf0*cosh_term))/(4.0*T*kT*kT*sinh_term*sinh_term);
    const Index nf = f_grid.nelem();
    // don't do this for the whole loop
    const Numeric f0_2 = f0 * f0;
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = dmafac_dT * (f_grid[i] * f_grid[i]) / f0_2;
    }
}
void lineshape_norm_quadratic_Rosenkranz_dF(Vector&         fac,
                                            const Numeric   f0,
                                            ConstVectorView f_grid,
                                            const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    const Numeric mafac = (PLANCK_CONST * f0) / (2.000e0 * BOLTZMAN_CONST * T) /
    sinh((PLANCK_CONST * f0) / (2.000e0 * BOLTZMAN_CONST * T));
    const Index nf = f_grid.nelem();
    // don't do this for the whole loop
    const Numeric f0_2 = f0 * f0;
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = mafac * (2.0 * f_grid[i]) / f0_2;
    }
}
void lineshape_norm_quadratic_Rosenkranz_dF0(Vector&         fac,
                                             const Numeric   f0,
                                             ConstVectorView f_grid,
                                             const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    const Numeric kT = 2.0*BOLTZMAN_CONST*T;
    const Index nf = f_grid.nelem();
    // don't do this for the whole loop
    const Numeric term1 = sinh((f0*PLANCK_CONST)/kT);
    
    const Numeric ddenom = -(PLANCK_CONST*(kT*term1 + f0*PLANCK_CONST*cosh((f0*PLANCK_CONST)/kT)))/(f0*f0*kT*kT*term1*term1);
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = (f_grid[i] * f_grid[i]) * ddenom;
    }
}

/*!  Van Vleck Huber normalization factor of the lineshape function
  with (f*tanh(h*f/(2*k*T))) / (f0*tanh(h*f0/(2*k*T))). The
  denominator is a result of catalogue intensities. See P. Rayer, The
  VVH and VVW Spectral Functions, Atmospheric Millimeter and
  Sub-Millimeter Wave Radiative Transfer Modeling II, Editors:
  P. Eriksson, S. Buehler, Berichte aus derm Institut fuer
  Umweltphysik, Band 4, 2001.

    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_grid The frequency grid.
    \param  T      Temperature

    \author Axel von Engeln 2003-07-28 */
void lineshape_norm_VVH(Vector&         fac,
                        const Numeric   f0,
                        ConstVectorView f_grid,
                        const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const Index nf = f_grid.nelem();
    
    // 2kT is constant for the loop
    const Numeric kT = 2.0 * BOLTZMAN_CONST * T;
    
    // denominator is constant for the loop
    const Numeric denom = abs(f0) * tanh( PLANCK_CONST * abs(f0) / kT );
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = f_grid[i] * tanh( PLANCK_CONST * f_grid[i] / kT ) /
        denom;
    }
}
void lineshape_norm_VVH_dT(Vector&         fac,
                           const Numeric   f0,
                           ConstVectorView f_grid,
                           const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const Index nf = f_grid.nelem();
    
    // Various constants for the loop
    const Numeric kT = BOLTZMAN_CONST * T;
    const Numeric hf0 = PLANCK_CONST * f0;
    const Numeric tanh_term2=tanh((hf0)/(2.0*kT));
    const Numeric tanh_term2_squared=tanh_term2*tanh_term2;
    
    // denominator is constant for the loop
    const Numeric denom1 = 2.0*T*kT*f0*tanh_term2,denom2=2.0*T*kT*tanh_term2_squared;
    
    for ( Index i=0; i<nf; ++i )
    {
        const Numeric hf = PLANCK_CONST * f_grid[i], tanh_term1=tanh((hf)/(2.0*kT)), tanh_term1_squared=tanh_term1*tanh_term1;
        
        fac[i] = (f_grid[i]*hf*(tanh_term1_squared - 1.0))/denom1 - (hf*tanh_term1*(tanh_term2_squared - 1.0))/denom2;
    }
}
void lineshape_norm_VVH_dF(Vector&         fac,
                           const Numeric   f0,
                           ConstVectorView f_grid,
                           const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const Index nf = f_grid.nelem();
    
    // Various constants for the loop
    const Numeric kT  = BOLTZMAN_CONST * T, 
                  hf0 = PLANCK_CONST*f0,
                  tanh_term2=tanh((hf0)/(2.0*kT)),
                  denom =(2*kT*f0*tanh_term2);
    
    for ( Index i=0; i<nf; ++i )
    {
        const Numeric hf = PLANCK_CONST*f_grid[i], 
                      tanh_term1=tanh((hf)/(2.0*kT)),
                      tanh_term1_squared=tanh_term1*tanh_term1;
        fac[i] = (hf + 2*kT*tanh_term1 - hf*tanh_term1_squared)/denom;
    }
}
void lineshape_norm_VVH_dF0(Vector&         fac,
                            const Numeric   f0,
                            ConstVectorView f_grid,
                            const Numeric   T)
{
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const Index nf = f_grid.nelem();
    
    // 2kT is constant for the loop
    const Numeric kT = 2.0 * BOLTZMAN_CONST * T, term1 = tanh((PLANCK_CONST*abs(f0))/kT);
    
    // denominator is constant for the loop
    const Numeric ddenom = (PLANCK_CONST*f0/abs(f0)*(term1*term1 - 1))/(kT*abs(f0)*term1*term1) - f0/abs(f0)/(abs(f0)*abs(f0)*term1);
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = f_grid[i] * tanh( PLANCK_CONST * f_grid[i] / kT ) * ddenom;
    }
}

/*!  Van Vleck Weiskopf normalization factor of the lineshape function
 * with (f*f) / (f0*f0). The
 * denominator is a result of catalogue intensities. See P. Rayer, The
 * VVH and VVW Spectral Functions, Atmospheric Millimeter and
 * Sub-Millimeter Wave Radiative Transfer Modeling II, Editors:
 * P. Eriksson, S. Buehler, Berichte aus derm Institut fuer
 * Umweltphysik, Band 4, 2001.
 * 
 *   \retval fac    Normalization factor to the lineshape function.
 *   \param  f0     Line center frequency.
 *   \param  f_grid The frequency grid.
 *   \param  T      Temperature
 * 
 *   \author Axel von Engeln 2003-07-28 */
void lineshape_norm_VVW(Vector&         fac,
                        const Numeric   f0,
                        ConstVectorView f_grid,
                        const Numeric)
{
    
    const Index nf = f_grid.nelem();
    
    // denominator is constant for the loop
    const Numeric denom = abs(f0) * abs(f0);
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = f_grid[i] * f_grid[i] /
        denom;
    }
}
void lineshape_norm_VVW_dT(Vector&         fac,
                           const Numeric,//   f0,
                           ConstVectorView,// f_grid,
                           const Numeric)
{
    fac = 0.0;
}
void lineshape_norm_VVW_dF(Vector&         fac,
                           const Numeric   f0,
                           ConstVectorView f_grid,
                           const Numeric)
{
    
    const Index nf = f_grid.nelem();
    
    // denominator is constant for the loop
    const Numeric denom = abs(f0) * abs(f0);
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = 2.0 * f_grid[i] /
        denom;
    }
}
void lineshape_norm_VVW_dF0(Vector&         fac,
                            const Numeric   f0,
                            ConstVectorView f_grid,
                            const Numeric)
{
    
    const Index nf = f_grid.nelem();
    
    // denominator is constant for the loop
    const Numeric ddenom = -2 * f0 / abs(f0) / abs(f0) / abs(f0) / abs(f0);
    
    for ( Index i=0; i<nf; ++i )
    {
        fac[i] = f_grid[i] * f_grid[i] * ddenom;
    }
}


//------------------------------------------------------------------------
// Available Lineshapes
//------------------------------------------------------------------------

/*! The lookup data for the different lineshapes. */
namespace global_data {
Array<LineshapeRecord> lineshape_data;
}

void define_lineshape_data()
{
  using global_data::lineshape_data;

  const bool PHASE = true;
  const bool NO_PHASE = false;
  const bool PARTIALS = true;
  const bool NO_PARTIALS = false;

  // Initialize to empty, just in case.
  lineshape_data.resize(0);

  lineshape_data.push_back
    (LineshapeRecord
     ("no_shape",
      "This lineshape does nothing. It only exists, because formally\n"
      "you have to specify a lineshape also for continuum tags.", 
      lineshape_no_shape, NO_PHASE, NO_PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("Lorentz",
      "The Lorentz line shape.",
      lineshape_lorentz, PHASE, NO_PARTIALS));
    
    lineshape_data.push_back
    (LineshapeRecord
    ("Mirrored Lorentz",
     "The mirrored Lorentz line shape.",
     lineshape_mirrored_lorentz, PHASE, NO_PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("Doppler",
      "The Doppler line shape.",
      lineshape_doppler, NO_PHASE, NO_PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("Voigt_Kuntz6",
      "The Voigt line shape. Approximation by Kuntz: Accuracy 2*10-6",
      lineshape_voigt_kuntz6,
      w_x_plus_iy_dT,
      w_x_plus_iy_dF,
      w_x_plus_iy_dF0,
      w_x_plus_iy_dgamma,
      w_x_plus_iy_dH,
      w_x_plus_iy_dDF,
      NO_PHASE, PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("Voigt_Kuntz3",
      "The Voigt line shape. Approximation by Kuntz: Accuracy 2*10-3",
      deprecated, NO_PHASE, NO_PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("Voigt_Kuntz4",
      "The Voigt line shape. Approximation by Kuntz: Accuracy 2*10-4",
      deprecated, NO_PHASE, NO_PARTIALS ));

  lineshape_data.push_back
    (LineshapeRecord
     ("Voigt_Drayson",
      "The Voigt line shape. Approximation by Drayson.",
      lineshape_voigt_drayson, NO_PHASE, NO_PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("CO2_Lorentz",
      "Special line shape for CO2 in the infrared, neglecting Doppler\n"
      "broadening and details of line mixing. The line shape can be\n"
      "expressed as\n"
      "   chi(f,f0) * Lorentz(f,f0) \n"
      "\n"
      "The chi-factor follows Cousin et al. 1985. Broadening by N2 and O2\n"
      "is considered, while self-broadening (CO2-CO2) is neglected."
      "\n"
      "NOTE: Temperature dependency is not yet included. The chi factor is\n"
      "valid for 238 K.",
      lineshape_CO2_lorentz, NO_PHASE, NO_PARTIALS));

  lineshape_data.push_back
    (LineshapeRecord
     ("CO2_Drayson",
      "Special line shape for CO2 in the infrared, neglecting details of\n"
      "line mixing. The line shape can be expressed as\n"
      "   chi(f,f0) * Drayson(f,f0) \n"
      "\n"
      "The chi-factor follows Cousin et al. 1985. Broadening by N2 and O2\n"
      "is considered, while self-broadening (CO2-CO2) is neglected.\n"
      "\n"
      "NOTE: Temperature dependency is not yet included. The chi factor is\n"
      "valid for 238 K.",
      lineshape_CO2_drayson, NO_PHASE, NO_PARTIALS));

    lineshape_data.push_back
    (LineshapeRecord
    ("Faddeeva_Algorithm_916",
    "Voigt and Faraday-Voigt function as per Faddeeva function solution by JPL.\n"
    "Line shape is considered from w(z)=exp(-z^2)*erfc(-iz) where z=v'+ia, and \n"
    "v' is a Doppler weighted freqeuncy parameter and a is a Doppler weighted  \n"
    "pressure parameter.",
     faddeeva_algorithm_916, 
     w_x_plus_iy_dT,
     w_x_plus_iy_dF,
     w_x_plus_iy_dF0,
     w_x_plus_iy_dgamma,
     w_x_plus_iy_dH,
     w_x_plus_iy_dDF,
     PHASE, PARTIALS));

    lineshape_data.push_back
    (LineshapeRecord
    ("Hartmann-Tran",
    "Line shape is considered as described by the Hartmann-Tran profile.",
     hartmann_tran_lineshape, PHASE, NO_PARTIALS));
    
    lineshape_data.push_back
    (LineshapeRecord
    ("Hui_etal_1978",
    "Classic line shape.  Solving the complex error function returns both parts\n"
    "of the refractive index.",
     hui_etal_1978_lineshape, PHASE, NO_PARTIALS));
    
    lineshape_data.push_back
    (LineshapeRecord
    ("O2NonResonant",
    "Special line shape.  Only use for non-resonant O2.",
     lineshape_o2nonresonant, NO_PHASE, NO_PARTIALS));
}

/*! The lookup data for the different normalization factors to the
  lineshapes. */
namespace global_data {
Array<LineshapeNormRecord> lineshape_norm_data;
}

void define_lineshape_norm_data()
{
  using global_data::lineshape_norm_data;

  // Initialize to empty, just in case.
  lineshape_norm_data.resize(0);

  lineshape_norm_data.push_back
    (LineshapeNormRecord
     ("no_norm",
      "No normalization of the lineshape.",
      lineshape_norm_no_norm,lineshape_norm_no_norm_dT,lineshape_norm_no_norm_dF,
      lineshape_norm_no_norm_dF0));

  lineshape_norm_data.push_back
    (LineshapeNormRecord
     ("Rosenkranz_quadratic",
      "Quadratic normalization of the lineshape with (f/f0)^2*h*f0/(2*k*T)/sinh(h*f0/(2*k*T)).",
      lineshape_norm_quadratic_Rosenkranz,lineshape_norm_quadratic_Rosenkranz_dT,
      lineshape_norm_quadratic_Rosenkranz_dF,lineshape_norm_quadratic_Rosenkranz_dF0));

  lineshape_norm_data.push_back
    (LineshapeNormRecord
     ("VVH",
      "Van Vleck Huber normalization of the lineshape with\n"
      "             (f*tanh(h*f/(2*k*T))) / (f0*tanh(h*f0/(2*k*T))).\n"
      "             The denominator is a result of catalogue intensities.",
      lineshape_norm_VVH,lineshape_norm_VVH_dT,lineshape_norm_VVH_dF,lineshape_norm_VVH_dF0));
    
    lineshape_norm_data.push_back
    (LineshapeNormRecord
    ("VVW",
     "Van Vleck Weiskopf normalization of the lineshape with (f*f) / (f0*f0).\n",
     lineshape_norm_VVW,lineshape_norm_VVW_dT,lineshape_norm_VVW_dF,lineshape_norm_VVW_dF0));
}
