/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>

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

  \author Stefan Buehler
  \date   2000-08-21
*/

#include "arts.h"
#include "vecmat.h"
#include "absorption.h"


/*! The Lorentz line shape. This is a quick and dirty implementation.

    \retval ls     The shape function.
    \param  f0     Line center frequency.
    \param  gamma  The pressure broadening parameter.
    \param  sigma  The Doppler broadening parameter. (Not used.)
    \param  f_mono The frequency grid.

    \author Stefan Buehler 16.06.2000 */
void lineshape_lorentz(VECTOR&       ls,
		       Numeric	      f0,
		       Numeric       gamma,
		       Numeric       sigma,
		       const VECTOR& f_mono)
{
  // FIXME: Maybe try if call by reference is faster for f0 and gamma?

  // 1/PI:
  extern const Numeric PI;
  static const Numeric invPI = 1. / PI;

  assert( ls.size() == f_mono.size() );

  for ( size_t i=0; i<f_mono.size(); ++i )
    {
      ls[i] = invPI * gamma / ( pow( f_mono[i]-f0, 2) + pow(gamma,2) );
    }
}



//------------------------------------------------------------------------

// help function for lineshape_voigt_kuntz1
int bfun(Numeric y,Numeric x)
{
      Numeric s;

      s = fabs(x)+y;
      if (s >= 15)
        return(1);
      else if (s >= 5.5)
        return(2); 
      else if (y >= (0.195*fabs(x))-0.176)
        return(3);
      else
        return(4);
}

void lineshape_voigt_kuntz1(VECTOR&       prb,
			    Numeric	  f0,
			    Numeric       gamma,
			    Numeric       sigma,
			    const VECTOR& f_mono)

/*! The Voigt line shape. Kuntz approximation of the Voigt line
  shape, extracted from the Skuld forward model.

    \retval ls     The shape function.
    \param  f0     Line center frequency.
    \param  gamma  The pressure broadening parameter.
    \param  sigma  The Doppler broadening parameter. (Not used.)
    \param  f_mono The frequency grid.

    Original skuld function call and documention:

    void voigt1(int nx, double *x,double y,double *prb,double fak)
    
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


    About 'voigt1' : The program was originally written by M. Kuntz
    in Fortran77 but has been translated into C by Frank Merino and 
    into C++ by Oliver Lemke and Axel von Engeln. fak is set to 1.0.
    

    \author Oliver Lemke and Axel von Engeln
    \date 2000-09-27 */ 
{

  // seems not necessary for Doppler correction
  //    extern const Numeric SQRT_NAT_LOG_2;

      int lauf[5][5], stack[20][4], stackp, bmin, bmax;
      Numeric static a8,b8,c8,d8,e8,f8,g8,h8,o8,p8,q8,r8,s8,t8,a7,b7,c7,d7,e7,f7, 
	g7,h7,o7,p7,q7,r7,s7,t7,a6,b6,c6,d6,e6,a5,b5,c5,d5,e5,a4,b4, 
	c4,d4,a3,b3,c3,d3,a2,a1,b1;

      Numeric yps0,yps1,yps2,yps3,yps4,y2,ym2,x2;
      Numeric c1,b2; // a0,b0;
      int i1,i2,imin,imax,imitte;

      // variables for C++ code FIXME: currently the x array is
      // calculated with each call to this function, it seems better
      // to do this in an intelligent way somewhere before the
      // function call. Note: x needs the Doppler broadening
      // parameter. As well, nx should be defined somewhere else.
      size_t nx;
      Numeric y;
      VECTOR x(f_mono.size());
      Numeric fak;

      // dimension of arrays
      nx = f_mono.size();

      // Ratio of the Lorentz halfwidth to the Doppler halfwidth
      y = gamma / sigma;

      // frequency in units of Doppler 
      // FIXME: see above
      {
	for (i1=0; i1< (int) nx; i1++)
	  {
	    x[i1] = (f_mono[i1] - f0) / sigma;
	  }
      }

      // dummy factor
      fak = 1.0;

      yps0 = -1.0;
      yps1 = -1.0;
      yps2 = -1.0;
      yps3 = -1.0;
      yps4 = -1.0;

      y2 = y*y;

      if ( (y >= 15.0) || (x[0] >= 15.0) || (x[nx-1] <= -15.0) )
      {
        lauf[1][1] = 1;
        lauf[1][2] = nx;
        lauf[1][3] = nx;
        lauf[1][4] = 0;
        goto l7;
      }

      for(i2=1;i2<=4;i2++)
        for(i1=1;i1<=4;i1++)
          lauf[i1][i2] = (i2 % 2)*(nx+1);

      stackp = 1;
      stack[stackp][1] = 1;
      stack[stackp][2] = nx;
      stack[stackp][3] = bfun(y,x[0]);
      stack[stackp][4] = bfun(y,x[nx-1]);

l2:   imin = stack[stackp][1];
      imax = stack[stackp][2];
      bmin = stack[stackp][3];
      bmax = stack[stackp][4];
      if (bmin == bmax)
      {
        if (x[imax-1] < 0.0)
	{
          lauf[bmin][1] = min(imin,lauf[bmin][1]);
          lauf[bmax][2] = max(imax,lauf[bmax][2]);
          stackp = stackp-1;
          goto l3;
	}
        else if (x[imin-1] >= 0.0)
	{
          lauf[bmin][3] = min(imin,lauf[bmin][3]);
          lauf[bmax][4] = max(imax,lauf[bmax][4]);
          stackp = stackp-1;
          goto l3;
	}
      }
      imitte = (imax+imin)/2;
      stack[stackp][1] = imitte+1;
      stack[stackp][2] = imax;
      stack[stackp][3] = bfun(y,x[imitte]);
      stack[stackp][4] = bmax;
      stackp = stackp+1;
      stack[stackp][1] = imin;
      stack[stackp][2] = imitte;
      stack[stackp][3] = bmin;
      stack[stackp][4] = bfun(y,x[imitte-1]);
l3:   if (stackp > 0) goto l2;

/*---- Region 4
--------------------------------------------------------------------*/
      if ( (lauf[4][2] >= lauf[4][1]) || (lauf[4][4] >= lauf[4][3]) )
      {
        if (fabs(y-yps4) > 1e-8)
	{
          yps4 = y;
          a7 = y*(1.16028e9+y2*(-9.86604e8+y2*(4.56662e8+y2* \
          (-1.53575e8+y2*(4.08168e7+y2*(- 9.69463e6+y2*(1.6841e6+y2* \
          (-320772.+y2*(40649.2+y2*(-5860.68+y2*(571.687+y2*(-72.9359 \
          +y2*(2.35944-y2*0.56419)))))))))))));
          b7 = y*(-5.60505e8+y2*(-9.85386e8+y2*(8.06985e8+y2* \
          (-2.91876e8+y2*(8.64829e7+y2*(-7.72359e6+y2*(3.59915e6+y2* \
          (-234417.+y2*(45251.3+y2*(-2269.19+y2*(-234.143+y2* \
          (23.0312-y2*7.33447))))))))))));
          c7 = y*(-6.51523e8+y2*(2.47157e8+y2*(2.94262e8+y2* \
          (-2.04467e8+y2*(2.29302e7+y2*(-2.3818e7+y2*(576054.+y2* \
          (98079.1+y2*(-25338.3+y2*(1097.77+y2* \
          (97.6203-y2*44.0068)))))))))));
          d7 = y*(-2.63894e8+y2*(2.70167e8+y2*(-9.96224e7+y2* \
          (-4.15013e7+y2*(3.83112e7+y2*(2.2404e6+y2*(-303569.+y2* \
          (-66431.2+y2*(8381.97+y2*(228.563-y2*161.358))))))))));
          e7 = y*(-6.31771e7+y2*(1.40677e8+y2*(5.56965e6+y2* \
          (2.46201e7+y2*(468142.+y2*(-1.003e6+y2*(-66212.1+y2* \
          (23507.6+y2*(296.38-y2*403.396)))))))));
          f7 = y*(-1.69846e7+y2*(4.07382e6+y2*(-3.32896e7+y2* \
          (-1.93114e6+y2*(-934717.+y2*(8820.94+y2*(37544.8+y2* \
          (125.591-y2*726.113))))))));
          g7 = y*(-1.23165e6+y2*(7.52883e6+y2*(-900010.+y2*(-186682.+ \
          y2*(79902.5+y2*(37371.9+y2*(-260.198-y2*968.15)))))));
          h7 = y*(-610622.+y2*(86407.6+y2*(153468.+y2*(72520.9+y2* \
          (23137.1+y2*(-571.645-y2*968.15))))));
          o7 = y*(-23586.5+y2*(49883.8+y2*(26538.5+y2*(8073.15+y2* \
          (-575.164-y2*726.113)))));
          p7 = y*(-8009.1+y2*(2198.86+y2*(953.655+y2* \
          (-352.467-y2*403.396))));
          q7 = y*(-622.056+y2*(-271.202+y2*(-134.792-y2*161.358)));
          r7 = y*(-77.0535+y2*(-29.7896-y2*44.0068));
          s7 = y*(-2.92264-y2*7.33447);
          t7 = y*(-0.56419);
          a8 = 1.02827e9+y2*(-1.5599e9+y2*(1.17022e9+y2*(-5.79099e8+y2* \
          (2.11107e8+y2*(-6.11148e7+y2*(1.44647e7+y2*(-2.85721e6+y2* \
          (483737.+y2*(-70946.1+y2*(9504.65+y2*(-955.194+y2*(126.532 \
          +y2*(-3.68288+y2)))))))))))));
          b8 = 1.5599e9+y2*(-2.28855e9+y2*(1.66421e9+y2*(-7.53828e8+y2* \
          (2.89676e8+y2*(-7.01358e7+y2*(1.39465e7+y2*(-2.84954e6+y2* \
          (498334.+y2*(-55600.+y2*(3058.26+y2*(533.254+y2* \
          (-40.5117+y2*14.))))))))))));
          c8 = 1.17022e9+y2*(-1.66421e9+y2*(1.06002e9+y2*(-6.60078e8+y2* \
          (6.33496e7+y2*(-4.60396e7+y2*(1.4841e7+y2*(-1.06352e6+y2* \
          (-217801.+y2*(48153.3+y2*(-1500.17+y2* \
          (-198.876+y2*91)))))))))));
          d8 = 5.79099e8+y2*(-7.53828e8+y2*(6.60078e8+y2*(5.40367e7+y2* \
          (1.99846e8+y2*(-6.87656e6+y2*(-6.89002e6+y2*(280428.+y2* \
          (161461.+y2*(-16493.7+y2*(-567.164+y2*364))))))))));
          e8 = 2.11107e8+y2*(-2.89676e8+y2*(6.33496e7+y2*(-1.99846e8+y2* \
          (-5.01017e7+y2*(-5.25722e6+y2*(1.9547e6+y2*(240373.+y2* \
          (-55582.+y2*(-1012.79+y2*1001)))))))));
          f8 = 6.11148e7+y2*(-7.01358e7+y2*(4.60396e7+y2*(-6.87656e6+y2* \
          (5.25722e6+y2*(3.04316e6+y2*(123052.+y2*(-106663.+y2* \
          (-1093.82+y2*2002))))))));
          g8 = 1.44647e7+y2*(-1.39465e7+y2*(1.4841e7+y2* \
          (6.89002e6+y2*(1.9547e6+y2*(-123052.+y2*(-131337.+y2* \
          (-486.14+y2*3003)))))));
          h8 = 2.85721e6+y2*(-2.84954e6+y2*(1.06352e6+y2*(280428.+y2* \
          (-240373.+y2*(-106663.+y2*(486.14+y2*3432))))));
          o8 = 483737.+y2*(-498334.+y2*(-217801.+y2*(-161461.+y2* \
          (-55582.+y2*(1093.82+y2*3003)))));
          p8 = 70946.1+y2*(-55600.+y2*(-48153.3+y2*(-16493.7+y2* \
          (1012.79+y2*2002))));
          q8 = 9504.65+y2*(-3058.26+y2*(-1500.17+y2*(567.164+y2*1001.)));
          r8 = 955.194+y2*(533.254+y2*(198.876+y2*364));
          s8 = 126.532+y2*(40.5117+y2*91.);
          t8 = 3.68288+y2*14.0;
	}
        ym2 = y*2;
        for(i2=1;i2<=3;i2+=2)
        {
	  for(i1=lauf[4][i2];i1<=lauf[4][i2+1];i1++)
	  {
            x2 = x[i1-1]*x[i1-1];
            prb[i1-1] = fak*(exp(y2-x2)*cos(x[i1-1]*ym2) - \
              (a7+x2*(b7+x2*(c7+x2*(d7+x2*(e7+x2*(f7+x2*(g7+x2*(h7+x2* \
              (o7+x2*(p7+x2*(q7+x2*(r7+x2*(s7+x2*t7))))))))))))) / \
              (a8+x2*(b8+x2*(c8+x2*(d8+x2*(e8+x2*(f8+x2*(g8+x2*(h8+x2* \
              (o8+x2*(p8+x2*(q8+x2*(r8+x2*(s8+x2*(t8+x2)))))))))))))));
          }
        }
      }


/*---- Region 3
--------------------------------------------------------------------*/
      if ( (lauf[3][2] >= lauf[3][1]) || (lauf[3][4] >= lauf[3][3]) )
      { 
        if (fabs(y-yps3) > 1e-8)
	{
          yps3 = y;
          a5 = (272.102+y*(973.778+y*(1629.76+y*(1678.33+y*(1174.8+y* \
          (581.746+y*(204.501+y*(49.5213+y*(7.55895+y*0.564224)))))))));
          b5 = (-60.5644+y*(-2.34403+y*(220.843+y*(336.364+y*(247.198 \
          +y*(100.705+y*(22.6778+y*2.25689)))))));
          c5 = (4.58029+y*(18.546+y*(42.5683+y*(52.8454+y*(22.6798+y* \
          3.38534)))));
          d5 = (-0.128922+y*(1.66203+y*(7.56186+y*2.25689)));
          e5 = (0.000971457+y*0.564224);
          a6 = 272.102+y*(1280.83+y*(2802.87+y*(3764.97+y*(3447.63+y* \
          (2256.98+y*(1074.41+y*(369.199+y*(88.2674+ \
          y*(13.3988+y)))))))));
          b6 = 211.678+y*(902.306+y*(1758.34+y*(2037.31+y*(1549.68+y* \
          (793.427+y*(266.299+y*(53.5952+y*5.)))))));
          c6 = 78.866+y*(308.186+y*(497.302+y*(479.258+y*(269.292+y* \
          (80.3928+y*10.)))));
          d6 = 22.0353+y*(55.0293+y*(92.7568+y*(53.5952+y*10.)));
          e6 = 1.49645+y*(13.3988+y*5.);
	}
	for(i2=1;i2<=3;i2+=2)
	{
	  for(i1=lauf[3][i2];i1<=lauf[3][i2+1];i1++)
	  {
            x2 = x[i1-1]*x[i1-1];
            prb[i1-1] = fak*(a5+x2*(b5+x2*(c5+x2*(d5+x2*e5))))/ \
                     (a6+x2*(b6+x2*(c6+x2*(d6+x2*(e6+x2)))));
	  }
	}
      }


/*---- Region 2
--------------------------------------------------------------------*/
      if ( (lauf[2][2] >= lauf[2][1]) || (lauf[2][4] >= lauf[2][3]) )
      {
        if (fabs(y-yps2) > 1e-8)
	{
          yps2 = y;
          a3 = y*(1.05786+y2*(4.65456+y2*(3.10304+y2*0.56419)));
          b3 = y*(2.962+y2*(0.56419+y2*1.69257));
          c3 = y*(1.69257*y2-2.53885);
          d3 = y*(0.56419);
          a4 = 0.5625+y2*(4.5+y2*(10.5+y2*(6.+y2)));
          b4 = -4.5+y2*(9.+y2*(6.+y2*4.));
          c4 = 10.5+y2*(-6.+y2*6.);
          d4 = -6.+y2*4.0;
	}
	for(i2=1;i2<=3;i2+=2)
	{
	  for(i1=lauf[2][i2];i1<=lauf[2][i2+1];i1++)
	  {
            x2 = x[i1-1]*x[i1-1];
            prb[i1-1] = fak*(a3+x2*(b3+x2*(c3+x2*d3)))/ 
                     (a4+x2*(b4+x2*(c4+x2*(d4+x2))));
	  }
	}
      }


/*---- Region 1
--------------------------------------------------------------------*/
l7:   if ( (lauf[1][2] >= lauf[1][1]) || (lauf[1][4] >= lauf[1][3]) )
      {

        if (fabs(y-yps1) > 1e-8)
	{
          yps1 = y;
          a1 = 0.5641896*y;
          b1 = 0.5+y2;
          a2 = 4*y2;
	}
        c1 = fak*a1;
	for(i2=1;i2<=3;i2+=2)
	{
	  for(i1=lauf[1][i2];i1<=lauf[1][i2+1];i1++)
	  {
            x2 = x[i1-1]*x[i1-1];
            b2 = b1-x2;
            prb[i1-1] = c1*(b1+x2)/(b2*b2+a2*x2);
	  }
	}
      }
}

//---------------------------------------------------------------------------------

/*!  No normalization of the lineshape function.

    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_mono The frequency grid.

    \author Axel von Engeln 30.11.2000 */
void lineshape_norm_no_norm(VECTOR&       fac,
			    Numeric	 f0,
			    const VECTOR& f_mono)
{

  assert( fac.size() == f_mono.size() );

  for ( size_t i=0; i<f_mono.size(); ++i )
    {
      fac[i] = 1.0;
    }
}



/*!  Linear normalization factor of the lineshape function with f/f0.

    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_mono The frequency grid.

    \author Axel von Engeln 30.11.2000 */
void lineshape_norm_linear(VECTOR&       fac,
			   Numeric	 f0,
			   const VECTOR& f_mono)
{
  // FIXME: Maybe try if call by reference is faster for f0?

  assert( fac.size() == f_mono.size() );

  for ( size_t i=0; i<f_mono.size(); ++i )
    {
      fac[i] = f_mono[i] / f0;
    }
}

/*!  Quadratic normalization factor of the lineshape function with (f/f0)^2.

    \retval fac    Normalization factor to the lineshape function.
    \param  f0     Line center frequency.
    \param  f_mono The frequency grid.

    \author Axel von Engeln 30.11.2000 */
void lineshape_norm_quadratic(VECTOR&       fac,
			      Numeric	 f0,
			      const VECTOR& f_mono)
{
  // FIXME: Maybe try if call by reference is faster for f0?

  assert( fac.size() == f_mono.size() );

  // don't do this for the whole loop
  Numeric f0_2 = f0 * f0;

  for ( size_t i=0; i<f_mono.size(); ++i )
    {
      fac[i] = f_mono[i] * f_mono[i] / f0_2;
    }
}




//---------------------------------------------------------------------------------




/*! The lookup data for the different lineshapes. */
ARRAY<LineshapeRecord> lineshape_data;

void define_lineshape_data()
{
  // Initialize to empty, just in case.
  lineshape_data.resize(0);

  lineshape_data.push_back
    (LineshapeRecord
     ("Lorentz",
      "The Lorentz line shape. This is a quick and dirty implementation.",
      -1,
      lineshape_lorentz));

  lineshape_data.push_back
    (LineshapeRecord
     ("Voigt_Kuntz1",
      "The Voigt line shape. Approximation by Kuntz.",
      -1,
      lineshape_voigt_kuntz1));
}

/*! The lookup data for the different normalization factors to the
  lineshapes. */
ARRAY<LineshapeNormRecord> lineshape_norm_data;

void define_lineshape_norm_data()
{
  // Initialize to empty, just in case.
  lineshape_norm_data.resize(0);

  lineshape_norm_data.push_back
    (LineshapeNormRecord
     ("no_norm",
      "No normalization of the lineshape.",
      lineshape_norm_no_norm));

  lineshape_norm_data.push_back
    (LineshapeNormRecord
     ("linear",
      "Linear normalization of the lineshape with f/f0.",
      lineshape_norm_linear));

  lineshape_norm_data.push_back
    (LineshapeNormRecord
     ("quadratic",
      "Quadratic normalization of the lineshape with (f/f0)^2.",
      lineshape_norm_quadratic));
}
