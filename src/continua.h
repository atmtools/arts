/* Copyright (C) 2001 Thomas Kuhn    <tkuhn@uni-bremen.de>
                      Stefan Buehler <sbuehler@uni-bremen.de>

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

/**
   \file   continua.h

   This header file contains all the declarations of the implemented 
   continua and full absorption (lines+continuum) models.

   \author Thomas Kuhn
   \date   2001-11-05
*/

#ifndef continua_h
#define continua_h

#include "matpackI.h"


//////////////////////////////////////////////////////////////////////////// 
// entry function to all continua and full model functions
//////////////////////////////////////////////////////////////////////////// 

void xsec_continuum_tag( MatrixView         xsec,       // calculated x-section
			 const String&      name,       // model name
			 ConstVectorView    parameters, // model 
			 const String&      model,      // model option
			 ConstVectorView    f_mono,     // frequency vector
			 ConstVectorView    p_abs,      // pressure vector
			 ConstVectorView    t_abs,      // temperature vector 
			 ConstVectorView    n2_abs,     // N2 vmr profile
			 ConstVectorView    h2o_abs,    // H2O vmr profile
			 ConstVectorView    vmr );      // species vmr profile

//////////////////////////////////////////////////////////////////////////// 
// check of consistency of all full and continua absorption models
//////////////////////////////////////////////////////////////////////////// 

void check_continuum_model(const String& name);


//////////////////////////////////////////////////////////////////////////// 
// water vapor line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM87H2OAbsModel( MatrixView        xsec,       // calculated x-section
                       const Numeric	 CC,         // continuum scale factor 
		       const Numeric	 CL,         // line strength scale factor
		       const Numeric	 CW,         // line broadening scale factor
                       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void MPM89H2OAbsModel( MatrixView        xsec,       // calculated x-section
                       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void MPM93H2OAbsModel( MatrixView        xsec,
                       const Numeric     CCin,       // continuum scale factor 
		       const Numeric	 CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void PWR98H2OAbsModel( MatrixView        xsec,       // calculated x-section
		       const Numeric	 CCin,       // continuum scale factor 
		       const Numeric     CLin,       // line strength scale factor
		       const Numeric	 CWin,       // line broadening scale factor
		       const String&     model,      // model option
		       ConstVectorView   f_mono,     // frequency vector
		       ConstVectorView   p_abs,      // pressure vector
		       ConstVectorView   t_abs,      // temperature vector
		       ConstVectorView   vmr );      // H2O vmr profile

void CP98H2OAbsModel( MatrixView        xsec,        // calculated x-section
                      const Numeric     CCin,        // continuum scale factor 
		      const Numeric     CLin,        // line strength scale factor
		      const Numeric     CWin,        // line broadening scale factor
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   vmr );       // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// water vapor continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void Pardo_ATM_H2O_ForeignContinuum( MatrixView          xsec,   // calculated x-section
				     const Numeric       Cin,    // model parameter
				     const String&       model,  // model option
				     ConstVectorView     f_mono, // frequency vector
				     ConstVectorView     p_abs,  // pressure vector
				     ConstVectorView     t_abs,  // temperature vector 
				     ConstVectorView     vmr);   // H2O vmr profile

void Standard_H2O_self_continuum( MatrixView        xsec,        // calculated x-section
				  const Numeric     C,           // model parameter
				  const Numeric     x,           // model parameter
				  const String&     model,       // model option
				  ConstVectorView   f_mono,      // frequency vector
				  ConstVectorView   p_abs,       // pressure vector
				  ConstVectorView   t_abs,       // temperature vector 
				  ConstVectorView   vmr);        // H2O vmr profile

void Standard_H2O_foreign_continuum( MatrixView        xsec,     // calculated x-section
				     const Numeric	 C,      // model parameter
				     const Numeric	 x,      // model parameter
				     const String&     model,    // model option
				     ConstVectorView   f_mono,   // frequency vector
				     ConstVectorView   p_abs,    // pressure vector
				     ConstVectorView   t_abs,    // temperature vector 
				     ConstVectorView   vmr);     // H2O vmr profile

void MaTipping_H2O_foreign_continuum( MatrixView        xsec,     // calculated x-section
				     const Numeric	 C,      // model parameter
				     const Numeric	 x,      // model parameter
				     const String&     model,    // model option
				     ConstVectorView   f_mono,   // frequency vector
				     ConstVectorView   p_abs,    // pressure vector
				     ConstVectorView   t_abs,    // temperature vector 
				     ConstVectorView   vmr);     // H2O vmr profile

void MPM93_H2O_continuum( MatrixView        xsec,                // calculated x-section
			  const Numeric	    fcenter,             // model parameter
			  const Numeric	    b1,                  // model parameter
			  const Numeric	    b2,                  // model parameter
			  const Numeric	    b3,                  // model parameter
			  const Numeric	    b4,                  // model parameter
			  const Numeric	    b5,                  // model parameter
			  const Numeric	    b6,                  // model parameter
			  const String&     model,               // model option
			  ConstVectorView   f_mono,              // frequency vector
			  ConstVectorView   p_abs,               // pressure vector
			  ConstVectorView   t_abs,               // temperature vector
			  ConstVectorView   vmr	 );              // H2O vmr profile

void CKD24_H20( MatrixView          xsec,      // calculated x-section
		int                 isf,       // flag if self or foreign cont.
		const Numeric       Cin,       // model scaling factor
		const String&       model,     // model option
		ConstVectorView     f_mono,    // frequency vector
		ConstVectorView     p_abs,     // pressure vector
		ConstVectorView     t_abs,     // temperature vector
		ConstVectorView     vmr,       // H2O vmr profile
                ConstVectorView     n2_abs );  // N2 vmr profile

void CKD_mt_self_h2o( MatrixView          xsec,      // calculated x-section
		      const Numeric       Cin,       // model scaling factor
		      const String&       model,     // model option
		      ConstVectorView     f_mono,    // frequency vector
		      ConstVectorView     p_abs,     // pressure vector
		      ConstVectorView     t_abs,     // temperature vector
		      ConstVectorView     vmr,       // H2O vmr profile
		      ConstVectorView     n2_abs );  // N2 vmr profile

void CKD_mt_foreign_h2o( MatrixView          xsec,      // calculated x-section
			 const Numeric       Cin,       // model scaling factor
			 const String&       model,     // model option
			 ConstVectorView     f_mono,    // frequency vector
			 ConstVectorView     p_abs,     // pressure vector
			 ConstVectorView     t_abs,     // temperature vector
			 ConstVectorView     vmr,       // H2O vmr profile
			 ConstVectorView     n2_abs );  // N2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// oxygen line+continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM85O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM87O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM89O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM92O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile


void MPM93O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile

void PWR93O2AbsModel( MatrixView        xsec,        // calculated x-section
		      const Numeric     CC,          // model parameter
		      const Numeric     CL,          // model parameter
		      const Numeric     CW,          // model parameter
		      const Numeric     CO,          // model parameter
		      const String&     model,       // model option
		      const String&     version,     // model version 1993 or 1988
		      ConstVectorView   f_mono,      // frequency vector
		      ConstVectorView   p_abs,       // pressure vector
		      ConstVectorView   t_abs,       // temperature vector
		      ConstVectorView   h2o_abs,     // H2O vmr profile
		      ConstVectorView   vmr );       // O2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// oxygen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void MPM93_O2_continuum( MatrixView        xsec,             // calculated x-section
			 const Numeric     S0in,             // model parameter
			 const Numeric     G0in,             // model parameter
			 const Numeric     XSOin,            // model parameter
			 const Numeric     XG0in,            // model parameter
			 const String&     model,            // model option
			 ConstVectorView   f_mono,           // frequency vector
			 ConstVectorView   p_abs,            // pressure vector
			 ConstVectorView   t_abs,            // temperature vector
			 ConstVectorView   h2o_abs,          // H2O vmr profile
			 ConstVectorView   vmr	 );          // O2 vmr profile

void Rosenkranz_O2_continuum( MatrixView        xsec,        // calculated x-section
			      const Numeric     S0in,        // model parameter
			      const Numeric     G0in,        // model parameter
			      const Numeric     XSOin,       // model parameter
			      const Numeric     XG0in,       // model parameter
			      const String&     model,       // model option
			      ConstVectorView  	f_mono,      // frequency vector
			      ConstVectorView  	p_abs,       // pressure vector
			      ConstVectorView  	t_abs,       // temperature vector
			      ConstVectorView   h2o_abs,     // H2O vmr profile
			      ConstVectorView   vmr);        // O2 vmr profile

void CKD_mt_CIAfun_o2( MatrixView         xsec,        // calculated x-section
		      const Numeric       Cin,         // scaling factor
		      const String&       model,       // model option
		      ConstVectorView     f_mono,      // frequency vector
		      ConstVectorView     p_abs,       // pressure vector
		      ConstVectorView     t_abs,       // temperature vector
		      ConstVectorView     vmr,         // O2 vmr profile
		      ConstVectorView     h2o_abs );   // H2O vmr profile

void CKD_mt_v0v0_o2( MatrixView          xsec,        // calculated x-section
		     const Numeric       Cin,         // scaling factor
		     const String&       model,       // model option
		     ConstVectorView     f_mono,      // frequency vector
		     ConstVectorView     p_abs,       // pressure vector
		     ConstVectorView     t_abs,       // temperature vector
		     ConstVectorView     vmr,         // O2 vmr profile
		     ConstVectorView     n2_abs,      // N2 vmr profile
		     ConstVectorView     h2o_abs );   // H2O vmr profile

void CKD_mt_v1v0_o2( MatrixView          xsec,        // calculated x-section
		     const Numeric       Cin,         // scaling factor
		     const String&       model,       // model option
		     ConstVectorView     f_mono,      // frequency vector
		     ConstVectorView     p_abs,       // pressure vector
		     ConstVectorView     t_abs,       // temperature vector
		     ConstVectorView     vmr,         // O2 vmr profile
		     ConstVectorView     n2_abs,      // N2 vmr profile
		     ConstVectorView     h2o_abs );   // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// nitrogen continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void CKD_mt_CIArot_n2( MatrixView         xsec,        // calculated x-section
		       const Numeric      Cin,         // scaling factor
		       const String&      model,       // model option
		       ConstVectorView    f_mono,      // frequency vector
		       ConstVectorView    p_abs,       // pressure vector
		       ConstVectorView    t_abs,       // temperature vector
		       ConstVectorView    vmr,         // N2 vmr profile
		       ConstVectorView    h2o_abs );   // H2O vmr profile

void CKD_mt_CIAfun_n2( MatrixView         xsec,        // calculated x-section
		      const Numeric       Cin,         // scaling factor
		      const String&       model,       // model option
		      ConstVectorView     f_mono,      // frequency vector
		      ConstVectorView     p_abs,       // pressure vector
		      ConstVectorView     t_abs,       // temperature vector
		      ConstVectorView     vmr,         // N2 vmr profile
		      ConstVectorView     h2o_abs );   // H2O vmr profile

void BF86_CIA_N2( MatrixView              xsec,        // calculated x-section
		  const Numeric           Cin,         // model parameter
		  const String&           model,       // model option 
		  ConstVectorView         f_mono,      // frequency vector
		  ConstVectorView         p_abs,       // pressure vector
		  ConstVectorView         t_abs,       // temperature vector
		  ConstVectorView         vmr   );     // N2 vmr profile 

void MPM93_N2_continuum( MatrixView       xsec,        // calculated x-section
			 const Numeric    Cin,         // model parameter
			 const Numeric    Gin,         // model parameter
			 const Numeric    xTin,        // model parameter
			 const Numeric    xfin,        // model parameter
			 const String&    model,       // model option
			 ConstVectorView  f_mono,      // frequency vector
			 ConstVectorView  p_abs,       // pressure vector
			 ConstVectorView  t_abs,       // temperature vector
			 ConstVectorView  h2o_abs,     // H2O vmr profile
			 ConstVectorView  vmr	 );    // N2 vmr profile

void Rosenkranz_N2_self_continuum( MatrixView        xsec,        // calculated x-section
				   const Numeric     Cin,         // model parameter
				   const Numeric     xin,         // model parameter
				   const String&     model,       // model option
				   ConstVectorView   f_mono,      // frequency vector
				   ConstVectorView   p_abs,       // pressure vector
				   ConstVectorView   t_abs,       // temperature vector
				   ConstVectorView   vmr );       // N2 vmr profile

void Standard_N2_self_continuum(   MatrixView        xsec,        // calculated x-section
                                   const Numeric     Cin,         // model parameter
                                   const Numeric     xfin,        // model parameter
                                   const Numeric     xtin,        // model parameter
                                   const Numeric     xpin,        // model parameter
				   const String&     model,       // model option
				   ConstVectorView   f_mono,      // frequency vector
				   ConstVectorView   p_abs,       // pressure vector
				   ConstVectorView   t_abs,       // temperature vector
				   ConstVectorView   vmr );       // N2 vmr profile

void Pardo_ATM_N2_dry_continuum( MatrixView        xsec,          // calculated x-section
				  const Numeric     Cin,          // model parameter
				  const String&     model,        // model option
				  ConstVectorView   f_mono,       // frequency vector
				  ConstVectorView   p_abs,        // pressure vector
				  ConstVectorView   t_abs,        // temperature vector
				  ConstVectorView   vmr,	  // N2 vmr profile
				  ConstVectorView   h2o_abs);     // H2O vmr profile

//////////////////////////////////////////////////////////////////////////// 
// carbon dioxide continuum absorption models
//////////////////////////////////////////////////////////////////////////// 

void CKD_mt_co2( MatrixView          xsec,         // calculated x-section
		 const Numeric       Cin,          // scaling factor
		 const String&       model,        // model option
		 ConstVectorView     f_mono,       // frequency vector
		 ConstVectorView     p_abs,        // pressure vector
		 ConstVectorView     t_abs,        // temperature vector
		 ConstVectorView     vmr,	   // CO2 vmr profile
                 ConstVectorView     h2o_abs);     // H2O vmr profile

void Rosenkranz_CO2_self_continuum( MatrixView        xsec,       // calculated x-section
				    const Numeric     C,          // model parameter
				    const Numeric     x,          // model parameter
				    const String&     model,      // model option
				    ConstVectorView   f_mono,     // frequency vector
				    ConstVectorView   p_abs,      // pressure vector
				    ConstVectorView   t_abs,      // temperature vector
				    ConstVectorView   vmr );      // CO2 vmr profile

void Rosenkranz_CO2_foreign_continuum( MatrixView        xsec,    // calculated x-section
				       const Numeric     C,       // model parameter
				       const Numeric     x,       // model parameter
				       const String&     model,   // model option
				       ConstVectorView   f_mono,  // frequency vector
				       ConstVectorView   p_abs,   // pressure vector
				       ConstVectorView   t_abs,   // temperature vector
				       ConstVectorView   n2_abs,  // N2 vmr profile
				       ConstVectorView   vmr );   // CO2 vmr profile

//////////////////////////////////////////////////////////////////////////// 
// water droplet and ice particle absorption (clouds)
//////////////////////////////////////////////////////////////////////////// 

void MPM93WaterDropletAbs( MatrixView        xsec,     // calculated x-section
			   const Numeric     CC,       // model parameter
			   const Numeric     CG,       // model parameter
			   const Numeric     CE,       // model parameter
			   const String&     model,    // model option
			   ConstVectorView   f_mono,   // frequency vector
			   ConstVectorView   p_abs,    // pressure vector
			   ConstVectorView   t_abs,    // temperature vector
			   ConstVectorView   vmr);     // suspended water droplet density vector

void MPM93IceCrystalAbs( MatrixView        xsec,       // calculated x-section
			 const Numeric     CC,         // model parameter
			 const Numeric     CA,         // model parameter
			 const Numeric     CB,         // model parameter
			 const String&     model,      // model option
			 ConstVectorView   f_mono,     // frequency vector
			 ConstVectorView   p_abs,      // pressure vector
			 ConstVectorView   t_abs,      // temperature vector
			 ConstVectorView   vmr	 );    // suspended ice particle density vector, 

void MPM93RainExt( MatrixView        xsec,       // calculated x-section
		   const Numeric     CE,         // model parameter
		   const Numeric     CA,         // model parameter
		   const Numeric     CB,         // model parameter
		   const String&     model,      // model option
		   ConstVectorView   f_mono,     // frequency vector
		   ConstVectorView   p_abs,      // pressure vector
		   ConstVectorView   t_abs,      // temperature vector
		   ConstVectorView   vmr   );    // rain rate vector, 

//////////////////////////////////////////////////////////////////////////// 
// help functions
//////////////////////////////////////////////////////////////////////////// 

Numeric MPMLineShapeFunction( const Numeric gamma,     // line width
			      const Numeric fl,        // line center frequency
			      const Numeric f);        // frequency

Numeric MPMLineShapeO2Function( const Numeric gamma,   // line width
				const Numeric fl,      // line center frequency
				const Numeric f,       // frequency
                                const Numeric delta);  // line coupling

Numeric WVSatPressureLiquidWater(const Numeric t);     // temperature

Numeric WVSatPressureIce(const Numeric t);             // temperature




//////////////////////////////////////////////////////////////////////////// 
// arrays of the CKD H2O, CO2, N2, O2 absorption models
//////////////////////////////////////////////////////////////////////////// 

/*

11 February 2003
 
Release Notes for MT_CKD_1.00
 
Prepared by S. A. Clough,
            AER Inc.,
            131 Harwell Avenue
            Lexington,  MA  02421
            clough@aer.com

This is the initial release of the MT_CKD water vapor continuum and
represents the  first recomputation of the entire self and foreign
broadened continuum since the original model was developed in the
1980s.  This version of the continuum is implemented in the line-by-line
model LBLRTM v7.01 and will be utilized in all related AER Radiative
Transfer models.

further information can be found under
http://www.rtweb.aer.com/continuum_frame.html

Transformation from original F77 code to C/C++ by T. Kuhn, iup Bremen, August 2003 
 */

// additional array fields due to different numbering schemes of F77 and C/C++
const int addF77fields = 1;

// H2O self continuum parameters at T=296 K
// date of last update: 11/18/02
// units of (CM**3/MOL)*1.E-20
const Numeric SL296_ckd_mt_100_v1  =   -20.0;
const Numeric SL296_ckd_mt_100_v2  = 20000.0;
const Numeric SL296_ckd_mt_100_dv  =    10.0;
const int     SL296_ckd_mt_100_npt =  2003;
const double  SL296_ckd_mt_100[SL296_ckd_mt_100_npt+addF77fields] = {
        0.000e0,    1.720e-01,  1.695e-01,
        1.700e-01,  1.695e-01,  1.720e-01,  1.680e-01,  1.687e-01,
        1.624e-01,  1.606e-01,  1.508e-01,  1.447e-01,  1.344e-01,
        1.214e-01,  1.133e-01,  1.009e-01,  9.217e-02,  8.297e-02,
        6.989e-02,  6.513e-02,  5.469e-02,  5.056e-02,  4.417e-02,
        3.779e-02,  3.484e-02,  2.994e-02,  2.720e-02,  2.325e-02,
        2.063e-02,  1.818e-02,  1.592e-02,  1.405e-02,  1.251e-02,
        1.080e-02,  9.647e-03,  8.424e-03,  7.519e-03,  6.555e-03,
        5.880e-03,  5.136e-03,  4.511e-03,  3.989e-03,  3.509e-03,
        3.114e-03,  2.740e-03,  2.446e-03,  2.144e-03,  1.895e-03,
        1.676e-03,  1.486e-03,  1.312e-03,  1.164e-03,  1.031e-03,
        9.129e-04,  8.106e-04,  7.213e-04,  6.400e-04,  5.687e-04,
        5.063e-04,  4.511e-04,  4.029e-04,  3.596e-04,  3.220e-04,
        2.889e-04,  2.597e-04,  2.337e-04,  2.108e-04,  1.907e-04,
        1.728e-04,  1.570e-04,  1.430e-04,  1.305e-04,  1.195e-04,
        1.097e-04,  1.009e-04,  9.307e-05,  8.604e-05,  7.971e-05,
        7.407e-05,  6.896e-05,  6.433e-05,  6.013e-05,  5.631e-05,
        5.283e-05,  4.963e-05,  4.669e-05,  4.398e-05,  4.148e-05,
        3.917e-05,  3.702e-05,  3.502e-05,  3.316e-05,  3.142e-05,
        2.978e-05,  2.825e-05,  2.681e-05,  2.546e-05,  2.419e-05,
        2.299e-05,  2.186e-05,  2.079e-05,  1.979e-05,  1.884e-05,
        1.795e-05,  1.711e-05,  1.633e-05,  1.559e-05,  1.490e-05,
        1.426e-05,  1.367e-05,  1.312e-05,  1.263e-05,  1.218e-05,
        1.178e-05,  1.143e-05,  1.112e-05,  1.088e-05,  1.070e-05,
        1.057e-05,  1.050e-05,  1.051e-05,  1.059e-05,  1.076e-05,
        1.100e-05,  1.133e-05,  1.180e-05,  1.237e-05,  1.308e-05,
        1.393e-05,  1.483e-05,  1.614e-05,  1.758e-05,  1.930e-05,
        2.123e-05,  2.346e-05,  2.647e-05,  2.930e-05,  3.279e-05,
        3.745e-05,  4.152e-05,  4.813e-05,  5.477e-05,  6.203e-05,
        7.331e-05,  8.056e-05,  9.882e-05,  1.050e-04,  1.210e-04,
        1.341e-04,  1.572e-04,  1.698e-04,  1.968e-04,  2.175e-04,
        2.431e-04,  2.735e-04,  2.867e-04,  3.190e-04,  3.371e-04,
        3.554e-04,  3.726e-04,  3.837e-04,  3.878e-04,  3.864e-04,
        3.858e-04,  3.841e-04,  3.852e-04,  3.815e-04,  3.762e-04,
        3.618e-04,  3.579e-04,  3.450e-04,  3.202e-04,  3.018e-04,
        2.785e-04,  2.602e-04,  2.416e-04,  2.097e-04,  1.939e-04,
        1.689e-04,  1.498e-04,  1.308e-04,  1.170e-04,  1.011e-04,
        9.237e-05,  7.909e-05,  7.006e-05,  6.112e-05,  5.401e-05,
        4.914e-05,  4.266e-05,  3.963e-05,  3.316e-05,  3.037e-05,
        2.598e-05,  2.294e-05,  2.066e-05,  1.813e-05,  1.583e-05,
        1.423e-05,  1.247e-05,  1.116e-05,  9.760e-06,  8.596e-06,
        7.720e-06,  6.825e-06,  6.108e-06,  5.366e-06,  4.733e-06,
        4.229e-06,  3.731e-06,  3.346e-06,  2.972e-06,  2.628e-06,
        2.356e-06,  2.102e-06,  1.878e-06,  1.678e-06,  1.507e-06,
        1.348e-06,  1.210e-06,  1.089e-06,  9.806e-07,  8.857e-07,
        8.004e-07,  7.261e-07,  6.599e-07,  6.005e-07,  5.479e-07,
        5.011e-07,  4.595e-07,  4.219e-07,  3.885e-07,  3.583e-07,
        3.314e-07,  3.071e-07,  2.852e-07,  2.654e-07,  2.474e-07,
        2.311e-07,  2.162e-07,  2.026e-07,  1.902e-07,  1.788e-07,
        1.683e-07,  1.587e-07,  1.497e-07,  1.415e-07,  1.338e-07,
        1.266e-07,  1.200e-07,  1.138e-07,  1.080e-07,  1.027e-07,
        9.764e-08,  9.296e-08,  8.862e-08,  8.458e-08,  8.087e-08,
        7.744e-08,  7.429e-08,  7.145e-08,  6.893e-08,  6.664e-08,
        6.468e-08,  6.322e-08,  6.162e-08,  6.070e-08,  5.992e-08,
        5.913e-08,  5.841e-08,  5.796e-08,  5.757e-08,  5.746e-08,
        5.731e-08,  5.679e-08,  5.577e-08,  5.671e-08,  5.656e-08,
        5.594e-08,  5.593e-08,  5.602e-08,  5.620e-08,  5.693e-08,
        5.725e-08,  5.858e-08,  6.037e-08,  6.249e-08,  6.535e-08,
        6.899e-08,  7.356e-08,  7.918e-08,  8.618e-08,  9.385e-08,
        1.039e-07,  1.158e-07,  1.290e-07,  1.437e-07,  1.650e-07,
        1.871e-07,  2.121e-07,  2.427e-07,  2.773e-07,  3.247e-07,
        3.677e-07,  4.037e-07,  4.776e-07,  5.101e-07,  6.214e-07,
        6.936e-07,  7.581e-07,  8.486e-07,  9.355e-07,  9.942e-07,
        1.063e-06,  1.123e-06,  1.191e-06,  1.215e-06,  1.247e-06,
        1.260e-06,  1.271e-06,  1.284e-06,  1.317e-06,  1.323e-06,
        1.349e-06,  1.353e-06,  1.362e-06,  1.344e-06,  1.329e-06,
        1.336e-06,  1.327e-06,  1.325e-06,  1.359e-06,  1.374e-06,
        1.415e-06,  1.462e-06,  1.526e-06,  1.619e-06,  1.735e-06,
        1.863e-06,  2.034e-06,  2.265e-06,  2.482e-06,  2.756e-06,
        3.103e-06,  3.466e-06,  3.832e-06,  4.378e-06,  4.913e-06,
        5.651e-06,  6.311e-06,  7.169e-06,  8.057e-06,  9.253e-06,
        1.047e-05,  1.212e-05,  1.360e-05,  1.569e-05,  1.776e-05,
        2.020e-05,  2.281e-05,  2.683e-05,  2.994e-05,  3.488e-05,
        3.896e-05,  4.499e-05,  5.175e-05,  6.035e-05,  6.340e-05,
        7.281e-05,  7.923e-05,  8.348e-05,  9.631e-05,  1.044e-04,
        1.102e-04,  1.176e-04,  1.244e-04,  1.283e-04,  1.326e-04,
        1.400e-04,  1.395e-04,  1.387e-04,  1.363e-04,  1.314e-04,
        1.241e-04,  1.228e-04,  1.148e-04,  1.086e-04,  1.018e-04,
        8.890e-05,  8.316e-05,  7.292e-05,  6.452e-05,  5.625e-05,
        5.045e-05,  4.380e-05,  3.762e-05,  3.290e-05,  2.836e-05,
        2.485e-05,  2.168e-05,  1.895e-05,  1.659e-05,  1.453e-05,
        1.282e-05,  1.132e-05,  1.001e-05,  8.836e-06,  7.804e-06,
        6.922e-06,  6.116e-06,  5.429e-06,  4.824e-06,  4.278e-06,
        3.788e-06,  3.371e-06,  2.985e-06,  2.649e-06,  2.357e-06,
        2.090e-06,  1.858e-06,  1.647e-06,  1.462e-06,  1.299e-06,
        1.155e-06,  1.028e-06,  9.142e-07,  8.132e-07,  7.246e-07,
        6.451e-07,  5.764e-07,  5.151e-07,  4.603e-07,  4.121e-07,
        3.694e-07,  3.318e-07,  2.985e-07,  2.690e-07,  2.428e-07,
        2.197e-07,  1.992e-07,  1.810e-07,  1.649e-07,  1.506e-07,
        1.378e-07,  1.265e-07,  1.163e-07,  1.073e-07,  9.918e-08,
        9.191e-08,  8.538e-08,  7.949e-08,  7.419e-08,  6.940e-08,
        6.508e-08,  6.114e-08,  5.761e-08,  5.437e-08,  5.146e-08,
        4.890e-08,  4.636e-08,  4.406e-08,  4.201e-08,  4.015e-08,
        3.840e-08,  3.661e-08,  3.510e-08,  3.377e-08,  3.242e-08,
        3.130e-08,  3.015e-08,  2.918e-08,  2.830e-08,  2.758e-08,
        2.707e-08,  2.656e-08,  2.619e-08,  2.609e-08,  2.615e-08,
        2.630e-08,  2.675e-08,  2.745e-08,  2.842e-08,  2.966e-08,
        3.125e-08,  3.318e-08,  3.565e-08,  3.850e-08,  4.191e-08,
        4.590e-08,  5.059e-08,  5.607e-08,  6.239e-08,  6.958e-08,
        7.796e-08,  8.773e-08,  9.880e-08,  1.114e-07,  1.258e-07,
        1.422e-07,  1.610e-07,  1.822e-07,  2.060e-07,  2.337e-07,
        2.645e-07,  2.996e-07,  3.393e-07,  3.843e-07,  4.363e-07,
        4.935e-07,  5.607e-07,  6.363e-07,  7.242e-07,  8.230e-07,
        9.411e-07,  1.071e-06,  1.232e-06,  1.402e-06,  1.600e-06,
        1.820e-06,  2.128e-06,  2.386e-06,  2.781e-06,  3.242e-06,
        3.653e-06,  4.323e-06,  4.747e-06,  5.321e-06,  5.919e-06,
        6.681e-06,  7.101e-06,  7.983e-06,  8.342e-06,  8.741e-06,
        9.431e-06,  9.952e-06,  1.026e-05,  1.055e-05,  1.095e-05,
        1.095e-05,  1.087e-05,  1.056e-05,  1.026e-05,  9.715e-06,
        9.252e-06,  8.452e-06,  7.958e-06,  7.268e-06,  6.295e-06,
        6.003e-06,  5.000e-06,  4.591e-06,  3.983e-06,  3.479e-06,
        3.058e-06,  2.667e-06,  2.293e-06,  1.995e-06,  1.747e-06,
        1.517e-06,  1.335e-06,  1.165e-06,  1.028e-06,  9.007e-07,
        7.956e-07,  7.015e-07,  6.192e-07,  5.491e-07,  4.859e-07,
        4.297e-07,  3.799e-07,  3.380e-07,  3.002e-07,  2.659e-07,
        2.366e-07,  2.103e-07,  1.861e-07,  1.655e-07,  1.469e-07,
        1.309e-07,  1.162e-07,  1.032e-07,  9.198e-08,  8.181e-08,
        7.294e-08,  6.516e-08,  5.787e-08,  5.163e-08,  4.612e-08,
        4.119e-08,  3.695e-08,  3.308e-08,  2.976e-08,  2.670e-08,
        2.407e-08,  2.171e-08,  1.965e-08,  1.780e-08,  1.617e-08,
        1.470e-08,  1.341e-08,  1.227e-08,  1.125e-08,  1.033e-08,
        9.524e-09,  8.797e-09,  8.162e-09,  7.565e-09,  7.040e-09,
        6.560e-09,  6.129e-09,  5.733e-09,  5.376e-09,  5.043e-09,
        4.750e-09,  4.466e-09,  4.211e-09,  3.977e-09,  3.759e-09,
        3.558e-09,  3.373e-09,  3.201e-09,  3.043e-09,  2.895e-09,
        2.760e-09,  2.635e-09,  2.518e-09,  2.411e-09,  2.314e-09,
        2.230e-09,  2.151e-09,  2.087e-09,  2.035e-09,  1.988e-09,
        1.946e-09,  1.927e-09,  1.916e-09,  1.916e-09,  1.933e-09,
        1.966e-09,  2.018e-09,  2.090e-09,  2.182e-09,  2.299e-09,
        2.442e-09,  2.623e-09,  2.832e-09,  3.079e-09,  3.368e-09,
        3.714e-09,  4.104e-09,  4.567e-09,  5.091e-09,  5.701e-09,
        6.398e-09,  7.194e-09,  8.127e-09,  9.141e-09,  1.035e-08,
        1.177e-08,  1.338e-08,  1.508e-08,  1.711e-08,  1.955e-08,
        2.216e-08,  2.534e-08,  2.871e-08,  3.291e-08,  3.711e-08,
        4.285e-08,  4.868e-08,  5.509e-08,  6.276e-08,  7.262e-08,
        8.252e-08,  9.400e-08,  1.064e-07,  1.247e-07,  1.411e-07,
        1.626e-07,  1.827e-07,  2.044e-07,  2.284e-07,  2.452e-07,
        2.854e-07,  3.026e-07,  3.278e-07,  3.474e-07,  3.693e-07,
        3.930e-07,  4.104e-07,  4.220e-07,  4.439e-07,  4.545e-07,
        4.778e-07,  4.812e-07,  5.018e-07,  4.899e-07,  5.075e-07,
        5.073e-07,  5.171e-07,  5.131e-07,  5.250e-07,  5.617e-07,
        5.846e-07,  6.239e-07,  6.696e-07,  7.398e-07,  8.073e-07,
        9.150e-07,  1.009e-06,  1.116e-06,  1.264e-06,  1.439e-06,
        1.644e-06,  1.856e-06,  2.147e-06,  2.317e-06,  2.713e-06,
        2.882e-06,  2.990e-06,  3.489e-06,  3.581e-06,  4.033e-06,
        4.260e-06,  4.543e-06,  4.840e-06,  4.826e-06,  5.013e-06,
        5.252e-06,  5.277e-06,  5.306e-06,  5.236e-06,  5.123e-06,
        5.171e-06,  4.843e-06,  4.615e-06,  4.385e-06,  3.970e-06,
        3.693e-06,  3.231e-06,  2.915e-06,  2.495e-06,  2.144e-06,
        1.910e-06,  1.639e-06,  1.417e-06,  1.226e-06,  1.065e-06,
        9.290e-07,  8.142e-07,  7.161e-07,  6.318e-07,  5.581e-07,
        4.943e-07,  4.376e-07,  3.884e-07,  3.449e-07,  3.060e-07,
        2.712e-07,  2.412e-07,  2.139e-07,  1.903e-07,  1.689e-07,
        1.499e-07,  1.331e-07,  1.183e-07,  1.050e-07,  9.362e-08,
        8.306e-08,  7.403e-08,  6.578e-08,  5.853e-08,  5.216e-08,
        4.632e-08,  4.127e-08,  3.678e-08,  3.279e-08,  2.923e-08,
        2.612e-08,  2.339e-08,  2.094e-08,  1.877e-08,  1.686e-08,
        1.516e-08,  1.366e-08,  1.234e-08,  1.114e-08,  1.012e-08,
        9.182e-09,  8.362e-09,  7.634e-09,  6.981e-09,  6.406e-09,
        5.888e-09,  5.428e-09,  5.021e-09,  4.650e-09,  4.326e-09,
        4.033e-09,  3.770e-09,  3.536e-09,  3.327e-09,  3.141e-09,
        2.974e-09,  2.825e-09,  2.697e-09,  2.584e-09,  2.488e-09,
        2.406e-09,  2.340e-09,  2.292e-09,  2.259e-09,  2.244e-09,
        2.243e-09,  2.272e-09,  2.310e-09,  2.378e-09,  2.454e-09,
        2.618e-09,  2.672e-09,  2.831e-09,  3.050e-09,  3.225e-09,
        3.425e-09,  3.677e-09,  3.968e-09,  4.221e-09,  4.639e-09,
        4.960e-09,  5.359e-09,  5.649e-09,  6.230e-09,  6.716e-09,
        7.218e-09,  7.746e-09,  7.988e-09,  8.627e-09,  8.999e-09,
        9.442e-09,  9.820e-09,  1.015e-08,  1.060e-08,  1.079e-08,
        1.109e-08,  1.137e-08,  1.186e-08,  1.180e-08,  1.187e-08,
        1.194e-08,  1.192e-08,  1.224e-08,  1.245e-08,  1.246e-08,
        1.318e-08,  1.377e-08,  1.471e-08,  1.582e-08,  1.713e-08,
        1.853e-08,  2.063e-08,  2.270e-08,  2.567e-08,  2.891e-08,
        3.264e-08,  3.744e-08,  4.286e-08,  4.915e-08,  5.623e-08,
        6.336e-08,  7.293e-08,  8.309e-08,  9.319e-08,  1.091e-07,
        1.243e-07,  1.348e-07,  1.449e-07,  1.620e-07,  1.846e-07,
        1.937e-07,  2.040e-07,  2.179e-07,  2.298e-07,  2.433e-07,
        2.439e-07,  2.464e-07,  2.611e-07,  2.617e-07,  2.582e-07,
        2.453e-07,  2.401e-07,  2.349e-07,  2.203e-07,  2.066e-07,
        1.939e-07,  1.780e-07,  1.558e-07,  1.391e-07,  1.203e-07,
        1.048e-07,  9.464e-08,  8.306e-08,  7.239e-08,  6.317e-08,
        5.520e-08,  4.847e-08,  4.282e-08,  3.796e-08,  3.377e-08,
        2.996e-08,  2.678e-08,  2.400e-08,  2.134e-08,  1.904e-08,
        1.705e-08,  1.523e-08,  1.350e-08,  1.204e-08,  1.070e-08,
        9.408e-09,  8.476e-09,  7.470e-09,  6.679e-09,  5.929e-09,
        5.267e-09,  4.711e-09,  4.172e-09,  3.761e-09,  3.288e-09,
        2.929e-09,  2.609e-09,  2.315e-09,  2.042e-09,  1.844e-09,
        1.640e-09,  1.470e-09,  1.310e-09,  1.176e-09,  1.049e-09,
        9.377e-10,  8.462e-10,  7.616e-10,  6.854e-10,  6.191e-10,
        5.596e-10,  5.078e-10,  4.611e-10,  4.197e-10,  3.830e-10,
        3.505e-10,  3.215e-10,  2.956e-10,  2.726e-10,  2.521e-10,
        2.338e-10,  2.173e-10,  2.026e-10,  1.895e-10,  1.777e-10,
        1.672e-10,  1.579e-10,  1.496e-10,  1.423e-10,  1.358e-10,
        1.302e-10,  1.254e-10,  1.216e-10,  1.187e-10,  1.163e-10,
        1.147e-10,  1.145e-10,  1.150e-10,  1.170e-10,  1.192e-10,
        1.250e-10,  1.298e-10,  1.345e-10,  1.405e-10,  1.538e-10,
        1.648e-10,  1.721e-10,  1.872e-10,  1.968e-10,  2.089e-10,
        2.172e-10,  2.317e-10,  2.389e-10,  2.503e-10,  2.585e-10,
        2.686e-10,  2.800e-10,  2.895e-10,  3.019e-10,  3.037e-10,
        3.076e-10,  3.146e-10,  3.198e-10,  3.332e-10,  3.397e-10,
        3.540e-10,  3.667e-10,  3.895e-10,  4.071e-10,  4.565e-10,
        4.983e-10,  5.439e-10,  5.968e-10,  6.676e-10,  7.456e-10,
        8.405e-10,  9.478e-10,  1.064e-09,  1.218e-09,  1.386e-09,
        1.581e-09,  1.787e-09,  2.032e-09,  2.347e-09,  2.677e-09,
        3.008e-09,  3.544e-09,  4.056e-09,  4.687e-09,  5.331e-09,
        6.227e-09,  6.854e-09,  8.139e-09,  8.945e-09,  9.865e-09,
        1.125e-08,  1.178e-08,  1.364e-08,  1.436e-08,  1.540e-08,
        1.672e-08,  1.793e-08,  1.906e-08,  2.036e-08,  2.144e-08,
        2.292e-08,  2.371e-08,  2.493e-08,  2.606e-08,  2.706e-08,
        2.866e-08,  3.036e-08,  3.136e-08,  3.405e-08,  3.665e-08,
        3.837e-08,  4.229e-08,  4.748e-08,  5.320e-08,  5.763e-08,
        6.677e-08,  7.216e-08,  7.716e-08,  8.958e-08,  9.419e-08,
        1.036e-07,  1.108e-07,  1.189e-07,  1.246e-07,  1.348e-07,
        1.310e-07,  1.361e-07,  1.364e-07,  1.363e-07,  1.343e-07,
        1.293e-07,  1.254e-07,  1.235e-07,  1.158e-07,  1.107e-07,
        9.961e-08,  9.011e-08,  7.910e-08,  6.916e-08,  6.338e-08,
        5.564e-08,  4.827e-08,  4.198e-08,  3.695e-08,  3.276e-08,
        2.929e-08,  2.633e-08,  2.391e-08,  2.192e-08,  2.021e-08,
        1.890e-08,  1.772e-08,  1.667e-08,  1.603e-08,  1.547e-08,
        1.537e-08,  1.492e-08,  1.515e-08,  1.479e-08,  1.450e-08,
        1.513e-08,  1.495e-08,  1.529e-08,  1.565e-08,  1.564e-08,
        1.553e-08,  1.569e-08,  1.584e-08,  1.570e-08,  1.538e-08,
        1.513e-08,  1.472e-08,  1.425e-08,  1.349e-08,  1.328e-08,
        1.249e-08,  1.170e-08,  1.077e-08,  9.514e-09,  8.614e-09,
        7.460e-09,  6.621e-09,  5.775e-09,  5.006e-09,  4.308e-09,
        3.747e-09,  3.240e-09,  2.840e-09,  2.481e-09,  2.184e-09,
        1.923e-09,  1.710e-09,  1.504e-09,  1.334e-09,  1.187e-09,
        1.053e-09,  9.367e-10,  8.306e-10,  7.419e-10,  6.630e-10,
        5.918e-10,  5.277e-10,  4.717e-10,  4.222e-10,  3.783e-10,
        3.390e-10,  3.036e-10,  2.729e-10,  2.455e-10,  2.211e-10,
        1.995e-10,  1.804e-10,  1.635e-10,  1.485e-10,  1.355e-10,
        1.240e-10,  1.139e-10,  1.051e-10,  9.757e-11,  9.114e-11,
        8.577e-11,  8.139e-11,  7.792e-11,  7.520e-11,  7.390e-11,
        7.311e-11,  7.277e-11,  7.482e-11,  7.698e-11,  8.162e-11,
        8.517e-11,  8.968e-11,  9.905e-11,  1.075e-10,  1.187e-10,
        1.291e-10,  1.426e-10,  1.573e-10,  1.734e-10,  1.905e-10,
        2.097e-10,  2.280e-10,  2.473e-10,  2.718e-10,  2.922e-10,
        3.128e-10,  3.361e-10,  3.641e-10,  3.910e-10,  4.196e-10,
        4.501e-10,  4.932e-10,  5.258e-10,  5.755e-10,  6.253e-10,
        6.664e-10,  7.344e-10,  7.985e-10,  8.877e-10,  1.005e-09,
        1.118e-09,  1.251e-09,  1.428e-09,  1.610e-09,  1.888e-09,
        2.077e-09,  2.331e-09,  2.751e-09,  3.061e-09,  3.522e-09,
        3.805e-09,  4.181e-09,  4.575e-09,  5.167e-09,  5.634e-09,
        6.007e-09,  6.501e-09,  6.829e-09,  7.211e-09,  7.262e-09,
        7.696e-09,  7.832e-09,  7.799e-09,  7.651e-09,  7.304e-09,
        7.150e-09,  6.977e-09,  6.603e-09,  6.209e-09,  5.690e-09,
        5.432e-09,  4.764e-09,  4.189e-09,  3.640e-09,  3.203e-09,
        2.848e-09,  2.510e-09,  2.194e-09,  1.946e-09,  1.750e-09,
        1.567e-09,  1.426e-09,  1.302e-09,  1.197e-09,  1.109e-09,
        1.035e-09,  9.719e-10,  9.207e-10,  8.957e-10,  8.578e-10,
        8.262e-10,  8.117e-10,  7.987e-10,  7.875e-10,  7.741e-10,
        7.762e-10,  7.537e-10,  7.424e-10,  7.474e-10,  7.294e-10,
        7.216e-10,  7.233e-10,  7.075e-10,  6.892e-10,  6.618e-10,
        6.314e-10,  6.208e-10,  5.689e-10,  5.550e-10,  4.984e-10,
        4.600e-10,  4.078e-10,  3.879e-10,  3.459e-10,  2.982e-10,
        2.626e-10,  2.329e-10,  1.988e-10,  1.735e-10,  1.487e-10,
        1.297e-10,  1.133e-10,  9.943e-11,  8.736e-11,  7.726e-11,
        6.836e-11,  6.053e-11,  5.384e-11,  4.789e-11,  4.267e-11,
        3.804e-11,  3.398e-11,  3.034e-11,  2.710e-11,  2.425e-11,
        2.173e-11,  1.950e-11,  1.752e-11,  1.574e-11,  1.418e-11,
        1.278e-11,  1.154e-11,  1.044e-11,  9.463e-12,  8.602e-12,
        7.841e-12,  7.171e-12,  6.584e-12,  6.073e-12,  5.631e-12,
        5.254e-12,  4.937e-12,  4.679e-12,  4.476e-12,  4.328e-12,
        4.233e-12,  4.194e-12,  4.211e-12,  4.286e-12,  4.424e-12,
        4.628e-12,  4.906e-12,  5.262e-12,  5.708e-12,  6.254e-12,
        6.914e-12,  7.714e-12,  8.677e-12,  9.747e-12,  1.101e-11,
        1.256e-11,  1.409e-11,  1.597e-11,  1.807e-11,  2.034e-11,
        2.316e-11,  2.622e-11,  2.962e-11,  3.369e-11,  3.819e-11,
        4.329e-11,  4.932e-11,  5.589e-11,  6.364e-11,  7.284e-11,
        8.236e-11,  9.447e-11,  1.078e-10,  1.229e-10,  1.417e-10,
        1.614e-10,  1.843e-10,  2.107e-10,  2.406e-10,  2.728e-10,
        3.195e-10,  3.595e-10,  4.153e-10,  4.736e-10,  5.410e-10,
        6.088e-10,  6.769e-10,  7.691e-10,  8.545e-10,  9.621e-10,
        1.047e-09,  1.161e-09,  1.296e-09,  1.424e-09,  1.576e-09,
        1.739e-09,  1.893e-09,  2.080e-09,  2.336e-09,  2.604e-09,
        2.760e-09,  3.001e-09,  3.365e-09,  3.550e-09,  3.895e-09,
        4.183e-09,  4.614e-09,  4.846e-09,  5.068e-09,  5.427e-09,
        5.541e-09,  5.864e-09,  5.997e-09,  5.997e-09,  6.061e-09,
        5.944e-09,  5.855e-09,  5.661e-09,  5.523e-09,  5.374e-09,
        4.940e-09,  4.688e-09,  4.170e-09,  3.913e-09,  3.423e-09,
        2.997e-09,  2.598e-09,  2.253e-09,  1.946e-09,  1.710e-09,
        1.507e-09,  1.336e-09,  1.190e-09,  1.068e-09,  9.623e-10,
        8.772e-10,  8.007e-10,  7.420e-10,  6.884e-10,  6.483e-10,
        6.162e-10,  5.922e-10,  5.688e-10,  5.654e-10,  5.637e-10,
        5.701e-10,  5.781e-10,  5.874e-10,  6.268e-10,  6.357e-10,
        6.525e-10,  7.137e-10,  7.441e-10,  8.024e-10,  8.485e-10,
        9.143e-10,  9.536e-10,  9.717e-10,  1.018e-09,  1.042e-09,
        1.054e-09,  1.092e-09,  1.079e-09,  1.064e-09,  1.043e-09,
        1.020e-09,  9.687e-10,  9.273e-10,  9.208e-10,  9.068e-10,
        7.687e-10,  7.385e-10,  6.595e-10,  5.870e-10,  5.144e-10,
        4.417e-10,  3.804e-10,  3.301e-10,  2.866e-10,  2.509e-10,
        2.202e-10,  1.947e-10,  1.719e-10,  1.525e-10,  1.361e-10,
        1.210e-10,  1.084e-10,  9.800e-11,  8.801e-11,  7.954e-11,
        7.124e-11,  6.335e-11,  5.760e-11,  5.132e-11,  4.601e-11,
        4.096e-11,  3.657e-11,  3.250e-11,  2.909e-11,  2.587e-11,
        2.297e-11,  2.050e-11,  1.828e-11,  1.632e-11,  1.462e-11,
        1.314e-11,  1.185e-11,  1.073e-11,  9.760e-12,  8.922e-12,
        8.206e-12,  7.602e-12,  7.100e-12,  6.694e-12,  6.378e-12,
        6.149e-12,  6.004e-12,  5.941e-12,  5.962e-12,  6.069e-12,
        6.265e-12,  6.551e-12,  6.935e-12,  7.457e-12,  8.074e-12,
        8.811e-12,  9.852e-12,  1.086e-11,  1.207e-11,  1.361e-11,
        1.553e-11,  1.737e-11,  1.930e-11,  2.175e-11,  2.410e-11,
        2.706e-11,  3.023e-11,  3.313e-11,  3.657e-11,  4.118e-11,
        4.569e-11,  5.025e-11,  5.660e-11,  6.231e-11,  6.881e-11,
        7.996e-11,  8.526e-11,  9.694e-11,  1.106e-10,  1.222e-10,
        1.355e-10,  1.525e-10,  1.775e-10,  1.924e-10,  2.181e-10,
        2.379e-10,  2.662e-10,  2.907e-10,  3.154e-10,  3.366e-10,
        3.579e-10,  3.858e-10,  4.046e-10,  4.196e-10,  4.166e-10,
        4.457e-10,  4.466e-10,  4.404e-10,  4.337e-10,  4.150e-10,
        4.083e-10,  3.910e-10,  3.723e-10,  3.514e-10,  3.303e-10,
        2.847e-10,  2.546e-10,  2.230e-10,  1.994e-10,  1.733e-10,
        1.488e-10,  1.297e-10,  1.144e-10,  1.004e-10,  8.741e-11,
        7.928e-11,  7.034e-11,  6.323e-11,  5.754e-11,  5.250e-11,
        4.850e-11,  4.502e-11,  4.286e-11,  4.028e-11,  3.899e-11,
        3.824e-11,  3.761e-11,  3.804e-11,  3.839e-11,  3.845e-11,
        4.244e-11,  4.382e-11,  4.582e-11,  4.847e-11,  5.209e-11,
        5.384e-11,  5.887e-11,  6.371e-11,  6.737e-11,  7.168e-11,
        7.415e-11,  7.827e-11,  8.037e-11,  8.120e-11,  8.071e-11,
        8.008e-11,  7.851e-11,  7.544e-11,  7.377e-11,  7.173e-11,
        6.801e-11,  6.267e-11,  5.727e-11,  5.288e-11,  4.853e-11,
        4.082e-11,  3.645e-11,  3.136e-11,  2.672e-11,  2.304e-11,
        1.986e-11,  1.725e-11,  1.503e-11,  1.315e-11,  1.153e-11,
        1.014e-11,  8.942e-12,  7.901e-12,  6.993e-12,  6.199e-12,
        5.502e-12,  4.890e-12,  4.351e-12,  3.878e-12,  3.461e-12,
        3.094e-12,  2.771e-12,  2.488e-12,  2.241e-12,  2.025e-12,
        1.838e-12,  1.677e-12,  1.541e-12,  1.427e-12,  1.335e-12,
        1.262e-12,  1.209e-12,  1.176e-12,  1.161e-12,  1.165e-12,
        1.189e-12,  1.234e-12,  1.300e-12,  1.389e-12,  1.503e-12,
        1.644e-12,  1.814e-12,  2.017e-12,  2.255e-12,  2.534e-12,
        2.858e-12,  3.231e-12,  3.661e-12,  4.153e-12,  4.717e-12,
        5.360e-12,  6.094e-12,  6.930e-12,  7.882e-12,  8.966e-12,
        1.020e-11,  1.162e-11,  1.324e-11,  1.510e-11,  1.720e-11,
        1.965e-11,  2.237e-11,  2.560e-11,  2.927e-11,  3.371e-11,
        3.842e-11,  4.429e-11,  5.139e-11,  5.798e-11,  6.697e-11,
        7.626e-11,  8.647e-11,  1.022e-10,  1.136e-10,  1.300e-10,
        1.481e-10,  1.672e-10,  1.871e-10,  2.126e-10,  2.357e-10,
        2.583e-10,  2.997e-10,  3.289e-10,  3.702e-10,  4.012e-10,
        4.319e-10,  4.527e-10,  5.001e-10,  5.448e-10,  5.611e-10,
        5.760e-10,  5.965e-10,  6.079e-10,  6.207e-10,  6.276e-10,
        6.222e-10,  6.137e-10,  6.000e-10,  5.814e-10,  5.393e-10,
        5.350e-10,  4.947e-10,  4.629e-10,  4.117e-10,  3.712e-10,
        3.372e-10,  2.923e-10,  2.550e-10,  2.232e-10,  1.929e-10,
        1.679e-10,  1.460e-10,  1.289e-10,  1.130e-10,  9.953e-11,
        8.763e-11,  7.760e-11,  6.900e-11,  6.160e-11,  5.525e-11,
        4.958e-11,  4.489e-11,  4.072e-11,  3.728e-11,  3.438e-11,
        3.205e-11,  3.006e-11,  2.848e-11,  2.766e-11,  2.688e-11,
        2.664e-11,  2.670e-11,  2.696e-11,  2.786e-11,  2.861e-11,
        3.009e-11,  3.178e-11,  3.389e-11,  3.587e-11,  3.819e-11,
        4.054e-11,  4.417e-11,  4.703e-11,  5.137e-11,  5.460e-11,
        6.055e-11,  6.333e-11,  6.773e-11,  7.219e-11,  7.717e-11,
        8.131e-11,  8.491e-11,  8.574e-11,  9.010e-11,  9.017e-11,
        8.999e-11,  8.959e-11,  8.838e-11,  8.579e-11,  8.162e-11,
        8.098e-11,  7.472e-11,  7.108e-11,  6.559e-11,  5.994e-11,
        5.172e-11,  4.424e-11,  3.951e-11,  3.340e-11,  2.902e-11,
        2.541e-11,  2.215e-11,  1.945e-11,  1.716e-11,  1.503e-11,
        1.339e-11,  1.185e-11,  1.050e-11,  9.336e-12,  8.307e-12,
        7.312e-12,  6.550e-12,  5.836e-12,  5.178e-12,  4.600e-12,
        4.086e-12,  3.639e-12,  3.247e-12,  2.904e-12,  2.604e-12,
        2.341e-12,  2.112e-12,  1.914e-12,  1.744e-12,  1.598e-12,
        1.476e-12,  1.374e-12,  1.293e-12,  1.230e-12,  1.185e-12,
        1.158e-12,  1.147e-12,  1.154e-12,  1.177e-12,  1.219e-12,
        1.280e-12,  1.360e-12,  1.463e-12,  1.591e-12,  1.750e-12,
        1.940e-12,  2.156e-12,  2.430e-12,  2.748e-12,  3.052e-12,
        3.533e-12,  3.967e-12,  4.471e-12,  5.041e-12,  5.860e-12,
        6.664e-12,  7.522e-12,  8.342e-12,  9.412e-12,  1.072e-11,
        1.213e-11,  1.343e-11,  1.496e-11,  1.664e-11,  1.822e-11,
        2.029e-11,  2.233e-11,  2.457e-11,  2.709e-11,  2.928e-11,
        3.115e-11,  3.356e-11,  3.592e-11,  3.818e-11,  3.936e-11,
        4.061e-11,  4.149e-11,  4.299e-11,  4.223e-11,  4.251e-11,
        4.287e-11,  4.177e-11,  4.094e-11,  3.942e-11,  3.772e-11,
        3.614e-11,  3.394e-11,  3.222e-11,  2.791e-11,  2.665e-11,
        2.309e-11,  2.032e-11,  1.740e-11,  1.535e-11,  1.323e-11,
        1.151e-11,  9.803e-12,  8.650e-12,  7.540e-12,  6.619e-12,
        5.832e-12,  5.113e-12,  4.503e-12,  3.975e-12,  3.520e-12,
        3.112e-12,  2.797e-12,  2.500e-12,  2.240e-12,  2.013e-12,
        1.819e-12,  1.653e-12,  1.513e-12,  1.395e-12,  1.299e-12,
        1.225e-12,  1.168e-12,  1.124e-12,  1.148e-12,  1.107e-12,
        1.128e-12,  1.169e-12,  1.233e-12,  1.307e-12,  1.359e-12,
        1.543e-12,  1.686e-12,  1.794e-12,  2.028e-12,  2.210e-12,
        2.441e-12,  2.653e-12,  2.828e-12,  3.093e-12,  3.280e-12,
        3.551e-12,  3.677e-12,  3.803e-12,  3.844e-12,  4.068e-12,
        4.093e-12,  4.002e-12,  3.904e-12,  3.624e-12,  3.633e-12,
        3.622e-12,  3.443e-12,  3.184e-12,  2.934e-12,  2.476e-12,
        2.212e-12,  1.867e-12,  1.594e-12,  1.370e-12,  1.192e-12,
        1.045e-12,  9.211e-13,  8.170e-13,  7.290e-13,  6.550e-13,
        5.929e-13,  5.415e-13,  4.995e-13,  4.661e-13,  4.406e-13,
        4.225e-13,  4.116e-13,  4.075e-13,  4.102e-13,  4.198e-13,
        4.365e-13,  4.606e-13,  4.925e-13,  5.326e-13,  5.818e-13,
        6.407e-13,  7.104e-13,  7.920e-13,  8.868e-13,  9.964e-13,
        1.123e-12,  1.268e-12,  1.434e-12,  1.626e-12,  1.848e-12,
        2.107e-12,  2.422e-12,  2.772e-12,  3.145e-12,  3.704e-12,
        4.270e-12,  4.721e-12,  5.361e-12,  6.083e-12,  7.095e-12,
        7.968e-12,  9.228e-12,  1.048e-11,  1.187e-11,  1.336e-11,
        1.577e-11,  1.772e-11,  2.017e-11,  2.250e-11,  2.630e-11,
        2.911e-11,  3.356e-11,  3.820e-11,  4.173e-11,  4.811e-11,
        5.254e-11,  5.839e-11,  6.187e-11,  6.805e-11,  7.118e-11,
        7.369e-11,  7.664e-11,  7.794e-11,  7.947e-11,  8.036e-11,
        7.954e-11,  7.849e-11,  7.518e-11,  7.462e-11,  6.926e-11,
        6.531e-11,  6.197e-11,  5.421e-11,  4.777e-11,  4.111e-11,
        3.679e-11,  3.166e-11,  2.786e-11,  2.436e-11,  2.144e-11,
        1.859e-11,  1.628e-11,  1.414e-11,  1.237e-11,  1.093e-11,
        9.558e-12 };

// H2O self continuum parameters at T=260 K
// date of last update: 11/18/02
// units of (CM**3/MOL)*1.E-20
const Numeric SL260_ckd_mt_100_v1  =   -20.0;
const Numeric SL260_ckd_mt_100_v2  = 20000.0;
const Numeric SL260_ckd_mt_100_dv  =    10.0;
const int     SL260_ckd_mt_100_npt =  2003;
const double  SL260_ckd_mt_100[SL260_ckd_mt_100_npt+addF77fields] = {
        0.000e0,    2.749e-01,  2.732e-01,
        2.752e-01,  2.732e-01,  2.749e-01,  2.676e-01,  2.667e-01,
        2.545e-01,  2.497e-01,  2.327e-01,  2.218e-01,  2.036e-01,
        1.825e-01,  1.694e-01,  1.497e-01,  1.353e-01,  1.210e-01,
        1.014e-01,  9.405e-02,  7.848e-02,  7.195e-02,  6.246e-02,
        5.306e-02,  4.853e-02,  4.138e-02,  3.735e-02,  3.171e-02,
        2.785e-02,  2.431e-02,  2.111e-02,  1.845e-02,  1.640e-02,
        1.405e-02,  1.255e-02,  1.098e-02,  9.797e-03,  8.646e-03,
        7.779e-03,  6.898e-03,  6.099e-03,  5.453e-03,  4.909e-03,
        4.413e-03,  3.959e-03,  3.581e-03,  3.199e-03,  2.871e-03,
        2.583e-03,  2.330e-03,  2.086e-03,  1.874e-03,  1.684e-03,
        1.512e-03,  1.361e-03,  1.225e-03,  1.100e-03,  9.890e-04,
        8.916e-04,  8.039e-04,  7.256e-04,  6.545e-04,  5.918e-04,
        5.359e-04,  4.867e-04,  4.426e-04,  4.033e-04,  3.682e-04,
        3.366e-04,  3.085e-04,  2.833e-04,  2.605e-04,  2.403e-04,
        2.221e-04,  2.055e-04,  1.908e-04,  1.774e-04,  1.653e-04,
        1.544e-04,  1.443e-04,  1.351e-04,  1.267e-04,  1.190e-04,
        1.119e-04,  1.053e-04,  9.922e-05,  9.355e-05,  8.831e-05,
        8.339e-05,  7.878e-05,  7.449e-05,  7.043e-05,  6.664e-05,
        6.307e-05,  5.969e-05,  5.654e-05,  5.357e-05,  5.075e-05,
        4.810e-05,  4.560e-05,  4.322e-05,  4.102e-05,  3.892e-05,
        3.696e-05,  3.511e-05,  3.339e-05,  3.177e-05,  3.026e-05,
        2.886e-05,  2.756e-05,  2.636e-05,  2.527e-05,  2.427e-05,
        2.337e-05,  2.257e-05,  2.185e-05,  2.127e-05,  2.080e-05,
        2.041e-05,  2.013e-05,  2.000e-05,  1.997e-05,  2.009e-05,
        2.031e-05,  2.068e-05,  2.124e-05,  2.189e-05,  2.267e-05,
        2.364e-05,  2.463e-05,  2.618e-05,  2.774e-05,  2.937e-05,
        3.144e-05,  3.359e-05,  3.695e-05,  4.002e-05,  4.374e-05,
        4.947e-05,  5.431e-05,  6.281e-05,  7.169e-05,  8.157e-05,
        9.728e-05,  1.079e-04,  1.337e-04,  1.442e-04,  1.683e-04,
        1.879e-04,  2.223e-04,  2.425e-04,  2.838e-04,  3.143e-04,
        3.527e-04,  4.012e-04,  4.237e-04,  4.747e-04,  5.057e-04,
        5.409e-04,  5.734e-04,  5.944e-04,  6.077e-04,  6.175e-04,
        6.238e-04,  6.226e-04,  6.248e-04,  6.192e-04,  6.098e-04,
        5.818e-04,  5.709e-04,  5.465e-04,  5.043e-04,  4.699e-04,
        4.294e-04,  3.984e-04,  3.672e-04,  3.152e-04,  2.883e-04,
        2.503e-04,  2.211e-04,  1.920e-04,  1.714e-04,  1.485e-04,
        1.358e-04,  1.156e-04,  1.021e-04,  8.887e-05,  7.842e-05,
        7.120e-05,  6.186e-05,  5.730e-05,  4.792e-05,  4.364e-05,
        3.720e-05,  3.280e-05,  2.946e-05,  2.591e-05,  2.261e-05,
        2.048e-05,  1.813e-05,  1.630e-05,  1.447e-05,  1.282e-05,
        1.167e-05,  1.041e-05,  9.449e-06,  8.510e-06,  7.596e-06,
        6.961e-06,  6.272e-06,  5.728e-06,  5.198e-06,  4.667e-06,
        4.288e-06,  3.897e-06,  3.551e-06,  3.235e-06,  2.952e-06,
        2.688e-06,  2.449e-06,  2.241e-06,  2.050e-06,  1.879e-06,
        1.722e-06,  1.582e-06,  1.456e-06,  1.339e-06,  1.236e-06,
        1.144e-06,  1.060e-06,  9.830e-07,  9.149e-07,  8.535e-07,
        7.973e-07,  7.466e-07,  6.999e-07,  6.574e-07,  6.180e-07,
        5.821e-07,  5.487e-07,  5.180e-07,  4.896e-07,  4.631e-07,
        4.386e-07,  4.160e-07,  3.945e-07,  3.748e-07,  3.562e-07,
        3.385e-07,  3.222e-07,  3.068e-07,  2.922e-07,  2.788e-07,
        2.659e-07,  2.539e-07,  2.425e-07,  2.318e-07,  2.219e-07,
        2.127e-07,  2.039e-07,  1.958e-07,  1.885e-07,  1.818e-07,
        1.758e-07,  1.711e-07,  1.662e-07,  1.630e-07,  1.605e-07,
        1.580e-07,  1.559e-07,  1.545e-07,  1.532e-07,  1.522e-07,
        1.510e-07,  1.495e-07,  1.465e-07,  1.483e-07,  1.469e-07,
        1.448e-07,  1.444e-07,  1.436e-07,  1.426e-07,  1.431e-07,
        1.425e-07,  1.445e-07,  1.477e-07,  1.515e-07,  1.567e-07,
        1.634e-07,  1.712e-07,  1.802e-07,  1.914e-07,  2.024e-07,
        2.159e-07,  2.295e-07,  2.461e-07,  2.621e-07,  2.868e-07,
        3.102e-07,  3.394e-07,  3.784e-07,  4.223e-07,  4.864e-07,
        5.501e-07,  6.039e-07,  7.193e-07,  7.728e-07,  9.514e-07,
        1.073e-06,  1.180e-06,  1.333e-06,  1.472e-06,  1.566e-06,
        1.677e-06,  1.784e-06,  1.904e-06,  1.953e-06,  2.020e-06,
        2.074e-06,  2.128e-06,  2.162e-06,  2.219e-06,  2.221e-06,
        2.249e-06,  2.239e-06,  2.235e-06,  2.185e-06,  2.141e-06,
        2.124e-06,  2.090e-06,  2.068e-06,  2.100e-06,  2.104e-06,
        2.142e-06,  2.181e-06,  2.257e-06,  2.362e-06,  2.500e-06,
        2.664e-06,  2.884e-06,  3.189e-06,  3.480e-06,  3.847e-06,
        4.313e-06,  4.790e-06,  5.250e-06,  5.989e-06,  6.692e-06,
        7.668e-06,  8.520e-06,  9.606e-06,  1.073e-05,  1.225e-05,
        1.377e-05,  1.582e-05,  1.761e-05,  2.029e-05,  2.284e-05,
        2.602e-05,  2.940e-05,  3.483e-05,  3.928e-05,  4.618e-05,
        5.240e-05,  6.132e-05,  7.183e-05,  8.521e-05,  9.111e-05,
        1.070e-04,  1.184e-04,  1.264e-04,  1.475e-04,  1.612e-04,
        1.704e-04,  1.818e-04,  1.924e-04,  1.994e-04,  2.061e-04,
        2.180e-04,  2.187e-04,  2.200e-04,  2.196e-04,  2.131e-04,
        2.015e-04,  1.988e-04,  1.847e-04,  1.729e-04,  1.597e-04,
        1.373e-04,  1.262e-04,  1.087e-04,  9.439e-05,  8.061e-05,
        7.093e-05,  6.049e-05,  5.120e-05,  4.435e-05,  3.817e-05,
        3.340e-05,  2.927e-05,  2.573e-05,  2.291e-05,  2.040e-05,
        1.827e-05,  1.636e-05,  1.463e-05,  1.309e-05,  1.170e-05,
        1.047e-05,  9.315e-06,  8.328e-06,  7.458e-06,  6.665e-06,
        5.940e-06,  5.316e-06,  4.752e-06,  4.252e-06,  3.825e-06,
        3.421e-06,  3.064e-06,  2.746e-06,  2.465e-06,  2.216e-06,
        1.990e-06,  1.790e-06,  1.609e-06,  1.449e-06,  1.306e-06,
        1.177e-06,  1.063e-06,  9.607e-07,  8.672e-07,  7.855e-07,
        7.118e-07,  6.460e-07,  5.871e-07,  5.340e-07,  4.868e-07,
        4.447e-07,  4.068e-07,  3.729e-07,  3.423e-07,  3.151e-07,
        2.905e-07,  2.686e-07,  2.484e-07,  2.306e-07,  2.142e-07,
        1.995e-07,  1.860e-07,  1.738e-07,  1.626e-07,  1.522e-07,
        1.427e-07,  1.338e-07,  1.258e-07,  1.183e-07,  1.116e-07,
        1.056e-07,  9.972e-08,  9.460e-08,  9.007e-08,  8.592e-08,
        8.195e-08,  7.816e-08,  7.483e-08,  7.193e-08,  6.892e-08,
        6.642e-08,  6.386e-08,  6.154e-08,  5.949e-08,  5.764e-08,
        5.622e-08,  5.479e-08,  5.364e-08,  5.301e-08,  5.267e-08,
        5.263e-08,  5.313e-08,  5.410e-08,  5.550e-08,  5.745e-08,
        6.003e-08,  6.311e-08,  6.713e-08,  7.173e-08,  7.724e-08,
        8.368e-08,  9.121e-08,  9.986e-08,  1.097e-07,  1.209e-07,
        1.338e-07,  1.486e-07,  1.651e-07,  1.837e-07,  2.048e-07,
        2.289e-07,  2.557e-07,  2.857e-07,  3.195e-07,  3.587e-07,
        4.015e-07,  4.497e-07,  5.049e-07,  5.665e-07,  6.366e-07,
        7.121e-07,  7.996e-07,  8.946e-07,  1.002e-06,  1.117e-06,
        1.262e-06,  1.416e-06,  1.611e-06,  1.807e-06,  2.056e-06,
        2.351e-06,  2.769e-06,  3.138e-06,  3.699e-06,  4.386e-06,
        5.041e-06,  6.074e-06,  6.812e-06,  7.790e-06,  8.855e-06,
        1.014e-05,  1.095e-05,  1.245e-05,  1.316e-05,  1.390e-05,
        1.504e-05,  1.583e-05,  1.617e-05,  1.652e-05,  1.713e-05,
        1.724e-05,  1.715e-05,  1.668e-05,  1.629e-05,  1.552e-05,
        1.478e-05,  1.340e-05,  1.245e-05,  1.121e-05,  9.575e-06,
        8.956e-06,  7.345e-06,  6.597e-06,  5.612e-06,  4.818e-06,
        4.165e-06,  3.579e-06,  3.041e-06,  2.623e-06,  2.290e-06,
        1.984e-06,  1.748e-06,  1.534e-06,  1.369e-06,  1.219e-06,
        1.092e-06,  9.800e-07,  8.762e-07,  7.896e-07,  7.104e-07,
        6.364e-07,  5.691e-07,  5.107e-07,  4.575e-07,  4.090e-07,
        3.667e-07,  3.287e-07,  2.931e-07,  2.633e-07,  2.356e-07,
        2.111e-07,  1.895e-07,  1.697e-07,  1.525e-07,  1.369e-07,
        1.233e-07,  1.114e-07,  9.988e-08,  9.004e-08,  8.149e-08,
        7.352e-08,  6.662e-08,  6.030e-08,  5.479e-08,  4.974e-08,
        4.532e-08,  4.129e-08,  3.781e-08,  3.462e-08,  3.176e-08,
        2.919e-08,  2.687e-08,  2.481e-08,  2.292e-08,  2.119e-08,
        1.967e-08,  1.828e-08,  1.706e-08,  1.589e-08,  1.487e-08,
        1.393e-08,  1.307e-08,  1.228e-08,  1.156e-08,  1.089e-08,
        1.028e-08,  9.696e-09,  9.159e-09,  8.658e-09,  8.187e-09,
        7.746e-09,  7.340e-09,  6.953e-09,  6.594e-09,  6.259e-09,
        5.948e-09,  5.660e-09,  5.386e-09,  5.135e-09,  4.903e-09,
        4.703e-09,  4.515e-09,  4.362e-09,  4.233e-09,  4.117e-09,
        4.017e-09,  3.962e-09,  3.924e-09,  3.905e-09,  3.922e-09,
        3.967e-09,  4.046e-09,  4.165e-09,  4.320e-09,  4.522e-09,
        4.769e-09,  5.083e-09,  5.443e-09,  5.872e-09,  6.366e-09,
        6.949e-09,  7.601e-09,  8.371e-09,  9.220e-09,  1.020e-08,
        1.129e-08,  1.251e-08,  1.393e-08,  1.542e-08,  1.720e-08,
        1.926e-08,  2.152e-08,  2.392e-08,  2.678e-08,  3.028e-08,
        3.390e-08,  3.836e-08,  4.309e-08,  4.900e-08,  5.481e-08,
        6.252e-08,  7.039e-08,  7.883e-08,  8.849e-08,  1.012e-07,
        1.142e-07,  1.300e-07,  1.475e-07,  1.732e-07,  1.978e-07,
        2.304e-07,  2.631e-07,  2.988e-07,  3.392e-07,  3.690e-07,
        4.355e-07,  4.672e-07,  5.110e-07,  5.461e-07,  5.828e-07,
        6.233e-07,  6.509e-07,  6.672e-07,  6.969e-07,  7.104e-07,
        7.439e-07,  7.463e-07,  7.708e-07,  7.466e-07,  7.668e-07,
        7.549e-07,  7.586e-07,  7.384e-07,  7.439e-07,  7.785e-07,
        7.915e-07,  8.310e-07,  8.745e-07,  9.558e-07,  1.038e-06,
        1.173e-06,  1.304e-06,  1.452e-06,  1.671e-06,  1.931e-06,
        2.239e-06,  2.578e-06,  3.032e-06,  3.334e-06,  3.980e-06,
        4.300e-06,  4.518e-06,  5.321e-06,  5.508e-06,  6.211e-06,
        6.590e-06,  7.046e-06,  7.555e-06,  7.558e-06,  7.875e-06,
        8.319e-06,  8.433e-06,  8.590e-06,  8.503e-06,  8.304e-06,
        8.336e-06,  7.739e-06,  7.301e-06,  6.827e-06,  6.078e-06,
        5.551e-06,  4.762e-06,  4.224e-06,  3.538e-06,  2.984e-06,
        2.619e-06,  2.227e-06,  1.923e-06,  1.669e-06,  1.462e-06,
        1.294e-06,  1.155e-06,  1.033e-06,  9.231e-07,  8.238e-07,
        7.360e-07,  6.564e-07,  5.869e-07,  5.236e-07,  4.673e-07,
        4.174e-07,  3.736e-07,  3.330e-07,  2.976e-07,  2.657e-07,
        2.367e-07,  2.106e-07,  1.877e-07,  1.671e-07,  1.494e-07,
        1.332e-07,  1.192e-07,  1.065e-07,  9.558e-08,  8.586e-08,
        7.717e-08,  6.958e-08,  6.278e-08,  5.666e-08,  5.121e-08,
        4.647e-08,  4.213e-08,  3.815e-08,  3.459e-08,  3.146e-08,
        2.862e-08,  2.604e-08,  2.375e-08,  2.162e-08,  1.981e-08,
        1.817e-08,  1.670e-08,  1.537e-08,  1.417e-08,  1.310e-08,
        1.215e-08,  1.128e-08,  1.050e-08,  9.793e-09,  9.158e-09,
        8.586e-09,  8.068e-09,  7.595e-09,  7.166e-09,  6.778e-09,
        6.427e-09,  6.108e-09,  5.826e-09,  5.571e-09,  5.347e-09,
        5.144e-09,  4.968e-09,  4.822e-09,  4.692e-09,  4.589e-09,
        4.506e-09,  4.467e-09,  4.440e-09,  4.466e-09,  4.515e-09,
        4.718e-09,  4.729e-09,  4.937e-09,  5.249e-09,  5.466e-09,
        5.713e-09,  6.030e-09,  6.436e-09,  6.741e-09,  7.330e-09,
        7.787e-09,  8.414e-09,  8.908e-09,  9.868e-09,  1.069e-08,
        1.158e-08,  1.253e-08,  1.300e-08,  1.409e-08,  1.470e-08,
        1.548e-08,  1.612e-08,  1.666e-08,  1.736e-08,  1.763e-08,
        1.812e-08,  1.852e-08,  1.923e-08,  1.897e-08,  1.893e-08,
        1.888e-08,  1.868e-08,  1.895e-08,  1.899e-08,  1.876e-08,
        1.960e-08,  2.020e-08,  2.121e-08,  2.239e-08,  2.379e-08,
        2.526e-08,  2.766e-08,  2.994e-08,  3.332e-08,  3.703e-08,
        4.158e-08,  4.774e-08,  5.499e-08,  6.355e-08,  7.349e-08,
        8.414e-08,  9.846e-08,  1.143e-07,  1.307e-07,  1.562e-07,
        1.817e-07,  2.011e-07,  2.192e-07,  2.485e-07,  2.867e-07,
        3.035e-07,  3.223e-07,  3.443e-07,  3.617e-07,  3.793e-07,
        3.793e-07,  3.839e-07,  4.081e-07,  4.117e-07,  4.085e-07,
        3.920e-07,  3.851e-07,  3.754e-07,  3.490e-07,  3.229e-07,
        2.978e-07,  2.691e-07,  2.312e-07,  2.029e-07,  1.721e-07,
        1.472e-07,  1.308e-07,  1.132e-07,  9.736e-08,  8.458e-08,
        7.402e-08,  6.534e-08,  5.811e-08,  5.235e-08,  4.762e-08,
        4.293e-08,  3.896e-08,  3.526e-08,  3.165e-08,  2.833e-08,
        2.551e-08,  2.288e-08,  2.036e-08,  1.820e-08,  1.626e-08,
        1.438e-08,  1.299e-08,  1.149e-08,  1.030e-08,  9.148e-09,
        8.122e-09,  7.264e-09,  6.425e-09,  5.777e-09,  5.060e-09,
        4.502e-09,  4.013e-09,  3.567e-09,  3.145e-09,  2.864e-09,
        2.553e-09,  2.311e-09,  2.087e-09,  1.886e-09,  1.716e-09,
        1.556e-09,  1.432e-09,  1.311e-09,  1.202e-09,  1.104e-09,
        1.013e-09,  9.293e-10,  8.493e-10,  7.790e-10,  7.185e-10,
        6.642e-10,  6.141e-10,  5.684e-10,  5.346e-10,  5.032e-10,
        4.725e-10,  4.439e-10,  4.176e-10,  3.930e-10,  3.714e-10,
        3.515e-10,  3.332e-10,  3.167e-10,  3.020e-10,  2.887e-10,
        2.769e-10,  2.665e-10,  2.578e-10,  2.503e-10,  2.436e-10,
        2.377e-10,  2.342e-10,  2.305e-10,  2.296e-10,  2.278e-10,
        2.321e-10,  2.355e-10,  2.402e-10,  2.478e-10,  2.670e-10,
        2.848e-10,  2.982e-10,  3.263e-10,  3.438e-10,  3.649e-10,
        3.829e-10,  4.115e-10,  4.264e-10,  4.473e-10,  4.630e-10,
        4.808e-10,  4.995e-10,  5.142e-10,  5.313e-10,  5.318e-10,
        5.358e-10,  5.452e-10,  5.507e-10,  5.698e-10,  5.782e-10,
        5.983e-10,  6.164e-10,  6.532e-10,  6.811e-10,  7.624e-10,
        8.302e-10,  9.067e-10,  9.937e-10,  1.104e-09,  1.221e-09,
        1.361e-09,  1.516e-09,  1.675e-09,  1.883e-09,  2.101e-09,
        2.349e-09,  2.614e-09,  2.920e-09,  3.305e-09,  3.724e-09,
        4.142e-09,  4.887e-09,  5.614e-09,  6.506e-09,  7.463e-09,
        8.817e-09,  9.849e-09,  1.187e-08,  1.321e-08,  1.474e-08,
        1.698e-08,  1.794e-08,  2.090e-08,  2.211e-08,  2.362e-08,
        2.556e-08,  2.729e-08,  2.880e-08,  3.046e-08,  3.167e-08,
        3.367e-08,  3.457e-08,  3.590e-08,  3.711e-08,  3.826e-08,
        4.001e-08,  4.211e-08,  4.315e-08,  4.661e-08,  5.010e-08,
        5.249e-08,  5.840e-08,  6.628e-08,  7.512e-08,  8.253e-08,
        9.722e-08,  1.067e-07,  1.153e-07,  1.347e-07,  1.428e-07,
        1.577e-07,  1.694e-07,  1.833e-07,  1.938e-07,  2.108e-07,
        2.059e-07,  2.157e-07,  2.185e-07,  2.208e-07,  2.182e-07,
        2.093e-07,  2.014e-07,  1.962e-07,  1.819e-07,  1.713e-07,
        1.510e-07,  1.340e-07,  1.154e-07,  9.890e-08,  8.880e-08,
        7.673e-08,  6.599e-08,  5.730e-08,  5.081e-08,  4.567e-08,
        4.147e-08,  3.773e-08,  3.460e-08,  3.194e-08,  2.953e-08,
        2.759e-08,  2.594e-08,  2.442e-08,  2.355e-08,  2.283e-08,
        2.279e-08,  2.231e-08,  2.279e-08,  2.239e-08,  2.210e-08,
        2.309e-08,  2.293e-08,  2.352e-08,  2.415e-08,  2.430e-08,
        2.426e-08,  2.465e-08,  2.500e-08,  2.496e-08,  2.465e-08,
        2.445e-08,  2.383e-08,  2.299e-08,  2.165e-08,  2.113e-08,
        1.968e-08,  1.819e-08,  1.644e-08,  1.427e-08,  1.270e-08,
        1.082e-08,  9.428e-09,  8.091e-09,  6.958e-09,  5.988e-09,
        5.246e-09,  4.601e-09,  4.098e-09,  3.664e-09,  3.287e-09,
        2.942e-09,  2.656e-09,  2.364e-09,  2.118e-09,  1.903e-09,
        1.703e-09,  1.525e-09,  1.365e-09,  1.229e-09,  1.107e-09,
        9.960e-10,  8.945e-10,  8.080e-10,  7.308e-10,  6.616e-10,
        5.994e-10,  5.422e-10,  4.929e-10,  4.478e-10,  4.070e-10,
        3.707e-10,  3.379e-10,  3.087e-10,  2.823e-10,  2.592e-10,
        2.385e-10,  2.201e-10,  2.038e-10,  1.897e-10,  1.774e-10,
        1.667e-10,  1.577e-10,  1.502e-10,  1.437e-10,  1.394e-10,
        1.358e-10,  1.324e-10,  1.329e-10,  1.324e-10,  1.360e-10,
        1.390e-10,  1.424e-10,  1.544e-10,  1.651e-10,  1.817e-10,
        1.984e-10,  2.195e-10,  2.438e-10,  2.700e-10,  2.991e-10,
        3.322e-10,  3.632e-10,  3.957e-10,  4.360e-10,  4.701e-10,
        5.030e-10,  5.381e-10,  5.793e-10,  6.190e-10,  6.596e-10,
        7.004e-10,  7.561e-10,  7.934e-10,  8.552e-10,  9.142e-10,
        9.570e-10,  1.027e-09,  1.097e-09,  1.193e-09,  1.334e-09,
        1.470e-09,  1.636e-09,  1.871e-09,  2.122e-09,  2.519e-09,
        2.806e-09,  3.203e-09,  3.846e-09,  4.362e-09,  5.114e-09,
        5.643e-09,  6.305e-09,  6.981e-09,  7.983e-09,  8.783e-09,
        9.419e-09,  1.017e-08,  1.063e-08,  1.121e-08,  1.130e-08,
        1.201e-08,  1.225e-08,  1.232e-08,  1.223e-08,  1.177e-08,
        1.151e-08,  1.116e-08,  1.047e-08,  9.698e-09,  8.734e-09,
        8.202e-09,  7.041e-09,  6.074e-09,  5.172e-09,  4.468e-09,
        3.913e-09,  3.414e-09,  2.975e-09,  2.650e-09,  2.406e-09,
        2.173e-09,  2.009e-09,  1.861e-09,  1.727e-09,  1.612e-09,
        1.514e-09,  1.430e-09,  1.362e-09,  1.333e-09,  1.288e-09,
        1.249e-09,  1.238e-09,  1.228e-09,  1.217e-09,  1.202e-09,
        1.209e-09,  1.177e-09,  1.157e-09,  1.165e-09,  1.142e-09,
        1.131e-09,  1.138e-09,  1.117e-09,  1.100e-09,  1.069e-09,
        1.023e-09,  1.005e-09,  9.159e-10,  8.863e-10,  7.865e-10,
        7.153e-10,  6.247e-10,  5.846e-10,  5.133e-10,  4.360e-10,
        3.789e-10,  3.335e-10,  2.833e-10,  2.483e-10,  2.155e-10,
        1.918e-10,  1.709e-10,  1.529e-10,  1.374e-10,  1.235e-10,
        1.108e-10,  9.933e-11,  8.932e-11,  8.022e-11,  7.224e-11,
        6.520e-11,  5.896e-11,  5.328e-11,  4.813e-11,  4.365e-11,
        3.961e-11,  3.594e-11,  3.266e-11,  2.967e-11,  2.701e-11,
        2.464e-11,  2.248e-11,  2.054e-11,  1.878e-11,  1.721e-11,
        1.579e-11,  1.453e-11,  1.341e-11,  1.241e-11,  1.154e-11,
        1.078e-11,  1.014e-11,  9.601e-12,  9.167e-12,  8.838e-12,
        8.614e-12,  8.493e-12,  8.481e-12,  8.581e-12,  8.795e-12,
        9.131e-12,  9.601e-12,  1.021e-11,  1.097e-11,  1.191e-11,
        1.303e-11,  1.439e-11,  1.601e-11,  1.778e-11,  1.984e-11,
        2.234e-11,  2.474e-11,  2.766e-11,  3.085e-11,  3.415e-11,
        3.821e-11,  4.261e-11,  4.748e-11,  5.323e-11,  5.935e-11,
        6.619e-11,  7.418e-11,  8.294e-11,  9.260e-11,  1.039e-10,
        1.156e-10,  1.297e-10,  1.460e-10,  1.641e-10,  1.858e-10,
        2.100e-10,  2.383e-10,  2.724e-10,  3.116e-10,  3.538e-10,
        4.173e-10,  4.727e-10,  5.503e-10,  6.337e-10,  7.320e-10,
        8.298e-10,  9.328e-10,  1.059e-09,  1.176e-09,  1.328e-09,
        1.445e-09,  1.593e-09,  1.770e-09,  1.954e-09,  2.175e-09,
        2.405e-09,  2.622e-09,  2.906e-09,  3.294e-09,  3.713e-09,
        3.980e-09,  4.384e-09,  4.987e-09,  5.311e-09,  5.874e-09,
        6.337e-09,  7.027e-09,  7.390e-09,  7.769e-09,  8.374e-09,
        8.605e-09,  9.165e-09,  9.415e-09,  9.511e-09,  9.704e-09,
        9.588e-09,  9.450e-09,  9.086e-09,  8.798e-09,  8.469e-09,
        7.697e-09,  7.168e-09,  6.255e-09,  5.772e-09,  4.970e-09,
        4.271e-09,  3.653e-09,  3.154e-09,  2.742e-09,  2.435e-09,
        2.166e-09,  1.936e-09,  1.731e-09,  1.556e-09,  1.399e-09,
        1.272e-09,  1.157e-09,  1.066e-09,  9.844e-10,  9.258e-10,
        8.787e-10,  8.421e-10,  8.083e-10,  8.046e-10,  8.067e-10,
        8.181e-10,  8.325e-10,  8.517e-10,  9.151e-10,  9.351e-10,
        9.677e-10,  1.071e-09,  1.126e-09,  1.219e-09,  1.297e-09,
        1.408e-09,  1.476e-09,  1.517e-09,  1.600e-09,  1.649e-09,
        1.678e-09,  1.746e-09,  1.742e-09,  1.728e-09,  1.699e-09,
        1.655e-09,  1.561e-09,  1.480e-09,  1.451e-09,  1.411e-09,
        1.171e-09,  1.106e-09,  9.714e-10,  8.523e-10,  7.346e-10,
        6.241e-10,  5.371e-10,  4.704e-10,  4.144e-10,  3.683e-10,
        3.292e-10,  2.942e-10,  2.620e-10,  2.341e-10,  2.104e-10,
        1.884e-10,  1.700e-10,  1.546e-10,  1.394e-10,  1.265e-10,
        1.140e-10,  1.019e-10,  9.279e-11,  8.283e-11,  7.458e-11,
        6.668e-11,  5.976e-11,  5.330e-11,  4.794e-11,  4.289e-11,
        3.841e-11,  3.467e-11,  3.130e-11,  2.832e-11,  2.582e-11,
        2.356e-11,  2.152e-11,  1.970e-11,  1.808e-11,  1.664e-11,
        1.539e-11,  1.434e-11,  1.344e-11,  1.269e-11,  1.209e-11,
        1.162e-11,  1.129e-11,  1.108e-11,  1.099e-11,  1.103e-11,
        1.119e-11,  1.148e-11,  1.193e-11,  1.252e-11,  1.329e-11,
        1.421e-11,  1.555e-11,  1.685e-11,  1.839e-11,  2.054e-11,
        2.317e-11,  2.571e-11,  2.839e-11,  3.171e-11,  3.490e-11,
        3.886e-11,  4.287e-11,  4.645e-11,  5.047e-11,  5.592e-11,
        6.109e-11,  6.628e-11,  7.381e-11,  8.088e-11,  8.966e-11,
        1.045e-10,  1.120e-10,  1.287e-10,  1.486e-10,  1.662e-10,
        1.866e-10,  2.133e-10,  2.524e-10,  2.776e-10,  3.204e-10,
        3.559e-10,  4.028e-10,  4.448e-10,  4.882e-10,  5.244e-10,
        5.605e-10,  6.018e-10,  6.328e-10,  6.579e-10,  6.541e-10,
        7.024e-10,  7.074e-10,  7.068e-10,  7.009e-10,  6.698e-10,
        6.545e-10,  6.209e-10,  5.834e-10,  5.412e-10,  5.001e-10,
        4.231e-10,  3.727e-10,  3.211e-10,  2.833e-10,  2.447e-10,
        2.097e-10,  1.843e-10,  1.639e-10,  1.449e-10,  1.270e-10,
        1.161e-10,  1.033e-10,  9.282e-11,  8.407e-11,  7.639e-11,
        7.023e-11,  6.474e-11,  6.142e-11,  5.760e-11,  5.568e-11,
        5.472e-11,  5.390e-11,  5.455e-11,  5.540e-11,  5.587e-11,
        6.230e-11,  6.490e-11,  6.868e-11,  7.382e-11,  8.022e-11,
        8.372e-11,  9.243e-11,  1.004e-10,  1.062e-10,  1.130e-10,
        1.176e-10,  1.244e-10,  1.279e-10,  1.298e-10,  1.302e-10,
        1.312e-10,  1.295e-10,  1.244e-10,  1.211e-10,  1.167e-10,
        1.098e-10,  9.927e-11,  8.854e-11,  8.011e-11,  7.182e-11,
        5.923e-11,  5.212e-11,  4.453e-11,  3.832e-11,  3.371e-11,
        2.987e-11,  2.651e-11,  2.354e-11,  2.093e-11,  1.863e-11,
        1.662e-11,  1.486e-11,  1.331e-11,  1.193e-11,  1.071e-11,
        9.628e-12,  8.660e-12,  7.801e-12,  7.031e-12,  6.347e-12,
        5.733e-12,  5.182e-12,  4.695e-12,  4.260e-12,  3.874e-12,
        3.533e-12,  3.235e-12,  2.979e-12,  2.760e-12,  2.579e-12,
        2.432e-12,  2.321e-12,  2.246e-12,  2.205e-12,  2.196e-12,
        2.223e-12,  2.288e-12,  2.387e-12,  2.525e-12,  2.704e-12,
        2.925e-12,  3.191e-12,  3.508e-12,  3.876e-12,  4.303e-12,
        4.793e-12,  5.347e-12,  5.978e-12,  6.682e-12,  7.467e-12,
        8.340e-12,  9.293e-12,  1.035e-11,  1.152e-11,  1.285e-11,
        1.428e-11,  1.586e-11,  1.764e-11,  1.972e-11,  2.214e-11,
        2.478e-11,  2.776e-11,  3.151e-11,  3.591e-11,  4.103e-11,
        4.660e-11,  5.395e-11,  6.306e-11,  7.172e-11,  8.358e-11,
        9.670e-11,  1.110e-10,  1.325e-10,  1.494e-10,  1.736e-10,
        2.007e-10,  2.296e-10,  2.608e-10,  3.004e-10,  3.361e-10,
        3.727e-10,  4.373e-10,  4.838e-10,  5.483e-10,  6.006e-10,
        6.535e-10,  6.899e-10,  7.687e-10,  8.444e-10,  8.798e-10,
        9.135e-10,  9.532e-10,  9.757e-10,  9.968e-10,  1.006e-09,
        9.949e-10,  9.789e-10,  9.564e-10,  9.215e-10,  8.510e-10,
        8.394e-10,  7.707e-10,  7.152e-10,  6.274e-10,  5.598e-10,
        5.028e-10,  4.300e-10,  3.710e-10,  3.245e-10,  2.809e-10,
        2.461e-10,  2.154e-10,  1.910e-10,  1.685e-10,  1.487e-10,
        1.313e-10,  1.163e-10,  1.031e-10,  9.172e-11,  8.221e-11,
        7.382e-11,  6.693e-11,  6.079e-11,  5.581e-11,  5.167e-11,
        4.811e-11,  4.506e-11,  4.255e-11,  4.083e-11,  3.949e-11,
        3.881e-11,  3.861e-11,  3.858e-11,  3.951e-11,  4.045e-11,
        4.240e-11,  4.487e-11,  4.806e-11,  5.133e-11,  5.518e-11,
        5.919e-11,  6.533e-11,  7.031e-11,  7.762e-11,  8.305e-11,
        9.252e-11,  9.727e-11,  1.045e-10,  1.117e-10,  1.200e-10,
        1.275e-10,  1.341e-10,  1.362e-10,  1.438e-10,  1.450e-10,
        1.455e-10,  1.455e-10,  1.434e-10,  1.381e-10,  1.301e-10,
        1.276e-10,  1.163e-10,  1.089e-10,  9.911e-11,  8.943e-11,
        7.618e-11,  6.424e-11,  5.717e-11,  4.866e-11,  4.257e-11,
        3.773e-11,  3.331e-11,  2.958e-11,  2.629e-11,  2.316e-11,
        2.073e-11,  1.841e-11,  1.635e-11,  1.464e-11,  1.310e-11,
        1.160e-11,  1.047e-11,  9.408e-12,  8.414e-12,  7.521e-12,
        6.705e-12,  5.993e-12,  5.371e-12,  4.815e-12,  4.338e-12,
        3.921e-12,  3.567e-12,  3.265e-12,  3.010e-12,  2.795e-12,
        2.613e-12,  2.464e-12,  2.346e-12,  2.256e-12,  2.195e-12,
        2.165e-12,  2.166e-12,  2.198e-12,  2.262e-12,  2.364e-12,
        2.502e-12,  2.682e-12,  2.908e-12,  3.187e-12,  3.533e-12,
        3.946e-12,  4.418e-12,  5.013e-12,  5.708e-12,  6.379e-12,
        7.430e-12,  8.390e-12,  9.510e-12,  1.078e-11,  1.259e-11,
        1.438e-11,  1.630e-11,  1.814e-11,  2.055e-11,  2.348e-11,
        2.664e-11,  2.956e-11,  3.300e-11,  3.677e-11,  4.032e-11,
        4.494e-11,  4.951e-11,  5.452e-11,  6.014e-11,  6.500e-11,
        6.915e-11,  7.450e-11,  7.971e-11,  8.468e-11,  8.726e-11,
        8.995e-11,  9.182e-11,  9.509e-11,  9.333e-11,  9.386e-11,
        9.457e-11,  9.210e-11,  9.019e-11,  8.680e-11,  8.298e-11,
        7.947e-11,  7.460e-11,  7.082e-11,  6.132e-11,  5.855e-11,
        5.073e-11,  4.464e-11,  3.825e-11,  3.375e-11,  2.911e-11,
        2.535e-11,  2.160e-11,  1.907e-11,  1.665e-11,  1.463e-11,
        1.291e-11,  1.133e-11,  9.997e-12,  8.836e-12,  7.839e-12,
        6.943e-12,  6.254e-12,  5.600e-12,  5.029e-12,  4.529e-12,
        4.102e-12,  3.737e-12,  3.428e-12,  3.169e-12,  2.959e-12,
        2.798e-12,  2.675e-12,  2.582e-12,  2.644e-12,  2.557e-12,
        2.614e-12,  2.717e-12,  2.874e-12,  3.056e-12,  3.187e-12,
        3.631e-12,  3.979e-12,  4.248e-12,  4.817e-12,  5.266e-12,
        5.836e-12,  6.365e-12,  6.807e-12,  7.470e-12,  7.951e-12,
        8.636e-12,  8.972e-12,  9.314e-12,  9.445e-12,  1.003e-11,
        1.013e-11,  9.937e-12,  9.729e-12,  9.064e-12,  9.119e-12,
        9.124e-12,  8.704e-12,  8.078e-12,  7.470e-12,  6.329e-12,
        5.674e-12,  4.808e-12,  4.119e-12,  3.554e-12,  3.103e-12,
        2.731e-12,  2.415e-12,  2.150e-12,  1.926e-12,  1.737e-12,
        1.578e-12,  1.447e-12,  1.340e-12,  1.255e-12,  1.191e-12,
        1.146e-12,  1.121e-12,  1.114e-12,  1.126e-12,  1.156e-12,
        1.207e-12,  1.278e-12,  1.372e-12,  1.490e-12,  1.633e-12,
        1.805e-12,  2.010e-12,  2.249e-12,  2.528e-12,  2.852e-12,
        3.228e-12,  3.658e-12,  4.153e-12,  4.728e-12,  5.394e-12,
        6.176e-12,  7.126e-12,  8.188e-12,  9.328e-12,  1.103e-11,
        1.276e-11,  1.417e-11,  1.615e-11,  1.840e-11,  2.155e-11,
        2.429e-11,  2.826e-11,  3.222e-11,  3.664e-11,  4.140e-11,
        4.906e-11,  5.536e-11,  6.327e-11,  7.088e-11,  8.316e-11,
        9.242e-11,  1.070e-10,  1.223e-10,  1.341e-10,  1.553e-10,
        1.703e-10,  1.900e-10,  2.022e-10,  2.233e-10,  2.345e-10,
        2.438e-10,  2.546e-10,  2.599e-10,  2.661e-10,  2.703e-10,
        2.686e-10,  2.662e-10,  2.560e-10,  2.552e-10,  2.378e-10,
        2.252e-10,  2.146e-10,  1.885e-10,  1.668e-10,  1.441e-10,
        1.295e-10,  1.119e-10,  9.893e-11,  8.687e-11,  7.678e-11,
        6.685e-11,  5.879e-11,  5.127e-11,  4.505e-11,  3.997e-11,
        3.511e-11};


// H2O foreign continuum parameters at all temperatures
// date of last update: 11/18/02
// units of (CM**3/MOL)*1.E-20
const Numeric FH2O_ckd_mt_100_v1  =   -20.0;
const Numeric FH2O_ckd_mt_100_v2  = 20000.0;
const Numeric FH2O_ckd_mt_100_dv  =    10.0;
const int     FH2O_ckd_mt_100_npt =  2003;
const double  FH2O_ckd_mt_100[FH2O_ckd_mt_100_npt+addF77fields] = {
        0.000e0,    1.205e-02,  1.126e-02,
        1.095e-02,  1.126e-02,  1.205e-02,  1.322e-02,  1.430e-02,
        1.506e-02,  1.548e-02,  1.534e-02,  1.486e-02,  1.373e-02,
        1.262e-02,  1.134e-02,  1.001e-02,  8.702e-03,  7.475e-03,
        6.481e-03,  5.480e-03,  4.600e-03,  3.833e-03,  3.110e-03,
        2.543e-03,  2.049e-03,  1.680e-03,  1.374e-03,  1.046e-03,
        8.193e-04,  6.267e-04,  4.968e-04,  3.924e-04,  2.983e-04,
        2.477e-04,  1.997e-04,  1.596e-04,  1.331e-04,  1.061e-04,
        8.942e-05,  7.168e-05,  5.887e-05,  4.848e-05,  3.817e-05,
        3.170e-05,  2.579e-05,  2.162e-05,  1.768e-05,  1.490e-05,
        1.231e-05,  1.013e-05,  8.555e-06,  7.328e-06,  6.148e-06,
        5.207e-06,  4.387e-06,  3.741e-06,  3.220e-06,  2.753e-06,
        2.346e-06,  1.985e-06,  1.716e-06,  1.475e-06,  1.286e-06,
        1.122e-06,  9.661e-07,  8.284e-07,  7.057e-07,  6.119e-07,
        5.290e-07,  4.571e-07,  3.948e-07,  3.432e-07,  2.983e-07,
        2.589e-07,  2.265e-07,  1.976e-07,  1.704e-07,  1.456e-07,
        1.260e-07,  1.101e-07,  9.648e-08,  8.415e-08,  7.340e-08,
        6.441e-08,  5.643e-08,  4.940e-08,  4.276e-08,  3.703e-08,
        3.227e-08,  2.825e-08,  2.478e-08,  2.174e-08,  1.898e-08,
        1.664e-08,  1.458e-08,  1.278e-08,  1.126e-08,  9.891e-09,
        8.709e-09,  7.652e-09,  6.759e-09,  5.975e-09,  5.310e-09,
        4.728e-09,  4.214e-09,  3.792e-09,  3.463e-09,  3.226e-09,
        2.992e-09,  2.813e-09,  2.749e-09,  2.809e-09,  2.913e-09,
        3.037e-09,  3.413e-09,  3.738e-09,  4.189e-09,  4.808e-09,
        5.978e-09,  7.088e-09,  8.071e-09,  9.610e-09,  1.210e-08,
        1.500e-08,  1.764e-08,  2.221e-08,  2.898e-08,  3.948e-08,
        5.068e-08,  6.227e-08,  7.898e-08,  1.033e-07,  1.437e-07,
        1.889e-07,  2.589e-07,  3.590e-07,  4.971e-07,  7.156e-07,
        9.983e-07,  1.381e-06,  1.929e-06,  2.591e-06,  3.453e-06,
        4.570e-06,  5.930e-06,  7.552e-06,  9.556e-06,  1.183e-05,
        1.425e-05,  1.681e-05,  1.978e-05,  2.335e-05,  2.668e-05,
        3.022e-05,  3.371e-05,  3.715e-05,  3.967e-05,  4.060e-05,
        4.010e-05,  3.809e-05,  3.491e-05,  3.155e-05,  2.848e-05,
        2.678e-05,  2.660e-05,  2.811e-05,  3.071e-05,  3.294e-05,
        3.459e-05,  3.569e-05,  3.560e-05,  3.434e-05,  3.186e-05,
        2.916e-05,  2.622e-05,  2.275e-05,  1.918e-05,  1.620e-05,
        1.373e-05,  1.182e-05,  1.006e-05,  8.556e-06,  7.260e-06,
        6.107e-06,  5.034e-06,  4.211e-06,  3.426e-06,  2.865e-06,
        2.446e-06,  1.998e-06,  1.628e-06,  1.242e-06,  1.005e-06,
        7.853e-07,  6.210e-07,  5.071e-07,  4.156e-07,  3.548e-07,
        2.825e-07,  2.261e-07,  1.916e-07,  1.510e-07,  1.279e-07,
        1.059e-07,  9.140e-08,  7.707e-08,  6.170e-08,  5.311e-08,
        4.263e-08,  3.518e-08,  2.961e-08,  2.457e-08,  2.119e-08,
        1.712e-08,  1.439e-08,  1.201e-08,  1.003e-08,  8.564e-09,
        7.199e-09,  6.184e-09,  5.206e-09,  4.376e-09,  3.708e-09,
        3.157e-09,  2.725e-09,  2.361e-09,  2.074e-09,  1.797e-09,
        1.562e-09,  1.364e-09,  1.196e-09,  1.042e-09,  8.862e-10,
        7.648e-10,  6.544e-10,  5.609e-10,  4.791e-10,  4.108e-10,
        3.531e-10,  3.038e-10,  2.618e-10,  2.268e-10,  1.969e-10,
        1.715e-10,  1.496e-10,  1.308e-10,  1.147e-10,  1.008e-10,
        8.894e-11,  7.885e-11,  7.031e-11,  6.355e-11,  5.854e-11,
        5.534e-11,  5.466e-11,  5.725e-11,  6.447e-11,  7.943e-11,
        1.038e-10,  1.437e-10,  2.040e-10,  2.901e-10,  4.051e-10,
        5.556e-10,  7.314e-10,  9.291e-10,  1.134e-09,  1.321e-09,
        1.482e-09,  1.596e-09,  1.669e-09,  1.715e-09,  1.762e-09,
        1.817e-09,  1.828e-09,  1.848e-09,  1.873e-09,  1.902e-09,
        1.894e-09,  1.864e-09,  1.841e-09,  1.797e-09,  1.704e-09,
        1.559e-09,  1.382e-09,  1.187e-09,  1.001e-09,  8.468e-10,
        7.265e-10,  6.521e-10,  6.381e-10,  6.660e-10,  7.637e-10,
        9.705e-10,  1.368e-09,  1.856e-09,  2.656e-09,  3.954e-09,
        5.960e-09,  8.720e-09,  1.247e-08,  1.781e-08,  2.491e-08,
        3.311e-08,  4.272e-08,  5.205e-08,  6.268e-08,  7.337e-08,
        8.277e-08,  9.185e-08,  1.004e-07,  1.091e-07,  1.159e-07,
        1.188e-07,  1.175e-07,  1.124e-07,  1.033e-07,  9.381e-08,
        8.501e-08,  7.956e-08,  7.894e-08,  8.331e-08,  9.102e-08,
        9.836e-08,  1.035e-07,  1.064e-07,  1.060e-07,  1.032e-07,
        9.808e-08,  9.139e-08,  8.442e-08,  7.641e-08,  6.881e-08,
        6.161e-08,  5.404e-08,  4.804e-08,  4.446e-08,  4.328e-08,
        4.259e-08,  4.421e-08,  4.673e-08,  4.985e-08,  5.335e-08,
        5.796e-08,  6.542e-08,  7.714e-08,  8.827e-08,  1.040e-07,
        1.238e-07,  1.499e-07,  1.829e-07,  2.222e-07,  2.689e-07,
        3.303e-07,  3.981e-07,  4.840e-07,  5.910e-07,  7.363e-07,
        9.087e-07,  1.139e-06,  1.455e-06,  1.866e-06,  2.440e-06,
        3.115e-06,  3.941e-06,  4.891e-06,  5.992e-06,  7.111e-06,
        8.296e-06,  9.210e-06,  9.987e-06,  1.044e-05,  1.073e-05,
        1.092e-05,  1.106e-05,  1.138e-05,  1.171e-05,  1.186e-05,
        1.186e-05,  1.179e-05,  1.166e-05,  1.151e-05,  1.160e-05,
        1.197e-05,  1.241e-05,  1.268e-05,  1.260e-05,  1.184e-05,
        1.063e-05,  9.204e-06,  7.584e-06,  6.053e-06,  4.482e-06,
        3.252e-06,  2.337e-06,  1.662e-06,  1.180e-06,  8.150e-07,
        5.950e-07,  4.354e-07,  3.302e-07,  2.494e-07,  1.930e-07,
        1.545e-07,  1.250e-07,  1.039e-07,  8.602e-08,  7.127e-08,
        5.897e-08,  4.838e-08,  4.018e-08,  3.280e-08,  2.720e-08,
        2.307e-08,  1.972e-08,  1.654e-08,  1.421e-08,  1.174e-08,
        1.004e-08,  8.739e-09,  7.358e-09,  6.242e-09,  5.303e-09,
        4.567e-09,  3.940e-09,  3.375e-09,  2.864e-09,  2.422e-09,
        2.057e-09,  1.750e-09,  1.505e-09,  1.294e-09,  1.101e-09,
        9.401e-10,  8.018e-10,  6.903e-10,  5.965e-10,  5.087e-10,
        4.364e-10,  3.759e-10,  3.247e-10,  2.809e-10,  2.438e-10,
        2.123e-10,  1.853e-10,  1.622e-10,  1.426e-10,  1.260e-10,
        1.125e-10,  1.022e-10,  9.582e-11,  9.388e-11,  9.801e-11,
        1.080e-10,  1.276e-10,  1.551e-10,  1.903e-10,  2.291e-10,
        2.724e-10,  3.117e-10,  3.400e-10,  3.562e-10,  3.625e-10,
        3.619e-10,  3.429e-10,  3.221e-10,  2.943e-10,  2.645e-10,
        2.338e-10,  2.062e-10,  1.901e-10,  1.814e-10,  1.827e-10,
        1.906e-10,  1.984e-10,  2.040e-10,  2.068e-10,  2.075e-10,
        2.018e-10,  1.959e-10,  1.897e-10,  1.852e-10,  1.791e-10,
        1.696e-10,  1.634e-10,  1.598e-10,  1.561e-10,  1.518e-10,
        1.443e-10,  1.377e-10,  1.346e-10,  1.342e-10,  1.375e-10,
        1.525e-10,  1.767e-10,  2.108e-10,  2.524e-10,  2.981e-10,
        3.477e-10,  4.262e-10,  5.326e-10,  6.646e-10,  8.321e-10,
        1.069e-09,  1.386e-09,  1.743e-09,  2.216e-09,  2.808e-09,
        3.585e-09,  4.552e-09,  5.907e-09,  7.611e-09,  9.774e-09,
        1.255e-08,  1.666e-08,  2.279e-08,  3.221e-08,  4.531e-08,
        6.400e-08,  9.187e-08,  1.295e-07,  1.825e-07,  2.431e-07,
        3.181e-07,  4.009e-07,  4.941e-07,  5.880e-07,  6.623e-07,
        7.155e-07,  7.451e-07,  7.594e-07,  7.541e-07,  7.467e-07,
        7.527e-07,  7.935e-07,  8.461e-07,  8.954e-07,  9.364e-07,
        9.843e-07,  1.024e-06,  1.050e-06,  1.059e-06,  1.074e-06,
        1.072e-06,  1.043e-06,  9.789e-07,  8.803e-07,  7.662e-07,
        6.378e-07,  5.133e-07,  3.958e-07,  2.914e-07,  2.144e-07,
        1.570e-07,  1.140e-07,  8.470e-08,  6.200e-08,  4.657e-08,
        3.559e-08,  2.813e-08,  2.222e-08,  1.769e-08,  1.391e-08,
        1.125e-08,  9.186e-09,  7.704e-09,  6.447e-09,  5.381e-09,
        4.442e-09,  3.669e-09,  3.057e-09,  2.564e-09,  2.153e-09,
        1.784e-09,  1.499e-09,  1.281e-09,  1.082e-09,  9.304e-10,
        8.169e-10,  6.856e-10,  5.866e-10,  5.043e-10,  4.336e-10,
        3.731e-10,  3.175e-10,  2.745e-10,  2.374e-10,  2.007e-10,
        1.737e-10,  1.508e-10,  1.302e-10,  1.130e-10,  9.672e-11,
        8.375e-11,  7.265e-11,  6.244e-11,  5.343e-11,  4.654e-11,
        3.975e-11,  3.488e-11,  3.097e-11,  2.834e-11,  2.649e-11,
        2.519e-11,  2.462e-11,  2.443e-11,  2.440e-11,  2.398e-11,
        2.306e-11,  2.183e-11,  2.021e-11,  1.821e-11,  1.599e-11,
        1.403e-11,  1.196e-11,  1.023e-11,  8.728e-12,  7.606e-12,
        6.941e-12,  6.545e-12,  6.484e-12,  6.600e-12,  6.718e-12,
        6.785e-12,  6.746e-12,  6.724e-12,  6.764e-12,  6.995e-12,
        7.144e-12,  7.320e-12,  7.330e-12,  7.208e-12,  6.789e-12,
        6.090e-12,  5.337e-12,  4.620e-12,  4.037e-12,  3.574e-12,
        3.311e-12,  3.346e-12,  3.566e-12,  3.836e-12,  4.076e-12,
        4.351e-12,  4.691e-12,  5.114e-12,  5.427e-12,  6.167e-12,
        7.436e-12,  8.842e-12,  1.038e-11,  1.249e-11,  1.540e-11,
        1.915e-11,  2.480e-11,  3.256e-11,  4.339e-11,  5.611e-11,
        7.519e-11,  1.037e-10,  1.409e-10,  1.883e-10,  2.503e-10,
        3.380e-10,  4.468e-10,  5.801e-10,  7.335e-10,  8.980e-10,
        1.110e-09,  1.363e-09,  1.677e-09,  2.104e-09,  2.681e-09,
        3.531e-09,  4.621e-09,  6.106e-09,  8.154e-09,  1.046e-08,
        1.312e-08,  1.607e-08,  1.948e-08,  2.266e-08,  2.495e-08,
        2.655e-08,  2.739e-08,  2.739e-08,  2.662e-08,  2.589e-08,
        2.590e-08,  2.664e-08,  2.833e-08,  3.023e-08,  3.305e-08,
        3.558e-08,  3.793e-08,  3.961e-08,  4.056e-08,  4.102e-08,
        4.025e-08,  3.917e-08,  3.706e-08,  3.493e-08,  3.249e-08,
        3.096e-08,  3.011e-08,  3.111e-08,  3.395e-08,  3.958e-08,
        4.875e-08,  6.066e-08,  7.915e-08,  1.011e-07,  1.300e-07,
        1.622e-07,  2.003e-07,  2.448e-07,  2.863e-07,  3.317e-07,
        3.655e-07,  3.960e-07,  4.098e-07,  4.168e-07,  4.198e-07,
        4.207e-07,  4.289e-07,  4.384e-07,  4.471e-07,  4.524e-07,
        4.574e-07,  4.633e-07,  4.785e-07,  5.028e-07,  5.371e-07,
        5.727e-07,  5.955e-07,  5.998e-07,  5.669e-07,  5.082e-07,
        4.397e-07,  3.596e-07,  2.814e-07,  2.074e-07,  1.486e-07,
        1.057e-07,  7.250e-08,  4.946e-08,  3.430e-08,  2.447e-08,
        1.793e-08,  1.375e-08,  1.096e-08,  9.091e-09,  7.709e-09,
        6.631e-09,  5.714e-09,  4.886e-09,  4.205e-09,  3.575e-09,
        3.070e-09,  2.631e-09,  2.284e-09,  2.002e-09,  1.745e-09,
        1.509e-09,  1.284e-09,  1.084e-09,  9.163e-10,  7.663e-10,
        6.346e-10,  5.283e-10,  4.354e-10,  3.590e-10,  2.982e-10,
        2.455e-10,  2.033e-10,  1.696e-10,  1.432e-10,  1.211e-10,
        1.020e-10,  8.702e-11,  7.380e-11,  6.293e-11,  5.343e-11,
        4.532e-11,  3.907e-11,  3.365e-11,  2.945e-11,  2.558e-11,
        2.192e-11,  1.895e-11,  1.636e-11,  1.420e-11,  1.228e-11,
        1.063e-11,  9.348e-12,  8.200e-12,  7.231e-12,  6.430e-12,
        5.702e-12,  5.052e-12,  4.469e-12,  4.000e-12,  3.679e-12,
        3.387e-12,  3.197e-12,  3.158e-12,  3.327e-12,  3.675e-12,
        4.292e-12,  5.437e-12,  7.197e-12,  1.008e-11,  1.437e-11,
        2.035e-11,  2.905e-11,  4.062e-11,  5.528e-11,  7.177e-11,
        9.064e-11,  1.109e-10,  1.297e-10,  1.473e-10,  1.652e-10,
        1.851e-10,  2.079e-10,  2.313e-10,  2.619e-10,  2.958e-10,
        3.352e-10,  3.796e-10,  4.295e-10,  4.923e-10,  5.490e-10,
        5.998e-10,  6.388e-10,  6.645e-10,  6.712e-10,  6.549e-10,
        6.380e-10,  6.255e-10,  6.253e-10,  6.459e-10,  6.977e-10,
        7.590e-10,  8.242e-10,  8.920e-10,  9.403e-10,  9.701e-10,
        9.483e-10,  9.135e-10,  8.617e-10,  7.921e-10,  7.168e-10,
        6.382e-10,  5.677e-10,  5.045e-10,  4.572e-10,  4.312e-10,
        4.145e-10,  4.192e-10,  4.541e-10,  5.368e-10,  6.771e-10,
        8.962e-10,  1.210e-09,  1.659e-09,  2.330e-09,  3.249e-09,
        4.495e-09,  5.923e-09,  7.642e-09,  9.607e-09,  1.178e-08,
        1.399e-08,  1.584e-08,  1.730e-08,  1.816e-08,  1.870e-08,
        1.868e-08,  1.870e-08,  1.884e-08,  1.990e-08,  2.150e-08,
        2.258e-08,  2.364e-08,  2.473e-08,  2.602e-08,  2.689e-08,
        2.731e-08,  2.816e-08,  2.859e-08,  2.839e-08,  2.703e-08,
        2.451e-08,  2.149e-08,  1.787e-08,  1.449e-08,  1.111e-08,
        8.282e-09,  6.121e-09,  4.494e-09,  3.367e-09,  2.487e-09,
        1.885e-09,  1.503e-09,  1.249e-09,  1.074e-09,  9.427e-10,
        8.439e-10,  7.563e-10,  6.772e-10,  6.002e-10,  5.254e-10,
        4.588e-10,  3.977e-10,  3.449e-10,  3.003e-10,  2.624e-10,
        2.335e-10,  2.040e-10,  1.771e-10,  1.534e-10,  1.296e-10,
        1.097e-10,  9.173e-11,  7.730e-11,  6.547e-11,  5.191e-11,
        4.198e-11,  3.361e-11,  2.732e-11,  2.244e-11,  1.791e-11,
        1.509e-11,  1.243e-11,  1.035e-11,  8.969e-12,  7.394e-12,
        6.323e-12,  5.282e-12,  4.543e-12,  3.752e-12,  3.140e-12,
        2.600e-12,  2.194e-12,  1.825e-12,  1.511e-12,  1.245e-12,
        1.024e-12,  8.539e-13,  7.227e-13,  6.102e-13,  5.189e-13,
        4.430e-13,  3.774e-13,  3.236e-13,  2.800e-13,  2.444e-13,
        2.156e-13,  1.932e-13,  1.775e-13,  1.695e-13,  1.672e-13,
        1.704e-13,  1.825e-13,  2.087e-13,  2.614e-13,  3.377e-13,
        4.817e-13,  6.989e-13,  1.062e-12,  1.562e-12,  2.288e-12,
        3.295e-12,  4.550e-12,  5.965e-12,  7.546e-12,  9.395e-12,
        1.103e-11,  1.228e-11,  1.318e-11,  1.380e-11,  1.421e-11,
        1.390e-11,  1.358e-11,  1.336e-11,  1.342e-11,  1.356e-11,
        1.424e-11,  1.552e-11,  1.730e-11,  1.951e-11,  2.128e-11,
        2.249e-11,  2.277e-11,  2.226e-11,  2.111e-11,  1.922e-11,
        1.775e-11,  1.661e-11,  1.547e-11,  1.446e-11,  1.323e-11,
        1.210e-11,  1.054e-11,  9.283e-12,  8.671e-12,  8.670e-12,
        9.429e-12,  1.062e-11,  1.255e-11,  1.506e-11,  1.818e-11,
        2.260e-11,  2.831e-11,  3.723e-11,  5.092e-11,  6.968e-11,
        9.826e-11,  1.349e-10,  1.870e-10,  2.580e-10,  3.430e-10,
        4.424e-10,  5.521e-10,  6.812e-10,  8.064e-10,  9.109e-10,
        9.839e-10,  1.028e-09,  1.044e-09,  1.029e-09,  1.005e-09,
        1.002e-09,  1.038e-09,  1.122e-09,  1.233e-09,  1.372e-09,
        1.524e-09,  1.665e-09,  1.804e-09,  1.908e-09,  2.015e-09,
        2.117e-09,  2.219e-09,  2.336e-09,  2.531e-09,  2.805e-09,
        3.189e-09,  3.617e-09,  4.208e-09,  4.911e-09,  5.619e-09,
        6.469e-09,  7.188e-09,  7.957e-09,  8.503e-09,  9.028e-09,
        9.571e-09,  9.990e-09,  1.055e-08,  1.102e-08,  1.132e-08,
        1.141e-08,  1.145e-08,  1.145e-08,  1.176e-08,  1.224e-08,
        1.304e-08,  1.388e-08,  1.445e-08,  1.453e-08,  1.368e-08,
        1.220e-08,  1.042e-08,  8.404e-09,  6.403e-09,  4.643e-09,
        3.325e-09,  2.335e-09,  1.638e-09,  1.190e-09,  9.161e-10,
        7.412e-10,  6.226e-10,  5.516e-10,  5.068e-10,  4.831e-10,
        4.856e-10,  5.162e-10,  5.785e-10,  6.539e-10,  7.485e-10,
        8.565e-10,  9.534e-10,  1.052e-09,  1.115e-09,  1.173e-09,
        1.203e-09,  1.224e-09,  1.243e-09,  1.248e-09,  1.261e-09,
        1.265e-09,  1.250e-09,  1.217e-09,  1.176e-09,  1.145e-09,
        1.153e-09,  1.199e-09,  1.278e-09,  1.366e-09,  1.426e-09,
        1.444e-09,  1.365e-09,  1.224e-09,  1.051e-09,  8.539e-10,
        6.564e-10,  4.751e-10,  3.404e-10,  2.377e-10,  1.631e-10,
        1.114e-10,  7.870e-11,  5.793e-11,  4.284e-11,  3.300e-11,
        2.620e-11,  2.152e-11,  1.777e-11,  1.496e-11,  1.242e-11,
        1.037e-11,  8.725e-12,  7.004e-12,  5.718e-12,  4.769e-12,
        3.952e-12,  3.336e-12,  2.712e-12,  2.213e-12,  1.803e-12,
        1.492e-12,  1.236e-12,  1.006e-12,  8.384e-13,  7.063e-13,
        5.879e-13,  4.930e-13,  4.171e-13,  3.569e-13,  3.083e-13,
        2.688e-13,  2.333e-13,  2.035e-13,  1.820e-13,  1.682e-13,
        1.635e-13,  1.674e-13,  1.769e-13,  2.022e-13,  2.485e-13,
        3.127e-13,  4.250e-13,  5.928e-13,  8.514e-13,  1.236e-12,
        1.701e-12,  2.392e-12,  3.231e-12,  4.350e-12,  5.559e-12,
        6.915e-12,  8.519e-12,  1.013e-11,  1.146e-11,  1.240e-11,
        1.305e-11,  1.333e-11,  1.318e-11,  1.263e-11,  1.238e-11,
        1.244e-11,  1.305e-11,  1.432e-11,  1.623e-11,  1.846e-11,
        2.090e-11,  2.328e-11,  2.526e-11,  2.637e-11,  2.702e-11,
        2.794e-11,  2.889e-11,  2.989e-11,  3.231e-11,  3.680e-11,
        4.375e-11,  5.504e-11,  7.159e-11,  9.502e-11,  1.279e-10,
        1.645e-10,  2.098e-10,  2.618e-10,  3.189e-10,  3.790e-10,
        4.303e-10,  4.753e-10,  5.027e-10,  5.221e-10,  5.293e-10,
        5.346e-10,  5.467e-10,  5.796e-10,  6.200e-10,  6.454e-10,
        6.705e-10,  6.925e-10,  7.233e-10,  7.350e-10,  7.538e-10,
        7.861e-10,  8.077e-10,  8.132e-10,  7.749e-10,  7.036e-10,
        6.143e-10,  5.093e-10,  4.089e-10,  3.092e-10,  2.299e-10,
        1.705e-10,  1.277e-10,  9.723e-11,  7.533e-11,  6.126e-11,
        5.154e-11,  4.428e-11,  3.913e-11,  3.521e-11,  3.297e-11,
        3.275e-11,  3.460e-11,  3.798e-11,  4.251e-11,  4.745e-11,
        5.232e-11,  5.606e-11,  5.820e-11,  5.880e-11,  5.790e-11,
        5.661e-11,  5.491e-11,  5.366e-11,  5.341e-11,  5.353e-11,
        5.336e-11,  5.293e-11,  5.248e-11,  5.235e-11,  5.208e-11,
        5.322e-11,  5.521e-11,  5.725e-11,  5.827e-11,  5.685e-11,
        5.245e-11,  4.612e-11,  3.884e-11,  3.129e-11,  2.404e-11,
        1.732e-11,  1.223e-11,  8.574e-12,  5.888e-12,  3.986e-12,
        2.732e-12,  1.948e-12,  1.414e-12,  1.061e-12,  8.298e-13,
        6.612e-13,  5.413e-13,  4.472e-13,  3.772e-13,  3.181e-13,
        2.645e-13,  2.171e-13,  1.778e-13,  1.464e-13,  1.183e-13,
        9.637e-14,  7.991e-14,  6.668e-14,  5.570e-14,  4.663e-14,
        3.848e-14,  3.233e-14,  2.706e-14,  2.284e-14,  1.944e-14,
        1.664e-14,  1.430e-14,  1.233e-14,  1.066e-14,  9.234e-15,
        8.023e-15,  6.993e-15,  6.119e-15,  5.384e-15,  4.774e-15,
        4.283e-15,  3.916e-15,  3.695e-15,  3.682e-15,  4.004e-15,
        4.912e-15,  6.853e-15,  1.056e-14,  1.712e-14,  2.804e-14,
        4.516e-14,  7.113e-14,  1.084e-13,  1.426e-13,  1.734e-13,
        1.978e-13,  2.194e-13,  2.388e-13,  2.489e-13,  2.626e-13,
        2.865e-13,  3.105e-13,  3.387e-13,  3.652e-13,  3.984e-13,
        4.398e-13,  4.906e-13,  5.550e-13,  6.517e-13,  7.813e-13,
        9.272e-13,  1.164e-12,  1.434e-12,  1.849e-12,  2.524e-12,
        3.328e-12,  4.523e-12,  6.108e-12,  8.207e-12,  1.122e-11,
        1.477e-11,  1.900e-11,  2.412e-11,  2.984e-11,  3.680e-11,
        4.353e-11,  4.963e-11,  5.478e-11,  5.903e-11,  6.233e-11,
        6.483e-11,  6.904e-11,  7.569e-11,  8.719e-11,  1.048e-10,
        1.278e-10,  1.557e-10,  1.869e-10,  2.218e-10,  2.610e-10,
        2.975e-10,  3.371e-10,  3.746e-10,  4.065e-10,  4.336e-10,
        4.503e-10,  4.701e-10,  4.800e-10,  4.917e-10,  5.038e-10,
        5.128e-10,  5.143e-10,  5.071e-10,  5.019e-10,  5.025e-10,
        5.183e-10,  5.496e-10,  5.877e-10,  6.235e-10,  6.420e-10,
        6.234e-10,  5.698e-10,  4.916e-10,  4.022e-10,  3.126e-10,
        2.282e-10,  1.639e-10,  1.142e-10,  7.919e-11,  5.690e-11,
        4.313e-11,  3.413e-11,  2.807e-11,  2.410e-11,  2.166e-11,
        2.024e-11,  1.946e-11,  1.929e-11,  1.963e-11,  2.035e-11,
        2.162e-11,  2.305e-11,  2.493e-11,  2.748e-11,  3.048e-11,
        3.413e-11,  3.754e-11,  4.155e-11,  4.635e-11,  5.110e-11,
        5.734e-11,  6.338e-11,  6.990e-11,  7.611e-11,  8.125e-11,
        8.654e-11,  8.951e-11,  9.182e-11,  9.310e-11,  9.273e-11,
        9.094e-11,  8.849e-11,  8.662e-11,  8.670e-11,  8.972e-11,
        9.566e-11,  1.025e-10,  1.083e-10,  1.111e-10,  1.074e-10,
        9.771e-11,  8.468e-11,  6.958e-11,  5.470e-11,  4.040e-11,
        2.940e-11,  2.075e-11,  1.442e-11,  1.010e-11,  7.281e-12,
        5.409e-12,  4.138e-12,  3.304e-12,  2.784e-12,  2.473e-12,
        2.273e-12,  2.186e-12,  2.118e-12,  2.066e-12,  1.958e-12,
        1.818e-12,  1.675e-12,  1.509e-12,  1.349e-12,  1.171e-12,
        9.838e-13,  8.213e-13,  6.765e-13,  5.378e-13,  4.161e-13,
        3.119e-13,  2.279e-13,  1.637e-13,  1.152e-13,  8.112e-14,
        5.919e-14,  4.470e-14,  3.492e-14,  2.811e-14,  2.319e-14,
        1.948e-14,  1.660e-14,  1.432e-14,  1.251e-14,  1.109e-14,
        1.006e-14,  9.450e-15,  9.384e-15,  1.012e-14,  1.216e-14,
        1.636e-14,  2.305e-14,  3.488e-14,  5.572e-14,  8.479e-14,
        1.265e-13,  1.905e-13,  2.730e-13,  3.809e-13,  4.955e-13,
        6.303e-13,  7.861e-13,  9.427e-13,  1.097e-12,  1.212e-12,
        1.328e-12,  1.415e-12,  1.463e-12,  1.495e-12,  1.571e-12,
        1.731e-12,  1.981e-12,  2.387e-12,  2.930e-12,  3.642e-12,
        4.584e-12,  5.822e-12,  7.278e-12,  9.193e-12,  1.135e-11,
        1.382e-11,  1.662e-11,  1.958e-11,  2.286e-11,  2.559e-11,
        2.805e-11,  2.988e-11,  3.106e-11,  3.182e-11,  3.200e-11,
        3.258e-11,  3.362e-11,  3.558e-11,  3.688e-11,  3.800e-11,
        3.929e-11,  4.062e-11,  4.186e-11,  4.293e-11,  4.480e-11,
        4.643e-11,  4.704e-11,  4.571e-11,  4.206e-11,  3.715e-11,
        3.131e-11,  2.541e-11,  1.978e-11,  1.508e-11,  1.146e-11,
        8.700e-12,  6.603e-12,  5.162e-12,  4.157e-12,  3.408e-12,
        2.829e-12,  2.405e-12,  2.071e-12,  1.826e-12,  1.648e-12,
        1.542e-12,  1.489e-12,  1.485e-12,  1.493e-12,  1.545e-12,
        1.637e-12,  1.814e-12,  2.061e-12,  2.312e-12,  2.651e-12,
        3.030e-12,  3.460e-12,  3.901e-12,  4.306e-12,  4.721e-12,
        5.008e-12,  5.281e-12,  5.541e-12,  5.791e-12,  6.115e-12,
        6.442e-12,  6.680e-12,  6.791e-12,  6.831e-12,  6.839e-12,
        6.946e-12,  7.128e-12,  7.537e-12,  8.036e-12,  8.392e-12,
        8.526e-12,  8.110e-12,  7.325e-12,  6.329e-12,  5.183e-12,
        4.081e-12,  2.985e-12,  2.141e-12,  1.492e-12,  1.015e-12,
        6.684e-13,  4.414e-13,  2.987e-13,  2.038e-13,  1.391e-13,
        9.860e-14,  7.240e-14,  5.493e-14,  4.288e-14,  3.427e-14,
        2.787e-14,  2.296e-14,  1.909e-14,  1.598e-14,  1.344e-14,
        1.135e-14,  9.616e-15,  8.169e-15,  6.957e-15,  5.938e-15,
        5.080e-15,  4.353e-15,  3.738e-15,  3.217e-15,  2.773e-15,
        2.397e-15,  2.077e-15,  1.805e-15,  1.575e-15,  1.382e-15,
        1.221e-15,  1.090e-15,  9.855e-16,  9.068e-16,  8.537e-16,
        8.270e-16,  8.290e-16,  8.634e-16,  9.359e-16,  1.055e-15,
        1.233e-15,  1.486e-15,  1.839e-15,  2.326e-15,  2.998e-15,
        3.934e-15,  5.256e-15,  7.164e-15,  9.984e-15,  1.427e-14,
        2.099e-14,  3.196e-14,  5.121e-14,  7.908e-14,  1.131e-13,
        1.602e-13,  2.239e-13,  3.075e-13,  4.134e-13,  5.749e-13,
        7.886e-13,  1.071e-12,  1.464e-12,  2.032e-12,  2.800e-12,
        3.732e-12,  4.996e-12,  6.483e-12,  8.143e-12,  1.006e-11,
        1.238e-11,  1.484e-11,  1.744e-11,  2.020e-11,  2.274e-11,
        2.562e-11,  2.848e-11,  3.191e-11,  3.617e-11,  4.081e-11,
        4.577e-11,  4.937e-11,  5.204e-11,  5.401e-11,  5.462e-11,
        5.507e-11,  5.510e-11,  5.605e-11,  5.686e-11,  5.739e-11,
        5.766e-11,  5.740e-11,  5.754e-11,  5.761e-11,  5.777e-11,
        5.712e-11,  5.510e-11,  5.088e-11,  4.438e-11,  3.728e-11,
        2.994e-11,  2.305e-11,  1.715e-11,  1.256e-11,  9.208e-12,
        6.745e-12,  5.014e-12,  3.785e-12,  2.900e-12,  2.239e-12,
        1.757e-12,  1.414e-12,  1.142e-12,  9.482e-13,  8.010e-13,
        6.961e-13,  6.253e-13,  5.735e-13,  5.433e-13,  5.352e-13,
        5.493e-13,  5.706e-13,  6.068e-13,  6.531e-13,  7.109e-13,
        7.767e-13,  8.590e-13,  9.792e-13,  1.142e-12,  1.371e-12,
        1.650e-12,  1.957e-12,  2.302e-12,  2.705e-12,  3.145e-12,
        3.608e-12,  4.071e-12,  4.602e-12,  5.133e-12,  5.572e-12,
        5.987e-12,  6.248e-12,  6.533e-12,  6.757e-12,  6.935e-12,
        7.224e-12,  7.422e-12,  7.538e-12,  7.547e-12,  7.495e-12,
        7.543e-12,  7.725e-12,  8.139e-12,  8.627e-12,  9.146e-12,
        9.443e-12,  9.318e-12,  8.649e-12,  7.512e-12,  6.261e-12,
        4.915e-12,  3.647e-12,  2.597e-12,  1.785e-12,  1.242e-12,
        8.660e-13,  6.207e-13,  4.610e-13,  3.444e-13,  2.634e-13,
        2.100e-13,  1.725e-13,  1.455e-13,  1.237e-13,  1.085e-13,
        9.513e-14,  7.978e-14,  6.603e-14,  5.288e-14,  4.084e-14,
        2.952e-14,  2.157e-14,  1.593e-14,  1.199e-14,  9.267e-15,
        7.365e-15,  6.004e-15,  4.995e-15,  4.218e-15,  3.601e-15,
        3.101e-15,  2.692e-15,  2.360e-15,  2.094e-15,  1.891e-15,
        1.755e-15,  1.699e-15,  1.755e-15,  1.987e-15,  2.506e-15,
        3.506e-15,  5.289e-15,  8.311e-15,  1.325e-14,  2.129e-14,
        3.237e-14,  4.595e-14,  6.441e-14,  8.433e-14,  1.074e-13,
        1.383e-13,  1.762e-13,  2.281e-13,  2.831e-13,  3.523e-13,
        4.380e-13,  5.304e-13,  6.290e-13,  7.142e-13,  8.032e-13,
        8.934e-13,  9.888e-13,  1.109e-12,  1.261e-12,  1.462e-12,
        1.740e-12,  2.099e-12,  2.535e-12,  3.008e-12,  3.462e-12,
        3.856e-12,  4.098e-12,  4.239e-12,  4.234e-12,  4.132e-12,
        3.986e-12,  3.866e-12,  3.829e-12,  3.742e-12,  3.705e-12,
        3.694e-12,  3.765e-12,  3.849e-12,  3.929e-12,  4.056e-12,
        4.092e-12,  4.047e-12,  3.792e-12,  3.407e-12,  2.953e-12,
        2.429e-12,  1.931e-12,  1.460e-12,  1.099e-12,  8.199e-13,
        6.077e-13,  4.449e-13,  3.359e-13,  2.524e-13,  1.881e-13,
        1.391e-13,  1.020e-13,  7.544e-14,  5.555e-14,  4.220e-14,
        3.321e-14,  2.686e-14,  2.212e-14,  1.780e-14,  1.369e-14,
        1.094e-14,  9.130e-15,  8.101e-15,  7.828e-15,  8.393e-15,
        1.012e-14,  1.259e-14,  1.538e-14,  1.961e-14,  2.619e-14,
        3.679e-14,  5.049e-14,  6.917e-14,  8.880e-14,  1.115e-13,
        1.373e-13,  1.619e-13,  1.878e-13,  2.111e-13,  2.330e-13,
        2.503e-13,  2.613e-13,  2.743e-13,  2.826e-13,  2.976e-13,
        3.162e-13,  3.360e-13,  3.491e-13,  3.541e-13,  3.595e-13,
        3.608e-13,  3.709e-13,  3.869e-13,  4.120e-13,  4.366e-13,
        4.504e-13,  4.379e-13,  3.955e-13,  3.385e-13,  2.741e-13,
        2.089e-13,  1.427e-13,  9.294e-14,  5.775e-14,  3.565e-14,
        2.210e-14,  1.398e-14,  9.194e-15,  6.363e-15,  4.644e-15,
        3.550e-15,  2.808e-15,  2.274e-15,  1.871e-15,  1.557e-15,
        1.308e-15,  1.108e-15,  9.488e-16,  8.222e-16,  7.238e-16,
        6.506e-16,  6.008e-16,  5.742e-16,  5.724e-16,  5.991e-16,
        6.625e-16,  7.775e-16,  9.734e-16,  1.306e-15,  1.880e-15,
        2.879e-15,  4.616e-15,  7.579e-15,  1.248e-14,  2.030e-14,
        3.244e-14,  5.171e-14,  7.394e-14,  9.676e-14,  1.199e-13,
        1.467e-13,  1.737e-13,  2.020e-13,  2.425e-13,  3.016e-13,
        3.700e-13,  4.617e-13,  5.949e-13,  7.473e-13,  9.378e-13,
        1.191e-12,  1.481e-12,  1.813e-12,  2.232e-12,  2.722e-12,
        3.254e-12,  3.845e-12,  4.458e-12,  5.048e-12,  5.511e-12,
        5.898e-12,  6.204e-12,  6.293e-12,  6.386e-12,  6.467e-12,
        6.507e-12,  6.466e-12,  6.443e-12,  6.598e-12,  6.873e-12,
        7.300e-12,  7.816e-12,  8.368e-12,  8.643e-12,  8.466e-12,
        7.871e-12,  6.853e-12,  5.714e-12,  4.482e-12,  3.392e-12,
        2.613e-12,  2.008e-12,  1.562e-12,  1.228e-12,  9.888e-13,
        7.646e-13,  5.769e-13,  4.368e-13,  3.324e-13,  2.508e-13,
        1.916e-13};


// CO2 continuum Ridgeway 1982, implementation of CKD_MT_1.00
// UNITS OF (CM**3/MOL)*1.E-20
const Numeric FCO2_ckd_mt_100_v1  =   -20.0;
const Numeric FCO2_ckd_mt_100_v2  = 10000.0;
const Numeric FCO2_ckd_mt_100_dv  =    10.0;
const int     FCO2_ckd_mt_100_npt =  1003;
const double  FCO2_ckd_mt_100[FCO2_ckd_mt_100_npt+addF77fields] = {
           0.000e0,    1.1110E-11, 1.0188E-11, 
           9.3516E-12, 1.0188E-11, 1.1110E-11, 1.2127E-11, 1.3251E-11,    // F17590
           1.4495E-11, 1.5872E-11, 1.7400E-11, 1.9097E-11, 2.0985E-11,    // F17600
           2.3087E-11, 2.5431E-11, 2.8051E-11, 3.0982E-11, 3.4268E-11,    // F17610
           3.7956E-11, 4.2105E-11, 4.6779E-11, 5.2056E-11, 5.8025E-11,    // F17620
           6.4791E-11, 7.2477E-11, 8.1226E-11, 9.1209E-11, 1.0263E-10,    // F17630
           1.1572E-10, 1.3078E-10, 1.4814E-10, 1.6821E-10, 1.9148E-10,    // F17640
           2.1857E-10, 2.5019E-10, 2.8723E-10, 3.3080E-10, 3.8223E-10,    // F17650
           4.4321E-10, 5.1583E-10, 6.0274E-10, 7.0725E-10, 8.3363E-10,    // F17660
           9.8735E-10, 1.1755E-09, 1.4074E-09, 1.6953E-09, 2.0557E-09,    // F17670
           2.5107E-09, 3.0909E-09, 3.8391E-09, 4.8165E-09, 6.1117E-09,    // F17680
           7.8550E-09, 1.0241E-08, 1.3593E-08, 1.8344E-08, 2.5408E-08,    // F17700
           3.6386E-08, 5.4251E-08, 8.4262E-08, 1.3273E-07, 2.1867E-07,    // F17710
           3.5007E-07, 6.0011E-07, 1.0797E-06, 1.8254E-06, 3.1621E-06,    // F17720
           4.0293E-06, 4.3683E-06, 4.4552E-06, 4.2684E-06, 3.9341E-06,    // F17730
           2.5972E-06, 1.5617E-06, 8.9063E-07, 5.0360E-07, 3.0616E-07,    // F17740
           1.9066E-07, 1.1904E-07, 7.6078E-08, 4.9304E-08, 3.3335E-08,    // F17750
           2.3494E-08, 1.7114E-08, 1.2742E-08, 9.6068E-09, 7.3706E-09,    // F17760
           5.7386E-09, 4.5302E-09, 3.6223E-09, 2.9309E-09, 2.4001E-09,    // F17770
           1.9927E-09, 1.6877E-09, 1.4602E-09, 1.2764E-09, 1.1317E-09,    // F17780
           1.0273E-09, 9.1943E-10, 8.0353E-10, 6.8746E-10, 5.9354E-10,    // F17790
           5.1722E-10, 4.4975E-10, 4.2350E-10, 4.2282E-10, 4.2610E-10,    // F17810
           4.5465E-10, 4.6166E-10, 4.3149E-10, 3.7615E-10, 3.1576E-10,    // F17820
           2.6490E-10, 1.9143E-10, 1.2885E-10, 9.4954E-11, 7.6499E-11,    // F17830
           6.4581E-11, 5.5923E-11, 4.9200E-11, 4.3813E-11, 3.9533E-11,    // F17840
           3.6338E-11, 3.4320E-11, 3.3329E-11, 3.2400E-11, 3.1700E-11,    // F17850
           3.1267E-11, 2.9940E-11, 2.7628E-11, 2.4496E-11, 2.1764E-11,    // F17860
           1.9306E-11, 1.7352E-11, 1.7292E-11, 1.8733E-11, 2.0224E-11,    // F17870
           2.2396E-11, 2.4225E-11, 2.4890E-11, 2.3513E-11, 2.0824E-11,    // F17880
           1.8642E-11, 1.5676E-11, 1.2882E-11, 1.1054E-11, 1.0074E-11,    // F17890
           9.6324E-12, 9.4910E-12, 9.5134E-12, 9.6427E-12, 9.8552E-12,    // F17900
           1.0140E-11, 1.0494E-11, 1.0915E-11, 1.1405E-11, 1.1965E-11,    // F17920
           1.2601E-11, 1.3316E-11, 1.4116E-11, 1.5006E-11, 1.5997E-11,    // F17930
           1.7092E-11, 1.8305E-11, 1.9641E-11, 2.1121E-11, 2.2744E-11,    // F17940
           2.4503E-11, 2.6419E-11, 2.8221E-11, 3.0609E-11, 3.3260E-11,    // F17950
           3.6247E-11, 3.9581E-11, 4.3279E-11, 4.7376E-11, 5.1932E-11,    // F17960
           5.7001E-11, 6.2654E-11, 6.8973E-11, 7.6058E-11, 8.4037E-11,    // F17970
           9.3081E-11, 1.0344E-10, 1.1547E-10, 1.2970E-10, 1.4659E-10,    // F17980
           1.6724E-10, 1.9481E-10, 2.3520E-10, 2.9424E-10, 3.6319E-10,    // F17990
           4.2279E-10, 4.8494E-10, 5.2296E-10, 5.6111E-10, 5.8935E-10,    // F18000
           6.0807E-10, 6.4204E-10, 6.8457E-10, 7.6709E-10, 8.7664E-10,    // F18010
           1.0183E-09, 1.2116E-09, 1.4874E-09, 1.8596E-09, 2.2742E-09,    // F18030
           2.7577E-09, 3.1932E-09, 3.6381E-09, 4.1207E-09, 4.6458E-09,    // F18040
           5.3065E-09, 6.0741E-09, 7.1942E-09, 8.7103E-09, 1.0713E-08,    // F18050
           1.3344E-08, 1.6831E-08, 2.1524E-08, 2.7967E-08, 3.7047E-08,    // F18060
           5.0312E-08, 7.0566E-08, 1.0275E-07, 1.5419E-07, 2.3309E-07,    // F18070
           3.4843E-07, 5.3194E-07, 8.7207E-07, 1.5075E-06, 2.7077E-06,    // F18080
           4.7125E-06, 7.1734E-06, 9.2381E-06, 1.1507E-05, 1.3737E-05,    // F18090
           1.4004E-05, 1.2679E-05, 1.0478E-05, 8.5684E-06, 6.1472E-06,    // F18100
           3.2424E-06, 1.5291E-06, 8.0390E-07, 4.6767E-07, 2.9170E-07,    // F18110
           1.9148E-07, 1.3076E-07, 9.2156E-08, 6.6652E-08, 4.9265E-08,    // F18120
           3.7094E-08, 2.8380E-08, 2.2019E-08, 1.7297E-08, 1.3738E-08,    // F18140
           1.1019E-08, 8.9178E-09, 7.2762E-09, 5.9810E-09, 4.9500E-09,    // F18150
           4.1226E-09, 3.4534E-09, 2.9082E-09, 2.4611E-09, 2.0922E-09,    // F18160
           1.7864E-09, 1.5313E-09, 1.3176E-09, 1.1379E-09, 9.8612E-10,    // F18170
           8.5741E-10, 7.4782E-10, 6.5416E-10, 5.7384E-10, 5.0471E-10,    // F18180
           4.4503E-10, 3.9334E-10, 3.4841E-10, 3.0927E-10, 2.7510E-10,    // F18190
           2.4519E-10, 2.1893E-10, 1.9587E-10, 1.7555E-10, 1.5762E-10,    // F18200
           1.4178E-10, 1.2772E-10, 1.1524E-10, 1.0414E-10, 9.4248E-11,    // F18210
           8.5421E-11, 7.7530E-11, 7.0466E-11, 6.4134E-11, 5.8450E-11,    // F18220
           5.3342E-11, 4.8746E-11, 4.4607E-11, 4.0874E-11, 3.7507E-11,    // F18230
           3.4466E-11, 3.1719E-11, 2.9237E-11, 2.6993E-11, 2.4968E-11,    // F18250
           2.3139E-11, 2.1494E-11, 2.0022E-11, 1.8709E-11, 1.7541E-11,    // F18260
           1.6533E-11, 1.5690E-11, 1.5027E-11, 1.4560E-11, 1.4169E-11,    // F18270
           1.3796E-11, 1.3553E-11, 1.3526E-11, 1.3567E-11, 1.3399E-11,    // F18280
           1.3149E-11, 1.3049E-11, 1.3078E-11, 1.3093E-11, 1.3168E-11,    // F18290
           1.3572E-11, 1.4383E-11, 1.5698E-11, 1.7658E-11, 2.0197E-11,    // F18300
           2.2845E-11, 2.5944E-11, 3.0250E-11, 3.5900E-11, 4.1482E-11,    // F18310
           4.6602E-11, 5.2453E-11, 5.9754E-11, 6.9308E-11, 8.0696E-11,    // F18320
           9.5737E-11, 1.1733E-10, 1.4793E-10, 1.9119E-10, 2.5355E-10,    // F18330
           3.4588E-10, 4.8343E-10, 6.9378E-10, 1.0212E-09, 1.4858E-09,    // F18340
           2.0906E-09, 3.0576E-09, 4.6318E-09, 7.1585E-09, 1.1259E-08,    // F18360
           1.7954E-08, 2.9760E-08, 4.6693E-08, 6.2035E-08, 7.4399E-08,    // F18370
           9.1705E-08, 9.9448E-08, 9.5181E-08, 8.3050E-08, 7.1756E-08,    // F18380
           6.6261E-08, 6.0357E-08, 6.6988E-08, 8.3419E-08, 9.8834E-08,    // F18390
           1.2385E-07, 1.3962E-07, 1.3651E-07, 1.1963E-07, 9.7731E-08,    // F18400
           8.0083E-08, 5.1660E-08, 2.5778E-08, 1.2600E-08, 6.8779E-09,    // F18410
           4.1161E-09, 2.6276E-09, 1.7595E-09, 1.2225E-09, 8.7493E-10,    // F18420
           6.4179E-10, 4.7987E-10, 3.6491E-10, 2.8191E-10, 2.2084E-10,    // F18430
           1.7507E-10, 1.4025E-10, 1.1344E-10, 9.2580E-11, 7.6170E-11,    // F18440
           6.3142E-11, 5.2694E-11, 4.4260E-11, 3.7421E-11, 3.1847E-11,    // F18450
           2.7263E-11, 2.3352E-11, 2.0081E-11, 1.7332E-11, 1.5000E-11,    // F18470
           1.2978E-11, 1.1204E-11, 9.7513E-12, 8.5300E-12, 7.4888E-12,    // F18480
           6.5947E-12, 5.8231E-12, 5.1548E-12, 4.5739E-12, 4.0675E-12,    // F18490
           3.6250E-12, 3.2371E-12, 2.8963E-12, 2.5964E-12, 2.3316E-12,    // F18500
           2.0975E-12, 1.8902E-12, 1.7061E-12, 1.5425E-12, 1.3967E-12,    // F18510
           1.2665E-12, 1.1503E-12, 1.0463E-12, 9.5319E-13, 8.6963E-13,    // F18520
           7.9461E-13, 7.2718E-13, 6.6654E-13, 6.1201E-13, 5.6296E-13,    // F18530
           5.1894E-13, 4.7969E-13, 4.4494E-13, 4.1320E-13, 3.8529E-13,    // F18540
           3.6202E-13, 3.4320E-13, 3.2546E-13, 3.0741E-13, 2.9156E-13,    // F18550
           2.7819E-13, 2.6576E-13, 2.5327E-13, 2.4319E-13, 2.3770E-13,    // F18560
           2.3645E-13, 2.3967E-13, 2.4960E-13, 2.6858E-13, 2.9679E-13,    // F18580
           3.3247E-13, 3.8487E-13, 4.7576E-13, 6.1833E-13, 8.0740E-13,    // F18590
           1.0267E-12, 1.2291E-12, 1.4710E-12, 1.7211E-12, 1.8251E-12,    // F18600
           1.8982E-12, 1.9768E-12, 2.1877E-12, 2.5008E-12, 3.0545E-12,    // F18610
           4.1513E-12, 5.7469E-12, 7.7913E-12, 1.0873E-11, 1.5538E-11,    // F18620
           2.2838E-11, 3.4153E-11, 4.9751E-11, 7.0591E-11, 1.0794E-10,    // F18630
           1.7287E-10, 2.6554E-10, 3.5250E-10, 4.1952E-10, 5.1979E-10,    // F18640
           5.7649E-10, 5.6168E-10, 5.0014E-10, 4.3670E-10, 4.0057E-10,    // F18650
           3.5169E-10, 3.7578E-10, 5.5054E-10, 8.8962E-10, 1.2940E-09,    // F18660
           1.6293E-09, 2.0553E-09, 2.3945E-09, 2.3926E-09, 2.1385E-09,    // F18670
           1.7637E-09, 1.4623E-09, 1.0150E-09, 5.5612E-10, 3.5162E-10,    // F18690
           3.4009E-10, 4.1744E-10, 5.0009E-10, 6.0748E-10, 7.3258E-10,    // F18700
           7.6553E-10, 7.2066E-10, 6.1317E-10, 5.1585E-10, 3.9136E-10,    // F18710
           2.2991E-10, 1.2590E-10, 6.9549E-11, 3.8699E-11, 2.2976E-11,    // F18720
           1.4702E-11, 9.9989E-12, 7.1233E-12, 5.2612E-12, 4.0298E-12,    // F18730
           3.2395E-12, 2.7932E-12, 2.6331E-12, 2.7835E-12, 3.3167E-12,    // F18740
           3.3581E-12, 3.3404E-12, 3.1243E-12, 2.8459E-12, 2.4092E-12,    // F18750
           1.5349E-12, 9.7039E-13, 5.8611E-13, 3.9686E-13, 2.9332E-13,    // F18760
           2.2795E-13, 1.8432E-13, 1.5287E-13, 1.2898E-13, 1.1019E-13,    // F18770
           9.5041E-14, 8.2617E-14, 7.2310E-14, 6.3711E-14, 5.6561E-14,    // F18780
           5.0763E-14, 4.6525E-14, 4.4418E-14, 4.4681E-14, 4.7199E-14,    // F18800
           5.0389E-14, 5.3620E-14, 6.0817E-14, 6.0192E-14, 5.5878E-14,    // F18810
           4.9874E-14, 4.3955E-14, 3.9854E-14, 3.1697E-14, 3.1135E-14,    // F18820
           3.4683E-14, 3.8789E-14, 4.6932E-14, 5.0213E-14, 4.7156E-14,    // F18830
           4.2130E-14, 3.5554E-14, 3.0465E-14, 1.9216E-14, 1.1378E-14,    // F18840
           8.2878E-15, 6.8260E-15, 6.0960E-15, 5.8135E-15, 5.9618E-15,    // F18850
           6.8295E-15, 9.2943E-15, 1.2572E-14, 1.4837E-14, 1.8595E-14,    // F18860
           2.1533E-14, 2.2008E-14, 2.1305E-14, 1.9743E-14, 2.0413E-14,    // F18870
           2.1131E-14, 2.5346E-14, 3.3709E-14, 4.3995E-14, 5.8911E-14,    // F18880
           7.8451E-14, 1.0537E-13, 1.4559E-13, 2.0405E-13, 2.6734E-13,    // F18890
           3.5029E-13, 4.9788E-13, 7.3207E-13, 1.0979E-12, 1.4960E-12,    // F18910
           1.7906E-12, 2.2171E-12, 2.5369E-12, 2.5873E-12, 2.3871E-12,    // F18920
           2.0730E-12, 1.9095E-12, 1.6227E-12, 1.3981E-12, 1.5228E-12,    // F18930
           2.0956E-12, 3.2493E-12, 5.2740E-12, 8.6666E-12, 1.2672E-11,    // F18940
           1.5725E-11, 1.9496E-11, 2.2858E-11, 2.2939E-11, 2.0597E-11,    // F18950
           1.7021E-11, 1.4456E-11, 1.0794E-11, 7.1327E-12, 6.5438E-12,    // F18960
           8.8057E-12, 1.2311E-11, 1.5284E-11, 1.9273E-11, 2.2796E-11,    // F18970
           2.3156E-11, 2.0914E-11, 1.7298E-11, 1.4424E-11, 1.0127E-11,    // F18980
           5.2952E-12, 2.5759E-12, 1.4304E-12, 9.4758E-13, 7.9895E-13,    // F18990
           9.1124E-13, 1.2297E-12, 1.5898E-12, 1.9056E-12, 2.3905E-12,    // F19000
           2.6695E-12, 2.6297E-12, 2.3467E-12, 2.0058E-12, 1.6773E-12,    // F19020
           1.1327E-12, 6.7331E-13, 4.0954E-13, 2.5152E-13, 1.4491E-13,    // F19030
           9.0916E-14, 6.6510E-14, 5.9022E-14, 6.4403E-14, 8.3126E-14,    // F19040
           1.2409E-13, 1.5153E-13, 1.6909E-13, 1.7938E-13, 1.9169E-13,    // F19050
           2.1173E-13, 2.1941E-13, 2.6360E-13, 3.5956E-13, 4.8369E-13,    // F19060
           5.9657E-13, 7.4062E-13, 8.9452E-13, 8.7899E-13, 8.2012E-13,    // F19070
           7.4109E-13, 6.9845E-13, 6.3130E-13, 5.6538E-13, 6.9516E-13,    // F19080
           9.9486E-13, 1.5226E-12, 2.4155E-12, 3.9119E-12, 6.3541E-12,    // F19090
           1.0075E-11, 1.5903E-11, 2.5091E-11, 3.6282E-11, 4.6076E-11,    // F19100
           5.6240E-11, 7.1126E-11, 7.0230E-11, 6.3642E-11, 5.3722E-11,    // F19110
           4.4651E-11, 3.4409E-11, 1.5287E-11, 7.2479E-12, 3.9218E-12,    // F19130
           2.3172E-12, 1.4585E-12, 9.6297E-13, 6.6017E-13, 4.6655E-13,    // F19140
           3.3814E-13, 2.5034E-13, 1.8874E-13, 1.4457E-13, 1.1228E-13,    // F19150
           8.8284E-14, 7.0188E-14, 5.6365E-14, 4.5685E-14, 3.7357E-14,    // F19160
           3.0817E-14, 2.5674E-14, 2.1679E-14, 1.8780E-14, 1.7243E-14,    // F19170
           1.6273E-14, 1.5201E-14, 1.5091E-14, 1.4725E-14, 1.3668E-14,    // F19180
           1.1940E-14, 1.0097E-14, 8.8905E-15, 7.1475E-15, 5.8080E-15,    // F19190
           5.5216E-15, 5.9338E-15, 7.1932E-15, 9.9780E-15, 1.6167E-14,    // F19200
           2.9100E-14, 5.2355E-14, 8.4889E-14, 1.1311E-13, 1.4192E-13,    // F19210
           1.7648E-13, 1.8657E-13, 1.7498E-13, 1.4877E-13, 1.2578E-13,    // F19220
           1.0051E-13, 6.7213E-14, 5.4750E-14, 7.0454E-14, 1.1351E-13,    // F19240
           1.8015E-13, 2.4825E-13, 3.0875E-13, 3.9200E-13, 4.2550E-13,    // F19250
           4.0067E-13, 3.4438E-13, 2.8204E-13, 2.2432E-13, 1.3172E-13,    // F19260
           6.2820E-14, 3.6474E-14, 2.9409E-14, 3.4164E-14, 4.8300E-14,    // F19270
           6.4140E-14, 7.7284E-14, 9.7973E-14, 1.0969E-13, 1.0580E-13,    // F19280
           9.2070E-14, 7.5008E-14, 6.1722E-14, 3.8874E-14, 1.9007E-14,    // F19290
           9.6765E-15, 5.5169E-15, 3.5254E-15, 2.5012E-15, 2.0013E-15,    // F19300
           1.8810E-15, 2.2143E-15, 3.5332E-15, 5.7552E-15, 7.3359E-15,    // F19310
           8.3292E-15, 9.9174E-15, 1.0930E-14, 1.1185E-14, 1.0884E-14,    // F19320
           1.0577E-14, 1.1048E-14, 1.1611E-14, 1.1128E-14, 1.0729E-14,    // F19330
           1.0248E-14, 1.0630E-14, 1.1793E-14, 1.3977E-14, 1.9857E-14,    // F19350
           2.9182E-14, 4.2229E-14, 6.2710E-14, 9.0717E-14, 1.2561E-13,    // F19360
           1.6951E-13, 2.2520E-13, 3.2470E-13, 4.5178E-13, 6.3104E-13,    // F19370
           8.7521E-13, 1.1073E-12, 1.3534E-12, 1.6954E-12, 1.7005E-12,    // F19380
           1.5993E-12, 1.4416E-12, 1.3280E-12, 1.2760E-12, 1.1076E-12,    // F19390
           1.2850E-12, 1.6208E-12, 1.9527E-12, 2.4941E-12, 2.5077E-12,    // F19400
           2.3156E-12, 2.0069E-12, 1.6301E-12, 1.2885E-12, 5.9863E-13,    // F19410
           2.8012E-13, 1.5065E-13, 8.8802E-14, 5.5888E-14, 3.6951E-14,    // F19420
           2.5393E-14, 1.8001E-14, 1.3093E-14, 9.7308E-15, 7.3665E-15,    // F19430
           5.6662E-15, 4.4194E-15, 3.4897E-15, 2.7857E-15, 2.2457E-15,    // F19440
           1.8264E-15, 1.4973E-15, 1.2365E-15, 1.0280E-15, 8.5996E-16,    // F19460
           7.2345E-16, 6.1182E-16, 5.1994E-16, 4.4388E-16, 3.8055E-16,    // F19470
           3.2756E-16, 2.8300E-16, 2.4537E-16, 2.1347E-16, 1.8630E-16,    // F19480
           1.6307E-16, 1.4314E-16, 1.2599E-16, 1.1117E-16, 9.8344E-17,    // F19490
           8.7197E-17, 7.7487E-17, 6.9004E-17, 6.1577E-17, 5.5060E-17,    // F19500
           4.9325E-17, 4.4271E-17, 3.9810E-17, 3.5861E-17, 3.2361E-17,    // F19510
           2.9252E-17, 2.6487E-17, 2.4023E-17, 2.1826E-17, 1.9862E-17,    // F19520
           1.8107E-17, 1.6536E-17, 1.5129E-17, 1.3869E-17, 1.2739E-17,    // F19530
           1.1726E-17, 1.0820E-17, 1.0009E-17, 9.2846E-18, 8.6398E-18,    // F19540
           8.0682E-18, 7.5641E-18, 7.1229E-18, 6.7411E-18, 6.4161E-18,    // F19550
           6.1455E-18, 5.9290E-18, 5.7662E-18, 5.6574E-18, 5.6049E-18,    // F19570
           5.6112E-18, 5.6811E-18, 5.8200E-18, 6.0364E-18, 6.3405E-18,    // F19580
           6.7450E-18, 7.2674E-18, 7.9298E-18, 8.7613E-18, 9.8010E-18,    // F19590
           1.1086E-17, 1.2686E-17, 1.4679E-17, 1.7177E-17, 2.0335E-17,    // F19600
           2.4384E-17, 2.9538E-17, 3.6416E-17, 4.5520E-17, 5.7788E-17,    // F19610
           7.4676E-17, 9.8513E-17, 1.3323E-16, 1.8570E-16, 2.6897E-16,    // F19620
           4.0958E-16, 6.6785E-16, 1.2064E-15, 2.4023E-15, 4.3240E-15,    // F19630
           6.6353E-15, 8.6393E-15, 1.1433E-14, 1.3946E-14, 1.3611E-14,    // F19640
           1.2557E-14, 1.0934E-14, 1.0039E-14, 8.5099E-15, 7.9557E-15,    // F19650
           1.1346E-14, 1.8512E-14, 2.9285E-14, 4.1585E-14, 5.2809E-14,    // F19660
           7.0377E-14, 7.8094E-14, 7.3735E-14, 6.5845E-14, 5.5023E-14,    // F19680
           4.6866E-14, 2.7430E-14, 1.5975E-14, 1.4522E-14, 1.7075E-14,    // F19690
           2.0408E-14, 2.5119E-14, 3.1194E-14, 3.0280E-14, 2.7676E-14,    // F19700
           2.3344E-14, 1.9466E-14, 1.4140E-14, 6.2087E-15, 3.0307E-15,    // F19710
           1.6815E-15, 1.0169E-15, 6.5448E-16, 4.4162E-16, 3.0928E-16,    // F19720
           2.2320E-16, 1.6511E-16, 1.2471E-16, 9.5881E-17, 7.4850E-17,    // F19730
           5.9216E-17, 4.7400E-17, 3.8338E-17, 3.1298E-17, 2.5765E-17,    // F19740
           2.1371E-17, 1.7848E-17, 1.5000E-17, 1.2679E-17, 1.0774E-17,    // F19750
           9.2002E-18, 7.8922E-18, 6.7987E-18, 5.8800E-18, 5.1042E-18,    // F19760
           4.4461E-18, 3.8855E-18, 3.4060E-18, 2.9944E-18, 2.6397E-18,    // F19770
           2.3331E-18};

// CKD_MT 1.00 implementation of N2-N2 model of 
// Borysow, A, and L. Frommhold, 
//  "Collision-induced rototranslational absorption spectra of N2-N2
//  pairs for temperatures from 50 to 300 K", The
//  Astrophysical Journal, 311, 1043-1057, 1986.
// absorption coefficient in units of [cm^-1 Amagat^-2] 
// these data are for T=296K
const Numeric N2N2_CT296_ckd_mt_100_v1  = -10.0;
const Numeric N2N2_CT296_ckd_mt_100_v2  = 350.0;
const Numeric N2N2_CT296_ckd_mt_100_dv  =   5.0;
const int     N2N2_CT296_ckd_mt_100_npt =  73;
const double  N2N2_CT296_ckd_mt_100[N2N2_CT296_ckd_mt_100_npt+addF77fields] = {
           0.0000e0, 
           0.4303E-06, 0.4850E-06, 0.4979E-06, 0.4850E-06, 0.4303E-06,
           0.3715E-06, 0.3292E-06, 0.3086E-06, 0.2920E-06, 0.2813E-06,
           0.2804E-06, 0.2738E-06, 0.2726E-06, 0.2724E-06, 0.2635E-06,
           0.2621E-06, 0.2547E-06, 0.2428E-06, 0.2371E-06, 0.2228E-06,
           0.2100E-06, 0.1991E-06, 0.1822E-06, 0.1697E-06, 0.1555E-06,
           0.1398E-06, 0.1281E-06, 0.1138E-06, 0.1012E-06, 0.9078E-07,
           0.7879E-07, 0.6944E-07, 0.6084E-07, 0.5207E-07, 0.4540E-07,
           0.3897E-07, 0.3313E-07, 0.2852E-07, 0.2413E-07, 0.2045E-07,
           0.1737E-07, 0.1458E-07, 0.1231E-07, 0.1031E-07, 0.8586E-08,
           0.7162E-08, 0.5963E-08, 0.4999E-08, 0.4226E-08, 0.3607E-08,
           0.3090E-08, 0.2669E-08, 0.2325E-08, 0.2024E-08, 0.1783E-08,
           0.1574E-08, 0.1387E-08, 0.1236E-08, 0.1098E-08, 0.9777E-09,
           0.8765E-09, 0.7833E-09, 0.7022E-09, 0.6317E-09, 0.5650E-09,
           0.5100E-09, 0.4572E-09, 0.4115E-09, 0.3721E-09, 0.3339E-09,
           0.3005E-09, 0.2715E-09, 0.2428E-09};


// CKD_MT 1.00 implementation of N2-N2 model of 
// Borysow, A, and L. Frommhold, 
//  "Collision-induced rototranslational absorption spectra of N2-N2
//  pairs for temperatures from 50 to 300 K", The
//  Astrophysical Journal, 311, 1043-1057, 1986.
// absorption coefficient in units of [cm^-1 Amagat^-2] 
// these data are for T=220K
const Numeric N2N2_CT220_ckd_mt_100_v1  = -10.0;
const Numeric N2N2_CT220_ckd_mt_100_v2  = 350.0;
const Numeric N2N2_CT220_ckd_mt_100_dv  =   5.0;
const int     N2N2_CT220_ckd_mt_100_npt =  73;
const double  N2N2_CT220_ckd_mt_100[N2N2_CT220_ckd_mt_100_npt+addF77fields] = {
           0.0000e0,
           0.4946E-06, 0.5756E-06, 0.5964E-06, 0.5756E-06, 0.4946E-06,
           0.4145E-06, 0.3641E-06, 0.3482E-06, 0.3340E-06, 0.3252E-06,
           0.3299E-06, 0.3206E-06, 0.3184E-06, 0.3167E-06, 0.2994E-06,
           0.2943E-06, 0.2794E-06, 0.2582E-06, 0.2468E-06, 0.2237E-06,
           0.2038E-06, 0.1873E-06, 0.1641E-06, 0.1474E-06, 0.1297E-06,
           0.1114E-06, 0.9813E-07, 0.8309E-07, 0.7059E-07, 0.6068E-07,
           0.5008E-07, 0.4221E-07, 0.3537E-07, 0.2885E-07, 0.2407E-07,
           0.1977E-07, 0.1605E-07, 0.1313E-07, 0.1057E-07, 0.8482E-08,
           0.6844E-08, 0.5595E-08, 0.4616E-08, 0.3854E-08, 0.3257E-08,
           0.2757E-08, 0.2372E-08, 0.2039E-08, 0.1767E-08, 0.1548E-08,
           0.1346E-08, 0.1181E-08, 0.1043E-08, 0.9110E-09, 0.8103E-09,
           0.7189E-09, 0.6314E-09, 0.5635E-09, 0.4976E-09, 0.4401E-09,
           0.3926E-09, 0.3477E-09, 0.3085E-09, 0.2745E-09, 0.2416E-09,
           0.2155E-09, 0.1895E-09, 0.1678E-09, 0.1493E-09, 0.1310E-09,
           0.1154E-09, 0.1019E-09, 0.8855E-10};


// CKD_MT 1.00 implementation of N2-N2 model of 
// Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M. Hartmann,
// Infrared collision-induced absorption by N2 near 4.3 microns for
// atmospheric applications: measurements and emprirical modeling, 
// Appl. Optics, 35, 5911-5917, (1996).
const Numeric N2N2_N2F_ckd_mt_100_v1  = 2085.000;
const Numeric N2N2_N2F_ckd_mt_100_v2  = 2670.000;
const Numeric N2N2_N2F_ckd_mt_100_dv  =    5.000;
const int     N2N2_N2F_ckd_mt_100_npt =  118;
const double  N2N2_N2F_ckd_mt_100[N2N2_N2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            0.000E+00,  2.000E-10,  5.200E-09,  1.020E-08,  1.520E-08,  
            2.020E-08,  2.520E-08,  3.020E-08,  4.450E-08,  5.220E-08,  
            6.460E-08,  7.750E-08,  9.030E-08,  1.060E-07,  1.210E-07,  
            1.370E-07,  1.570E-07,  1.750E-07,  2.010E-07,  2.300E-07,  
            2.590E-07,  2.950E-07,  3.260E-07,  3.660E-07,  4.050E-07,  
            4.470E-07,  4.920E-07,  5.340E-07,  5.840E-07,  6.240E-07,  
            6.670E-07,  7.140E-07,  7.260E-07,  7.540E-07,  7.840E-07,  
            8.090E-07,  8.420E-07,  8.620E-07,  8.870E-07,  9.110E-07,  
            9.360E-07,  9.760E-07,  1.030E-06,  1.110E-06,  1.230E-06,  
            1.390E-06,  1.610E-06,  1.760E-06,  1.940E-06,  1.970E-06,  
            1.870E-06,  1.750E-06,  1.560E-06,  1.420E-06,  1.350E-06,  
            1.320E-06,  1.290E-06,  1.290E-06,  1.290E-06,  1.300E-06,  
            1.320E-06,  1.330E-06,  1.340E-06,  1.350E-06,  1.330E-06,  
            1.310E-06,  1.290E-06,  1.240E-06,  1.200E-06,  1.160E-06,  
            1.100E-06,  1.040E-06,  9.960E-07,  9.380E-07,  8.630E-07,  
            7.980E-07,  7.260E-07,  6.550E-07,  5.940E-07,  5.350E-07,  
            4.740E-07,  4.240E-07,  3.770E-07,  3.330E-07,  2.960E-07,  
            2.630E-07,  2.340E-07,  2.080E-07,  1.850E-07,  1.670E-07,  
            1.470E-07,  1.320E-07,  1.200E-07,  1.090E-07,  9.850E-08,  
            9.080E-08,  8.180E-08,  7.560E-08,  6.850E-08,  6.140E-08,  
            5.830E-08,  5.770E-08,  5.000E-08,  4.320E-08,  3.140E-08,  
            2.890E-08,  2.640E-08,  2.390E-08,  2.140E-08,  1.890E-08,  
            1.640E-08,  1.390E-08,  1.140E-08,  8.900E-09,  6.400E-09,  
            3.900E-09,  1.400E-09,  0.000E+00};

//     temperature coefficients:
const double  N2N2_N2Ft_ckd_mt_100[N2N2_N2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            1.040E+03,  1.010E+03,  9.800E+02,  9.500E+02,  9.200E+02,  
            8.900E+02,  8.600E+02,  8.300E+02,  8.020E+02,  7.610E+02,  
            7.220E+02,  6.790E+02,  6.460E+02,  6.090E+02,  5.620E+02,  
            5.110E+02,  4.720E+02,  4.360E+02,  4.060E+02,  3.770E+02,  
            3.550E+02,  3.380E+02,  3.190E+02,  2.990E+02,  2.780E+02,  
            2.550E+02,  2.330E+02,  2.080E+02,  1.840E+02,  1.490E+02,  
            1.070E+02,  6.600E+01,  2.500E+01, -1.300E+01, -4.900E+01,  
           -8.200E+01, -1.040E+02, -1.190E+02, -1.300E+02, -1.390E+02,  
           -1.440E+02, -1.460E+02, -1.460E+02, -1.470E+02, -1.480E+02,  
           -1.500E+02, -1.530E+02, -1.600E+02, -1.690E+02, -1.810E+02,  
           -1.890E+02, -1.950E+02, -2.000E+02, -2.050E+02, -2.090E+02,  
           -2.110E+02, -2.100E+02, -2.100E+02, -2.090E+02, -2.050E+02,  
           -1.990E+02, -1.900E+02, -1.800E+02, -1.680E+02, -1.570E+02,  
           -1.430E+02, -1.260E+02, -1.080E+02, -8.900E+01, -6.300E+01,  
           -3.200E+01,  1.000E+00,  3.500E+01,  6.500E+01,  9.500E+01,  
            1.210E+02,  1.410E+02,  1.520E+02,  1.610E+02,  1.640E+02,  
            1.640E+02,  1.610E+02,  1.550E+02,  1.480E+02,  1.430E+02,  
            1.370E+02,  1.330E+02,  1.310E+02,  1.330E+02,  1.390E+02,  
            1.500E+02,  1.650E+02,  1.870E+02,  2.130E+02,  2.480E+02,  
            2.840E+02,  3.210E+02,  3.720E+02,  4.490E+02,  5.140E+02,  
            5.690E+02,  6.090E+02,  6.420E+02,  6.730E+02,  7.000E+02,  
            7.300E+02,  7.600E+02,  7.900E+02,  8.200E+02,  8.500E+02,  
            8.800E+02,  9.100E+02,  9.400E+02,  9.700E+02,  1.000E+03,  
            1.030E+03,  1.060E+03,  1.090E+03};


// CKD_MT 1.00 implementation of oxygen collision induced fundamental model of 
// F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, 
// J.-M. Hartmann, Ch. Boulet, 
// "Infrared collision-induced absorption by O2 near 6.4 microns for
// atmospheric applications: measurements and emprirical modeling", 
// Appl. Optics, 35, 5911-5917, (1996).
const Numeric O2O2_O2F_ckd_mt_100_v1  = 1340.000;
const Numeric O2O2_O2F_ckd_mt_100_v2  = 1850.000;
const Numeric O2O2_O2F_ckd_mt_100_dv  =    5.000;
const int     O2O2_O2F_ckd_mt_100_npt =  103;
const double  O2O2_O2Fo_ckd_mt_100[O2O2_O2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            0.000E+00,  9.744E-09,  2.256E-08,  3.538E-08,  4.820E-08,  
            6.100E-08,  7.400E-08,  8.400E-08,  9.600E-08,  1.200E-07,  
            1.620E-07,  2.080E-07,  2.460E-07,  2.850E-07,  3.140E-07,  
            3.800E-07,  4.440E-07,  5.000E-07,  5.710E-07,  6.730E-07,  
            7.680E-07,  8.530E-07,  9.660E-07,  1.100E-06,  1.210E-06,  
            1.330E-06,  1.470E-06,  1.590E-06,  1.690E-06,  1.800E-06,  
            1.920E-06,  2.040E-06,  2.150E-06,  2.260E-06,  2.370E-06,  
            2.510E-06,  2.670E-06,  2.850E-06,  3.070E-06,  3.420E-06,  
            3.830E-06,  4.200E-06,  4.450E-06,  4.600E-06,  4.530E-06,  
            4.280E-06,  3.960E-06,  3.680E-06,  3.480E-06,  3.350E-06,  
            3.290E-06,  3.250E-06,  3.230E-06,  3.230E-06,  3.210E-06,  
            3.190E-06,  3.110E-06,  3.030E-06,  2.910E-06,  2.800E-06,  
            2.650E-06,  2.510E-06,  2.320E-06,  2.130E-06,  1.930E-06,  
            1.760E-06,  1.590E-06,  1.420E-06,  1.250E-06,  1.110E-06,  
            9.900E-07,  8.880E-07,  7.910E-07,  6.780E-07,  5.870E-07,  
            5.240E-07,  4.640E-07,  4.030E-07,  3.570E-07,  3.200E-07,  
            2.900E-07,  2.670E-07,  2.420E-07,  2.150E-07,  1.820E-07,  
            1.600E-07,  1.460E-07,  1.280E-07,  1.030E-07,  8.700E-08,  
            8.100E-08,  7.100E-08,  6.400E-08,  5.807E-08,  5.139E-08,  
            4.496E-08,  3.854E-08,  3.212E-08,  2.569E-08,  1.927E-08,  
            1.285E-08,  6.423E-09,  0.000E+00};
                          
const double  O2O2_O2Ft_ckd_mt_100[O2O2_O2F_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  
            4.670E+02,  4.000E+02,  3.150E+02,  3.790E+02,  3.680E+02,  
            4.750E+02,  5.210E+02,  5.310E+02,  5.120E+02,  4.420E+02,  
            4.440E+02,  4.300E+02,  3.810E+02,  3.350E+02,  3.240E+02,  
            2.960E+02,  2.480E+02,  2.150E+02,  1.930E+02,  1.580E+02,  
            1.270E+02,  1.010E+02,  7.100E+01,  3.100E+01, -6.000E+00,  
           -2.600E+01, -4.700E+01, -6.300E+01, -7.900E+01, -8.800E+01,  
           -8.800E+01, -8.700E+01, -9.000E+01, -9.800E+01, -9.900E+01,  
           -1.090E+02, -1.340E+02, -1.600E+02, -1.670E+02, -1.640E+02,  
           -1.580E+02, -1.530E+02, -1.510E+02, -1.560E+02, -1.660E+02,  
           -1.680E+02, -1.730E+02, -1.700E+02, -1.610E+02, -1.450E+02,  
           -1.260E+02, -1.080E+02, -8.400E+01, -5.900E+01, -2.900E+01,  
            4.000E+00,  4.100E+01,  7.300E+01,  9.700E+01,  1.230E+02,  
            1.590E+02,  1.980E+02,  2.200E+02,  2.420E+02,  2.560E+02,  
            2.810E+02,  3.110E+02,  3.340E+02,  3.190E+02,  3.130E+02,  
            3.210E+02,  3.230E+02,  3.100E+02,  3.150E+02,  3.200E+02,  
            3.350E+02,  3.610E+02,  3.780E+02,  3.730E+02,  3.380E+02,  
            3.190E+02,  3.460E+02,  3.220E+02,  2.910E+02,  2.900E+02,  
            3.500E+02,  3.710E+02,  5.040E+02,  4.000E+02,  4.000E+02,  
            4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  4.000E+02,  
            4.000E+02,  4.000E+02,  4.000E+02};



// CKD_MT 1.00 implementation of oxygen v0<-v0 band model of 
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
// Journal of Geophysical Research (1998).
const Numeric O2_00_ckd_mt_100_v1  = 7536.000e0;
const Numeric O2_00_ckd_mt_100_v2  = 8500.000e0;
const Numeric O2_00_ckd_mt_100_dv  =    2.000e0;
const int     O2_00_ckd_mt_100_npt =  483;
const double  O2_00_ckd_mt_100[O2_00_ckd_mt_100_npt+addF77fields] = {
            0.000E+00, 
            0.000E+00,  4.355E-11,  8.709E-11,  1.742E-10,  3.484E-10,  
            6.968E-10,  1.394E-09,  2.787E-09,  3.561E-09,  3.314E-09,  
            3.368E-09,  3.435E-09,  2.855E-09,  3.244E-09,  3.447E-09,  
            3.891E-09,  4.355E-09,  3.709E-09,  4.265E-09,  4.772E-09,  
            4.541E-09,  4.557E-09,  4.915E-09,  4.688E-09,  5.282E-09,  
            5.755E-09,  5.096E-09,  5.027E-09,  4.860E-09,  4.724E-09,  
            5.048E-09,  5.248E-09,  5.473E-09,  4.852E-09,  5.362E-09,  
            6.157E-09,  6.150E-09,  6.347E-09,  6.388E-09,  6.213E-09,  
            6.521E-09,  8.470E-09,  8.236E-09,  8.269E-09,  8.776E-09,  
            9.122E-09,  9.189E-09,  9.778E-09,  8.433E-09,  9.964E-09,  
            9.827E-09,  1.064E-08,  1.063E-08,  1.031E-08,  1.098E-08,  
            1.156E-08,  1.295E-08,  1.326E-08,  1.467E-08,  1.427E-08,  
            1.452E-08,  1.456E-08,  1.554E-08,  1.605E-08,  1.659E-08,  
            1.754E-08,  1.757E-08,  1.876E-08,  1.903E-08,  1.876E-08,  
            1.869E-08,  2.036E-08,  2.203E-08,  2.221E-08,  2.284E-08,  
            2.288E-08,  2.394E-08,  2.509E-08,  2.663E-08,  2.720E-08,  
            2.839E-08,  2.923E-08,  2.893E-08,  2.949E-08,  2.962E-08,  
            3.057E-08,  3.056E-08,  3.364E-08,  3.563E-08,  3.743E-08,  
            3.813E-08,  3.946E-08,  4.082E-08,  4.201E-08,  4.297E-08,  
            4.528E-08,  4.587E-08,  4.704E-08,  4.962E-08,  5.115E-08,  
            5.341E-08,  5.365E-08,  5.557E-08,  5.891E-08,  6.084E-08,  
            6.270E-08,  6.448E-08,  6.622E-08,  6.939E-08,  7.233E-08,  
            7.498E-08,  7.749E-08,  8.027E-08,  8.387E-08,  8.605E-08,  
            8.888E-08,  9.277E-08,  9.523E-08,  9.880E-08,  1.037E-07,  
            1.076E-07,  1.114E-07,  1.151E-07,  1.203E-07,  1.246E-07,  
            1.285E-07,  1.345E-07,  1.408E-07,  1.465E-07,  1.519E-07,  
            1.578E-07,  1.628E-07,  1.685E-07,  1.760E-07,  1.847E-07,  
            1.929E-07,  2.002E-07,  2.070E-07,  2.177E-07,  2.262E-07,  
            2.365E-07,  2.482E-07,  2.587E-07,  2.655E-07,  2.789E-07,  
            2.925E-07,  3.023E-07,  3.153E-07,  3.296E-07,  3.409E-07,  
            3.532E-07,  3.680E-07,  3.859E-07,  3.951E-07,  4.074E-07,  
            4.210E-07,  4.381E-07,  4.588E-07,  4.792E-07,  4.958E-07,  
            5.104E-07,  5.271E-07,  5.501E-07,  5.674E-07,  5.913E-07,  
            6.243E-07,  6.471E-07,  6.622E-07,  6.831E-07,  6.987E-07,  
            7.159E-07,  7.412E-07,  7.698E-07,  7.599E-07,  7.600E-07,  
            7.918E-07,  8.026E-07,  8.051E-07,  8.049E-07,  7.914E-07,  
            7.968E-07,  7.945E-07,  7.861E-07,  7.864E-07,  7.741E-07,  
            7.675E-07,  7.592E-07,  7.400E-07,  7.362E-07,  7.285E-07,  
            7.173E-07,  6.966E-07,  6.744E-07,  6.597E-07,  6.413E-07,  
            6.265E-07,  6.110E-07,  5.929E-07,  5.717E-07,  5.592E-07,  
            5.411E-07,  5.235E-07,  5.061E-07,  4.845E-07,  4.732E-07,  
            4.593E-07,  4.467E-07,  4.328E-07,  4.161E-07,  4.035E-07,  
            3.922E-07,  3.820E-07,  3.707E-07,  3.585E-07,  3.475E-07,  
            3.407E-07,  3.317E-07,  3.226E-07,  3.134E-07,  3.016E-07,  
            2.969E-07,  2.894E-07,  2.814E-07,  2.749E-07,  2.657E-07,  
            2.610E-07,  2.536E-07,  2.467E-07,  2.394E-07,  2.337E-07,  
            2.302E-07,  2.241E-07,  2.191E-07,  2.140E-07,  2.093E-07,  
            2.052E-07,  1.998E-07,  1.963E-07,  1.920E-07,  1.862E-07,  
            1.834E-07,  1.795E-07,  1.745E-07,  1.723E-07,  1.686E-07,  
            1.658E-07,  1.629E-07,  1.595E-07,  1.558E-07,  1.523E-07,  
            1.498E-07,  1.466E-07,  1.452E-07,  1.431E-07,  1.408E-07,  
            1.381E-07,  1.362E-07,  1.320E-07,  1.298E-07,  1.262E-07,  
            1.247E-07,  1.234E-07,  1.221E-07,  1.197E-07,  1.176E-07,  
            1.142E-07,  1.121E-07,  1.099E-07,  1.081E-07,  1.073E-07,  
            1.061E-07,  1.041E-07,  1.019E-07,  9.969E-08,  9.727E-08,  
            9.642E-08,  9.487E-08,  9.318E-08,  9.116E-08,  9.046E-08,  
            8.827E-08,  8.689E-08,  8.433E-08,  8.324E-08,  8.204E-08,  
            8.036E-08,  7.951E-08,  7.804E-08,  7.524E-08,  7.392E-08,  
            7.227E-08,  7.176E-08,  6.975E-08,  6.914E-08,  6.859E-08,  
            6.664E-08,  6.506E-08,  6.368E-08,  6.262E-08,  6.026E-08,  
            6.002E-08,  5.866E-08,  5.867E-08,  5.641E-08,  5.589E-08,  
            5.499E-08,  5.309E-08,  5.188E-08,  5.139E-08,  4.991E-08,  
            4.951E-08,  4.833E-08,  4.640E-08,  4.524E-08,  4.479E-08,  
            4.304E-08,  4.228E-08,  4.251E-08,  4.130E-08,  3.984E-08,  
            3.894E-08,  3.815E-08,  3.732E-08,  3.664E-08,  3.512E-08,  
            3.463E-08,  3.503E-08,  3.218E-08,  3.253E-08,  3.107E-08,  
            2.964E-08,  2.920E-08,  2.888E-08,  2.981E-08,  2.830E-08,  
            2.750E-08,  2.580E-08,  2.528E-08,  2.444E-08,  2.378E-08,  
            2.413E-08,  2.234E-08,  2.316E-08,  2.199E-08,  2.088E-08,  
            1.998E-08,  1.920E-08,  1.942E-08,  1.859E-08,  1.954E-08,  
            1.955E-08,  1.749E-08,  1.720E-08,  1.702E-08,  1.521E-08,  
            1.589E-08,  1.469E-08,  1.471E-08,  1.543E-08,  1.433E-08,  
            1.298E-08,  1.274E-08,  1.226E-08,  1.204E-08,  1.201E-08,  
            1.298E-08,  1.220E-08,  1.220E-08,  1.096E-08,  1.080E-08,  
            9.868E-09,  9.701E-09,  1.130E-08,  9.874E-09,  9.754E-09,  
            9.651E-09,  9.725E-09,  8.413E-09,  7.705E-09,  7.846E-09,  
            8.037E-09,  9.163E-09,  8.098E-09,  8.160E-09,  7.511E-09,  
            7.011E-09,  6.281E-09,  6.502E-09,  7.323E-09,  7.569E-09,  
            5.941E-09,  5.867E-09,  5.676E-09,  4.840E-09,  5.063E-09,  
            5.207E-09,  4.917E-09,  5.033E-09,  5.356E-09,  3.795E-09,  
            4.983E-09,  4.600E-09,  3.635E-09,  3.099E-09,  2.502E-09,  
            3.823E-09,  3.464E-09,  4.332E-09,  3.612E-09,  3.682E-09,  
            3.709E-09,  3.043E-09,  3.593E-09,  3.995E-09,  4.460E-09,  
            3.583E-09,  3.290E-09,  3.132E-09,  2.812E-09,  3.109E-09,  
            3.874E-09,  3.802E-09,  4.024E-09,  3.901E-09,  2.370E-09,  
            1.821E-09,  2.519E-09,  4.701E-09,  3.855E-09,  4.685E-09,  
            5.170E-09,  4.387E-09,  4.148E-09,  4.043E-09,  3.545E-09,  
            3.392E-09,  3.609E-09,  4.635E-09,  3.467E-09,  2.558E-09,  
            3.389E-09,  2.672E-09,  2.468E-09,  1.989E-09,  2.816E-09,  
            4.023E-09,  2.664E-09,  2.219E-09,  3.169E-09,  1.654E-09,  
            3.189E-09,  2.535E-09,  2.618E-09,  3.265E-09,  2.138E-09,  
            1.822E-09,  2.920E-09,  2.002E-09,  1.300E-09,  3.764E-09,  
            3.212E-09,  3.222E-09,  2.961E-09,  2.108E-09,  1.708E-09,  
            2.636E-09,  2.937E-09,  2.939E-09,  2.732E-09,  2.218E-09,  
            1.046E-09,  6.419E-10,  1.842E-09,  1.112E-09,  1.265E-09,  
            4.087E-09,  2.044E-09,  1.022E-09,  5.109E-10,  2.554E-10,  
            1.277E-10,  6.386E-11,  0.000E+00};

#endif // continua_h
